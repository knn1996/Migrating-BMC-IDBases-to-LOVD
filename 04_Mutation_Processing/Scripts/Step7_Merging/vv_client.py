"""
vv_client.py — Consolidated VariantValidator REST client
=========================================================
Submits each variant's hgvs_input to the VV REST API, updates the
standard 23-column merge_tracks schema in-place, and appends responses
to a JSONL cache so repeated runs are free.

Env vars (all required unless noted)
--------------------------------------
IN_TSV       input TSV (merge_tracks 23-col schema; extra cols from
             classify_unresolved output are passed through unchanged)
OUT_TSV      output TSV (same cols as input; status/HGVS cols updated)
CACHE_JSONL  append-only JSONL cache  {"query": ..., "response": ...}
LOG_PATH     plain-text run log
TRANSCRIPTS  VV transcript selector: "mane_select" (default) or "all"
             mane_select — strict: rejects silent HGVS rewrites
             all         — permissive: accepts whichever transcript VV picks

Endpoint verified live 2026-07-21:
    https://rest.variantvalidator.org/VariantValidator/variantvalidator/
    GRCh38/{variant}/{transcripts}
"""

import csv
import json
import os
import re
import sys
import time
import urllib.parse
from pathlib import Path

import requests

IN_TSV      = os.environ["IN_TSV"]
OUT_TSV     = os.environ["OUT_TSV"]
CACHE_JSONL = os.environ["CACHE_JSONL"]
LOG_PATH    = os.environ["LOG_PATH"]
TRANSCRIPTS = os.environ.get("TRANSCRIPTS", "mane_select")

VV_BASE      = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"
VV_BUILD     = "GRCh38"
THROTTLE_SEC = 1.1
TIMEOUT_SEC  = 60
MAX_RETRIES  = 3
BACKOFF_BASE = 5

OUT_FIELDS = [
    "gene", "accession", "allele_num", "sysname", "mut_type", "ref", "alt",
    "chrom", "pos_hg38", "strand", "status", "source_track",
    "nc_hgvs", "c_hgvs", "g_hgvs", "r_hgvs", "p_hgvs",
    "protein_pos_first", "mane_select", "mutalyzer_gene",
    "hgvs_input", "normalized", "errors",
]


def make_logger(log_path):
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    lf = open(log_path, "w", encoding="utf-8")

    def log(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        lf.write(line + "\n")
        lf.flush()

    return log, lf


def load_cache(path):
    cache = {}
    if not Path(path).exists():
        return cache
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                cache[obj["query"]] = obj["response"]
            except (json.JSONDecodeError, KeyError):
                continue
    return cache


def append_cache(path, query, response):
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps({"query": query, "response": response}) + "\n")


def strip_nc_wrapper(hgvs_input):
    m = re.match(
        r"^NC_[0-9.]+\((NM_[0-9.]+)\):(c\..+)$",
        (hgvs_input or "").strip(),
    )
    if m:
        return f"{m.group(1)}:{m.group(2)}"
    return (hgvs_input or "").strip()


def query_vv(hgvs, session, log):
    encoded = urllib.parse.quote(hgvs, safe="")
    url = f"{VV_BASE}/{VV_BUILD}/{encoded}/{TRANSCRIPTS}"
    for attempt in range(MAX_RETRIES):
        try:
            r = session.get(
                url,
                params={"content-type": "application/json"},
                timeout=TIMEOUT_SEC,
            )
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 502, 503, 504):
                wait = int(r.headers.get("Retry-After", BACKOFF_BASE * (2 ** attempt)))
                log(f"    HTTP {r.status_code} — sleeping {wait}s (retry {attempt + 1}/{MAX_RETRIES})")
                time.sleep(wait)
                continue
            return {"_http_error": r.status_code, "_text": r.text[:500]}
        except (requests.RequestException, ValueError) as exc:
            log(f"    Attempt {attempt + 1}/{MAX_RETRIES} exception: {exc}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(BACKOFF_BASE * (attempt + 1))
    return {"_exception": "max retries exceeded"}


def parse_vv_response(resp, vv_query):
    if not isinstance(resp, dict) or "_http_error" in resp or "_exception" in resp:
        return {"status": "invalid", "errors": json.dumps(resp)[:500]}

    flag = resp.get("flag")
    if flag in (None, "warning", "intergenic"):
        warnings = []
        for k, v in resp.items():
            if k in ("flag", "metadata"):
                continue
            if isinstance(v, dict):
                warnings.extend(v.get("validation_warnings", []) or [])
        if not warnings:
            warnings = [str(flag)]
        return {"status": "invalid", "errors": "; ".join(str(w) for w in warnings)[:500]}

    keys = [
        k for k in resp
        if k not in ("flag", "metadata") and isinstance(resp[k], dict)
    ]
    if not keys:
        return {"status": "invalid", "errors": "no_result_block"}

    expected_base = vv_query.split(":", 1)[0].split(".")[0]
    chosen = None
    for k in keys:
        if k.split(".")[0].split(":")[0] == expected_base:
            chosen = resp[k]
            break
    if chosen is None:
        chosen = resp[keys[0]]

    warnings = chosen.get("validation_warnings", []) or []

    if any(
        "Variation has not been validated" in str(w) or "Cannot validate" in str(w)
        for w in warnings
    ):
        return {"status": "invalid", "errors": "; ".join(str(w) for w in warnings)[:500]}

    for w in warnings:
        ws = str(w)
        if "does not agree with reference sequence" in ws:
            return {"status": "invalid", "errors": f"VV_ref_mismatch: {ws[:300]}"}

    hgvs_t = chosen.get("hgvs_transcript_variant", "")

    if TRANSCRIPTS == "mane_select" and hgvs_t:
        submitted_desc = vv_query.split(":", 1)[-1] if ":" in vv_query else vv_query
        returned_desc  = hgvs_t.split(":", 1)[-1]   if ":" in hgvs_t  else hgvs_t
        is_ivs_query   = bool(re.search(r"\bIVS\d+", vv_query, re.IGNORECASE))
        if not is_ivs_query and submitted_desc != returned_desc:
            return {
                "status": "invalid",
                "errors": f"VV_silent_rewrite: submitted={vv_query} returned={hgvs_t}",
            }

    pc      = chosen.get("hgvs_predicted_protein_consequence", {}) or {}
    hgvs_p  = (
        (pc.get("tlr", "") if isinstance(pc, dict) else "")
        or (pc.get("slr", "") if isinstance(pc, dict) else "")
    )
    hgvs_r  = chosen.get("hgvs_rna_variant", "")
    ann     = chosen.get("annotations", {}) or {}
    mane_blk = ann.get("mane_select", {}) if isinstance(ann, dict) else {}
    if not isinstance(mane_blk, dict):
        mane_blk = {}
    mane_tx = mane_blk.get("hgvs_transcript_variant", "") or hgvs_t

    primary = (
        chosen.get("primary_assembly_loci", {}).get("grch38", {})
        or chosen.get("primary_assembly_loci", {}).get("hg38", {})
        or {}
    )
    hgvs_g = primary.get("hgvs_genomic_description", "")
    vcf    = primary.get("vcf", {}) or {}
    chrom  = vcf.get("chr", "")
    pos    = vcf.get("pos", "")

    nc_hgvs = hgvs_g if hgvs_g.startswith("NC_") else ""

    if not hgvs_t and not hgvs_g:
        return {"status": "invalid", "errors": "no_hgvs_returned"}

    return {
        "status":      "ok",
        "c_hgvs":      hgvs_t,
        "g_hgvs":      hgvs_g,
        "r_hgvs":      hgvs_r,
        "p_hgvs":      hgvs_p,
        "nc_hgvs":     nc_hgvs,
        "mane_select": mane_tx,
        "vv_chrom":    f"chr{chrom}" if chrom and not str(chrom).startswith("chr") else str(chrom),
        "vv_pos":      str(pos),
        "errors":      "",
    }


def ping_endpoint(log):
    encoded = urllib.parse.quote("NM_000022.4:c.500del", safe="")
    url = f"{VV_BASE}/{VV_BUILD}/{encoded}/mane_select"
    try:
        r = requests.get(
            url,
            params={"content-type": "application/json"},
            timeout=20,
        )
        if r.status_code == 200 and isinstance(r.json(), dict):
            flag = r.json().get("flag", "?")
            log(f"VV endpoint alive — HTTP 200, flag={flag!r}")
            return True
        log(f"VV endpoint unexpected response — HTTP {r.status_code}")
        return False
    except requests.RequestException as exc:
        log(f"VV endpoint unreachable: {exc}")
        return False


def main():
    for path in (OUT_TSV, CACHE_JSONL, LOG_PATH):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log, log_fh = make_logger(LOG_PATH)
    log("=" * 65)
    log("VariantValidator consolidated client")
    log(f"  endpoint    : {VV_BASE}/{VV_BUILD}/{{variant}}/{TRANSCRIPTS}")
    log(f"  input       : {IN_TSV}")
    log(f"  output      : {OUT_TSV}")
    log(f"  cache       : {CACHE_JSONL}")
    log(f"  throttle    : {THROTTLE_SEC}s")
    log("=" * 65)

    if not ping_endpoint(log):
        log("FATAL: VV endpoint unreachable — aborting")
        log_fh.close()
        sys.exit(1)

    with open(IN_TSV, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    if not rows:
        log("Input is empty — writing empty output and exiting")
        Path(OUT_TSV).write_text("\t".join(OUT_FIELDS) + "\n", encoding="utf-8")
        log_fh.close()
        return

    log(f"Loaded {len(rows)} rows")

    cache = load_cache(CACHE_JSONL)
    log(f"Cache: {len(cache)} existing entries")

    _filter = os.environ.get("FILTER_DISPOSITIONS", "")
    _allowed = {d.strip() for d in _filter.split(",") if d.strip()}
    if _allowed:
        before = len(rows)
        rows = [r for r in rows if r.get("disposition", "") in _allowed]
        log(f"FILTER_DISPOSITIONS={_filter!r}: kept {len(rows)}/{before} rows")

    extra_cols = [c for c in rows[0] if c not in OUT_FIELDS]
    out_fieldnames = OUT_FIELDS + extra_cols

    session = requests.Session()
    session.headers.update({
        "Accept": "application/json",
        "User-Agent": "IDbases-LOVD-migration/1.1",
    })

    out_rows = []
    n_ok = n_fail = n_cached = n_skip = 0
    stopped_early = False

    for i, row in enumerate(rows, 1):
        original_input = (row.get("hgvs_input") or "").strip()
        out_row = {c: row.get(c, "") for c in out_fieldnames}

        if not original_input:
            out_row["status"] = "skip_no_hgvs_input"
            out_rows.append(out_row)
            n_skip += 1
            continue

        vv_query = strip_nc_wrapper(original_input)
        out_row["hgvs_input"] = vv_query

        if vv_query in cache:
            resp = cache[vv_query]
            n_cached += 1
        else:
            log(f"  [{i}/{len(rows)}] {row.get('gene', ''):10s}  {vv_query}")
            resp = query_vv(vv_query, session, log)
            if resp.get("_exception") == "max retries exceeded":
                log("FATAL: max retries reached — flushing remaining rows as 'unreachable'")
                out_row["status"] = "unreachable"
                out_row["errors"] = "VV unreachable after max retries"
                out_rows.append(out_row)
                for remaining in rows[i:]:
                    r2 = {c: remaining.get(c, "") for c in out_fieldnames}
                    r2["hgvs_input"] = strip_nc_wrapper(remaining.get("hgvs_input") or "")
                    r2["status"] = "unreachable"
                    r2["errors"] = "VV unreachable after max retries"
                    out_rows.append(r2)
                stopped_early = True
                break
            append_cache(CACHE_JSONL, vv_query, resp)
            time.sleep(THROTTLE_SEC)

        parsed = parse_vv_response(resp, vv_query)

        if parsed["status"] == "ok":
            n_ok += 1
            out_row["status"]         = "ok"
            out_row["source_track"]   = "VariantValidator"
            out_row["c_hgvs"]         = parsed["c_hgvs"]
            out_row["g_hgvs"]         = parsed["g_hgvs"]
            out_row["r_hgvs"]         = parsed["r_hgvs"]
            out_row["p_hgvs"]         = parsed["p_hgvs"]
            out_row["nc_hgvs"]        = parsed["nc_hgvs"]
            out_row["mane_select"]    = parsed["mane_select"]
            out_row["normalized"]     = parsed["c_hgvs"]
            out_row["errors"]         = ""
            if parsed["vv_chrom"] and not out_row.get("chrom"):
                out_row["chrom"]    = parsed["vv_chrom"]
            if parsed["vv_pos"] and not out_row.get("pos_hg38"):
                out_row["pos_hg38"] = parsed["vv_pos"]
        else:
            n_fail += 1
            out_row["status"] = "invalid"
            out_row["errors"] = parsed.get("errors", "unknown")

        out_rows.append(out_row)

    with open(OUT_TSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    log("")
    log("=== Summary ===")
    log(f"  Input rows          : {len(rows)}")
    log(f"  Skipped (no input)  : {n_skip}")
    log(f"  Served from cache   : {n_cached}")
    log(f"  Resolved (ok)       : {n_ok}")
    log(f"  Still failed        : {n_fail}")
    if stopped_early:
        log("  WARNING: stopped early due to VV unreachability")
    log(f"  Output              : {OUT_TSV}")
    log(f"  Cache               : {CACHE_JSONL}")
    log("=" * 65)
    log_fh.close()


if __name__ == "__main__":
    main()
