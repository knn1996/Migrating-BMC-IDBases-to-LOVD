import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.parse
from pathlib import Path

import requests

THESIS_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
STEP7_OUT  = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output", "Step8_Merging")
LOG_DIR    = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Logs")

VV_BASE = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"
VV_BUILD = "GRCh38"
VV_TRANSCRIPTS = "mane_select"

THROTTLE_SEC = 1.1
TIMEOUT_SEC  = 60
MAX_RETRIES  = 3

OUT_FIELDS = [
    "gene", "accession", "allele_num", "sysname", "mut_type", "ref", "alt",
    "chrom", "pos_hg38", "strand", "status", "source_track",
    "nc_hgvs", "c_hgvs", "g_hgvs", "r_hgvs", "p_hgvs",
    "protein_pos_first", "mane_select", "mutalyzer_gene",
    "hgvs_input", "normalized", "errors",
]


def strip_nc_wrapper(hgvs_input: str) -> str:
    m = re.match(r"^NC_[0-9.]+\((NM_[0-9.]+)\):(c\..+)$", hgvs_input.strip())
    if m:
        return f"{m.group(1)}:{m.group(2)}"
    return hgvs_input.strip()


def query_vv(hgvs: str, session: requests.Session) -> dict:
    encoded = urllib.parse.quote(hgvs, safe="")
    url = f"{VV_BASE}/{VV_BUILD}/{encoded}/{VV_TRANSCRIPTS}"
    params = {"content-type": "application/json"}
    last_exc = None
    for attempt in range(MAX_RETRIES):
        try:
            r = session.get(url, params=params, timeout=TIMEOUT_SEC)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 502, 503, 504):
                time.sleep(2 ** attempt + THROTTLE_SEC)
                continue
            return {"_http_error": r.status_code, "_text": r.text[:500]}
        except (requests.RequestException, ValueError) as e:
            last_exc = e
            time.sleep(2 ** attempt)
    return {"_exception": str(last_exc) if last_exc else "unknown"}


def parse_vv_response(resp: dict, expected_input: str) -> dict:
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

    keys = [k for k in resp if k not in ("flag", "metadata") and isinstance(resp[k], dict)]
    if not keys:
        return {"status": "invalid", "errors": "no_result_block"}

    expected_nm = expected_input.split(":", 1)[0].split(".")[0]
    chosen = None
    for k in keys:
        if k.split(".")[0].split(":")[0] == expected_nm:
            chosen = resp[k]
            break
    if chosen is None:
        chosen = resp[keys[0]]

    submitted = chosen.get("submitted_variant", "")
    if submitted and submitted.split(".")[0] != expected_input.split(".")[0]:
        return {
            "status": "invalid",
            "errors": f"VV_mismatch_bug expected={expected_input} got={submitted}",
        }

    warnings = chosen.get("validation_warnings", []) or []
    if any("Variation has not been validated" in str(w) or "Cannot validate" in str(w)
           for w in warnings):
        return {"status": "invalid", "errors": "; ".join(str(w) for w in warnings)[:500]}

    for w in warnings:
        ws = str(w)
        if "automapped to" in ws or "does not agree with reference sequence" in ws:
            return {"status": "invalid",
                    "errors": f"VV_automap_or_ref_mismatch: {ws[:300]}"}

    hgvs_t      = chosen.get("hgvs_transcript_variant", "")
    if hgvs_t and hgvs_t.split(":", 1)[-1] != expected_input.split(":", 1)[-1]:
        return {"status": "invalid",
                "errors": f"VV_silent_rewrite: submitted={expected_input} returned={hgvs_t}"}

    pc = chosen.get("hgvs_predicted_protein_consequence", {}) or {}
    hgvs_p      = (pc.get("tlr", "") if isinstance(pc, dict) else "") or \
                  (pc.get("slr", "") if isinstance(pc, dict) else "")
    hgvs_r      = chosen.get("hgvs_rna_variant", "")
    annotations = chosen.get("annotations", {}) or {}
    mane_block  = annotations.get("mane_select", {}) if isinstance(annotations, dict) else {}
    if not isinstance(mane_block, dict):
        mane_block = {}
    mane_tx     = mane_block.get("hgvs_transcript_variant", "") or hgvs_t

    primary = chosen.get("primary_assembly_loci", {}).get("grch38", {}) or \
              chosen.get("primary_assembly_loci", {}).get("hg38", {})
    hgvs_g  = primary.get("hgvs_genomic_description", "")
    vcf     = primary.get("vcf", {}) or {}
    chrom   = vcf.get("chr", "")
    pos     = vcf.get("pos", "")
    ref     = vcf.get("ref", "")
    alt     = vcf.get("alt", "")

    nc_hgvs = ""
    if hgvs_g.startswith("NC_"):
        nc_hgvs = hgvs_g

    if not hgvs_t and not hgvs_g:
        return {"status": "invalid", "errors": "no_hgvs_returned"}

    return {
        "status": "ok",
        "c_hgvs": hgvs_t,
        "g_hgvs": hgvs_g,
        "r_hgvs": hgvs_r,
        "p_hgvs": hgvs_p,
        "nc_hgvs": nc_hgvs,
        "mane_select": mane_tx,
        "vv_chrom": chrom,
        "vv_pos": pos,
        "vv_ref": ref,
        "vv_alt": alt,
        "errors": "",
    }


def load_cache(cache_path: Path) -> dict:
    cache = {}
    if cache_path.exists():
        with cache_path.open(encoding="utf-8") as f:
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


def append_cache(cache_path: Path, query: str, response: dict) -> None:
    with cache_path.open("a", encoding="utf-8") as f:
        f.write(json.dumps({"query": query, "response": response}) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unresolved",
                    default=os.path.join(STEP7_OUT, "unresolved_variants.tsv"))
    ap.add_argument("--out",
                    default=os.path.join(STEP7_OUT, "eintronic_recovery.tsv"))
    ap.add_argument("--cache",
                    default=os.path.join(STEP7_OUT, "vv_cache.jsonl"))
    ap.add_argument("--log",
                    default=os.path.join(LOG_DIR, "recover_eintronic.log"))
    ap.add_argument("--limit", type=int, default=0,
                    help="Process only N rows (0=all). For testing.")
    args = ap.parse_args()

    Path(os.path.dirname(args.out)).mkdir(parents=True, exist_ok=True)
    Path(os.path.dirname(args.log)).mkdir(parents=True, exist_ok=True)
    Path(os.path.dirname(args.cache)).mkdir(parents=True, exist_ok=True)

    log_f = open(args.log, "w", encoding="utf-8")
    def log(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        log_f.write(line + "\n")
        log_f.flush()

    with open(args.unresolved, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    eintronic_rows = [r for r in rows if "EINTRONIC" in (r.get("errors") or "")]
    log(f"Found {len(eintronic_rows)} EINTRONIC rows out of {len(rows)} unresolved")

    if args.limit:
        eintronic_rows = eintronic_rows[: args.limit]
        log(f"Limiting to first {len(eintronic_rows)} for testing")

    cache_path = Path(args.cache)
    cache = load_cache(cache_path)
    log(f"Loaded {len(cache)} cached responses from {cache_path}")

    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "IDbases-LOVD-migration/1.0"})

    out_rows = []
    n_ok = 0
    n_fail = 0
    n_cached = 0

    for i, row in enumerate(eintronic_rows, 1):
        original_input = (row.get("hgvs_input") or "").strip()
        vv_query = strip_nc_wrapper(original_input)

        if vv_query in cache:
            resp = cache[vv_query]
            n_cached += 1
        else:
            resp = query_vv(vv_query, session)
            append_cache(cache_path, vv_query, resp)
            time.sleep(THROTTLE_SEC)

        parsed = parse_vv_response(resp, vv_query)

        new_row = {k: row.get(k, "") for k in OUT_FIELDS}
        new_row["source_track"] = "VariantValidator_EINTRONIC_recovery"
        new_row["hgvs_input"]   = vv_query

        if parsed["status"] == "ok":
            n_ok += 1
            new_row["status"]      = "ok"
            new_row["c_hgvs"]      = parsed["c_hgvs"]
            new_row["g_hgvs"]      = parsed["g_hgvs"]
            new_row["r_hgvs"]      = parsed["r_hgvs"]
            new_row["p_hgvs"]      = parsed["p_hgvs"]
            new_row["nc_hgvs"]     = parsed["nc_hgvs"]
            new_row["mane_select"] = parsed["mane_select"]
            new_row["normalized"]  = parsed["c_hgvs"]
            new_row["errors"]      = ""
            if parsed["vv_chrom"] and not new_row.get("chrom"):
                new_row["chrom"] = f"chr{parsed['vv_chrom']}"
            if parsed["vv_pos"] and not new_row.get("pos_hg38"):
                new_row["pos_hg38"] = str(parsed["vv_pos"])
            log(f"  [{i}/{len(eintronic_rows)}] OK  {row['gene']:8s} {vv_query} -> {parsed['c_hgvs']}")
        else:
            n_fail += 1
            new_row["status"] = "invalid"
            new_row["errors"] = parsed.get("errors", "unknown")
            log(f"  [{i}/{len(eintronic_rows)}] FAIL {row['gene']:8s} {vv_query} -> {parsed['errors'][:120]}")

        out_rows.append(new_row)

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    log("")
    log(f"=== Summary ===")
    log(f"Total submitted:  {len(eintronic_rows)}")
    log(f"From cache:       {n_cached}")
    log(f"Recovered (ok):   {n_ok}")
    log(f"Still failed:     {n_fail}")
    log(f"Output:           {args.out}")
    log(f"Cache:            {args.cache}")
    log(f"Log:              {args.log}")
    log_f.close()


if __name__ == "__main__":
    sys.exit(main())
