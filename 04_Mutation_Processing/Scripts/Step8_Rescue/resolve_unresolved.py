"""
resolve_unresolved.py
=====================
Aggregates mutalyzer_rescue + vv_rescue results and applies a ±1
c.-position offset fix for ESEQUENCEMISMATCH/positional_offset cases.

Resolution priority per variant
--------------------------------
ENOSELECTORFOUND   mutalyzer_rescue ok  →  mane_remap
                   vv_rescue ok         →  variantvalidator
EINTRONIC          vv_rescue ok         →  variantvalidator
positional_offset  vv_rescue ok         →  variantvalidator
                   ±1 c. shift + VV ok  →  offset_fix
RESCUABLE_MANUAL   not attempted (needs manual HGVS rewrite first)
UNRESCUABLE        not attempted

Env vars
--------
IN_DISPOSITION   unresolved_disposition.tsv  (classify_unresolved output)
IN_MUTALYZER     rescue_mutalyzer_results.tsv (mutalyzer_rescue output)
IN_VV            vv_rescue_results.tsv        (vv_rescue output)
OFFSET_CSV       lrg_offset_results.csv       (refcheck_offset output)
MERGED_TSV       merged_variants.tsv          (pre-rescue resolved count)
OUT_RESOLVED     resolved_unresolved.tsv
OUT_FUNNEL       resolution_funnel.tsv
OUT_NARRATIVE    resolution_narrative.md
VV_CACHE_JSONL   shared VV JSONL cache (same file used by vv_rescue)
LOG_PATH         plain-text run log
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hgvs_norm import normalize_hgvs_query, describe_normalization

IN_DISPOSITION = os.environ["IN_DISPOSITION"]
IN_MUTALYZER   = os.environ["IN_MUTALYZER"]
IN_VV          = os.environ["IN_VV"]
OFFSET_CSV     = os.environ["OFFSET_CSV"]
MERGED_TSV     = os.environ["MERGED_TSV"]
OUT_RESOLVED   = os.environ["OUT_RESOLVED"]
OUT_FUNNEL     = os.environ["OUT_FUNNEL"]
OUT_NARRATIVE  = os.environ["OUT_NARRATIVE"]
VV_CACHE_JSONL = os.environ["VV_CACHE_JSONL"]
LOG_PATH       = os.environ["LOG_PATH"]

for _p in (OUT_RESOLVED, OUT_FUNNEL, OUT_NARRATIVE, VV_CACHE_JSONL, LOG_PATH):
    Path(_p).parent.mkdir(parents=True, exist_ok=True)

VV_BASE      = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"
VV_BUILD     = "GRCh38"
VV_TX        = "mane_select"
THROTTLE_SEC = 1.1
TIMEOUT_SEC  = 60
MAX_RETRIES  = 3
BACKOFF_BASE = 5

OUT_FIELDS  = [
    "gene", "accession", "allele_num", "sysname", "mut_type", "ref", "alt",
    "chrom", "pos_hg38", "strand", "status", "source_track",
    "nc_hgvs", "c_hgvs", "g_hgvs", "r_hgvs", "p_hgvs",
    "protein_pos_first", "mane_select", "mutalyzer_gene",
    "hgvs_input", "normalized", "errors",
]
EXTRA_FIELDS = ["validated_hgvs", "transcript_used", "resolving_tool", "accession_list",
                "normalization_applied", "needs_curator_review"]

_SIMPLE_SUB = re.compile(r"^c\.(-?\d+)([ACGT]+>[ACGT]+)$", re.IGNORECASE)


def make_logger(path):
    lf = open(path, "w", encoding="utf-8")
    def log(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        lf.write(line + "\n"); lf.flush()
    return log, lf


def load_cache(path):
    cache = {}
    p = Path(path)
    if not p.exists():
        return cache
    with p.open(encoding="utf-8") as f:
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


def _parse_vv(resp, vv_query):
    if not isinstance(resp, dict) or "_http_error" in resp or "_exception" in resp:
        return None
    flag = resp.get("flag")
    if flag in (None, "warning", "intergenic"):
        return None
    keys = [k for k in resp if k not in ("flag", "metadata") and isinstance(resp[k], dict)]
    if not keys:
        return None
    base = vv_query.split(":", 1)[0].split(".")[0]
    chosen = next((resp[k] for k in keys if k.split(".")[0].split(":")[0] == base), resp[keys[0]])
    warns = chosen.get("validation_warnings", []) or []
    if any("Variation has not been validated" in str(w) or "Cannot validate" in str(w) for w in warns):
        return None
    for w in warns:
        if "does not agree with reference sequence" in str(w):
            return None
    hgvs_t = chosen.get("hgvs_transcript_variant", "")
    if hgvs_t:
        sub_d = vv_query.split(":", 1)[-1]
        ret_d = hgvs_t.split(":", 1)[-1]
        is_ivs = bool(re.search(r"\bIVS\d+", vv_query, re.IGNORECASE))
        if not is_ivs and sub_d != ret_d:
            return None
    pc     = chosen.get("hgvs_predicted_protein_consequence", {}) or {}
    hgvs_p = (pc.get("tlr", "") if isinstance(pc, dict) else "") or \
              (pc.get("slr", "") if isinstance(pc, dict) else "")
    hgvs_r = chosen.get("hgvs_rna_variant", "")
    ann    = chosen.get("annotations", {}) or {}
    mane_b = ann.get("mane_select", {}) if isinstance(ann, dict) else {}
    mane_t = (mane_b.get("hgvs_transcript_variant", "") if isinstance(mane_b, dict) else "") or hgvs_t
    prim   = chosen.get("primary_assembly_loci", {})
    loci   = prim.get("grch38", {}) or prim.get("hg38", {}) or {}
    hgvs_g = loci.get("hgvs_genomic_description", "")
    vcf    = loci.get("vcf", {}) or {}
    chrom  = vcf.get("chr", "")
    pos    = vcf.get("pos", "")
    nc     = hgvs_g if hgvs_g.startswith("NC_") else ""
    if not hgvs_t and not hgvs_g:
        return None
    return {
        "c_hgvs": hgvs_t, "g_hgvs": hgvs_g, "r_hgvs": hgvs_r,
        "p_hgvs": hgvs_p, "nc_hgvs": nc, "mane_select": mane_t,
        "vv_chrom": f"chr{chrom}" if chrom and not str(chrom).startswith("chr") else str(chrom),
        "vv_pos": str(pos),
    }


def call_vv(raw_hgvs, session, cache, log):
    hgvs = normalize_hgvs_query(raw_hgvs)
    if not hgvs:
        return None
    if hgvs != (raw_hgvs or "").strip():
        log(f"    normalised: {raw_hgvs} -> {hgvs}")
    if hgvs in cache:
        return _parse_vv(cache[hgvs], hgvs)
    encoded = urllib.parse.quote(hgvs, safe="")
    url = f"{VV_BASE}/{VV_BUILD}/{encoded}/{VV_TX}"
    for attempt in range(MAX_RETRIES):
        try:
            r = session.get(url, params={"content-type": "application/json"}, timeout=TIMEOUT_SEC)
            if r.status_code == 200:
                resp = r.json()
                append_cache(VV_CACHE_JSONL, hgvs, resp)
                time.sleep(THROTTLE_SEC)
                return _parse_vv(resp, hgvs)
            if r.status_code in (429, 502, 503, 504):
                wait = BACKOFF_BASE * (2 ** attempt)
                log(f"    HTTP {r.status_code} — sleeping {wait}s")
                time.sleep(wait)
                continue
            resp = {"_http_error": r.status_code}
            append_cache(VV_CACHE_JSONL, hgvs, resp)
            return None
        except requests.RequestException as exc:
            log(f"    VV request failed: {exc}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(BACKOFF_BASE * (attempt + 1))
    return None


def try_offset_shifts(hgvs_input, session, cache, log):
    hgvs_input = normalize_hgvs_query(hgvs_input)
    if ":" not in hgvs_input:
        return None, None
    prefix, c_part = hgvs_input.rsplit(":", 1)
    m = _SIMPLE_SUB.match(c_part.strip())
    if not m:
        return None, None
    pos  = int(m.group(1))
    rest = m.group(2)
    for delta in (-1, +1):
        new_pos = pos + delta
        if new_pos == 0:
            continue
        shifted = f"{prefix}:c.{new_pos}{rest}"
        parsed  = call_vv(shifted, session, cache, log)
        if parsed is not None:
            return parsed, shifted
    return None, None


def _build_row(base, vv, resolving_tool, new_hgvs_input=None):
    row = {c: base.get(c, "") for c in OUT_FIELDS + EXTRA_FIELDS}
    row["status"]       = "ok"
    row["source_track"] = resolving_tool
    for f in ("c_hgvs", "g_hgvs", "r_hgvs", "p_hgvs", "nc_hgvs", "mane_select"):
        row[f] = vv.get(f, "") or row[f]
    row["normalized"]      = vv.get("c_hgvs", "") or row.get("normalized", "")
    row["errors"]          = ""
    if vv.get("vv_chrom") and not row.get("chrom"):
        row["chrom"]    = vv["vv_chrom"]
    if vv.get("vv_pos") and not row.get("pos_hg38"):
        row["pos_hg38"] = vv["vv_pos"]
    if new_hgvs_input:
        row["hgvs_input"] = new_hgvs_input
    validated = row["c_hgvs"] or row["g_hgvs"] or row["normalized"]
    tx = ""
    if row["c_hgvs"] and ":" in row["c_hgvs"]:
        tx = row["c_hgvs"].split(":", 1)[0]
    elif row["mane_select"] and ":" in row["mane_select"]:
        tx = row["mane_select"].split(":", 1)[0]
    row["validated_hgvs"]  = validated
    row["transcript_used"] = tx
    row["resolving_tool"]  = resolving_tool
    norm = describe_normalization(base.get("hgvs_input", ""))
    row["normalization_applied"] = norm
    row["needs_curator_review"]  = "yes" if "deldup_shortened" in norm else ""
    return row


def main():
    log, log_fh = make_logger(LOG_PATH)
    log("=" * 65)
    log("resolve_unresolved.py")
    log(f"  disposition  : {IN_DISPOSITION}")
    log(f"  mutalyzer    : {IN_MUTALYZER}")
    log(f"  vv results   : {IN_VV}")
    log(f"  offset csv   : {OFFSET_CSV}")
    log(f"  vv cache     : {VV_CACHE_JSONL}")
    log("=" * 65)

    def _rows(path):
        with open(path, encoding="utf-8") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    disp_rows = _rows(IN_DISPOSITION)
    mut_rows  = _rows(IN_MUTALYZER)
    vv_rows   = _rows(IN_VV)

    mut_ok = {
        (r["gene"].upper(), r["sysname"].strip()): r
        for r in mut_rows
        if r.get("status") == "ok"
    }
    vv_ok = {
        (r["gene"].upper(), r["sysname"].strip()): r
        for r in vv_rows
        if r.get("status") == "ok"
    }

    offset_index = {}
    with open(OFFSET_CSV, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            gene = r.get("gene", "").strip().upper()
            if r.get("status", "").startswith("ok"):
                offset_index[gene] = {
                    "sstart": r.get("sstart", ""),
                    "strand": r.get("strand", ""),
                }

    n_merged_rows = sum(1 for _ in open(MERGED_TSV, encoding="utf-8")) - 1

    vv_cache = load_cache(VV_CACHE_JSONL)
    log(f"VV cache pre-loaded: {len(vv_cache)} entries")

    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "IDbases-LOVD-migration/1.1"})

    resolved = []
    residual = []
    counters  = {"mane_remap": 0, "variantvalidator": 0, "offset_fix": 0, "vv_mismatch": 0}
    VV_PROBE_SUBS = {"single_base_wrong", "single_base_complement",
                     "multibase_rc_match", "completely_different", "length_mismatch"}
    residual_cats = {}

    for row in disp_rows:
        key  = (row["gene"].upper(), row.get("sysname", "").strip())
        cat  = row.get("category", "")
        sub  = row.get("subcategory", "")
        disp = row.get("disposition", "")

        if disp in ("UNRESCUABLE", "UNRESCUABLE_DATA_ERROR") and not (
            cat == "ESEQUENCEMISMATCH" and sub in VV_PROBE_SUBS
        ):
            residual_cats[f"{cat}/{sub}"] = residual_cats.get(f"{cat}/{sub}", 0) + 1
            residual.append(row)
            continue

        if disp == "RESCUABLE_MANUAL":
            residual_cats[f"{cat}/{sub}"] = residual_cats.get(f"{cat}/{sub}", 0) + 1
            residual.append(row)
            continue

        if disp == "UNKNOWN":
            residual_cats["OTHER/uncategorised"] = residual_cats.get("OTHER/uncategorised", 0) + 1
            residual.append(row)
            continue

        if cat == "NO_REF":
            residual_cats["NO_REF/no_usable_reference"] = residual_cats.get("NO_REF/no_usable_reference", 0) + 1
            residual.append(row)
            continue

        resolved_row = None

        if cat == "ENOSELECTORFOUND":
            if key in mut_ok:
                r = mut_ok[key]
                vv_data = {
                    "c_hgvs": r.get("c_hgvs", "") or r.get("normalized", ""),
                    "g_hgvs": r.get("g_hgvs", ""),
                    "r_hgvs": r.get("r_hgvs", ""), "p_hgvs": r.get("p_hgvs", ""),
                    "nc_hgvs": r.get("nc_hgvs", ""),
                    "mane_select": r.get("mane_select", "") or r.get("normalized", ""),
                    "vv_chrom": r.get("chrom", ""), "vv_pos": r.get("pos_hg38", ""),
                }
                anchor_q = r.get("normalized", "") or r.get("c_hgvs", "")
                if anchor_q and not vv_data["g_hgvs"] and not vv_data["vv_pos"]:
                    anc = call_vv(anchor_q, session, vv_cache, log)
                    if anc:
                        vv_data["g_hgvs"]   = vv_data["g_hgvs"]   or anc.get("g_hgvs", "")
                        vv_data["nc_hgvs"]  = vv_data["nc_hgvs"]  or anc.get("nc_hgvs", "")
                        vv_data["vv_chrom"] = vv_data["vv_chrom"] or anc.get("vv_chrom", "")
                        vv_data["vv_pos"]   = vv_data["vv_pos"]   or anc.get("vv_pos", "")
                        vv_data["r_hgvs"]   = vv_data["r_hgvs"]   or anc.get("r_hgvs", "")
                        vv_data["p_hgvs"]   = vv_data["p_hgvs"]   or anc.get("p_hgvs", "")
                        vv_data["mane_select"] = vv_data["mane_select"] or anc.get("mane_select", "")
                resolved_row = _build_row(row, vv_data, "mane_remap",
                                          new_hgvs_input=r.get("hgvs_input", ""))
                counters["mane_remap"] += 1
            elif key in vv_ok:
                r = vv_ok[key]
                vv_data = {
                    "c_hgvs": r.get("c_hgvs", ""), "g_hgvs": r.get("g_hgvs", ""),
                    "r_hgvs": r.get("r_hgvs", ""), "p_hgvs": r.get("p_hgvs", ""),
                    "nc_hgvs": r.get("nc_hgvs", ""), "mane_select": r.get("mane_select", ""),
                    "vv_chrom": r.get("chrom", ""), "vv_pos": r.get("pos_hg38", ""),
                }
                resolved_row = _build_row(row, vv_data, "variantvalidator")
                counters["variantvalidator"] += 1

        elif cat == "EINTRONIC":
            if key in vv_ok:
                r = vv_ok[key]
                vv_data = {
                    "c_hgvs": r.get("c_hgvs", ""), "g_hgvs": r.get("g_hgvs", ""),
                    "r_hgvs": r.get("r_hgvs", ""), "p_hgvs": r.get("p_hgvs", ""),
                    "nc_hgvs": r.get("nc_hgvs", ""), "mane_select": r.get("mane_select", ""),
                    "vv_chrom": r.get("chrom", ""), "vv_pos": r.get("pos_hg38", ""),
                }
                resolved_row = _build_row(row, vv_data, "variantvalidator")
                counters["variantvalidator"] += 1

        elif cat == "ESEQUENCEMISMATCH" and sub == "positional_offset":
            if key in vv_ok:
                r = vv_ok[key]
                vv_data = {
                    "c_hgvs": r.get("c_hgvs", ""), "g_hgvs": r.get("g_hgvs", ""),
                    "r_hgvs": r.get("r_hgvs", ""), "p_hgvs": r.get("p_hgvs", ""),
                    "nc_hgvs": r.get("nc_hgvs", ""), "mane_select": r.get("mane_select", ""),
                    "vv_chrom": r.get("chrom", ""), "vv_pos": r.get("pos_hg38", ""),
                }
                resolved_row = _build_row(row, vv_data, "variantvalidator")
                counters["variantvalidator"] += 1
            else:
                hgvs_in   = row.get("hgvs_input", "").strip()
                gene_name = row["gene"].upper()
                sstart_ev = offset_index.get(gene_name, {}).get("sstart", "")
                log(f"  offset_fix attempt: {gene_name} {hgvs_in}  sstart={sstart_ev}")
                parsed, shifted = try_offset_shifts(hgvs_in, session, vv_cache, log)
                if parsed is not None:
                    resolved_row = _build_row(row, parsed, "offset_fix",
                                              new_hgvs_input=shifted)
                    counters["offset_fix"] += 1
                    log(f"    offset_fix OK: {shifted} -> {parsed.get('c_hgvs', '')}")
                else:
                    log(f"    offset_fix FAIL: {hgvs_in}")

        elif cat == "ESEQUENCEMISMATCH" and sub in VV_PROBE_SUBS:
            r = vv_ok.get(key)
            if r is not None:
                vv_data = {
                    "c_hgvs": r.get("c_hgvs", ""), "g_hgvs": r.get("g_hgvs", ""),
                    "r_hgvs": r.get("r_hgvs", ""), "p_hgvs": r.get("p_hgvs", ""),
                    "nc_hgvs": r.get("nc_hgvs", ""), "mane_select": r.get("mane_select", ""),
                    "vv_chrom": r.get("chrom", ""), "vv_pos": r.get("pos_hg38", ""),
                }
                resolved_row = _build_row(row, vv_data, "variantvalidator")
                counters["vv_mismatch"] += 1
            else:
                hgvs_in = normalize_hgvs_query(row.get("hgvs_input", ""))
                parsed = call_vv(hgvs_in, session, vv_cache, log) if hgvs_in else None
                if parsed is not None and parsed.get("c_hgvs"):
                    resolved_row = _build_row(row, parsed, "variantvalidator",
                                              new_hgvs_input=hgvs_in)
                    counters["vv_mismatch"] += 1
                    log(f"    VV mismatch-probe OK: {hgvs_in} -> {parsed.get('c_hgvs','')}")
                else:
                    log(f"    VV mismatch-probe FAIL: {hgvs_in}")

        if resolved_row is not None:
            resolved.append(resolved_row)
            log(f"  RESOLVED [{resolved_row['resolving_tool']}] {row['gene']} {row.get('sysname','')} "
                f"-> {resolved_row.get('validated_hgvs', '')[:60]}")
        else:
            residual_cats[f"{cat}/{sub}"] = residual_cats.get(f"{cat}/{sub}", 0) + 1
            residual.append(row)

    with open(OUT_RESOLVED, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUT_FIELDS + EXTRA_FIELDS,
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(resolved)

    n_unresolved  = len(disp_rows)
    n_rescued     = len(resolved)
    n_remaining   = len(residual)
    n_mane        = counters["mane_remap"]
    n_vv          = counters["variantvalidator"]
    n_vv_mismatch = counters["vv_mismatch"]
    n_offset      = counters["offset_fix"]
    n_review      = sum(1 for r in resolved if r.get("needs_curator_review") == "yes")

    pct_rescued = n_rescued / n_unresolved * 100 if n_unresolved else 0

    REF_MISMATCH_SUBS = {
        "single_base_wrong", "single_base_complement", "completely_different",
        "multibase_rc_match", "length_mismatch", "positional_offset",
    }
    TX_STRUCTURE_SUBS = {"nm_not_in_ng", "intronic_c_notation"}

    n_ref_mismatch = n_tx_structure = n_curation_loss = n_manual_only = 0
    for r in residual:
        sub  = r.get("subcategory", "")
        disp = r.get("disposition", "")
        if sub in REF_MISMATCH_SUBS:
            n_ref_mismatch += 1
        elif sub in TX_STRUCTURE_SUBS:
            n_tx_structure += 1
        elif disp == "RESCUABLE_MANUAL":
            n_manual_only += 1
        else:
            n_curation_loss += 1

    with open(OUT_FUNNEL, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["stage", "count"])
        w.writerow(["unresolved_in",                   n_unresolved])
        w.writerow(["resolved_by_MANE_remap",          n_mane])
        w.writerow(["resolved_by_VariantValidator",    n_vv])
        w.writerow(["resolved_by_VV_mismatch_probe",   n_vv_mismatch])
        w.writerow(["resolved_by_offset_fix",          n_offset])
        w.writerow(["total_rescued",                   n_rescued])
        w.writerow(["  of which need curator review",  n_review])
        w.writerow(["remaining_unresolved",            n_remaining])
        w.writerow(["  ref_mismatch",                  n_ref_mismatch])
        w.writerow(["  tx_structure",                  n_tx_structure])
        w.writerow(["  curation_loss",                 n_curation_loss])
        w.writerow(["  manual_rewrite",                n_manual_only])
        w.writerow(["", ""])
        w.writerow(["remaining_by_category", "count"])
        for cat_sub, n in sorted(residual_cats.items(), key=lambda x: -x[1]):
            w.writerow([cat_sub, n])

    narrative = f"""\
## Resolution of Unresolved Variants

Of the {n_unresolved} distinct variants that failed all three Mutalyzer
processing tracks (NG\\_IDRefseq, NM\\_MANE, and NM\\_IDRefseq), **{n_rescued}
({pct_rescued:.1f}%)** were automatically resolved by the rescue layer:
{n_mane} via MANE Select NM\\_ re-mapping (correcting obsolete transcript
version references), {n_vv} via VariantValidator against the current GRCh38
build (intronic and IVS-notation variants), {n_vv_mismatch} via a
VariantValidator probe of reference-mismatch candidates, and
{n_offset} via empirical coordinate-offset correction (±1 c.\\ position
shift confirmed against GRCh38).

Of these, **{n_review} require curator review before submission.** These are
descriptions that spelled out the deleted or duplicated sequence in the
deprecated pre-2016 style (e.g. c.1293delC rather than c.1293del). Shortening
the description discards the source record's assertion about which base was
deleted while retaining its assertion about the position, so where that base
disagreed with GRCh38 the variant now validates against the reference base
instead. Position is being trusted over allele; the affected rows are flagged
`needs_curator_review` in `resolved_unresolved.tsv` and should be confirmed
against the primary literature before deposition.

Resolution rates are reported on a distinct-variant basis in
`rescue_metrics_comparison.csv`; patient-level row counts before and after
augmentation are reported in `augment_stats.txt` and `dedup_rescued_stats.txt`.
The pre-rescue pool contained {n_merged_rows} patient-variant rows, which is
not directly comparable with the {n_unresolved} distinct unresolved variants
above and is therefore not combined into a single ratio here.

The residual **{n_remaining}** variants fall into three mechanistically
distinct groups:

- **{n_ref_mismatch} reference-sequence mismatches.** The base asserted by the
  source record disagrees with GRCh38 at the stated position, either through
  strand-orientation error in the original curation or through a legacy
  reference sequence that diverges from the current assembly. These cannot be
  corrected without fabricating sequence information and represent a genuine
  limit of automated rescue.
- **{n_tx_structure} legacy transcript-structure incompatibilities.** The
  historical c.\\ numbering anchors positions — particularly intronic offsets —
  to exon boundaries that do not exist in the MANE Select transcript, yielding
  EEXONBOUNDARY on re-validation. These are not source-data errors; resolving
  them requires reconstructing genomic coordinates from the legacy RefSeqGene
  reference rather than re-mapping the c.\\ description, which is outside the
  scope of the present rescue layer.
- **{n_curation_loss} curation losses and {n_manual_only} manual-rewrite cases.**
  Descriptions that lost their inserted sequence during original database entry,
  carry reversed or non-consecutive ranges, or otherwise require a targeted
  human HGVS rewrite before automated re-verification is meaningful.
"""

    with open(OUT_NARRATIVE, "w", encoding="utf-8") as f:
        f.write(narrative)

    log("")
    log("=== Resolution summary ===")
    log(f"  Unresolved input       : {n_unresolved}")
    log(f"  Resolved by mane_remap : {n_mane}")
    log(f"  Resolved by VV         : {n_vv}")
    log(f"  Resolved by VV mismatch: {n_vv_mismatch}")
    log(f"  Resolved by offset_fix : {n_offset}")
    log(f"  Total rescued          : {n_rescued}  ({pct_rescued:.1f}% of unresolved)")
    log(f"  Remaining              : {n_remaining}")
    log(f"    ref mismatch         : {n_ref_mismatch}")
    log(f"    tx structure         : {n_tx_structure}")
    log(f"    curation loss        : {n_curation_loss}")
    log(f"    manual rewrite       : {n_manual_only}")
    log(f"  Pre-rescue merged rows : {n_merged_rows}  (patient-variant rows, not variants)")
    log("  Resolution rates       : see rescue_metrics_comparison.csv")
    log(f"  OUT_RESOLVED           : {OUT_RESOLVED}")
    log(f"  OUT_FUNNEL             : {OUT_FUNNEL}")
    log(f"  OUT_NARRATIVE          : {OUT_NARRATIVE}")
    log("=" * 65)
    log_fh.close()


if __name__ == "__main__":
    main()
