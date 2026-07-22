"""
build_final_lovd.py
===================
Adds a provenance column to the rescue-deduped flat, emits lovd_flat_final.tsv,
computes published vs final metrics, and enforces four acceptance checks.

Env vars
--------
RESCUED_DEDUP_TSV     Step8_Rescue/dedup_merged_variants.rescued.tsv
PRIMARY_DEDUP_TSV     Step7_Merging/dedup_merged_variants.tsv
MERGED_TSV            Step7_Merging/merged_variants.tsv
UNRESOLVED_DISP_TSV   Step8_Rescue/unresolved_disposition.tsv
SCRIPTS_DIR           Scripts/ root (used to locate dedup script and git cwd)
OUT_FLAT              Step8_Rescue/lovd_flat_final.tsv
OUT_METRICS           Step8_Rescue/rescue_metrics_comparison.csv
"""

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

RESCUED_DEDUP_TSV   = os.environ["RESCUED_DEDUP_TSV"]
PRIMARY_DEDUP_TSV   = os.environ["PRIMARY_DEDUP_TSV"]
MERGED_TSV          = os.environ["MERGED_TSV"]
UNRESOLVED_DISP_TSV = os.environ["UNRESOLVED_DISP_TSV"]
SCRIPTS_DIR         = os.environ["SCRIPTS_DIR"]
OUT_FLAT            = os.environ["OUT_FLAT"]
OUT_METRICS         = os.environ["OUT_METRICS"]

for _p in (OUT_FLAT, OUT_METRICS):
    Path(_p).parent.mkdir(parents=True, exist_ok=True)


def _read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _provenance(source_track):
    st = (source_track or "").strip()
    if st == "mane_remap":
        return "rescued-mutalyzer"
    if st in ("variantvalidator", "offset_fix"):
        return "rescued-vv"
    return "published"


def _unique_patients(rows):
    patients = set()
    for r in rows:
        for a in str(r.get("accession_list", "")).split(";"):
            a = a.strip()
            if a:
                patients.add(a)
    return patients


def _recompute_dedup_count():
    dedup_script = os.path.join(SCRIPTS_DIR, "Step7_Merging", "dedup_merged_variants.py")
    with tempfile.TemporaryDirectory() as tmp:
        tmp_out   = os.path.join(tmp, "recomputed.tsv")
        tmp_stats = os.path.join(tmp, "recomputed_stats.txt")
        env = {
            **os.environ,
            "INPUT_TSV":  MERGED_TSV,
            "OUTPUT_TSV": tmp_out,
            "STATS_TXT":  tmp_stats,
        }
        subprocess.run(
            [sys.executable, dedup_script],
            env=env, check=True, capture_output=True,
        )
        rows = _read_tsv(tmp_out)
        distinct = len(rows)
        pairs    = sum(int(r.get("patient_count", 0) or 0) for r in rows)
        return distinct, pairs


def main():
    rescued_rows = _read_tsv(RESCUED_DEDUP_TSV)
    primary_rows = _read_tsv(PRIMARY_DEDUP_TSV)
    unresolved   = _read_tsv(UNRESOLVED_DISP_TSV)

    for row in rescued_rows:
        row["provenance"] = _provenance(row.get("source_track", ""))

    null_prov = [r for r in rescued_rows if not r.get("provenance")]
    if null_prov:
        sys.exit(f"ACCEPTANCE FAIL: {len(null_prov)} rows have null provenance")

    if rescued_rows:
        existing     = list(rescued_rows[0].keys())
        final_fields = ["provenance"] + [f for f in existing if f != "provenance"]
    else:
        final_fields = ["provenance"]

    with open(OUT_FLAT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=final_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rescued_rows)

    n_unresolved = len(unresolved)
    pub_distinct = len(primary_rows)
    pub_pairs    = sum(int(r.get("patient_count", 0) or 0) for r in primary_rows)
    pub_patients = len(_unique_patients(primary_rows))
    pub_denom    = pub_distinct + n_unresolved
    pub_rate     = pub_distinct / pub_denom if pub_denom else 0.0

    fin_distinct = len(rescued_rows)
    fin_pairs    = sum(int(r.get("patient_count", 0) or 0) for r in rescued_rows)
    fin_patients = len(_unique_patients(rescued_rows))
    fin_rate     = fin_distinct / pub_denom if pub_denom else 0.0

    recomputed, recomputed_pairs = _recompute_dedup_count()
    if recomputed != pub_distinct:
        sys.exit(
            f"ACCEPTANCE FAIL: recomputed published distinct {recomputed} "
            f"!= primary dedup row count {pub_distinct}"
        )
    if recomputed_pairs != pub_pairs:
        sys.exit(
            f"ACCEPTANCE FAIL: recomputed published pairs {recomputed_pairs} "
            f"!= primary dedup pair count {pub_pairs}"
        )

    if fin_distinct < pub_distinct:
        sys.exit(
            f"ACCEPTANCE FAIL: final distinct {fin_distinct} < published {pub_distinct}"
        )

    git_result = subprocess.run(
        ["git", "status", "--porcelain", "--", PRIMARY_DEDUP_TSV],
        capture_output=True, text=True, cwd=SCRIPTS_DIR,
    )
    if git_result.stdout.strip():
        sys.exit(
            f"ACCEPTANCE FAIL: {PRIMARY_DEDUP_TSV} is not git-clean: "
            f"{git_result.stdout.strip()!r}"
        )

    metrics = [
        ["metric",                "published",         "final"],
        ["distinct_variants",     pub_distinct,        fin_distinct],
        ["patient_variant_pairs", pub_pairs,           fin_pairs],
        ["unique_patients",       pub_patients,        fin_patients],
        ["resolution_rate",       f"{pub_rate:.4f}",   f"{fin_rate:.4f}"],
    ]
    with open(OUT_METRICS, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(metrics)

    print(f"lovd_flat_final          : {OUT_FLAT}  ({fin_distinct} rows)")
    print(f"rescue_metrics_comparison: {OUT_METRICS}")
    print(f"published rate           : {pub_rate:.1%}  ({pub_distinct}/{pub_denom})")
    print(f"final rate               : {fin_rate:.1%}  ({fin_distinct}/{pub_denom})")
    print("All acceptance checks passed.")


if __name__ == "__main__":
    main()
