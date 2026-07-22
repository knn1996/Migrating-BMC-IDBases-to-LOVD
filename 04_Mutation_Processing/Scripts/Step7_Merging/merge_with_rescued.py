"""
merge_with_rescued.py
=====================
Concatenates merged_variants.tsv with resolved variants from
resolve_unresolved.py, expanding each de-duplicated resolved row back
into one row per patient accession so that dedup_merged_variants.py
correctly accumulates patient counts.

Env vars
--------
MERGED_TSV    merged_variants.tsv  (23-col; one row per patient-variant)
RESOLVED_TSV  resolved_unresolved.tsv (23+4 col; one row per distinct variant)
OUT_TSV       augmented_merged_variants.tsv (23-col; expanded, concatenated)
STATS_TXT     plain-text before/after row counts
"""

import csv
import os
from pathlib import Path

MERGED_TSV   = os.environ["MERGED_TSV"]
RESOLVED_TSV = os.environ["RESOLVED_TSV"]
OUT_TSV      = os.environ["OUT_TSV"]
STATS_TXT    = os.environ["STATS_TXT"]

for _p in (OUT_TSV, STATS_TXT):
    Path(_p).parent.mkdir(parents=True, exist_ok=True)

OUT_FIELDS = [
    "gene", "accession", "allele_num", "sysname", "mut_type", "ref", "alt",
    "chrom", "pos_hg38", "strand", "status", "source_track",
    "nc_hgvs", "c_hgvs", "g_hgvs", "r_hgvs", "p_hgvs",
    "protein_pos_first", "mane_select", "mutalyzer_gene",
    "hgvs_input", "normalized", "errors",
]


def main():
    with open(MERGED_TSV, encoding="utf-8") as f:
        merged_rows = list(csv.DictReader(f, delimiter="\t"))

    with open(RESOLVED_TSV, encoding="utf-8") as f:
        resolved_rows = list(csv.DictReader(f, delimiter="\t"))

    n_merged = len(merged_rows)

    rescued = [r for r in resolved_rows if r.get("resolving_tool", "").strip()]

    expanded = []
    for row in rescued:
        acc_list_raw = row.get("accession_list", "").strip()
        accessions = [a.strip() for a in acc_list_raw.split(";") if a.strip()]
        if not accessions:
            accessions = [row.get("accession", "")]
        for acc in accessions:
            new_row = {c: row.get(c, "") for c in OUT_FIELDS}
            new_row["accession"] = acc
            expanded.append(new_row)

    n_rescued_variants = len(rescued)
    n_expanded = len(expanded)

    out_rows = [{c: r.get(c, "") for c in OUT_FIELDS} for r in merged_rows] + expanded

    with open(OUT_TSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    stats = (
        f"merge_with_rescued summary\n"
        f"==========================\n"
        f"  merged_variants rows (pre-rescue)    : {n_merged}\n"
        f"  resolved distinct variants           : {n_rescued_variants}\n"
        f"  expanded patient-variant rows added  : {n_expanded}\n"
        f"  total rows in augmented output       : {len(out_rows)}\n"
        f"  output                               : {OUT_TSV}\n"
    )
    Path(STATS_TXT).write_text(stats, encoding="utf-8")
    print(stats, end="")


if __name__ == "__main__":
    main()
