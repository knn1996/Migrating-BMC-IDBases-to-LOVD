"""
update_hg_version.py
====================
Updates the hg_version column in IDBases_Summary.csv based on
blast_first_hits_source_IDRefseq.tsv.

Selection rule: highest pct_identity; ties broken by highest build
(hg18 > hg17 > hg16).
"""

import os
import csv

THESIS_DIR   = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
PIPELINE_DIR = os.path.join(THESIS_DIR, "04_Mutation_Processing")

SUMMARY_CSV  = os.path.join(PIPELINE_DIR, "Output", "IDBases_Summary.csv")
BLAST_TSV    = os.path.join(PIPELINE_DIR, "Output", "Step3_BLAST", "blast_first_hits_source_IDRefseq.tsv")

BUILD_RANK = {"hg16": 0, "hg17": 1, "hg18": 2}


def load_best_build(blast_path):
    best = {}
    with open(blast_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene  = row["query_id"].strip()
            build = row["build"].strip().lower()
            pct   = float(row["pct_identity"])
            rank  = BUILD_RANK.get(build, -1)

            if gene not in best:
                best[gene] = (pct, rank, build)
            else:
                cur_pct, cur_rank, _ = best[gene]
                if (pct, rank) > (cur_pct, cur_rank):
                    best[gene] = (pct, rank, build)

    return {gene: build for gene, (_, _, build) in best.items()}


def main():
    best_build = load_best_build(BLAST_TSV)

    with open(SUMMARY_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)

    updated = 0
    for row in rows:
        gene = row["gene_name"].strip()
        build = best_build.get(gene, "")
        if build:
            row["hg_version"] = build
            updated += 1

    with open(SUMMARY_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Updated hg_version for {updated} genes.")
    print(f"No BLAST hit for {len(rows) - updated} genes.")


if __name__ == "__main__":
    main()
