"""
update_liftover_status.py
=========================
Updates the liftover_to_hg38 column in IDBases_Summary.csv.

Sets "yes" if a BED file exists for the gene in BED_HG38_DIR,
otherwise "no".
"""

import os
import csv

THESIS_DIR   = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
PIPELINE_DIR = os.path.join(THESIS_DIR, "04_Mutation_Processing")

SUMMARY_CSV  = os.path.join(PIPELINE_DIR, "Output", "IDBases_Summary.csv")
BED_HG38_DIR = os.path.join(THESIS_DIR, "03_BED_Files", "BED_IDRefseq", "hg38")


def get_genes_with_bed(bed_dir):
    if not os.path.isdir(bed_dir):
        return set()
    return {
        fname.split("_")[0]
        for fname in os.listdir(bed_dir)
        if fname.lower().endswith(".bed")
    }


def main():
    genes_with_bed = get_genes_with_bed(BED_HG38_DIR)

    with open(SUMMARY_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        rows = list(reader)

    yes_count = no_count = 0
    for row in rows:
        gene = row["gene_name"].strip()
        if gene in genes_with_bed:
            row["liftover_to_hg38"] = "yes"
            yes_count += 1
        else:
            row["liftover_to_hg38"] = "no"
            no_count += 1

    with open(SUMMARY_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"liftover_to_hg38 = yes: {yes_count}")
    print(f"liftover_to_hg38 = no : {no_count}")


if __name__ == "__main__":
    main()
