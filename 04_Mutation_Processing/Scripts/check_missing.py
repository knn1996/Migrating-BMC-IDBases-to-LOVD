import csv
from pathlib import Path

IDBASE_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\idbase_pub_html_only"
SUMMARY_CSV = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\IDBases_Summary.csv"

with open(SUMMARY_CSV, encoding="utf-8") as f:
    summary_genes = {row["gene_name"].upper() for row in csv.DictReader(f) if row["gene_name"]}

def gene_from_filename(path):
    stem = Path(path).stem.lower()
    for suffix in ("_pub", "pub", "_base", "base"):
        while stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    return stem.upper()

folder_genes = {gene_from_filename(fp): fp.name
                for fp in Path(IDBASE_DIR).glob("*pub.html")}

in_folder_not_summary = sorted(set(folder_genes) - summary_genes)
in_summary_not_folder = sorted(summary_genes - set(folder_genes))

print(f"Files in folder: {len(folder_genes)}")
print(f"Genes in summary: {len(summary_genes)}")
print()
print(f"In folder but NOT in summary ({len(in_folder_not_summary)}):")
for g in in_folder_not_summary:
    print(f"  - {g}  (file: {folder_genes[g]})")
print()
print(f"In summary but NOT in folder ({len(in_summary_not_folder)}):")
for g in in_summary_not_folder:
    print(f"  - {g}")