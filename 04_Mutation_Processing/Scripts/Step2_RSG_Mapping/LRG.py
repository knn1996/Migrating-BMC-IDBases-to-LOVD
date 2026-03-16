import pandas as pd
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
THESIS_DIR  = _SCRIPT_DIR.parent.parent.parent

excel_file  = THESIS_DIR / "02_Source_Database" / "IDBases_Summary.xlsx"
lrg_file    = THESIS_DIR / "04_Mutation_Processing" / "Output" / "LRG_RefSeqGene.txt"
output_file = THESIS_DIR / "04_Mutation_Processing" / "Output" / "LRG.csv"

# -----------------------------
# Read Excel file
# -----------------------------
df_excel = pd.read_excel(excel_file)

# Get gene names (change column name here if needed)
gene_names = df_excel['name'].dropna().unique()

# -----------------------------
# Read LRG reference file
# -----------------------------
df_lrg = pd.read_csv(
    lrg_file,
    sep="\t",
    comment="#",
    header=None,
    names=[
        "tax_id","GeneID","Symbol","RSG","LRG",
        "RNA","t","Protein","p","Category"
    ]
)

# -----------------------------
# Match gene name -> RSG
# -----------------------------
results = []

for gene in gene_names:
    match = df_lrg[df_lrg["Symbol"] == gene]

    if not match.empty:
        rsg = match.iloc[0]["RSG"]
    else:
        rsg = None

    results.append({"name": gene, "RSG": rsg})

# -----------------------------
# Save result
# -----------------------------
df_out = pd.DataFrame(results)
df_out.to_csv(output_file, index=False)

print(f"File saved to: {output_file}")