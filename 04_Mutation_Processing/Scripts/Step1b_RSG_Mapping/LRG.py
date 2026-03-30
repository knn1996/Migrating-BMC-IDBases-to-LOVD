import pandas as pd
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
THESIS_DIR  = _SCRIPT_DIR.parent.parent.parent

summary_file = THESIS_DIR / "04_Mutation_Processing" / "Output" / "IDBases_Summary.csv"
lrg_file     = THESIS_DIR / "04_Mutation_Processing" / "Output" / "Step1b_RSG_Mapping" / "LRG_RefSeqGene.txt"
alias_file   = THESIS_DIR / "04_Mutation_Processing" / "Output" / "alias.csv"
output_file  = THESIS_DIR / "04_Mutation_Processing" / "Output" / "Step1b_RSG_Mapping" / "LRG.csv"

df_summary = pd.read_csv(summary_file)
gene_names = df_summary['gene_name'].unique()

df_lrg = pd.read_csv(
    lrg_file,
    sep="\t",
    comment="#",
    header=None,
    names=["tax_id","GeneID","Symbol","RSG","LRG","RNA","t","Protein","p","Category"]
)

df_alias = pd.read_csv(alias_file)
alias_map = dict(zip(df_alias["gene"], df_alias["alias"]))

results = []
for gene in gene_names:
    match = df_lrg[df_lrg["Symbol"] == gene]
    if match.empty and gene in alias_map:
        match = df_lrg[df_lrg["Symbol"] == alias_map[gene]]
    rsg = match.iloc[0]["RSG"] if not match.empty else None
    results.append({"name": gene, "RSG": rsg})

df_out = pd.DataFrame(results)
df_out.to_csv(output_file, index=False)
print(f"File saved to: {output_file}")
