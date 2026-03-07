"""
add_nm_to_lrg.py
================
Adds a NM column to LRG.csv by looking up each gene NG_ accession
in LRG_RefSeqGene.txt.

Priority for NM_ selection:
    1. Category == "reference standard"  (preferred transcript)
    2. Any other NM_ aligned to the same NG_ (fallback)

Output: LRG_with_NM.csv  (same folder as LRG.csv)
"""

import pandas as pd

LRG_CSV     = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\LRG.csv"
REFSEQ_TXT  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\LRG_RefSeqGene.txt"
OUTPUT_CSV  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\LRG_with_NM.csv"

# Load files
lrg = pd.read_csv(LRG_CSV)
ref = pd.read_csv(REFSEQ_TXT, sep="\t", comment="#",
                  names=["tax_id","GeneID","Symbol","RSG","LRG","RNA",
                         "t","Protein","p","Category"])

# Keep only rows with an NM_ transcript
ref = ref[ref["RNA"].str.startswith("NM_", na=False)].copy()

# Build NG_ -> NM_ lookup
# Priority 1: reference standard
ref_standard = (
    ref[ref["Category"].str.strip() == "reference standard"]
    .drop_duplicates(subset="RSG")
    .set_index("RSG")["RNA"]
)

# Priority 2: any NM_ aligned to that NG_ (first one found)
ref_any = (
    ref.drop_duplicates(subset="RSG")
    .set_index("RSG")["RNA"]
)

def pick_nm(ng):
    if pd.isna(ng) or ng == "":
        return ""
    if ng in ref_standard:
        return ref_standard[ng]
    if ng in ref_any:
        return ref_any[ng]
    return ""

# Apply and save
lrg["NM"] = lrg["RSG"].apply(pick_nm)
lrg.to_csv(OUTPUT_CSV, index=False)

# Report
total  = len(lrg)
has_ng = lrg["RSG"].notna() & (lrg["RSG"] != "")
has_nm = lrg["NM"] != ""
no_ng  = lrg[~has_ng]["name"].tolist()
no_nm  = lrg[has_ng & ~has_nm]["name"].tolist()

print(f"Total genes                     : {total}")
print(f"Genes with NG_                  : {has_ng.sum()}")
print(f"Genes with NM_                  : {has_nm.sum()}")
print(f"Genes without NG_               : {len(no_ng)}  -> {no_ng}")
print(f"Genes with NG_ but no NM_ found : {len(no_nm)}  -> {no_nm}")
print(f"\nOutput saved to: {OUTPUT_CSV}")
print("\nSample output:")
print(lrg.head(10).to_string(index=False))
