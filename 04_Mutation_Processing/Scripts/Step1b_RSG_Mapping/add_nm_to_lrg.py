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

import os
import pandas as pd
from pathlib import Path

_SCRIPT_DIR   = Path(__file__).parent
_MUTATION_DIR = _SCRIPT_DIR.parent.parent

LRG_CSV     = str(_MUTATION_DIR / "Output" / "Step1b_RSG_Mapping" / "LRG.csv")
REFSEQ_TXT  = str(_MUTATION_DIR / "Output" / "Step1b_RSG_Mapping" / "LRG_RefSeqGene.txt")
OUTPUT_CSV  = str(_MUTATION_DIR / "Output" / "Step1b_RSG_Mapping" / "LRG_with_NM.csv")

lrg = pd.read_csv(LRG_CSV)
ref = pd.read_csv(REFSEQ_TXT, sep="\t", comment="#",
                  names=["tax_id","GeneID","Symbol","RSG","LRG","RNA",
                         "t","Protein","p","Category"])

ref = ref[ref["RNA"].str.startswith("NM_", na=False)].copy()

ref_standard = (
    ref[ref["Category"].str.strip() == "reference standard"]
    .drop_duplicates(subset="RSG")
    .set_index("RSG")["RNA"]
)

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

lrg["NM"] = lrg["RSG"].apply(pick_nm)
lrg.to_csv(OUTPUT_CSV, index=False)

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
