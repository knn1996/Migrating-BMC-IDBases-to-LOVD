import os
import re
import pandas as pd
import requests
from pathlib import Path
import time

IDBASE_DIR   = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\idbase"
REF_SUMMARY  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step2_RefCheck\reference_summary.csv"
OUT_CSV      = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\IDBases_Summary.csv"

SLEEP = 0.3

genes = sorted(
    re.sub(r"base$", "", d, flags=re.IGNORECASE)
    for d in os.listdir(IDBASE_DIR)
    if os.path.isdir(os.path.join(IDBASE_DIR, d))
    and d.lower().endswith("base")
    and d.lower() != "immunomebase"
)

ref_df = pd.read_csv(REF_SUMMARY)
ref_df["match_pct"] = pd.to_numeric(ref_df["match_pct"], errors="coerce")

def best_ref(gene):
    rows = ref_df[ref_df["gene"].str.upper() == gene.upper()]
    if rows.empty:
        return "", ""
    best = rows.loc[rows["match_pct"].idxmax()]
    return best["ref"], best["match_pct"]

def get_uniprot(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_id:9606+AND+reviewed:true&fields=accession&format=json&size=1"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        results = r.json().get("results", [])
        if results:
            acc = results[0]["primaryAccession"]
            return acc, f"https://www.uniprot.org/uniprot/{acc}"
    except Exception:
        pass
    return "", ""

rows = []
for gene in genes:
    ref, pct = best_ref(gene)
    uniprot_id, uniprot_url = get_uniprot(gene)
    print(f"  {gene:<15} UniProt={uniprot_id}")
    rows.append({
        "gene_name":          gene,
        "reference_sequence": ref,
        "mutation_match_pct": pct,
        "uniprot_id":         uniprot_id,
        "uniprot_link":       uniprot_url,
        "hg_version":         "",
        "bed_generated":      "",
        "liftover_to_hg38":   "",
        "converted_to_vcf":   "",
        "mutalyzered":        "",
    })
    time.sleep(SLEEP)

pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f"\nWritten {len(rows)} genes to {OUT_CSV}")
