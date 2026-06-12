import os
import re
import pandas as pd

THESIS_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
STEP8_DIR  = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output", "Step8_Merging")
IN_PATH    = os.path.join(STEP8_DIR, "dedup_merged_variants.tsv")
OUT_PATH   = os.path.join(STEP8_DIR, "lovd_flat.tsv")

NM_IN_PARENS = re.compile(r"\((N[MPR]_[^)]+)\)")
NM_BARE      = re.compile(r"^(N[MPR]_[^:()]+)")

OUTPUT_COLUMNS = [
    "Individual/Gender",
    "Individual/Origin/Geographic",
    "Individual/Age",
    "Individual/Remarks",
    "Phenotype/Disease",
    "Phenotype/Additional",
    "Screening/Technique",
    "chromosome",
    "VariantOnGenome/DNA",
    "VariantOnGenome/Remarks",
    "Reference/PubMed_ID",
    "Reference/Authors",
    "Reference/Title",
    "Reference/Location",
    "transcriptid",
    "VariantOnTranscript/DNA",
    "VariantOnTranscript/RNA",
    "VariantOnTranscript/Protein",
    "_gene",
    "_dedup_key",
    "_patient_count",
    "_accession_list",
    "_mane_select",
    "_source_track",
    "_status",
    "_variant_class",
]

PATIENT_PLACEHOLDERS = [
    "Individual/Gender",
    "Individual/Origin/Geographic",
    "Individual/Age",
    "Individual/Remarks",
    "Phenotype/Disease",
    "Phenotype/Additional",
    "Screening/Technique",
    "VariantOnGenome/Remarks",
    "Reference/PubMed_ID",
    "Reference/Authors",
    "Reference/Title",
    "Reference/Location",
]


def split_body(s):
    if not isinstance(s, str) or ":" not in s:
        return ""
    return s.split(":", 1)[1].strip()


def extract_transcript_any(s):
    if not isinstance(s, str) or not s:
        return ""
    m = NM_IN_PARENS.search(s)
    if m:
        return m.group(1)
    m = NM_BARE.search(s)
    return m.group(1) if m else ""


def clean_chr(s):
    if not isinstance(s, str):
        return ""
    return s.replace("chr", "").strip()


RE_REPEAT = re.compile(r"[ACGT]\[\d+\]")

def classify_variant(*notations):
    for n in notations:
        if not isinstance(n, str) or not n:
            continue
        if n.endswith(":c.=") or n.endswith(":g.=") or n in ("c.=", "g.="):
            return "no_change"
        if RE_REPEAT.search(n):
            return "repeat"
        if "delins" in n: return "delins"
        if "dup"    in n: return "duplication"
        if "ins"    in n: return "insertion"
        if "del"    in n: return "deletion"
        if "inv"    in n: return "inversion"
        if ">"      in n: return "substitution"
    return "other"


def main():
    df = pd.read_csv(IN_PATH, sep="\t", dtype=str).fillna("")

    is_genomic = df["normalized"].str.startswith("NC_")

    out = pd.DataFrame(index=df.index, columns=OUTPUT_COLUMNS)
    out[:] = ""

    out["chromosome"]            = df["chrom"].apply(clean_chr)
    out["VariantOnGenome/DNA"]   = df["normalized"].where(is_genomic, "")

    c_source = df["c_hgvs"].where(df["c_hgvs"] != "", df["normalized"].where(~is_genomic, ""))

    out["transcriptid"]                  = c_source.apply(extract_transcript_any)
    out["VariantOnTranscript/DNA"]       = c_source.apply(split_body)
    out["VariantOnTranscript/RNA"]       = df["r_hgvs"].apply(split_body)
    out["VariantOnTranscript/Protein"]   = df["p_hgvs"].apply(split_body)

    out["_gene"]            = df["gene"]
    out["_dedup_key"]       = df["dedup_key"]
    out["_patient_count"]   = df["patient_count"]
    out["_accession_list"]  = df["accession_list"]
    out["_mane_select"]     = df.get("mane_select", "")
    out["_source_track"]    = df.get("source_track", "")
    out["_status"]          = df.get("status", "")

    out["_variant_class"] = [
        classify_variant(c, g)
        for c, g in zip(out["VariantOnTranscript/DNA"], out["VariantOnGenome/DNA"])
    ]

    for col in PATIENT_PLACEHOLDERS:
        if col not in out.columns:
            out[col] = ""

    out = out[OUTPUT_COLUMNS]
    out.to_csv(OUT_PATH, sep="\t", index=False)

    n_genomic = is_genomic.sum()
    n_tx_only = (~is_genomic).sum()

    print(f"Wrote {len(out)} rows -> {OUT_PATH}")
    print(f"  genomic-anchored rows (NC_)   : {n_genomic}")
    print(f"  transcript-only rows (NM_)    : {n_tx_only}")
    print(f"  rows with chromosome filled   : {(out['chromosome'] != '').sum()}")
    print(f"  rows with VariantOnGenome     : {(out['VariantOnGenome/DNA'] != '').sum()}")
    print(f"  rows with transcriptid        : {(out['transcriptid'] != '').sum()}")
    print(f"  rows with c. body             : {(out['VariantOnTranscript/DNA'] != '').sum()}")
    print(f"  rows with p. body             : {(out['VariantOnTranscript/Protein'] != '').sum()}")
    print(f"  total patient links (sum)     : {pd.to_numeric(out['_patient_count'], errors='coerce').sum():.0f}")
    print()
    print(f"=== Variant class distribution ===")
    print(out['_variant_class'].value_counts().to_string())


if __name__ == "__main__":
    main()
