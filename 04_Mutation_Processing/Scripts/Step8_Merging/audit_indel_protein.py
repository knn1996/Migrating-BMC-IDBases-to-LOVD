import os
import re
import pandas as pd

THESIS_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
STEP8_DIR  = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output", "Step8_Merging")
LOG_DIR    = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Logs")
IN_PATH    = os.path.join(STEP8_DIR, "lovd_flat.tsv")
OUT_PATH   = os.path.join(STEP8_DIR, "indel_audit.tsv")
LOG_PATH   = os.path.join(LOG_DIR,   "indel_audit_summary.txt")

RE_RANGE       = re.compile(r"c\.(-?\d+)(?:[+\-]\d+)?_(-?\d+)(?:[+\-]\d+)?")
RE_SINGLE_POS  = re.compile(r"c\.(-?\d+)(?:[+\-]\d+)?(?:del|dup|ins)")
RE_INS_BASES   = re.compile(r"ins([ACGT]+)")
RE_DELINS_BASES = re.compile(r"delins([ACGT]+)")
RE_DEL_BASES   = re.compile(r"del([ACGT]+)")
RE_DUP_BASES   = re.compile(r"dup([ACGT]+)")


def length_from_range(c_dna):
    m = RE_RANGE.search(c_dna)
    if m:
        return abs(int(m.group(2)) - int(m.group(1))) + 1
    if RE_SINGLE_POS.search(c_dna):
        return 1
    return None


def length_explicit_bases(c_dna):
    m = RE_DELINS_BASES.search(c_dna)
    if m: return None
    m = RE_DEL_BASES.search(c_dna)
    if m: return len(m.group(1))
    m = RE_DUP_BASES.search(c_dna)
    if m: return len(m.group(1))
    m = RE_INS_BASES.search(c_dna)
    if m: return len(m.group(1))
    return None


def compute_length_change(c_dna, vclass):
    if not c_dna:
        return None
    if vclass == "substitution":
        return 0
    if vclass in ("deletion", "duplication", "insertion"):
        n_explicit = length_explicit_bases(c_dna)
        if n_explicit is not None:
            return n_explicit
        return length_from_range(c_dna)
    if vclass == "delins":
        n_del = length_from_range(c_dna)
        m = RE_DELINS_BASES.search(c_dna)
        n_ins = len(m.group(1)) if m else None
        if n_del is None or n_ins is None:
            return None
        return abs(n_ins - n_del)
    return None


def classify_protein(p_body):
    if not p_body:
        return "empty"
    p = p_body.strip()
    if p in ("p.(=)", "p.=", "p.?"):
        return "wildtype"
    if "fs" in p or "Ter" in p or "*" in p:
        return "frameshift_or_stop"
    if "del" in p or "dup" in p or "ins" in p:
        return "inframe_indel"
    if "ext" in p:
        return "extension"
    if re.search(r"p\.\(?[A-Z][a-z][a-z]\d+[A-Z][a-z][a-z]", p):
        return "missense"
    return "other"


def audit_row(row):
    vclass = row["_variant_class"]
    c_dna  = row["VariantOnTranscript/DNA"]
    p_body = row["VariantOnTranscript/Protein"]
    p_state = classify_protein(p_body)
    length_change = compute_length_change(c_dna, vclass)

    flags = []

    if vclass == "no_change":
        flags.append("no_change_exclude_from_lovd")

    elif vclass == "substitution":
        if p_state == "empty":
            flags.append("substitution_protein_empty")
        elif p_state == "frameshift_or_stop" and "Ter" not in p_body and "*" not in p_body:
            flags.append("substitution_protein_says_frameshift")

    elif vclass in ("deletion", "duplication", "insertion"):
        if length_change is None:
            flags.append(f"{vclass}_length_unparseable")
        else:
            expected_inframe = (length_change % 3 == 0)
            if p_state == "empty":
                flags.append(f"{vclass}_protein_empty")
            elif p_state == "wildtype":
                flags.append(f"{vclass}_protein_wildtype")
            elif expected_inframe and p_state == "frameshift_or_stop":
                flags.append(f"{vclass}_expected_inframe_got_frameshift")
            elif not expected_inframe and p_state == "inframe_indel":
                flags.append(f"{vclass}_expected_frameshift_got_inframe")
            elif not expected_inframe and p_state == "missense":
                flags.append(f"{vclass}_expected_frameshift_got_missense")
            elif p_state == "missense":
                flags.append(f"{vclass}_protein_says_missense")

    elif vclass == "delins":
        if length_change is None:
            flags.append("delins_length_unparseable")
        else:
            expected_inframe = (length_change % 3 == 0)
            if p_state == "empty":
                flags.append("delins_protein_empty")
            elif p_state == "wildtype":
                flags.append("delins_protein_wildtype")
            elif expected_inframe and p_state == "frameshift_or_stop":
                flags.append("delins_expected_inframe_got_frameshift")
            elif not expected_inframe and p_state == "inframe_indel":
                flags.append("delins_expected_frameshift_got_inframe")

    elif vclass == "repeat":
        flags.append("repeat_manual_review_required")

    elif vclass == "inversion":
        if p_state == "empty":
            flags.append("inversion_protein_empty")

    return pd.Series({
        "length_change_bp": length_change,
        "protein_state":    p_state,
        "flag":              ";".join(flags) if flags else "",
    })


def main():
    df = pd.read_csv(IN_PATH, sep="\t", dtype=str).fillna("")
    print(f"Loaded {len(df)} variants from {IN_PATH}")

    audit = df.apply(audit_row, axis=1)
    df_out = pd.concat([df, audit], axis=1)

    flagged = df_out[df_out["flag"] != ""].copy()

    out_cols = [
        "_gene", "_dedup_key", "_variant_class",
        "VariantOnGenome/DNA", "VariantOnTranscript/DNA",
        "VariantOnTranscript/RNA", "VariantOnTranscript/Protein",
        "length_change_bp", "protein_state", "flag",
        "_source_track", "_mane_select", "_patient_count",
    ]
    flagged[out_cols].to_csv(OUT_PATH, sep="\t", index=False)

    print(f"\nWrote {len(flagged)} flagged variants -> {OUT_PATH}")
    print(f"  ({100 * len(flagged) / len(df):.1f}% of all variants)")

    print("\n=== Flag distribution ===")
    flag_counts = flagged["flag"].value_counts()
    for flag, n in flag_counts.items():
        print(f"  {n:>5}  {flag}")

    print("\n=== Per variant_class breakdown ===")
    class_breakdown = flagged.groupby("_variant_class").size().sort_values(ascending=False)
    for klass, n in class_breakdown.items():
        total_in_class = (df["_variant_class"] == klass).sum()
        print(f"  {klass:<15}  {n:>4}/{total_in_class:<4}  ({100*n/total_in_class:.1f}% flagged)")

    print("\n=== Top 10 genes by flagged-variant count ===")
    gene_breakdown = flagged.groupby("_gene").size().sort_values(ascending=False).head(10)
    for gene, n in gene_breakdown.items():
        print(f"  {gene:<15}  {n}")

    with open(LOG_PATH, "w", encoding="utf-8") as f:
        f.write(f"Indel/Protein Audit Summary\n")
        f.write(f"Total variants audited: {len(df)}\n")
        f.write(f"Flagged variants:       {len(flagged)} ({100*len(flagged)/len(df):.1f}%)\n\n")
        f.write("Flag distribution:\n")
        for flag, n in flag_counts.items():
            f.write(f"  {n:>5}  {flag}\n")
    print(f"\nSummary log -> {LOG_PATH}")


if __name__ == "__main__":
    main()
