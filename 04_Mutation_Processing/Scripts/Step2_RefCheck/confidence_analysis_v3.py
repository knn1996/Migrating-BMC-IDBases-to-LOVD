import os
import argparse
import pandas as pd
from scipy.stats import binomtest

THESIS_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
PROC_DIR   = os.path.join(THESIS_DIR, "04_Mutation_Processing")

ALL_MUT_TSV = os.path.join(PROC_DIR, "Output", "Step1_Extraction", "all_mutations.tsv")
OUT_DIR     = os.path.join(PROC_DIR, "Output", "Step2_RefCheck")

OFFSET_FILES_BASE = {
    "NG_IDRefseq": "lrg_offset_results.csv",
    "NG_MANE":     "MANE_offset_results.csv",
}

NULL_P     = 0.25
SUB_ONLY   = {"substitution"}
ALL_ANCHOR = {"substitution", "deletion", "indel", "duplication"}


def analyze(all_mut, unique_vars, offset_df, anchor_classes):
    rows = []
    for _, off in offset_df.iterrows():
        gene   = off["gene"]
        status = off["status"]
        strand = off.get("strand", "")
        acc    = off.get("ng_accession", "")

        gene_uniq = unique_vars[unique_vars["gene"] == gene]
        mc_counts = gene_uniq["mut_class"].value_counts().to_dict()

        n_substitution = mc_counts.get("substitution", 0)
        n_deletion     = mc_counts.get("deletion",     0)
        n_duplication  = mc_counts.get("duplication",  0)
        n_indel        = mc_counts.get("indel",        0)
        n_insertion    = mc_counts.get("insertion",    0)
        n_unknown      = mc_counts.get("unknown",      0)
        n_total        = len(gene_uniq)

        n_anchorable    = n_substitution + n_deletion + n_duplication + n_indel
        n_unanchorable  = n_insertion
        frac_anchorable = n_anchorable / n_total if n_total else 0.0

        n_matched  = int(off["matched"]) if pd.notna(off["matched"]) else 0
        n_in_test  = int(off["total"])   if pd.notna(off["total"])   else 0
        match_rate = n_matched / n_in_test if n_in_test else 0.0

        binomial_p  = ""
        wilson_low  = ""
        wilson_high = ""
        if n_in_test > 0:
            bt          = binomtest(k=n_matched, n=n_in_test, p=NULL_P, alternative="greater")
            binomial_p  = bt.pvalue
            ci          = bt.proportion_ci(confidence_level=0.95, method="wilson")
            wilson_low  = round(ci.low,  4)
            wilson_high = round(ci.high, 4)

        gene_span_bp   = 0
        anchor_span_bp = 0
        coverage_frac  = 0.0
        max_gap_bp     = 0
        anchors_per_kb = 0.0

        if pd.notna(off["sstart"]) and pd.notna(off["send"]):
            gene_span_bp = int(off["send"]) - int(off["sstart"]) + 1

            anchor_rows = all_mut[
                (all_mut["gene"] == gene) &
                (all_mut["variant_type"] == "genomic") &
                (all_mut["mut_class"].isin(anchor_classes))
            ].dropna(subset=["pos_start"])

            anchor_rows = anchor_rows.drop_duplicates(subset=["accession", "allele_num"])
            positions   = sorted(anchor_rows["pos_start"].astype(int).tolist())

            n_positions = len(positions)
            if n_positions > 0 and gene_span_bp > 0:
                anchor_span_bp = positions[-1] - positions[0]
                if anchor_span_bp > gene_span_bp:
                    anchor_span_bp = gene_span_bp
                coverage_frac  = anchor_span_bp / gene_span_bp
                anchors_per_kb = n_positions / (gene_span_bp / 1000)
                if n_positions >= 2:
                    gaps       = [positions[i+1] - positions[i] for i in range(n_positions - 1)]
                    max_gap_bp = max(gaps)

        rows.append({
            "gene":             gene,
            "accession":        acc,
            "status":           status,
            "strand":           strand,
            "n_total":          n_total,
            "n_substitution":   n_substitution,
            "n_deletion":       n_deletion,
            "n_duplication":    n_duplication,
            "n_indel":          n_indel,
            "n_insertion":      n_insertion,
            "n_unknown":        n_unknown,
            "n_anchorable":     n_anchorable,
            "n_unanchorable":   n_unanchorable,
            "frac_anchorable":  round(frac_anchorable, 4),
            "n_matched":        n_matched,
            "n_in_test":        n_in_test,
            "match_rate":       round(match_rate, 4),
            "binomial_p":       binomial_p,
            "wilson_ci_low":    wilson_low,
            "wilson_ci_high":   wilson_high,
            "gene_span_bp":     gene_span_bp,
            "anchor_span_bp":   anchor_span_bp,
            "coverage_frac":    round(coverage_frac, 4),
            "max_gap_bp":       max_gap_bp,
            "anchors_per_kb":   round(anchors_per_kb, 3),
        })

    return pd.DataFrame(rows).sort_values("gene")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--include-indels", action="store_true",
                        help="Read _with_indels offset files; spatial coverage uses all anchorable types")
    args = parser.parse_args()

    suffix         = "_with_indels" if args.include_indels else ""
    anchor_classes = ALL_ANCHOR if args.include_indels else SUB_ONLY

    print(f"Anchor types for spatial: {sorted(anchor_classes)}")
    print(f"Offset file suffix:       '{suffix}'")
    print()

    all_mut = pd.read_csv(ALL_MUT_TSV, sep="\t", dtype=str).fillna("")
    all_mut["gene"]         = all_mut["gene"].str.strip().str.upper()
    all_mut["mut_class"]    = all_mut["mut_class"].str.strip().str.lower()
    all_mut["variant_type"] = all_mut["variant_type"].str.strip().str.lower()
    all_mut["pos_start"]    = pd.to_numeric(all_mut["pos_start"], errors="coerce")

    unique_keys = ["gene", "accession", "allele_num", "mut_class"]
    unique_vars = all_mut[unique_keys].drop_duplicates()

    for track_name, base_filename in OFFSET_FILES_BASE.items():
        stem, ext       = os.path.splitext(base_filename)
        offset_filename = f"{stem}{suffix}{ext}"
        offset_path     = os.path.join(OUT_DIR, offset_filename)

        if not os.path.exists(offset_path):
            print(f"[{track_name}] SKIP - file not found: {offset_path}\n")
            continue

        offset = pd.read_csv(offset_path, dtype=str).fillna("")
        offset["gene"]    = offset["gene"].str.strip().str.upper()
        offset["matched"] = pd.to_numeric(offset["matched"], errors="coerce")
        offset["total"]   = pd.to_numeric(offset["total"],   errors="coerce")
        offset["sstart"]  = pd.to_numeric(offset["sstart"],  errors="coerce")
        offset["send"]    = pd.to_numeric(offset["send"],    errors="coerce")

        out_df   = analyze(all_mut, unique_vars, offset, anchor_classes)
        out_path = os.path.join(OUT_DIR, f"gene_confidence_{track_name}{suffix}.csv")
        out_df.to_csv(out_path, index=False)
        print(f"[{track_name}] Wrote {len(out_df)} rows to {out_path}")

        ok   = out_df[out_df["status"].isin(["ok", "below_threshold", "ok_reverse"])]
        ok_p = ok[ok["binomial_p"] != ""].copy()
        ok_p["binomial_p"] = pd.to_numeric(ok_p["binomial_p"])
        print(f"  status counts: {out_df['status'].value_counts().to_dict()}")
        print(f"  match_rate mean: {ok['match_rate'].mean():.4f}")
        print(f"  median n_in_test: {int(ok['n_in_test'].median())}")
        print(f"  p<0.05: {(ok_p['binomial_p']<0.05).sum()}/{len(ok_p)}  "
              f"p<0.001: {(ok_p['binomial_p']<0.001).sum()}/{len(ok_p)}  "
              f"underpowered (p>=0.05): {(ok_p['binomial_p']>=0.05).sum()}/{len(ok_p)}")
        print(f"  coverage_frac>1: {(ok['coverage_frac']>1).sum()} (should be 0)")
        print()


if __name__ == "__main__":
    main()
