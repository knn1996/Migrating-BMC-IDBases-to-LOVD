import os
import re
import argparse
import pandas as pd
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
_THESIS_DIR = (_SCRIPT_DIR / ".." / ".." / "..").resolve()
_PROC_DIR   = _THESIS_DIR / "04_Mutation_Processing"

SOURCE_DIRS = {
    "LRG":  _PROC_DIR / "DNA sequences" / "LRG FASTA file (Source 3)" / "NG",
    "MANE": _PROC_DIR / "DNA sequences" / "Mane_Select_NG",
}

OUTPUT_CSVS = {
    "LRG":  _PROC_DIR / "Output" / "Step2_RefCheck" / "lrg_offset_results.csv",
    "MANE": _PROC_DIR / "Output" / "Step2_RefCheck" / "MANE_offset_results.csv",
}

MUTATIONS_TSV      = _PROC_DIR / "Output" / "Step1_Extraction" / "all_mutations.tsv"
MAX_SEQ_LEN        = 400_000
MATCH_THRESHOLD    = 0.90
TARGET_MUT_CLASSES = {"substitution"}


def parse_ref_nucleotides(notation: str) -> str:
    if ">" in notation:
        m = re.search(r"([ACGTN]+)>[ACGTN]+$", notation, re.IGNORECASE)
        return m.group(1).upper() if m else ""
    return ""


def load_fasta(path):
    with open(path, encoding="utf-8") as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


def find_best_offset(ng_seq, mutations):
    seq_len   = len(ng_seq)
    positions = [(r["accession"], int(r["pos_start"]), int(r["pos_end"]), r["ref_nucleotides"].upper())
                 for r in mutations]
    total = len(positions)

    candidate_offsets = []
    for seed_acc, seed_start, seed_end, seed_ref in sorted(positions, key=lambda x: len(x[3]), reverse=True):
        candidate_offsets = []
        idx = 0
        while True:
            idx = ng_seq.find(seed_ref, idx)
            if idx == -1:
                break
            offset = idx - (seed_start - 1)
            if offset >= 0:
                candidate_offsets.append(offset)
            idx += 1
        if candidate_offsets:
            break

    if not candidate_offsets:
        return None, 0, total, []

    best_offset, best_count, best_non_matching = None, 0, []
    for offset in candidate_offsets:
        matched, non_matching = 0, []
        for acc, pos_start, pos_end, ref in positions:
            i0, i1 = offset + pos_start - 1, offset + pos_end
            if i1 > seq_len or ng_seq[i0:i1] != ref:
                non_matching.append(acc)
            else:
                matched += 1
        if matched > best_count:
            best_count, best_offset, best_non_matching = matched, offset, non_matching
        if best_count == total:
            break

    return best_offset, best_count, total, best_non_matching


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", choices=["LRG", "MANE"], default="LRG")
    args = parser.parse_args()

    ng_dir  = SOURCE_DIRS[args.source]
    out_csv = OUTPUT_CSVS[args.source]
    log_path = _PROC_DIR / "Logs" / f"find_lrg_offset_{args.source}.log"

    os.makedirs(str(out_csv.parent), exist_ok=True)
    os.makedirs(str(log_path.parent), exist_ok=True)

    print(f"Source: {args.source} -> {ng_dir}")

    with open(log_path, "w", encoding="utf-8") as logf:
        def log(msg):
            print(msg)
            logf.write(msg + "\n")

        df = pd.read_csv(MUTATIONS_TSV, sep="\t", dtype=str).fillna("")
        df = df[
            (df["variant_type"].str.strip().str.lower() == "genomic") &
            (df["mut_class"].str.strip().str.lower().isin(TARGET_MUT_CLASSES))
        ].copy()
        df["ref_nucleotides"] = df["notation"].apply(parse_ref_nucleotides)
        df["pos_start"] = pd.to_numeric(df["pos_start"], errors="coerce")
        df["pos_end"]   = pd.to_numeric(df["pos_end"],   errors="coerce")
        df = df.dropna(subset=["pos_start", "pos_end"])
        df = df[df["ref_nucleotides"].str.strip() != ""]

        log(f"Filtered genomic substitution rows with ref: {len(df)}")

        ng_files = {}
        for fname in os.listdir(ng_dir):
            if not fname.lower().endswith((".fasta", ".fa")):
                continue
            gene = fname.split("_")[0].upper()
            if gene not in ng_files:
                ng_files[gene] = Path(ng_dir) / fname

        genes = sorted(df["gene"].str.strip().str.upper().unique())
        results = []

        for gene in genes:
            def skip(status):
                results.append({"gene": gene, "ng_accession": "", "status": status,
                                 "sstart": "", "send": "", "matched": "", "total": "",
                                 "match_pct": "", "non_matching_accessions": ""})

            fasta_path = ng_files.get(gene)
            if fasta_path is None:
                log(f"  SKIP {gene}: no NG FASTA found"); skip("no_fasta"); continue

            ng_accession = ""
            m = re.search(r"(NG_\d+\.\d+)", fasta_path.name)
            if m:
                ng_accession = m.group(1)

            if fasta_path.stat().st_size > MAX_SEQ_LEN * 5:
                log(f"  SKIP {gene}: FASTA too large"); skip("too_large"); continue

            ng_seq = load_fasta(fasta_path)
            if len(ng_seq) > MAX_SEQ_LEN:
                log(f"  SKIP {gene}: sequence too long ({len(ng_seq)} bp)"); skip("too_large"); continue

            gene_muts = df[df["gene"].str.upper() == gene].to_dict("records")
            if not gene_muts:
                log(f"  SKIP {gene}: no mutations after filtering"); skip("no_mutations"); continue

            log(f"  {gene} ({ng_accession}): {len(gene_muts)} mutations, NG seq {len(ng_seq)} bp")
            best_offset, best_count, total, non_matching = find_best_offset(ng_seq, gene_muts)

            if best_offset is None:
                log(f"    no seed found after trying all mutations")
                results.append({"gene": gene, "ng_accession": ng_accession, "status": "seed_not_found",
                                 "sstart": "", "send": "", "matched": 0, "total": total,
                                 "match_pct": 0.0, "non_matching_accessions": ""})
                continue

            max_pos_end      = max(int(r["pos_end"]) for r in gene_muts)
            sstart           = best_offset + 1
            send             = best_offset + max_pos_end
            match_pct        = best_count / total
            n_non            = len(set(non_matching))
            status           = "ok" if (n_non <= 2 or match_pct >= MATCH_THRESHOLD) else "below_threshold"
            non_matching_str = ";".join(sorted(set(non_matching)))

            log(f"    offset={best_offset}  sstart={sstart}  send={send}  "
                f"match={best_count}/{total} ({match_pct:.1%})  [{status}]")
            if non_matching_str:
                log(f"    non-matching: {non_matching_str}")

            results.append({"gene": gene, "ng_accession": ng_accession, "status": status,
                             "sstart": sstart, "send": send, "matched": best_count, "total": total,
                             "match_pct": round(match_pct, 4), "non_matching_accessions": non_matching_str})

        pd.DataFrame(results).to_csv(out_csv, index=False)
        log(f"\nResults written to {out_csv}")


if __name__ == "__main__":
    main()
