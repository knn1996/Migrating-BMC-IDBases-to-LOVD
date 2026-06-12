import os
import re
import argparse
import pandas as pd
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
_THESIS_DIR = (_SCRIPT_DIR / ".." / ".." / "..").resolve()
_PROC_DIR   = _THESIS_DIR / "04_Mutation_Processing"
_SEQ_DIR    = _PROC_DIR / "DNA sequences"

SOURCE_DIRS = {
    "IDRefseq": _SEQ_DIR / "IDRefseq",
    "Mane_NG":  _SEQ_DIR / "Mane_Select_NG",
}

OUTPUT_CSVS = {
    "IDRefseq": _PROC_DIR / "Output" / "Step2_RefCheck" / "lrg_offset_results.csv",
    "Mane_NG":  _PROC_DIR / "Output" / "Step2_RefCheck" / "MANE_offset_results.csv",
}

MUTATIONS_TSV   = _PROC_DIR / "Output" / "Step1_Extraction" / "all_mutations.tsv"
MAX_SEQ_LEN     = 400_000
MATCH_THRESHOLD = 0.90

SUB_ONLY_CLASSES = {"substitution"}
ALL_ANCHOR_CLASSES = {"substitution", "deletion", "indel", "duplication"}

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def parse_ref_nucleotides(notation):
    notation = notation.strip()
    m = re.search(r"([ACGTN]+)>[ACGTN]+$", notation, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.search(r"del([ACGTN]+)ins[ACGTN]+$", notation, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.search(r"del([ACGTN]+)$", notation, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    m = re.search(r"dup([ACGTN]+)$", notation, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    return ""


def load_fasta(path):
    with open(path, encoding="utf-8") as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def find_best_offset(seq, mutations):
    seq_len   = len(seq)
    positions = [(r["accession"], int(r["pos_start"]), int(r["pos_end"]), r["ref_nucleotides"].upper())
                 for r in mutations]
    total = len(positions)

    candidate_offsets = []
    for _, seed_start, _, seed_ref in sorted(positions, key=lambda x: len(x[3]), reverse=True):
        candidate_offsets = []
        idx = 0
        while True:
            idx = seq.find(seed_ref, idx)
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

    best_offset, best_count, best_non = None, 0, []
    for offset in candidate_offsets:
        matched, non_matching = 0, []
        for acc, pos_start, pos_end, ref in positions:
            i0, i1 = offset + pos_start - 1, offset + pos_end
            if i1 > seq_len or seq[i0:i1] != ref:
                non_matching.append(acc)
            else:
                matched += 1
        if matched > best_count:
            best_count, best_offset, best_non = matched, offset, non_matching
        if best_count == total:
            break

    return best_offset, best_count, total, best_non


def is_ok(count, total, non_matching):
    return len(set(non_matching)) <= 2 or (count / total) >= MATCH_THRESHOLD


def shift_mutations(mutations, delta):
    return [{**r, "pos_start": r["pos_start"] + delta, "pos_end": r["pos_end"] + delta}
            for r in mutations]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", choices=["IDRefseq", "Mane_NG"], default="IDRefseq")
    parser.add_argument("--include-indels", action="store_true",
                        help="Anchor on substitutions + deletions + indels + duplications (default: subs only)")
    args = parser.parse_args()

    target_classes = ALL_ANCHOR_CLASSES if args.include_indels else SUB_ONLY_CLASSES
    suffix         = "_with_indels" if args.include_indels else ""

    ng_dir       = SOURCE_DIRS[args.source]
    base_out     = OUTPUT_CSVS[args.source]
    out_csv      = base_out.with_name(base_out.stem + suffix + base_out.suffix)
    log_path     = _PROC_DIR / "Logs" / f"find_lrg_offset_{args.source}{suffix}.log"
    allow_rc     = args.source != "IDRefseq"

    os.makedirs(str(out_csv.parent), exist_ok=True)
    os.makedirs(str(log_path.parent), exist_ok=True)

    print(f"Source: {args.source} -> {ng_dir}")
    print(f"Anchor types: {sorted(target_classes)}")
    print(f"Reverse complement pass: {'enabled' if allow_rc else 'disabled (IDRefseq pre-extracted)'}")
    print(f"Output: {out_csv}")

    with open(log_path, "w", encoding="utf-8") as logf:
        def log(msg):
            print(msg)
            logf.write(msg + "\n")

        df = pd.read_csv(MUTATIONS_TSV, sep="\t", dtype=str).fillna("")
        df = df[
            (df["variant_type"].str.strip().str.lower() == "genomic") &
            (df["mut_class"].str.strip().str.lower().isin(target_classes))
        ].copy()
        df["ref_nucleotides"] = df["notation"].apply(parse_ref_nucleotides)
        df["pos_start"] = pd.to_numeric(df["pos_start"], errors="coerce")
        df["pos_end"]   = pd.to_numeric(df["pos_end"],   errors="coerce")
        df = df.dropna(subset=["pos_start", "pos_end"])
        df = df[df["ref_nucleotides"].str.strip() != ""]
        log(f"Filtered genomic anchorable rows with parseable ref: {len(df)}")

        if args.include_indels:
            class_counts = df["mut_class"].str.lower().value_counts().to_dict()
            log(f"  Anchor type breakdown: {class_counts}")

        ng_files = {}
        for fname in os.listdir(ng_dir):
            if not fname.lower().endswith((".fasta", ".fa")):
                continue
            gene = fname.split("_")[0].upper()
            if gene not in ng_files:
                ng_files[gene] = Path(ng_dir) / fname

        genes   = sorted(df["gene"].str.strip().str.upper().unique())
        results = []

        for gene in genes:
            def skip(status, gene=gene):
                results.append({"gene": gene, "ng_accession": "", "status": status,
                                 "strand": "", "sstart": "", "send": "",
                                 "matched": "", "total": "", "match_pct": "",
                                 "non_matching_accessions": ""})

            fasta_path = ng_files.get(gene)
            if fasta_path is None:
                log(f"  SKIP {gene}: no NG FASTA found"); skip("no_fasta"); continue

            m = re.search(r"(NG_\d+\.\d+)", fasta_path.name)
            ng_accession = m.group(1) if m else ""

            if fasta_path.stat().st_size > MAX_SEQ_LEN * 5:
                log(f"  SKIP {gene}: FASTA too large"); skip("too_large"); continue

            ng_seq = load_fasta(fasta_path)
            if len(ng_seq) > MAX_SEQ_LEN:
                log(f"  SKIP {gene}: sequence too long ({len(ng_seq)} bp)"); skip("too_large"); continue

            gene_muts = df[df["gene"].str.upper() == gene].to_dict("records")
            if not gene_muts:
                log(f"  SKIP {gene}: no mutations after filtering"); skip("no_mutations"); continue

            total = len(gene_muts)
            log(f"  {gene} ({ng_accession}): {total} anchors, NG seq {len(ng_seq)} bp")

            best_offset = None
            best_count  = 0
            best_strand = "+"
            best_non    = []

            strands = [("+", ng_seq)]
            if allow_rc:
                strands.append(("-", reverse_complement(ng_seq)))

            for strand, seq in strands:
                if best_count == total:
                    break
                offset, count, _, non = find_best_offset(seq, gene_muts)
                if offset is not None and count > best_count:
                    best_offset, best_count, best_strand, best_non = offset, count, strand, non

            if best_offset is None:
                log(f"    seed_not_found after all strand passes")
                results.append({"gene": gene, "ng_accession": ng_accession, "status": "seed_not_found",
                                 "strand": "", "sstart": "", "send": "",
                                 "matched": 0, "total": total, "match_pct": 0.0,
                                 "non_matching_accessions": ""})
                continue

            if not is_ok(best_count, total, best_non):
                best_seq = reverse_complement(ng_seq) if best_strand == "-" else ng_seq
                for delta in (+1, -1):
                    offset, count, _, non = find_best_offset(best_seq, shift_mutations(gene_muts, delta))
                    if offset is not None and count > best_count:
                        best_offset, best_count, best_non = offset, count, non
                        log(f"    shift ({delta:+d}) on {best_strand} strand improved: {count}/{total} ({count/total:.1%})")
                        if is_ok(count, total, non):
                            break

            match_pct        = best_count / total
            ok               = is_ok(best_count, total, best_non)
            status           = ("ok_reverse" if best_strand == "-" else "ok") if ok else "below_threshold"
            max_pos_end      = max(int(r["pos_end"]) for r in gene_muts)
            sstart           = best_offset + 1
            send             = best_offset + max_pos_end
            non_matching_str = ";".join(sorted(set(best_non)))

            log(f"    offset={best_offset}  sstart={sstart}  send={send}  "
                f"strand={best_strand}  match={best_count}/{total} ({match_pct:.1%})  [{status}]")
            if non_matching_str:
                log(f"    non-matching: {non_matching_str}")

            results.append({"gene": gene, "ng_accession": ng_accession, "status": status,
                             "strand": best_strand, "sstart": sstart, "send": send,
                             "matched": best_count, "total": total,
                             "match_pct": round(match_pct, 4),
                             "non_matching_accessions": non_matching_str})

        pd.DataFrame(results).to_csv(out_csv, index=False)
        log(f"\nResults written to {out_csv}")


if __name__ == "__main__":
    main()
