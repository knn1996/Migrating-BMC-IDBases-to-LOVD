import os
import re
import argparse
import pandas as pd
from pathlib import Path

SOURCE_DIRS = {
    "IDRefseq_NM": Path(os.environ["SEQ_DIR"]),
    "Mane_NM":     Path(os.environ["SEQ_DIR"]),
}

OUTPUT_CSVS = {
    "IDRefseq_NM": Path(os.environ["OUT_CSV"]),
    "Mane_NM":     Path(os.environ["OUT_CSV"]),
}

MUTATIONS_TSV      = Path(os.environ["MUTATIONS_TSV"])
MAX_SEQ_LEN        = 100_000
MATCH_THRESHOLD    = 0.90
TARGET_MUT_CLASSES = {"substitution"}


def parse_ref_nucleotides(notation):
    if ">" in notation:
        m = re.search(r"([ACGTN]+)>[ACGTN]+$", notation, re.IGNORECASE)
        return m.group(1).upper() if m else ""
    return ""


def load_fasta(path):
    with open(path, encoding="utf-8") as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


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
    parser.add_argument("--source", choices=["IDRefseq_NM", "Mane_NM"], default="Mane_NM")
    args = parser.parse_args()

    nm_dir   = SOURCE_DIRS[args.source]
    out_csv  = OUTPUT_CSVS[args.source]
    log_path = Path(os.environ["LOG_PATH"])

    os.makedirs(str(out_csv.parent), exist_ok=True)
    os.makedirs(str(log_path.parent), exist_ok=True)

    print(f"Source: {args.source} -> {nm_dir}")

    with open(log_path, "w", encoding="utf-8") as logf:
        def log(msg):
            print(msg)
            logf.write(msg + "\n")

        df = pd.read_csv(MUTATIONS_TSV, sep="\t", dtype=str).fillna("")
        df = df[
            (df["variant_type"].str.strip().str.lower() == "coding") &
            (df["mut_class"].str.strip().str.lower().isin(TARGET_MUT_CLASSES))
        ].copy()
        df["ref_nucleotides"] = df["notation"].apply(parse_ref_nucleotides)
        df["pos_start"] = pd.to_numeric(df["pos_start"], errors="coerce")
        df["pos_end"]   = pd.to_numeric(df["pos_end"],   errors="coerce")
        df = df.dropna(subset=["pos_start", "pos_end"])
        df = df[df["ref_nucleotides"].str.strip() != ""]
        log(f"Filtered coding substitution rows with ref: {len(df)}")

        nm_files = {}
        for fname in os.listdir(nm_dir):
            if not fname.lower().endswith(".fasta"):
                continue
            m = re.match(r"^([A-Z0-9]+)_NM_", fname, re.IGNORECASE)
            if m:
                gene = m.group(1).upper()
                if gene not in nm_files:
                    nm_files[gene] = Path(nm_dir) / fname

        genes   = sorted(df["gene"].str.strip().str.upper().unique())
        results = []

        for gene in genes:
            def skip(status, gene=gene):
                results.append({"gene": gene, "nm_accession": "", "status": status,
                                 "sstart": "", "send": "", "matched": "", "total": "",
                                 "match_pct": "", "non_matching_accessions": ""})

            fasta_path = nm_files.get(gene)
            if fasta_path is None:
                log(f"  SKIP {gene}: no NM FASTA found"); skip("no_fasta"); continue

            m = re.search(r"(NM_\d+\.\d+)", fasta_path.name)
            nm_accession = m.group(1) if m else ""

            if fasta_path.stat().st_size > MAX_SEQ_LEN * 5:
                log(f"  SKIP {gene}: FASTA too large"); skip("too_large"); continue

            nm_seq = load_fasta(fasta_path)
            if len(nm_seq) > MAX_SEQ_LEN:
                log(f"  SKIP {gene}: sequence too long ({len(nm_seq)} bp)"); skip("too_large"); continue

            gene_muts = df[df["gene"].str.upper() == gene].to_dict("records")
            if not gene_muts:
                log(f"  SKIP {gene}: no mutations after filtering"); skip("no_mutations"); continue

            total = len(gene_muts)
            log(f"  {gene} ({nm_accession}): {total} mutations, NM seq {len(nm_seq)} bp")

            best_offset, best_count, _, best_non = find_best_offset(nm_seq, gene_muts)

            if best_offset is None:
                log(f"    seed_not_found")
                results.append({"gene": gene, "nm_accession": nm_accession, "status": "seed_not_found",
                                 "sstart": "", "send": "", "matched": 0, "total": total,
                                 "match_pct": 0.0, "non_matching_accessions": ""})
                continue

            if not is_ok(best_count, total, best_non):
                for delta in (+1, -1):
                    offset, count, _, non = find_best_offset(nm_seq, shift_mutations(gene_muts, delta))
                    if offset is not None and count > best_count:
                        best_offset, best_count, best_non = offset, count, non
                        log(f"    shift ({delta:+d}) improved: {count}/{total} ({count/total:.1%})")
                        if is_ok(count, total, non):
                            break

            match_pct        = best_count / total
            ok               = is_ok(best_count, total, best_non)
            status           = "ok" if ok else "below_threshold"
            sstart           = best_offset + 1
            send             = len(nm_seq)
            non_matching_str = ";".join(sorted(set(best_non)))

            log(f"    offset={best_offset}  sstart={sstart}  send={send}  "
                f"match={best_count}/{total} ({match_pct:.1%})  [{status}]")
            if non_matching_str:
                log(f"    non-matching: {non_matching_str}")

            results.append({"gene": gene, "nm_accession": nm_accession, "status": status,
                             "sstart": sstart, "send": send, "matched": best_count, "total": total,
                             "match_pct": round(match_pct, 4),
                             "non_matching_accessions": non_matching_str})

        pd.DataFrame(results).to_csv(out_csv, index=False)
        log(f"\nResults written to {out_csv}")


if __name__ == "__main__":
    main()