import os
import pandas as pd
from pathlib import Path

CROSS_CHECK_TSV = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\all_mutations\cross_check_mutations.tsv"
LRG_NG_DIR      = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\LRG FASTA file (Source 3)\NG"
OUT_CSV         = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\lrg_offset_results.csv"

MAX_SEQ_LEN     = 400_000
MATCH_THRESHOLD = 0.90


def load_fasta(path):
    seq_lines = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                seq_lines.append(line.strip().upper())
    return "".join(seq_lines)


def find_best_offset(lrg_seq, mutations):
    seq_len   = len(lrg_seq)
    positions = [(int(r["pos_start"]), int(r["pos_end"]), r["ref_nucleotides"].upper())
                 for r in mutations]
    total     = len(positions)

    seed_pos_start, seed_pos_end, seed_ref = max(positions, key=lambda x: len(x[2]))

    candidate_offsets = []
    search_start = 0
    while True:
        idx = lrg_seq.find(seed_ref, search_start)
        if idx == -1:
            break
        offset = idx - (seed_pos_start - 1)
        if offset >= 0:
            candidate_offsets.append(offset)
        search_start = idx + 1

    print(f"    seed='{seed_ref[:20]}{'...' if len(seed_ref)>20 else ''}' "
          f"len={len(seed_ref)}  candidates={len(candidate_offsets)}")

    if not candidate_offsets:
        return None, 0, total

    best_offset = None
    best_count  = 0

    for offset in candidate_offsets:
        matched = 0
        for pos_start, pos_end, ref in positions:
            idx_start = offset + pos_start - 1
            idx_end   = offset + pos_end
            if idx_end > seq_len:
                continue
            if lrg_seq[idx_start:idx_end] == ref:
                matched += 1

        if matched > best_count:
            best_count  = matched
            best_offset = offset

        if best_count == total:
            break

    return best_offset, best_count, total


def main():
    df = pd.read_csv(CROSS_CHECK_TSV, sep="\t", dtype=str).fillna("")
    df["pos_start"] = pd.to_numeric(df["pos_start"], errors="coerce")
    df["pos_end"]   = pd.to_numeric(df["pos_end"],   errors="coerce")
    df = df.dropna(subset=["pos_start", "pos_end", "ref_nucleotides"])
    df = df[df["ref_nucleotides"].str.strip() != ""]

    fasta_files = {
        Path(p).stem.upper(): Path(LRG_NG_DIR) / p
        for p in os.listdir(LRG_NG_DIR)
        if p.lower().endswith((".fasta", ".fa"))
    }

    genes   = df["gene"].str.strip().str.upper().unique()
    results = []

    for gene in sorted(genes):
        fasta_path = fasta_files.get(gene)
        if fasta_path is None:
            print(f"  SKIP {gene}: no FASTA found")
            results.append({"gene": gene, "status": "no_fasta",
                            "sstart": "", "send": "", "matched": "", "total": "", "match_pct": ""})
            continue

        fasta_size = fasta_path.stat().st_size
        if fasta_size > MAX_SEQ_LEN * 1.1:
            print(f"  SKIP {gene}: FASTA too large ({fasta_size/1000:.0f} kb > {MAX_SEQ_LEN//1000} kb limit)")
            results.append({"gene": gene, "status": "too_large",
                            "sstart": "", "send": "", "matched": "", "total": "", "match_pct": ""})
            continue

        lrg_seq = load_fasta(fasta_path)
        if len(lrg_seq) > MAX_SEQ_LEN:
            print(f"  SKIP {gene}: sequence length {len(lrg_seq)} > {MAX_SEQ_LEN} limit")
            results.append({"gene": gene, "status": "too_large",
                            "sstart": "", "send": "", "matched": "", "total": "", "match_pct": ""})
            continue

        gene_muts = df[df["gene"].str.upper() == gene].to_dict("records")
        if not gene_muts:
            print(f"  SKIP {gene}: no mutations after filtering")
            results.append({"gene": gene, "status": "no_mutations",
                            "sstart": "", "send": "", "matched": "", "total": "", "match_pct": ""})
            continue

        print(f"  {gene}: {len(gene_muts)} mutations, LRG seq {len(lrg_seq)} bp")

        best_offset, best_count, total = find_best_offset(lrg_seq, gene_muts)

        if best_offset is None:
            print(f"    seed not found in LRG sequence")
            results.append({"gene": gene, "status": "seed_not_found",
                            "sstart": "", "send": "", "matched": 0, "total": total, "match_pct": 0.0})
            continue

        max_pos_end = max(int(r["pos_end"]) for r in gene_muts)
        sstart      = best_offset + 1
        send        = best_offset + max_pos_end
        match_pct   = best_count / total

        status = "ok" if match_pct >= MATCH_THRESHOLD else "below_threshold"
        print(f"    offset={best_offset}  sstart={sstart}  send={send}  "
              f"match={best_count}/{total} ({match_pct:.1%})  [{status}]")

        results.append({
            "gene":      gene,
            "status":    status,
            "sstart":    sstart,
            "send":      send,
            "matched":   best_count,
            "total":     total,
            "match_pct": round(match_pct, 4),
        })

    pd.DataFrame(results).to_csv(OUT_CSV, index=False)
    print(f"\nResults written to {OUT_CSV}")


if __name__ == "__main__":
    main()
