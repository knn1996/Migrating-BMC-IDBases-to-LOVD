import csv
import os
from pathlib import Path

OFFSET_CSV = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step2_RefCheck\lrg_offset_results_NM.csv"
SOURCE3_NM = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\LRG FASTA file (Source 3)\NM"
OUT_DIR    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\IDRefseq_NM"

MATCH_PCT_THRESHOLD = 0.90
MAX_MISMATCHES      = 2

os.makedirs(OUT_DIR, exist_ok=True)


def load_fasta_seq(path):
    with open(path, encoding="utf-8") as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


def make_header(gene, nm_acc, sstart, send, match_pct, seq_len):
    return (f">{gene} ref={nm_acc} sstart={sstart} send={send} "
            f"match_pct={match_pct} length={seq_len}")


def write_fasta(path, header, seq, width=60):
    with open(path, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i+width] + "\n")


def build_fasta_index(directory):
    index = {}
    for p in Path(directory).iterdir():
        if p.suffix.lower() not in (".fasta", ".fa"):
            continue
        parts = p.stem.split("_", 1)
        if len(parts) == 2:
            gene, nm_acc = parts[0].upper(), parts[1]
            index[gene] = (p, nm_acc)
    return index


def main():
    fasta_index = build_fasta_index(SOURCE3_NM)

    with open(OFFSET_CSV, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        gene      = row["gene"].strip().upper()
        status    = row.get("status", "").strip()
        sstart    = row.get("sstart", "").strip()
        send      = row.get("send", "").strip()
        match_pct = row.get("match_pct", "").strip()

        if status not in ("ok", "below_threshold") or not sstart or not send:
            print(f"  SKIP {gene}: status={status}")
            continue

        try:
            matched    = int(row.get("matched", 0))
            total      = int(row.get("total", 0))
            pct        = float(match_pct) if match_pct else 0.0
            mismatches = total - matched
        except ValueError:
            print(f"  SKIP {gene}: cannot parse matched/total/match_pct")
            continue

        if pct < MATCH_PCT_THRESHOLD and mismatches > MAX_MISMATCHES:
            print(f"  SKIP {gene}: match_pct={pct:.2%} mismatches={mismatches} (below threshold)")
            continue

        if gene not in fasta_index:
            print(f"  SKIP {gene}: FASTA not found in {SOURCE3_NM}")
            continue

        fasta_path, nm_acc = fasta_index[gene]
        full_seq = load_fasta_seq(fasta_path)
        s, e     = int(sstart), int(send)
        subseq   = full_seq[s - 1 : e]

        out_name = f"{gene}_{nm_acc}.fasta"
        header   = make_header(gene, nm_acc, sstart, send, match_pct, len(subseq))
        write_fasta(Path(OUT_DIR) / out_name, header, subseq)
        print(f"  WROTE {out_name}  ({len(subseq)} bp, positions {s}-{e} of {len(full_seq)})  "
              f"match={pct:.2%} mismatches={mismatches}")


if __name__ == "__main__":
    main()
