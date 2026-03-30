import csv
import os
from pathlib import Path
from collections import defaultdict

_SCRIPT_DIR       = Path(__file__).parent
_THESIS_DIR       = (_SCRIPT_DIR / ".." / ".." / "..").resolve()
_PROC_DIR         = _THESIS_DIR / "04_Mutation_Processing"

REFERENCE_SUMMARY = _PROC_DIR / "Output" / "Step2_RefCheck" / "reference_summary.csv"
OUT_DIR           = _PROC_DIR / "DNA sequences" / "IDRefseq"

os.makedirs(OUT_DIR, exist_ok=True)


def load_fasta_seq(path):
    with open(path, encoding="utf-8") as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


def make_header(gene, ref, sstart, send, match_pct, non_matching, seq_len, source):
    nm = non_matching if non_matching else "none"
    return (f">{gene} ref={ref} source={source} sstart={sstart} send={send} "
            f"match_pct={match_pct} non_matching={nm} length={seq_len}")


def write_fasta(path, header, seq, width=60):
    with open(path, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i+width] + "\n")


def pick_row(rows):
    if len(rows) == 1:
        return rows[0]
    def sort_key(r):
        pct = float(r["match_pct"]) if r["match_pct"] else 0.0
        src_priority = 1 if r["source"] == "source_3" else 0
        return (pct, src_priority)
    return max(rows, key=sort_key)


def main():
    with open(REFERENCE_SUMMARY, encoding="utf-8") as f:
        all_rows = list(csv.DictReader(f))

    by_gene = defaultdict(list)
    for row in all_rows:
        by_gene[row["gene"]].append(row)

    for gene, rows in sorted(by_gene.items()):
        row          = pick_row(rows)
        source       = row["source"]
        ref          = row["ref"]
        sstart       = row["sstart"]
        send         = row["send"]
        match_pct    = row["match_pct"]
        non_matching = row["non_matching"]

        if len(rows) > 1:
            skipped = [r["source"] for r in rows if r is not row]
            print(f"  {gene}: picked {source} (match_pct={match_pct}), skipped {skipped}")

        if source == "source_1":
            src_path = Path(SOURCE1_DIR) / f"{gene}.FASTA"
            if not src_path.exists():
                src_path = Path(SOURCE1_DIR) / f"{gene}.fasta"
            if not src_path.exists():
                print(f"  SKIP {gene} source_1: FASTA not found"); continue

            seq      = load_fasta_seq(src_path)
            label    = ref if ref != "not_available_yet" else "NA"
            out_name = f"{gene}_{label}_source_1.fasta"
            header   = make_header(gene, ref, sstart, send, match_pct, non_matching, len(seq), "source_1")
            write_fasta(Path(OUT_DIR) / out_name, header, seq)
            print(f"  WROTE {out_name}  ({len(seq)} bp)")

        elif source == "source_3":
            src_path = Path(SOURCE3_DIR) / f"{gene}.fasta"
            if not src_path.exists():
                src_path = Path(SOURCE3_DIR) / f"{gene}.FASTA"
            if not src_path.exists():
                print(f"  SKIP {gene} source_3: FASTA not found"); continue

            full_seq = load_fasta_seq(src_path)
            s, e     = int(sstart), int(send)
            subseq   = full_seq[s - 1 : e]

            out_name = f"{gene}_{ref}_source_3.fasta"
            header   = make_header(gene, ref, sstart, send, match_pct, non_matching, len(subseq), "source_3")
            write_fasta(Path(OUT_DIR) / out_name, header, subseq)
            print(f"  WROTE {out_name}  ({len(subseq)} bp, positions {s}-{e} of {len(full_seq)})")


if __name__ == "__main__":
    main()
