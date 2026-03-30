import os
import csv
from pathlib import Path
from collections import defaultdict

_SCRIPT_DIR  = Path(__file__).parent
_THESIS_DIR  = (_SCRIPT_DIR / ".." / ".." / "..").resolve()
_PROC_DIR    = _THESIS_DIR / "04_Mutation_Processing"

IDREFSEQ_DIR = _PROC_DIR / "DNA sequences" / "IDRefseq"
OUT_CSV      = _PROC_DIR / "Output" / "Step2_RefCheck" / "idrefseq_pair_comparison.csv"


def load_seq(path):
    with open(path, encoding="utf-8") as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


def compare(s1, s3):
    if s1 == s3:
        return "identical"
    if s1 in s3:
        return "s1_inside_s3"
    if s3 in s1:
        return "s3_inside_s1"
    return "no_containment"


def main():
    files = [f for f in os.listdir(IDREFSEQ_DIR) if f.endswith(".fasta")]

    by_gene = defaultdict(dict)
    for fname in files:
        if fname.endswith("_source_1.fasta"):
            gene = fname.split("_")[0]
            by_gene[gene]["source_1"] = fname
        elif fname.endswith("_source_3.fasta"):
            gene = fname.split("_")[0]
            by_gene[gene]["source_3"] = fname

    pairs = {g: v for g, v in by_gene.items() if "source_1" in v and "source_3" in v}
    print(f"Found {len(pairs)} gene pairs to compare\n")

    rows = []
    for gene in sorted(pairs):
        f1 = pairs[gene]["source_1"]
        f3 = pairs[gene]["source_3"]

        s1 = load_seq(Path(IDREFSEQ_DIR) / f1)
        s3 = load_seq(Path(IDREFSEQ_DIR) / f3)

        result = compare(s1, s3)
        offset = s3.index(s1) if result == "s1_inside_s3" else \
                 s1.index(s3) if result == "s3_inside_s1" else ""

        print(f"  {gene:<12}  s1={len(s1)} bp  s3={len(s3)} bp  → {result}"
              + (f"  offset={offset}" if offset != "" else ""))

        rows.append({
            "gene":              gene,
            "source_1":          f1,
            "source_3":          f3,
            "len_s1":            len(s1),
            "len_s3":            len(s3),
            "result":            result,
            "offset_in_longer":  offset,
        })

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "source_1", "source_3", "len_s1", "len_s3", "result", "offset_in_longer"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nResults written to {OUT_CSV}")
    counts = {}
    for r in rows:
        counts[r["result"]] = counts.get(r["result"], 0) + 1
    for k, v in sorted(counts.items()):
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
