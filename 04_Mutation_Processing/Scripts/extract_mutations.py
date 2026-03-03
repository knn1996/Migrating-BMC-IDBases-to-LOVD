"""
extract_mutations.py
====================
Extracts all DNA mutations from every *pub.html file in the idbase directory
and writes them to a single TSV: all_mutations.tsv

Output columns:
  gene        — gene symbol (derived from folder name, e.g. ADAbase → ADA)
  accession   — patient/allele accession (e.g. A0022)
  sysname     — Systematic name if present (e.g. g.4037C>T, c.967C>T, p.S291L)
  pos_start   — IDRefSeq DNA position start (1-based)
  pos_end     — IDRefSeq DNA position end (1-based; same as pos_start for SNVs)
  feat_name   — mutation type from /name field (point, deletion, insertion, ...)
  allele_num  — allele number within the entry (1, 2, ...)

Usage:
  python3 extract_mutations.py
"""

import os
import re
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────────────────────
THESIS_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
IDBASE_DIR = os.path.join(THESIS_DIR, "02_Source_Database", "idbase")
OUT_PATH   = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output", "all_mutations.tsv")

RE_STRIP = re.compile(r"<[^>]+>")
RE_LOC   = re.compile(r"/loc:\s+IDRefSeq:\s+\w+:\s*(\d+)(?:\.\.(\d+))?", re.IGNORECASE)

# ── Parser ─────────────────────────────────────────────────────────────────────
def parse_mutations(pub_html_path, gene):
    """Parse all Feature dna /loc entries from one pub.html.
    Yields dicts with gene, accession, sysname, pos_start, pos_end,
    feat_name, allele_num.
    """
    lines = RE_STRIP.sub(
        "", Path(pub_html_path).read_text(encoding="utf-8", errors="replace")
    ).splitlines()

    acc = sys_name = feat_name = ""
    in_dna = False
    allele_num = 0
    seen = set()

    for line in lines:
        s = line.strip()

        if s.startswith("Accession"):
            m = re.match(r"Accession\s+(\S+)", s)
            if m:
                acc, sys_name, feat_name, in_dna, allele_num = m.group(1), "", "", False, 0

        elif s.startswith("Systematic name"):
            m = re.match(r"Systematic name\s+(.+)", s)
            if m:
                # Append additional allele names (some entries have multiple lines)
                part = m.group(1).strip()
                sys_name = f"{sys_name}; {part}" if sys_name else part

        elif re.match(r"FeatureHeader\s+allele;", s):
            m = re.match(r"FeatureHeader\s+allele;\s*(\d+)", s)
            if m:
                allele_num = int(m.group(1))

        elif re.match(r"Feature\s+dna;", s):
            in_dna, feat_name = True, ""

        elif re.match(r"Feature\s+(rna|aa);", s):
            in_dna = False

        elif in_dna:
            m = re.match(r"/name:\s+(\S+)", s)
            if m:
                feat_name = m.group(1)

            m = RE_LOC.search(s)
            if m:
                p1 = int(m.group(1))
                p2 = int(m.group(2)) if m.group(2) else p1
                key = (acc, allele_num, p1, p2)
                if key not in seen:
                    seen.add(key)
                    yield {
                        "gene":       gene,
                        "accession":  acc,
                        "sysname":    sys_name,
                        "pos_start":  p1,
                        "pos_end":    p2,
                        "feat_name":  feat_name,
                        "allele_num": allele_num,
                    }
                # in_dna intentionally NOT reset — captures all /loc per dna block


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    folders = sorted(
        d for d in os.listdir(IDBASE_DIR)
        if os.path.isdir(os.path.join(IDBASE_DIR, d))
        and d.lower().endswith("base")
        and d.lower() != "immunomebase"
    )

    rows = []
    for folder in folders:
        gene = re.sub(r"base$", "", folder, flags=re.IGNORECASE)
        folder_path = os.path.join(IDBASE_DIR, folder)
        pub = next(
            (f for f in os.listdir(folder_path) if f.lower().endswith("pub.html")),
            None
        )
        if not pub:
            print(f"  SKIP {gene}: no pub.html")
            continue

        pub_path = os.path.join(folder_path, pub)
        gene_rows = list(parse_mutations(pub_path, gene))
        rows.extend(gene_rows)
        print(f"  {gene:<15}  {len(gene_rows):4d} mutations  ({pub})")

    # Write TSV
    header = "gene\taccession\tsysname\tpos_start\tpos_end\tfeat_name\tallele_num\n"
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        f.write(header)
        for r in rows:
            f.write(
                f"{r['gene']}\t{r['accession']}\t{r['sysname']}\t"
                f"{r['pos_start']}\t{r['pos_end']}\t{r['feat_name']}\t{r['allele_num']}\n"
            )

    print(f"\nTotal mutations : {len(rows)}")
    print(f"Output written  : {OUT_PATH}")


if __name__ == "__main__":
    main()
