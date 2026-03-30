"""
extract_mutations.py
====================
Extracts all DNA mutations from every *pub.html in the idbase directory.

Output TSV columns:
  gene         — gene symbol (ADAbase → ADA)
  accession    — patient accession (e.g. A0022)
  allele_num   — 0 (unknown/single), 1, or 2
  variant_type — genomic / coding / rna / protein / mitochondrial / non_coding / unknown
  notation     — individual HGVS notation for this row
  pos_start    — extracted from notation; for genomic rows falls back to IDRefSeq /loc
  pos_end      — same logic as pos_start
  mut_class    — substitution / deletion / insertion / indel / duplication / unknown
                 derived from the genomic (g.) notation and applied to all levels of
                 the same accession so c. and p. rows share the same classification
  missing_flag — comma-separated missing levels among g./c./p., empty if all present
"""

import os
import re
import sys
import logging
from pathlib import Path
from datetime import datetime

THESIS_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
IDBASE_DIR = os.path.join(THESIS_DIR, "02_Source_Database", "idbase")
OUT_PATH   = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output", "Step1_Extraction", "all_mutations.tsv")
LOG_PATH   = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Logs", "extract_mutations_log.txt")

RE_STRIP = re.compile(r"<[^>]+>")
RE_LOC   = re.compile(r"/loc:.*?:\s*(\d+)(?:\.\.(\d+))?", re.IGNORECASE)

VARIANT_PREFIXES = [
    ("g.", "genomic"),
    ("c.", "coding"),
    ("r.", "rna"),
    ("p.", "protein"),
    ("m.", "mitochondrial"),
    ("n.", "non_coding"),
]

REQUIRED_LEVELS = {"genomic", "coding", "protein"}


def setup_logging():
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(LOG_PATH, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger(__name__)


def detect_variant_type(notation):
    for prefix, vtype in VARIANT_PREFIXES:
        if notation.startswith(prefix):
            return vtype
    return "unknown"


def classify_mutation(notation):
    """Classify based on the notation string. Applied from g. notation when available."""
    n = notation.lower()
    if "delins" in n:
        return "indel"
    if "dup" in n:
        return "duplication"
    if "ins" in n:
        return "insertion"
    if "del" in n:
        return "deletion"
    if ">" in n:
        return "substitution"
    return "unknown"


def extract_pos_from_notation(notation):
    """Extract (pos_start, pos_end) from a notation string. Returns (None, None) on failure."""
    if notation.startswith("p."):
        m = re.search(r"[A-Za-z*]+(\d+)", notation[2:])
        if m:
            pos = int(m.group(1))
            m2  = re.search(r"_[A-Za-z*]*(\d+)", notation[2:])
            return pos, int(m2.group(1)) if m2 else pos
        return None, None
    m = re.search(r"[gcrnm]\.[-*]?(\d+)", notation)
    if m:
        pos = int(m.group(1))
        m2  = re.search(r"_(\d+)", notation)
        return pos, int(m2.group(1)) if m2 else pos
    return None, None


def parse_allele_label(sysname):
    m = re.match(r"Allele\s+1\s+and\s+2\s*:\s*(.+)", sysname, re.IGNORECASE)
    if m:
        return "both", m.group(1).strip()
    m = re.match(r"Allele\s+(\d+)\s*:\s*(.+)", sysname, re.IGNORECASE)
    if m:
        return m.group(1), m.group(2).strip()
    return None, sysname


def split_by_level(sysname_clean):
    """Parse a cleaned sysname into {variant_type: [token, ...]} dict."""
    levels = {}
    for seg in sysname_clean.split(","):
        tokens = [t.strip() for t in re.split(r";+", seg) if t.strip()]
        if not tokens:
            continue
        vtype = detect_variant_type(tokens[0])
        if vtype != "unknown":
            levels[vtype] = tokens
    return levels


def compute_missing_flag(levels):
    missing = REQUIRED_LEVELS - set(levels.keys())
    return ",".join(sorted(missing)) if missing else ""


def expand_rows(gene, accession, sysname, genomic_loc_start, genomic_loc_end):
    """
    Yield one output row per (allele, variant_type) combination.

    Position logic:
      - genomic rows: extract from g. notation, fall back to IDRefSeq /loc
      - all other rows (coding, protein, rna, ...): extract from their own notation only;
        do NOT fall back to genomic /loc (those coordinates are meaningless for c./p.)

    mut_class logic:
      - derived from the genomic (g.) notation when present, shared across all levels
      - falls back to per-notation classification only when no g. notation exists
    """
    allele_mode, clean = parse_allele_label(sysname)
    levels       = split_by_level(clean) if clean else {}
    missing_flag = compute_missing_flag(levels)

    # Determine shared mut_class from g. notation if available
    genomic_tokens = levels.get("genomic", [])
    shared_mut_class = classify_mutation(genomic_tokens[0]) if genomic_tokens else None

    if not levels:
        ps, pe = extract_pos_from_notation(clean)
        yield {
            "gene": gene, "accession": accession, "allele_num": 0,
            "variant_type": "unknown", "notation": clean,
            "pos_start": ps if ps is not None else genomic_loc_start,
            "pos_end":   pe if pe is not None else genomic_loc_end,
            "mut_class": classify_mutation(clean),
            "missing_flag": "genomic,coding,protein",
        }
        return

    if allele_mode == "both":
        allele_pairs = [(1, 0), (2, 1)]
    elif allele_mode is not None:
        allele_pairs = [(int(allele_mode), 0)]
    else:
        allele_pairs = [(0, 0)]

    for vtype, tokens in levels.items():
        for allele_num, token_idx in allele_pairs:
            notation = tokens[token_idx] if token_idx < len(tokens) else tokens[-1]
            ps, pe   = extract_pos_from_notation(notation)

            # Only genomic rows may fall back to IDRefSeq /loc coordinates
            if vtype == "genomic":
                ps = ps if ps is not None else genomic_loc_start
                pe = pe if pe is not None else genomic_loc_end

            mut_class = shared_mut_class if shared_mut_class else classify_mutation(notation)

            yield {
                "gene": gene, "accession": accession, "allele_num": allele_num,
                "variant_type": vtype, "notation": notation,
                "pos_start": ps,
                "pos_end":   pe,
                "mut_class": mut_class,
                "missing_flag": missing_flag,
            }


def parse_pub_html(pub_html_path, gene, log):
    """
    Parse one pub.html and yield output rows.

    Emission strategy: one emission per accession, triggered at the '//' record
    separator. IDRefSeq /loc lines are collected during parsing and the FIRST
    genomic DNA /loc is used as the positional fallback for that accession.
    This avoids duplicates caused by emitting once per /loc line.
    """
    lines = RE_STRIP.sub(
        "", Path(pub_html_path).read_text(encoding="utf-8", errors="replace")
    ).splitlines()

    acc      = ""
    sys_name = ""
    loc_start = None
    loc_end   = None
    emitted   = set()
    total     = 0

    for line in lines:
        s = line.strip()

        if s.startswith("Accession"):
            m = re.match(r"Accession\s+(\S+)", s)
            if m:
                acc       = m.group(1)
                sys_name  = ""
                loc_start = None
                loc_end   = None

        elif s.startswith("Systematic name"):
            m = re.match(r"Systematic name\s+(.+)", s)
            if m:
                part = m.group(1).strip()
                sys_name = f"{sys_name}; {part}" if sys_name else part

        elif RE_LOC.search(s):
            # Capture only the first /loc per accession as the genomic position fallback
            if loc_start is None:
                m = RE_LOC.search(s)
                loc_start = int(m.group(1))
                loc_end   = int(m.group(2)) if m.group(2) else loc_start

        elif s == "//":
            if acc and sys_name and acc not in emitted:
                emitted.add(acc)
                total += 1
                yield from expand_rows(gene, acc, sys_name, loc_start, loc_end)

    if total == 0:
        log.warning(f"  {gene}: 0 entries emitted — check pub.html structure")


def main():
    log   = setup_logging()
    start = datetime.now()

    log.info("=" * 60)
    log.info("extract_mutations.py — starting")
    log.info(f"IDBASE_DIR : {IDBASE_DIR}")
    log.info(f"OUT_PATH   : {OUT_PATH}")
    log.info(f"Started at : {start.strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("=" * 60)

    folders = sorted(
        d for d in os.listdir(IDBASE_DIR)
        if os.path.isdir(os.path.join(IDBASE_DIR, d))
        and d.lower().endswith("base")
        and d.lower() != "immunomebase"
    )

    log.info(f"Found {len(folders)} gene folders\n")

    rows           = []
    skipped_genes  = []
    zero_row_genes = []

    for folder in folders:
        gene        = re.sub(r"base$", "", folder, flags=re.IGNORECASE)
        folder_path = os.path.join(IDBASE_DIR, folder)
        pub         = next((f for f in os.listdir(folder_path) if f.lower().endswith("pub.html")), None)

        if not pub:
            log.warning(f"  SKIP {gene}: no pub.html found")
            skipped_genes.append(gene)
            continue

        gene_rows = list(parse_pub_html(os.path.join(folder_path, pub), gene, log))
        rows.extend(gene_rows)

        missing_flags = sum(1 for r in gene_rows if r["missing_flag"])
        unknown_types = sum(1 for r in gene_rows if r["variant_type"] == "unknown")

        log.info(
            f"  {gene:<15}  {len(gene_rows):5d} rows"
            + (f"  [{missing_flags} missing_flag]" if missing_flags else "")
            + (f"  [{unknown_types} unknown_type]"  if unknown_types  else "")
        )

        if not gene_rows:
            zero_row_genes.append(gene)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    header = "gene\taccession\tallele_num\tvariant_type\tnotation\tpos_start\tpos_end\tmut_class\tmissing_flag\n"
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        f.write(header)
        for r in rows:
            f.write(
                f"{r['gene']}\t{r['accession']}\t{r['allele_num']}\t{r['variant_type']}\t"
                f"{r['notation']}\t{r['pos_start']}\t{r['pos_end']}\t{r['mut_class']}\t{r['missing_flag']}\n"
            )

    elapsed = datetime.now() - start
    log.info("")
    log.info("=" * 60)
    log.info("SUMMARY")
    log.info("=" * 60)
    log.info(f"  Genes processed   : {len(folders) - len(skipped_genes)}")
    log.info(f"  Genes skipped     : {len(skipped_genes)}  {skipped_genes or ''}")
    log.info(f"  Genes with 0 rows : {len(zero_row_genes)}  {zero_row_genes or ''}")
    log.info(f"  Total rows        : {len(rows)}")
    log.info(f"  Output            : {OUT_PATH}")
    log.info(f"  Log               : {LOG_PATH}")
    log.info(f"  Elapsed           : {str(elapsed).split('.')[0]}")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
