"""
extract_mutations.py (v2)
=========================
Extracts all DNA mutations from every *pub.html file in the idbase directory
and writes them to a single TSV: all_mutations.tsv

Additionally classifies and separates three known data quality issues:

  Problem 1 – Incomplete substitution variants (all_mutations_no_subs.tsv):
      Sysname has a DNA/RNA change ending with '>' but no target nucleotide,
      e.g. 'g.206889>', 'c.13242+1>', 'r.13242+1>'. These are ambiguous –
      they may be deletions whose sequence is recorded in the Feature /change
      field rather than in the sysname.  An extra column 'is_sequence_in_feature'
      records the concatenated /change string from the Feature dna block when it
      contains a nucleotide sequence of length > 3; otherwise "None".

  Problem 2 – Long-sequence mutations in Feature (all_mutations_long_seq.tsv):
      Sysname is exactly 'g.c.r.' (placeholder), meaning the full mutation
      sequence was too long for the sysname field and was placed in the Feature
      /change lines instead.  Rows already in Problem 1 are NOT duplicated here.

  Problem 3 – Sysname semicolon / multi-allele splitting (all_mutations.tsv):
      Many sysnames encode both alleles in one string separated by '; Allele N:'.
      The cleaned output splits the sysname per-allele into 'sysname_per_allele':
        - 'Allele 1: X; Allele 2: Y' → two rows with 'Allele 1: X' / 'Allele 2: Y'
        - 'Allele 1 and 2: X' → two rows with 'Allele 1: X' / 'Allele 2: X'
      Irregular cases (e.g. 'Allele 1 and 2:' that actually encodes two different
      interleaved mutations) are flagged in 'irregular_sysname_flag'.
      Problem-1 and Problem-2 rows are excluded from this clean output.
      The raw sysname column is NOT included in the clean output.

Output columns (base):
  gene             — gene symbol (derived from folder name)
  accession        — patient/allele accession
  sysname          — raw Systematic name as recorded in the source file
  pos_start        — IDRefSeq DNA position start (1-based)
  pos_end          — IDRefSeq DNA position end (same as pos_start for SNVs)
  feat_name        — mutation type from /name field
  allele_num       — allele number within the entry (0 = single, 1/2 = multi)

Extra columns in all_mutations_no_subs.tsv:
  is_sequence_in_feature — concatenated /change text if it contains >3 nt chars,
                           else "None"

Columns in all_mutations.tsv (Output/all_mutations/) — sysname is NOT included:
  sysname_per_allele    — per-allele portion of sysname after splitting
  irregular_sysname_flag — True if sysname was flagged as ambiguously encoded

Usage:
  python3 extract_mutations.py
"""

import os
import re
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
THESIS_DIR  = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", ".."))
IDBASE_DIR  = os.path.join(THESIS_DIR, "02_Source_Database", "idbase")

OUT_DIR     = os.path.join(_SCRIPT_DIR, "..", "..", "Output")
OUT_SUB_DIR = os.path.join(OUT_DIR, "all_mutations")

# Original output path (kept for backward compatibility)
PATH_MAIN      = os.path.join(OUT_DIR,     "all_mutations.tsv")
# New output paths in the subfolder
PATH_CLEAN     = os.path.join(OUT_SUB_DIR, "all_mutations.tsv")
PATH_NO_SUBS   = os.path.join(OUT_SUB_DIR, "all_mutations_no_subs.tsv")
PATH_LONG_SEQ  = os.path.join(OUT_SUB_DIR, "all_mutations_long_seq.tsv")

# Minimum nucleotide-character count in /change text to be considered "long"
LONG_SEQ_NT_THRESHOLD = 3

RE_STRIP = re.compile(r"<[^>]+>")
RE_LOC   = re.compile(r"/loc:\s+IDRefSeq:\s+\w+:\s*(\d+)(?:\.\.(\d+))?", re.IGNORECASE)


# ── Problem detection helpers ──────────────────────────────────────────────────

def is_problem1(sysname: str) -> bool:
    """
    Return True if the sysname contains an incomplete substitution variant –
    a DNA/RNA change notation that ends with '>' without a target nucleotide.
    Examples: g.206889>  c.13242+1>  r.13242+1>  g.IVS5-25GGTATTGCCC>
    """
    # Match [gcr]. followed by any non-whitespace/non-separator chars ending in '>'
    # The '>' must be followed by a comma, semicolon, whitespace, or end-of-string
    return bool(re.search(
        r'[gcr]\.[^\s,;]*[A-Za-z0-9]>(?=[,;\s]|\Z)',
        sysname,
        re.IGNORECASE
    ))


def is_problem2(sysname: str) -> bool:
    """
    Return True if sysname is the 'g.c.r.' placeholder indicating the full
    mutation sequence was placed in the Feature /change lines.
    """
    return bool(re.match(r'^g\.\s*c\.\s*r\.?\s*$', sysname.strip(), re.IGNORECASE))


def compute_is_sequence_in_feature(dna_change: str) -> str:
    """
    Examine the concatenated Feature dna /change text.
    If it contains more than LONG_SEQ_NT_THRESHOLD nucleotide characters,
    return the change string as-is (trimmed).  Otherwise return "None".
    """
    if not dna_change:
        return "None"
    # Count nucleotide chars (ACGTU, case-insensitive)
    nt_count = len(re.sub(r'[^ACGTUacgtu]', '', dna_change))
    return dna_change.strip() if nt_count > LONG_SEQ_NT_THRESHOLD else "None"


# ── Sysname per-allele splitting (Problem 3) ──────────────────────────────────

def split_sysname_by_allele(sysname: str, allele_num: int, pos_start: int | None = None):
    """
    Split a combined sysname into the portion relevant to the given allele_num.

    Returns (per_allele_sysname: str, is_irregular: bool).

    is_irregular is True for cases where the 'Allele 1 and 2:' label was used
    but the content actually encodes two *different* interleaved mutations
    (e.g. SMARCAL1 S0020, AIRE A0211).  In those cases the original sysname is
    returned unchanged and the row is flagged.
    """
    # No allele marker at all → return as-is, not irregular
    if not re.search(r'Allele\s+\d', sysname, re.IGNORECASE):
        return sysname, False

    # ── Case A: "Allele X and Y: [content]" ───────────────────────────────────
    m_and = re.match(
        r'^Allele\s+(\d+)\s+and\s+(\d+)\s*:\s*(.*)',
        sysname,
        re.IGNORECASE | re.DOTALL
    )
    if m_and:
        content = m_and.group(3)

        # Irregular check 1: content also contains an explicit "Allele N:" label
        if re.search(r';\s*Allele\s+\d+\s*:', content, re.IGNORECASE):
            return sysname, True

        # Irregular check 2: multiple distinct g. positions separated by ';'
        # (interleaved allele notation)
        g_entries = re.findall(r'g\.[^\s,;]+', content, re.IGNORECASE)
        if len(set(g_entries)) > 1:
            return sysname, True

        # Normal "Allele 1 and 2" – same mutation for both alleles.
        # Label each row with its specific allele number so the sysname_per_allele
        # column is clean: "Allele 1: <content>" / "Allele 2: <content>".
        label = allele_num if allele_num > 0 else int(m_and.group(1))
        return f"Allele {label}: {content}", False

    # ── Case B: "Allele 1: [data]; Allele 2: [data]" ─────────────────────────
    if re.search(r'Allele\s+\d+\s*:', sysname, re.IGNORECASE):
        # Split at positions where a ';' is immediately followed by 'Allele N:'
        # Use lookahead so the delimiter content is included in the second part
        parts = re.split(r';\s*(?=Allele\s+\d+\s*:)', sysname, flags=re.IGNORECASE)

        allele_entries = {}
        for part in parts:
            part_clean = part.strip().lstrip(';').strip()
            m = re.match(r'Allele\s+(\d+)\s*:\s*(.*)', part_clean,
                         re.IGNORECASE | re.DOTALL)
            if m:
                num = int(m.group(1))
                body = m.group(2).strip().rstrip(',; ')
                entry = f"Allele {num}: {body}"
                allele_entries.setdefault(num, []).append(entry)

        if len(allele_entries) > 1 and allele_num in allele_entries:
            candidates = allele_entries[allele_num]
            if len(candidates) == 1:
                return candidates[0], False
            if pos_start is not None:
                for cand in candidates:
                    if re.search(rf'g\.{pos_start}(?![0-9_])', cand):
                        return cand, False
            body_parts = [re.sub(r'^Allele\s+\d+\s*:\s*', '', c) for c in candidates]
            return f"Allele {allele_num}: {'; '.join(body_parts)}", False

    # Could not split meaningfully → return original
    return sysname, False


# ── Block-based parser ─────────────────────────────────────────────────────────

def _emit_entry(acc, sys_name, allele_blocks, gene, seen):
    """
    Yield row dicts for a fully-parsed entry.
    allele_blocks: dict  allele_num → list of dna-block dicts
                         {name, loc:(p1,p2), changes:[str]}
    """
    rows = []
    for an in sorted(allele_blocks):
        for blk in allele_blocks[an]:
            if blk.get('loc'):
                p1, p2 = blk['loc']
                key = (acc, an, p1, p2)
                if key not in seen:
                    seen.add(key)
                    rows.append({
                        'gene':       gene,
                        'accession':  acc,
                        'sysname':    sys_name,
                        'pos_start':  p1,
                        'pos_end':    p2,
                        'feat_name':  blk.get('name', ''),
                        'allele_num': an,
                        'dna_change': "".join(blk.get('changes', [])),
                    })
    return rows


def parse_mutations(pub_html_path, gene):
    """
    Parse all Feature dna /loc entries from one pub.html.
    Yields dicts with gene, accession, sysname, pos_start, pos_end,
    feat_name, allele_num, dna_change (internal – used for is_sequence_in_feature).
    """
    lines = RE_STRIP.sub(
        "", Path(pub_html_path).read_text(encoding="utf-8", errors="replace")
    ).splitlines()

    seen = set()

    # Entry-level state
    acc       = ""
    sys_name  = ""
    allele_num    = 0
    allele_blocks = {}   # allele_num → [blk dict]

    # Feature-block state
    in_dna  = False
    cur_dna = None       # {'name': str, 'loc': (p1,p2)|None, 'changes': [str]}

    def _close_dna_block():
        """Finalise the current dna block and store it."""
        nonlocal in_dna, cur_dna
        if in_dna and cur_dna is not None and cur_dna.get('loc'):
            allele_blocks.setdefault(allele_num, []).append(cur_dna)
        in_dna  = False
        cur_dna = None

    for line in lines:
        s = line.strip()

        # ── New accession → flush previous entry ──────────────────────────────
        if s.startswith("Accession"):
            m = re.match(r"Accession\s+(\S+)", s)
            if m:
                _close_dna_block()
                if acc:
                    yield from _emit_entry(acc, sys_name, allele_blocks, gene, seen)
                # Reset for new entry
                acc           = m.group(1)
                sys_name      = ""
                allele_num    = 0
                allele_blocks = {}
                in_dna        = False
                cur_dna       = None

        elif s.startswith("Systematic name"):
            m = re.match(r"Systematic name\s+(.+)", s)
            if m:
                part = m.group(1).strip()
                sys_name = f"{sys_name}; {part}" if sys_name else part

        elif re.match(r"FeatureHeader\s+allele;", s):
            m = re.match(r"FeatureHeader\s+allele;\s*(\d+)", s)
            if m:
                _close_dna_block()
                allele_num = int(m.group(1))

        elif re.match(r"Feature\s+dna;", s):
            # Close any open dna block before starting a new one
            _close_dna_block()
            in_dna  = True
            cur_dna = {'name': '', 'loc': None, 'changes': []}

        elif re.match(r"Feature\s+(rna|aa);", s):
            _close_dna_block()

        elif in_dna and cur_dna is not None:
            # NOTE: lines look like "Feature           /name: point" after stripping,
            # so we need re.search (not re.match) for /name: and /change: patterns.
            nm = re.search(r"/name:\s+(\S+)", s)
            if nm:
                cur_dna['name'] = nm.group(1)

            lm = RE_LOC.search(s)
            if lm:
                p1 = int(lm.group(1))
                p2 = int(lm.group(2)) if lm.group(2) else p1
                cur_dna['loc'] = (p1, p2)

            cm = re.search(r"/change:\s*(.*)", s)
            if cm:
                cur_dna['changes'].append(cm.group(1).strip())

    # Flush last entry
    _close_dna_block()
    if acc:
        yield from _emit_entry(acc, sys_name, allele_blocks, gene, seen)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUT_SUB_DIR, exist_ok=True)

    folders = sorted(
        d for d in os.listdir(IDBASE_DIR)
        if os.path.isdir(os.path.join(IDBASE_DIR, d))
        and d.lower().endswith("base")
        and d.lower() != "immunomebase"
    )

    all_rows = []
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

        pub_path  = os.path.join(folder_path, pub)
        gene_rows = list(parse_mutations(pub_path, gene))
        all_rows.extend(gene_rows)
        print(f"  {gene:<15}  {len(gene_rows):4d} mutations  ({pub})")

    # ── Classify rows ──────────────────────────────────────────────────────────
    prob1_rows     = []   # Problem 1: incomplete substitution
    prob2_rows     = []   # Problem 2: g.c.r. placeholder (not already in prob1)
    prob1_keys     = set()
    clean_rows     = []   # Normal rows (problems excluded)

    for r in all_rows:
        key = (r['gene'], r['accession'], r['allele_num'], r['pos_start'], r['pos_end'])

        if is_problem1(r['sysname']):
            r['is_sequence_in_feature'] = compute_is_sequence_in_feature(r['dna_change'])
            prob1_rows.append(r)
            prob1_keys.add(key)

        elif is_problem2(r['sysname']):
            prob2_rows.append(r)

        else:
            clean_rows.append(r)

    # ── Apply sysname splitting (Problem 3) to clean rows ─────────────────────
    irregular_count = 0
    for r in clean_rows:
        per_allele, is_irr = split_sysname_by_allele(r['sysname'], r['allele_num'], r['pos_start'])
        r['sysname_per_allele']     = per_allele
        r['irregular_sysname_flag'] = is_irr
        if is_irr:
            irregular_count += 1

    # ── Write original all_mutations.tsv (unchanged format, all rows) ─────────
    base_header = "gene\taccession\tsysname\tpos_start\tpos_end\tfeat_name\tallele_num\n"
    with open(PATH_MAIN, "w", encoding="utf-8") as f:
        f.write(base_header)
        for r in all_rows:
            f.write(
                f"{r['gene']}\t{r['accession']}\t{r['sysname']}\t"
                f"{r['pos_start']}\t{r['pos_end']}\t{r['feat_name']}\t{r['allele_num']}\n"
            )

    # ── Write all_mutations_no_subs.tsv (Problem 1) ───────────────────────────
    nosubs_header = (
        "gene\taccession\tsysname\tpos_start\tpos_end\t"
        "feat_name\tallele_num\tis_sequence_in_feature\n"
    )
    with open(PATH_NO_SUBS, "w", encoding="utf-8") as f:
        f.write(nosubs_header)
        for r in prob1_rows:
            f.write(
                f"{r['gene']}\t{r['accession']}\t{r['sysname']}\t"
                f"{r['pos_start']}\t{r['pos_end']}\t{r['feat_name']}\t"
                f"{r['allele_num']}\t{r['is_sequence_in_feature']}\n"
            )

    # ── Write all_mutations_long_seq.tsv (Problem 2, no overlap with Problem 1) ─
    with open(PATH_LONG_SEQ, "w", encoding="utf-8") as f:
        f.write(base_header)
        for r in prob2_rows:
            f.write(
                f"{r['gene']}\t{r['accession']}\t{r['sysname']}\t"
                f"{r['pos_start']}\t{r['pos_end']}\t{r['feat_name']}\t{r['allele_num']}\n"
            )

    # ── Write clean all_mutations.tsv (Output/all_mutations/) ─────────────────
    # sysname is intentionally omitted – sysname_per_allele is the clean field.
    clean_header = (
        "gene\taccession\tsysname_per_allele\t"
        "pos_start\tpos_end\tfeat_name\tallele_num\tirregular_sysname_flag\n"
    )
    with open(PATH_CLEAN, "w", encoding="utf-8") as f:
        f.write(clean_header)
        for r in clean_rows:
            f.write(
                f"{r['gene']}\t{r['accession']}\t"
                f"{r['sysname_per_allele']}\t"
                f"{r['pos_start']}\t{r['pos_end']}\t{r['feat_name']}\t"
                f"{r['allele_num']}\t{r['irregular_sysname_flag']}\n"
            )

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print(f"Total mutations extracted  : {len(all_rows)}")
    print()
    print(f"Problem 1 – Incomplete substitution (trailing '>') :")
    print(f"  Rows extracted           : {len(prob1_rows)}")
    p1_seq = sum(1 for r in prob1_rows if r['is_sequence_in_feature'] != 'None')
    p1_none = len(prob1_rows) - p1_seq
    print(f"  → with sequence in Feature : {p1_seq}")
    print(f"  → no sequence in Feature   : {p1_none}")
    print(f"  → Output: {PATH_NO_SUBS}")
    print()
    print(f"Problem 2 – g.c.r. placeholder (long seq in Feature):")
    print(f"  Rows extracted           : {len(prob2_rows)}")
    print(f"  → Output: {PATH_LONG_SEQ}")
    print()
    print(f"Problem 3 – Sysname semicolon / allele splitting:")
    semi_total = sum(1 for r in clean_rows if ';' in r['sysname'] and
                     re.search(r'Allele', r['sysname'], re.IGNORECASE))
    print(f"  Rows with allele-split sysname : {semi_total}")
    print(f"  Irregular sysname flags        : {irregular_count}")
    print()
    print(f"Clean all_mutations.tsv (problems excluded):")
    print(f"  Rows written             : {len(clean_rows)}")
    print(f"  → Output: {PATH_CLEAN}")
    print()
    print(f"Original all_mutations.tsv (all rows, unchanged format):")
    print(f"  → Output: {PATH_MAIN}")
    print("=" * 60)


if __name__ == "__main__":
    main()
