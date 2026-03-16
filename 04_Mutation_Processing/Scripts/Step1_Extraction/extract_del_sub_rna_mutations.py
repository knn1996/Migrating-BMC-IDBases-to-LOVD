import os
import re
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'all_mutations', 'all_mutations.tsv')
OUTPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'all_mutations', 'cross_check_rna_mutations.tsv')

TARGET_FEAT_TYPES = {'point', 'deletion', 'indel', 'complex'}
SKIP_MUT_TYPES = {'delins', 'unknown'}


def extract_c_notation(sysname: str) -> str | None:
    match = re.search(r'c\.([^\s,;]+)', sysname)
    if match:
        return 'c.' + match.group(1).rstrip(',;')
    return None


def is_non_exonic(c_notation: str) -> bool:
    """Skip intronic (+/- offset), 5'UTR (negative), 3'UTR (*), and EX-style notations."""
    body = c_notation[2:]
    if re.match(r'-\d', body):
        return True
    if re.match(r'\*', body):
        return True
    if re.search(r'\d[+-]\d', body):
        return True
    if re.match(r'EX', body, re.IGNORECASE):
        return True
    return False


def parse_c_positions(c_notation: str) -> tuple[int, int] | None:
    """
    Extract (pos_start, pos_end) from c. notation (1-based, exonic only).
    Returns None if position cannot be parsed.
    """
    body = c_notation[2:]
    range_match = re.match(r'(\d+)_(\d+)', body)
    if range_match:
        return int(range_match.group(1)), int(range_match.group(2))
    single_match = re.match(r'(\d+)', body)
    if single_match:
        pos = int(single_match.group(1))
        return pos, pos
    return None


def parse_ref_nucleotides(c_notation: str) -> tuple[str | None, str]:
    if not c_notation:
        return None, 'unknown'

    if 'delins' in c_notation:
        return None, 'delins'

    if '>' in c_notation:
        match = re.search(r'([ACGTN]+)>[ACGTN]+$', c_notation, re.IGNORECASE)
        if match:
            return match.group(1).upper(), 'substitution'
        return None, 'substitution'

    if 'del' in c_notation:
        match = re.search(r'del([ACGTN]+)$', c_notation, re.IGNORECASE)
        if match:
            return match.group(1).upper(), 'deletion'
        return None, 'deletion'

    return None, 'unknown'


def main():
    df = pd.read_csv(INPUT_PATH, sep='\t')
    df = df[df['feat_name'].isin(TARGET_FEAT_TYPES)].copy()

    rows = []
    seen = set()
    skipped_non_exonic = 0
    skipped_no_pos = 0

    for _, row in df.iterrows():
        c_notation = extract_c_notation(str(row['sysname_per_allele']))
        if not c_notation:
            continue

        dedup_key = (row['gene'], c_notation)
        if dedup_key in seen:
            continue
        seen.add(dedup_key)

        if is_non_exonic(c_notation):
            skipped_non_exonic += 1
            continue

        ref_nucs, mut_type = parse_ref_nucleotides(c_notation)

        if mut_type in SKIP_MUT_TYPES:
            continue

        positions = parse_c_positions(c_notation)
        if positions is None:
            skipped_no_pos += 1
            continue

        c_pos_start, c_pos_end = positions

        rows.append({
            'gene': row['gene'],
            'accession': row['accession'],
            'c_notation': c_notation,
            'c_pos_start': c_pos_start,
            'c_pos_end': c_pos_end,
            'ref_nucleotides': ref_nucs if ref_nucs else '',
            'mut_type': mut_type,
        })

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    out_df.to_csv(OUTPUT_PATH, sep='\t', index=False)

    total = len(out_df)
    checkable = (out_df['ref_nucleotides'] != '').sum()
    print(f"Total exonic mutations extracted : {total}")
    print(f"  With ref_nucleotides           : {checkable}")
    print(f"  Without ref_nucleotides        : {total - checkable}  (large deletions with no nucleotide specified)")
    print(f"Skipped (intronic/UTR/EX)        : {skipped_non_exonic}")
    print(f"Skipped (unparseable position)   : {skipped_no_pos}")
    print(f"Saved to: {OUTPUT_PATH}")


if __name__ == '__main__':
    main()
