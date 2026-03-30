import os
import re
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'Step1_Extraction', 'all_mutations.tsv')
OUTPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'Step1_Extraction', 'cross_check_mutations.tsv')

TARGET_MUT_CLASSES = {'substitution', 'deletion', 'insertion', 'indel'}
SKIP_MUT_TYPES = {'delins', 'unknown'}


def parse_positions_from_g_notation(g_notation: str, fallback_start: int, fallback_end: int) -> tuple[int, int]:
    range_match = re.search(r'g\.(\d+)_(\d+)', g_notation)
    if range_match:
        return int(range_match.group(1)), int(range_match.group(2))
    single_match = re.search(r'g\.(\d+)', g_notation)
    if single_match:
        pos = int(single_match.group(1))
        return pos, pos
    return fallback_start, fallback_end


def parse_ref_nucleotides(g_notation: str) -> tuple[str | None, str]:
    if not g_notation:
        return None, 'unknown'

    if 'delins' in g_notation:
        return None, 'delins'

    if '>' in g_notation:
        match = re.search(r'([ACGTN]+)>[ACGTN]+$', g_notation, re.IGNORECASE)
        if match:
            return match.group(1).upper(), 'substitution'
        return None, 'substitution'

    if 'del' in g_notation:
        match = re.search(r'del([ACGTN]+)$', g_notation, re.IGNORECASE)
        if match:
            return match.group(1).upper(), 'deletion'
        return None, 'deletion'

    return None, 'unknown'


def main():
    df = pd.read_csv(INPUT_PATH, sep='\t')
    df = df[df['mut_class'].isin(TARGET_MUT_CLASSES)].copy()
    df = df[df['variant_type'] == 'genomic'].copy()

    rows = []
    seen = set()

    for _, row in df.iterrows():
        g_notation = str(row['notation']).strip()
        if not g_notation or not g_notation.startswith('g.'):
            continue

        dedup_key = (row['gene'], g_notation)
        if dedup_key in seen:
            continue
        seen.add(dedup_key)

        ref_nucs, mut_type = parse_ref_nucleotides(g_notation)

        if mut_type in SKIP_MUT_TYPES:
            continue

        pos_start, pos_end = parse_positions_from_g_notation(
            g_notation, int(row['pos_start']), int(row['pos_end'])
        )

        rows.append({
            'gene': row['gene'],
            'accession': row['accession'],
            'g_notation': g_notation,
            'pos_start': pos_start,
            'pos_end': pos_end,
            'ref_nucleotides': ref_nucs if ref_nucs else '',
            'mut_type': mut_type,
        })

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    out_df.to_csv(OUTPUT_PATH, sep='\t', index=False)

    total = len(out_df)
    checkable = (out_df['ref_nucleotides'] != '').sum()
    print(f"Total mutations extracted: {total}")
    print(f"  With ref_nucleotides   : {checkable}")
    print(f"  Without ref_nucleotides: {total - checkable}  (large deletions with no nucleotide specified)")
    print(f"Saved to: {OUTPUT_PATH}")


if __name__ == '__main__':
    main()
