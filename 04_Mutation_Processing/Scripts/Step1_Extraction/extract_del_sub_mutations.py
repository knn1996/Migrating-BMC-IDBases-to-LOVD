import os
import re
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'all_mutations', 'all_mutations.tsv')
OUTPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'all_mutations', 'cross_check_mutations.tsv')

TARGET_FEAT_TYPES = {'point', 'deletion', 'indel', 'complex'}
SKIP_MUT_TYPES = {'delins', 'unknown'}


def extract_all_g_notations(sysname: str) -> list[str]:
    return ['g.' + m.rstrip(',;') for m in re.findall(r'g\.([^\s,;]+)', sysname)]


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
    df = df[df['feat_name'].isin(TARGET_FEAT_TYPES)].copy()

    rows = []
    seen = set()

    for _, row in df.iterrows():
        g_notations = extract_all_g_notations(str(row['sysname_per_allele']))
        if not g_notations:
            continue

        for g_notation in g_notations:
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
