import os
import re
import argparse
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SOURCE_DIRS = {
    1: os.path.join(SCRIPT_DIR, '..', '..', 'DNA sequences', 'Processed FASTA file (Source 1)'),
    3: os.path.join(SCRIPT_DIR, '..', '..', 'DNA sequences', 'LRG FASTA file (Source 3)', 'NG'),
}

INPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'all_mutations', 'cross_check_mutations.tsv')
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'reference_check')


def load_fasta_sequence(fasta_path: str) -> str:
    seq_parts = []
    with open(fasta_path, 'r') as f:
        recording = False
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if recording:
                    break
                recording = True
            elif recording:
                seq_parts.append(re.sub(r'[^ACGTN]', '', line.upper()))
    return ''.join(seq_parts)


def find_fasta_file(gene: str, base_dir: str) -> str | None:
    for ext in (f'{gene}.fasta', f'{gene}.FASTA'):
        path = os.path.join(base_dir, ext)
        if os.path.exists(path):
            return path
    return None


def sequence_matches(seq: str, pos_start: int, pos_end: int, ref_nucs: str) -> bool:
    start_idx = pos_start - 1
    end_idx = pos_end
    if start_idx < 0 or end_idx > len(seq):
        return False
    return seq[start_idx:end_idx] == ref_nucs.upper()


def check_gene(gene: str, gene_df: pd.DataFrame, base_dir: str, seq_cache: dict) -> dict:
    fasta_path = find_fasta_file(gene, base_dir)
    total = len(gene_df)

    if not fasta_path:
        return {
            'gene': gene,
            'total_mutations': total,
            'matching_count': 0,
            'non_matching_accessions': 'NO_FASTA',
            'percentage_match': 0.0,
        }

    if fasta_path not in seq_cache:
        seq_cache[fasta_path] = load_fasta_sequence(fasta_path)
    seq = seq_cache[fasta_path]

    match_count = 0
    non_matching = []

    for _, row in gene_df.iterrows():
        try:
            if sequence_matches(seq, int(row['pos_start']), int(row['pos_end']),
                                str(row['ref_nucleotides'])):
                match_count += 1
            else:
                non_matching.append(str(row['accession']))
        except Exception:
            non_matching.append(str(row['accession']))

    assert match_count <= total, (
        f"Bug: match_count ({match_count}) > total ({total}) for gene {gene}. "
        "Regenerate cross_check_mutations.tsv by re-running extract_del_sub_mutations.py."
    )

    pct = round(match_count / total * 100, 2) if total > 0 else 0.0

    return {
        'gene': gene,
        'total_mutations': total,
        'matching_count': match_count,
        'non_matching_accessions': ';'.join(sorted(set(non_matching))) if non_matching else '',
        'percentage_match': pct,
    }


def main():
    parser = argparse.ArgumentParser(description='Check genomic FASTA reference sequences (source 1 or 3).')
    parser.add_argument('--source', type=int, required=True, choices=[1, 3],
                        help='Source number: 1=Processed FASTA, 3=LRG NG')
    args = parser.parse_args()

    base_dir = SOURCE_DIRS[args.source]
    output_path = os.path.join(OUTPUT_DIR, f'reference_check_source_{args.source}.tsv')
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_PATH, sep='\t')
    df_checkable = df[df['ref_nucleotides'].notna() & (df['ref_nucleotides'].astype(str) != '')].copy()

    print(f"Source {args.source}: checking {len(df_checkable)} mutations across {df_checkable['gene'].nunique()} genes")

    seq_cache = {}
    results = []

    for gene, gene_df in df_checkable.groupby('gene'):
        results.append(check_gene(gene, gene_df, base_dir, seq_cache))

    out_df = pd.DataFrame(results, columns=[
        'gene', 'total_mutations', 'matching_count',
        'non_matching_accessions', 'percentage_match'
    ])
    out_df.to_csv(output_path, sep='\t', index=False)

    total_mut = out_df['total_mutations'].sum()
    total_match = out_df['matching_count'].sum()
    overall_pct = round(total_match / total_mut * 100, 2) if total_mut > 0 else 0.0

    print(f"Overall match rate: {total_match}/{total_mut} ({overall_pct}%)")
    print(f"Saved to: {output_path}")


if __name__ == '__main__':
    main()
