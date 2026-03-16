import os
import re
import glob
import argparse
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SOURCE_2_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'DNA sequences', 'Reference sequence (Source 2)')
INPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'all_mutations', 'cross_check_rna_mutations.tsv')
OUTPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'reference_check', 'reference_check_source_2.tsv')

# c. positions are relative to CDS start (c.1 = first nucleotide of ATG).
# This script assumes the FASTA sequence also starts at the CDS (position 1).
# If your FASTA includes a 5'UTR, pass --utr-offset <length> to shift accordingly.


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


def find_fasta_files(gene: str) -> list[str]:
    found = set(
        glob.glob(os.path.join(SOURCE_2_DIR, f'*_{gene}.fasta'))
        + glob.glob(os.path.join(SOURCE_2_DIR, f'*_{gene}.FASTA'))
    )
    return sorted(found)


def sequence_matches(seq: str, pos_start: int, pos_end: int, ref_nucs: str,
                     utr_offset: int) -> bool:
    start_idx = pos_start - 1 + utr_offset
    end_idx = pos_end + utr_offset
    if start_idx < 0 or end_idx > len(seq):
        return False
    return seq[start_idx:end_idx] == ref_nucs.upper()


def check_gene(gene: str, gene_df: pd.DataFrame, utr_offset: int,
               seq_cache: dict) -> dict:
    fasta_paths = find_fasta_files(gene)
    total = len(gene_df)

    if not fasta_paths:
        return {
            'gene': gene,
            'total_mutations': total,
            'matching_count': 0,
            'non_matching_accessions': 'NO_FASTA',
            'percentage_match': 0.0,
        }

    best_match_count = -1
    best_non_matching = []

    for fasta_path in fasta_paths:
        if fasta_path not in seq_cache:
            seq_cache[fasta_path] = load_fasta_sequence(fasta_path)
        seq = seq_cache[fasta_path]

        match_count = 0
        non_matching = []

        for _, row in gene_df.iterrows():
            try:
                if sequence_matches(seq, int(row['c_pos_start']), int(row['c_pos_end']),
                                    str(row['ref_nucleotides']), utr_offset):
                    match_count += 1
                else:
                    non_matching.append(str(row['accession']))
            except Exception:
                non_matching.append(str(row['accession']))

        if match_count > best_match_count:
            best_match_count = match_count
            best_non_matching = non_matching

    assert best_match_count <= total, (
        f"Bug: match_count ({best_match_count}) > total ({total}) for gene {gene}. "
        "Regenerate cross_check_rna_mutations.tsv by re-running extract_del_sub_rna_mutations.py."
    )

    pct = round(best_match_count / total * 100, 2) if total > 0 else 0.0

    return {
        'gene': gene,
        'total_mutations': total,
        'matching_count': best_match_count,
        'non_matching_accessions': ';'.join(sorted(set(best_non_matching))) if best_non_matching else '',
        'percentage_match': pct,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Check mRNA FASTA reference sequences against c. mutations (source 2).')
    parser.add_argument('--utr-offset', type=int, default=0,
                        help='Length of 5\'UTR in the FASTA (default: 0, i.e. FASTA starts at CDS)')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

    df = pd.read_csv(INPUT_PATH, sep='\t')
    df_checkable = df[df['ref_nucleotides'].notna() & (df['ref_nucleotides'].astype(str) != '')].copy()

    print(f"Source 2: checking {len(df_checkable)} mutations across {df_checkable['gene'].nunique()} genes")
    if args.utr_offset:
        print(f"Applying 5'UTR offset: {args.utr_offset} bp")

    seq_cache = {}
    results = []

    for gene, gene_df in df_checkable.groupby('gene'):
        results.append(check_gene(gene, gene_df, args.utr_offset, seq_cache))

    out_df = pd.DataFrame(results, columns=[
        'gene', 'total_mutations', 'matching_count',
        'non_matching_accessions', 'percentage_match'
    ])
    out_df.to_csv(OUTPUT_PATH, sep='\t', index=False)

    total_mut = out_df['total_mutations'].sum()
    total_match = out_df['matching_count'].sum()
    overall_pct = round(total_match / total_mut * 100, 2) if total_mut > 0 else 0.0

    print(f"Overall match rate: {total_match}/{total_mut} ({overall_pct}%)")
    print(f"Saved to: {OUTPUT_PATH}")


if __name__ == '__main__':
    main()
