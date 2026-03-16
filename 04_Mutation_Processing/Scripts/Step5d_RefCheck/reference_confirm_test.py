import os
import re
import glob
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

FASTA_DIR = os.path.join(SCRIPT_DIR, '..', '..', '..', '05_Testing', 'reference_check')
INPUT_PATH = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'all_mutations', 'cross_check_mutations.tsv')
OUTPUT_PATH = os.path.join(FASTA_DIR, 'test_reference_check.tsv')


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


def find_fasta_files(gene: str) -> list[tuple[str, str]]:
    """Return list of (sequence_accession, fasta_path) for all versions of a gene."""
    pattern = os.path.join(FASTA_DIR, f'{gene}_*.fasta')
    paths = sorted(glob.glob(pattern))
    results = []
    for path in paths:
        basename = os.path.splitext(os.path.basename(path))[0]
        seq_accession = basename[len(gene) + 1:]
        results.append((seq_accession, path))
    return results


def sequence_matches(seq: str, pos_start: int, pos_end: int, ref_nucs: str) -> bool:
    start_idx = pos_start - 1
    end_idx = pos_end
    if start_idx < 0 or end_idx > len(seq):
        return False
    return seq[start_idx:end_idx] == ref_nucs.upper()


def check_gene(gene: str, gene_df: pd.DataFrame, seq_cache: dict) -> list[dict]:
    fasta_files = find_fasta_files(gene)
    total = len(gene_df)

    if not fasta_files:
        return [{
            'gene': gene,
            'sequence_accession': 'NO_FASTA',
            'total_mutations': total,
            'matching_count': 0,
            'non_matching_accessions': '',
            'percentage_match': 0.0,
        }]

    rows = []
    for seq_accession, fasta_path in fasta_files:
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
            f"Bug: match_count ({match_count}) > total ({total}) for {gene} / {seq_accession}. "
            "Regenerate cross_check_mutations.tsv by re-running extract_del_sub_mutations.py."
        )

        pct = round(match_count / total * 100, 2) if total > 0 else 0.0

        rows.append({
            'gene': gene,
            'sequence_accession': seq_accession,
            'total_mutations': total,
            'matching_count': match_count,
            'non_matching_accessions': ';'.join(non_matching) if non_matching else '',
            'percentage_match': pct,
        })

    return rows


def main():
    df = pd.read_csv(INPUT_PATH, sep='\t')
    df_checkable = df[df['ref_nucleotides'].notna() & (df['ref_nucleotides'].astype(str) != '')].copy()

    print(f"Checking {len(df_checkable)} mutations across {df_checkable['gene'].nunique()} genes")
    print(f"FASTA dir : {os.path.normpath(FASTA_DIR)}")

    seq_cache = {}
    results = []

    for gene, gene_df in df_checkable.groupby('gene'):
        results.extend(check_gene(gene, gene_df, seq_cache))

    out_df = pd.DataFrame(results, columns=[
        'gene', 'sequence_accession', 'total_mutations', 'matching_count',
        'non_matching_accessions', 'percentage_match'
    ])
    out_df.to_csv(OUTPUT_PATH, sep='\t', index=False)

    print(f"\nResults preview (>=90% match):")
    high_match = out_df[out_df['percentage_match'] >= 90]
    if high_match.empty:
        print("  None found.")
    else:
        for _, r in high_match.iterrows():
            print(f"  {r['gene']} / {r['sequence_accession']}: {r['matching_count']}/{r['total_mutations']} ({r['percentage_match']}%)")

    print(f"\nSaved to: {os.path.normpath(OUTPUT_PATH)}")


if __name__ == '__main__':
    main()
