import os
import re
import glob
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SOURCE_2_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'DNA sequences', 'Reference sequences (Source 2)')
INPUT_PATH   = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'Step1_Extraction', 'all_mutations.tsv')
OUTPUT_PATH  = os.path.join(SCRIPT_DIR, '..', '..', 'Output', 'Step2_RefCheck', 'reference_check_source_2.tsv')


def fasta_length(fasta_path: str) -> int:
    length = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                length += len(re.sub(r'[^ACGTN]', '', line.strip().upper()))
    return length


def find_fasta_files(gene: str) -> list[str]:
    found = set(
        glob.glob(os.path.join(SOURCE_2_DIR, f'{gene}_*.fasta'))
        + glob.glob(os.path.join(SOURCE_2_DIR, f'{gene}_*.FASTA'))
    )
    return sorted(found)


def check_gene(gene: str, gene_df: pd.DataFrame, seq_len_cache: dict) -> dict:
    fasta_paths = find_fasta_files(gene)
    total = len(gene_df)

    if not fasta_paths:
        return {
            'gene': gene,
            'total_mutations': total,
            'in_bounds_count': 0,
            'out_of_bounds_accessions': 'NO_FASTA',
            'best_fasta': '',
            'fasta_length': '',
            'percentage_in_bounds': 0.0,
        }

    best_in_bounds = -1
    best_oob = []
    best_path = ''
    best_len = 0

    for fasta_path in fasta_paths:
        if fasta_path not in seq_len_cache:
            seq_len_cache[fasta_path] = fasta_length(fasta_path)
        seq_len = seq_len_cache[fasta_path]

        in_bounds_count = 0
        oob_accessions = []

        for _, row in gene_df.iterrows():
            try:
                pos_start = int(row['pos_start'])
                pos_end   = int(row['pos_end'])
                if pos_start >= 1 and pos_end <= seq_len:
                    in_bounds_count += 1
                else:
                    oob_accessions.append(str(row['accession']))
            except Exception:
                oob_accessions.append(str(row['accession']))

        if in_bounds_count > best_in_bounds:
            best_in_bounds = in_bounds_count
            best_oob       = oob_accessions
            best_path      = os.path.basename(fasta_path)
            best_len       = seq_len

    pct = round(best_in_bounds / total * 100, 2) if total > 0 else 0.0

    return {
        'gene': gene,
        'total_mutations': total,
        'in_bounds_count': best_in_bounds,
        'out_of_bounds_accessions': ';'.join(sorted(set(best_oob))) if best_oob else '',
        'best_fasta': best_path,
        'fasta_length': best_len,
        'percentage_in_bounds': pct,
    }


def main():
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

    df = pd.read_csv(INPUT_PATH, sep='\t', dtype=str)
    df = df[df['variant_type'] == 'genomic'].copy()
    df['pos_start'] = pd.to_numeric(df['pos_start'], errors='coerce')
    df['pos_end']   = pd.to_numeric(df['pos_end'],   errors='coerce')
    df = df.dropna(subset=['pos_start', 'pos_end'])

    print(f"Source 2: checking {len(df)} mutations across {df['gene'].nunique()} genes")

    seq_len_cache = {}
    results = []

    for gene, gene_df in df.groupby('gene'):
        result = check_gene(gene, gene_df, seq_len_cache)
        results.append(result)
        print(f"  {gene}: {result['in_bounds_count']}/{result['total_mutations']} in-bounds "
              f"({result['percentage_in_bounds']}%)  [{result['best_fasta']}  len={result['fasta_length']}]")

    out_df = pd.DataFrame(results, columns=[
        'gene', 'total_mutations', 'in_bounds_count',
        'out_of_bounds_accessions', 'best_fasta', 'fasta_length', 'percentage_in_bounds'
    ])
    out_df.to_csv(OUTPUT_PATH, sep='\t', index=False)

    total_mut   = out_df['total_mutations'].sum()
    total_match = out_df['in_bounds_count'].sum()
    overall_pct = round(total_match / total_mut * 100, 2) if total_mut > 0 else 0.0

    print(f"\nOverall in-bounds rate: {total_match}/{total_mut} ({overall_pct}%)")
    print(f"Saved to: {OUTPUT_PATH}")


if __name__ == '__main__':
    main()


