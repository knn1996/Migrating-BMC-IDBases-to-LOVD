#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
F="$HERE/fixtures"

rm -rf "$F"
mkdir -p "$F/thesis/idbase"
mkdir -p "$F/inputs"
mkdir -p "$F/genomes/HG16" "$F/genomes/HG17" "$F/genomes/HG18"
mkdir -p "$F/dna/Mane_Select_NG" "$F/dna/Mane_Select_NM" "$F/dna/IDRefseq_NM"
mkdir -p "$F/dna/IDRefseq" "$F/dna/Processed_FASTA_Source1" "$F/dna/Reference_Source2"

: > "$F/inputs/IDBases_Summary.csv"
: > "$F/inputs/alias.csv"
: > "$F/inputs/LRG_RefSeqGene.txt"
printf '>stub\nACGT\n' > "$F/inputs/MANE.GRCh38.refseq_rna.fna"

for b in HG16 HG17 HG18; do
  printf '>chr1\nACGT\n' > "$F/genomes/$b/chr1.fa"
done

for g in RAB27A UNC13D ADA; do
  printf '>%s_NG\nACGT\n' "$g" > "$F/dna/Mane_Select_NG/${g}_NG.fasta"
done

echo "fixtures created under $F"
