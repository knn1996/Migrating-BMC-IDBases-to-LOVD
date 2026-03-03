# Migrating-BMC-IDBases-to-LOVD
# Mutation Extraction and Matching

A bioinformatics pipeline for extracting DNA mutations from IDBases 
and cross-checking them against reference mRNA sequences from EMBL/ENA.

## Pipeline Overview
1. `extract_mutations.py` — Parses IDBases HTML files to extract mutation data
2. `download_embl_sequences.py` — Downloads mRNA EMBL flat-files via UniProt + ENA APIs
3. `process_embl_mutations.py` — Converts EMBL → FASTA and validates mutations against sequences
4. `copy_reference_sequences.py` — Collects confirmed reference sequences

## Dependencies
pip install requests pandas openpyxl biopython tqdm

## Usage
Run scripts in order 1 → 4. Edit the CONSTANTS block in each script to set your paths.
