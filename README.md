# Immunodeficiency Disease Mutation Pipeline

A bioinformatics pipeline for extracting, processing, and mapping DNA mutations
from the IDBases immunodeficiency mutation databases to human reference genomes.

## Project Structure

| Folder | Contents |
|--------|----------|
| `02_Source_Database/` | IDBases gene/UniProt summary spreadsheets |
| `03_BED_Files/` | Gene locus BED files for hg16–hg18 and hg38 |
| `04_Mutation_Processing/Scripts/` | Core pipeline scripts |
| `04_Mutation_Processing/Mutation Extraction and Matching/` | EMBL download & mutation matching scripts |
| `04_Mutation_Processing/Output/` | Pipeline outputs (TSV, VCF) |

## Pipeline Overview

### Stage 1 – Extract Mutations
`extract_mutations.py` — Parses IDBases pub.html files and extracts all DNA
mutations into `all_mutations.tsv` (gene, accession, position, type).

### Stage 2 – Generate BED Files
`generate_bed.py` / `generate_bed.ps1` — Generates BED coordinate files for
each gene across genome assemblies hg16, hg17, hg18.

### Stage 3 – Download Reference Sequences
`download_embl_sequences.py` — Queries UniProt REST API for mRNA EMBL
accessions and downloads EMBL flat-files from ENA for each gene.

### Stage 4 – Match Mutations to Sequences
`process_embl_mutations.py` — Converts EMBL → FASTA, extracts CDS offsets,
and cross-checks each mutation against downloaded sequences.

### Stage 5 – Copy Reference Sequences
`copy_reference_sequences.py` — Collects confirmed reference FASTA sequences
(genes with 1–3 perfect mutation matches) into a single folder.

### Stage 6 – Convert to VCF
`bed_to_vcf.py` — Converts BED mutation coordinates to hg38 VCF format.

### Stage 7 – LRG Mapping
`LRG.py` — Maps mutations to Locus Reference Genomic (LRG) records.

## Dependencies
```bash
pip install requests pandas openpyxl biopython tqdm
```

## Usage

1. Place `IDBases_Summary_with_UniProt.xlsx` in the working directory
2. Run scripts in order (Stage 1 → 7)
3. Edit the `CONSTANTS` block at the top of each script to set your local paths

## Data Sources

- [IDBases](http://structure.bmc.lu.se/idbase/) — Immunodeficiency mutation databases
- [UniProt REST API](https://rest.uniprot.org/) — Protein/gene cross-references
- [ENA](https://www.ebi.ac.uk/ena/) — EMBL nucleotide flat-files
- [UCSC Genome Browser](https://genome.ucsc.edu/) — Reference genome assemblies
