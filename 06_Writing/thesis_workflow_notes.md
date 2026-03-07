# Thesis Project Workflow Notes
*Last updated: 2026-03-04*

---

## Overview

The original reference sequences were lost during maintenance and had to be recovered from two sources.

---

## Step 1: Reference Sequence Recovery

### Source 1 — DNA from `idbases` HTML files
- Extracted DNA sequences from `*dna.html` files for each ID in the `idbases` folder
- Raw sequences saved to:
  `C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\Raw file`
- Processed to remove:
  - Tabs at start of lines
  - Numbering at end (or start) of lines
- Cleaned sequences saved to:
  `... \DNA sequences\Processed FASTA file\`
- **~60 sequences** recovered this way

### Source 2 — mRNA from ENA via UniProt IDs
- Scraped UniProt IDs from `*pub.html` files in `idbases`
- Downloaded mRNA sequences from **ENA (European Nucleotide Archive)**
- Confirmed via **Mutation Extraction and Matching** folder:
  - Extracted all **deletions** and **substitutions** from idbases
  - Retrieved original mutated positions and nucleotides
  - Verified RNA mutations matched the reference sequences

### Source 3 — Refseq project from National Library of Medicine, NCBI
- List of genes in this project was downloaded (LRG_RefseqGene.txt)



---

## Step 2: BLAST (run on WSL / Linux environment)

- Ran `blast_all.sh` across all sequences
- Some runs failed with **OOM (Out of Memory)** errors → moved to `BLAST fail due to memory/` folder
- OOM jobs resubmitted on **LUNARC** (HPC cluster)
- Results stored in: `BLAST results/`

---

## Step 3: Summarize BLAST Results

- Ran `summarize_blast_first_hits.sh` to extract first hits from BLAST results
- Selected hits with **100% match**

---

## Step 4: Generate BED Files

- Script: `generate_bed.py`
- Input: 100% match BLAST hits
- Output: BED files

### Approach by source

**Sources 1 & 2 — All three genome builds (cross-check mode):**
- BED files generated for **all three** BLAST targets: hg16, hg17, and hg18
- This was intentional — to cross-check whether lifting over from different source assemblies yields different results (e.g. coordinate loss or out-of-bounds after liftover)
- Results in up to 3 BED files per gene (e.g. `ADA_vs_hg16.bed`, `ADA_vs_hg17.bed`, `ADA_vs_hg18.bed`) and correspondingly up to 3 lifted-over hg38 outputs

**Source 3 — Most recent build only (optimised mode):**
- Since Sources 1 & 2 confirmed that liftover results are consistent across builds, Source 3 was processed more efficiently
- For each gene, only the **first 100% match hit from the most recent available build** was used (priority: hg18 > hg17 > hg16)
- Only that single best hit was lifted over to hg38

---

## Step 5: LiftOver to hg38

- Ran **LiftOver** on BED files
- Mapped coordinates to most recent human genome build: **hg38 (GRCh38)**
- Sources 1 & 2: up to 3 liftover outputs per gene (one per source assembly)
- Source 3: single liftover output per gene (most recent 100% hit only)

---

## Step 6: BED → VCF Conversion

- Script: `bed_to_vcf.py`
- Converted lifted-over BED files to VCF format

---

## Step 7: Variant Validation *(Next Step)*

- Script: `validate_variants.py`
- Calls the **validate_variants API**
- **⬅ Paused here — resume from this step**

---

## Notes
- WSL used for local Linux environment
- LUNARC used for jobs that exceeded local memory limits
