# From Legacy Databases to Modern Standards

**A bioinformatics pipeline for migrating immunodeficiency variant data to LOVD 3.0 with GRCh38 standardisation**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.10+-3776AB.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io/)
[![Reference](https://img.shields.io/badge/Genome-GRCh38-2C6E91.svg)](https://www.ncbi.nlm.nih.gov/grc)
[![Nomenclature](https://img.shields.io/badge/Variants-HGVS-12303D.svg)](https://hgvs-nomenclature.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20886949.svg)](https://doi.org/10.5281/zenodo.20886949)

A reproducible Snakemake pipeline that migrates three decades of curated primary
immunodeficiency variant data from the legacy **IDbases / MUTbase** collection into
**LOVD 3.0**, converting locally-anchored coordinates into GRCh38 genomic positions
and validated HGVS nomenclature.

> **134 genes in scope · 8,790 source variants · 2,240 distinct variants resolved (76.2%)
> across 115 genes, 5,773 patient–variant pairs and 2,729 patients.**

Post-thesis (2026, Run 13): 6,439 of 7,019 in-scope source units in the final table (91.74%), and
98.93% of patient–variant pairs consistent with the source. See
[Post-thesis correctness work](#post-thesis-correctness-work-2026).

📄 **Thesis (open access):** [LUP record 9242187](https://lup.lub.lu.se/student-papers/record/9242187) ·
[full text PDF](https://lup.lub.lu.se/student-papers/record/9242187/file/9242188.pdf) ·
[popular-science summary](https://lup.lub.lu.se/student-papers/record/9242187/file/9242190.pdf)
Master's thesis KMBM01, Biotechnology (MSc), Lund University, 2026.
Supervisor: Prof. Mauno Vihinen. Examiner: Ghjuvan Micaelu Grimaud.

---

## Why this was hard: the nested-reference problem

IDbases variant coordinates are anchored to **IDRefSeq** sequences — partial ENA clone
fragments chosen as local references in the late 1990s and early 2000s, before
genome-assembly-anchored annotation became standard. These fragments are not independent
references: each is **nested inside** the modern full-length LRG / `NG_` RefSeqGene record
for the same gene, at a gene-specific and previously undocumented offset.

Naïvely treating an IDRefSeq coordinate as an `NG_` coordinate therefore fails **silently
and systematically** — the position is wrong by a constant that differs per gene.

The **empirical offset algorithm** recovers that constant by BLAST-aligning each IDRefSeq
fragment to its LRG / `NG_` parent, then confirming the recovered offset against
substitution variants: a substitution asserts a reference base at a position, so a correct
offset reproduces that base and an incorrect one does not.

![Dot plot of anchor REF-base match rate for eight representative genes with 95% Wilson confidence intervals. AIRE, CTSC, PRF1, ITGB2, CYBB and CD40L all sit at or near 100% and are accepted into Track A; RAG1 (89%) and RAG2 (79%) fall below the 90% threshold and are rejected to Track B, each with a significant binomial p-value.](docs/figures/fig01_anchor_match_rates_slide.png)

*Per-gene anchor confirmation. Across all **94 testable genes**, **92** reach a ≥90% anchor
match rate, with **99.29% aggregate anchor concordance** (3,799 / 3,826 anchors). Only
`RAG1` and `RAG2` fail, and they fall through to the transcript-based track rather than
being silently mis-placed.*

**Anchor-transferred validation.** This confirmation is only directly available for
substitutions. Insertions and duplications carry no REF-base assertion and pass HGVS
validation *silently* even under a coordinate frame mismatch. Because all variant classes
within a gene share one offset, substitution-anchored confirmation transfers to the
length-changing classes. This validates the **positional placement** of insertions and
duplications — it does **not** independently validate inserted sequence content.

---

## Results (as defended)

The figures in this section are the thesis results as defended and are left unchanged.

### Read the hierarchy, not a single number

| Level | Count |
|---|---|
| Source patient records (IDbases) | 6,259 |
| Source variants | 8,790 (15 homozygous / 5,236 heterozygous / 789 hemizygous) |
| Resolved patient–variant records (row-level) | 7,776 |
| Patient–variant pairs | 5,773 |
| **Distinct variants** | **2,240** |
| Genes with ≥1 resolved variant | 115 of 134 (85.8%) |
| Patients linked | 2,729 |
| Automated resolution rate | **76.2%** |
| Addressable ceiling (incl. rescuable failures) | ~94% |

The merged table holds 7,776 resolved *records*; these collapse to 5,773
*patient–variant pairs* once track duplication is removed; those represent 2,240
*distinct variants* in 2,729 patients. Quoting the record count as a variant count
overstates the dataset; quoting only the distinct count understates its clinical reach.
The row-level figures are merge artifacts and are **not** headline counts.

### Three-track resolution

Variants are resolved in parallel down three independent annotation routes, then merged
and deduplicated:

| Track | Reference route | Resolution rate | Representative for |
|---|---|---|---|
| **A** | Genomic `NG_` (offset-corrected IDRefSeq) | 99.3% (3,582 / 3,606) | 951 variants |
| **B** | MANE Select `NM_` | 94.1% (4,496 / 4,775) | 1,187 variants |
| **C** | Legacy `NM_` (IDRefSeq-era transcript) | 59.5% (2,862 / 4,811) | 102 variants |

![Stacked bar chart of the 2,240 distinct variants by representative resolution track and variant class. Track A contributes 951 variants (579 substitutions, 268 deletions, 61 duplications, 43 insertions), Track B contributes 1,187 (845 substitutions, 236 deletions, 58 duplications, 48 insertions) and Track C contributes 102 (75 substitutions plus a small deletion remainder).](docs/figures/fig03_class_by_track.png)

The tracks are **complementary, not nested**: 863 variants (38.5%) were independently
resolved by exactly two of the three, and 1,377 by a single track. Track A carries
proportionally more of the rarer, length-changing and legacy-annotated variants that the
MANE route cannot resolve natively.

### Unresolved variants

698 distinct variants (23.8%) did not resolve automatically.

![Horizontal bar chart of the 698 unresolved distinct variants by Mutalyzer error code, coloured by disposition. ENOSELECTORFOUND dominates with 512 variants and is automatically rescuable; ESEQUENCEMISMATCH accounts for 171 and is classified as an unrescuable data error; EINTRONIC (12), ESYNTAXUEOF, ERANGEREVERSED and EINSERTIONRANGE (1 each) make up the remainder.](docs/figures/fig08_failure_taxonomy.png)

The dominant failure mode is `ENOSELECTORFOUND` (512 variants, **73.4%** of failures),
caused by RefSeqGene `NG_` records pinning transcript versions incompatible with current
MANE releases. That is a mechanical, recoverable failure rather than a data-quality one —
**75.8%** of unresolved variants are classified rescuable (525 automated, 4 manual),
which is why the addressable ceiling sits near 94%.

---

## Post-thesis correctness work (2026)

The thesis results above are reported as defended. After the defence, the work moved from *how much* of IDbases can be migrated to *whether what is migrated is right*. A LOVD submission is only useful if every record is correct, so precision now takes priority over recovery.

**Counting definitions.** The thesis figures above count *distinct variants* (and patient–variant pairs). The audit below counts *source units*: one (gene, accession, allele) entry of the published IDbases pages. The two are not interchangeable.

### An independent audit, not the pipeline's own validation

Mutalyzer and VariantValidator check that an HGVS description is valid on its reference. They cannot tell whether it describes the variant the original authors reported. To answer that, a separate audit (`Step10_Audit/audit_against_source.py`) re-parses the published IDbases pages (`pub.html`). It shares no code with the extraction step, and it traces every source unit to its fate in the final table.

- **Fixed ruler.** The key is (gene, accession, allele) from the source, never the remapped HGVS, which can change between runs. The audit covers 7,047 source units, 7,019 of them in scope.
- **Three axes, always reported together:**
  - Recovery: share of source units present in the final table.
  - Best-record accuracy: the protein consequence matches the source p. description.
  - Precision: each patient–variant pair is consistent with its own source change, and no individual has excess pairs.
- **c./p. coherence.** The p. description must be the product of the transcript used for the c. description.
- **Run ledger** (`Step10_Audit/audit_ledger.py`). Each fix is recorded as one run against committed outputs. A run lists gains, regressions, GOOD→BAD and BAD→GOOD transitions, and pair transitions against the previous run. When the audit itself changes, the whole history is re-scored (`--rescore-all`), and every existing column must reproduce exactly.

### Coverage gains watched against consistency

A higher mapping rate is only accepted if consistency holds. The first fix shows why this matters:

| Run | Fix | In final | Protein good | Pair consistency | Excess pairs |
|---|---|---|---|---|---|
| 0 | Pre-fix baseline | 83.33% | 90.79% | 87.04% | 653 |
| 1 | Extraction fix | 89.84% | 92.74% | 85.50% | 1,067 |

In Run 1, recovery rose by 6.5 points, but pair consistency fell and excess pairs grew from 653 to 1,067. There were 10 regressions and 18 GOOD→BAD transitions. Looking at recovery alone would have hidden this. Most of the later fixes went into precision.

### Biological inconsistencies between source and output protein

For each unit, the audit compares the protein consequence in the final record with the source p. description. It classifies every mismatch: different change, shifted position, frameshift start or stop differences, or type difference at the same position. Tracing these mismatches back to their cause found most of the bugs listed below.

### Mixed historical conventions for c. positions

IDbases stores c. positions as published. Many pre-HGVS papers counted from the first base of their own cDNA record instead of the A of the start codon. The coding sequence is the same; only the origin differs, by the length of that record's 5′ UTR. A full scan of 7,920 source changes in 133 genes found 28 genes with a verified constant shift, such as RAG1 −112, AICDA −76 or AIRE −127. Eighteen of these genes mix both conventions in one database. The style of the p. description usually tells them apart: one-letter (old) versus three-letter (HGVS).

A shift is only accepted when at least two independent changes reconcile with the MANE Select sequence. The remap is applied per source change from a committed table (`refs/legacy_remap.tsv`), so a change that already follows HGVS numbering is never shifted. The thesis-era global ADA offset of −95 was replaced in this way, because it also shifted ADA changes that were already correct.

Extending the remap to every verified gene (Run 10) recovered exactly the 110 source units the scan had predicted, with no regressions and no new inconsistent pairs. Seventeen existing pairs moved from wrong keys to the correct coding change; for example, 7 RAG1 pairs had been keyed on an intronic position and are now keyed on c.256_257del. Three RFXAP pairs that VariantValidator had placed at the legacy position went from a mismatching protein to a match.

The remaining mismatches are increasingly errors in the source rather than in the pipeline. For example, AIRE A0082 reports p.S64del for a 36-nt deletion, and CD40L C0047 gives a c. range that crosses an exon–intron boundary. The audit now flags 24 suspected source errors (29 pairs, pending review by the original curator) in separate columns, without changing any existing measure. Excluding them, pair consistency in Run 13 is 99.39%.

### Fixes, measured on both right and wrong records

Each fix is checked against units that were already correct (regressions, GOOD→BAD) and against those it was meant to repair (gains, BAD→GOOD). Before any live run, it is simulated offline against cached results, and the real outputs must then match the simulation byte for byte.

| Run | Fix | In final | Protein good | Pair consistency | Excess | c./p. coherence |
|---|---|---|---|---|---|---|
| 0 | Pre-fix baseline | 83.33% | 90.79% | 87.04% | 653 | 75.34% |
| 1 | Extraction: compound heterozygous alleles, glued records, brackets | 89.84% | 92.74% | 85.50% | 1,067 | 74.55% |
| 2 | Dedup: position-only XM/NM merge; `individual_id = GENE_accession` | 90.03% | 93.06% | 86.51% | 987 | 74.14% |
| 3 | BED → Mutalyzer: delins, dup, minus strand | 89.98% | 93.11% | 89.93% | 733 | 73.18% |
| 4 | IVS positions inherited from neighbours | 89.98% | 93.11% | 89.93% | 731 | 73.17% |
| 5 | Track resolution: one primary record per source change | 89.98% | 96.20% | 96.41% | 23 | 73.69% |
| 6 | Dedup transcript split, stage-2 links, accession fix | 89.98% | 96.20% | 96.65% | 8 | 73.69% |
| 7 | Track B c. taken from the same MANE transcript as p. | 89.98% | 96.20% | 96.65% | 8 | 100% |
| 8 | Legacy c. numbering remap | 89.98% | 98.19% | 98.66% | 8 | 100% |
| 9 | ADA per-change remap | 90.14% | 98.33% | 98.80% | 8 | 100% |
| 10 | Legacy remap extended to the full numbering scan | 91.71% | 98.42% | 98.88% | 8 | 100% |
| 11 | Old-numbering records of remapped changes not eligible | 91.74% | 98.42% | 98.91% | 8 | 100% |
| 12 | MANE Select normalization; transcript-picker fix | 91.74% | 98.43% | 98.93% | 8 | 100% |
| 13 | Known-wrong and unverifiable records demoted to Review | 91.74% | 98.43% | 98.93% | 8 | 100% |

From Run 0 to Run 13:
- Recovery rose from 5,849 to 6,439 source units (83.33% → 91.74%).
- Inconsistent pairs fell from 745 to 60.
- 6,518 of 6,562 final pairs (99.3%) are on the exact MANE Select transcript; the rest are flagged records or have no MANE Select (IGHG2).
- Excess pairs fell from 653 to 8.
- Protein BAD→GOOD totalled 460 units, with all GOOD→BAD (18) in Run 1.
- c./p. coherence rose from 75.34% to 100% of checkable pairs.

### Knowing when to stop

Each run fixed less than the one before (from Run 8 onwards, about 0.1 point of consistency per run), and later candidates increasingly needed rules tailored to single records. To avoid fitting the pipeline to its own audit, a stopping rule was set:
- A run must fix a whole class of records, covering at least 0.5% of in-scope source units (35), or remove records known to be wrong.
- Smaller or one-off cases go to a Review workbook for the curator instead of into the code.

The pipeline was frozen after Run 13. Run 13 changes no record content. It moves records the audit cannot vouch for to the Review tier, with a stated reason: genomic positions that land deep in an intron or outside the gene while the source describes an exonic change (41 pairs), records still in legacy numbering (6), and source entries that describe no change (4).

### Traceability and resolution policy

- `source_change_id` links every record to the source change it came from.
- `individual_id = GENE_accession`, because accessions are only unique within one gene database.
- When the tracks disagree, one primary record is kept per source change. Other descriptions of the same genomic change go to `alt_descriptions`, and every final record carries a `confidence_tier` (A, B or Review).
- In-cis alleles (`c.[A;B]`) legitimately give several records per source unit and are never collapsed.

### Reproducibility discipline

- Mutalyzer, VariantValidator and GenBank responses are cached append-only.
- Before each live run, the uncached queries are forecast and reviewed.
- Each fix is committed as code, then primary outputs, then rescue outputs, then the ledger entry, so any recorded run can be re-scored from git.

---

## Pipeline

![Three-panel pipeline diagram. Steps 1 to 5, Reference Harmonisation: extract legacy variants from 136 IDbases (~8,790 variants), confirm reference sequences by resolving IDRefSeq clones within LRG/NG sequences using the offset algorithm, and convert coordinates via BLAST alignment, BED mapping and LiftOver to GRCh38, outputting variants on GRCh38 with strand information. Step 6, HGVS Normalisation with Mutalyzer: three complementary annotation strategies run in parallel — Track A normalises at genomic level using NG_ references (position validated on GRCh38, independent of transcript version), Track B at coding level using MANE Select transcripts (current clinical standard), Track C at coding level using legacy NM_ transcripts (backward compatibility and traceability to original records); Mutalyzer 3 validates and returns HGVS at g., c., r. and p. levels. Step 7, Priority Merge and LOVD 3.0 Import: representative records are chosen in priority order MANE (Track B), Genomic (Track A), Legacy (Track C), and all variants are imported into LOVD 3.0 as HGVS-compliant entries.](docs/figures/pipeline_overview.png)

*Steps 1–5 harmonise references and coordinates, Step 6 runs the three annotation tracks
through Mutalyzer 3, Step 7 merges and deduplicates into the variant table, and
`rescue.smk` (Step 8) recovers a further tranche of failures. Outcome: 136 IDbases
processed, 2,240 distinct variants migrated across **5,773 patient–variant pairs** in
**2,729 patients**.*

### Two-stage deduplication

Row counts overstate biological content, so deduplication runs in two stages:

1. **Version-drift collapse** — strip the transcript accession from the normalised coding
   key, collapsing records that differ only by transcript version (`NM_….2` vs `.3`).
   7,776 input rows → 2,661.
2. **XM/NM resolution** — per gene + chromosome + genomic position, union patient lists
   into the curated `NM_` record and drop the predicted `XM_` duplicate. This removes 421
   rows, of which `CYBB` accounts for 417. Final: **2,240** distinct variants.

---

## Repository structure

```
.
├── README.md
├── CITATION.cff
├── LICENSE
├── docs/figures/               # README figures
├── 01_Reference_Genomes/
├── 02_Source_Database/         # IDbases gene / UniProt summary spreadsheets
├── 03_BED_Files/               # gene locus BED files (hg16–hg18, hg38)
├── 04_Mutation_Processing/     # core pipeline
│   ├── Scripts/
│   │   ├── Snakefile           # primary workflow (Steps 1–7)
│   │   ├── rescue.smk          # rescue workflow — reads primary outputs as fixed inputs
│   │   ├── config/             # config.example.yaml (template) + config.yaml (gitignored)
│   │   ├── envs/               # per-rule conda environments
│   │   ├── R for figures/      # ggplot2 scripts for every figure
│   │   ├── refs/               # VariO ontology, legacy_remap.tsv, legacy_shift_evidence.tsv, source_errors.tsv
│   │   ├── Step1a_Extraction/  Step1b_RSG_Mapping/  Step1c_SeqDownload/
│   │   ├── Step2_RefCheck/     Step3_BLAST/  Step4_BED/  Step5_LiftOver/
│   │   ├── Step6_Mutalyzer/    Step7_Merging/  Step8_Rescue/
│   │   ├── Step9_NMD/          Step10_Audit/  Step11_Export/
│   │   ├── Utilities/          # one-off checks and maintenance scripts
│   │   └── tests/
│   ├── Output/                 # mirrors the Step* layout above, plus cache/ and Step10_Audit/history/
│   └── Logs/                   # run logs (gitignored)
└── 05_Testing/                 # validation and test scripts
```

The two workflows are deliberately separate. `Snakefile` runs Steps 1–7. `rescue.smk` declares
the primary outputs as **external inputs with no producing rules**, so it never regenerates or
overwrites them, and writes its results from `Output/Step8_Rescue/` onwards.

The as-defended thesis results (see [Results (as defended)](#results-as-defended)) are
reproducible from tag [`v1.1.1`](https://github.com/knn1996/Migrating-BMC-IDBases-to-LOVD/tree/v1.1.1),
which predates the rescue workflow. From v2.0.0 on, the outputs reflect the
[post-thesis correctness work](#post-thesis-correctness-work-2026).

---

## Reproducing

The pipeline is a config-driven [Snakemake](https://snakemake.readthedocs.io/) workflow.
You need `snakemake` and `conda`/`mamba`; every rule pulls its own dependencies from
`envs/*.yaml` via `--use-conda`, so no manual `pip install` is required. Network access to
Mutalyzer 3, VariantValidator and NCBI E-utilities is required.

Every relative path in the config resolves against the `Scripts` directory, so always
launch Snakemake from there:

```bash
cd 04_Mutation_Processing/Scripts

# 1. Create your local config from the committed template, then edit the paths
cp config/config.example.yaml config/config.yaml

# 2. Dry run — resolves the DAG without executing anything
snakemake --use-conda --cores 8 -n

# 3. Primary workflow
snakemake --use-conda --cores 8

# 4. Rescue workflow — reads primary outputs, writes to Output/Step8_Rescue/,
#    then NMD annotation (Output/Step9_NMD/) and the source audit (Output/Step10_Audit/)
snakemake -s rescue.smk --use-conda --cores 8
```

**Outputs.** The final deliverable is the canonical export
(`Step11_Export/canonical_export.tsv` under the output root); the primary workflow's LOVD 3.0
review workbook is `lovd_review.xlsx`. Intermediates land under `results/`: `results/seq` (downloaded sequences),
`results/blast` (BLAST DBs and hits), `results/bed` (`raw` → `filtered` → `hg38`),
`results/cache` (Mutalyzer response cache) and `results/logs`.

`Step11_Export/build_canonical_export.py` is a standalone offline script, not a Snakemake rule; run it
from `Scripts` after the rescue workflow.

Step 3 (BLAST) was executed on the LUNARC COSMOS cluster; the rule is written for SLURM
submission but falls back to local execution. Mutalyzer results are cached by exact HGVS
string — **do not clear the cache**, as roughly 5,000 valid entries are reusable and
re-querying is slow.

### Choosing the gene set (`run_set`)

**The default config runs a 3-gene smoke test** (`RAB27A`, `UNC13D`, `ADA`) so you can
validate the install quickly. To reproduce the full thesis dataset, set:

```yaml
run_set: "all"
```

in `config/config.yaml`. With `run_set: all`, the gene list is auto-discovered from the
IDbases source directories (`idbase_subdir`).

### Required external data

These inputs are third-party or large and are **not** stored in git. Download or deposit
them and point the matching config key at their location:

| Input | Source | Config key |
| --- | --- | --- |
| IDbases source directories (per-gene MUTbase records) | [BMC IDbases](http://structure.bmc.lu.se/idbase/) | `idbase_subdir` |
| hg16 / hg17 / hg18 chromosome FASTA (`chr*.fa`) | [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html) | `genome_dir` + `build_dirs` |
| MANE RNA FASTA (`MANE.GRCh38.*.refseq_rna.fna`) | [NCBI MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) | `inputs.mane_rna_fasta` |
| MANE GFF (`MANE.GRCh38.v1.5.ensembl_genomic.gff.gz`), read by the rescue workflow | [NCBI MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/) | `thesis_dir`/`01_Reference_Genomes/GRCh38_MANE/` |
| `LRG_RefSeqGene.txt` | [NCBI RefSeqGene / LRG](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/) | `inputs.lrg_refseqgene_txt` |
| Gene alias table (`alias.csv`) | Project-provided | `inputs.alias_csv` |
| Curated IDRefSeq sequences (NG / NM / processed FASTA) | Zenodo / project deposit | `ref_seq.*`, `dna_seq.*` |

---

## Key outputs

| File | Contents |
|---|---|
| `Output/Step1a_Extraction/all_mutations.tsv` | Extracted patient-level variant records, source notation preserved |
| `Output/Step2_RefCheck/lrg_offset_results.csv` | Per-gene IDRefSeq → LRG / `NG_` offset, alignment span, anchor match rate |
| `Output/Step2_RefCheck/gene_confidence_NG_IDRefseq.csv` | Per-gene anchor statistics, Wilson CI, acceptance status (optional rule) |
| `Output/Step7_Merging/merged_variants.tsv` | Three-track merged, resolved records (row-level) |
| `Output/Step7_Merging/dedup_merged_variants.tsv` | Primary deduplicated variant set (fixed input to the rescue workflow) |
| `Output/Step7_Merging/lovd_review.xlsx` | LOVD 3.0 review workbook |
| `Output/Step8_Rescue/unresolved_disposition.tsv` | Failures with error code, subcategory and rescuability class |
| `Output/Step8_Rescue/lovd_flat_final.tsv` | Rescue-inclusive LOVD flat table |
| `Output/Step8_Rescue/lovd_flat_with_patients_final.tsv` | Rescue-inclusive table with VariO terms and patient metadata |
| `Output/Step8_Rescue/lovd_review_final.xlsx` | Rescue-inclusive LOVD 3.0 review workbook |
| `Output/Step8_Rescue/rescue_metrics_comparison.csv` | Primary vs rescue-inclusive resolution metrics |
| `Output/Step9_NMD/nmd_annotation.tsv` | Premature-stop variants with NMD prediction and proposed `r.`/`p.` descriptions |
| `Output/Step10_Audit/audit_summary.tsv` | Source-to-output audit: fate of every IDbases record, metadata and protein checks |
| `Output/Step10_Audit/history/` | Audit ledger: one row per recorded run (`audit_history.tsv`), per-run snapshots and change lists |
| `Output/Step11_Export/canonical_export.tsv` | Canonical export: one row per patient–variant pair with allele, phase, provenance, confidence tier and metadata |
| `Output/Step11_Export/review_workbook_for_mauno.xlsx` | Curator review workbook: Review-tier pairs and suspected source errors (regenerable, untracked) |
| `Scripts/refs/legacy_remap.tsv` | Per-source-change legacy c. numbering remap, with its evidence in `legacy_shift_evidence.tsv` |
| `Scripts/refs/source_errors.tsv` | Suspected IDbases source errors (candidates pending curator review) |

## Limitations

- Mutalyzer validates position *and* REF base for substitutions, but position range only
  for deletions, insertions, duplications and indels. Positional confidence for
  length-changing classes is inherited through anchor transfer, not asserted directly.
- Several source databases are structurally non-standard and are excluded or handled
  specially: `BLNK` (single unparseable IVS entry), `LRRC8A` (chromosomal deletion only),
  `TAPBP` (curation error), `SH2` (a domain subset of `BTK`), and `IGHG2` / `IGHM` (no LRG
  records exist, owing to somatic recombination at these loci — handled as `NG_`
  RefSeqGene directly). `BTK` and `SH2` are curated separately by Prof. Vihinen.
- Legacy c. numbering is handled per source change from `refs/legacy_remap.tsv`: 26 genes
  remapped, out of the 28 with a verified shift in the full numbering scan (CYBB and NCF1 were
  left out); the old global `ADA` offset of −95 was removed.
- Phenotype and disease free-text fields are sparsely and inconsistently populated in the
  source (20–58% coverage, versus 87–99% for variant-level fields) and are migrated as-is.
- 286 remaining genuine errors still need characterisation by error code and gene.

## Data sources

- [IDbases](http://structure.bmc.lu.se/idbase/) — immunodeficiency mutation databases
- [LOVD](https://www.lovd.nl/) — Leiden Open Variation Database (migration target)
- [UniProt REST API](https://rest.uniprot.org/) — protein / gene cross-references
- [ENA](https://www.ebi.ac.uk/ena/) — EMBL nucleotide flat-files
- [Mutalyzer](https://mutalyzer.nl/) and [VariantValidator](https://variantvalidator.org/) — HGVS validation
- [UCSC Genome Browser](https://genome.ucsc.edu/) — reference genome assemblies

## Citation

If you use this pipeline or the migrated dataset, please cite the thesis:

> Nguyen, Ngoc Khoi (2026). *From Legacy Databases to Modern Standards: A Bioinformatics
> Pipeline for Migrating Immunodeficiency Variant Data to LOVD 3.0 with GRCh38
> Standardisation.* Master's thesis, Lund University.
> <https://lup.lub.lu.se/student-papers/record/9242187>

and the software:

> Nguyen, Ngoc Khoi (2026). *Migrating BMC IDbases to LOVD.* Zenodo.
> <https://doi.org/10.5281/zenodo.20886949>

A machine-readable `CITATION.cff` is provided; GitHub renders a **Cite this repository**
button from it.

## Acknowledgements

Supervised by **Prof. Mauno Vihinen** (Lund University), who developed the original
MUTbase / IDbases system. Examined by **Ghjuvan Micaelu Grimaud**. Thanks to the IDbases
curators for three decades of work, and to the LOVD, ENA, UniProt, Mutalyzer and
VariantValidator teams for open infrastructure. Computation performed on the LUNARC COSMOS
cluster (project `lu2025-2-88`).

## License

Released under the [MIT License](LICENSE).
