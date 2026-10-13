# Release notes — v2.0.0

**Type:** major release. **Outputs differ from v1.1.1.**
The thesis results as defended are reproducible from tag `v1.1.1`; v2.0.0 reflects the
post-thesis correctness work (2026).

v1.1.1 measured *how much* of IDbases the pipeline migrates. v2.0.0 asks *whether what is
migrated is right*: an independent audit traces every published IDbases source unit
(gene, accession, allele) to its fate in the final table, and the work is recorded as a baseline
(Run 0) and 13 fix runs. The pipeline is frozen after Run 13; the remaining open cases
go to a curator review workbook instead of into the code.

## Highlights

- **Independent source audit and run ledger** (`Step10_Audit/`). `audit_against_source.py`
  re-parses the published IDbases pages, sharing no code with the extraction.
  `audit_ledger.py` records each fix as a run against committed outputs, reports gains,
  regressions, GOOD→BAD / BAD→GOOD and pair transitions, and re-scores the whole history
  (`--rescore-all`) whenever the ruler changes, with every existing column required to
  reproduce. History: `Output/Step10_Audit/history/`.
- **Extraction, dedup, BED and IVS fixes** (Runs 1–4): compound-heterozygous alleles, glued
  records and bracket notation in extraction; position-only XM/NM merge and
  `individual_id = GENE_accession` in dedup; delins/dup/minus-strand descriptions in the
  BED → Mutalyzer step; IVS positions no longer inherited from neighbouring features.
- **Track resolution with confidence tiers** (Runs 5–6). One primary record per source
  change (`source_change_id`), decided by protein consistency, then c.-region match, MANE
  Select, gene prior and track priority. Other descriptions of the same change go to
  `alt_descriptions`; every pair carries `confidence_tier` A, B or Review.
- **c./p. coherence** (Run 7). The Track B c. is taken from the same MANE transcript as its
  p.; the audit checks that every p. is the product of its c.'s transcript.
- **Legacy c. numbering remap** (Runs 8–11). A full scan of the source changes found genes
  whose published c. counts from the cDNA record start rather than the start codon. The
  remap is applied per source change from `refs/legacy_remap.tsv` (evidence in
  `refs/legacy_shift_evidence.tsv`, built by `Utilities/build_legacy_remap.py`), only where
  at least two independent changes reconcile with MANE Select. The thesis-era global `ADA`
  offset (−95) is removed; literal records of remapped changes are no longer eligible.
- **MANE Select normalization and `pick_mane` fix** (Run 12). Records on an older version of
  the MANE transcript with an identical CDS are relabelled; records on other transcripts are
  re-queried on MANE Select and accepted only for the same genomic change. `pick_mane` now
  reads Mutalyzer's `tag.details == "MANE Select"`.
- **Review demotions** (Run 13). Pairs the audit cannot vouch for move to tier Review with a
  `review_reason`: genomic positions outside the gene or deep intronic while the source
  change is exonic, records kept in legacy numbering, and source entries describing no
  change. No record content changes.
- **Suspected source errors** (`Utilities/build_source_errors.py`): published notations that
  contradict themselves, listed as candidates for the original curator and reported in separate
  ledger columns without changing any existing measure. The record-level candidates are kept
  privately (`private/`, git-ignored) pending curator review.
- **Canonical export and review workbook** (`Step11_Export/`, standalone offline script):
  `Output/Step11_Export/canonical_export.tsv`, one row per patient–variant pair with allele,
  zygosity, phase, provenance, tier and metadata; and `review_workbook_for_mauno.xlsx`
  (regenerable, not tracked) with the Review-tier pairs and the suspected source errors.
- **Step-based folder reorganization.** The tree now follows the pipeline steps
  (`Step1a_Extraction/` … `Step11_Export/`, mirrored under `Output/`); superseded one-off
  material moved to an untracked archive.

## Before / after (audit ledger, Run 0 → Run 13)

From `Output/Step10_Audit/history/audit_history.tsv`; exact-MANE pairs counted from each
run's committed final table.

| Measure | Run 0 (pre-fix baseline) | Run 13 |
|---|---|---|
| Source units in the final table (of 7,019 in scope) | 5,849 (83.33%) | 6,439 (91.74%) |
| Protein matches source (best record) | 4,660 / 5,133 (90.79%) | 5,471 / 5,558 (98.43%) |
| Pairs consistent with their own source change | 5,002 / 5,747 (87.04%) | 5,541 / 5,601 (98.93%) |
| — excluding suspected source errors | 87.36% | 99.39% |
| Inconsistent pairs | 745 | 60 |
| Excess pairs (individuals) | 653 (652) | 8 (7) |
| c./p. coherence (of known) | 75.34% (4,072 / 5,405) | 100% (5,863 / 5,863) |
| Pairs on the exact MANE Select transcript | 2,777 / 6,589 (42.1%) | 6,518 / 6,562 (99.3%) |
| Final patient–variant pairs | 6,589 | 6,562 |

## Breaking changes

- **Output schema.** The final table (`lovd_flat_with_patients_final.tsv`) gains
  `confidence_tier`, `alt_descriptions` and `review_reason`; `resolution_pairs.tsv` carries
  the per-pair `source_change_ids`. New intermediate files in `Output/Step8_Rescue/`
  (`*.tracks.tsv`, `*.mane.tsv`, `mane_norm_*.tsv`, `review_log.tsv`). The canonical
  deliverable is now `Output/Step11_Export/canonical_export.tsv`.
- **Transcripts and descriptions moved to MANE Select.** Transcript ids, c. and p. differ from
  v1.1.1 wherever records were remapped or normalized; earlier descriptions are kept in
  `alt_descriptions`. Identity is `individual_id = GENE_accession`.
- **New required reference tables** under `Scripts/refs/`: `legacy_remap.tsv`,
  `legacy_shift_evidence.tsv` (and the builder input `legacy_remap_manual.tsv`). The
  source-error candidates are kept privately pending curator review; without them the ledger's
  source-error columns stay empty. The rescue workflow also reads the MANE v1.5 GFF
  (`MANE.GRCh38.v1.5.ensembl_genomic.gff.gz`).
- **Changed counts.** All headline counts change (table above). The thesis figures count
  distinct variants; the audit counts source units (gene, accession, allele). They are not
  interchangeable.

## Validation

From Run 6 on, every run was simulated offline against cached Mutalyzer, VariantValidator
and GenBank responses before the live run, and the live outputs matched the simulation byte
for byte. After each change to the audit itself, the ledger re-scored the full recorded
history, with every pre-existing column required to reproduce exactly. Caches are append-only; uncached queries were forecast and reviewed
before each live run.

## Not in this release

- The LOVD submission file itself: the format is pending with the LOVD curators, and the
  canonical export is designed to map onto it.
- Curator decisions on the Review tier and the suspected source errors.
- Snakemake rule for Step11 (the export is a standalone script for now).
