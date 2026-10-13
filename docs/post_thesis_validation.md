# Post-thesis validation (2026)

The thesis results in the [README](../README.md#results-as-defended) are reported as defended. After the
defence, the work moved from *how much* of IDbases can be migrated to *whether what is migrated is
right*: a LOVD submission is only useful if every record is correct, so precision takes priority over
recovery. How the audit measures this is described in [audit_methodology.md](audit_methodology.md).

## Coverage gains watched against consistency

A higher mapping rate is only accepted if consistency holds. The first fix shows why this matters:

| Run | Fix | In final | Protein good | Pair consistency | Excess pairs |
|---|---|---|---|---|---|
| 0 | Pre-fix baseline | 83.33% | 90.79% | 87.04% | 653 |
| 1 | Extraction fix | 89.84% | 92.74% | 85.50% | 1,067 |

In Run 1, recovery rose by 6.5 points, but pair consistency fell and excess pairs grew from 653 to 1,067. There were 10 regressions and 18 GOOD→BAD transitions. Looking at recovery alone would have hidden this. Most of the later fixes went into precision.

## Biological inconsistencies between source and output protein

For each unit, the audit compares the protein consequence in the final record with the source p. description. It classifies every mismatch: different change, shifted position, frameshift start or stop differences, or type difference at the same position. Tracing these mismatches back to their cause found most of the bugs listed below.

## Mixed historical conventions for c. positions

IDbases stores c. positions as published. Many pre-HGVS papers counted from the first base of their own cDNA record instead of the A of the start codon. The coding sequence is the same; only the origin differs, by the length of that record's 5′ UTR. A full scan of 7,920 source changes in 133 genes found 28 genes with a verified constant shift, such as RAG1 −112, AICDA −76 or AIRE −127. Eighteen of these genes mix both conventions in one database. The style of the p. description usually tells them apart: one-letter (old) versus three-letter (HGVS).

A shift is only accepted when at least two independent changes reconcile with the MANE Select sequence. The remap is applied per source change from a committed table (`refs/legacy_remap.tsv`), so a change that already follows HGVS numbering is never shifted. The thesis-era global ADA offset of −95 was replaced in this way, because it also shifted ADA changes that were already correct.

Extending the remap to every verified gene (Run 10) recovered exactly the 110 source units the scan had predicted, with no regressions and no new inconsistent pairs. Seventeen existing pairs moved from wrong keys to the correct coding change; for example, 7 RAG1 pairs had been keyed on an intronic position and are now keyed on c.256_257del. Three RFXAP pairs that VariantValidator had placed at the legacy position went from a mismatching protein to a match.

The remaining mismatches are increasingly errors in the source rather than in the pipeline. The audit now flags 24 candidate source changes (29 pairs) whose own c. and p., or g. and c., contradict each other; they are pending review by the original curator and are reported in separate columns, without changing any existing measure (criteria and categories in [audit_methodology.md](audit_methodology.md#suspected-source-errors)). Excluding them, pair consistency in Run 13 is 99.39%.

## Fixes, measured on both right and wrong records

The work is recorded as a baseline (Run 0) and 13 fix runs. Each fix is checked against units that were already correct (regressions, GOOD→BAD) and against those it was meant to repair (gains, BAD→GOOD). Before any live run, it is simulated offline against cached results, and the real outputs must then match the simulation byte for byte.

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

## Knowing when to stop

Each run fixed less than the one before (each measure reached its plateau at a different run — coherence at Run 7, consistency after the legacy remap in Run 8, recovery after Run 10 — and from Run 11 on no measure moved by more than 0.03 points), and later candidates increasingly needed rules tailored to single records. To avoid fitting the pipeline to its own audit, a stopping rule was set:
- A run must fix a whole class of records, covering at least 0.5% of in-scope source units (35), or remove records known to be wrong.
- Smaller or one-off cases go to a Review workbook for the curator instead of into the code.

![Heatmap of each fix run's share of the total Run 0 to Run 13 improvement, for five measures (rows) and Runs 1 to 13 (columns). Gains are blue and losses orange. Each measure improves in a few runs: recovery mostly in Run 1 (77%) and Run 10 (19%); protein matches in Runs 1, 5 and 8 (26%, 40%, 26%); pair consistency in Runs 3, 5 and 8 (29%, 54%, 17%) after a loss in Run 1 (−13%); excess pairs in Run 5 (110%) after Run 1 added excess (−64%); c./p. coherence in Run 7 (107%). Runs 11 to 13 are outlined as a plateau where every share is below 1%.](figures/fig09_gain_heatmap.png)

*Where each measure's improvement came from: every measure has its own decisive runs, and from Run 11
on none moves.*

The pipeline was frozen after Run 13. Run 13 changes no record content. It moves records the audit cannot vouch for to the Review tier, with a stated reason: genomic positions that land deep in an intron or outside the gene while the source describes an exonic change (41 pairs), records still in legacy numbering (6), and source entries that describe no change (4). The pipeline stays frozen: a new fix run is made only if a major problem is found that meets the stopping rule, for example from the curator's review. Curator-given annotation rules, such as the NMD rule below, are applied outside the audit runs.

**NMD rule applied after the freeze.** The curator's rule marks a premature stop more than 50 nt
before the junction of the last protein-coding exon as predicted nonsense-mediated decay. Only
protein-coding exons count, so a last exon that is entirely 3′ UTR is ignored. In the canonical
export the 2,378 pairs with a predicted-NMD record are annotated `r.0?` / `p.0?`, and the Mutalyzer
r. and p. descriptions are kept in `alt_descriptions`. The audit still compares the Mutalyzer
protein with the source, so it keeps testing the c. mapping itself.

## Traceability and resolution policy

- `source_change_id` links every record to the source change it came from.
- `individual_id = GENE_accession`, because accessions are only unique within one gene database.
- When the tracks disagree, one primary record is kept per source change. Other descriptions of the same genomic change go to `alt_descriptions`, and every final record carries a `confidence_tier` (A, B or Review).
- In-cis alleles (`c.[A;B]`) legitimately give several records per source unit and are never collapsed.

## Reproducibility discipline

- Mutalyzer, VariantValidator and GenBank responses are cached append-only.
- Before each live run, the uncached queries are forecast and reviewed.
- Each fix is committed as code, then primary outputs, then rescue outputs, then the ledger entry, so any recorded run can be re-scored from git.
