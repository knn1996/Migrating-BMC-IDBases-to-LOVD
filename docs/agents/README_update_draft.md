# Draft: README section "Post-thesis correctness work"

Status: living draft. The planner extends it after every recorded run; it is merged into the root README.md once the user judges the work at a satisfactory stopping point.
Figures come from Output/Step10_Audit/history/audit_history.tsv (currently Runs 0–11). Recompute them from the history after each run; never hard-code.

---

## Post-thesis correctness work (2026)

The thesis results above are reported as defended. After the defence, the work moved from *how much* of IDbases can be migrated to *whether what is migrated is right*. A LOVD submission is only useful if every record is correct, so precision now takes priority over recovery.

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

The remaining mismatches are increasingly errors in the source rather than in the pipeline. For example, AIRE A0082 reports p.S64del for a 36-nt deletion, and CD40L C0047 gives a c. range that crosses an exon–intron boundary. The audit now flags 24 suspected source errors (29 pairs, pending review by the original curator) in separate columns, without changing any existing measure. Excluding them, pair consistency in Run 11 is 99.37%.

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

From Run 0 to Run 11:
- Recovery rose from 5,849 to 6,439 source units (83.33% → 91.74%).
- Inconsistent pairs fell from 745 to 61.
- Excess pairs fell from 653 to 8.
- Protein BAD→GOOD totalled 459 units, with all GOOD→BAD (18) in Run 1.
- c./p. coherence rose from 75.34% to 100% of checkable pairs.

### Traceability and resolution policy

- `source_change_id` links every record to the source change it came from.
- `individual_id = GENE_accession`, because accessions are only unique within one gene database.
- When the tracks disagree, one primary record is kept per source change. Other descriptions of the same genomic change go to `alt_descriptions`, and every final record carries a `confidence_tier` (A, B or Review).
- In-cis alleles (`c.[A;B]`) legitimately give several records per source unit and are never collapsed.

### Reproducibility discipline

- Mutalyzer, VariantValidator and GenBank responses are cached append-only.
- Before each live run, the uncached queries are forecast and reviewed.
- Each fix is committed as code, then primary outputs, then rescue outputs, then the ledger entry, so any recorded run can be re-scored from git.

### README statements now out of date

- Limitations: "`ADA` requires a documented coding-coordinate offset of −95". This was replaced by the per-change remap (Run 9).
- The headline counts (76.2%, 2,240 distinct variants) are thesis-era and use distinct-variant counting. The audit figures above use source units (7,019 in scope). State both definitions, or add a current-state line.
