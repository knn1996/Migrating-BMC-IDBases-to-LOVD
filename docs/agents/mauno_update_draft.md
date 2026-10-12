# Update for Mauno (draft, 2026-10-12)

Figures from Output/Step10_Audit/history/audit_history.tsv (Run 0 → Run 13) and Output/Step11_Export/canonical_export.tsv.
Attachment: Output/Step11_Export/review_workbook_for_mauno.xlsx

## Email

Since the defence, I have focused on making sure that every record going into LOVD is correct, rather than on resolving more variants.

1. **Independent audit.** I wrote an audit that re-reads the published IDbases pages separately from the pipeline. It traces each of the 7,019 in-scope source variant units (gene, accession, allele) to the final table. For each unit it checks three things: whether it is present, whether our protein matches the source p., and whether each patient–variant pair matches its own source change. It also checks that c. and p. come from the same transcript. Every fix is scored against the previous run, counting gains, regressions, and correct→incorrect and incorrect→correct protein changes.
2. **Coverage is never accepted at the cost of consistency.** The first fix raised coverage, but consistency fell. That pattern shaped how every later fix was judged.
3. **Mixed historical numbering.** Many pre-HGVS papers number c. from the start of their own cDNA record, which gives a constant per-gene shift (e.g. RAG1 −112, AIRE −127, AICDA −76). 28 genes have a verified shift, and 18 of them mix both conventions in one database. A shift is applied per source change, and only when at least two independent changes confirm it. This replaced the old global ADA −95 offset.
4. **Pipeline bugs fixed**, each measured against records that were already right as well as those that were wrong.
5. **A stopping point.** Once each fix improved less than 0.5% of the source units, I froze the pipeline. Records we cannot vouch for are kept at a Review tier with a stated reason, rather than patched one by one.

| | Before | Now |
|---|---|---|
| Source units in the final table | 83.33% | 91.74% |
| Protein matches source (best record) | 90.79% | 98.43% |
| Patient–variant pairs consistent with source | 87.04% | 98.93% (99.39% excluding suspected source errors) |
| Excess (wrong extra) patient–variant pairs | 653 | 8 |
| c./p. from the same transcript | 75.34% | 100% |
| Pairs on the exact MANE Select transcript | – | 6,518 of 6,562 |

These numbers count source units, so they are not directly comparable with the thesis figure of 76.2%, which counts distinct variants.

The export table now has 6,558 patient–variant pairs, in 5,712 individuals and 128 genes.
- Each row carries the GRCh38 g., the MANE c./r./p., the allele, phase and zygosity, the notation as published, a confidence tier, and the patient metadata.
- Gene symbols are updated to current names (e.g. ELA2 → ELANE, CD40L → CD40LG, BIRC4 → XIAP).

**Attached workbook, for your review:**
- **Review_tier:** 76 patient–variant pairs that we will not submit as confirmed, each with the reason. In most of them our genomic mapping lands outside the coding region, while the published c. is exonic.
- **Suspected_source_errors:** 27 published entries whose notations contradict each other. For example, AIRE A0082 gives p.S64del for a 36-nt deletion. A few others describe no change at all, such as CYBB A1432 c.998C>C.

Each sheet has a Mauno_decision column (keep / correct / drop). Your answers can be read back in directly.

## Questions
- Has the LOVD director said anything about the submission format? The export can be mapped to it once it is known.
- IGHG2 and SH2: I plan to drop IGHG2 (the heavy-chain reference from the BTKbase migration), as we discussed, and to leave SH2 to the BTKbase curation. Is that still right?
- The Screening, Genetic origin and Effect fields are empty in the source. Should they stay empty, or is there a default you would like?
