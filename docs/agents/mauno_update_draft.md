# Update for Mauno (draft, 2026-10-12)

Figures from Output/Step10_Audit/history/audit_history.tsv, Run 0 → Run 10.

## Short version

Since the defence, I have focused on making sure that every record going into LOVD is correct, rather than on resolving more variants.

1. **Independent audit.** I wrote an audit that re-reads the published IDbases pages separately from the pipeline. It traces each of the 7,019 in-scope source variant units (gene, accession, allele) to the final table. For each record it checks three things: whether the unit is present, whether our protein consequence matches the source p., and whether each patient–variant pair matches its own source change. It also checks that the c. and p. come from the same transcript. Every fix is scored as one run against the previous one, counting gains, regressions, and correct→incorrect and incorrect→correct protein changes.
2. **Coverage is never accepted at the cost of consistency.** The first fix raised coverage, but consistency fell. That pattern shaped how every later fix was judged.
3. **Mixed historical numbering.** Many pre-HGVS papers number c. from the start of their own cDNA record. The result is a constant per-gene shift equal to that record's 5′ UTR length (e.g. RAG1 −112, AIRE −127, AICDA −76). 28 genes have a verified shift, and 18 of them mix both conventions in one database. A shift is applied per source change, and only when at least two independent changes confirm it. This replaced the old global ADA −95 offset, which also moved ADA changes that were already correct.
4. **Pipeline bugs fixed**, each measured against records that were already right as well as those that were wrong: extraction (compound heterozygous alleles, merged records), deduplication, BED-to-HGVS conversion, IVS positions, choosing one primary record per source change across tracks, and taking c. and p. from the same transcript.

| | Before (Run 0) | Now (Run 10) |
|---|---|---|
| Source units in final table | 83.33% | 91.71% |
| Protein matches source (best record) | 90.79% | 98.42% |
| Patient–variant pairs consistent with source | 87.04% | 98.88% |
| Excess (wrong extra) patient–variant pairs | 653 | 8 |
| c./p. from the same transcript | 75.34% | 100% |

Across the 10 runs, 459 units went from an incorrect to a correct protein. 18 went the other way, all in the first run, and all were found and fixed later.

These numbers count source units, so they are not directly comparable with the thesis figure of 76.2%, which counts distinct variants.

## What's next
- **Suspected source errors.** I am compiling a list of entries where the source itself looks internally inconsistent: for example, AIRE A0082 gives p.S64del for a 36-nt deletion, and CD40L C0047 has a c. range that crosses an exon–intron boundary. I will send it as a CSV for you to check.
- **MANE Select.** I will move the remaining records that are not on MANE Select onto it, keeping the original descriptions as alternatives.
- **Export.** I will then build one canonical export table, so that it can be mapped to whatever format the LOVD director specifies.

## Questions
- Has the LOVD director said anything about the submission format?
- IGHG2 and SH2: I plan to drop IGHG2 (the heavy-chain reference from the BTKbase migration), as we discussed, and to leave SH2 to the BTKbase curation. Is that still right?
- Review-tier records: I will send them as an xlsx review workbook, as during the thesis. Is that format fine for you?
