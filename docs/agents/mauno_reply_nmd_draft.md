# Reply to Mauno: NMD (draft, 2026-10-13)

Attachments:
- Output/Step11_Export/nmd_sample_for_mauno.csv (20 variants)
- Output/Step11_Export/review_workbook_for_mauno.xlsx (updated)
- Optional: Output/Step11_Export/nmd_review_for_mauno.csv (all 1,091 premature-stop variants)

## Email

Hi Mauno,

Thank you, and thanks for the NMD rule. I have now applied it:
- The boundary is 50 nt before the end of the penultimate protein-coding exon. 3′UTR-only exons are ignored; 7 of our MANE transcripts end in one.
- If the stop lies upstream of that boundary, the variant is annotated as predicted missing RNA and protein (r.0?, p.0?). The Mutalyzer protein description is kept only as an alternative description.

Of the 1,091 distinct variants with a premature stop:
- 872 (2,378 patient–variant pairs) are predicted to undergo NMD and are now r.0?/p.0?;
- 186 have the stop in the last coding exon;
- 33 have it within 50 nt of the last junction.

Since this is too many to check by hand, I attach a sample of 20. It has 10 predicted and 10 escaping variants, half frameshifts and half nonsense, all from different genes. It includes your JAK3 example and the borderline cases on both sides of 50 nt (TNFRSF13B at 52 nt, ITGB2 at 48 nt). The rule is applied in the same way to all of them; the full list is available if you want it.

The suspected source-error list in the workbook is now shorter (27 → 21 entries). The frameshift cases whose stop triggers NMD have been removed, as you pointed out.

One question: for variants that escape NMD under this rule, should we keep the truncated-protein description, e.g. CXCR4 p.(Arg334Ter) or IFNGR1 p.(Asn274HisfsTer2)? Or should those also be annotated as p.0?

Best,
Khoi
