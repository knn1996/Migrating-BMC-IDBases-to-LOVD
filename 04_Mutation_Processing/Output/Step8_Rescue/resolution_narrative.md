## Resolution of Unresolved Variants

Of the 848 distinct variants that failed all three Mutalyzer
processing tracks (NG\_IDRefseq, NM\_MANE, and NM\_IDRefseq), **566
(66.7%)** were automatically resolved by the rescue layer:
451 via MANE Select NM\_ re-mapping (correcting obsolete transcript
version references), 77 via VariantValidator against the current GRCh38
build (intronic and IVS-notation variants), 38 via a
VariantValidator probe of reference-mismatch candidates, and
0 via empirical coordinate-offset correction (±1 c.\ position
shift confirmed against GRCh38).

Of these, **149 require curator review before submission.** These are
descriptions that spelled out the deleted or duplicated sequence in the
deprecated pre-2016 style (e.g. c.1293delC rather than c.1293del). Shortening
the description discards the source record's assertion about which base was
deleted while retaining its assertion about the position, so where that base
disagreed with GRCh38 the variant now validates against the reference base
instead. Position is being trusted over allele; the affected rows are flagged
`needs_curator_review` in `resolved_unresolved.tsv` and should be confirmed
against the primary literature before deposition.

Resolution rates are reported on a distinct-variant basis in
`rescue_metrics_comparison.csv`; patient-level row counts before and after
augmentation are reported in `augment_stats.txt` and `dedup_rescued_stats.txt`.
The pre-rescue pool contained 9052 patient-variant rows, which is
not directly comparable with the 848 distinct unresolved variants
above and is therefore not combined into a single ratio here.

The residual **282** variants fall into three mechanistically
distinct groups:

- **210 reference-sequence mismatches.** The base asserted by the
  source record disagrees with GRCh38 at the stated position, either through
  strand-orientation error in the original curation or through a legacy
  reference sequence that diverges from the current assembly. These cannot be
  corrected without fabricating sequence information and represent a genuine
  limit of automated rescue.
- **67 legacy transcript-structure incompatibilities.** The
  historical c.\ numbering anchors positions — particularly intronic offsets —
  to exon boundaries that do not exist in the MANE Select transcript, yielding
  EEXONBOUNDARY on re-validation. These are not source-data errors; resolving
  them requires reconstructing genomic coordinates from the legacy RefSeqGene
  reference rather than re-mapping the c.\ description, which is outside the
  scope of the present rescue layer.
- **2 curation losses and 3 manual-rewrite cases.**
  Descriptions that lost their inserted sequence during original database entry,
  carry reversed or non-consecutive ranges, or otherwise require a targeted
  human HGVS rewrite before automated re-verification is meaningful.
