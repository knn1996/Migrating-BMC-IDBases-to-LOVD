# Audit methodology

How the post-thesis work measures whether the migrated records are right. Results per run are in
[post_thesis_validation.md](post_thesis_validation.md); the code is in
`04_Mutation_Processing/Scripts/Step10_Audit/`.

## An independent audit, not the pipeline's own validation

Mutalyzer and VariantValidator check that an HGVS description is valid on its reference. They cannot
tell whether it describes the variant the original authors reported. A separate audit
(`audit_against_source.py`) therefore re-parses the published IDbases pages (`pub.html`). It shares no
code with the extraction step, and it traces every source unit to its fate in the final table.

## The ruler

- **Source unit**: one (gene, accession, allele) entry of the published pages. The key never uses the
  remapped HGVS, which can change between runs.
- The audit covers **7,047 source units, 7,019 of them in scope**.
- `individual_id = GENE_accession`, because accessions are only unique within one gene database.

### Counting definitions

The thesis figures count *distinct variants* (and patient–variant pairs). The audit counts *source
units*. The two are not interchangeable: one distinct variant can carry many source units, and one
source unit can need several records (in-cis alleles, `c.[A;B]`).

## Three axes, always reported together

- **Recovery**: share of in-scope source units present in the final table.
- **Best-record accuracy**: the protein consequence of the final record matches the source p.
  description (`protein_good`).
- **Precision**: each final patient–variant pair is consistent with the protein of *its own* source
  change (`pair_consistency`), and no individual has more final pairs than source changes
  (`excess_pairs`).

For each mismatch the audit records the kind: a different change, a shifted position, a different
frameshift start or stop, or a different change type at the same position.

## c./p. coherence

A pair is coherent when the protein accession of its p. is the product of the transcript of its c.
Products come from the GenBank cache (`/protein_id`) and, for MANE Select transcripts, from cached
Mutalyzer responses. A pair whose transcript product is not known offline counts as unknown, and
coherence is reported over the known pairs.

## Notation equivalences

Two published notations are accepted as equivalent to the HGVS form:

- **E1, immediate stop**: a source `XnfsX1`, a bare `XnfsX`, or `XnfsX<n>` with n equal to the start
  position (absolute stop numbering) equals the final `p.(XnTer)`.
- **E2, insertion as duplication**: an insertion whose right flank is the duplication start, with the
  duplication length and the duplication end residues as its first and last inserted residues.

## The run ledger

`audit_ledger.py` records each fix as one run against committed outputs: a baseline (Run 0) and 13 fix
runs (`Output/Step10_Audit/history/`). For every run it lists gains, regressions, GOOD→BAD and
BAD→GOOD transitions, and pair transitions against the previous run.

Every pair is linked to its own source change (`source_change_id`). The link is rebuilt by replaying
the run's own deduplication step on its committed inputs, and it must reproduce the stored links.

When the audit itself changes, `--rescore-all` re-scores the whole recorded history from git, and every
existing column must reproduce exactly before anything is written.

## Suspected source errors

Some remaining mismatches are errors in the published record rather than in the pipeline. A source
change is a **candidate** when its own c. and p. (or g. and c.) contradict each other:

- tested on the MANE Select sequence at shift 0 and at the gene's verified legacy shift, and on every
  cached non-MANE transcript;
- and not explained by anything on our side (an isoform, an unverified shift, a non-GenBank coding
  sequence, an intronic record or a record kept in legacy numbering), unless the conflict is intrinsic
  (for example a frameshift-length indel given as a missense change).

Today there are **24 candidate source changes (29 pairs)**:

| Category | Changes |
|---|---|
| Notation: right position, wrong change type, extent or stop | 19 |
| c./p. mismatch: the c. cannot give the stated p. | 2 |
| g./c. mismatch: the g. and c. describe different changes | 2 |
| Exon boundary: a c. range crosses an exon–intron junction | 1 |

16 are high severity (the protein consequence differs) and 8 low (only the frameshift stop count or the
g. differs). They are pending review by the original curator. The record-level list is kept out of the
public repository until then; the ledger reports the candidates in separate columns
(`pairs_source_error*`, `pair_consistency_excl_source_error_pct`) without changing any other measure,
and the canonical export flags only source errors the curator has confirmed (none yet).
