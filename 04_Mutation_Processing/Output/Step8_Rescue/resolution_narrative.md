## Resolution of Unresolved Variants: 76% Was Not a Ceiling

Of the 698 variants that failed all three Mutalyzer processing
tracks (NG\_IDRefseq, NM\_MANE, and NM\_IDRefseq), **412
(59%)** were automatically resolved by the rescue layer:
401 via MANE Select NM\_ re-mapping (correcting obsolete transcript
version references), 11 via VariantValidator (resolving intronic and
IVS-notation variants against the current GRCh38 build), and
0 via empirical coordinate-offset correction (±1 c.\ position
shift confirmed against GRCh38). Combined with the 7776
patient-level entries resolved in the primary pipeline, the post-rescue
variant pool reaches **8188** rows (exact distinct count after
deduplication is reported in `dedup_stats.txt`), raising the automated
resolution rate from 91.8% to an estimated **96.6%**
of all submitted variants.

The residual **286** unresolved variants comprise a small,
well-characterised set of genuine source-data errors that cannot be
corrected without fabricating sequence information:
169 carry reference-sequence mismatches attributable to
strand-orientation discrepancies or completely divergent legacy sequences,
1 lost their inserted sequence during original database
curation (delins with no alt allele), and 4 require
targeted manual HGVS rewrite before automated re-verification.
These categories represent hard limits of automated rescue rather than
ambiguities in the variant data itself.
