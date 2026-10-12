# Handoff — IDbases → LOVD pipeline (source of truth; updated 2026-10-12, after Run 13; pipeline frozen)

This file is the single source of truth. The planner (`thesis-plan`) updates it after each recorded run. The older copy in the claude.ai project (claude/handoff_2026-10-11.md) is frozen and may be stale.

## Context
MSc thesis (Lund, supervisor Mauno Vihinen) is DONE. Goal: finish migrating IDbases from MUTbase to LOVD 3.0 (GRCh38, MANE Select, HGVS) so the data can be submitted. Priority: correctness of what goes into LOVD (precision) over recovery. Repo: D:\Assignment\Thesis = /mnt/d/Assignment/Thesis in WSL2 Ubuntu.

Submission format: Mauno is contacting the current LOVD director; for a database this size LOVD likely defines a custom submission format. Do NOT build a generic LOVD import file. Do format-independent correctness first, then one canonical export table (one row per patient–variant pair with everything any format could need); the final format becomes a mapping job once known.

Sessions: executor `thesis-exec` (Claude Code, repo root, Remote Control on), planner `thesis-plan` (Claude Code, repo root, see docs/agents/PLANNER.md). Shared rules in CLAUDE.md.

Shell setup for the user: `tmux new -A -s thesis`, then `thesis` (bash function: cd to Scripts + conda activate snakemake).
Opening routine: (1) tmux session `thesis` (`tmux new -A -s thesis`) holds the Claude sessions: one window for the executor, Claude Code session named `thesis-exec`, Remote Control on; one window for the planner, Claude Code session named `thesis-plan` (`/rename thesis-plan` if needed), told "Read docs/agents/PLANNER.md and follow its start-of-session steps"; (2) tmux session `pipeline` (`tmux new -A -s pipeline`, then `thesis`) for Snakemake runs. Both sessions run Opus for everything (no per-task model switching).

Pipeline (Snakemake 9.27, env `snakemake`, run from 04_Mutation_Processing/Scripts):
- `Snakefile`: Step1a extraction → … → Step7 merge/dedup
- `rescue.smk`: Step8 rescue, resolve_tracks, build_final_lovd, VariO; Step9 NMD (fetches missing GenBank transcripts from NCBI); Step10 audit
- Always `--cores 1`. Force code changes with `-R <rule>`. Commit primary outputs before the rescue run (build_final_lovd refuses otherwise).
- audit_pairs.tsv is produced by the ledger (throwaway ledger run before committing rescue outputs).
- Forecast tools don't cover NMD GenBank fetches (backlog). vv_client.py pings VariantValidator at startup, so every rescue run needs network even with 0 queries.
- matplotlib missing in the snakemake env → audit_trend.png stale: `conda install -n snakemake -c conda-forge matplotlib` (network; user runs).

## Audit ledger
- Scripts/Step10_Audit/audit_ledger.py (--rescore-all, --allow-move); history in Output/Step10_Audit/history/.
- Key = (gene, accession, allele) from pub.html; 7,047 source units, 7,019 in scope.
- Axes: recovery (in_final), best-record accuracy (protein_good), precision (pair consistency vs own source change; excess pairs), c./p. coherence (p.'s protein is the product of c.'s transcript; GenBank /protein_id + cached Mutalyzer for MANE v1.5; ambiguous = unknown).
- Notation equivalences: E1 immediate stop (fsX1 / bare fsX), E2 ins≡dup, parsing fixes. No FS_STOP_DIFF rule.
- Runs 1–5 are rescored on a reconstructed pre-Run-6 extraction (99_Archive/…/all_mutations_pre_run6.tsv).
- source_error flag: mechanism present, SOURCE_ERRORS empty. Parked until Run 8.

| Run | Label | In final | Protein good | Pair consistency | Excess | Pairs | Coherence |
|---|---|---|---|---|---|---|---|
| 0 | pre-fix baseline | 83.33% | 90.79% | 87.04% | 653 | 6,589 | 75.34% |
| 1 | extraction fix | 89.84% | 92.74% | 85.50% | – | 7,508 | 74.55% |
| 2 | dedup fix + individual_id | 90.03% | 93.06% | 86.51% | – | 7,428 | 74.14% |
| 3 | bed_to_mutalyzer delins/dup fix | 89.98% | 93.11% | 89.93% | 733 | 7,171 | 73.18% |
| 4 | IVS position fix | 89.98% | 93.11% | 89.93% | 731 | 7,169 | 73.17% |
| 5 | R2 track resolution | 89.98% | 96.20% | 96.41% | 23 | 6,456 | 73.69% |
| 6 | dedup transcript split + stage-2 links + accession fix | 89.98% | 96.20% | 96.65% | 8 | 6,441 | 73.69% |
| 7 | Track B c./p. coherence | 89.98% | 96.20% | 96.65% | 8 | 6,441 | 100.00% |
| 8 | legacy c. numbering remap | 89.98% | 98.19% | 98.66% | 8 | 6,441 | 100.00% |
| 9 | ADA per-change remap | 90.14% | 98.33% | 98.8% | 8 | 6,452 | 100.00% |
| 10 | full-scan legacy remap | 91.71% | 98.42% | 98.88% | 8 | 6,562 | 100.00% |
| 11 | literal remapped records not eligible | 91.74% | 98.42% | 98.91% | 8 | 6,562 | 100.00% |
| 12 | MANE normalization | 91.74% | 98.43% | 98.93% | 8 | 6,562 | 100.00% |
| 13 | Review demotions | 91.74% | 98.43% | 98.93% | 8 | 6,562 | 100.00% |
Coherence ruler: 3fc3b8a (code), 4c7c6b1 (history). Run 6: dcde390, 9e0e2b7, f6de49d, d82a5c5, 7543cbb (rescue), aa6c043 (history). Run 7: 0c255fe (code), 1e12fd6 (primary), d6b09bb (rescue), 3b4cc34 (history, run_id 20261012_d6b09bb). Run 8: 12945a1 (rescue grouping fix), f786738 (code), 50833c3 (primary), f64272c (rescue), 3de036a (history, run_id 20261012_f64272c). Run 9: bdc9ac4 (code), d1879a9 (primary), 30fdf51 (rescue), 9c4c1f9 (history, run_id 20261012_30fdf51). Run 10: 1bb3eda (code), 78b539b (primary), a3fe960 (rescue), a216ed4 (history, run_id 20261012_a3fe960). Run 11: 2d9d79c (code), 98bdc7d (primary), 62c8f01 (rescue), 4ee52a4 (history, run_id 20261012_62c8f01). Run 12: 195a769 (code), 7cad191 (rescue; no primary change), 25e50c0 (history, run_id 20261012_7cad191). Run 13: 57f1196 (code), bd661f1 (rescue), 1df4b49 (history, run_id 20261012_bd661f1), 607e00a (docs). Tiers A 6,482 / B 4 / Review 76. Nothing pushed.

## Fixed so far
1. Extraction: compound-het alleles, glued records, brackets (Run 1)
2. Position-only XM/NM merge (Run 2)
3. individual_id = GENE_accession (Run 2)
4. bed_to_mutalyzer delins/dup/minus-strand (Run 3)
5. IVS positions inherited from neighbours (Run 4)
6. Track A allele labels + source_change_id (Run 5 prep)
7. Track priority never fired → R2 track resolution with alt_descriptions and confidence_tier A/B/Review (Run 5)
8. Transcript-blind stage-1 dedup, resolve_tracks follows stage-2 merges, trailing ";" in accessions (Run 6)
9. Track B c. taken from the normalized MANE transcript (same as p.); old c. kept as alt_c_hgvs → alt_descriptions (Run 7)
Measurement: count definitions, precision metric, own-source-change ruler, c./p. coherence axis, --rescore-all.

## Run 7 — DONE (Track B c./p. coherence)
- Cause: run_mutalyzer_NM.py took p. from the normalized MANE transcript but c. from pick_mane's fallback (another transcript).
- Result: coherence 73.69% → 100% of known (coherent 3,926 → 5,703, incoherent 1,402 → 0, unknown 1,113 → 738). in_final, protein_good, excess, consistency, pairs unchanged. Gains 0, regressions 0, GOOD↔BAD 0.
- Pair transitions: 9 removed + 9 inconsistent = the same 9 pairs re-keyed to MANE c. strings, DIFFERENT before and after (C2_C0002 ×2, C2_C0003 ×2, LYST_C0027, LYST_C0028, RFX5_R0005/R0007/R0008). They belong to Run 8 (legacy numbering).
- Track B on exact MANE 1,829/1,829; all final records 2,464/2,663. Dedup identical_change 449 → 4; rescue skip_different_change_same_unit 23 → 74 (primary 5 → 56).
- NMD fetched 59 GenBank records (146 → 205 files). Real outputs = offline replay byte for byte.
- Note: lovd_review_final.xlsx needs openpyxl (missing in the base env; only affects offline replays).

## Run 8 — DONE (legacy c. numbering remap)
- Result: BAD→GOOD 107, regressions 0, inconsistent 183 → 73, 114 pairs re-keyed (0 lost). Live: 56 Track B + 16 rescue Mutalyzer, 32 VV.
- Accepted residuals: 3 source_error candidates (long in-frame indels) and 1 FS_STOP_DIFF residual (detail private); WAS W0009 consistent by E1 at pair level but unit-level protein_check ignores equivalences (ledger gap).

### Run 8 scoping notes (report-only, 2026-10-12; scratch /tmp/claude-1000/run8_*)
- 114 of 183 inconsistent pairs (80 changes) fixable by a verified gene-wide shift (≥2 independent changes; MANE pos = source c. + k): AICDA −76, RAG1 −112, TCIRG1 −60, WAS −34, C9 −4, RAG2 −1201, CIITA −138, RFXAP −116, C2 −36, RFX5 −161, CFH −73, C3 −60, C5 −12, LYST −189, IL7R −22, CD40L −59, IGLL1 −118. 72 pairs carry the literal legacy c.; 41 come from Track A errors.
- ADA (6 pairs): over-remap. Source c. is already MANE; the pipeline's offset path shifts it by −95. A different bug.
- Excluded: JAK3 isoform (−167 codons), C1QB (no common shift on MANE), CYBB band artifact, 4 k=0 stop/type cases.
- Forecast: about 61 uncached Mutalyzer queries (upper bound), 0 VV, 0 GenBank for the fix genes.
- source_error candidates: 29 changes / 41 pairs (13 irreconcilable, 15 single-change, 1 untestable).
- Design: shared remap module + committed per-gene shift table; targeted hooks in generate_mutalyzer_input_NM.py, the rescue input generator and resolve_unresolved (ADA).
- Ceiling if all 114 become consistent: 5,399/5,468 = 98.74%.
- Implemented (uncommitted): legacy_remap.py + build_legacy_remap.py; refs/legacy_remap.tsv keyed per source change (gene, accession, source c.), 114 entries; 15 genes kept, RFXAP/CD40L/SERPING1 dropped by the evidence check. Hooks: generate_mutalyzer_input_NM.py (MANE only), generate_rescue_input.py.
- Cached-only sim: protein_good 96.83%, consistency 97.27%, BAD→GOOD 34, 0 regressions, excess 8; exposed a rescue grouping bug (classify_unresolved groups by query c., merge_with_rescued expands with the first row's sysname) → 13 unlinked.
- User decisions: grouping fix is a separate code commit first (byte-identical on Run 7); queries staged: 56 Track B Mutalyzer in the primary run, then re-forecast rescue (≤16 Mutalyzer, ≤16 VV) for separate approval.
- ADA −95 comes from hard-coded CODING_OFFSET in generate_mutalyzer_input_NM.py (for the ADA run).
- Token hygiene: /clear both sessions after each recorded run.
- Committed: 12945a1 (rescue grouping fix; 77 augmented allele labels corrected, finals unchanged), f786738 (Run 8 code). Primary: 56 live Track B Mutalyzer queries approved.

## Full numbering scan (2026-10-12, scratch/run8/scan/scan_changes.tsv, scan_genes.tsv)
- Cause: IDbases stores c. as published; pre-HGVS papers counted from the cDNA record start (offset = that record's 5' UTR length). Same coding sequence, different origin. p. style tracks it (one-letter = old, three-letter = HGVS) in every mixed gene except CYBB.
- 7,920 source changes / 133 genes: k0 4,216, shift 701, both 3, isoform 63, unverified 226, neither 360, untestable 2,351. 28 genes with a verified shift (18 mixed, 10 pure).
- New mixed genes: AIRE −127 (188/157), C1S −202, CD40L −56 (85/98), CYBB −14 (doubtful), IKBKG −185, PRF1 −101, RFXANK −417, TAP2 −28. New pure-shift: CD40 −47, CFP −242, NCF1 −12 (doubtful).
- 110 recoverable not-in-final units (CFH 20, RAG1 19, AICDA 14, IKBKG 13, RAG2 13, WAS 11, TCIRG1 9, ...). 1 wrong-kept Track A record (AIRE). 14 shift-class changes still inconsistent. Isoform PTPN11 −71 codons: all consistent, no action.
- Most shift-class units in mixed genes are already right via Track A genomic mapping.

## Next after Run 8
0. (done, Run 9) ADA: CODING_OFFSET removed; 81 −95 entries + 1 manual + 11 excluded in legacy_remap. Gains 11, BAD→GOOD 7, regressions 0, consistency 98.80%. Live: 3+8 Mutalyzer, 16 VV.
0b. (done, Run 10) Extend the legacy remap with the full scan. Result: +110 in_final (exactly the 110 scan-recoverable units: CFH 20, RAG1 19, AICDA 14, IKBKG 13, RAG2 13, WAS 11, TCIRG1 9, CD40 3, TAP2 3, AIRE 2, RFXAP 2, CIITA 1), regressions 0, GOOD→BAD 0, BAD→GOOD 3 (RFXAP c.484del → c.368del), new inconsistent 0, 17 pairs re-keyed (0 lost), CD40 −47 confirmed (3/3 MATCH). Live: 170 Track B + 36 rescue Mutalyzer, VV, 1 GenBank (CD40). Procedure: throwaway ledger dirties audit_pairs.tsv; restore it to HEAD before the recording ledger.
   Original plan: extend the legacy remap with the full scan (110 recoverable units, new mixed genes, AIRE Track A record).
   Scoping (2026-10-12, partial; executor hit usage limit): scratch 99_Archive/04_Mutation_Processing/scratch/run10/ (table_main = 694 rows, keeps all 207 existing rows; table_all adds CYBB/NCF1). Sim copy: planner job tmp sim_main (temporary).
   - 509 shift units not yet in remap, 19 genes, all verified. Builder gotcha: ADA must be in --all-changes.
   - Cached-only sim (table_main): in_final 6327 → 6367, protein_good 98.39%, consistency 98.86%, excess 8, coherence 100%. Gains 44 (RAG1 12, AICDA 8, CFH 8, RAG2 8, WAS 6, AIRE 2), GOOD→BAD 0, BAD→GOOD 0.
   - 4 regressions (RFXAP R0001/R0003/R0004, IKBKG I0018) = uncached Track B queries (no_hgvs offline). Planner verified offline: remapped c. reproduces source p. exactly (RFXAP c.368delG → p.S123fsX137, c.163C>T → p.Q55X; IKBKG c.1217_1218insA → p.D406fsX419). The 3 RFXAP pairs removed were wrong in Run 9 (VV c.484del, DIFFERENT).
   - Forecast: 171 distinct uncached Track B Mutalyzer queries (CD40L 51, LYST 19, AIRE 16, PRF1 13, CFP 11, IKBKG 10, RAG1 9, TCIRG1 7, CFH 6, RAG2 5, WAS 5, AICDA 4, RFXANK 3, RFXAP 3, C1S 2, CD40 2, CIITA 2, TAP2 2, C3 1); Track C 0. Rescue re-forecast after the primary run. NMD GenBank: ≤1 (CD40 NM_001250.6 not cached).
   - CYBB/NCF1: exclude (no gain at cached level; CYBB doubtful, NCF1 residual is FS_STOP_DIFF).
   - AIRE (1 unit): Track A record is right at the −127 shift; the source p. disagrees → source_error candidate (detail private), stays inconsistent. The remap removes the wrong Track B alternative c.319_354del.
   - Committed 1bb3eda (refs/legacy_remap.tsv 694 rows = 207 old identical + 487 new; legacy_shift_evidence.tsv). Installed-table replay = sim_main byte for byte.
   - Forecast tools gap: check_mutalyzer_cache.py / check_rescue_cache.py ignore the remap unless LEGACY_REMAP_TSV=$PWD/refs/legacy_remap.tsv is exported (otherwise 21 false ADA hits, 171 missed). Fix as a separate commit after Run 10.
   - Primary: `snakemake -s Snakefile --use-conda --cores 1 -R mutalyzer_input_NM` (171 live Track B). Then commit primary, re-forecast rescue (≤36 Mutalyzer, VV 0 + ≤36 anchors) for separate approval.
   - Expected: 4 regressions back GOOD, in_final ≥ 6,371, GOOD→BAD 0, excess ≤ 8, coherence 100%, no new inconsistent pairs. Gate: new inconsistency or ESEQUENCEMISMATCH in RFXAP/TAP2/C3/CD40 → drop that gene's rows and rerun before committing primary. CD40 CDS from longest ORF (no GenBank): review its records.
   - translation_test.py stays uncommitted until source_error work.
   - Primary live run done; committed 78b539b. 170 new Track B queries: 133 ok, 36 ENOSELECTORFOUND (to rescue), 1 ESEQUENCEMISMATCH (CD40L (1 unit): source c. range crosses an exon/intron boundary; the Track A record is correct and kept; source-notation candidate, detail private).
   - Cached-only rescue projection: in_final 6383, protein_good 98.40%, consistency 98.86%, excess 8, pairs 6508, coherence 100%. Gains 60 (RAG1 19, RAG2 13, AICDA 8, CFH 8, WAS 6, TAP2 3, AIRE 2, CIITA 1), regressions 4 (pending uncached rescue), GOOD→BAD 0, new inconsistent 0; 13 pairs re-keyed (bogus intronic keys → coding keys).
   - Rescue forecast: 36 Mutalyzer, 36 VV probe, 0–37 VV resolve, ≤1 GenBank (CD40). Expected: in_final ≥ 6387, 4 pending back GOOD, GOOD→BAD 0, excess ≤ 8.
1. (done) Run 8 — Legacy c. numbering remap (~70 pairs; RAG1 +112, AICDA +76, CFH +73, TCIRG1 +60, WAS 34–36, LYST 190, RFXAP 116) + ADA over-remap. Needs Utilities/translation_test.py (draft, uncommitted): forecast missing transcripts first (fetch = network step), per-gene shift verified when ≥2 independent changes reconcile (C1QB −2 codons, RFX5 −40 codons are isoform numbering), include small frameshifts/in-frame indels, one shared remap function for c. tracks and rescue. ~25 queries + transcripts. Also yields source_error flags.
2. Run 9 — MANE normalization of the ~200 records still off exact MANE after Run 7 (Track A on pick_mane fallback transcripts, Track C non-MANE NMs, VariantValidator/rescue). Recompute scope and forecast. First find why pick_mane never matches MANE_SELECT (needs 1 live raw Mutalyzer response; user runs it). Only IGHG2 has no MANE Select (ELA2 = ELANE NM_001972.4). Keep originals in alt_descriptions.
3. Patient metadata: extract_patient_metadata.py glued-record bug (metadata flags ~3,147).
4. Canonical export table; then the LOVD format mapping once the director's format is known.

## PIPELINE FROZEN after Run 13 (2026-10-12)
No further Step1–Step10 code changes without an explicit user decision. Next phase: canonical export table (one row per patient–variant pair; display VariantOnTranscript/DNA, never _dedup_key; exclude review_reason "no change described (c.=)"; allele/zygosity from augmented_merged_variants via source_change_id), Review xlsx for Mauno (the 76 Review pairs with reasons), tier semantics (921 uncheckable single-track tier-A pairs), then merge the README section.

## Export phase (started 2026-10-12)
- User: build the canonical export + Review workbook first, merge the README after.
- Step11_Export (standalone offline script, outside the frozen Steps 1–10): canonical_export.tsv (one row per pair, c.= excluded, expected 6,558 rows; allele/phase per pair via source_change_id; HGNC symbols; published_as; tier, review_reason, source_error flag; Screening/Genetic_origin/Effect left empty for Mauno) and review_workbook_for_mauno.xlsx (README, 76 Review pairs, suspected source errors, counts; Mauno_decision column).

- Export built (Step11_Export/build_canonical_export.py): 6,558 rows, 43 columns; 0 copy mismatches, 0 unlinked. Zygosity hom 1,713 (source "Allele 1 and 2"), het 1,635, hemi 782, unknown 2,428 (X with sex not M 1,094; allele 0 931; no other allele 403). Gene symbols offline from alias.csv + MANE v1.5: BIRC4→XIAP, CD40L→CD40LG, ELA2→ELANE, MAPBPIP→LAMTOR2, MRE11A→MRE11, NP→PNP, TAZ→TAFAZZIN (IGHG2 unchecked). Workbook: README, Review_tier 76 (genomic 41, literal 6, no change 4, track conflict 25), Suspected_source_errors 27, Counts. xlsx stays untracked (regenerable). Garbled multi-individual sex strings flagged (metadata_flag), not interpreted.
- Backlog (frozen pipeline): 7 truncated source_change_id c. tokens (c.501 for c.501+1G>C; CYBB A1416).

- Export committed 0071f6a (44 columns; metadata_flag on 46 garbled-sex rows). Workbook at Output/Step11_Export/review_workbook_for_mauno.xlsx (untracked). Mauno email draft updated to Run 13 (docs/agents/mauno_update_draft.md; export 6,558 pairs, 5,712 individuals, 128 genes). README merge sent to the executor.

## Stopping rule (user decision, 2026-10-12)
Guard against overfitting the pipeline to the audit ruler.
1. Known-wrong records in final must be resolved: repaired by a general rule, or demoted to Review tier / dropped. Prefer demotion when a repair needs special cases.
2. A pipeline run needs a general fix (whole class: verified gene-wide shift, step bug) that fixes ≥ 0.5% of in-scope units (≥ 35 of 7,019) or removes known-wrong records. Smaller or one-off cases go to the Review workbook for Mauno.
3. No per-unit patches beyond the evidence-backed legacy_remap table (≥2 independent changes per gene).
4. FREEZE after Run 13. Remaining items (35 possibly-ours inconsistent pairs, C1QB isoform, CYBB, VV rescue quality, residual excess, AIRE Bug 2, NMD validation) go to the Review workbook unless they pass rule 2.
Remaining runs: audit option A (DONE: e921821 code+refs, e078d72 history; all pre-existing columns byte-identical across 11 runs; Run 10 pairs_source_error 29, consistency excl. source errors 99.34% = 5538/5575; our_record_wrong pairs stay inconsistent) → Run 11 DONE (RFXAP R0006/R0007 → NM_000538.4:c.28_34dup, consistent; gains 2, regressions 0, inconsistent 61, excl-SE consistency 99.37%; procedure note: rescue refuses until primary is committed — give the rescue command only after the primary commit is confirmed) (scoped: a general positional remap would fix ≤9 units, below the 35 threshold → rejected. Run 11 = rule (a): a literal record of a legacy_remap change is not eligible, so the unit goes to the remapped rescue; hits RFXAP R0006/R0007 only, 1 rescue Mutalyzer + ≤2 VV. Rule (b): demote the 6 other our_record_wrong changes (IKBKG I0007, TCIRG1 T0004/2, T0005, T0013, WAS W0005, RAG1 R0006) to Review, folded into Run 13) → Run 12 DONE (107 Track C relabelled, 8 normalized live, 16 kept for Run 13; exact-MANE pairs 6,104 → 6,518 of 6,562; LYST C0017 BAD→GOOD: it was XM numbering, so drop it from the pipeline ticket list; coherence unknown 704 → 699; pick_mane fixed) → Run 13 DONE and recorded (code 57f1196; live run matches the sim: tiers A 6,482 / B 4 / Review 76; 51 pairs with review_reason: genomic 41, literal 6, no change 4; all rates unchanged) (blind-spot check 2026-10-12, scratch …/mane_norm/blind.py: 43 tier-A pairs whose g. projects deep intronic (19), flank (14) or outside the gene (10) while the source c. is exonic; 34 invisible to the consistency ruler; cause resolve_tracks tier.get(pk, "A") for single-track records. Demote these 43 to Review, never remap; drop the 4 c.=/g.= records (source REF = ALT: 4 allele units, added to the private Mauno CSV); scope a global single-track-uncheckable → B change as report only) → FREEZE → canonical export table + Review xlsx → merge README section.

## Run 13 decisions (2026-10-12)
- Genomic demotion rule uses the pair's own source c.: 41 pairs (not the 43 accession-level count; JAK3 J0034 and UNC13D U0069 agree with their own intronic/UTR source c.).
- c.= records stay in the pipeline at Review ("no change described (c.=)"); the canonical export excludes them. Reason: source units are allele-level, so excluding them would lose 4 units as ledger regressions (1 CYBB and 3 CFH allele units; detail private).
- A global single-track-uncheckable → B default would move 921 tier-A pairs (NM_MANE 713, VV 166, NG 24, mane_remap 17, NM_IDRefseq 1). Not changed; decide at export (tier semantics for Mauno).

## Next after Run 10
- (a) Small separate commit: forecast tools default LEGACY_REMAP_TSV to refs/legacy_remap.tsv.
- (b) source_error flag (Run 11 candidate): residuals (detail private), plus the 29-change / 41-pair candidate list from Run 8 scoping. Utilities/translation_test.py belongs here.
   Scoping done (2026-10-12, scratch 99_Archive/04_Mutation_Processing/scratch/source_error/, classify.py). 63 inconsistent pairs → 28 explained by source errors (24 changes; Mauno file = private/source_errors_for_mauno.csv, 16 high / 8 low). 7 changes / 8 pairs have a wrong record of ours (literal numbering in a verified-shift gene): IKBKG I0007, TCIRG1 T0004/2, T0005, T0013, WAS W0005, RAG1 R0006, RFXAP R0006/R0007 → Run 11 candidate: the remap does not reach Track C / VV / surviving literal records. Audit change: option A approved by the user (new columns, existing columns frozen; source_errors.tsv status = candidate until Mauno reviews). translation_test.py: archive approved. Mauno CSV is source-only (gene, accessions, allele, source g/c/p, issue, severity, why_flagged).
- (c) MANE normalization of records off exact MANE (item 2 below), then patient metadata, then the canonical export table.
   MANE scoping done (2026-10-12, scratch …/scratch/mane_norm/). 2,509/2,698 final records (5,949/6,562 pairs) already on exact MANE; ELA2 ×52 are on ELANE MANE (alias only). Off-MANE 137 records / 460 pairs:
   - Track C old version of same NM: 109 / 405 pairs, 10 genes; CDS identical (planner re-verified all 10) → offline relabel, 1 FOXP3 3′UTR query.
   - Track A on pick_mane fallback: 7 CDS/splice records worth normalizing (≤7 queries); 16 project deep intronic/UTR = Track A mapping errors.
   - No c.: C1QA ×6, CD79B, CD8A ×3 project 2–21 kb upstream yet tier A (precision blind spot: uncheckable in the audit); IGHG2 ×2 (to drop); one c.= record in final (source-error candidate, detail private).
   - pick_mane ROOT CAUSE (raw AK2 response, 99_Archive/ak2_raw.json): the API marks MANE as entry["tag"]["details"] == "MANE Select"; both pick_mane functions test "MANE_SELECT" in entry["selector"]["tags"], which never exists → always entries[0]. Fix folded into Run 12's code commit (byte-identical on cached data; cache stores parsed fields). User approved the recommended Run 12 scope (Track C relabel + FOXP3 + 7 Track A CDS records).
   - Forecast: ≤8 Mutalyzer, 0 VV, 6 GenBank (CARD9, DNMT3B, PTPN11, NRAS, CFD, FASLG MANE NMs), + 1 raw query.
   - Order: audit option A → Run 11 (remap reaches Track C/VV literal records; must precede the relabel, RFXAP R0006/R0007 is Track C) → Run 12 MANE normalization → Run 13 Track A mapping errors (26 records: demote tier / drop; audit should treat upstream/intronic Track A as a flag).

## README update (living draft)
- docs/agents/README_update_draft.md: post-thesis correctness section for the root README.md. After each recorded run the planner adds the run row and any new finding (figures from audit_history.tsv). Merge into README.md only when the user calls the stopping point; at merge also fix the outdated ADA −95 limitation and the thesis-era vs source-unit counting note.

## Backlog
- Track A offset errors (C8B, NHEJ1, C1QB, C5, IL7R, C2, CASP8, RAG2, CASP10, IGLL1): ~21 pairs left after R2.
- Gene symbols: current HGNC (ELA2 → ELANE; check all). IGHG2 reference decision.
- VariantValidator rescue route quality (12/14 checkable Run 1 gains protein-bad); vv_client.py aborts offline on its startup ping.
- Residual 8 excess pairs.
- Bug 2: per-allele exclusion in filter_bed (AIRE A0006/A0008).
- NMD: validate 50-nt-rule calls against aenmd (offline R) or VEP NMD plugin; Mutalyzer fs*N as stop source of truth; GenBank fetches in the forecast.
- Track C (NM_IDRefseq): KEEP — only track for 10 genes, 403 units at ~99% consistency.
- Infrastructure: scripts as rule inputs; audit_pairs.tsv in the rescue audit rule.
- VariO frameshift term; remove hard-coded "131 IDbases".
- Export key: _dedup_key is an internal record id; after Run 12 accepted Track A records keep their old token there (e.g. LYST c.2623del) while VariantOnTranscript/DNA shows the MANE c. The export must display VariantOnTranscript/DNA, never _dedup_key.
- Export allele: final tables carry no allele (dedup keeps the representative's allele_num; build_lovd_flat drops it). Take per-pair allele from augmented_merged_variants via source_change_id.
- Export fields: Screening (?/DNA), Genetic_origin, Effect, r.(?), Published as (source bracket notation), allele field from phase_group + phase-known (cis/trans/unknown; LOVD stores in-cis changes as separate variants linked by allele), individual_id, alt_descriptions, confidence_tier, OMIM.

## User decisions (2026-10-12)
- Accessions are unique only within a gene database (data agrees: 1,865 of 3,516 accessions occur in more than one gene, e.g. A0052 in ADA, AICDA, AIRE, CYBB, ITGB2, PTPN11, SH2). GENE_accession stands; dropped from the Mauno questions.
- Review-tier records go out as an xlsx review workbook, like the thesis-era lovd_review(_final).xlsx.
- IGHG2: the heavy-chain reference Mauno brought over from the BTKbase migration; he said dropping it is OK. User confirms IGHG2 and SH2 with Mauno this week.

## Mauno feedback 1 (2026-10-13): NMD rule
- Rule: if the PTC lies more than 50 nt before the end of the penultimate CODING exon → annotate predicted missing RNA and protein (r.0?, p.0?). Use protein-coding exons only (ignore 3′UTR-only last exons). Mutalyzer/VV frameshift p. is wrong for these.
- Before: Step9 measured from the transcript's last exon (7 MANE transcripts end in a 3′UTR-only exon: CD3G, CIITA, CD19, IL12B, NRAS, PTPN11, RAC2); the export kept the Mutalyzer fs p.; the suspect list ignored NMD (7 NMD-predicted candidate rows; detail private).
- User approved unfreezing Step9 for this rule. Plan: Step9 last-coding-exon boundary → export r.0?/p.0? (Mutalyzer p. to alt_descriptions) → suspect list drops NMD rows → nmd_review_for_mauno.csv (all PTC variants, per distinct variant). The audit keeps comparing the Mutalyzer protein (it tests the c. mapping).
- Open question for Mauno: NMD-escaping PTCs keep the Mutalyzer fs/nonsense p.?
- Done: Step9 rule df09bd4 (changes no call in current data; only 7 distances in CD19/CD3G/CIITA/IL12B). NMD CSV 5e81713: Output/Step11_Export/nmd_review_for_mauno.csv, 1,091 distinct PTC variants (predicted 872, escape last coding exon 186, escape within 50 nt 33; borderline 2; start-proximal 91), built from scratch Step9; tier/review columns being added. Pending: user's `-R annotate_nmd` run, then export r.0?/p.0? + workbook without the NMD suspect rows.

- NMD done (2026-10-13): Step9 outputs 6111487 (real run = scratch, 0 fetches); full CSV unchanged (f0cc89e); stratified sample 16f5f30 (nmd_sample_for_mauno.csv, 20 variants, seed 20261013; JAK3 example, borderline TNFRSF13B 52 nt / ITGB2 48 nt, CXCR4, IFNGR1); export 83e85f4: 2,378 NMD_predicted pairs → r.0?/p.0?, Mutalyzer r./p. in alt_descriptions, nmd_rule_applied column (45 columns); workbook suspect list 27 → 21. Reply draft: docs/agents/mauno_reply_nmd_draft.md. User view: all PTC → p.0? was considered; kept Mauno's rule as stated, pending his answer on escapes.

## Note for a future pipeline optimization (user, 2026-10-12)
- The lovd_* builders (Step7 build_lovd_flat.py, build_lovd_review_xlsx.py; Step8 build_final_lovd.py) sit mid-pipeline only because the LOVD table was the end of the thesis pipeline; post-thesis steps (VariO, patient metadata, NMD, audit) were added after it. If the pipeline is ever reworked: make the audited final variant table the end of the core pipeline and keep every export (LOVD tables, review workbooks, the canonical export) at the end, in Step11_Export. Consider reusing build_lovd_review_xlsx.py's writer there.
- Step11_Export stays as is (user decision); README 9d3a90e stands.

## Questions for Mauno
- Submission format from the LOVD director (pending).
- Accession IDs per gene database (confirms GENE_accession)?
- Screening placeholder, Effect "+", merged reference fields, NMD 50-nt rule, Review-tier handling, non-MANE isoform variants (CTSC), IGHG2 reference.
- Final numbers differ from the thesis (relevant only if published/deposited).
