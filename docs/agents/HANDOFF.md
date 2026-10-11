# Handoff — IDbases → LOVD pipeline (source of truth; updated 2026-10-12 01:10)

This file is the single source of truth. The planner (`thesis-plan`) updates it after each recorded run. The older copy in the claude.ai project (claude/handoff_2026-10-11.md) is frozen and may be stale.

## Context
MSc thesis (Lund, supervisor Mauno Vihinen) is DONE. Goal: finish migrating IDbases from MUTbase to LOVD 3.0 (GRCh38, MANE Select, HGVS) so the data can be submitted. Priority: correctness of what goes into LOVD (precision) over recovery. Repo: D:\Assignment\Thesis = /mnt/d/Assignment/Thesis in WSL2 Ubuntu.

Submission format: Mauno is contacting the current LOVD director; for a database this size LOVD likely defines a custom submission format. Do NOT build a generic LOVD import file. Do format-independent correctness first, then one canonical export table (one row per patient–variant pair with everything any format could need); the final format becomes a mapping job once known.

Sessions: executor `thesis-exec` (Claude Code, repo root, Remote Control on), planner `thesis-plan` (Claude Code, repo root, see docs/agents/PLANNER.md). Shared rules in CLAUDE.md.

Shell setup for the user: `tmux new -A -s thesis`, then `thesis` (bash function: cd to Scripts + conda activate snakemake).

Pipeline (Snakemake 9.27, env `snakemake`, run from 04_Mutation_Processing/Scripts):
- `Snakefile`: Step1a extraction → … → Step7 merge/dedup
- `rescue.smk`: Step8 rescue, resolve_tracks, build_final_lovd, VariO; Step9 NMD (fetches missing GenBank transcripts from NCBI); Step10 audit
- Always `--cores 1`. Force code changes with `-R <rule>`. Commit primary outputs before the rescue run (build_final_lovd refuses otherwise).
- audit_pairs.tsv is produced by the ledger (throwaway ledger run before committing rescue outputs).
- Forecast tools don't cover NMD GenBank fetches (backlog).
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
Coherence ruler: 3fc3b8a (code), 4c7c6b1 (history). Run 6: dcde390, 9e0e2b7, f6de49d, d82a5c5, 7543cbb (rescue), aa6c043 (history). Nothing pushed.

## Fixed so far
1. Extraction: compound-het alleles, glued records, brackets (Run 1)
2. Position-only XM/NM merge (Run 2)
3. individual_id = GENE_accession (Run 2)
4. bed_to_mutalyzer delins/dup/minus-strand (Run 3)
5. IVS positions inherited from neighbours (Run 4)
6. Track A allele labels + source_change_id (Run 5 prep)
7. Track priority never fired → R2 track resolution with alt_descriptions and confidence_tier A/B/Review (Run 5)
8. Transcript-blind stage-1 dedup, resolve_tracks follows stage-2 merges, trailing ";" in accessions (Run 6)
Measurement: count definitions, precision metric, own-source-change ruler, c./p. coherence axis, --rescore-all.

## Run 7 — CURRENT (Track B c./p. coherence)
- Cause: run_mutalyzer_NM.py call_normalize took p. from the normalized (MANE) transcript but c. from pick_mane(equivalent_descriptions["c"]), which falls back to the first entry (another transcript).
- Code edit DONE, UNCOMMITTED in the working tree: run_mutalyzer_NM.py emits c. = cached normalized when c_hgvs is on another transcript; old c. → new alt_c_hgvs column (own p. from the Track A c.-pass cache where available); passed through merge_tracks.py, merge_with_rescued.py, prepended into alt_descriptions by resolve_tracks.py; Snakefile passes ALT_P_CACHE. Cache reinterpreted, never rewritten.
- Offline simulation (NM tracks → audit; NMD not simulated): in_final 6,316, 0 regressions, excess 8, pair consistency 96.65%, inconsistent 183; coherence 73.69% → 100% of known (5,703 coherent, 0 incoherent, 738 unknown); Track B on exact current MANE 64/1,322 → 1,829/1,829; all final records on exact current MANE 1,206 → 2,464 of 2,663. 676 pairs renamed (c. string changed); stage-2 identical_change merges 449 → 4; skip_different_change_same_unit 23 → 74; R2 conflicts 515 unchanged, 198 gone, 199 new.
- Queries: 0 Mutalyzer, 0 VariantValidator; 59 NCBI GenBank fetches in the NMD step (new MANE transcripts, e.g. NM_000057.4, NM_001625.4, NM_003664.5).
- DECISION: user chose option 1 — the real rescue run fetches the 59 GenBank records (network only in the user's own run, then cached).
- Remaining steps: executor commits the code → user runs `snakemake -s Snakefile --use-conda --cores 1 -R mutalyzer_track_NM` → executor commits primary outputs → user runs rescue.smk (with the fetches) → executor commits rescue outputs and records Run 7 ("Track B c./p. coherence"). Expect: in_final 6,316, 0 regressions, excess 8, consistency 96.65%, coherence ~100% of known, real outputs = simulation except NMD.

## Next after Run 7
1. Run 8 — Legacy c. numbering remap (~70 pairs; RAG1 +112, AICDA +76, CFH +73, TCIRG1 +60, WAS 34–36, LYST 190, RFXAP 116) + ADA over-remap. Needs Utilities/translation_test.py (draft, uncommitted): forecast missing transcripts first (fetch = network step), per-gene shift verified when ≥2 independent changes reconcile (C1QB −2 codons, RFX5 −40 codons are isoform numbering), include small frameshifts/in-frame indels, one shared remap function for c. tracks and rescue. ~25 queries + transcripts. Also yields source_error flags.
2. Run 9 — MANE normalization of the ~200 records still off exact MANE after Run 7 (Track A on pick_mane fallback transcripts, Track C non-MANE NMs, VariantValidator/rescue). Recompute scope and forecast. First find why pick_mane never matches MANE_SELECT (needs 1 live raw Mutalyzer response; user runs it). Only IGHG2 has no MANE Select (ELA2 = ELANE NM_001972.4). Keep originals in alt_descriptions.
3. Patient metadata: extract_patient_metadata.py glued-record bug (metadata flags ~3,147).
4. Canonical export table; then the LOVD format mapping once the director's format is known.

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
- Export fields: Screening (?/DNA), Genetic_origin, Effect, r.(?), Published as (source bracket notation), allele field from phase_group + phase-known (cis/trans/unknown; LOVD stores in-cis changes as separate variants linked by allele), individual_id, alt_descriptions, confidence_tier, OMIM.

## Questions for Mauno
- Submission format from the LOVD director (pending).
- Accession IDs per gene database (confirms GENE_accession)?
- Screening placeholder, Effect "+", merged reference fields, NMD 50-nt rule, Review-tier handling, non-MANE isoform variants (CTSC), IGHG2 reference.
- Final numbers differ from the thesis (relevant only if published/deposited).
