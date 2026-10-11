# Handoff — IDbases → LOVD pipeline (source of truth; updated 2026-10-12, after Run 7)

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
Coherence ruler: 3fc3b8a (code), 4c7c6b1 (history). Run 6: dcde390, 9e0e2b7, f6de49d, d82a5c5, 7543cbb (rescue), aa6c043 (history). Run 7: 0c255fe (code), 1e12fd6 (primary), d6b09bb (rescue), 3b4cc34 (history, run_id 20261012_d6b09bb). Nothing pushed.

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
