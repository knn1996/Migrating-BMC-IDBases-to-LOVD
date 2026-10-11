# Handoff — continue IDbases → LOVD pipeline work (updated 2026-10-11 21:15)

## Context
MSc thesis (Lund, supervisor Mauno Vihinen) is DONE. Current goal: finish the migration of IDbases from MUTbase to LOVD 3.0 (GRCh38, MANE Select, HGVS) so the data can be submitted. Priority is correctness of what goes into LOVD (precision), not thesis figures. Repo: D:\Assignment\Thesis = /mnt/d/Assignment/Thesis in WSL2 Ubuntu. Hands-on work runs through Claude Code in WSL; the chat is for planning and review.

Claude Code: launch from the repo root (`cd /mnt/d/Assignment/Thesis && claude`). CLAUDE.md (working rules, run order, ledger and data conventions) drafted for the repo root; commit it on its own. Use Opus for investigative tasks; Sonnet (`/model`) for mechanical chores.

Shell setup: `tmux new -A -s thesis`, then `thesis` (bash function in ~/.bashrc: cd to Scripts + conda activate snakemake).

Pipeline = Snakemake (conda env `snakemake`, Snakemake 9.27), run from 04_Mutation_Processing/Scripts:
- `Snakefile`: primary workflow (Step1a extraction → … → Step7 merge/dedup)
- `rescue.smk`: Step8 rescue (mutalyzer_rescue, vv_rescue, annotate_vario, resolve_tracks, build_final_lovd), Step9 NMD, Step10 audit
- Run pattern: `snakemake [-s rescue.smk] --use-conda --cores 1 2>&1 | tee ../Logs/<name>_$(date +%Y%m%d_%H%M).log`. Always `--cores 1`.
- Snakemake does NOT detect Python script changes: force with `-R <rule>`.
- Commit primary outputs (Step6/Step7 plus any tracked BED/RefCheck files) between the primary and rescue runs; build_final_lovd refuses otherwise.
- audit_pairs.tsv is written by the ledger, not the rescue audit rule: Claude Code generates it with a throwaway ledger run before committing rescue outputs (backlog: move pair labelling into audit_against_source.py).

## Audit ledger
- Scripts/Step10_Audit/audit_ledger.py (+ --rescore-all), Scripts/Utilities/plot_audit_history.py. History in Output/Step10_Audit/history/ (tracked).
- Key = (gene, accession, allele) from pub.html. Ruler: 7,047 source units, 7,019 in scope.
- Axes: recovery (in_final %), best-record accuracy (protein_good %), precision (excess_pairs, pair_consistency_pct over every final pair).

| Run | Commit | Label | In final | Protein good | Pair consistency | Excess pairs | Variants / individuals / pairs |
|---|---|---|---|---|---|---|---|
| 0 | 7c7f6fc | pre-fix baseline | 83.33% | 90.79% | 86.56% | 654 | 2,651 / 5,326 / 6,589 |
| 1 | 0641835 | extraction fix | 89.84% | 92.74% | 84.58% | 1,068 | 3,003 / 5,619 / 7,508 |
| 2 | d8abc09 | dedup fix + individual_id | 90.03% | 93.06% | 85.58% | 988 | 3,018 / 5,619 / 7,428 |
| 3 | c28efe5 | bed_to_mutalyzer delins/dup fix | 89.98% | 93.11% | 88.94% | 734 | 2,967 / 5,616 / 7,171 |
| 4 | fb548f7 | IVS position fix | 89.98% | 93.11% | 88.94% | 732 | 2,966 / 5,616 / 7,169 |
| 5 | f51bc07 | R2 track resolution | 89.98% | 96.20% | 95.38% | 24 | 2,663 / – / 6,456 |

Run 5 detail: inconsistent pairs 680 → 255 (DIFFERENT 175, FS_START 6, FS_STOP 6, POS_SHIFT 13, TYPE_DIFF 55); uncheckable 931; 161 BAD→GOOD; 0 gains/regressions. Confidence tiers: A 6,208, B 199, Review 49 (Step8_Rescue/resolution_log.tsv, tier = Review).

Current precision estimate of the final build: ~95% (range set by 931 uncheckable pairs).

## Done (local commits, NOT pushed; origin at 84f25ba)
- Runs 0–3: see table; 010f8eb, 6380e17, 0641835, 9005d31, dd183da, d8abc09, c8cb1f2, 2014d12, b400dda, c28efe5, 396eb1d.
- Precision metrics + --rescore-all: a0d5cc4 (code), 7968f8a (history). Validation: all old columns reproduced exactly; Run 2→3 excess fell by exactly 254.
- Run 4 (bug 3, IVS position inheritance; option A = leave position empty): da1f05f code, c4ef82d primary (amended to include 12 BED/RefCheck files; hash may have changed), fb548f7 rescue, a36f153 history.
- Bug 5: 581b4cb — allele_num and source_change_id through merge_tracks (Track A rows had no allele; no output change).
- Run 5 (R2 track resolution): resolve_tracks rule after dedup_rescued; winner by protein consistency → c.-region match → MANE Select → gene prior → track priority; same-genomic-change alternates kept in alt_descriptions; confidence_tier column. Code commit (hash in Claude Code report), Run 5 primary commit, f51bc07 rescue, b7b1f59 history.
- Root cause found: merge_tracks keyed on sysname (genomic for Track A, coding for others), so track priority never fired and nearly every source change got two final records.

## Next (in order)
1. Report-only characterization of the 255 inconsistent pairs + 49 Review: group by gene/BAD_P label; constant per-gene protein offset (numbering convention → reclassify) vs pipeline fix (legacy c. numbering, Track A offset; forecast queries) vs likely source error. Prompt drafted in chat; use Opus.
2. Run 6 based on that report.
3. Patient metadata: fix extract_patient_metadata.py glued-record bug (metadata flags 3,141).
4. LOVD import builder.

## Backlog
- Root causes behind resolved conflicts: Track A offsets (C8B 68 groups, NHEJ1), legacy c. numbering applied literally by the MANE track and VariantValidator (AIRE, CD40L, PRF1, NCF1, CFP, F1, LYST), transcript choice (CTSC, HAX1, CXCR4 — source change in an exon absent from MANE).
- Residual 24 excess pairs: 16 with no source-change mapping, 7 individuals whose expected count differs.
- Uncheckable pairs (931): non-protein check via c. notation comparison within the CDS / intronic offsets.
- NMD validation: cross-check 50-nt-rule calls on ~1,326 PTC variants against aenmd (R, offline) or VEP NMD plugin/LOFTEE; make Mutalyzer's p.(…fs*N) the frameshift-stop source of truth (53 review rows disagree).
- Bug 2 (small): per-allele exclusion in filter_bed (AIRE A0006/A0008).
- Track A RAG1/RAG2: source writes cDNA numbers in the genomic field.
- Add each rule's script as a rule input; move audit_pairs.tsv generation into the rescue audit rule; vv_client.py startup ping.
- VariO frameshift term; VariO in xlsx; remove hard-coded "131 IDbases" in build_lovd_review_xlsx.py.
- LOVD import builder: Screening (?/DNA), Genetic_origin, Effect, r.(?), Published as (source bracket notation), allele field from phase_group + phase-known (cis/trans/unknown), individual_id, alt_descriptions as extra transcript annotations, confidence_tier, OMIM linking, test import.

## Questions for Mauno
- Final migrated numbers differ from the thesis (extraction, dedup, NG delins/dup, track-resolution fixes) — relevant if published or deposited.
- Are accession IDs per gene database (confirms GENE_accession)?
- Which LOVD instance; Screening placeholder; Effect "+"; merged reference fields; NM_MANE offset issue; NMD 50-nt rule; how to handle the Review tier and variants on non-MANE isoforms (CTSC).

## Working rules
Brief, logic-focused answers; no code comments; ask when direction is unclear; I run code myself; full-script rewrites; never `git add -A`; plan table + approval before reorganizing; no network git through the bridge; forecast before any API run; one fix per ledger run; `--cores 1`.
