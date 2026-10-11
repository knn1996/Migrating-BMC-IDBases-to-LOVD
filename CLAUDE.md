# IDbases → LOVD 3.0 migration

Migrating the IDbases immunodeficiency variant databases (MUTbase) to LOVD 3.0: GRCh38, MANE Select, HGVS. The goal is a correct LOVD submission, so precision (no wrong records) matters more than recovery.

## Repo
- Root: `/mnt/d/Assignment/Thesis` (Windows `D:\Assignment\Thesis`), run under WSL2 Ubuntu.
- Pipeline: `04_Mutation_Processing/Scripts`
  - `Snakefile`: primary workflow (Step1a extraction → Step7 merge/dedup)
  - `rescue.smk`: Step8 rescue + track resolution, Step9 NMD, Step10 audit
  - Config: `Scripts/config/config.yaml` (gitignored)
- Conda env `snakemake` (Snakemake 9.27). Logs in `04_Mutation_Processing/Logs`.
- Caches in `Output/cache/` (Mutalyzer, VariantValidator, GenBank) are append-only: never delete or rewrite entries.
- `99_Archive/` is gitignored and untracked.

## Hard rules
- Stage by explicit path. Never `git add -A` or `git add .`.
- No network git (fetch, pull, push). The user pushes.
- No live API calls (Mutalyzer, VariantValidator, NCBI) unless the user explicitly asks. Work offline (`unshare -rn` where possible).
- Before any run that may query an API, forecast uncached queries with `Utilities/check_mutalyzer_cache.py` and `Utilities/check_rescue_cache.py`, then stop for review.
- Do not run Snakemake for real unless asked; dry runs (`-n`) are fine. The user runs the pipeline in tmux.
- Offline scripts (audit, ledger, `--rescore-all`, cache forecasts, scratch verifications) you may run yourself. Real Snakemake runs and live API calls are the user's.
- Always use `--cores 1` in commands you give the user.
- Code: no comments. Default to targeted edits. Do a full rewrite only when the bug is structural (wrong key, data model or data flow), the change touches about a third of the script or more, or the script has already needed repeated patches. Say which you chose and why. Either way, verify in scratch that unaffected outputs are byte-identical.
- Compute figures from data; never hard-code counts.
- Present a plan table and get approval before reorganizing files.
- If a verification differs from the expected numbers, stop and report. Do not commit.
- Stop at the stop point the task gives.

## Pipeline mechanics
- Snakemake does not notice changes to Python scripts in shell rules. Force reruns with `-R <rule>`.
- Run order after a fix:
  1. `snakemake -s Snakefile --use-conda --cores 1 [-R <rule>]`
  2. Commit the changed primary outputs (Step6/Step7 and any other tracked primary outputs, e.g. BED and RefCheck files).
  3. `snakemake -s rescue.smk --use-conda --cores 1`
  4. Commit the rescue outputs.
  5. Record the ledger run, commit `history/`.
- `build_final_lovd` refuses to run until the primary dedup output is committed.
- Commit messages: `Run N primary: <fix>`, `Run N rescue: <fix>`, `Audit history: add Run N (<fix>)`.

## Audit ledger
- `Scripts/Step10_Audit/audit_ledger.py --label "<fix>"`; history in `Output/Step10_Audit/history/` (tracked).
- One fix per ledger run. Record only against committed, clean outputs.
- Ruler: source units parsed from pub.html (7,047; 7,019 in scope). Key = (gene, accession, allele). Never key on remapped HGVS.
- Three axes: recovery (`in_final_pct`), best-record accuracy (`protein_good_pct`), precision (`excess_pairs`, `pair_consistency_pct`).
- If the audit script changes, use `--rescore-all`. Every pre-existing column must reproduce exactly.
- After each run, report: previous vs current table including the precision columns, gains, regressions (list all), GOOD→BAD, BAD→GOOD, pair transitions.

## Data conventions
- `individual_id = GENE_accession`, keeping the letter prefix (e.g. `CYBB_A1416`). Accessions are only unique per gene database. The raw accession column stays; the audit matches on gene + raw accession.
- `source_change_id` links every track record to its source change.
- Track resolution keeps one primary record per source change (protein consistency → c.-region match → MANE Select → gene prior → track priority). Other descriptions of the same genomic change go to `alt_descriptions`. Every final record has `confidence_tier` A, B or Review.
- In-cis alleles (`c.[A;B]`) legitimately produce several records for one source unit. Never collapse two NM records of the same source unit.

## Reports
The user reviews reports in a separate planning chat. Lead with the result, then the numbers, then decisions needed. Keep it brief.
