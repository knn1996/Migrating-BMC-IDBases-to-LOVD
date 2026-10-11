# Plan — Audit ledger (original design, 2026-10-11; current state in HANDOFF.md)

## Goal
After each pipeline fix, score the final output against pub.html (the independent ground truth) and keep a history, so we can show how much each fix improved the pipeline and catch regressions.

## Principles
- Fixed ruler: the denominator is the set of source variant units parsed from pub.html by audit_against_source.py. It must stay constant between runs. If the audit script changes, re-score past runs with the new version.
- Coverage is always shown with accuracy: recovery % is never reported without accuracy/precision, so a gain can't come from loosening the gates.
- Stable variant key: (gene, accession, allele) from pub.html. Never the remapped HGVS, which can change between runs.
- Score only committed outputs: each ledger row records the pipeline commit, so any past run can be re-scored from git.

## Files
- `Scripts/Step10_Audit/audit_ledger.py` (args: --label, --audit-dir, --pipeline-commit, --run-id, --allow-dirty, --rescore-all, --allow-move)
- `Output/Step10_Audit/history/` (tracked): audit_history.tsv (one wide row per run; same run_id replaces), snapshots/<run_id>.tsv.gz and .pairs.tsv.gz, changes/<run_id>_vs_<prev>.tsv, audit_trend.png
- `Scripts/Utilities/plot_audit_history.py`

## Metrics per run
- run_id, date, label, pipeline_commit, audit_script_sha
- source_units, out_of_scope, in_final(_pct), known_unresolved, not_registered (+ split by fate_detail)
- protein_checkable / good / good_pct / bad (+ per BAD_P label)
- precision: individuals_with_excess, excess_pairs, pairs_checkable/consistent/inconsistent (+ per label), pair_consistency_pct, pairs_uncheckable, pair transitions
- c./p. coherence: coherent / incoherent / unknown, coherence_pct of known
- metadata_flags (excl. SUBHEADERS_STRIPPED, EXTRA_REFS_NOT_CARRIED)
- final_unique_variants, final_individuals (by individual_id), final_pairs
- vs previous run: gains, regressions, protein GOOD→BAD, BAD→GOOD

## Guardrails
- If source_units ≠ previous run's value: stop with "ruler changed — re-score history".
- --rescore-all: every pre-existing column must reproduce exactly; only new or explicitly allowed columns may move.
- Not a Snakemake rule (needs a human label). Run after each committed audit.
