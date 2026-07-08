# Release notes — v1.1.1

**Type:** documentation, portability, and reproducibility patch.
**No scientific logic, DAG shape, or reference-selection decisions were changed.**

v1.1.1 builds on the v1.1.0 config-driven Snakemake refactor. It makes the
repository portable and reproducible for anyone cloning it, and replaces the stale
"edit the CONSTANTS block" instructions with the real Snakemake workflow.

## Highlights

- **README Quickstart rewritten** around the actual Snakemake workflow
  (`cp config/config.example.yaml config/config.yaml` → `snakemake --use-conda
  --cores 8 -n` → real run). Final output documented as
  `results/lovd_review.xlsx`, with the intermediate layout under `results/`.
- **Portable config layout.** Personal absolute Windows paths (which exposed a
  username) are removed from version control. The pipeline now uses
  `config/config.example.yaml` (committed, neutral placeholder paths) plus
  `config/config.yaml` (gitignored, your real paths). The Snakefile `configfile:`
  directive points at `config/config.yaml` and resolves relative to the Snakemake
  working directory — always launch from `04_Mutation_Processing/Scripts`.
- **"Required external data" table** added to the README, mapping each third-party
  or large input (IDbases dirs, hg16/17/18 FASTA, MANE RNA FASTA,
  `LRG_RefSeqGene.txt`, alias table, curated IDRefSeq sequences) to its source and
  config key.
- **`run_set` documented prominently.** The default config is a **3-gene smoke
  test** (`RAB27A`, `UNC13D`, `ADA`); set `run_set: all` for the full thesis
  dataset (auto-discovered from the IDbases directories).
- **`.gitignore` fixed.** The concatenated `Mutation_extraction.ipynb` /
  `Scripts.rar` line is split into two separate rules; both are confirmed ignored.
- **`CITATION.cff` bumped** to 1.1.1 with an updated release date and the Zenodo
  version DOI. Thesis URL/DOI remain `[TBA]`.
- **Conda envs pinned to minor versions** (`pandas=2.2`, `numpy=2.0`, `scipy=1.14`,
  `openpyxl=3.1`, `requests=2.32`, `blast=2.16`, `ucsc-liftover=469`;
  `python=3.11` unchanged). Full lock files are a later improvement.
- **GitHub Actions CI added** (`.github/workflows/ci.yml`) that runs
  `snakemake --lint` and `snakemake -n` against a stub-fixture config
  (`tests/config.ci.yaml` + `tests/make_fixtures.sh`). The full scientific run
  cannot execute on CI (large reference data + Mutalyzer network).

## Validation status — dry-run + lint only

A full 3-gene execution was **not** performed for this release: the environment
lacks the reference genomes, the IDbases source data, the curated IDRefSeq
sequences, and Mutalyzer network access. Validation was therefore limited to DAG
resolution and linting.

What was verified:

- **`snakemake -n` resolves cleanly** with the stub-fixture config and the default
  `run_set: test` — **33 core jobs**, matching the v1.1.0 dry-run
  (`blast_gene` ×9, `make_blast_db` ×3, all others ×1), terminating at
  `results/lovd_review.xlsx`. This confirms the `config/` move and the
  `configfile:` directive change did not break the workflow.
- **`snakemake --lint` runs.** Remaining warnings are non-blocking and deferred:
  28 × "no `log:` directive" (the wired scripts already emit their own logs via
  `LOG_PATH`/`LOG_DIR`), 2 × "hardcoded prefix param" (BLAST DB prefixes in
  `make_blast_db`/`blast_gene`), and 1 × "mixed rules and functions" (standard for
  a single-file Snakefile). CI reports lint but does not fail on it; the hard gate
  is that `snakemake -n` must resolve. Adding `log:` directives across all rules is
  a v2.0 workflow-hardening item.

Environment used for validation: Ubuntu 22.04.5 LTS, Snakemake 7.32.4, Python 3.11
(conda envs) — dry run under sandbox Python 3.10.

### Ready-to-run command for a real 3-gene validation

On a machine/HPC with the reference data in place and `config/config.yaml`
pointing at it:

```bash
cd 04_Mutation_Processing/Scripts
cp config/config.example.yaml config/config.yaml   # then edit paths (run_set: test)
snakemake --use-conda --cores 8 -n                 # expect 33 jobs
snakemake --use-conda --cores 8                     # produces results/lovd_review.xlsx
```

After a successful run, record: exact command, OS/Snakemake version, runtime,
completed job count, and row counts / checksums of the major outputs
(`all_mutations.tsv`, `merged_variants.tsv`, `lovd_flat_with_patients.tsv`,
`results/lovd_review.xlsx`). The release can then be re-labelled
"execution-validated (3-gene test)".

## Not in this release (deferred to v2.0)

- Converting scripts from `os.environ` to argparse/CLI.
- Replacing `directory()` rule outputs with per-file outputs.
- Adding `log:` directives to every rule / conda lock files.
- Any change to scientific logic, the DAG, or reference-selection decisions.
