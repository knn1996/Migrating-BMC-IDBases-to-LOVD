#!/usr/bin/env bash
set -euo pipefail
cd "${THESIS_ROOT:-/mnt/d/Assignment/Thesis}"

S="04_Mutation_Processing/Scripts"
FILES=(
  "$S/Snakefile"
  "$S/rescue.smk"
  "$S/Step1a_Extraction/extract_mutations.py"
  "$S/Step7_Merging/build_lovd_flat.py"
  "$S/Step7_Merging/recover_eintronic.py"
  "$S/Step7_Merging/resolve_unresolved.py"
  "$S/Step7_Merging/vv_client.py"
  "$S/Step7_Merging/annotate_vario.py"
  "$S/Step7_Merging/hgvs_norm.py"
  "$S/Step9_Audit/audit_against_source.py"
  "$S/Step9_NMD/annotate_nmd.py"
  "$S/refs/VariO_1.09.owx"
  "$S/R for figures/fig01_anchor_match_rates_slide.R"
  "$S/Utilities/reorganize_thesis.sh"
  "$S/Utilities/git_snapshot.sh"
)

missing=0
for f in "${FILES[@]}"; do
  if [ ! -e "$f" ]; then echo "MISSING: $f"; missing=1; fi
done
[ $missing -eq 0 ] || { echo "Aborting: fix missing paths first."; exit 1; }

for f in "${FILES[@]}"; do git add -- "$f"; done

echo "Staged for commit:"
git diff --cached --name-status

git commit -m "Snapshot before reorg: extraction fix, VariO, rescue helpers, NMD and audit steps"
echo
git log --oneline -1
