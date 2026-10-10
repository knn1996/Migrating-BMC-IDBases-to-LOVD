#!/usr/bin/env bash
set -uo pipefail

ROOT="${THESIS_ROOT:-/mnt/d/Assignment/Thesis}"
APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

P="04_Mutation_Processing"
S="$P/Scripts"
O="$P/Output"
A="99_Archive"
LOGDIR="$P/Logs"
STAMP="$(date +%Y%m%d_%H%M%S)"
MANIFEST="$ROOT/$LOGDIR/reorganize_${STAMP}.log"

cd "$ROOT" || { echo "Cannot cd to $ROOT"; exit 1; }
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || { echo "Not a git repo: $ROOT"; exit 1; }

n_ok=0; n_skip=0; n_warn=0

say() { echo "$*"; [ $APPLY -eq 1 ] && echo "$*" >> "$MANIFEST"; }

is_tracked() { [ -n "$(git ls-files -- "$1" 2>/dev/null | head -n1)" ]; }

do_move() {
    local src="$1" dst="$2"
    if [ $APPLY -eq 0 ]; then
        say "  [dry] $( is_tracked "$src" && echo 'git mv' || echo 'mv    ' )  $src  ->  $dst"
        return 0
    fi
    mkdir -p "$(dirname "$dst")"
    if is_tracked "$src"; then
        git mv -- "$src" "$dst" 2>/dev/null || mv -- "$src" "$dst"
    else
        mv -- "$src" "$dst"
    fi
}

move() {
    local src="$1" dst="$2"
    if [ ! -e "$src" ]; then say "  SKIP (missing)  $src"; n_skip=$((n_skip+1)); return; fi
    if [ -e "$dst" ]; then say "  WARN (target exists, not moved)  $dst"; n_warn=$((n_warn+1)); return; fi
    do_move "$src" "$dst"
    [ $APPLY -eq 1 ] && say "  moved  $src  ->  $dst"
    n_ok=$((n_ok+1))
}

case_rename() {
    local src="$1" dst="$2" tmp="$1.__tmp_rename__"
    if [ ! -e "$src" ]; then say "  SKIP (missing)  $src"; n_skip=$((n_skip+1)); return; fi
    if [ $APPLY -eq 0 ]; then say "  [dry] rename (via temp)  $src  ->  $dst"; n_ok=$((n_ok+1)); return; fi
    do_move "$src" "$tmp" && do_move "$tmp" "$dst"
    say "  renamed  $src  ->  $dst"
    n_ok=$((n_ok+1))
}

archive() { move "$1" "$A/$1"; }

copy_new() {
    local src="$1" dst="$2"
    if [ ! -e "$src" ]; then say "  SKIP (missing)  $src"; n_skip=$((n_skip+1)); return; fi
    if [ -e "$dst" ]; then say "  WARN (target exists, not copied)  $dst"; n_warn=$((n_warn+1)); return; fi
    if [ $APPLY -eq 0 ]; then say "  [dry] cp  $src  ->  $dst"; n_ok=$((n_ok+1)); return; fi
    mkdir -p "$(dirname "$dst")"; cp -p -- "$src" "$dst"
    say "  copied  $src  ->  $dst"; n_ok=$((n_ok+1))
}

if [ $APPLY -eq 1 ]; then
    mkdir -p "$LOGDIR"
    echo "reorganize_thesis.sh  applied $STAMP  root=$ROOT" > "$MANIFEST"
else
    echo "DRY RUN — nothing will be changed. Re-run with --apply to execute."
fi

say ""
say "== 1. Mutalyzer / VariantValidator caches -> names Snakemake expects =="
move "$O/Step6_Mutalyzer/mutalyzer_cache.jsonl"             "$O/cache/mutalyzer/track_a_g.jsonl"
move "$O/Step6_Mutalyzer/mutalyzer_cache_c.jsonl"           "$O/cache/mutalyzer/track_a_c.jsonl"
move "$O/Step6_Mutalyzer/mutalyzer_cache_NM_MANE.jsonl"     "$O/cache/mutalyzer/track_b_MANE.jsonl"
move "$O/Step6_Mutalyzer/mutalyzer_cache_NM_IDRefseq.jsonl" "$O/cache/mutalyzer/track_c_IDRefseq.jsonl"
copy_new "$O/cache/variantvalidator/vv_rescue.jsonl"        "$O/cache/variantvalidator/vv_rescue_probe.jsonl"
move "$O/cache/variantvalidator/vv_rescue.jsonl"            "$O/cache/variantvalidator/vv_resolve.jsonl"

say ""
say ""
say "== 2. Single log folder: Output/logs -> $LOGDIR =="
if [ -d "$O/logs" ]; then
    for f in "$O/logs"/*; do
        [ -e "$f" ] || continue
        move "$f" "$LOGDIR/$(basename "$f")"
    done
    [ $APPLY -eq 1 ] && rmdir "$O/logs" 2>/dev/null && say "  removed empty $O/logs"
else
    say "  SKIP (missing)  $O/logs"
fi

say ""
say "== 3. Scripts: Step numbering and rescue scripts =="
move "$S/Step9_Audit" "$S/Step10_Audit"
for f in classify_unresolved.py generate_rescue_input.py vv_client.py resolve_unresolved.py \
         merge_with_rescued.py hgvs_norm.py recover_eintronic.py; do
    move "$S/Step7_Merging/$f" "$S/Step8_Rescue/$f"
done
for f in check_missing.py count_patients_variants.py reconcile_source_to_output.py; do
    move "$S/$f" "$S/Utilities/$f"
done

say ""
say "== 4. Archive: superseded scripts and junk in Scripts =="
archive "$S/Step1a_Extraction/extract_mutations_fixed.py"
archive "$S/Step1a_Extraction/extract_mutations_old.py"
archive "$S/Scripts.rar"
archive "$S/results"
for d in "$S"/C*Users*BornLoser*; do
    [ -e "$d" ] && archive "$d"
done

say ""
say "== 5. Archive: Output clutter =="
for f in all_mutations.tsv all_mutations_v0.tsv lrg_offset_results.csv debug.txt \
         "for low matching rate cases.txt" Mutation_extraction_refined.pdf strand.csv \
         Thesis_Number_Inventory.xlsx; do
    archive "$O/$f"
done
for f in all_mutations_long_seq.tsv all_mutations_no_subs.tsv all_mutations_pre_fix.tsv \
         all_mutations_v0.tsv all_mutations_v1.tsv batch_mutation_extraction_summary.csv \
         cross_check_mutations.tsv cross_check_rna_mutations.tsv; do
    archive "$O/Step1_Extraction/$f"
done
for f in mutalyzer_cache_NM.jsonl mutalyzer_results_NM.tsv mutalyzer_results_rc.tsv; do
    archive "$O/Step6_Mutalyzer/$f"
done
archive "$O/Step7_Merging/back_up"

say ""
say "== 6. Output step folders aligned to Scripts names =="
move        "$O/Step1_Extraction" "$O/Step1a_Extraction"
case_rename "$O/Step5_Liftover"   "$O/Step5_LiftOver"

say ""
say "== 7. Archive: legacy pre-Snakemake material =="
archive "$P/BLAST_Results"
archive "$P/Mutalyzer"
archive "$P/Mutation Extraction and Matching"
archive "$P/blast_unconfirmed.sh"
archive "03_BED_Files/BED_source_1"
archive "03_BED_Files/BED_LRG_source_3"
archive "protein_only_variants.csv"
archive "blast-first-hits.skill"

say ""
say "== 8. .gitignore: ignore the archive =="
if grep -qxF "99_Archive/" .gitignore; then
    say "  SKIP (already present)  99_Archive/ in .gitignore"
elif [ $APPLY -eq 1 ]; then
    printf '\n# ---- Local archive of superseded / legacy material ----\n99_Archive/\n' >> .gitignore
    say "  added 99_Archive/ to .gitignore"
else
    say "  [dry] append 99_Archive/ to .gitignore"
fi

say ""
say "Summary: $n_ok actions, $n_skip skipped (missing), $n_warn warnings (target exists)"
[ $APPLY -eq 1 ] && say "Manifest: $MANIFEST" && say "Review with: git status --short | head -50"
