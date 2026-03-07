#!/usr/bin/env bash
# summarize_blast_first_hits.sh
# ==============================
# Extracts the first hit from each BLAST outfmt7 result file matching
#   *_vs_hg16*, *_vs_hg17*, *_vs_hg18*
# and writes a tab-separated summary CSV.
#
# Usage
# -----
#   bash summarize_blast_first_hits.sh [RESULTS_DIR] [-o OUTPUT_FILE]
#
# Examples
#   bash summarize_blast_first_hits.sh
#   bash summarize_blast_first_hits.sh /home/khoi1996/Documents/BLAST_Results
#   bash summarize_blast_first_hits.sh /path/to/results -o /path/to/summary.tsv

# ── Defaults ──────────────────────────────────────────────────────────────────
RESULTS_DIR="/home/khoi1996/Documents/BLAST_Results"
OUTPUT_FILE=""

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output) OUTPUT_FILE="$2"; shift 2 ;;
        -*) echo "Unknown option: $1"; exit 1 ;;
        *)  RESULTS_DIR="$1"; shift ;;
    esac
done

OUTPUT_FILE="${OUTPUT_FILE:-${RESULTS_DIR}/blast_first_hits_summary.tsv}"

# ── Sanity check ──────────────────────────────────────────────────────────────
if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "ERROR: Directory not found: $RESULTS_DIR"
    exit 1
fi

# ── Write header ──────────────────────────────────────────────────────────────
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "gene" "build" "file" "total_hits" \
    "query_id" "subject_id" "pct_identity" "alignment_length" \
    "mismatches" "gap_opens" \
    "q_start" "q_end" "s_start" "s_end" \
    "evalue" \
    > "$OUTPUT_FILE"

# ── Process each matching file ────────────────────────────────────────────────
count=0

for filepath in "$RESULTS_DIR"/*_vs_hg1[678]*; do
    [[ -f "$filepath" ]] || continue

    filename=$(basename "$filepath")

    # Extract build (hg16 / hg17 / hg18)
    if [[ "$filename" =~ (hg1[678]) ]]; then
        build="${BASH_REMATCH[1]}"
    else
        build="unknown"
    fi

    # Extract gene name: everything before _vs_hg1X
    gene="${filename%%_vs_${build}*}"

    # Total hits from header comment
    total_hits=$(grep -m1 "hits found" "$filepath" \
                 | grep -oP '[\d,]+(?=\s+hits found)' \
                 | tr -d ',')
    total_hits="${total_hits:-0}"

    # First non-comment, non-empty data line
    first_hit=$(grep -v '^#' "$filepath" | grep -m1 '.')

    if [[ -z "$first_hit" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$gene" "$build" "$filename" "$total_hits" \
            "NO_HITS" "" "" "" "" "" "" "" "" "" "" \
            >> "$OUTPUT_FILE"
    else
        # Parse the 12 tab-separated outfmt7 fields
        IFS=$'\t' read -r qid sid pct alen mis gap qs qe ss se ev bs <<< "$first_hit"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$gene" "$build" "$filename" "$total_hits" \
            "$qid" "$sid" "$pct" "$alen" "$mis" "$gap" \
            "$qs" "$qe" "$ss" "$se" "$ev" \
            >> "$OUTPUT_FILE"
    fi

    echo "  Processed: $filename  (hits: $total_hits)"
    (( count++ ))
done

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "Done. $count file(s) processed."
echo "Summary written to: $OUTPUT_FILE"