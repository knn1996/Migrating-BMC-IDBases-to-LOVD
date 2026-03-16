#!/bin/bash
#SBATCH -A lu2025-2-88
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH -J BLAST_crosscheck
#SBATCH --array=0-60%20
#SBATCH --output=logs/%x_%j_%a.out
#SBATCH --error=logs/%x_%j_%a.err

# ==============================================================================
# BLAST Crosscheck: Processed FASTA (Source 1) vs LRG NM & NG (Source 3)
#
# Each gene's processed sequence is compared pairwise against its own NM and NG
# LRG counterpart. Uses -subject (no makeblastdb needed) since comparisons are
# gene-to-gene, not against a full genome database.
# ==============================================================================

module load GCC/13.3.0
module load OpenMPI/5.0.3
module load BLAST+/2.16.0

# ── Directory paths (update these when moving to HPC) ─────────────────────────
SOURCE1_DIR="/home/khoi1996/Documents/DNA_sequences/Processed_FASTA_Source1"
LRG_NM_DIR="/home/khoi1996/Documents/DNA_sequences/LRG_FASTA_Source3/NM"
LRG_NG_DIR="/home/khoi1996/Documents/DNA_sequences/LRG_FASTA_Source3/NG"
OUTPUT_DIR="/home/khoi1996/Documents/BLAST_Results/crosscheck"
# ──────────────────────────────────────────────────────────────────────────────

mkdir -p "$OUTPUT_DIR/logs"

echo "=============================================="
echo "BLAST Crosscheck Job"
echo "Array task ID : $SLURM_ARRAY_TASK_ID"
echo "Job ID        : $SLURM_JOB_ID"
echo "Node          : $(hostname)"
echo "Started       : $(date)"
echo "=============================================="

# Load all Processed FASTA files into an indexed array
FASTA_FILES=("$SOURCE1_DIR"/*.FASTA)
TOTAL=${#FASTA_FILES[@]}

echo "Total Source 1 FASTA files found : $TOTAL"

# Bounds check
if [[ $SLURM_ARRAY_TASK_ID -ge $TOTAL ]]; then
    echo "ERROR: Task ID $SLURM_ARRAY_TASK_ID exceeds file count ($TOTAL). Exiting."
    exit 1
fi

query_file="${FASTA_FILES[$SLURM_ARRAY_TASK_ID]}"

if [[ ! -f "$query_file" ]]; then
    echo "ERROR: Query file not found: $query_file"
    exit 1
fi

# Derive gene name from filename (strip path and .FASTA extension)
gene_name=$(basename "$query_file" .FASTA)

echo ""
echo "Processing gene  : $gene_name"
echo "Query file       : $query_file"
echo ""

# ── Common BLAST parameters ───────────────────────────────────────────────────
BLAST_OUTFMT="7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
BLAST_EVALUE="1e-10"
BLAST_THREADS=8
BLAST_MAX_HITS=100
# ──────────────────────────────────────────────────────────────────────────────

run_blast() {
    local label="$1"
    local subject="$2"
    local outfile="$3"

    echo "  Running BLAST : $gene_name vs $label"
    echo "  Subject       : $subject"
    echo "  Output        : $outfile"

    timeout 4h blastn \
        -query   "$query_file" \
        -subject "$subject" \
        -out     "$outfile" \
        -outfmt  "$BLAST_OUTFMT" \
        -num_threads "$BLAST_THREADS" \
        -evalue  "$BLAST_EVALUE" \
        -max_target_seqs "$BLAST_MAX_HITS"

    local exit_code=$?
    if [[ $exit_code -eq 0 ]]; then
        echo "  ✓ Done : $label"
    elif [[ $exit_code -eq 124 ]]; then
        echo "  ✗ TIMEOUT after 4h : $label — result file may be incomplete."
    else
        echo "  ✗ BLAST failed (exit $exit_code) : $label"
    fi
    echo ""
    return $exit_code
}

overall_status=0

# ── 1. Compare against LRG NM counterpart ─────────────────────────────────────
# NM files are named:  GENE_NM_XXXXXX.X.fasta
# Use glob to find the matching file for this gene
nm_files=("$LRG_NM_DIR/${gene_name}_NM_"*.fasta)

if [[ ${#nm_files[@]} -eq 1 && -f "${nm_files[0]}" ]]; then
    nm_subject="${nm_files[0]}"
    nm_outfile="$OUTPUT_DIR/${gene_name}_vs_NM.txt"
    run_blast "NM (LRG Source 3)" "$nm_subject" "$nm_outfile"
    [[ $? -ne 0 ]] && overall_status=1
elif [[ ${#nm_files[@]} -gt 1 ]]; then
    echo "  WARNING: Multiple NM files found for $gene_name — using the first match."
    echo "  Files found: ${nm_files[*]}"
    nm_subject="${nm_files[0]}"
    nm_outfile="$OUTPUT_DIR/${gene_name}_vs_NM.txt"
    run_blast "NM (LRG Source 3)" "$nm_subject" "$nm_outfile"
    [[ $? -ne 0 ]] && overall_status=1
else
    echo "  SKIP: No NM LRG file found for gene '$gene_name' in $LRG_NM_DIR"
    echo "        (Expected pattern: ${gene_name}_NM_*.fasta)"
    echo ""
fi

# ── 2. Compare against LRG NG counterpart ─────────────────────────────────────
# NG files are named:  GENE.fasta
ng_subject="$LRG_NG_DIR/${gene_name}.fasta"

if [[ -f "$ng_subject" ]]; then
    ng_outfile="$OUTPUT_DIR/${gene_name}_vs_NG.txt"
    run_blast "NG (LRG Source 3)" "$ng_subject" "$ng_outfile"
    [[ $? -ne 0 ]] && overall_status=1
else
    echo "  SKIP: No NG LRG file found for gene '$gene_name' in $LRG_NG_DIR"
    echo "        (Expected: ${gene_name}.fasta)"
    echo ""
fi

echo "=============================================="
echo "Task $SLURM_ARRAY_TASK_ID ($gene_name) finished : $(date)"
echo "Overall status : $([ $overall_status -eq 0 ] && echo 'OK' || echo 'ERRORS — check log')"
echo "=============================================="

exit $overall_status
