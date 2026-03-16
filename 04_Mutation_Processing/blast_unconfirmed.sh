#!/bin/bash
#SBATCH -A lu2025-2-88
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH -J BLAST_ALL
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load GCC/13.3.0
module load OpenMPI/5.0.3
module load BLAST+/2.16.0
module load blast_databases/v5

FASTA_DIR="/home/khoi1996/Documents/DNA_sequences/confirmed_reference_sequence"
OUTPUT_DIR="/home/khoi1996/Documents/DNA_sequences/blast_results"
BLAST_DB="${BLASTDB}/nt"
THREADS=8
OUTFMT="7 qseqid sacc stitle pident length qcovs evalue bitscore"
SUMMARY="${OUTPUT_DIR}/blast_summary.tsv"

mkdir -p logs "${OUTPUT_DIR}"

echo "BLASTDB=$BLASTDB"
ls $BLASTDB/nt* 2>/dev/null | head -5

if ! command -v blastn &>/dev/null; then
    echo "ERROR: blastn not found in PATH after module load." >&2
    exit 1
fi

mapfile -t QUERIES < <(find "${FASTA_DIR}" -maxdepth 1 -name "*.fasta" | sort)
echo "Found ${#QUERIES[@]} files to BLAST"

printf "gene\ttop_accession\ttop_description\tpident\tcoverage\tevalue\tbitscore\n" > "${SUMMARY}"

for FASTA in "${QUERIES[@]}"; do
    GENE=$(basename "${FASTA}" .fasta)
    OUT_TSV="${OUTPUT_DIR}/${GENE}_blast.tsv"
    echo "  BLAST: $(basename "${FASTA}") -> $(basename "${OUT_TSV}")"

    blastn \
        -query           "${FASTA}" \
        -db              "${BLAST_DB}" \
        -out             "${OUT_TSV}" \
        -outfmt          "${OUTFMT}" \
        -num_threads     "${THREADS}" \
        -max_target_seqs 5 \
        -evalue          1e-5

    if [[ $? -ne 0 ]]; then
        echo "  ERROR on ${GENE}" >&2
        printf "%s\tERROR\t-\t-\t-\t-\t-\n" "${GENE}" >> "${SUMMARY}"
        continue
    fi

    TOP=$(grep -v '^#' "${OUT_TSV}" | head -1)
    if [[ -z "${TOP}" ]]; then
        printf "%s\tNO_HIT\t-\t-\t-\t-\t-\n" "${GENE}" >> "${SUMMARY}"
    else
        ACC=$(echo "${TOP}"    | cut -f2)
        DESC=$(echo "${TOP}"   | cut -f3)
        PIDENT=$(echo "${TOP}" | cut -f4)
        COV=$(echo "${TOP}"    | cut -f6)
        EVAL=$(echo "${TOP}"   | cut -f7)
        BITS=$(echo "${TOP}"   | cut -f8)
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${GENE}" "${ACC}" "${DESC}" "${PIDENT}" "${COV}" "${EVAL}" "${BITS}" >> "${SUMMARY}"
    fi
done

echo "Done. Summary -> ${SUMMARY}"