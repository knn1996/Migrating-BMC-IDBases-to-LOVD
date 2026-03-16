#!/bin/bash
#SBATCH --job-name=generate_bed
#SBATCH --output=/home/khoi1996/Documents/BED/LRG/logs/generate_bed_%j.out
#SBATCH --error=/home/khoi1996/Documents/BED/LRG/logs/generate_bed_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

# ── Setup ──────────────────────────────────────────────────────────────────────

mkdir -p /home/khoi1996/Documents/BED/LRG/logs
mkdir -p /home/khoi1996/Documents/BED/LRG/hg38

echo "Job ID     : $SLURM_JOB_ID"
echo "Node       : $SLURMD_NODENAME"
echo "Started at : $(date)"
echo "──────────────────────────────────────────"

# ── Run ────────────────────────────────────────────────────────────────────────

python3 /home/khoi1996/Documents/BED/LRG/generate_bed.py

# ── Done ───────────────────────────────────────────────────────────────────────

echo "──────────────────────────────────────────"
echo "Finished at : $(date)"
