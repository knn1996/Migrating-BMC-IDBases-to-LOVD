"""
copy_reference_sequences.py
============================
Author  : Generated for Assignment/Thesis – Mutation Extraction and Matching
Purpose : 1. Creates a "Reference sequence" folder under the thesis directory.
          2. Reads process_summary.csv to find genes with matched accessions.
          3. For genes with 1–3 matched accessions, copies each matching FASTA
             file into the Reference sequence folder, renaming it:
                 <accession>_<GENE>.fasta
             e.g.  X02994_ADA.fasta
          4. Genes with 0 matches or more than 3 matches are skipped and logged.
          5. Writes a copy_log.txt and copy_report.csv in the Reference sequence
             folder summarising every action taken.

Usage
-----
    python copy_reference_sequences.py

Configuration
-------------
    Edit the CONSTANTS block below if your paths differ.

Dependencies
------------
    pip install pandas openpyxl
"""

# ──────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ──────────────────────────────────────────────────────────────────────────────
import os
import sys
import shutil
import logging
import pandas as pd
from pathlib import Path
from datetime import datetime

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS  – edit these as needed
# ──────────────────────────────────────────────────────────────────────────────

_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
THESIS_DIR   = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
MATCHING_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', 'Mutation Extraction and Matching'))
SUMMARY_CSV  = os.path.join(MATCHING_DIR, "process_summary.csv")
OUTPUT_DIR   = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Reference sequence")

# Genes with more than this many matches are skipped (too ambiguous)
MAX_MATCHES  = 3

# ──────────────────────────────────────────────────────────────────────────────
# LOGGING SETUP
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)
log_path = os.path.join(OUTPUT_DIR, "copy_log.txt")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(log_path, encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    start_time = datetime.now()
    log.info("=" * 70)
    log.info("Reference Sequence Copier – starting")
    log.info(f"Summary CSV  : {SUMMARY_CSV}")
    log.info(f"Source root  : {MATCHING_DIR}")
    log.info(f"Output folder: {OUTPUT_DIR}")
    log.info(f"Started at   : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("=" * 70)

    # ── 1. Load summary CSV ───────────────────────────────────────────────────
    if not os.path.isfile(SUMMARY_CSV):
        log.error(f"Summary CSV not found: {SUMMARY_CSV}")
        sys.exit(1)

    df = pd.read_csv(SUMMARY_CSV)
    df.columns = df.columns.str.strip()
    log.info(f"Loaded {len(df)} rows from process_summary.csv")
    log.info("")

    report_rows = []

    # ── 2. Process each gene row ──────────────────────────────────────────────
    for _, row in df.iterrows():
        gene         = str(row.get("Gene", "")).strip()
        full_matches = int(row.get("Full matches", 0))
        matched_str  = str(row.get("Matched accessions", "")).strip()

        # ── Skip: no matches ──────────────────────────────────────────────────
        if full_matches == 0 or matched_str in ("", "nan"):
            log.info(f"  {gene:<15} – SKIP  (0 matches)")
            report_rows.append({
                "Gene": gene, "Action": "skipped – no matches",
                "Accession": "", "Destination file": "",
            })
            continue

        # Parse comma-separated accessions
        accessions = [a.strip() for a in matched_str.split(",") if a.strip()]

        # ── Skip: too many matches ────────────────────────────────────────────
        if len(accessions) > MAX_MATCHES:
            log.info(f"  {gene:<15} – SKIP  ({len(accessions)} matches > {MAX_MATCHES} limit): "
                     f"{accessions}")
            report_rows.append({
                "Gene": gene,
                "Action": f"skipped – {len(accessions)} matches exceeds limit of {MAX_MATCHES}",
                "Accession": matched_str, "Destination file": "",
            })
            continue

        log.info(f"  {gene:<15} – {len(accessions)} match(es): {accessions}")

        # ── Copy each matched FASTA ───────────────────────────────────────────
        gene_folder = Path(MATCHING_DIR) / gene.lower()

        for accession in accessions:
            src = gene_folder / f"{accession}.fasta"
            dst_name = f"{accession}_{gene}.fasta"
            dst = Path(OUTPUT_DIR) / dst_name

            if not src.exists():
                log.warning(f"    {accession}: source FASTA not found: {src}")
                report_rows.append({
                    "Gene": gene, "Action": "ERROR – source fasta missing",
                    "Accession": accession, "Destination file": str(dst),
                })
                continue

            shutil.copy2(src, dst)
            log.info(f"    Copied  {src.name}  →  {dst_name}")
            report_rows.append({
                "Gene": gene, "Action": "copied",
                "Accession": accession, "Destination file": dst_name,
            })

    # ── 3. Write copy report CSV ──────────────────────────────────────────────
    df_report = pd.DataFrame(report_rows)
    report_path = os.path.join(OUTPUT_DIR, "copy_report.csv")
    df_report.to_csv(report_path, index=False)

    n_copied  = (df_report["Action"] == "copied").sum()
    n_skipped = df_report["Action"].str.startswith("skipped").sum()
    n_errors  = df_report["Action"].str.startswith("ERROR").sum()

    log.info("")
    log.info("=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)
    log.info(f"  Files copied             : {n_copied}")
    log.info(f"  Genes skipped            : {n_skipped}")
    log.info(f"  Errors (missing FASTA)   : {n_errors}")
    log.info(f"  Report saved to          : {report_path}")
    log.info(f"  Log saved to             : {log_path}")
    elapsed = datetime.now() - start_time
    log.info(f"  Elapsed time             : {str(elapsed).split('.')[0]}")
    log.info("=" * 70)
    log.info("Done.")


if __name__ == "__main__":
    main()
