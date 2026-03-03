"""
download_embl_sequences.py
==========================
Author  : Generated for Assignment/Thesis – Mutation Extraction and Matching
Purpose : For every gene in IDBases_Summary_with_UniProt.xlsx, this script:
            1. Queries the UniProt REST API to retrieve all EMBL cross-references
               that correspond to mRNA nucleotide sequences.
            2. For each EMBL accession found, downloads the EMBL flat-file from
               the European Nucleotide Archive (ENA).
            3. Saves each flat-file into a gene-specific sub-folder under the
               user-supplied base directory.
            4. Writes a detailed log file (download_log.txt) and a summary
               report (download_summary.csv) alongside the downloaded files.

Usage
-----
    python download_embl_sequences.py

Configuration
-------------
    Edit the CONSTANTS block below to change the base directory or input file.

Dependencies
------------
    pip install requests pandas openpyxl tqdm

APIs used
---------
    UniProt REST API : https://rest.uniprot.org/uniprotkb/{UniProtID}.json
    ENA flat-file    : https://www.ebi.ac.uk/ena/browser/api/embl/{AccessionID}?download=true
"""

# ──────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ──────────────────────────────────────────────────────────────────────────────
import os
import sys
import time
import logging
import requests
import pandas as pd
from pathlib import Path
from datetime import datetime

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS  – edit these as needed
# ──────────────────────────────────────────────────────────────────────────────

# Path to the merged Excel file produced in the previous step
INPUT_XLSX = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\IDBases_Summary_with_UniProt.xlsx"

# Root folder that already contains one sub-folder per gene (lower-case names)
BASE_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Mutation Extraction and Matching"

# Seconds to wait between consecutive API requests (be polite to public APIs)
SLEEP_BETWEEN_ACCESSIONS = 0.5   # between each EMBL accession download
SLEEP_BETWEEN_GENES      = 0.3   # between UniProt lookups

# Maximum number of EMBL accessions to download per gene (None = unlimited)
# Set to e.g. 20 if you only want a representative subset per gene
MAX_ACCESSIONS_PER_GENE = None

# ──────────────────────────────────────────────────────────────────────────────
# LOGGING SETUP
# ──────────────────────────────────────────────────────────────────────────────

log_path = os.path.join(BASE_DIR, "download_log.txt")

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
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def get_mrna_embl_accessions(uniprot_id: str) -> list[str]:
    """
    Query the UniProt REST API and return a list of EMBL accession numbers
    whose moleculeType is 'mRNA' (case-insensitive match).

    The current UniProt REST API (2024+) JSON structure uses:
        "uniProtKBCrossReferences" (NOT the old "dbReferences" key)

    Each EMBL cross-reference entry looks like:
        {
          "database": "EMBL",                      ← filter on this  (NOT "type")
          "id": "M13792",                          ← nucleotide accession to download
          "properties": [
            {"key": "ProteinSequenceID", "value": "AAA51772.1"},
            {"key": "MoleculeType",      "value": "mRNA"}        ← filter on this
          ]
        }

    Note: "properties" is a LIST of {key, value} dicts in the new API,
    not a plain dict as in the legacy API.

    We match MoleculeType case-insensitively and log all distinct values seen
    so any future API changes are immediately visible in the log.

    Parameters
    ----------
    uniprot_id : str
        A UniProt accession such as 'P00813'.

    Returns
    -------
    list[str]
        Possibly-empty list of EMBL nucleotide accessions, e.g. ['M13792', …]
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
    except requests.RequestException as exc:
        log.error(f"    UniProt fetch failed for {uniprot_id}: {exc}")
        return []

    data = resp.json()

    # ── Collect all EMBL cross-references ─────────────────────────────────────
    # Field is "uniProtKBCrossReferences" in the current API (was "dbReferences")
    # "database" key identifies the DB type (was "type" in legacy API)
    all_xrefs = data.get("uniProtKBCrossReferences", [])
    embl_refs = [ref for ref in all_xrefs if ref.get("database") == "EMBL"]

    if not embl_refs:
        # Fallback: try legacy field name in case API version differs
        all_xrefs_legacy = data.get("dbReferences", [])
        embl_refs = [ref for ref in all_xrefs_legacy if ref.get("type") == "EMBL"]
        if embl_refs:
            log.info(f"  Using legacy 'dbReferences' field for {uniprot_id}")

    if not embl_refs:
        log.info(f"  No EMBL cross-references found for {uniprot_id}")
        return []

    # ── Helper: extract a property value from the list-of-dicts format ────────
    def get_prop(properties, key: str) -> str:
        """
        Properties in the new API are a list of {"key": ..., "value": ...} dicts.
        In the legacy API they were a plain dict.
        This helper handles both formats transparently.
        """
        if isinstance(properties, dict):
            # Legacy flat-dict format: {"moleculeType": "mRNA", ...}
            return properties.get(key, properties.get(key.lower(), ""))
        elif isinstance(properties, list):
            # Current list-of-dicts format
            for prop in properties:
                if prop.get("key", "").lower() == key.lower():
                    return prop.get("value", "")
        return ""

    # ── Log all distinct MoleculeType values seen (aids diagnosis) ────────────
    seen_types = sorted({
        get_prop(ref.get("properties", []), "MoleculeType") or "<missing>"
        for ref in embl_refs
    })
    log.info(f"  EMBL cross-refs found: {len(embl_refs)}  |  MoleculeType values seen: {seen_types}")

    # ── Filter for mRNA entries only ───────────────────────────────────────────
    accessions = []
    for ref in embl_refs:
        mol_type = get_prop(ref.get("properties", []), "MoleculeType")
        if mol_type.lower() == "mrna":
            accessions.append(ref["id"])

    return accessions


def download_embl_file(accession: str, dest_path: Path) -> bool:
    """
    Download the EMBL flat-file for a given accession from ENA and write it
    to *dest_path*.

    Parameters
    ----------
    accession : str
        ENA / EMBL nucleotide accession, e.g. 'AK123456'.
    dest_path : Path
        Full path (including filename) where the file should be saved.

    Returns
    -------
    bool
        True if the download succeeded and the file was written, False otherwise.
    """
    url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=true"
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException as exc:
        log.error(f"      ENA fetch failed for {accession}: {exc}")
        return False

    # ENA returns a plain-text error body (not HTTP 4xx) for unknown accessions
    if resp.text.strip().startswith("ERROR"):
        log.warning(f"      ENA returned error body for {accession}: {resp.text[:120]}")
        return False

    dest_path.write_text(resp.text, encoding="utf-8")
    return True


def gene_folder(base: str, gene_name: str) -> Path:
    """
    Resolve the on-disk folder for a gene.  The user stores folders in
    lower-case, so we convert the gene name accordingly.

    Parameters
    ----------
    base : str
        Root directory that contains all gene sub-folders.
    gene_name : str
        Gene symbol as it appears in the Excel file (e.g. 'ADA').

    Returns
    -------
    Path
        Path object pointing to  <base>/<gene_name_lower>.
        The folder is NOT created here – we only locate it.
    """
    return Path(base) / gene_name.lower()


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    start_time = datetime.now()
    log.info("=" * 70)
    log.info("EMBL mRNA Sequence Downloader – starting")
    log.info(f"Input file  : {INPUT_XLSX}")
    log.info(f"Output root : {BASE_DIR}")
    log.info(f"Started at  : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("=" * 70)

    # ── 1. Load the merged Excel file ─────────────────────────────────────────
    if not os.path.isfile(INPUT_XLSX):
        log.error(f"Input file not found: {INPUT_XLSX}")
        sys.exit(1)

    df = pd.read_excel(INPUT_XLSX)
    required_cols = {"name", "UniProt ID"}
    if not required_cols.issubset(df.columns):
        log.error(f"Expected columns {required_cols} in the Excel file. Found: {list(df.columns)}")
        sys.exit(1)

    # Drop duplicates and rows with missing UniProt IDs
    df_valid   = df.dropna(subset=["UniProt ID"]).drop_duplicates(subset=["name"])
    df_missing = df[df["UniProt ID"].isna()][["name"]].drop_duplicates()

    log.info(f"Total genes in file  : {len(df.drop_duplicates(subset=['name']))}")
    log.info(f"Genes with UniProt ID: {len(df_valid)}")
    log.info(f"Genes WITHOUT UniProt: {len(df_missing)}  → {df_missing['name'].tolist()}")
    log.info("")

    # ── 2. Per-gene processing ─────────────────────────────────────────────────
    summary_rows = []          # collected for the CSV report at the end
    iterator = df_valid.iterrows()
    if HAS_TQDM:
        iterator = tqdm(list(iterator), desc="Genes", unit="gene")

    for _, row in iterator:
        gene       = str(row["name"]).strip()
        uniprot_id = str(row["UniProt ID"]).strip()

        log.info(f"──── {gene}  ({uniprot_id}) ────")

        # Locate / validate the gene folder
        folder = gene_folder(BASE_DIR, gene)
        if not folder.exists():
            log.warning(f"  Folder not found, creating: {folder}")
            folder.mkdir(parents=True, exist_ok=True)

        # ── 2a. Get mRNA EMBL accessions from UniProt ─────────────────────────
        accessions = get_mrna_embl_accessions(uniprot_id)
        if not accessions:
            log.info(f"  No mRNA EMBL accessions found for {uniprot_id}")
            summary_rows.append({
                "Gene": gene, "UniProt ID": uniprot_id,
                "mRNA Accessions Found": 0, "Downloaded": 0,
                "Failed": 0, "Skipped (already existed)": 0,
            })
            time.sleep(SLEEP_BETWEEN_GENES)
            continue

        log.info(f"  Found {len(accessions)} mRNA accession(s)")

        # Optionally cap the number of accessions per gene
        if MAX_ACCESSIONS_PER_GENE:
            accessions = accessions[:MAX_ACCESSIONS_PER_GENE]
            log.info(f"  (Capped at {MAX_ACCESSIONS_PER_GENE})")

        # ── 2b. Download each EMBL flat-file ──────────────────────────────────
        n_downloaded = 0
        n_failed     = 0
        n_skipped    = 0

        for acc in accessions:
            dest = folder / f"{acc}.embl"

            # Skip files that were already downloaded in a previous run
            if dest.exists() and dest.stat().st_size > 0:
                log.info(f"    {acc}.embl  – already exists, skipping")
                n_skipped += 1
                continue

            log.info(f"    Downloading {acc} → {dest.name} …")
            success = download_embl_file(acc, dest)

            if success:
                size_kb = dest.stat().st_size / 1024
                log.info(f"    {acc}.embl  saved  ({size_kb:.1f} KB)")
                n_downloaded += 1
            else:
                n_failed += 1

            time.sleep(SLEEP_BETWEEN_ACCESSIONS)

        log.info(
            f"  {gene} done – downloaded: {n_downloaded}, "
            f"skipped: {n_skipped}, failed: {n_failed}"
        )
        summary_rows.append({
            "Gene": gene, "UniProt ID": uniprot_id,
            "mRNA Accessions Found": len(accessions),
            "Downloaded": n_downloaded,
            "Failed": n_failed,
            "Skipped (already existed)": n_skipped,
        })

        time.sleep(SLEEP_BETWEEN_GENES)

    # ── 3. Summary report ─────────────────────────────────────────────────────
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(BASE_DIR, "download_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    log.info("")
    log.info("=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)
    log.info(f"  Genes processed        : {len(summary_rows)}")
    log.info(f"  Total files downloaded : {summary_df['Downloaded'].sum()}")
    log.info(f"  Total files skipped    : {summary_df['Skipped (already existed)'].sum()}")
    log.info(f"  Total failures         : {summary_df['Failed'].sum()}")
    log.info(f"  Summary CSV saved to   : {summary_path}")
    log.info(f"  Log saved to           : {log_path}")

    elapsed = datetime.now() - start_time
    log.info(f"  Elapsed time           : {str(elapsed).split('.')[0]}")
    log.info("=" * 70)
    log.info("Done.")


if __name__ == "__main__":
    main()