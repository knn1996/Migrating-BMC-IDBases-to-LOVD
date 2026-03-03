"""
process_embl_mutations.py
=========================
Author  : Generated for Assignment/Thesis – Mutation Extraction and Matching
Purpose : For every gene that has downloaded EMBL files and a mutation CSV:
            1. Parses each EMBL flat-file and converts it to a FASTA file.
            2. Extracts the CDS start position automatically from the EMBL
               feature table — no manual offset needed.
            3. Cross-checks every mutation (position + reference base) in
               <gene>_mut.csv against each FASTA sequence using the CDS offset.
            4. Declares a sequence as the "reference" if ALL mutations match.
            5. Writes per-gene CSV results and a master summary CSV.

Usage
-----
    python process_embl_mutations.py

Configuration
-------------
    Edit the CONSTANTS block below to change paths.

Dependencies
------------
    pip install biopython pandas openpyxl tqdm

Output files (written inside each gene sub-folder)
---------------------------------------------------
    <gene>/<accession>.fasta          – FASTA converted from EMBL
    <gene>/<gene>_mutation_check.csv  – per-mutation match table for all sequences
    <gene>/<gene>_reference_match.csv – which sequence(s) matched all mutations

Master outputs (written to BASE_DIR)
-------------------------------------
    process_log.txt       – timestamped log of every action
    process_summary.csv   – one row per gene: how many sequences checked,
                            how many matched all mutations
"""

# ──────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ──────────────────────────────────────────────────────────────────────────────
import os
import sys
import logging
import pandas as pd
from pathlib import Path
from datetime import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS  – edit these as needed
# ──────────────────────────────────────────────────────────────────────────────

# Root folder containing one sub-folder per gene (lower-case)
BASE_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Mutation Extraction and Matching"

# Merged Excel file – used only to get the canonical gene name list
INPUT_XLSX = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\IDBases_Summary_with_UniProt.xlsx"

# Name pattern for mutation CSVs inside each gene folder
# e.g.  ada\ada_mut.csv,  aicda\aicda_mut.csv, …
MUT_CSV_SUFFIX = "_mut.csv"

# ──────────────────────────────────────────────────────────────────────────────
# LOGGING SETUP
# ──────────────────────────────────────────────────────────────────────────────

log_path = os.path.join(BASE_DIR, "process_log.txt")

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

def extract_cds_offset(embl_path: Path) -> int:
    """
    Parse an EMBL flat-file and return the 0-based start position of the
    first CDS feature.  This is the offset that must be added to a
    CDS-relative mutation position to get the FASTA/sequence coordinate.

    How it works
    ------------
    EMBL flat-files list sequence features (gene, CDS, exon, …) in the FT
    block.  BioPython parses these into SeqFeature objects.  Each feature
    has a .location with a .start attribute that is already 0-based (Python
    convention).

    Example
    -------
    If the CDS annotation reads  "join(77..965)"  in the EMBL file,
    BioPython gives location.start == 76  (0-based).
    A mutation at CDS position 1 therefore sits at sequence index 76,
    i.e.  seq[76 + 1 - 1]  =  seq[76].

    Parameters
    ----------
    embl_path : Path
        Path to the .embl flat-file.

    Returns
    -------
    int
        0-based offset of the CDS start.  Returns 0 if no CDS feature is
        found (offset-free assumption).
    """
    try:
        record = SeqIO.read(str(embl_path), "embl")
    except Exception as exc:
        log.warning(f"    Could not parse EMBL file {embl_path.name}: {exc}")
        return 0

    for feature in record.features:
        if feature.type == "CDS":
            cds_start = int(feature.location.start)   # 0-based
            log.info(f"    CDS found at 0-based position {cds_start} in {embl_path.name}")
            return cds_start

    log.warning(f"    No CDS feature found in {embl_path.name} — using offset 0")
    return 0


def embl_to_fasta(embl_path: Path, fasta_path: Path) -> SeqRecord | None:
    """
    Convert an EMBL flat-file to a FASTA file and return the SeqRecord.

    Parameters
    ----------
    embl_path : Path
        Source .embl file.
    fasta_path : Path
        Destination .fasta file (created or overwritten).

    Returns
    -------
    SeqRecord | None
        The parsed record, or None if conversion failed.
    """
    try:
        record = SeqIO.read(str(embl_path), "embl")
    except Exception as exc:
        log.warning(f"    Cannot read {embl_path.name}: {exc}")
        return None

    # Write FASTA (overwrite if already exists from a previous run)
    SeqIO.write(record, str(fasta_path), "fasta")
    return record


def check_mutations(seq: str, df_mut: pd.DataFrame, cds_offset: int,
                    accession: str) -> list[dict]:
    """
    Check each mutation row against the sequence.

    Mutation positions in the CSV are CDS-relative (1-based).
    The actual 0-based index into the full sequence is:
        index = cds_offset + position - 1

    For range positions like "314-315", both endpoints are shifted the
    same way.

    Parameters
    ----------
    seq        : str   – full nucleotide sequence (uppercase)
    df_mut     : DataFrame with columns [mutation type, position,
                                         reference base/sequence]
    cds_offset : int   – 0-based index where the CDS starts in seq
    accession  : str   – label for results (e.g. 'X02994')

    Returns
    -------
    list[dict]
        One dict per mutation row.
    """
    results = []
    seq_len = len(seq)

    for _, row in df_mut.iterrows():
        pos_str  = str(row["position"]).strip()
        ref      = str(row["reference base/sequence"]).strip().upper()
        mut_type = str(row["mutation type"]).strip()

        try:
            if "-" in pos_str:
                # Range e.g. "314-315"  →  two 1-based CDS positions
                start_cds, end_cds = (int(x) for x in pos_str.split("-"))
                idx_start = cds_offset + start_cds - 1   # 0-based inclusive
                idx_end   = cds_offset + end_cds          # 0-based exclusive (slice)
                if idx_end > seq_len:
                    fasta_base = "OUT_OF_BOUNDS"
                else:
                    fasta_base = seq[idx_start:idx_end]
            else:
                pos       = int(pos_str)
                idx       = cds_offset + pos - 1          # 0-based
                if idx >= seq_len:
                    fasta_base = "OUT_OF_BOUNDS"
                else:
                    fasta_base = seq[idx]
        except (ValueError, IndexError):
            fasta_base = "PARSE_ERROR"

        results.append({
            "accession":         accession,
            "mutation_type":     mut_type,
            "position":          pos_str,
            "reference_base_df": ref,
            "fasta_base_at_pos": fasta_base,
            "cds_offset":        cds_offset,
            "is_match":          fasta_base == ref,
        })

    return results


def process_gene(gene: str, base_dir: str) -> dict:
    """
    Full pipeline for a single gene:
        1. Locate its folder and mutation CSV.
        2. For each .embl file:  convert → FASTA, extract CDS offset,
           check all mutations.
        3. Write per-gene CSV outputs.

    Parameters
    ----------
    gene     : str  – gene symbol as in Excel (e.g. 'ADA')
    base_dir : str  – root directory containing gene sub-folders

    Returns
    -------
    dict  – summary row for the master CSV
    """
    folder   = Path(base_dir) / gene.lower()
    mut_csv  = folder / f"{gene.lower()}{MUT_CSV_SUFFIX}"

    log.info(f"──── {gene} ────")

    # ── Sanity checks ─────────────────────────────────────────────────────────
    if not folder.exists():
        log.warning(f"  Folder not found: {folder}  – skipping")
        return {"Gene": gene, "EMBL files": 0, "Sequences checked": 0,
                "Full matches": 0, "No mutation CSV": True, "Error": "folder missing"}

    embl_files = sorted(folder.glob("*.embl"))
    if not embl_files:
        log.info(f"  No .embl files found in {folder}  – skipping")
        return {"Gene": gene, "EMBL files": 0, "Sequences checked": 0,
                "Full matches": 0, "No mutation CSV": not mut_csv.exists(), "Error": "no embl files"}

    if not mut_csv.exists():
        log.warning(f"  Mutation CSV not found: {mut_csv}  – FASTA conversion only")
        # Still convert EMBL → FASTA even without mutation data
        for embl_path in embl_files:
            fasta_path = embl_path.with_suffix(".fasta")
            embl_to_fasta(embl_path, fasta_path)
        return {"Gene": gene, "EMBL files": len(embl_files), "Sequences checked": 0,
                "Full matches": 0, "No mutation CSV": True, "Error": "no mutation csv"}

    # ── Load mutation table ────────────────────────────────────────────────────
    try:
        df_mut = pd.read_csv(mut_csv)
        df_mut.columns = df_mut.columns.str.strip()
        log.info(f"  Loaded {len(df_mut)} mutations from {mut_csv.name}")
    except Exception as exc:
        log.error(f"  Cannot read {mut_csv}: {exc}")
        return {"Gene": gene, "EMBL files": len(embl_files), "Sequences checked": 0,
                "Full matches": 0, "No mutation CSV": False, "Error": str(exc)}

    required = {"mutation type", "position", "reference base/sequence"}
    if not required.issubset(df_mut.columns):
        log.error(f"  Missing columns in {mut_csv.name}. Found: {df_mut.columns.tolist()}")
        return {"Gene": gene, "EMBL files": len(embl_files), "Sequences checked": 0,
                "Full matches": 0, "No mutation CSV": False,
                "Error": f"bad columns: {df_mut.columns.tolist()}"}

    # ── Process each EMBL file ────────────────────────────────────────────────
    all_results   = []
    n_checked     = 0

    for embl_path in embl_files:
        accession  = embl_path.stem                        # e.g. 'X02994'
        fasta_path = embl_path.with_suffix(".fasta")

        log.info(f"  Processing {embl_path.name} …")

        # Convert EMBL → FASTA
        record = embl_to_fasta(embl_path, fasta_path)
        if record is None:
            continue

        seq = str(record.seq).upper()

        # Auto-extract CDS offset from EMBL features
        cds_offset = extract_cds_offset(embl_path)

        # Check all mutations
        results = check_mutations(seq, df_mut, cds_offset, accession)
        all_results.extend(results)
        n_checked += 1

        n_match = sum(r["is_match"] for r in results)
        log.info(f"    {accession}: {n_match}/{len(results)} mutations matched "
                 f"(offset={cds_offset})")

    if not all_results:
        return {"Gene": gene, "EMBL files": len(embl_files), "Sequences checked": 0,
                "Full matches": 0, "No mutation CSV": False, "Error": "all embl parse failed"}

    # ── Build results DataFrame ────────────────────────────────────────────────
    df_results = pd.DataFrame(all_results)

    # Per-accession: did ALL mutations match?
    match_summary = (
        df_results.groupby("accession")["is_match"]
        .agg(total_mutations="count", matched="sum")
        .assign(all_match=lambda x: x["total_mutations"] == x["matched"])
        .reset_index()
    )
    match_summary.insert(0, "gene", gene)

    full_matches = match_summary[match_summary["all_match"]]
    n_full = len(full_matches)

    if n_full > 0:
        log.info(f"  ✓ {n_full} sequence(s) matched ALL mutations: "
                 f"{full_matches['accession'].tolist()}")
    else:
        best = match_summary.loc[match_summary["matched"].idxmax()]
        log.info(f"  ✗ No perfect match. Best: {best['accession']} "
                 f"({best['matched']}/{best['total_mutations']} mutations)")

    # ── Write per-gene output CSVs ─────────────────────────────────────────────
    detail_path  = folder / f"{gene.lower()}_mutation_check.csv"
    summary_path = folder / f"{gene.lower()}_reference_match.csv"

    df_results.to_csv(detail_path,  index=False)
    match_summary.to_csv(summary_path, index=False)
    log.info(f"  Saved: {detail_path.name}, {summary_path.name}")

    return {
        "Gene":              gene,
        "EMBL files":        len(embl_files),
        "Sequences checked": n_checked,
        "Full matches":      n_full,
        "Matched accessions": ", ".join(full_matches["accession"].tolist()) if n_full else "",
        "No mutation CSV":   False,
        "Error":             "",
    }


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    start_time = datetime.now()
    log.info("=" * 70)
    log.info("EMBL → FASTA Converter & Mutation Cross-Checker – starting")
    log.info(f"Base directory : {BASE_DIR}")
    log.info(f"Started at     : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("=" * 70)

    # ── Collect gene list from Excel ──────────────────────────────────────────
    if os.path.isfile(INPUT_XLSX):
        df = pd.read_excel(INPUT_XLSX)
        genes = df["name"].dropna().drop_duplicates().str.strip().tolist()
        log.info(f"Gene list loaded from {INPUT_XLSX}: {len(genes)} genes")
    else:
        # Fallback: discover genes from sub-folder names
        log.warning(f"{INPUT_XLSX} not found – discovering genes from sub-folders")
        genes = [
            d.name.upper()
            for d in Path(BASE_DIR).iterdir()
            if d.is_dir() and not d.name.startswith(".")
        ]
        log.info(f"Discovered {len(genes)} gene folders")

    log.info("")

    # ── Process each gene ─────────────────────────────────────────────────────
    summary_rows = []
    iterator = genes if not HAS_TQDM else tqdm(genes, desc="Genes", unit="gene")

    for gene in iterator:
        row = process_gene(gene, BASE_DIR)
        summary_rows.append(row)

    # ── Master summary ────────────────────────────────────────────────────────
    df_summary = pd.DataFrame(summary_rows)
    master_path = os.path.join(BASE_DIR, "process_summary.csv")
    df_summary.to_csv(master_path, index=False)

    log.info("")
    log.info("=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)
    log.info(f"  Genes processed          : {len(summary_rows)}")
    log.info(f"  Genes with full match    : {df_summary['Full matches'].gt(0).sum()}")
    log.info(f"  Genes without mut CSV    : {df_summary['No mutation CSV'].sum()}")
    log.info(f"  Genes with errors        : {df_summary['Error'].ne('').sum()}")
    log.info(f"  Master summary saved to  : {master_path}")
    log.info(f"  Log saved to             : {log_path}")
    elapsed = datetime.now() - start_time
    log.info(f"  Elapsed time             : {str(elapsed).split('.')[0]}")
    log.info("=" * 70)
    log.info("Done.")


if __name__ == "__main__":
    main()
