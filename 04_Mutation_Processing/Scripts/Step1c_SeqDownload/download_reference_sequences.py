"""
download_reference_sequences.py
================================
Reads origin_paper_with_refseq.tsv, collects all unique
(idbase, reference_sequence) pairs, and downloads a FASTA file
for each reference sequence accession from ENA (EMBL flat-file → FASTA)
or NCBI (efetch FASTA), saving as:

    <reference_sequence>_<IDBASE>.fasta
    e.g. M13792_ADA.fasta

Output folder: 05_Testing\reference_check

Logic
-----
- Tries ENA first (EMBL flat-file, converted to FASTA inline).
- Falls back to NCBI efetch if ENA returns an error body.
- Skips accessions that are "Unknown" or contain spaces (inferred/notes).
- Skips files already downloaded.

Dependencies: pip install requests pandas biopython
"""

import os
import re
import csv
import time
import requests
from pathlib import Path

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_TSV  = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..', '05_Testing', 'reference_check', 'origin_paper_with_refseq.tsv'))
OUT_DIR    = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..', '05_Testing', 'reference_check'))

ENA_FASTA_URL  = "https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}?download=true"
NCBI_FASTA_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

SLEEP = 0.5

RE_ACC = re.compile(r"^[A-Z]{1,2}\d{5,6}(\.\d+)?$")


def clean_accession(raw):
    """
    Extract the bare accession from strings like:
      "M13792"
      "M13792 (inferred)"
      "M13792, X02994"
    Returns a list of clean accession strings.
    """
    parts = re.split(r"[,;]", raw)
    accessions = []
    for p in parts:
        p = p.strip()
        m = re.match(r"([A-Za-z]{1,2}\d{5,8}(?:\.\d+)?)", p)
        if m:
            accessions.append(m.group(1).upper())
    return accessions


def download_fasta_ena(acc):
    url = ENA_FASTA_URL.format(acc=acc)
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        if r.text.strip().startswith("ERROR") or r.text.strip().startswith("<?xml"):
            return None
        if r.text.strip().startswith(">"):
            return r.text
    except requests.RequestException as e:
        print(f"      ENA error for {acc}: {e}")
    return None


def download_fasta_ncbi(acc):
    params = {
        "db": "nuccore", "id": acc,
        "rettype": "fasta", "retmode": "text",
    }
    try:
        r = requests.get(NCBI_FASTA_URL, params=params, timeout=60)
        r.raise_for_status()
        if r.text.strip().startswith(">"):
            return r.text
    except requests.RequestException as e:
        print(f"      NCBI error for {acc}: {e}")
    return None


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    rows = []
    with open(INPUT_TSV, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append(r)

    seen = set()
    pairs = []
    for r in rows:
        idbase = r["idbase"]
        raw    = r.get("reference_sequence", "")
        if not raw or raw.lower() in ("unknown", "none", ""):
            continue
        for acc in clean_accession(raw):
            key = (idbase, acc)
            if key not in seen:
                seen.add(key)
                pairs.append((idbase, acc))

    print(f"Unique (idbase, accession) pairs to download: {len(pairs)}")

    n_ok = n_skip = n_fail = 0

    for idbase, acc in pairs:
        fname = f"{acc}_{idbase}.fasta"
        dest  = Path(OUT_DIR) / fname

        if dest.exists() and dest.stat().st_size > 0:
            print(f"  SKIP  {fname}  (already exists)")
            n_skip += 1
            continue

        print(f"  Downloading {acc} for {idbase} → {fname}")

        fasta = download_fasta_ena(acc)
        if not fasta:
            time.sleep(SLEEP)
            fasta = download_fasta_ncbi(acc)

        if fasta:
            dest.write_text(fasta, encoding="utf-8")
            kb = dest.stat().st_size / 1024
            print(f"    saved  ({kb:.1f} KB)")
            n_ok += 1
        else:
            print(f"    FAILED to download {acc}")
            n_fail += 1

        time.sleep(SLEEP)

    print(f"\nDownloaded : {n_ok}")
    print(f"Skipped    : {n_skip}")
    print(f"Failed     : {n_fail}")
    print(f"Output dir : {OUT_DIR}")


if __name__ == "__main__":
    main()
