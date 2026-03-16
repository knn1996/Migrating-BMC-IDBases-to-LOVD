import os
import sys
import time
import logging
import argparse
import requests
import pandas as pd
from datetime import datetime

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

INPUT_XLSX = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\IDBases_Summary_with_UniProt.xlsx"
OUT_DIR    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\Reference sequences (Source 2)"

SLEEP_BETWEEN_ACCESSIONS = 0.5
SLEEP_BETWEEN_GENES      = 0.3
EXCLUDE_PREFIXES         = ("CH47",)

os.makedirs(OUT_DIR, exist_ok=True)
log_path = os.path.join(OUT_DIR, "download_log.txt")

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


def get_genomic_accessions(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
    except requests.RequestException as exc:
        log.error(f"  UniProt fetch failed for {uniprot_id}: {exc}")
        return []

    data   = resp.json()
    xrefs  = data.get("uniProtKBCrossReferences", []) or data.get("dbReferences", [])
    db_key = "database" if "uniProtKBCrossReferences" in data else "type"

    def get_prop(props, key):
        if isinstance(props, dict):
            return props.get(key, "")
        for p in (props or []):
            if p.get("key", "").lower() == key.lower():
                return p.get("value", "")
        return ""

    return [r["id"] for r in xrefs
            if r.get(db_key) == "EMBL"
            and get_prop(r.get("properties", []), "MoleculeType") == "Genomic_DNA"
            and not r["id"].startswith(EXCLUDE_PREFIXES)]


def get_ena_length(accession):
    url = f"https://www.ebi.ac.uk/ena/browser/api/summary/{accession}?format=json"
    try:
        resp = requests.get(url, timeout=20)
        resp.raise_for_status()
        data = resp.json()
        summaries = data.get("summaries", [])
        return summaries[0].get("sequenceLength", 0) if summaries else 0
    except Exception:
        return 0


def download_fasta(accession, dest_path):
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException as exc:
        log.error(f"    ENA fetch failed for {accession}: {exc}")
        return False

    if not resp.text.strip().startswith(">"):
        log.warning(f"    ENA did not return FASTA for {accession}: {resp.text[:120]}")
        return False

    with open(dest_path, "w", encoding="utf-8") as f:
        f.write(resp.text)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                        help="Query UniProt + ENA lengths only; no downloads")
    parser.add_argument("--min-length", type=int, default=0,
                        help="Only keep accessions with sequence length >= this value (bp)")
    args = parser.parse_args()

    start_time = datetime.now()
    mode_label = "DRY RUN" if args.dry_run else "DOWNLOAD"
    log.info("=" * 70)
    log.info(f"Genomic DNA Downloader (Source 2) – {mode_label}")
    if args.min_length:
        log.info(f"Min length filter  : {args.min_length:,} bp")
    log.info(f"Excluded prefixes  : {EXCLUDE_PREFIXES}")
    log.info(f"Output dir : {OUT_DIR}")
    log.info("=" * 70)

    df       = pd.read_excel(INPUT_XLSX)
    df_valid = df.dropna(subset=["UniProt ID"]).drop_duplicates(subset=["name"])
    log.info(f"Genes to process: {len(df_valid)}\n")

    summary_rows = []
    iterator = tqdm(list(df_valid.iterrows()), desc="Genes", unit="gene") if HAS_TQDM else df_valid.iterrows()

    for _, row in iterator:
        gene       = str(row["name"]).strip()
        uniprot_id = str(row["UniProt ID"]).strip()
        log.info(f"──── {gene}  ({uniprot_id}) ────")

        accessions = list(dict.fromkeys(get_genomic_accessions(uniprot_id)))

        if args.dry_run:
            acc_info = []
            for acc in accessions:
                length = get_ena_length(acc)
                passes = length >= args.min_length if args.min_length else True
                acc_info.append({"accession": acc, "length": length, "passes": passes})
                time.sleep(0.2)

            passing = [x for x in acc_info if x["passes"]]
            log.info(f"  {len(accessions)} total  |  {len(passing)} pass filter")
            for x in acc_info:
                flag = "✓" if x["passes"] else "✗"
                log.info(f"    {flag} {x['accession']}  {x['length']:,} bp")

            summary_rows.append({
                "Gene":               gene,
                "UniProt ID":         uniprot_id,
                "Total accessions":   len(accessions),
                "Pass filter":        len(passing),
                "Passing accessions": ", ".join(x["accession"] for x in passing),
            })
            time.sleep(SLEEP_BETWEEN_GENES)
            continue

        if args.min_length:
            filtered = []
            for acc in accessions:
                if get_ena_length(acc) >= args.min_length:
                    filtered.append(acc)
                time.sleep(0.2)
            accessions = filtered

        if not accessions:
            log.info(f"  No accessions pass filter")
            summary_rows.append({"Gene": gene, "Found": 0, "Downloaded": 0, "Failed": 0, "Skipped": 0})
            time.sleep(SLEEP_BETWEEN_GENES)
            continue

        n_downloaded = n_failed = n_skipped = 0

        for acc in accessions:
            dest = os.path.join(OUT_DIR, f"{gene}_{acc}.fasta")

            if os.path.exists(dest) and os.path.getsize(dest) > 0:
                log.info(f"    {gene}_{acc}.fasta – exists, skipping")
                n_skipped += 1
                continue

            log.info(f"    Downloading {acc} …")
            if download_fasta(acc, dest):
                log.info(f"    {gene}_{acc}.fasta saved ({os.path.getsize(dest)/1024:.1f} KB)")
                n_downloaded += 1
            else:
                n_failed += 1

            time.sleep(SLEEP_BETWEEN_ACCESSIONS)

        log.info(f"  {gene}: downloaded={n_downloaded} skipped={n_skipped} failed={n_failed}")
        summary_rows.append({"Gene": gene, "Found": len(accessions),
                              "Downloaded": n_downloaded, "Failed": n_failed, "Skipped": n_skipped})
        time.sleep(SLEEP_BETWEEN_GENES)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(OUT_DIR, "download_summary.csv"), index=False)

    log.info("\n" + "=" * 70)
    if args.dry_run:
        log.info(f"  Genes surveyed : {len(summary_rows)}")
        log.info(f"  Summary saved  : {os.path.join(OUT_DIR, 'download_summary.csv')}")
    else:
        log.info(f"  Total downloaded : {summary_df['Downloaded'].sum()}")
        log.info(f"  Total skipped    : {summary_df['Skipped'].sum()}")
        log.info(f"  Total failed     : {summary_df['Failed'].sum()}")
    log.info(f"  Elapsed : {str(datetime.now() - start_time).split('.')[0]}")
    log.info("=" * 70)


if __name__ == "__main__":
    main()
