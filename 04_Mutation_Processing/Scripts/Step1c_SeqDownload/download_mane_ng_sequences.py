import os
import time
import requests
from pathlib import Path

IDBASE_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\idbase"
OUT_DIR    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\Mane_Select_NG"
LOG_DIR    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Logs"
OUT_CSV    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step2_RefCheck\mane_ng_accessions.csv"

EFETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ELINK_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
ESUMMARY_URL= "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
SLEEP       = 0.4
MAX_RETRIES = 3

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

LOG_PATH = os.path.join(LOG_DIR, "download_mane_ng_sequences.log")


def get_gene_id(gene_symbol):
    params = {
        "db": "gene",
        "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
        "retmode": "json",
        "retmax": 5,
    }
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(ESEARCH_URL, params=params, timeout=30)
            r.raise_for_status()
            ids = r.json()["esearchresult"]["idlist"]
            if ids:
                return ids[0]
            return None
        except Exception as e:
            if attempt < MAX_RETRIES:
                time.sleep(SLEEP * attempt * 3)
    return None


def get_ng_accession(gene_id):
    params = {
        "dbfrom": "gene",
        "db": "nuccore",
        "id": gene_id,
        "linkname": "gene_nuccore_refseqgene",
        "retmode": "json",
    }
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(ELINK_URL, params=params, timeout=30)
            r.raise_for_status()
            data = r.json()
            linksetdbs = data["linksets"][0].get("linksetdbs", [])
            if not linksetdbs:
                return None
            nuccore_ids = linksetdbs[0]["links"]
            if not nuccore_ids:
                return None
            time.sleep(SLEEP)
            s = requests.get(ESUMMARY_URL, params={
                "db": "nuccore", "id": nuccore_ids[0], "retmode": "json"
            }, timeout=30)
            s.raise_for_status()
            acc = s.json()["result"][str(nuccore_ids[0])]["accessionversion"]
            return acc
        except Exception as e:
            if attempt < MAX_RETRIES:
                time.sleep(SLEEP * attempt * 3)
    return None


def fetch_fasta(ng_acc):
    params = {"db": "nuccore", "id": ng_acc, "rettype": "fasta", "retmode": "text"}
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(EFETCH_URL, params=params, timeout=120)
            if r.status_code == 200 and r.text.strip().startswith(">"):
                return r.text
            print(f"  [{attempt}] Unexpected response for {ng_acc}: HTTP {r.status_code}")
        except Exception as e:
            print(f"  [{attempt}] Request error for {ng_acc}: {e}")
        if attempt < MAX_RETRIES:
            time.sleep(SLEEP * attempt * 3)
    return None


def discover_genes(idbase_dir):
    import re
    genes = []
    for d in sorted(os.listdir(idbase_dir)):
        full = os.path.join(idbase_dir, d)
        if os.path.isdir(full) and d.lower().endswith("base") and d.lower() != "immunomebase":
            genes.append(re.sub(r"base$", "", d, flags=re.IGNORECASE))
    return genes


def main():
    genes = discover_genes(IDBASE_DIR)
    print(f"Genes to process: {len(genes)}")

    csv_rows = []
    downloaded = []
    already_exist = []
    failed = []
    no_ng = []

    with open(LOG_PATH, "w", encoding="utf-8") as log:
        for gene in genes:
            print(f"\n--- {gene} ---")

            gene_id = get_gene_id(gene)
            time.sleep(SLEEP)

            if not gene_id:
                msg = f"SKIP {gene}: no NCBI Gene ID found"
                print(f"  {msg}"); log.write(msg + "\n")
                no_ng.append(gene)
                csv_rows.append({"gene": gene, "gene_id": "", "ng_accession": "", "status": "no_gene_id"})
                continue

            ng_acc = get_ng_accession(gene_id)
            time.sleep(SLEEP)

            if not ng_acc or not ng_acc.startswith("NG_"):
                msg = f"SKIP {gene} (GeneID={gene_id}): no NG_ accession found"
                print(f"  {msg}"); log.write(msg + "\n")
                no_ng.append(gene)
                csv_rows.append({"gene": gene, "gene_id": gene_id, "ng_accession": ng_acc or "", "status": "no_ng"})
                continue

            print(f"  GeneID={gene_id}  NG={ng_acc}")

            out_path = Path(OUT_DIR) / f"{gene}_{ng_acc}.fasta"

            if out_path.exists() and out_path.stat().st_size > 0:
                msg = f"EXISTS {gene} ({ng_acc})"
                print(f"  {msg}"); log.write(msg + "\n")
                already_exist.append(gene)
                csv_rows.append({"gene": gene, "gene_id": gene_id, "ng_accession": ng_acc, "status": "already_exists"})
                continue

            print(f"  Downloading {ng_acc} ...")
            seq = fetch_fasta(ng_acc)
            time.sleep(SLEEP)

            if seq:
                out_path.write_text(seq, encoding="utf-8")
                size_kb = out_path.stat().st_size / 1024
                msg = f"OK {gene} ({ng_acc}): {size_kb:.1f} KB"
                print(f"  Saved {out_path.name} ({size_kb:.1f} KB)")
                log.write(msg + "\n")
                downloaded.append(gene)
                csv_rows.append({"gene": gene, "gene_id": gene_id, "ng_accession": ng_acc, "status": "downloaded"})
            else:
                msg = f"FAIL {gene} ({ng_acc}): download failed"
                print(f"  FAILED"); log.write(msg + "\n")
                failed.append(gene)
                csv_rows.append({"gene": gene, "gene_id": gene_id, "ng_accession": ng_acc, "status": "download_failed"})

    import csv
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "gene_id", "ng_accession", "status"])
        writer.writeheader()
        writer.writerows(csv_rows)

    print(f"\n=== Summary ===")
    print(f"Downloaded      : {len(downloaded)}")
    print(f"Already existed : {len(already_exist)}")
    print(f"No NG_ found    : {len(no_ng)}  — {no_ng}")
    print(f"Failed          : {len(failed)}  — {failed}")
    print(f"CSV saved to    : {OUT_CSV}")
    print(f"Log saved to    : {LOG_PATH}")


if __name__ == "__main__":
    main()
