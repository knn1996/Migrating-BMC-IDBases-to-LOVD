import os
import re
import csv
import time
import requests
from pathlib import Path

IDBASE_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\idbase"
OUT_DIR    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\Mane_Select_NG"
LOG_DIR    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Logs"
OUT_CSV    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step2_RefCheck\mane_ng_accessions.csv"
ALIAS_CSV  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\alias.csv"

ESEARCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ELINK_URL    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

SLEEP       = 0.4
MAX_RETRIES = 3
SKIP_GENES  = {"IGHG2", "IGHM", "SH2", "BTK"}

LOG_PATH = os.path.join(LOG_DIR, "download_mane_ng_sequences.log")

for d in [OUT_DIR, LOG_DIR, os.path.dirname(OUT_CSV)]:
    os.makedirs(d, exist_ok=True)


def load_aliases(path):
    if not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8") as f:
        return {r["gene"].strip().upper(): r["alias"].strip() for r in csv.DictReader(f)}


def discover_genes(idbase_dir):
    return [
        re.sub(r"base$", "", d, flags=re.IGNORECASE)
        for d in sorted(os.listdir(idbase_dir))
        if os.path.isdir(os.path.join(idbase_dir, d))
        and d.lower().endswith("base")
        and d.lower() != "immunomebase"
    ]


def ncbi_get(url, params):
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(url, params=params, timeout=30)
            r.raise_for_status()
            return r
        except Exception:
            if attempt < MAX_RETRIES:
                time.sleep(SLEEP * attempt * 3)
    return None


def get_gene_id(symbol):
    r = ncbi_get(ESEARCH_URL, {
        "db": "gene",
        "term": f"{symbol}[Gene/Protein Symbol] AND Homo sapiens[Organism]",
        "retmode": "json",
        "retmax": 5,
    })
    if not r:
        return None
    for gene_id in r.json()["esearchresult"]["idlist"]:
        time.sleep(SLEEP)
        s = ncbi_get(ESUMMARY_URL, {"db": "gene", "id": gene_id, "retmode": "json"})
        if s and s.json()["result"].get(str(gene_id), {}).get("name", "").upper() == symbol.upper():
            return gene_id
    return None


def get_ng_accession(gene_id):
    r = ncbi_get(ELINK_URL, {
        "dbfrom": "gene",
        "db": "nuccore",
        "id": gene_id,
        "linkname": "gene_nuccore_refseqgene",
        "retmode": "json",
    })
    if not r:
        return None
    linksetdbs = r.json()["linksets"][0].get("linksetdbs", [])
    if not linksetdbs:
        return None
    nuccore_ids = linksetdbs[0].get("links", [])
    if not nuccore_ids:
        return None
    time.sleep(SLEEP)
    s = ncbi_get(ESUMMARY_URL, {"db": "nuccore", "id": nuccore_ids[0], "retmode": "json"})
    if not s:
        return None
    result = s.json()["result"]
    entry = result.get(str(nuccore_ids[0])) or next((v for k, v in result.items() if k != "uids"), None)
    if not entry:
        return None
    acc = entry.get("accessionversion", "")
    return acc if acc.startswith("NG_") else None


def fetch_fasta(ng_acc):
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(EFETCH_URL, params={
                "db": "nuccore", "id": ng_acc, "rettype": "fasta", "retmode": "text"
            }, timeout=120)
            if r.status_code == 200 and r.text.strip().startswith(">"):
                return r.text
            print(f"  [{attempt}] Unexpected response for {ng_acc}: HTTP {r.status_code}")
        except Exception as e:
            print(f"  [{attempt}] Request error: {e}")
        if attempt < MAX_RETRIES:
            time.sleep(SLEEP * attempt * 3)
    return None


def resolve_gene(gene, aliases, log):
    for symbol in [gene, aliases.get(gene.upper())]:
        if not symbol:
            continue
        gene_id = get_gene_id(symbol)
        time.sleep(SLEEP)
        if not gene_id:
            continue
        log.write(f"  GeneID={gene_id} (searched as '{symbol}')\n")
        ng_acc = get_ng_accession(gene_id)
        time.sleep(SLEEP)
        if ng_acc:
            return gene_id, ng_acc, symbol
    return None, None, None


def csv_row(gene, searched_as="", gene_id="", ng_accession="", status=""):
    return {"gene": gene, "searched_as": searched_as, "gene_id": gene_id,
            "ng_accession": ng_accession, "status": status}


def main():
    aliases = load_aliases(ALIAS_CSV)
    genes   = discover_genes(IDBASE_DIR)
    print(f"Genes to process: {len(genes)}")
    print(f"Aliases loaded  : {aliases}")
    print(f"Skipping        : {SKIP_GENES}\n")

    rows          = []
    downloaded    = []
    already_exist = []
    failed        = []
    skipped       = []
    no_ng         = []

    with open(LOG_PATH, "w", encoding="utf-8") as log:
        for gene in genes:
            print(f"\n--- {gene} ---")

            if gene.upper() in SKIP_GENES:
                log.write(f"SKIP {gene}: in skip list\n")
                skipped.append(gene)
                rows.append(csv_row(gene, status="skipped"))
                continue

            gene_id, ng_acc, searched_as = resolve_gene(gene, aliases, log)

            if not gene_id:
                log.write(f"SKIP {gene}: no NCBI Gene ID found\n")
                no_ng.append(gene)
                rows.append(csv_row(gene, status="no_gene_id"))
                continue

            if not ng_acc:
                log.write(f"SKIP {gene} (GeneID={gene_id}): no NG_ accession\n")
                no_ng.append(gene)
                rows.append(csv_row(gene, searched_as, gene_id, status="no_ng"))
                continue

            print(f"  NG={ng_acc}")
            out_path = Path(OUT_DIR) / f"{gene}_{ng_acc}.fasta"

            if out_path.exists() and out_path.stat().st_size > 0:
                log.write(f"EXISTS {gene} ({ng_acc})\n")
                already_exist.append(gene)
                rows.append(csv_row(gene, searched_as, gene_id, ng_acc, "already_exists"))
                continue

            print(f"  Downloading {ng_acc} ...")
            seq = fetch_fasta(ng_acc)
            time.sleep(SLEEP)

            if seq:
                out_path.write_text(seq, encoding="utf-8")
                size_kb = out_path.stat().st_size / 1024
                print(f"  Saved {out_path.name} ({size_kb:.1f} KB)")
                log.write(f"OK {gene} ({ng_acc}): {size_kb:.1f} KB\n")
                downloaded.append(gene)
                rows.append(csv_row(gene, searched_as, gene_id, ng_acc, "downloaded"))
            else:
                print(f"  FAILED")
                log.write(f"FAIL {gene} ({ng_acc}): download failed\n")
                failed.append(gene)
                rows.append(csv_row(gene, searched_as, gene_id, ng_acc, "download_failed"))

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "searched_as", "gene_id", "ng_accession", "status"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n=== Summary ===")
    print(f"Downloaded        : {len(downloaded)}")
    print(f"Already existed   : {len(already_exist)}")
    print(f"Skipped           : {len(skipped)}  -- {skipped}")
    print(f"No NG_ found      : {len(no_ng)}  -- {no_ng}")
    print(f"Failed            : {len(failed)}  -- {failed}")
    print(f"CSV  : {OUT_CSV}")
    print(f"Log  : {LOG_PATH}")


if __name__ == "__main__":
    main()
