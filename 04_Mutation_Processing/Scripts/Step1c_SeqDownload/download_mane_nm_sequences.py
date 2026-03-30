import os
import re
import time
import csv
import requests

IDBASE_SUMMARY = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\IDBases_Summary.csv"
ALIAS_CSV      = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\alias.csv"
MANE_FASTA     = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\MANE.GRCh38.v1.5.refseq_rna.fna"
OUT_DIR        = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\Mane_Select_NM"
LOG_DIR        = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Logs"

NCBI_ESEARCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_EFETCH   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

SLEEP = 0.4

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)


def load_genes(path):
    with open(path, encoding="utf-8") as f:
        return [row["gene_name"].strip() for row in csv.DictReader(f) if row["gene_name"].strip()]


def load_aliases(path):
    aliases = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            aliases[row["gene"].strip().upper()] = row["alias"].strip()
    return aliases


def parse_mane_fasta(path):
    print("Parsing local MANE FASTA ...")
    gene_to_nm  = {}
    nm_to_entry = {}
    re_header   = re.compile(r"^>(NM_\S+)\s+.*\(([^)]+)\)")

    current_nm = current_header = None
    current_seq = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_nm:
                    nm_to_entry[current_nm] = (current_header, "".join(current_seq))
                current_seq = []
                m = re_header.match(line)
                if m:
                    current_nm     = m.group(1)
                    current_header = line
                    gene_symbol    = m.group(2).strip().upper()
                    gene_to_nm.setdefault(gene_symbol, current_nm)
                else:
                    current_nm = current_header = None
            elif current_nm:
                current_seq.append(line)

    if current_nm:
        nm_to_entry[current_nm] = (current_header, "".join(current_seq))

    print(f"  {len(gene_to_nm)} genes, {len(nm_to_entry)} NM_ entries.")
    return gene_to_nm, nm_to_entry


def write_fasta(dest_path, header, seq, width=60):
    with open(dest_path, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i+width] + "\n")


def get_latest_nm_api(gene_symbol):
    r = requests.get(NCBI_ESEARCH, params={
        "db": "gene", "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]",
        "retmode": "json", "retmax": 1,
    }, timeout=30)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    time.sleep(SLEEP)

    r2 = requests.get(NCBI_ESUMMARY, params={"db": "gene", "id": ids[0], "retmode": "json"}, timeout=30)
    r2.raise_for_status()
    result = r2.json().get("result", {}).get(str(ids[0]), {})
    transcript_list = result.get("transcriptlist", []) or result.get("transcripts", [])

    latest_nm, latest_ver = None, -1
    for t in transcript_list:
        acc = t.get("accessionversion", "") or t.get("accver", "")
        if not acc.startswith("NM_"):
            continue
        try:
            ver = int(acc.split(".")[1])
            if ver > latest_ver:
                latest_ver, latest_nm = ver, acc
        except (IndexError, ValueError):
            pass
    if latest_nm:
        return latest_nm

    time.sleep(SLEEP)
    r3 = requests.get(NCBI_ESEARCH, params={
        "db": "nuccore", "retmode": "json", "retmax": 10, "sort": "relevance",
        "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism] AND mRNA[Filter]",
    }, timeout=30)
    r3.raise_for_status()
    ids3 = r3.json().get("esearchresult", {}).get("idlist", [])
    if not ids3:
        return None

    time.sleep(SLEEP)
    r4 = requests.get(NCBI_ESUMMARY, params={"db": "nuccore", "id": ",".join(ids3), "retmode": "json"}, timeout=30)
    r4.raise_for_status()
    result4 = r4.json().get("result", {})

    candidates = []
    for uid in ids3:
        acc = result4.get(uid, {}).get("accessionversion", "")
        if acc.startswith("NM_"):
            try:
                ver = int(acc.split(".")[1])
            except (IndexError, ValueError):
                ver = 0
            candidates.append((ver, acc))

    return candidates[0][1] if sorted(candidates, reverse=True) else None


def download_fasta_api(nm, dest_path):
    r = requests.get(NCBI_EFETCH, params={
        "db": "nuccore", "id": nm, "rettype": "fasta", "retmode": "text",
    }, timeout=60)
    r.raise_for_status()
    if not r.text.strip().startswith(">"):
        return False
    with open(dest_path, "w", encoding="utf-8") as f:
        f.write(r.text)
    return True


def resolve_and_save(gene, aliases, gene_to_nm, nm_to_entry):
    lookups = [(gene, ""), (aliases.get(gene.upper(), ""), "_alias")]

    for lookup, suffix in lookups:
        if not lookup:
            continue
        nm = gene_to_nm.get(lookup.upper())
        if nm and nm in nm_to_entry:
            header, seq = nm_to_entry[nm]
            dest = os.path.join(OUT_DIR, f"{gene}_{nm}.fasta")
            write_fasta(dest, header, seq)
            return nm, f"mane_local{suffix}", True

    for lookup, suffix in lookups:
        if not lookup:
            continue
        nm = get_latest_nm_api(lookup)
        time.sleep(SLEEP)
        if nm:
            dest = os.path.join(OUT_DIR, f"{gene}_{nm}.fasta")
            ok = download_fasta_api(nm, dest)
            time.sleep(SLEEP)
            return nm, f"latest_nm_api{suffix}", ok

    return None, None, False


def main():
    genes              = load_genes(IDBASE_SUMMARY)
    aliases            = load_aliases(ALIAS_CSV)
    gene_to_nm, nm_to_entry = parse_mane_fasta(MANE_FASTA)

    log_path     = os.path.join(LOG_DIR, "download_mane_nm.log")
    summary_rows = []

    print(f"Total genes: {len(genes)}")

    with open(log_path, "w", encoding="utf-8") as log:
        for gene in genes:
            existing = [f for f in os.listdir(OUT_DIR) if f.upper().startswith(gene.upper() + "_NM_")]
            if existing:
                print(f"  SKIP {gene}: {existing[0]}")
                log.write(f"SKIP\t{gene}\talready_exists\t{existing[0]}\n")
                summary_rows.append({"gene": gene, "status": "skipped", "nm": existing[0], "source": ""})
                continue

            nm, source, ok = resolve_and_save(gene, aliases, gene_to_nm, nm_to_entry)

            if nm is None:
                print(f"  FAIL {gene}: no NM_ found")
                log.write(f"FAIL\t{gene}\tno_nm_found\t\n")
                summary_rows.append({"gene": gene, "status": "no_nm_found", "nm": "", "source": ""})
            elif not ok:
                print(f"  FAIL {gene}: download failed for {nm}")
                log.write(f"FAIL\t{gene}\tdownload_failed\t{nm}\n")
                summary_rows.append({"gene": gene, "status": "download_failed", "nm": nm, "source": source})
            else:
                print(f"  SAVED {gene}_{nm}.fasta [{source}]")
                log.write(f"OK\t{gene}\t{nm}\t{source}\n")
                summary_rows.append({"gene": gene, "status": "ok", "nm": nm, "source": source})

    summary_path = os.path.join(LOG_DIR, "download_mane_nm_summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "status", "nm", "source"])
        writer.writeheader()
        writer.writerows(summary_rows)

    ok_count   = sum(1 for r in summary_rows if r["status"] == "ok")
    fail_count = sum(1 for r in summary_rows if r["status"] not in ("ok", "skipped"))
    skip_count = sum(1 for r in summary_rows if r["status"] == "skipped")
    print(f"\nDone. OK={ok_count}  Skipped={skip_count}  Failed={fail_count}")
    print(f"Log    : {log_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()

