import os
import sys
import time
import logging
import xml.etree.ElementTree as ET
import requests
import pandas as pd
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
_THESIS_DIR = (_SCRIPT_DIR / ".." / ".." / "..").resolve()
_PROC_DIR   = _THESIS_DIR / "04_Mutation_Processing"

INPUT_CSV   = _THESIS_DIR / "04_Mutation_Processing" / "Output" / "IDBases_Summary.csv"
OUT_DIR     = _PROC_DIR / "DNA sequences" / "Reference sequences (Source 2)"
CSV_OUT     = _SCRIPT_DIR / "ncbi_genomic_sequences.csv"

NCBI_BASE   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_FASTA   = "https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}?download=true"
ORGANISM    = "Homo sapiens"
MIN_BP      = 5000
SLEEP_API   = 0.4
SLEEP_GENE  = 0.3

# RefSeq prefixes that ENA does not host — skip these outright
REFSEQ_PREFIXES = ("NC_", "NG_", "NM_", "NR_", "NT_", "NW_", "NZ_",
                   "XM_", "XR_", "AC_", "NP_", "XP_", "WP_")

os.makedirs(OUT_DIR, exist_ok=True)
log_path = os.path.join(OUT_DIR, "ncbi_download_log.txt")

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


def ncbi_esearch(gene_symbol):
    params = {
        "db": "gene",
        "term": f'"{gene_symbol}"[Gene Name] AND "{ORGANISM}"[Organism]',
        "retmode": "json",
        "retmax": "1",
    }
    try:
        r = requests.get(f"{NCBI_BASE}/esearch.fcgi", params=params, timeout=30)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        return ids[0] if ids else None
    except Exception as exc:
        log.error(f"  esearch failed for {gene_symbol}: {exc}")
        return None


def ncbi_gene_xml(gene_id):
    params = {
        "db": "gene",
        "id": gene_id,
        "rettype": "xml",
        "retmode": "xml",
    }
    try:
        r = requests.get(f"{NCBI_BASE}/efetch.fcgi", params=params, timeout=60)
        r.raise_for_status()
        return r.text
    except Exception as exc:
        log.error(f"  efetch failed for gene_id={gene_id}: {exc}")
        return None


def parse_related_sequences(xml_text):
    """
    Extract accession.version for entries under:
      Gene-commentary[heading='Related Sequences'] > Gene-commentary_products >
        Gene-commentary[type=genomic, heading='Genomic']

    Skips RefSeq accessions (NC_, NG_, etc.) that ENA does not host.
    Returns list of accession strings like ['AL139352.16', 'HH769637.1']
    """
    accessions = []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        log.error(f"  XML parse error: {exc}")
        return accessions

    for gc in root.iter("Gene-commentary"):
        if gc.findtext("Gene-commentary_heading") != "Related Sequences":
            continue
        products = gc.find("Gene-commentary_products")
        if products is None:
            continue
        for child in products.findall("Gene-commentary"):
            ctype = child.find("Gene-commentary_type")
            if ctype is None or ctype.attrib.get("value") != "genomic":
                continue
            if child.findtext("Gene-commentary_heading") != "Genomic":
                continue
            acc = child.findtext("Gene-commentary_accession")
            ver = child.findtext("Gene-commentary_version")
            if not acc:
                continue
            full = f"{acc}.{ver}" if ver else acc
            if any(full.startswith(p) for p in REFSEQ_PREFIXES):
                continue
            if full not in accessions:
                accessions.append(full)

    return accessions


def existing_accessions_for_gene(gene, out_dir):
    existing = set()
    prefix = f"{gene}_"
    for fname in os.listdir(out_dir):
        if fname.startswith(prefix) and fname.endswith(".fasta"):
            acc_part = fname[len(prefix):-len(".fasta")]
            existing.add(acc_part)
    return existing


def accession_exists(acc_full, existing_set):
    base = acc_full.split(".")[0]
    return acc_full in existing_set or base in existing_set


def fasta_length(fasta_text):
    seq = "".join(
        line for line in fasta_text.splitlines()
        if not line.startswith(">")
    ).replace(" ", "")
    return len(seq)


def download_fasta(acc_full, dest_path):
    """
    Download from ENA. Returns (success, bp_count).
    Returns (False, 0) if download fails or sequence < MIN_BP.
    """
    bare = acc_full.split(".")[0]
    for acc in [acc_full, bare]:
        url = ENA_FASTA.format(acc=acc)
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            if r.text.strip().startswith(">"):
                bp = fasta_length(r.text)
                if bp < MIN_BP:
                    log.info(f"    SKIPPED {acc_full}: only {bp} bp (< {MIN_BP} bp threshold)")
                    return False, bp
                with open(dest_path, "w", encoding="utf-8") as f:
                    f.write(r.text)
                return True, bp
        except Exception:
            pass
        time.sleep(0.3)
    return False, 0


def main():
    log.info("=" * 70)
    log.info("NCBI Related Sequences – fetch, compare & download")
    log.info(f"Output dir  : {OUT_DIR}")
    log.info(f"Min bp      : {MIN_BP}")
    log.info("=" * 70)

    df = pd.read_csv(INPUT_CSV)
    df_valid = df.drop_duplicates(subset=["gene_name"])
    genes = df_valid["gene_name"].str.strip().tolist()
    log.info(f"Genes loaded: {len(genes)}")

    rows = []

    for gene in genes:
        log.info(f"\n──── {gene} ────")

        gene_id = ncbi_esearch(gene)
        time.sleep(SLEEP_API)

        if not gene_id:
            log.warning(f"  No NCBI Gene ID found for {gene}")
            rows.append({
                "Gene": gene,
                "NCBI_Link": "",
                "All_Related_Sequences": "",
                "New_Sequences": "",
                "Status": "Gene ID not found",
            })
            continue

        ncbi_url = f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"
        log.info(f"  Gene ID: {gene_id}  →  {ncbi_url}")

        xml_text = ncbi_gene_xml(gene_id)
        time.sleep(SLEEP_API)

        if not xml_text:
            rows.append({
                "Gene": gene,
                "NCBI_Link": ncbi_url,
                "All_Related_Sequences": "",
                "New_Sequences": "",
                "Status": "XML fetch failed",
            })
            continue

        accessions = parse_related_sequences(xml_text)
        log.info(f"  Related sequences found: {len(accessions)}  {accessions}")

        existing = existing_accessions_for_gene(gene, OUT_DIR)
        new_accs = [a for a in accessions if not accession_exists(a, existing)]
        skipped = len(accessions) - len(new_accs)
        log.info(f"  Already present: {skipped}   To download: {len(new_accs)}")

        downloaded_accs = []
        n_ok = n_fail = n_short = 0

        for acc in new_accs:
            dest = os.path.join(OUT_DIR, f"{gene}_{acc.split('.')[0]}.fasta")
            log.info(f"    Downloading {acc} …")
            ok, bp = download_fasta(acc, dest)
            if ok:
                kb = os.path.getsize(dest) / 1024
                log.info(f"    Saved {os.path.basename(dest)}  ({bp} bp / {kb:.1f} KB)")
                downloaded_accs.append(acc)
                n_ok += 1
            elif bp > 0:
                n_short += 1
            else:
                log.warning(f"    FAILED: {acc}")
                n_fail += 1
            time.sleep(SLEEP_API)

        rows.append({
            "Gene": gene,
            "NCBI_Link": ncbi_url,
            "All_Related_Sequences": ", ".join(accessions),
            "New_Sequences": ", ".join(downloaded_accs),
            "Status": f"downloaded={n_ok} too_short={n_short} failed={n_fail} already_present={skipped}",
        })
        time.sleep(SLEEP_GENE)

    out_df = pd.DataFrame(rows, columns=["Gene", "NCBI_Link", "All_Related_Sequences", "New_Sequences", "Status"])
    out_df.to_csv(CSV_OUT, index=False)
    log.info(f"\nCSV saved: {CSV_OUT}")
    log.info("Done.")


if __name__ == "__main__":
    main()
