import os
import re
import time
import json
import csv
import winreg
import requests

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_TSV  = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..', '05_Testing', 'reference_check', 'origin_paper.tsv'))
OUT_DIR    = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..', '05_Testing', 'reference_check'))
OUT_TSV    = os.path.join(OUT_DIR, "origin_paper_with_refseq.tsv")
CACHE_FILE = os.path.join(OUT_DIR, "refseq_cache.json")

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLAUDE_URL  = "https://api.anthropic.com/v1/messages"
SLEEP_API   = 0.4

RE_PMID = re.compile(r"list_uids=(\d+)", re.IGNORECASE)


def get_api_key():
    key = os.environ.get("ANTHROPIC_API_KEY", "")
    if key:
        return key
    try:
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as reg:
            key, _ = winreg.QueryValueEx(reg, "ANTHROPIC_API_KEY")
            return key
    except Exception:
        return ""


def load_cache(path):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_cache(path, cache):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)


def fetch_pubmed_abstract(pmid):
    params = {"db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text"}
    try:
        r = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params, timeout=30)
        if r.ok:
            return r.text[:4000]
    except Exception as e:
        print(f"    Entrez error for {pmid}: {e}")
    return ""


def ask_claude_for_refseq(api_key, idbase, title, abstract):
    prompt = (
        f"You are a bioinformatics expert. Given the article title and abstract below, "
        f"identify the SPECIFIC reference sequence accession (GenBank/EMBL/ENA accession number, "
        f"e.g. M13792, X02994, U15905) used as the reference for numbering mutations in the {idbase} gene.\n\n"
        f"Article title: {title}\n\nAbstract:\n{abstract}\n\n"
        f"Rules:\n"
        f"- Return ONLY the accession number(s), comma-separated if multiple.\n"
        f"- If strongly inferable from field conventions, state it with a note like \"M13792 (inferred)\".\n"
        f"- If truly unknown, return exactly: Unknown\n"
        f"- No other text."
    )

    headers = {
        "x-api-key": api_key,
        "anthropic-version": "2023-06-01",
        "content-type": "application/json",
    }
    payload = {
        "model": "claude-sonnet-4-20250514",
        "max_tokens": 100,
        "messages": [{"role": "user", "content": prompt}],
    }

    try:
        r = requests.post(CLAUDE_URL, headers=headers, json=payload, timeout=30)
        r.raise_for_status()
        return r.json()["content"][0]["text"].strip()
    except Exception as e:
        print(f"    Claude API error: {e}")
        return "Unknown"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    api_key = get_api_key()
    if not api_key:
        print("ERROR: ANTHROPIC_API_KEY not found in environment or registry")
        return
    print(f"API key loaded (length {len(api_key)})")

    cache = load_cache(CACHE_FILE)

    rows = []
    with open(INPUT_TSV, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)

    print(f"Loaded {len(rows)} rows")

    unique_articles = {}
    for r in rows:
        key = (r["idbase"], r["article_title"])
        if key not in unique_articles:
            unique_articles[key] = r["link"]

    print(f"Unique (idbase, article) pairs: {len(unique_articles)}")

    for idx, ((idbase, title), link) in enumerate(unique_articles.items(), 1):
        cache_key = f"{idbase}||{title}"
        if cache_key in cache:
            continue

        pmid_m = RE_PMID.search(link) if link and link != "None" else None
        pmid   = pmid_m.group(1) if pmid_m else None

        print(f"  [{idx}/{len(unique_articles)}] {idbase}: {title[:60]}...")

        abstract = ""
        if pmid:
            abstract = fetch_pubmed_abstract(pmid)
            time.sleep(SLEEP_API)

        refseq = ask_claude_for_refseq(api_key, idbase, title, abstract)
        print(f"    -> {refseq}")
        cache[cache_key] = refseq
        save_cache(CACHE_FILE, cache)
        time.sleep(SLEEP_API)

    with open(OUT_TSV, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["idbase", "accessions", "article_title", "link", "reference_sequence"])
        for r in rows:
            cache_key = f"{r['idbase']}||{r['article_title']}"
            refseq = cache.get(cache_key, "Unknown")
            writer.writerow([r["idbase"], r["accessions"], r["article_title"], r["link"], refseq])

    print(f"\nOutput : {OUT_TSV}")
    print(f"Cache  : {CACHE_FILE}")


if __name__ == "__main__":
    main()
