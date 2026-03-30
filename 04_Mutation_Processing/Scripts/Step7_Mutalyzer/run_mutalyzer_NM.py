import argparse
import os
import csv
import json
import time
import requests
import urllib.parse
from pathlib import Path
from datetime import datetime

INPUT_TSV_TEMPLATE  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_input_NM_{source}.tsv"
OUT_TSV_TEMPLATE    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_results_NM_{source}.tsv"
CACHE_TEMPLATE      = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_cache_NM_{source}.jsonl"
LOG_TEMPLATE        = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Logs\mutalyzer_run_NM_{source}.log"

SOURCE_FASTA_DIRS = {
    "MANE":     r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\Mane_Select_NM",
    "IDRefseq": r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\IDRefseq_NM",
}

os.makedirs(r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer", exist_ok=True)
os.makedirs(r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Logs", exist_ok=True)

API_BASE   = "https://mutalyzer.nl/api"
SLEEP      = 0.5
MAX_RETRY  = 3
RETRY_WAIT = 10


def make_logger(log_path):
    def log(msg):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"{ts}  {msg}"
        print(line)
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    return log


def load_cache(path):
    cache = {}
    if not os.path.exists(path):
        return cache
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                cache[obj["key"]] = obj["value"]
            except Exception:
                pass
    return cache


def save_cache(path, key, value):
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps({"key": key, "value": value}, ensure_ascii=False) + "\n")


def build_refseq_index(directory):
    index = {}
    for p in Path(directory).iterdir():
        if p.suffix.lower() not in (".fasta", ".fa"):
            continue
        parts = p.stem.split("_", 1)
        if len(parts) == 2:
            gene = parts[0].upper()
            with open(p, encoding="utf-8") as f:
                header = f.readline().strip()
            index[gene] = header.lstrip(">")
    return index


def pick_mane(entries):
    for entry in entries:
        if "MANE_SELECT" in entry.get("selector", {}).get("tags", []):
            return entry.get("description", "")
    return entries[0].get("description", "") if entries else ""


def call_normalize(hgvs, log):
    encoded = urllib.parse.quote(hgvs, safe="")
    url = f"{API_BASE}/normalize/{encoded}"

    for attempt in range(1, MAX_RETRY + 1):
        try:
            r = requests.get(url, timeout=30, headers={"Accept": "application/json"})
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", RETRY_WAIT * attempt))
                log(f"  429 rate limit — sleeping {wait}s")
                time.sleep(wait)
                continue
            if r.status_code == 422:
                return {"_status": "invalid", "_errors": str(r.json())}
            if not r.ok:
                return {"_status": "error", "_errors": f"HTTP {r.status_code}: {r.text[:200]}"}

            data    = r.json()
            errors  = data.get("errors", [])
            infos   = data.get("infos", [])
            equiv   = data.get("equivalent_descriptions", {})
            protein = data.get("protein") or {}
            rna     = data.get("rna") or {}

            nc_hgvs = next(
                (e["description"] for e in equiv.get("g", []) if e.get("description", "").startswith("NC_")),
                ""
            )

            return {
                "_status":            "ok" if not errors else "warning",
                "_normalized":        data.get("normalized_description", ""),
                "_nc_hgvs":           nc_hgvs,
                "_c_hgvs":            pick_mane(equiv.get("c", [])),
                "_g_hgvs":            pick_mane(equiv.get("g", [])),
                "_r_hgvs":            rna.get("description", ""),
                "_p_hgvs":            protein.get("description", ""),
                "_protein_pos_first": str(protein["position_first"]) if "position_first" in protein else "",
                "_mane_select":       (data.get("tag") or {}).get("details", ""),
                "_mutalyzer_gene":    data.get("gene_id", ""),
                "_errors":            "; ".join(str(e) for e in errors),
                "_infos":             "; ".join(str(i) for i in infos),
            }
        except requests.RequestException as exc:
            log(f"  Attempt {attempt}/{MAX_RETRY} failed: {exc}")
            if attempt < MAX_RETRY:
                time.sleep(RETRY_WAIT * attempt)

    return {"_status": "unreachable", "_errors": "max retries exceeded"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", choices=["MANE", "IDRefseq"], default="MANE",
                        help="Which input TSV to run: MANE or IDRefseq")
    parser.add_argument("--clear-cache", action="store_true",
                        help="Delete the cache file before running to start fresh")
    args = parser.parse_args()

    input_tsv  = INPUT_TSV_TEMPLATE.format(source=args.source)
    out_tsv    = OUT_TSV_TEMPLATE.format(source=args.source)
    cache_path = CACHE_TEMPLATE.format(source=args.source)
    log_path   = LOG_TEMPLATE.format(source=args.source)
    fasta_dir  = SOURCE_FASTA_DIRS[args.source]
    log        = make_logger(log_path)

    if args.clear_cache and os.path.exists(cache_path):
        os.remove(cache_path)
        log(f"Cache cleared: {cache_path}")

    log("=" * 60)
    log(f"Mutalyzer 3 NM batch normalizer — source: {args.source}")
    log(f"Input : {input_tsv}")
    log("=" * 60)

    with open(input_tsv, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    log(f"Loaded {len(rows)} rows")

    refseq_index = build_refseq_index(fasta_dir)
    log(f"FASTA headers loaded for {len(refseq_index)} genes")

    cache    = load_cache(cache_path)
    log(f"Cache: {len(cache)} entries")

    unique_hgvs = list(dict.fromkeys(r["hgvs_input"] for r in rows if r.get("hgvs_input")))
    todo        = [h for h in unique_hgvs if h not in cache]
    log(f"Unique HGVS: {len(unique_hgvs)}  |  to fetch: {len(todo)}")

    for i, hgvs in enumerate(todo, 1):
        log(f"  [{i}/{len(todo)}] {hgvs}")
        result = call_normalize(hgvs, log)
        entry = {
            "status":            result.get("_status", "error"),
            "normalized":        result.get("_normalized", ""),
            "nc_hgvs":           result.get("_nc_hgvs", ""),
            "c_hgvs":            result.get("_c_hgvs", ""),
            "g_hgvs":            result.get("_g_hgvs", ""),
            "r_hgvs":            result.get("_r_hgvs", ""),
            "p_hgvs":            result.get("_p_hgvs", ""),
            "protein_pos_first": result.get("_protein_pos_first", ""),
            "mane_select":       result.get("_mane_select", ""),
            "mutalyzer_gene":    result.get("_mutalyzer_gene", ""),
            "errors":            result.get("_errors", ""),
            "infos":             result.get("_infos", ""),
        }
        cache[hgvs] = entry
        save_cache(cache_path, hgvs, entry)
        if entry["status"] == "unreachable":
            log("  API unreachable — stopping early")
            break
        time.sleep(SLEEP)

    out_rows = []
    for row in rows:
        hgvs  = row.get("hgvs_input", "")
        entry = cache.get(hgvs, {})
        gene  = row.get("gene", "").upper()
        out_rows.append({
            **row,
            "reference_coding":  refseq_index.get(gene, ""),
            "status":            entry.get("status",             "no_hgvs"),
            "normalized":        entry.get("normalized",         ""),
            "nc_hgvs":           entry.get("nc_hgvs",           ""),
            "c_hgvs":            entry.get("c_hgvs",            ""),
            "g_hgvs":            entry.get("g_hgvs",            ""),
            "r_hgvs":            entry.get("r_hgvs",            ""),
            "p_hgvs":            entry.get("p_hgvs",            ""),
            "protein_pos_first": entry.get("protein_pos_first", ""),
            "mane_select":       entry.get("mane_select",       ""),
            "mutalyzer_gene":    entry.get("mutalyzer_gene",    ""),
            "errors":            entry.get("errors",            ""),
            "infos":             entry.get("infos",             ""),
        })

    result_fields = ["reference_coding", "status", "normalized", "nc_hgvs", "c_hgvs",
                     "g_hgvs", "r_hgvs", "p_hgvs", "protein_pos_first",
                     "mane_select", "mutalyzer_gene", "errors", "infos"]
    fieldnames = list(rows[0].keys()) + result_fields

    with open(out_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)

    n_ok      = sum(1 for r in out_rows if r["status"] == "ok")
    n_warning = sum(1 for r in out_rows if r["status"] == "warning")
    n_invalid = sum(1 for r in out_rows if r["status"] in ("error", "invalid"))
    n_unreach = sum(1 for r in out_rows if r["status"] == "unreachable")
    n_mane    = sum(1 for r in out_rows if r.get("mane_select"))
    n_p       = sum(1 for r in out_rows if r.get("p_hgvs"))
    n_r       = sum(1 for r in out_rows if r.get("r_hgvs"))

    log("=" * 60)
    log(f"Results written : {out_tsv}")
    log(f"  ok            : {n_ok}")
    log(f"  warning       : {n_warning}")
    log(f"  error/invalid : {n_invalid}")
    log(f"  unreachable   : {n_unreach}")
    log(f"  MANE Select   : {n_mane}")
    log(f"  with p_hgvs   : {n_p}")
    log(f"  with r_hgvs   : {n_r}")
    log("=" * 60)
    log("Done.")


if __name__ == "__main__":
    main()
