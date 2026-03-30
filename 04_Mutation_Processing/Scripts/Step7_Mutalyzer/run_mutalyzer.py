import os
import csv
import json
import time
import requests
import urllib.parse
from pathlib import Path
from datetime import datetime

INPUT_TSV    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_input.tsv"
OUT_TSV      = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_results.tsv"
CACHE_G_JSONL = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_cache.jsonl"
CACHE_C_JSONL = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_cache_c.jsonl"
LOG_PATH     = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Logs\mutalyzer_run.log"

os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)

API_BASE   = "https://mutalyzer.nl/api"
SLEEP      = 0.5
MAX_RETRY  = 3
RETRY_WAIT = 10


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"{ts}  {msg}"
    print(line)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line + "\n")


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


def pick_mane(entries):
    """Return MANE Select description from equivalent_descriptions list.
    Falls back to first entry if no MANE Select present."""
    for entry in entries:
        tags = entry.get("selector", {}).get("tags", [])
        if "MANE_SELECT" in tags:
            return entry.get("description", "")
    if entries:
        return entries[0].get("description", "")
    return ""


def call_normalize(hgvs):
    """
    Call /api/normalize with any HGVS string.
    When input is g. format: returns normalized, nc_hgvs, c_hgvs, status.
    When input is c. format: additionally returns p_hgvs, r_hgvs,
      protein_pos_first, mane_select, mutalyzer_gene.
    Returns a flat dict with all fields (empty string for unavailable ones).
    """
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
            tag     = data.get("tag") or {}

            nc_hgvs = ""
            for entry in equiv.get("g", []):
                if entry.get("description", "").startswith("NC_"):
                    nc_hgvs = entry["description"]
                    break

            return {
                "_status":           "ok" if not errors else "warning",
                "_normalized":       data.get("normalized_description", ""),
                "_nc_hgvs":          nc_hgvs,
                "_c_hgvs":           pick_mane(equiv.get("c", [])),
                "_r_hgvs":           rna.get("description", ""),
                "_p_hgvs":           protein.get("description", ""),
                "_protein_pos_first": str(protein["position_first"]) if "position_first" in protein else "",
                "_mane_select":      tag.get("details", ""),
                "_mutalyzer_gene":   data.get("gene_id", ""),
                "_errors":           "; ".join(str(e) for e in errors) if errors else "",
                "_infos":            "; ".join(str(i) for i in infos) if infos else "",
            }
        except requests.RequestException as exc:
            log(f"  Attempt {attempt}/{MAX_RETRY} failed: {exc}")
            if attempt < MAX_RETRY:
                time.sleep(RETRY_WAIT * attempt)

    return {"_status": "unreachable", "_errors": "max retries exceeded"}


EMPTY = {
    "status": "no_hgvs", "normalized": "", "nc_hgvs": "", "c_hgvs": "",
    "r_hgvs": "", "p_hgvs": "", "protein_pos_first": "",
    "mane_select": "", "mutalyzer_gene": "", "errors": "", "infos": "",
}


def main():
    log("=" * 60)
    log("Mutalyzer 3 batch normalizer — starting")
    log(f"Input : {INPUT_TSV}")
    log("=" * 60)

    with open(INPUT_TSV, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    log(f"Loaded {len(rows)} rows")

    # ── Pass 1: g. → c. (normalize genomic coordinates) ─────────────────────
    cache_g = load_cache(CACHE_G_JSONL)
    log(f"Pass 1 cache: {len(cache_g)} entries")

    unique_g = list(dict.fromkeys(r["hgvs_input"] for r in rows if r.get("hgvs_input")))
    todo_g   = [h for h in unique_g if h not in cache_g]
    log(f"Pass 1 — unique g. HGVS: {len(unique_g)}  |  to fetch: {len(todo_g)}")

    for i, hgvs in enumerate(todo_g, 1):
        log(f"  P1 [{i}/{len(todo_g)}] {hgvs}")
        result = call_normalize(hgvs)
        entry  = {
            "status":    result.get("_status", "error"),
            "normalized": result.get("_normalized", ""),
            "nc_hgvs":   result.get("_nc_hgvs", ""),
            "c_hgvs":    result.get("_c_hgvs", ""),
            "errors":    result.get("_errors", ""),
            "infos":     result.get("_infos", ""),
        }
        cache_g[hgvs] = entry
        save_cache(CACHE_G_JSONL, hgvs, entry)
        if entry["status"] == "unreachable":
            log("  API unreachable — stopping early")
            break
        time.sleep(SLEEP)

    # ── Pass 2: c. → p./r./gene/mane (re-submit c_hgvs to get protein info) ─
    cache_c = load_cache(CACHE_C_JSONL)
    log(f"Pass 2 cache: {len(cache_c)} entries")

    unique_c = list(dict.fromkeys(
        cache_g[h]["c_hgvs"]
        for h in unique_g
        if h in cache_g and cache_g[h].get("c_hgvs")
    ))
    todo_c = [c for c in unique_c if c not in cache_c]
    log(f"Pass 2 — unique c. HGVS: {len(unique_c)}  |  to fetch: {len(todo_c)}")

    for i, c_hgvs in enumerate(todo_c, 1):
        log(f"  P2 [{i}/{len(todo_c)}] {c_hgvs}")
        result = call_normalize(c_hgvs)
        entry  = {
            "r_hgvs":            result.get("_r_hgvs", ""),
            "p_hgvs":            result.get("_p_hgvs", ""),
            "protein_pos_first": result.get("_protein_pos_first", ""),
            "mane_select":       result.get("_mane_select", ""),
            "mutalyzer_gene":    result.get("_mutalyzer_gene", ""),
        }
        cache_c[c_hgvs] = entry
        save_cache(CACHE_C_JSONL, c_hgvs, entry)
        if result.get("_status") == "unreachable":
            log("  API unreachable — stopping early")
            break
        time.sleep(SLEEP)

    # ── Merge and write output ────────────────────────────────────────────────
    empty_c = {"r_hgvs": "", "p_hgvs": "", "protein_pos_first": "", "mane_select": "", "mutalyzer_gene": ""}

    out_rows = []
    for row in rows:
        g_hgvs  = row.get("hgvs_input", "")
        g_entry = cache_g.get(g_hgvs, {}) if g_hgvs else {}
        c_hgvs  = g_entry.get("c_hgvs", "")
        c_entry = cache_c.get(c_hgvs, empty_c) if c_hgvs else empty_c

        merged = {
            "status":            g_entry.get("status",     EMPTY["status"]),
            "normalized":        g_entry.get("normalized", ""),
            "nc_hgvs":           g_entry.get("nc_hgvs",   ""),
            "c_hgvs":            c_hgvs,
            "r_hgvs":            c_entry.get("r_hgvs",            ""),
            "p_hgvs":            c_entry.get("p_hgvs",            ""),
            "protein_pos_first": c_entry.get("protein_pos_first", ""),
            "mane_select":       c_entry.get("mane_select",       ""),
            "mutalyzer_gene":    c_entry.get("mutalyzer_gene",    ""),
            "errors":            g_entry.get("errors",    ""),
            "infos":             g_entry.get("infos",     ""),
        }
        out_rows.append({**row, **merged})

    result_fields = ["status", "normalized", "nc_hgvs", "c_hgvs", "r_hgvs", "p_hgvs",
                     "protein_pos_first", "mane_select", "mutalyzer_gene", "errors", "infos"]
    fieldnames = list(rows[0].keys()) + result_fields

    with open(OUT_TSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)

    n_ok      = sum(1 for r in out_rows if r["status"] == "ok")
    n_warning = sum(1 for r in out_rows if r["status"] == "warning")
    n_invalid = sum(1 for r in out_rows if r["status"] in ("error", "invalid"))
    n_unreach = sum(1 for r in out_rows if r["status"] == "unreachable")
    n_mane    = sum(1 for r in out_rows if r.get("mane_select") == "MANE Select")
    n_p       = sum(1 for r in out_rows if r.get("p_hgvs"))
    n_r       = sum(1 for r in out_rows if r.get("r_hgvs"))

    log("=" * 60)
    log(f"Results written : {OUT_TSV}")
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
