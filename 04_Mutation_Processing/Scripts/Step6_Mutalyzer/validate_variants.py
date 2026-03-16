"""
validate_variants.py
====================
Author  : Generated for Assignment/Thesis – Mutation Extraction and Matching
Purpose : Takes the all_mutations.tsv produced by extract_mutations.py and
          validates / normalises every HGVS c. notation variant using two
          complementary REST APIs:

            1. Variant Validator  – https://rest.variantvalidator.org
               • Primary validator; returns normalised HGVS, RefSeq transcript,
                 genomic coordinates, and error messages.
               • Called one variant at a time (no official batch endpoint).
               • Rate limit: ~1 req/s recommended (no key required).

            2. Ensembl VEP        – https://rest.ensembl.org/vep/human/hgvs
               • Accepts POST with up to 200 variants per request  ← fast batch
               • Returns consequence terms, gene symbol, transcript impact.
               • No key required; rate limit: 15 req/s.

          The two results are merged into a single output CSV.

Pipeline
--------
    extract_mutations.py  →  all_mutations.tsv
                                    ↓
                         validate_variants.py
                                    ↓
                      validated_variants.csv   (main output)
                      vv_raw_cache.jsonl       (raw VV responses – resumable)
                      vep_raw_cache.jsonl      (raw VEP responses – resumable)
                      validation_log.txt

Usage
-----
    python validate_variants.py

Configuration
-------------
    Edit the CONSTANTS block below.

Dependencies
------------
    pip install requests pandas tqdm

Notes
-----
  • Only rows whose sysname contains a c. HGVS notation are sent to the APIs.
    Rows without a usable c. notation are still written to the output with
    vv_status = "no_hgvs".
  • Both caches (vv_raw_cache.jsonl and vep_raw_cache.jsonl) are written
    incrementally – the script is fully resumable after interruption.
  • VEP batches up to VEP_BATCH_SIZE variants per POST request.
  • If VariantValidator is unreachable the script falls back to VEP-only mode.
"""

# ──────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ──────────────────────────────────────────────────────────────────────────────
import os
import re
import sys
import json
import time
import logging
import requests
import pandas as pd
from pathlib import Path
from datetime import datetime

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS  – edit these as needed
# ──────────────────────────────────────────────────────────────────────────────

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
THESIS_DIR  = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
INPUT_TSV   = os.path.join(THESIS_DIR, '04_Mutation_Processing', 'Output', 'all_mutations.tsv')
OUTPUT_DIR  = os.path.join(THESIS_DIR, '04_Mutation_Processing', 'Output')

# Gene → primary RefSeq transcript mapping.
# Variant Validator needs a transcript to anchor c. coordinates.
# Add every gene in your dataset. If a gene is absent, VV will try to
# resolve it from the c. notation alone (less reliable).
TRANSCRIPT_MAP: dict[str, str] = {
    "ADA":    "NM_000022.4",
    "AICDA":  "NM_020661.3",
    "ATM":    "NM_000051.4",
    "BRCA1":  "NM_007294.4",
    "BRCA2":  "NM_000059.4",
    "CFTR":   "NM_000492.4",
    "MLH1":   "NM_000249.4",
    "MSH2":   "NM_000251.3",
    "MSH6":   "NM_000179.3",
    "PMS2":   "NM_000535.7",
    "TP53":   "NM_000546.6",
    # ← add more genes here
}

# Genome build used by both APIs ("GRCh38" or "GRCh37")
GENOME_BUILD = "GRCh38"

# VEP batch size – max 200 per POST request
VEP_BATCH_SIZE = 200

# Seconds to sleep between VariantValidator calls (respect 1 req/s limit)
VV_SLEEP = 1.1

# Seconds to sleep between VEP batch POSTs
VEP_SLEEP = 0.5

# HTTP retry settings
MAX_RETRIES    = 3
RETRY_BACKOFF  = 5   # seconds

# ──────────────────────────────────────────────────────────────────────────────
# DERIVED PATHS  (auto-set from OUTPUT_DIR)
# ──────────────────────────────────────────────────────────────────────────────

LOG_PATH       = os.path.join(OUTPUT_DIR, "validation_log.txt")
OUT_CSV        = os.path.join(OUTPUT_DIR, "validated_variants.csv")
VV_CACHE_PATH  = os.path.join(OUTPUT_DIR, "vv_raw_cache.jsonl")
VEP_CACHE_PATH = os.path.join(OUTPUT_DIR, "vep_raw_cache.jsonl")

# ──────────────────────────────────────────────────────────────────────────────
# LOGGING
# ──────────────────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(LOG_PATH, encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# REGEX
# ──────────────────────────────────────────────────────────────────────────────

# Matches c. HGVS: c.976C>T  c.1235delA  c.100_101insAT  c.-15C>T  c.*3A>G
RE_HGVS_C = re.compile(r"(c\.[-*]?[0-9_+\-*]+[A-Za-z>_]+[A-Za-z0-9]*)", re.IGNORECASE)

# ──────────────────────────────────────────────────────────────────────────────
# CACHE HELPERS  (JSONL – one JSON object per line, append-only)
# ──────────────────────────────────────────────────────────────────────────────

def load_cache(path: str) -> dict:
    """Load a JSONL cache file into a dict keyed by the 'key' field."""
    cache = {}
    if not os.path.exists(path):
        return cache
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                cache[obj["key"]] = obj["value"]
            except (json.JSONDecodeError, KeyError):
                pass
    return cache


def append_cache(path: str, key: str, value) -> None:
    """Append one entry to a JSONL cache file (atomic-ish, no rewrite)."""
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps({"key": key, "value": value}, ensure_ascii=False) + "\n")


# ──────────────────────────────────────────────────────────────────────────────
# HTTP HELPERS  (with retry + back-off)
# ──────────────────────────────────────────────────────────────────────────────

def http_get(url: str, params: dict | None = None,
             headers: dict | None = None) -> requests.Response | None:
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=30)
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", RETRY_BACKOFF * attempt))
                log.warning(f"  429 rate-limit – sleeping {wait}s")
                time.sleep(wait)
                continue
            return r
        except requests.RequestException as exc:
            log.warning(f"  GET attempt {attempt}/{MAX_RETRIES} failed: {exc}")
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_BACKOFF * attempt)
    return None


def http_post(url: str, payload: dict,
              headers: dict | None = None) -> requests.Response | None:
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.post(url, json=payload, headers=headers, timeout=60)
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", RETRY_BACKOFF * attempt))
                log.warning(f"  429 rate-limit – sleeping {wait}s")
                time.sleep(wait)
                continue
            return r
        except requests.RequestException as exc:
            log.warning(f"  POST attempt {attempt}/{MAX_RETRIES} failed: {exc}")
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_BACKOFF * attempt)
    return None


# ──────────────────────────────────────────────────────────────────────────────
# HGVS HELPERS
# ──────────────────────────────────────────────────────────────────────────────

def extract_hgvs_c(sysname: str) -> str | None:
    """
    Pull the first c. notation from a sysname string.
    "g.4037C>T; c.967C>T; p.S291L"  →  "c.967C>T"
    """
    m = RE_HGVS_C.search(sysname)
    return m.group(1) if m else None


def build_full_hgvs(gene: str, hgvs_c: str) -> str:
    """
    Prepend the configured RefSeq transcript to the c. notation.
    Falls back to bare c. string if no transcript is mapped.
    e.g.  "NM_000022.4:c.976C>T"
    """
    tx = TRANSCRIPT_MAP.get(gene.upper(), "")
    return f"{tx}:{hgvs_c}" if tx else hgvs_c


# ──────────────────────────────────────────────────────────────────────────────
# VARIANT VALIDATOR  (single-variant GET)
# ──────────────────────────────────────────────────────────────────────────────

VV_BASE = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"

def call_variant_validator(full_hgvs: str, build: str = GENOME_BUILD) -> dict:
    """
    Validate one HGVS variant via VariantValidator REST API.

    Endpoint:
        GET /VariantValidator/variantvalidator/{build}/{variant}/all

    Returns dict with:
        vv_status       : "ok" | "warning" | "error" | "unreachable"
        vv_hgvs_norm    : normalised HGVS c. notation
        vv_gene_symbol  : HGNC gene symbol
        vv_transcript   : RefSeq transcript used
        vv_g_hgvs       : genomic g. HGVS notation
        vv_error        : validation warning / error text
    """
    empty = {
        "vv_status": "unreachable", "vv_hgvs_norm": "", "vv_gene_symbol": "",
        "vv_transcript": "", "vv_g_hgvs": "", "vv_error": "API unreachable",
    }

    # URL-encode colons and > signs in the HGVS string
    encoded = requests.utils.quote(full_hgvs, safe="")
    url = f"{VV_BASE}/{build}/{encoded}/all"

    resp = http_get(url, headers={"Accept": "application/json"})
    if resp is None:
        return empty
    if not resp.ok:
        return {**empty, "vv_status": "error",
                "vv_error": f"HTTP {resp.status_code}: {resp.text[:200]}"}

    try:
        data = resp.json()
    except ValueError:
        return {**empty, "vv_status": "error", "vv_error": "JSON parse error"}

    # VV response: top-level keys are variant IDs + "flag"
    result_keys = [k for k in data if k != "flag"]
    if not result_keys:
        return {**empty, "vv_status": "error",
                "vv_error": f"empty result (flag={data.get('flag','')})"}

    vdata = data[result_keys[0]]

    warnings  = vdata.get("validation_warnings", [])
    error_msg = "; ".join(warnings) if warnings else ""

    # Normalised transcript variant
    hgvs_norm  = vdata.get("hgvs_transcript_variant", "")
    gene_info  = vdata.get("gene_ids", {})
    gene_sym   = gene_info.get("hgnc_symbol", "")
    transcript = hgvs_norm.split(":")[0] if ":" in hgvs_norm else ""

    # Genomic HGVS from primary_assembly_loci
    g_hgvs = ""
    for build_key in [build, build.lower(), "grch38", "grch37"]:
        loci = vdata.get("primary_assembly_loci", {}).get(build_key, {})
        if loci:
            g_hgvs = loci.get("hgvs_genomic_description", "")
            break

    return {
        "vv_status":      "warning" if error_msg else "ok",
        "vv_hgvs_norm":   hgvs_norm,
        "vv_gene_symbol": gene_sym,
        "vv_transcript":  transcript,
        "vv_g_hgvs":      g_hgvs,
        "vv_error":       error_msg,
    }


# ──────────────────────────────────────────────────────────────────────────────
# ENSEMBL VEP  (batch POST, up to 200 variants per request)
# ──────────────────────────────────────────────────────────────────────────────

VEP_URL = "https://rest.ensembl.org/vep/human/hgvs"

def call_vep_batch(hgvs_list: list[str]) -> dict[str, dict]:
    """
    Annotate a batch of HGVS notations via Ensembl VEP POST endpoint.

    Returns dict mapping input HGVS → result dict with keys:
        vep_consequence     : most severe consequence term(s)
        vep_gene_symbol     : gene symbol
        vep_transcript_id   : Ensembl transcript ID
        vep_impact          : HIGH / MODERATE / LOW / MODIFIER
        vep_biotype         : transcript biotype
        vep_sift            : SIFT prediction(score)
        vep_polyphen        : PolyPhen prediction(score)
        vep_error           : empty if ok
    """
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = {"hgvs_notations": hgvs_list}

    empty = {
        "vep_consequence": "", "vep_gene_symbol": "", "vep_transcript_id": "",
        "vep_impact": "", "vep_biotype": "", "vep_sift": "",
        "vep_polyphen": "", "vep_error": "no response",
    }
    results = {h: dict(empty) for h in hgvs_list}

    resp = http_post(VEP_URL, payload, headers=headers)
    if resp is None:
        return results

    if not resp.ok:
        err = f"HTTP {resp.status_code}: {resp.text[:200]}"
        log.warning(f"  VEP batch error: {err}")
        for h in hgvs_list:
            results[h]["vep_error"] = err
        return results

    try:
        data = resp.json()
    except ValueError:
        for h in hgvs_list:
            results[h]["vep_error"] = "JSON parse error"
        return results

    for entry in data:
        input_hgvs = entry.get("input", "")
        if not input_hgvs:
            continue

        # Pick most severe transcript consequence (first in list)
        tc   = entry.get("transcript_consequences") or \
               entry.get("intergenic_consequences") or []
        best = tc[0] if tc else {}

        csq  = ", ".join(best.get("consequence_terms", []))
        sift = best.get("sift_prediction", "")
        sift_score = best.get("sift_score", "")
        pp   = best.get("polyphen_prediction", "")
        pp_score = best.get("polyphen_score", "")

        results[input_hgvs] = {
            "vep_consequence":   csq,
            "vep_gene_symbol":   best.get("gene_symbol", entry.get("gene_symbol", "")),
            "vep_transcript_id": best.get("transcript_id", ""),
            "vep_impact":        best.get("impact", ""),
            "vep_biotype":       best.get("biotype", ""),
            "vep_sift":          f"{sift}({sift_score})" if sift else "",
            "vep_polyphen":      f"{pp}({pp_score})" if pp else "",
            "vep_error":         "",
        }

    return results


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    start_time = datetime.now()
    log.info("=" * 70)
    log.info("Variant Validator + Ensembl VEP Batch Annotator – starting")
    log.info(f"Input TSV    : {INPUT_TSV}")
    log.info(f"Output dir   : {OUTPUT_DIR}")
    log.info(f"Genome build : {GENOME_BUILD}")
    log.info(f"Started at   : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("=" * 70)

    # ── 1. Load mutations TSV ─────────────────────────────────────────────────
    if not os.path.isfile(INPUT_TSV):
        log.error(f"Input TSV not found: {INPUT_TSV}")
        sys.exit(1)

    df = pd.read_csv(INPUT_TSV, sep="\t", dtype=str).fillna("")
    log.info(f"Loaded {len(df)} rows")

    # Extract c. HGVS from sysname column
    df["hgvs_c"]    = df["sysname"].apply(extract_hgvs_c)
    df["full_hgvs"] = df.apply(
        lambda r: build_full_hgvs(r["gene"], r["hgvs_c"])
        if r["hgvs_c"] else "",
        axis=1,
    )

    n_with    = (df["hgvs_c"] != "").sum()
    n_without = (df["hgvs_c"] == "").sum()
    log.info(f"Rows with c. HGVS    : {n_with}")
    log.info(f"Rows without c. HGVS : {n_without}  (passed through with vv_status=no_hgvs)")
    log.info("")

    # ── 2. Load caches ────────────────────────────────────────────────────────
    vv_cache  = load_cache(VV_CACHE_PATH)
    vep_cache = load_cache(VEP_CACHE_PATH)
    log.info(f"Cache: {len(vv_cache)} VV, {len(vep_cache)} VEP entries already stored")

    # ── 3. Variant Validator – one call per unique full_hgvs ──────────────────
    unique_hgvs    = [h for h in df["full_hgvs"].unique() if h]
    todo_vv        = [h for h in unique_hgvs if h not in vv_cache]
    vv_unreachable = False

    log.info(f"VV: {len(unique_hgvs)} unique HGVS  |  {len(todo_vv)} to fetch")

    if todo_vv:
        it = tqdm(todo_vv, desc="VariantValidator", unit="var") if HAS_TQDM else todo_vv
        for full_hgvs in it:
            log.info(f"  VV ← {full_hgvs}")
            result = call_variant_validator(full_hgvs)
            vv_cache[full_hgvs] = result
            append_cache(VV_CACHE_PATH, full_hgvs, result)

            if result["vv_status"] == "unreachable":
                log.warning("  VariantValidator unreachable – switching to VEP-only mode")
                vv_unreachable = True
                break

            time.sleep(VV_SLEEP)

    # ── 4. Ensembl VEP – batched POST ─────────────────────────────────────────
    todo_vep = [h for h in unique_hgvs if h not in vep_cache]
    log.info(f"VEP: {len(unique_hgvs)} unique HGVS  |  {len(todo_vep)} to fetch")

    if todo_vep:
        batches = [todo_vep[i:i + VEP_BATCH_SIZE]
                   for i in range(0, len(todo_vep), VEP_BATCH_SIZE)]
        it = tqdm(batches, desc="VEP batches", unit="batch") if HAS_TQDM else batches

        for batch in it:
            log.info(f"  VEP POST: {len(batch)} variants")
            batch_res = call_vep_batch(batch)
            for hgvs, res in batch_res.items():
                vep_cache[hgvs] = res
                append_cache(VEP_CACHE_PATH, hgvs, res)
            time.sleep(VEP_SLEEP)

    # ── 5. Merge results into DataFrame ──────────────────────────────────────
    VV_EMPTY = {
        "vv_status": "no_hgvs", "vv_hgvs_norm": "", "vv_gene_symbol": "",
        "vv_transcript": "", "vv_g_hgvs": "", "vv_error": "",
    }
    VEP_EMPTY = {
        "vep_consequence": "", "vep_gene_symbol": "", "vep_transcript_id": "",
        "vep_impact": "", "vep_biotype": "", "vep_sift": "",
        "vep_polyphen": "", "vep_error": "",
    }

    vv_rows  = [vv_cache.get(r["full_hgvs"],  VV_EMPTY)  if r["full_hgvs"] else VV_EMPTY
                for _, r in df.iterrows()]
    vep_rows = [vep_cache.get(r["full_hgvs"], VEP_EMPTY) if r["full_hgvs"] else VEP_EMPTY
                for _, r in df.iterrows()]

    df_out = pd.concat([
        df.reset_index(drop=True),
        pd.DataFrame(vv_rows),
        pd.DataFrame(vep_rows),
    ], axis=1)

    # ── 6. Write output CSV ───────────────────────────────────────────────────
    df_out.to_csv(OUT_CSV, index=False)

    log.info("")
    log.info("=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)
    log.info(f"  Total rows              : {len(df_out)}")
    log.info(f"  Rows with HGVS          : {n_with}")
    log.info(f"  VV ok                   : {(df_out['vv_status'] == 'ok').sum()}")
    log.info(f"  VV warnings             : {(df_out['vv_status'] == 'warning').sum()}")
    log.info(f"  VV errors               : {(df_out['vv_status'] == 'error').sum()}")
    log.info(f"  VV unreachable/skipped  : {(df_out['vv_status'] == 'unreachable').sum()}")
    log.info(f"  VEP annotated           : {(df_out['vep_consequence'] != '').sum()}")
    log.info(f"  Output CSV              : {OUT_CSV}")
    log.info(f"  VV cache                : {VV_CACHE_PATH}")
    log.info(f"  VEP cache               : {VEP_CACHE_PATH}")
    log.info(f"  Log                     : {LOG_PATH}")
    elapsed = datetime.now() - start_time
    log.info(f"  Elapsed time            : {str(elapsed).split('.')[0]}")
    log.info("=" * 70)
    log.info("Done.")


if __name__ == "__main__":
    main()
