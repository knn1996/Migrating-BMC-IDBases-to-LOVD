"""
vcf_to_mutalyzer.py
===================
Reads the VCF produced by bed_to_vcf.py, converts each variant to an HGVS
genomic notation, and queries the Mutalyzer 3 REST API
(https://mutalyzer.nl/api/normalize/{description}) for each one.

Results are written to two TSV files:
  mutalyzer_results.tsv  — full results (g. / c. / p. / errors)
  mutalyzer_errors.tsv   — only rows where Mutalyzer returned errors/warnings

Usage:
  python3 vcf_to_mutalyzer.py [--workers N] [--delay S] [--no-progress]

Options:
  --workers N    parallel threads for API calls (default: 4, max: 8)
  --delay S      seconds to sleep between requests per thread (default: 0.5)
  --no-progress  suppress the progress bar
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

# ── Configuration ──────────────────────────────────────────────────────────────

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
VCF_PATH    = os.path.join(_SCRIPT_DIR, "all_variants_hg38.vcf")
OUT_PATH    = os.path.join(_SCRIPT_DIR, "mutalyzer_results.tsv")
ERR_PATH    = os.path.join(_SCRIPT_DIR, "mutalyzer_errors.tsv")

API_BASE    = "https://mutalyzer.nl/api"
TIMEOUT_SEC = 30   # per-request HTTP timeout (seconds)
MAX_RETRIES = 3    # retries on transient network errors
RETRY_DELAY = 5.0  # seconds to wait between retries

# ── hg38 chromosome → RefSeq NC_ accession (GRCh38.p14) ──────────────────────

CHR_TO_NC = {
    "chr1":  "NC_000001.11", "chr2":  "NC_000002.12", "chr3":  "NC_000003.12",
    "chr4":  "NC_000004.12", "chr5":  "NC_000005.10", "chr6":  "NC_000006.12",
    "chr7":  "NC_000007.14", "chr8":  "NC_000008.11", "chr9":  "NC_000009.12",
    "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9",  "chr15": "NC_000015.10",
    "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
    "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
    "chr22": "NC_000022.11", "chrX":  "NC_000023.11", "chrY":  "NC_000024.10",
    "chrM":  "NC_012920.1",
}

# Matches a genomic range in a BED name field, e.g. "g.12345_12350del"
RE_G_RANGE = re.compile(r"g\.(\d+)_(\d+)(?:del|dup)", re.IGNORECASE)

# Extracts the GENE value from the VCF INFO column
RE_GENE = re.compile(r"GENE=([^;]+)")

# ── TSV output columns ─────────────────────────────────────────────────────────

TSV_HEADER = "\t".join([
    "gene", "vcf_name", "hgvs_g",
    "normalized_g", "mane_c", "mane_p",
    "all_c", "all_p", "errors",
])

# ── HGVS builder ──────────────────────────────────────────────────────────────

def vcf_to_hgvs(chrom, pos, ref, alt, name):
    """
    Convert a VCF row to an HGVS genomic string anchored on the NC_ accession
    for the given chromosome.  Returns None if the chromosome is not mapped.
    """
    nc = CHR_TO_NC.get(chrom)
    if not nc:
        return None

    # SNV: single-base substitution
    if len(ref) == 1 and len(alt) == 1 and ref not in "N." and alt not in ("N", ".", "DUP"):
        return f"{nc}:g.{pos}{ref}>{alt}"

    # Range deletion / duplication encoded with sentinel ALT values
    if alt in (".", "DUP"):
        m = RE_G_RANGE.search(name)
        span = int(m.group(2)) - int(m.group(1)) if m else 0
        op   = "del" if alt == "." else "dup"
        return f"{nc}:g.{pos}_{pos + span}{op}" if span else f"{nc}:g.{pos}{op}"

    # Insertion / duplication with sequence: REF=N, ALT=Nxxx
    if ref == "N" and alt.startswith("N") and len(alt) > 1:
        ins = alt[1:]
        if "dup" in name.lower():
            return (
                f"{nc}:g.{pos}_{pos + len(ins) - 1}dup"
                if len(ins) > 1
                else f"{nc}:g.{pos}dup"
            )
        return f"{nc}:g.{pos - 1}_{pos}ins{ins}"

    # Deletion with sequence: REF=Nxxx, ALT=N
    if ref.startswith("N") and len(ref) > 1 and alt == "N":
        d = ref[1:]
        return (
            f"{nc}:g.{pos + 1}del"
            if len(d) == 1
            else f"{nc}:g.{pos + 1}_{pos + len(d)}del"
        )

    # Delins: REF=N, ALT=plain sequence
    if ref == "N" and alt not in ("N", ".", "DUP"):
        return f"{nc}:g.{pos}delins{alt}"

    # Fallback: express as delins spanning the reference allele
    return f"{nc}:g.{pos}_{pos + max(len(ref), 1) - 1}del{ref}ins{alt}"


# ── Mutalyzer 3 API ───────────────────────────────────────────────────────────

def call_mutalyzer(hgvs_desc, delay=0.5):
    """
    Call GET /api/normalize/{description} and return a normalised result dict.

    Keys returned:
      normalized   (str | None)  canonical g. description from Mutalyzer
      c_notations  (list[str])   all equivalent c. descriptions
      p_notations  (list[str])   all equivalent p. descriptions
      mane_c       (str | None)  MANE Select c. description
      mane_p       (str | None)  MANE Select p. description
      errors       (list[str])   error / warning messages
    """
    result = {
        "normalized":  None,
        "c_notations": [],
        "p_notations": [],
        "mane_c":      None,
        "mane_p":      None,
        "errors":      [],
    }

    encoded = urllib.parse.quote(hgvs_desc, safe="")
    url     = f"{API_BASE}/normalize/{encoded}"
    req     = urllib.request.Request(url, headers={"Accept": "application/json"})

    data = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            time.sleep(delay)
            with urllib.request.urlopen(req, timeout=TIMEOUT_SEC) as resp:
                data = json.loads(resp.read().decode())
            break  # success — exit retry loop
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                result["errors"].append("HTTP 404 – description not recognised by Mutalyzer")
                return result
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY)
                continue
            body = exc.read().decode(errors="replace")
            result["errors"].append(f"HTTP {exc.code}: {body[:200]}")
            return result
        except Exception as exc:
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY)
                continue
            result["errors"].append(f"Request error: {exc}")
            return result

    if data is None:
        return result

    # Normalised g. description (fall back to corrected if normalised absent)
    result["normalized"] = (
        data.get("normalized_description") or data.get("corrected_description")
    )

    # Collect any errors / warnings reported by Mutalyzer
    for severity in ("errors", "infos", "warnings"):
        for item in data.get(severity, []):
            msg = (
                item.get("details") or item.get("description") or json.dumps(item)
                if isinstance(item, dict)
                else str(item)
            )
            result["errors"].append(f"[{severity}] {msg}")

    # Equivalent c. and p. descriptions; flag the MANE Select entry
    equiv = data.get("equivalent_descriptions", {})

    for entry in equiv.get("c", []):
        desc = entry.get("description", "")
        result["c_notations"].append(desc)
        tag = entry.get("tag", {})
        if isinstance(tag, dict) and "MANE Select" in tag.get("details", ""):
            result["mane_c"] = desc

    for entry in equiv.get("p", []):
        desc = entry.get("description", "")
        result["p_notations"].append(desc)
        tag = entry.get("tag", {})
        if isinstance(tag, dict) and "MANE Select" in tag.get("details", ""):
            result["mane_p"] = desc

    return result


# ── VCF reader ────────────────────────────────────────────────────────────────

def read_vcf(vcf_path):
    """
    Parse the VCF and return (rows, n_skipped).

    rows      — list of (gene, vcf_name, hgvs_g) tuples
    n_skipped — rows dropped because the chromosome had no NC_ mapping
    """
    rows    = []
    skipped = 0

    with open(vcf_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, name, ref, alt, info = (
                parts[0], int(parts[1]), parts[2], parts[3], parts[4], parts[7]
            )

            m    = RE_GENE.search(info)
            gene = m.group(1) if m else ""
            hgvs = vcf_to_hgvs(chrom, pos, ref, alt, name)

            if hgvs:
                rows.append((gene, name, hgvs))
            else:
                skipped += 1

    return rows, skipped


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="VCF → Mutalyzer 3 API")
    parser.add_argument("--workers",     type=int,   default=4,
                        help="Parallel threads, 1–8 (default: 4)")
    parser.add_argument("--delay",       type=float, default=0.5,
                        help="Seconds between requests per thread (default: 0.5)")
    parser.add_argument("--no-progress", action="store_true",
                        help="Suppress the progress bar")
    args = parser.parse_args()

    workers       = min(max(args.workers, 1), 8)
    delay         = max(args.delay, 0.1)
    show_progress = not args.no_progress

    print(f"VCF → Mutalyzer 3 API  |  {VCF_PATH}")
    print(f"API    : {API_BASE}")
    print(f"Workers: {workers}   Delay: {delay}s/thread\n")

    rows, n_skipped = read_vcf(VCF_PATH)
    total = len(rows)
    print(f"Variants loaded : {total}  ({n_skipped} skipped – unmapped chromosome)\n")

    n_errors = 0
    done     = 0

    def process(item):
        gene, name, hgvs = item
        return gene, name, hgvs, call_mutalyzer(hgvs, delay=delay)

    with open(OUT_PATH, "w", encoding="utf-8") as out_fh, \
         open(ERR_PATH, "w", encoding="utf-8") as err_fh:

        out_fh.write(TSV_HEADER + "\n")
        err_fh.write(TSV_HEADER + "\n")

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(process, row): row for row in rows}
            for future in as_completed(futures):
                gene, name, hgvs, res = future.result()
                done += 1

                row_str = "\t".join([
                    gene, name, hgvs,
                    res["normalized"]          or "",
                    res["mane_c"]              or "",
                    res["mane_p"]              or "",
                    "; ".join(res["c_notations"]),
                    "; ".join(res["p_notations"]),
                    "; ".join(res["errors"]),
                ])
                out_fh.write(row_str + "\n")

                if res["errors"]:
                    err_fh.write(row_str + "\n")
                    n_errors += 1

                if show_progress:
                    pct = done / total * 100
                    bar = "#" * int(pct / 2)
                    sys.stdout.write(f"\r  [{bar:<50}] {done}/{total}  ({pct:.1f}%)  ")
                    sys.stdout.flush()

    if show_progress:
        print()

    print(f"\nDone.")
    print(f"  Results → {OUT_PATH}  ({total} rows)")
    print(f"  Errors  → {ERR_PATH}  ({n_errors} rows)")
    print()
    print("Output columns:")
    print("  gene         gene symbol from VCF INFO")
    print("  vcf_name     variant ID from the VCF NAME field")
    print("  hgvs_g       input genomic HGVS sent to Mutalyzer (NC_:g.)")
    print("  normalized_g Mutalyzer-normalised g. description")
    print("  mane_c       MANE Select transcript c. notation")
    print("  mane_p       MANE Select protein p. notation")
    print("  all_c        all equivalent c. notations (semicolon-separated)")
    print("  all_p        all equivalent p. notations (semicolon-separated)")
    print("  errors       Mutalyzer errors / warnings")


if __name__ == "__main__":
    main()