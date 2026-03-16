"""
generate_bed.py
===============
Generates LiftOver-ready BED files from pre-computed BLAST first-hits.

Inputs
------
  blast_first_hits_source_<N>.tsv   Pre-computed first 100%-identity BLAST
                                    hits (one row per gene per genome build).
                                    Located in Output/blast_first_hits/.

  all_mutations.tsv                 Mutation positions in IDRefSeq coordinates.
                                    Located in Output/all_mutations/.

Logic
-----
  For each gene in all_mutations.tsv:
    1. Find the row in the hits TSV where pct_identity == 100.0.
    2. If the gene has hits in multiple builds, prefer the latest
       (hg18 > hg17 > hg16).
    3. Convert every mutation's (pos_start, pos_end) from IDRefSeq
       coordinates to chromosomal coordinates using the BLAST alignment
       endpoints (linear interpolation).
    4. Write one BED file (UCSC 0-based half-open, LiftOver format).

  Validation: the total number of BED lines written for matched genes
  must equal the number of rows in all_mutations.tsv for those genes.

Usage
-----
    python generate_bed.py --source 1   # source 1 BLAST results
    python generate_bed.py --source 3   # LRG / source 3 BLAST results
"""

import argparse
import csv
import os
import re
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_THESIS_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", ".."))

TSV_PATH = os.path.join(
    _SCRIPT_DIR, "..", "..", "Output", "all_mutations", "all_mutations.tsv"
)
HITS_DIR = os.path.join(_SCRIPT_DIR, "..", "..", "Output", "blast_first_hits")

# Output directory for pre-LiftOver BED files, keyed by source number.
BED_DIRS = {
    1: os.path.join(_THESIS_DIR, "03_BED_Files", "BED_source_1"),
    3: os.path.join(_THESIS_DIR, "03_BED_Files", "BED_LRG_source_3"),
}

# Build priority: higher number = preferred.
BUILD_PRIORITY = {"hg18": 3, "hg17": 2, "hg16": 1}


# ── BLAST hit loading ───────────────────────────────────────────────────────────

def _parse_source1_row(row):
    """Extract fields from a source-1 TSV row (columns from BLAST tabular output)."""
    fname = row["file"]
    gene  = re.sub(r"_DNA$", "", fname.split("_vs_")[0], flags=re.I)
    hg_m  = re.search(r"(hg\d+)", fname, re.I)
    return {
        "gene":         gene,
        "build":        hg_m.group(1).lower() if hg_m else "hg_unknown",
        "chrom":        row["sseqid"].strip(),
        "pct_identity": float(row["pident"]),
        "qstart":       int(row["qstart"]),
        "qend":         int(row["qend"]),
        "sstart":       int(row["sstart"]),
        "send":         int(row["send"]),
    }


def _parse_source3_row(row):
    """Extract fields from a source-3 TSV row (pre-labelled columns)."""
    return {
        "gene":         row["gene"].strip(),
        "build":        row["build"].strip().lower(),
        "chrom":        row["subject_id"].strip(),
        "pct_identity": float(row["pct_identity"]),
        "qstart":       int(row["q_start"]),
        "qend":         int(row["q_end"]),
        "sstart":       int(row["s_start"]),
        "send":         int(row["s_end"]),
    }


def load_blast_hits(source_num):
    """
    Load blast_first_hits_source_<N>.tsv.

    Returns a dict mapping GENE (upper-case) -> best-build entry dict.
    Only rows with pct_identity == 100.0 are kept.
    When a gene appears in multiple builds, the latest build is chosen
    (hg18 > hg17 > hg16).
    """
    fpath = os.path.join(HITS_DIR, f"blast_first_hits_source_{source_num}.tsv")
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"Hits file not found: {fpath}")

    parse_row = _parse_source1_row if source_num == 1 else _parse_source3_row

    hits        = {}   # gene_upper -> entry
    rows_read   = 0
    rows_kept   = 0
    rows_sub100 = 0

    with open(fpath, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows_read += 1
            entry = parse_row(row)

            if entry["pct_identity"] < 100.0:
                rows_sub100 += 1
                continue            # strict: 100% identity only

            gene_key = entry["gene"].upper()
            current  = hits.get(gene_key)
            priority = BUILD_PRIORITY.get(entry["build"], 0)

            if current is None or priority > BUILD_PRIORITY.get(current["build"], 0):
                hits[gene_key] = entry
                rows_kept += 1

    print(f"[HITS] {os.path.basename(fpath)}")
    print(f"       {rows_read} rows read  |  "
          f"{rows_sub100} below 100% skipped  |  "
          f"{len(hits)} genes with a 100%-identity hit")
    return hits


# ── Mutation loading ────────────────────────────────────────────────────────────

def load_mutations(tsv_path):
    """
    Read all_mutations.tsv and group rows by gene (upper-case key).

    Returns:
        mutations  dict[gene -> list of mutation dicts]
        total      total number of rows across all genes
    """
    mutations = defaultdict(list)
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["gene"].strip().upper()
            mutations[gene].append({
                "accession": row["accession"].strip(),
                "sysname":   row["sysname_per_allele"].strip(),
                "pos_start": int(row["pos_start"]),
                "pos_end":   int(row["pos_end"]),
            })

    total = sum(len(v) for v in mutations.values())
    print(f"[TSV]  {os.path.basename(tsv_path)}")
    print(f"       {total} mutations across {len(mutations)} genes")
    return mutations, total


# ── Coordinate conversion ───────────────────────────────────────────────────────

def idref_to_genomic(pos, hit):
    """
    Convert a 1-based IDRefSeq position to a chromosomal coordinate.

    Uses linear interpolation from the BLAST alignment endpoints.
    Handles both forward (sstart <= send) and reverse strands.
    """
    offset = pos - hit["qstart"]
    if hit["sstart"] <= hit["send"]:
        return hit["sstart"] + offset   # forward strand
    return hit["sstart"] - offset       # reverse strand


# ── BED writing ─────────────────────────────────────────────────────────────────

def _make_name(accession, sysname):
    """Return a BED-safe name: accession_sysname with special chars replaced."""
    safe = re.sub(r"[,\s;]+", "_", sysname)
    safe = re.sub(r"[^A-Za-z0-9_.>\-]", "", safe)
    return f"{accession}_{safe}" if safe else accession


def write_bed(gene, hit, mutations, bed_dir):
    """
    Write one BED file for gene using the given BLAST hit entry.

    BED format: 0-based half-open coordinates (UCSC convention).
    SNVs (pos_start == pos_end) are written as exactly 1 bp wide.

    Returns the number of data lines written (== len(mutations)).
    """
    strand   = "+" if hit["sstart"] <= hit["send"] else "-"
    build    = hit["build"]
    out_path = os.path.join(bed_dir, f"{gene}_{build}.BED")

    os.makedirs(bed_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(f"# BED file for {gene} ({build}) — LiftOver format\n")
        fh.write(
            f"# BLAST: {hit['chrom']}  "
            f"qstart={hit['qstart']} qend={hit['qend']}  "
            f"sstart={hit['sstart']} send={hit['send']}\n"
        )
        for mut in mutations:
            g1 = idref_to_genomic(mut["pos_start"], hit)
            g2 = idref_to_genomic(mut["pos_end"],   hit)
            cs = min(g1, g2) - 1    # 0-based start
            ce = max(g1, g2)        # half-open end
            if mut["pos_start"] == mut["pos_end"]:
                ce = cs + 1         # SNV: exactly 1 bp
            fh.write(
                f"{hit['chrom']}\t{cs}\t{ce}\t"
                f"{_make_name(mut['accession'], mut['sysname'])}\t0\t{strand}\n"
            )

    return len(mutations)


# ── Main ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate LiftOver-ready BED files from pre-computed BLAST hits."
    )
    parser.add_argument(
        "--source", type=int, choices=[1, 3], required=True,
        help="BLAST source number (1 = source1, 3 = LRG/source3)"
    )
    args = parser.parse_args()

    bed_dir = BED_DIRS[args.source]
    print(f"Source   : {args.source}")
    print(f"BED dir  : {bed_dir}")
    print()

    hits                 = load_blast_hits(args.source)
    mutations, total_tsv = load_mutations(TSV_PATH)
    print()

    written_files = 0
    written_muts  = 0
    skipped       = []      # (gene, n_muts) — in TSV but no 100% hit

    for gene in sorted(mutations):
        hit = hits.get(gene)
        if hit is None:
            skipped.append((gene, len(mutations[gene])))
            continue

        n = write_bed(gene, hit, mutations[gene], bed_dir)
        print(
            f"  WROTE  {gene}_{hit['build']}.BED  "
            f"({n} mutations  {hit['chrom']}  "
            f"{'fwd' if hit['sstart'] <= hit['send'] else 'rev'})"
        )
        written_files += 1
        written_muts  += n

    skipped_muts = sum(n for _, n in skipped)
    expected     = written_muts + skipped_muts
    match_symbol = "✓" if expected == total_tsv else "✗"

    print()
    print("=" * 62)
    print(f"  BED files written        : {written_files}")
    print(f"  Mutations in TSV         : {total_tsv}")
    print(f"  BED lines written        : {written_muts}")
    print(f"  Skipped (no 100% hit)    : {skipped_muts} "
          f"mutations across {len(skipped)} genes")
    print(f"  Count check              : {match_symbol} "
          f"{expected} accounted / {total_tsv} in TSV")
    print("=" * 62)

    if skipped:
        print("\nSkipped genes (no 100%-identity BLAST hit):")
        for gene, n in skipped:
            print(f"  {gene:20s}  {n} mutations")


if __name__ == "__main__":
    main()
