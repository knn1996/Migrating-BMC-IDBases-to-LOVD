"""
generate_bed.py
===============
Generates LiftOver-ready BED files from all_mutations.tsv.

For each gene in all_mutations.tsv:
  - Reads DNA mutation positions (pos_start, pos_end) directly from the TSV
  - Maps IDRefSeq positions to hg16/17/18 chromosomal coordinates via
    100%-identity BLAST hits across all *_vs_hg*.txt files in BASE_DIR
  - Writes GENE_hgXX.BED files into BED_DIR, ready for UCSC LiftOver

Usage:
    python generate_bed.py
"""

import csv
import glob
import os
import re
from collections import defaultdict
from typing import Optional

# -- Configuration --------------------------------------------------------------

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR    = os.path.join(_SCRIPT_DIR, "..", "BLAST_Results", "source1")
TSV_PATH    = os.path.join(_SCRIPT_DIR, "..", "Output", "all_mutations.tsv")
BED_DIR     = os.path.join(_SCRIPT_DIR, "..", "BED", "hg38")

# Gene-name aliases: maps gene symbol -> list of alternative keys in the BLAST files.
ALIASES = {}

# -- Helpers --------------------------------------------------------------------

def load_blast_map(base_dir):
    """
    Glob all *_vs_hg*.txt files in base_dir, load each one, and merge
    results into a single blast_map.

    Per-build deduplication
    -----------------------
    Pass 1 - exact-duplicate removal:
        Rows sharing the same (hg, chrom, sstart, send) are collapsed to one.
        These arise when a BLAST file contains the same alignment twice
        (e.g. forward/reverse tie-breaking artefacts).

    Pass 2 - best-hit-per-build selection:
        Of all remaining 100%-identity hits for a given (gene, hg) pair,
        only the one with the longest query alignment (qend - qstart) is kept.
        This eliminates spurious hits against pseudogenes, paralogs and
        repeat regions. The longest match is the true full-length locus.

    Returns dict keyed by GENE_KEY_UPPER; each gene has at most one entry per build.
    """
    pattern   = os.path.join(base_dir, "*_vs_hg*.txt")
    txt_files = sorted(glob.glob(pattern))

    if not txt_files:
        raise FileNotFoundError(f"No *_vs_hg*.txt files found in {base_dir}")

    blast_map  = defaultdict(list)
    total_skip = 0

    for path in txt_files:
        fname    = os.path.basename(path)
        hg_match = re.search(r"(hg\d+)", fname, re.IGNORECASE)
        hg       = hg_match.group(1).lower() if hg_match else "hg_unknown"
        gene_key = fname.split("_vs_")[0].upper()

        skipped = 0
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 10:
                    continue
                try:
                    if float(cols[2]) != 100.0:
                        skipped += 1
                        continue
                except ValueError:
                    continue

                blast_map[gene_key].append({
                    "hg":     hg,
                    "chrom":  cols[1].strip(),
                    "qstart": int(cols[6]),
                    "qend":   int(cols[7]),
                    "sstart": int(cols[8]),
                    "send":   int(cols[9]),
                })

        total_skip += skipped
        print(f"  [BLAST] {fname}  ->  gene={gene_key}  build={hg}"
              f"  ({skipped} partial-identity rows skipped)")

    # Pass 1: remove exact duplicates (same genomic span)
    for key in blast_map:
        seen    = set()
        deduped = []
        for e in blast_map[key]:
            sig = (e["hg"], e["chrom"], e["sstart"], e["send"])
            if sig not in seen:
                seen.add(sig)
                deduped.append(e)
        n_before = len(blast_map[key])
        blast_map[key] = deduped
        if len(deduped) < n_before:
            print(f"  [DEDUP] {key}: {n_before} -> {len(deduped)} entries"
                  f" after exact-duplicate removal")

    # Pass 2: keep only the best hit per (gene, build)
    # Best = longest query alignment, which identifies the true full-length locus.
    for key in blast_map:
        best_by_hg = {}
        for e in blast_map[key]:
            hg      = e["hg"]
            aln_len = e["qend"] - e["qstart"]
            if hg not in best_by_hg or aln_len > (best_by_hg[hg]["qend"] - best_by_hg[hg]["qstart"]):
                best_by_hg[hg] = e
        n_before = len(blast_map[key])
        blast_map[key] = list(best_by_hg.values())
        if len(blast_map[key]) < n_before:
            kept = {e["hg"]: f"{e['chrom']}:{e['sstart']}-{e['send']}" for e in blast_map[key]}
            print(f"  [BEST]  {key}: {n_before} -> {len(blast_map[key])} entries"
                  f" after best-hit-per-build selection  |  kept: {kept}")

    print(f"[BLAST] {len(txt_files)} file(s) loaded  |  "
          f"{len(blast_map)} unique genes  |  "
          f"{total_skip} partial-identity rows skipped total")
    return blast_map


def load_mutations(tsv_path):
    """
    Read all_mutations.tsv and group mutations by gene.
    Expected columns: gene, accession, sysname, pos_start, pos_end, feat_name, allele_num
    """
    mutations = defaultdict(list)
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = row["gene"].strip().upper()
            mutations[gene].append({
                "accession": row["accession"].strip(),
                "sysname":   row["sysname"].strip(),
                "pos_start": int(row["pos_start"]),
                "pos_end":   int(row["pos_end"]),
            })
    print(f"[TSV] {tsv_path} loaded  |  "
          f"{sum(len(v) for v in mutations.values())} mutations across "
          f"{len(mutations)} genes")
    return mutations


def blast_key(gene, blast_map):
    """Return the blast_map key for gene, or None."""
    candidates = [f"{gene}_DNA", gene] + ALIASES.get(gene, [])
    for candidate in candidates:
        if candidate.upper() in blast_map:
            return candidate.upper()
    return None


def idref_to_genomic(pos, entry):
    """
    Convert a 1-based IDRefSeq position to a chromosomal coordinate via
    linear extrapolation from the BLAST alignment endpoints.
    """
    offset = pos - entry["qstart"]
    if entry["sstart"] <= entry["send"]:
        return entry["sstart"] + offset   # forward strand
    return entry["sstart"] - offset       # reverse strand


def make_bed_name(accession, sysname):
    """Build a BED-safe name field from accession + sysname."""
    safe = re.sub(r"[,\s;]+", "_", sysname)
    safe = re.sub(r"[^A-Za-z0-9_.>\-]", "", safe)
    return f"{accession}_{safe}" if safe else accession


def write_bed(gene, hg, entry, mutations, source_path):
    """
    Write a BED file for one gene / genome-build combination.
    BED coordinates are 0-based half-open (UCSC convention).
    SNVs (pos_start == pos_end) are written as exactly 1 bp wide.
    Returns number of data lines written.
    """
    strand = "+" if entry["sstart"] <= entry["send"] else "-"

    header_lines = [
        f"# BED file for {gene} ({hg}) -- LiftOver format",
        f"# Source: {source_path}",
        f"# Chrom: {entry['chrom']}  |  BLAST: "
        f"qstart={entry['qstart']} qend={entry['qend']} "
        f"sstart={entry['sstart']} send={entry['send']}",
    ]

    bed_lines = []
    for mut in mutations:
        g1 = idref_to_genomic(mut["pos_start"], entry)
        g2 = idref_to_genomic(mut["pos_end"],   entry)
        cs = min(g1, g2) - 1
        ce = max(g1, g2)
        if mut["pos_start"] == mut["pos_end"]:
            ce = cs + 1

        name = make_bed_name(mut["accession"], mut["sysname"])
        bed_lines.append(f"{entry['chrom']}\t{cs}\t{ce}\t{name}\t0\t{strand}")

    out_path = os.path.join(BED_DIR, f"{gene}_{hg}.BED")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(header_lines + bed_lines) + "\n")

    return len(bed_lines)


# -- Main -----------------------------------------------------------------------

def main():
    os.makedirs(BED_DIR, exist_ok=True)

    blast_map = load_blast_map(BASE_DIR)
    mutations = load_mutations(TSV_PATH)

    total_files = total_muts = 0
    skipped = []

    for gene, gene_mutations in sorted(mutations.items()):
        key = blast_key(gene, blast_map)
        if not key:
            skipped.append((gene, "no 100% BLAST match"))
            continue

        for entry in blast_map[key]:
            n = write_bed(gene, entry["hg"], entry, gene_mutations, TSV_PATH)
            print(f"  WROTE  {gene}_{entry['hg']}.BED  ({n} mutations)")
            total_files += 1
            total_muts  += n

    print(f"\nBED files written : {total_files}")
    print(f"Total entries     : {total_muts}")
    print(f"Skipped genes     : {len(skipped)}")
    for gene, reason in skipped:
        print(f"  SKIP  {gene}: {reason}")


if __name__ == "__main__":
    main()