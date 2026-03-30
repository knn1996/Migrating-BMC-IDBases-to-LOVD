"""
filter_bed.py
=============
Filters non-matching accessions out of IDRefseq BED files based on
reference_summary.csv, then writes clean BED files to Output directories.

Usage:
    python filter_bed.py

Input:
    03_BED_Files/BED_IDRefseq/*.BED          -> Output/Step4_BED/
    03_BED_Files/BED_IDRefseq/hg38/*.BED     -> Output/Step5_Liftover/

Logic:
    - For each gene, pick the IDRefseq source (higher match_pct; source_3 if tied)
    - Accessions listed in non_matching for that source are excluded from the BED
    - Genes with 100% match are copied as-is (no row filtering needed)
    - Excluded rows are logged to a TSV for manual curation
"""

import os
import re
import csv
import glob
import shutil

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_THESIS_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", ".."))
_OUT_BASE   = os.path.join(_SCRIPT_DIR, "..", "..", "Output")

REF_SUMMARY  = os.path.join(_OUT_BASE, "Step2_RefCheck", "reference_summary.csv")
BED_ORIGINAL = os.path.join(_THESIS_DIR, "03_BED_Files", "BED_IDRefseq")
BED_HG38     = os.path.join(_THESIS_DIR, "03_BED_Files", "BED_IDRefseq", "hg38")
OUT_STEP4    = os.path.join(_OUT_BASE, "Step4_BED")
OUT_STEP5    = os.path.join(_OUT_BASE, "Step5_Liftover")
OUT_EXCLUDED = os.path.join(_OUT_BASE, "Step4_BED", "excluded_accessions.tsv")


def load_excluded_accessions():
    """Return {gene: set(accession_ids)} for accessions that failed reference matching."""
    best = {}
    with open(REF_SUMMARY, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            gene = row["gene"]
            pct  = float(row["match_pct"])
            src  = row["source"]
            if gene not in best:
                best[gene] = row
            else:
                cur_pct = float(best[gene]["match_pct"])
                cur_src = best[gene]["source"]
                if pct > cur_pct or (pct == cur_pct and src == "source_3"):
                    best[gene] = row

    excluded = {}
    for gene, row in best.items():
        nm = row.get("non_matching", "").strip()
        if nm:
            excluded[gene] = set(nm.split(";"))
    return excluded


def gene_from_filename(fname):
    m = re.match(r"^(.+?)_hg\d+", fname, re.I)
    return m.group(1) if m else None


def filter_bed(src_path, dst_path, excluded_accs, gene, excluded_log):
    """
    Filter one BED file. Returns (kept, removed) counts.
    If gene has no excluded accessions, copies the file directly.
    """
    if gene not in excluded_accs:
        shutil.copy2(src_path, dst_path)
        total = sum(1 for l in open(src_path, encoding="utf-8")
                    if l.strip() and not l.startswith("#"))
        return total, 0

    to_exclude = excluded_accs[gene]
    kept = removed = 0
    with open(src_path, encoding="utf-8") as fin, \
         open(dst_path, "w", encoding="utf-8") as fout:
        for line in fin:
            stripped = line.rstrip()
            if not stripped or stripped.startswith("#"):
                fout.write(line)
                continue
            parts = stripped.split("\t")
            if len(parts) < 4:
                fout.write(line)
                continue
            acc = parts[3].split("_")[0]
            if acc in to_exclude:
                excluded_log.append({
                    "gene": gene,
                    "file": os.path.basename(src_path),
                    "chrom": parts[0],
                    "pos": parts[1],
                    "name": parts[3],
                    "accession": acc,
                })
                removed += 1
            else:
                fout.write(line)
                kept += 1
    return kept, removed


def process_dir(src_dir, out_dir, excluded_accs, label, excluded_log):
    files = sorted(glob.glob(os.path.join(src_dir, "*.BED")))
    print(f"\n{label}: {len(files)} files -> {out_dir}")
    os.makedirs(out_dir, exist_ok=True)

    total_kept = total_removed = total_copied = 0
    for fpath in files:
        fname = os.path.basename(fpath)
        gene  = gene_from_filename(fname)
        dst   = os.path.join(out_dir, fname)
        kept, removed = filter_bed(fpath, dst, excluded_accs, gene, excluded_log)
        if removed == 0:
            total_copied += 1
        else:
            total_removed += removed
            print(f"  {gene}: removed {removed} rows ({kept} kept)")
        total_kept += kept

    print(f"  Copied as-is : {total_copied} files")
    print(f"  Rows removed : {total_removed}")
    print(f"  Rows kept    : {total_kept}")


def write_excluded_log(excluded_log):
    with open(OUT_EXCLUDED, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "file", "chrom", "pos", "name", "accession"], delimiter="\t")
        writer.writeheader()
        writer.writerows(excluded_log)
    print(f"\nExclusion log: {OUT_EXCLUDED} ({len(excluded_log)} rows)")


def main():
    print("filter_bed.py  |  IDRefseq BED filtering\n")
    excluded_accs = load_excluded_accessions()
    print(f"Genes with excluded accessions: {len(excluded_accs)}")
    for gene, accs in sorted(excluded_accs.items()):
        print(f"  {gene}: {', '.join(sorted(accs))}")

    excluded_log = []
    process_dir(BED_ORIGINAL, OUT_STEP4, excluded_accs, "Original BED (hg18)", excluded_log)
    process_dir(BED_HG38,     OUT_STEP5, excluded_accs, "hg38 BED",            excluded_log)
    write_excluded_log(excluded_log)


if __name__ == "__main__":
    main()
