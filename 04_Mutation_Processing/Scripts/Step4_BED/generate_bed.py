"""
generate_bed.py
===============
Generates LiftOver-ready BED files from all_mutations.tsv + BLAST first-hit
coordinates.

Usage:
    python generate_bed.py --source IDRefseq

Selection rule per gene:
    Best = highest pct_identity; ties broken by highest build (hg18 > hg17 > hg16).
"""

import os
import re
import argparse
import csv

THESIS_DIR    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
PIPELINE_DIR  = os.path.join(THESIS_DIR, "04_Mutation_Processing")
MUTATIONS_TSV = os.path.join(PIPELINE_DIR, "Output", "Step1_Extraction", "all_mutations.tsv")
BED_DIR       = os.path.join(THESIS_DIR, "03_BED_Files", "BED_IDRefseq")
LOG_DIR       = os.path.join(PIPELINE_DIR, "Logs")

SOURCE_BLAST = {
    "IDRefseq": os.path.join(PIPELINE_DIR, "Output", "Step3_BLAST", "blast_first_hits_source_IDRefseq.tsv"),
}

BUILD_RANK = {"hg16": 0, "hg17": 1, "hg18": 2}


def load_blast_map(blast_path):
    """
    Returns dict: gene_symbol (query_id) -> best hit dict.
    Best = highest pct_identity; ties broken by highest build rank.
    """
    best = {}

    with open(blast_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene  = row["query_id"].strip()
            build = row["build"].strip().lower()
            pct   = float(row["pct_identity"])
            rank  = BUILD_RANK.get(build, -1)

            entry = {
                "build":        build,
                "chrom":        row["subject_id"].strip(),
                "q_start":      int(row["q_start"]),
                "q_end":        int(row["q_end"]),
                "s_start":      int(row["s_start"]),
                "s_end":        int(row["s_end"]),
                "pct_identity": pct,
                "rank":         rank,
            }

            if gene not in best:
                best[gene] = entry
            else:
                cur = best[gene]
                if (pct, rank) > (cur["pct_identity"], cur["rank"]):
                    best[gene] = entry

    return best


def idref_to_genomic(pos, entry):
    offset = pos - entry["q_start"]
    if entry["s_start"] <= entry["s_end"]:
        return entry["s_start"] + offset
    return entry["s_start"] - offset


def make_bed_name(accession, notation):
    safe = re.sub(r"[,\s;]+", "_", notation)
    safe = re.sub(r"[^A-Za-z0-9_.>\-]", "", safe)
    return f"{accession}_{safe}" if safe else accession


def load_mutations(tsv_path):
    """
    Returns dict: gene -> list of mutation dicts.
    Skips rows with missing pos_start/pos_end or non-empty missing_flag.
    """
    by_gene = {}
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("missing_flag", "").strip():
                continue
            try:
                ps = int(row["pos_start"])
                pe = int(row["pos_end"])
            except (ValueError, KeyError):
                continue

            gene = row["gene"].strip()
            by_gene.setdefault(gene, []).append({
                "accession": row["accession"].strip(),
                "notation":  row["notation"].strip(),
                "pos_start": ps,
                "pos_end":   pe,
            })
    return by_gene


def write_bed(gene, entry, mutations, out_dir):
    strand = "+" if entry["s_start"] <= entry["s_end"] else "-"
    lines = [
        f"# BED for {gene} — build={entry['build']} chrom={entry['chrom']} "
        f"pct_id={entry['pct_identity']} "
        f"q={entry['q_start']}-{entry['q_end']} "
        f"s={entry['s_start']}-{entry['s_end']}"
    ]

    seen = set()
    for mut in mutations:
        g1 = idref_to_genomic(mut["pos_start"], entry)
        g2 = idref_to_genomic(mut["pos_end"],   entry)
        cs = min(g1, g2) - 1
        ce = max(g1, g2)
        if mut["pos_start"] == mut["pos_end"]:
            ce = cs + 1

        name = make_bed_name(mut["accession"], mut["notation"])
        key  = (entry["chrom"], cs, ce, name)
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"{entry['chrom']}\t{cs}\t{ce}\t{name}\t0\t{strand}")

    out_path = os.path.join(out_dir, f"{gene}_{entry['build']}.BED")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    return len(lines) - 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, choices=list(SOURCE_BLAST.keys()))
    args = parser.parse_args()

    blast_path = SOURCE_BLAST[args.source]
    os.makedirs(BED_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    log_path  = os.path.join(LOG_DIR, f"generate_bed_{args.source}.log")
    log_lines = []

    def log(msg):
        print(msg)
        log_lines.append(msg)

    log(f"Source     : {args.source}")
    log(f"BLAST file : {blast_path}")
    log(f"Mutations  : {MUTATIONS_TSV}")
    log(f"BED output : {BED_DIR}")
    log("")

    blast_map = load_blast_map(blast_path)
    by_gene   = load_mutations(MUTATIONS_TSV)

    log(f"BLAST entries loaded : {len(blast_map)} genes")
    log(f"Mutation genes       : {len(by_gene)}")
    log("")

    n_written = n_skipped = 0
    skipped = []

    for gene in sorted(by_gene.keys()):
        entry = blast_map.get(gene)
        if entry is None:
            skipped.append((gene, "no BLAST hit"))
            n_skipped += 1
            continue

        n = write_bed(gene, entry, by_gene[gene], BED_DIR)
        log(f"  WROTE  {gene}_{entry['build']}.BED  ({n} entries, pct_id={entry['pct_identity']}, build={entry['build']})")
        n_written += 1

    log("")
    log(f"BED files written : {n_written}")
    log(f"Genes skipped     : {n_skipped}")
    for gene, reason in skipped:
        log(f"  SKIP  {gene}: {reason}")

    with open(log_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()


