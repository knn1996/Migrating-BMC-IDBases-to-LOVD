"""
generate_rescue_input.py
========================
Builds a Mutalyzer NM input TSV for ENOSELECTORFOUND variants by replacing
the obsolete or mismatched NM_ accession with the gene's current MANE Select
NM_ transcript.

NM lookup priority
------------------
1. MANE_NM_DIR FASTA file index (same source as Track B; GENE_NM_xxx.x.fasta)
2. LRG_with_NM.csv NM column (rule add_nm_to_lrg, RefSeqGene reference standard)

Output format is identical to generate_mutalyzer_input_NM.py --source MANE
so that run_mutalyzer_NM.py --source MANE can process it without modification.

Env vars
--------
IN_DISPOSITION   unresolved_disposition.tsv (output of classify_unresolved.py)
LRG_NM_CSV       results/LRG_with_NM.csv    (rule add_nm_to_lrg)
MANE_NM_DIR      dna_seq.mane_nm directory  (GENE_NM_xxx.fasta files)
OUT_TSV          rescue Mutalyzer input TSV  (fed to run_mutalyzer_NM.py)
LOG_TSV          tab-separated skip log
"""

import csv
import os
from pathlib import Path

IN_DISPOSITION = os.environ["IN_DISPOSITION"]
LRG_NM_CSV     = os.environ["LRG_NM_CSV"]
MANE_NM_DIR    = os.environ["MANE_NM_DIR"]
OUT_TSV        = os.environ["OUT_TSV"]
LOG_TSV        = os.environ["LOG_TSV"]

for _p in (OUT_TSV, LOG_TSV):
    Path(_p).parent.mkdir(parents=True, exist_ok=True)

CHR_TO_NC = {
    "chr1":  "NC_000001.11", "chr2":  "NC_000002.12", "chr3":  "NC_000003.12",
    "chr4":  "NC_000004.12", "chr5":  "NC_000005.10", "chr6":  "NC_000006.12",
    "chr7":  "NC_000007.14", "chr8":  "NC_000008.11", "chr9":  "NC_000009.12",
    "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9",  "chr15": "NC_000015.10",
    "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
    "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
    "chr22": "NC_000022.11", "chrX":  "NC_000023.11", "chrY":  "NC_000024.10",
}

OUT_COLS = [
    "gene", "accession", "allele_num", "sysname", "hgvs_input",
    "mut_type", "ref", "alt", "chrom", "pos_hg38", "strand",
]


def build_mane_nm_index(mane_nm_dir):
    index = {}
    for p in Path(mane_nm_dir).iterdir():
        if p.suffix.lower() not in (".fasta", ".fa"):
            continue
        parts = p.stem.split("_", 1)
        if len(parts) == 2:
            index[parts[0].upper()] = parts[1]
    return index


def build_lrg_nm_index(csv_path):
    index = {}
    with open(csv_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            gene = (row.get("name") or "").strip().upper()
            nm   = (row.get("NM") or "").strip()
            if gene and nm:
                index[gene] = nm
    return index


def extract_c_notation(hgvs_input):
    s = (hgvs_input or "").strip()
    if ":" in s:
        return s.rsplit(":", 1)[1].strip()
    return s


def main():
    mane_index = build_mane_nm_index(MANE_NM_DIR)
    lrg_index  = build_lrg_nm_index(LRG_NM_CSV)

    print(f"MANE NM index : {len(mane_index)} genes from {MANE_NM_DIR}")
    print(f"LRG NM index  : {len(lrg_index)} genes from {LRG_NM_CSV}")

    with open(IN_DISPOSITION, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    eligible = [
        r for r in rows
        if r.get("category") == "ENOSELECTORFOUND"
        and r.get("disposition") == "RESCUABLE_AUTO"
    ]

    print(f"Total unresolved (distinct)  : {len(rows)}")
    print(f"ENOSELECTORFOUND / RESCUABLE : {len(eligible)}")

    out_rows  = []
    skip_rows = []

    for row in eligible:
        gene = row.get("gene", "").strip().upper()

        nm_acc = mane_index.get(gene) or lrg_index.get(gene)
        if not nm_acc:
            skip_rows.append({
                "gene":       gene,
                "hgvs_input": row.get("hgvs_input", ""),
                "reason":     "no_MANE_NM_found_in_index_or_LRG",
            })
            continue

        c_notation = extract_c_notation(row.get("hgvs_input", ""))
        if not c_notation or not c_notation.startswith("c."):
            skip_rows.append({
                "gene":       gene,
                "hgvs_input": row.get("hgvs_input", ""),
                "reason":     f"cannot_extract_c_notation:{c_notation!r}",
            })
            continue

        chrom  = (row.get("chrom") or "").strip()
        nc_acc = CHR_TO_NC.get(chrom, "")

        if nc_acc:
            new_hgvs = f"{nc_acc}({nm_acc}):{c_notation}"
        else:
            new_hgvs = f"{nm_acc}:{c_notation}"

        out_rows.append({
            "gene":       gene,
            "accession":  row.get("accession", ""),
            "allele_num": row.get("allele_num", ""),
            "sysname":    row.get("sysname", ""),
            "hgvs_input": new_hgvs,
            "mut_type":   row.get("mut_type", ""),
            "ref":        row.get("ref", ""),
            "alt":        row.get("alt", ""),
            "chrom":      chrom,
            "pos_hg38":   row.get("pos_hg38", ""),
            "strand":     row.get("strand", ""),
        })

    with open(OUT_TSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    with open(LOG_TSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "hgvs_input", "reason"], delimiter="\t")
        w.writeheader()
        w.writerows(skip_rows)

    print(f"Rescue input written : {len(out_rows)} rows  -> {OUT_TSV}")
    print(f"Skipped              : {len(skip_rows)} rows  -> {LOG_TSV}")
    if out_rows:
        print("Sample new hgvs_input constructions:")
        for r in out_rows[:5]:
            print(f"  {r['gene']:10s}  {r['hgvs_input']}")


if __name__ == "__main__":
    main()
