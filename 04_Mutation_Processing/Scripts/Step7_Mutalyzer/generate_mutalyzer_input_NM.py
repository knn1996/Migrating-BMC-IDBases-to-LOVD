import argparse
import csv
import os
import re
from pathlib import Path

ALL_MUTATIONS_TSV = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step1_Extraction\all_mutations.tsv"
BED_DIR           = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step5_Liftover"
OUT_TSV_TEMPLATE  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_input_NM_{source}.tsv"
LOG_TSV_TEMPLATE  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_input_NM_skipped_{source}.tsv"

SOURCE_DIRS = {
    "MANE":     r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\Mane_Select_NM",
    "IDRefseq": r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\IDRefseq_NM",
}

OFFSET_CSVS = {
    "MANE":     r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step2_RefCheck\lrg_offset_results_NM_MANE.csv",
    "IDRefseq": r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step2_RefCheck\lrg_offset_results_NM.csv",
}

MATCH_THRESHOLD = 0.9

os.makedirs(r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer", exist_ok=True)

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


def load_valid_genes(offset_csv):
    valid = set()
    with open(offset_csv, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                matched = int(row["matched"])
                total   = int(row["total"])
                pct     = float(row["match_pct"])
            except (ValueError, KeyError):
                continue
            if pct >= MATCH_THRESHOLD or (total - matched) <= 2:
                valid.add(row["gene"].strip().upper())
    return valid


def build_nm_index(directory):
    index = {}
    for p in Path(directory).iterdir():
        if p.suffix.lower() not in (".fasta", ".fa"):
            continue
        parts = p.stem.split("_", 1)
        if len(parts) == 2:
            index[parts[0].upper()] = parts[1]
    return index


def build_bed_index(bed_dir):
    index = {}
    for p in Path(bed_dir).glob("*.BED"):
        gene = p.stem.split("_")[0].upper()
        if gene in index:
            continue
        with open(p, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 6:
                    index[gene] = {"chrom": parts[0], "strand": parts[5]}
                    break
    return index


def parse_c_notation(notation):
    notation = notation.strip()

    m = re.match(r"^c\.([-*]?\d+[+-]?\d*)([A-Z])>([A-Z])$", notation, re.IGNORECASE)
    if m:
        return notation, m.group(2).upper(), m.group(3).upper(), "substitution"

    m = re.match(r"^c\.([-*]?\d+[+-]?\d*)(?:_([-*]?\d+[+-]?\d*))?del([A-Z]*)$", notation, re.IGNORECASE)
    if m:
        return notation, m.group(3).upper() if m.group(3) else "", "", "deletion"

    m = re.match(r"^c\.([-*]?\d+[+-]?\d*)(?:_([-*]?\d+[+-]?\d*))?ins([A-Z]+)$", notation, re.IGNORECASE)
    if m:
        return notation, "", m.group(3).upper(), "insertion"

    m = re.match(r"^c\.([-*]?\d+[+-]?\d*)(?:_([-*]?\d+[+-]?\d*))?delins([A-Z]+)$", notation, re.IGNORECASE)
    if m:
        return notation, "", m.group(3).upper(), "delins"

    m = re.match(r"^c\.([-*]?\d+[+-]?\d*)(?:_([-*]?\d+[+-]?\d*))?dup([A-Z]*)$", notation, re.IGNORECASE)
    if m:
        duped = m.group(3).upper() if m.group(3) else ""
        return notation, duped, duped, "duplication"

    return None, None, None, "unparseable"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", choices=["MANE", "IDRefseq"], default="MANE",
                        help="NM FASTA source: MANE (Mane_Select_NM) or IDRefseq (IDRefseq_NM)")
    args = parser.parse_args()

    fasta_dir   = SOURCE_DIRS[args.source]
    out_tsv     = OUT_TSV_TEMPLATE.format(source=args.source)
    log_tsv     = LOG_TSV_TEMPLATE.format(source=args.source)
    valid_genes = load_valid_genes(OFFSET_CSVS[args.source])

    print(f"Source      : {args.source} -> {fasta_dir}")
    print(f"Valid genes : {len(valid_genes)} (match_pct >= {MATCH_THRESHOLD} or total-matched <= 2)")

    nm_index  = build_nm_index(fasta_dir)
    bed_index = build_bed_index(BED_DIR)

    out_rows     = []
    skipped_rows = []

    with open(ALL_MUTATIONS_TSV, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("variant_type", "").strip() != "coding":
                continue

            gene      = row["gene"].strip().upper()
            accession = row["accession"].strip()
            notation  = row["notation"].strip()

            def skip(reason):
                skipped_rows.append({"gene": gene, "accession": accession,
                                     "notation": notation, "reason": reason})

            if gene not in valid_genes:
                skip("below match threshold"); continue
            if gene not in nm_index:
                skip("no NM fasta found"); continue

            c_notation, ref, alt, mut_type = parse_c_notation(notation)
            if c_notation is None:
                skip("unparseable c. notation"); continue

            bed_info   = bed_index.get(gene, {})
            chrom      = bed_info.get("chrom", "")
            strand     = bed_info.get("strand", "") or "+"
            nc_acc     = CHR_TO_NC.get(chrom, "")
            nm_acc     = nm_index[gene]
            hgvs_input = f"{nc_acc}({nm_acc}):{c_notation}" if nc_acc else f"{nm_acc}:{c_notation}"

            out_rows.append({
                "gene":       gene,
                "accession":  accession,
                "allele_num": row["allele_num"].strip(),
                "sysname":    notation,
                "hgvs_input": hgvs_input,
                "mut_type":   mut_type,
                "ref":        ref,
                "alt":        alt,
                "chrom":      chrom,
                "pos_hg38":   "",
                "strand":     strand,
            })

    with open(out_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "gene", "accession", "allele_num", "sysname", "hgvs_input",
            "mut_type", "ref", "alt", "chrom", "pos_hg38", "strand",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    with open(log_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "accession", "notation", "reason"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(skipped_rows)

    no_nc = sum(1 for r in out_rows if not r["chrom"])
    print(f"Written : {len(out_rows)} rows  -> {out_tsv}")
    print(f"Skipped : {len(skipped_rows)} rows  -> {log_tsv}")
    print(f"No NC_/chrom (NM-only fallback): {no_nc} rows")


if __name__ == "__main__":
    main()
