import argparse
import csv
import os
import re
from pathlib import Path

ALL_MUTATIONS_TSV = os.environ["ALL_MUTATIONS_TSV"]
BED_DIR           = os.environ["BED_DIR"]
OUT_TSV  = os.environ["OUT_TSV"]
LOG_TSV  = os.environ["LOG_TSV"]
SOURCE_DIRS = {
    "MANE":     os.environ.get("MANE_NM_DIR", ""),
    "IDRefseq": os.environ.get("IDREFSEQ_NM_DIR", ""),
}

SOURCE_NG_DIRS = {
    "MANE":     os.environ.get("MANE_NG_DIR", ""),
    "IDRefseq": None,
}

EXCLUDE_GENES = {"BTK", "SH2"}

CODING_OFFSET = {"ADA": -95}

NM_ONLY_GENES = {"ORAI1"}

os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
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


def build_acc_index(directory):
    index = {}
    if not directory:
        return index
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


def apply_coding_offset(c_notation, offset):
    if not c_notation.startswith("c."):
        return c_notation

    def repl(m):
        sign = m.group("sign")
        if sign in ("-", "*"):
            return m.group(0)
        new = int(m.group("num")) + offset
        if new <= 0:
            new -= 1
        return f"{new}{m.group('intron') or ''}"

    pattern = re.compile(r"(?P<sign>[-*]?)(?P<num>\d+)(?P<intron>[+-]\d+)?")
    return "c." + pattern.sub(repl, c_notation[2:])


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

    fasta_dir = SOURCE_DIRS[args.source]
    ng_dir    = SOURCE_NG_DIRS.get(args.source)
    out_tsv   = OUT_TSV
    log_tsv   = LOG_TSV

    nm_index  = build_acc_index(fasta_dir)
    ng_index  = build_acc_index(ng_dir)
    bed_index = build_bed_index(BED_DIR)

    apply_offsets = CODING_OFFSET if args.source == "MANE" else {}

    print(f"Source      : {args.source} -> {fasta_dir}")
    print(f"Scope       : {len(nm_index)} genes with an NM FASTA, minus excluded {sorted(EXCLUDE_GENES)}")
    print(f"NG reference: {len(ng_index)} genes ({ng_dir})")
    print(f"Coding offset corrections applied: {apply_offsets}")

    out_rows     = []
    skipped_rows = []
    n_offset_applied = 0

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

            if gene in EXCLUDE_GENES:
                skip("out of scope"); continue
            if gene not in nm_index:
                skip("no NM fasta found"); continue

            c_notation, ref, alt, mut_type = parse_c_notation(notation)
            if c_notation is None:
                skip("unparseable c. notation"); continue

            if gene in apply_offsets:
                c_notation = apply_coding_offset(c_notation, apply_offsets[gene])
                n_offset_applied += 1

            bed_info = bed_index.get(gene, {})
            chrom    = bed_info.get("chrom", "")
            strand   = bed_info.get("strand", "") or "+"
            nc_acc   = CHR_TO_NC.get(chrom, "")
            nm_acc   = nm_index[gene]
            ng_acc   = ng_index.get(gene, "")

            if gene in NM_ONLY_GENES:
                hgvs_input = f"{nm_acc}:{c_notation}"
            elif nc_acc:
                hgvs_input = f"{nc_acc}({nm_acc}):{c_notation}"
            elif ng_acc:
                hgvs_input = f"{ng_acc}({nm_acc}):{c_notation}"
            else:
                hgvs_input = f"{nm_acc}:{c_notation}"

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

    via_nc = sum(1 for r in out_rows if r["chrom"])
    via_ng = sum(1 for r in out_rows if not r["chrom"] and ng_index.get(r["gene"]) and r["gene"] not in NM_ONLY_GENES)
    via_nm = len(out_rows) - via_nc - via_ng
    print(f"Written : {len(out_rows)} rows  -> {out_tsv}")
    print(f"Skipped : {len(skipped_rows)} rows  -> {log_tsv}")
    print(f"Offset-corrected rows: {n_offset_applied}")
    print(f"  reference used: NC_(NM_)={via_nc}  NG_(NM_)={via_ng}  NM_only={via_nm}")


if __name__ == "__main__":
    main()
