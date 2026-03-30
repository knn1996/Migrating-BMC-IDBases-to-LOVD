import os
import re
import csv
from pathlib import Path

BED_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step5_Liftover"
OUT_TSV = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_input.tsv"
LOG_TSV = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer\mutalyzer_skipped.tsv"

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

COMP = {"A":"T","T":"A","G":"C","C":"G","a":"t","t":"a","g":"c","c":"g"}

RE_SUBST   = re.compile(r"^g\.(?:IVS\d+[+-]?\d*)?(\d+)([A-Z])>([A-Z])$", re.IGNORECASE)
RE_DEL     = re.compile(r"^g\.(?:IVS\d+[+-]?\d*)?(\d+)(?:[._-]+(\d+))?del([A-Z]*)$", re.IGNORECASE)
RE_INS     = re.compile(r"^g\.(?:IVS\d+[+-]?\d*)?(\d+)(?:[._-]+(\d+))?ins([A-Z]+)$", re.IGNORECASE)
RE_DELINS  = re.compile(r"^g\.(?:IVS\d+[+-]?\d*)?(\d+)(?:[._-]+(\d+))?delins([A-Z]+)$", re.IGNORECASE)
RE_DUP     = re.compile(r"^g\.(?:IVS\d+[+-]?\d*)?(\d+)(?:[._-]+(\d+))?dup([A-Z]*)$", re.IGNORECASE)

RE_IVS_SUBST  = re.compile(r"^g\.IVS\d+[+-]\d+([A-Z])>([A-Z])$", re.IGNORECASE)
RE_IVS_DEL    = re.compile(r"^g\.IVS\d+[+-]\d+([A-Z]+)>$", re.IGNORECASE)
RE_IVS_BARE   = re.compile(r"^g\.IVS\d+[+-]?\d*([A-Z]*)>([A-Z]*)$", re.IGNORECASE)


def complement(seq):
    return "".join(COMP.get(b, b) for b in seq)


def strip_ivs(sysname):
    m = re.match(r"^(g\.)IVS\d+[+-]?\d*(.*)$", sysname, re.IGNORECASE)
    if m:
        return m.group(1) + m.group(2)
    return sysname


def parse_sysname(sysname, chrom_pos_0based, strand):
    pos    = chrom_pos_0based + 1
    is_rev = (strand == "-")
    is_ivs = bool(re.match(r"^g\.IVS", sysname, re.IGNORECASE))

    normalized = strip_ivs(sysname) if is_ivs else sysname

    m = RE_SUBST.match(normalized)
    if m:
        ref = m.group(2).upper()
        alt = m.group(3).upper()
        if is_rev:
            ref, alt = complement(ref), complement(alt)
        return f"g.{pos}{ref}>{alt}", ref, alt, "substitution", is_ivs

    m = RE_DEL.match(normalized)
    if m:
        deleted = m.group(3).upper() if m.group(3) else "N"
        end_raw = m.group(2)
        span = (int(end_raw) - int(m.group(1)) + 1) if end_raw else (len(deleted) if deleted != "N" else 1)
        if deleted == "N" or len(deleted) != span:
            deleted = "N" * span
        elif is_rev:
            deleted = complement(deleted)
        end_pos  = pos + span - 1
        notation = f"g.{pos}_{end_pos}del" if span > 1 else f"g.{pos}del"
        return notation, deleted, "", "deletion", is_ivs

    m = RE_INS.match(normalized)
    if m:
        inserted = m.group(3).upper()
        if is_rev:
            inserted = complement(inserted)
        return f"g.{pos}_{pos+1}ins{inserted}", "", inserted, "insertion", is_ivs

    m = RE_DELINS.match(normalized)
    if m:
        inserted = m.group(3).upper()
        if is_rev:
            inserted = complement(inserted)
        end_raw = m.group(2)
        end_pos = int(end_raw) if end_raw else pos
        return f"g.{pos}_{end_pos}delins{inserted}", "N", inserted, "delins", is_ivs

    m = RE_DUP.match(normalized)
    if m:
        duped = m.group(3).upper() if m.group(3) else "N"
        if is_rev and duped != "N":
            duped = complement(duped)
        return f"g.{pos}dup", duped, duped, "duplication", is_ivs

    if is_ivs:
        m = RE_IVS_SUBST.match(sysname)
        if m:
            ref, alt = m.group(1).upper(), m.group(2).upper()
            if is_rev:
                ref, alt = complement(ref), complement(alt)
            return f"g.{pos}{ref}>{alt}", ref, alt, "substitution", True

        m = RE_IVS_DEL.match(sysname)
        if m:
            deleted = m.group(1).upper()
            if is_rev:
                deleted = complement(deleted)
            span    = len(deleted)
            end_pos = pos + span - 1
            notation = f"g.{pos}_{end_pos}del" if span > 1 else f"g.{pos}del"
            return notation, deleted, "", "deletion", True

        m = RE_IVS_BARE.match(sysname)
        if m:
            ref, alt = m.group(1).upper(), m.group(2).upper()
            if ref and not alt:
                if is_rev:
                    ref = complement(ref)
                span    = len(ref)
                end_pos = pos + span - 1
                notation = f"g.{pos}_{end_pos}del" if span > 1 else f"g.{pos}del"
                return notation, ref, "", "deletion", True
            if ref and alt:
                if is_rev:
                    ref, alt = complement(ref), complement(alt)
                return f"g.{pos}{ref}>{alt}", ref, alt, "substitution", True

    return None, None, None, "unparseable", is_ivs


def main():
    out_rows     = []
    skipped_rows = []

    for bed_file in sorted(Path(BED_DIR).glob("*.BED")):
        gene = bed_file.stem.split("_")[0]

        with open(bed_file, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue

                chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                strand = parts[5] if len(parts) > 5 else "+"

                if "_" not in name:
                    skipped_rows.append({"gene": gene, "accession": "", "sysname": name, "reason": "no underscore in name"})
                    continue

                accession, sysname = name.split("_", 1)

                if not sysname.lower().startswith("g."):
                    continue

                nc = CHR_TO_NC.get(chrom)
                if not nc:
                    skipped_rows.append({"gene": gene, "accession": accession, "sysname": sysname, "reason": f"unknown chrom: {chrom}"})
                    continue

                g_notation, ref, alt, mut_type, was_ivs = parse_sysname(sysname, start, strand)

                if g_notation is None:
                    skipped_rows.append({"gene": gene, "accession": accession, "sysname": sysname, "reason": "unparseable sysname"})
                    continue

                hgvs = f"{nc}:{g_notation}"
                out_rows.append({
                    "gene":       gene,
                    "accession":  accession,
                    "sysname":    sysname,
                    "hgvs_input": hgvs,
                    "mut_type":   mut_type,
                    "ref":        ref,
                    "alt":        alt,
                    "chrom":      chrom,
                    "pos_hg38":   start + 1,
                    "strand":     strand,
                    "was_ivs":    was_ivs,
                })

    with open(OUT_TSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene","accession","sysname","hgvs_input","mut_type","ref","alt","chrom","pos_hg38","strand","was_ivs"], delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    with open(LOG_TSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene","accession","sysname","reason"], delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(skipped_rows)

    print(f"Written  {len(out_rows)} rows to {OUT_TSV}")
    print(f"Skipped  {len(skipped_rows)} rows to {LOG_TSV}")


if __name__ == "__main__":
    main()
