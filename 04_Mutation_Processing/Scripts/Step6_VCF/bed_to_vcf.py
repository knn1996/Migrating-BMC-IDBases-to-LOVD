"""
bed_to_vcf.py
=============
Converts IDRefseq hg38 BED files into a VCF 4.2 file.

Usage:
    python bed_to_vcf.py

Input:
    Output/Step5_Liftover/*.BED   (filtered hg38 BED files from filter_bed.py)

Output (written to Output/Step6_VCF/):
    all_variants_hg38_IDRefseq.vcf      one row per unique (chrom, pos, ref, alt, gene)
    variant_accessions_IDRefseq.tsv     all accessions sharing each variant
    bed_to_vcf_skipped_IDRefseq.tsv     genuinely unparseable g. entries

Notes:
    - Only g. (genomic) rows are processed — Mutalyzer re-derives c., r., p. from g.
    - Verified: every c., r., p., and bare-notation accession has a g. companion in the BED files
    - Non-matching accessions already excluded upstream by filter_bed.py
    - Dedup key is (chrom, pos, ref, alt, gene) — gene prevents cross-gene collapse
    - All accessions sharing a variant are preserved in the accession mapping file
"""

import os
import re
import glob
from collections import Counter

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_OUT_DIR    = os.path.join(_SCRIPT_DIR, "..", "..", "Output", "Step6_VCF")

BED_DIR  = os.path.join(_SCRIPT_DIR, "..", "..", "Output", "Step5_Liftover")
BED_GLOB = "*.BED"

OUT_VCF        = os.path.join(_OUT_DIR, "all_variants_hg38_IDRefseq.vcf")
OUT_ACCESSIONS = os.path.join(_OUT_DIR, "variant_accessions_IDRefseq.tsv")
OUT_SKIPPED    = os.path.join(_OUT_DIR, "bed_to_vcf_skipped_IDRefseq.tsv")

RE_SUBST_G     = re.compile(r"g\.(?:IVS[\d_+\-]+|[-+*\d]+)([ACGT])>([ACGT])", re.I)
RE_IVS_DEL     = re.compile(r"g\.IVS[\d_+\-]+([ACGT]+)>\s*$", re.I)
RE_SUBST_NOALT = re.compile(r"g\.[^\t\n]+[ACGT]>(?:_|\s*$)", re.I)
RE_DELINS      = re.compile(r"g\.\d+(?:_\d+)?delins([ACGT]+)", re.I)
RE_INTRON_DEL  = re.compile(r"g\.\d+[-+]\d+_\d+del([ACGTacgt]+)", re.I)
RE_INTRON_DINS = re.compile(r"g\.\d+[-+]\d+([ACGT]{2,})>([ACGT]*)", re.I)
RE_DEL_SEQ     = re.compile(r"g\.\d+(?:_\d+)?del([ACGT]+)", re.I)
RE_DEL_NOSEQ   = re.compile(r"g\.\d+(?:_\d+)?del(?![ACGT])", re.I)
RE_DOTDOT_DEL  = re.compile(r"g\.(\d+)\.\.(\d+)del", re.I)
RE_DUP_SEQ     = re.compile(r"g\.\d+(?:_\d+)?dup([ACGT]+)", re.I)
RE_DUP_NOSEQ   = re.compile(r"g\.\d+(?:_\d+)?dup(?![ACGT])", re.I)
RE_INS         = re.compile(r"g\.\d+_\d+ins([ACGT]+)", re.I)


def parse_ref_alt(name):
    for pat, fn in [
        (RE_DELINS,      lambda m: ("N", m.group(1).upper())),
        (RE_INTRON_DEL,  lambda m: ("N" + m.group(1).upper(), "N")),
        (RE_INTRON_DINS, lambda m: ("N" + m.group(1).upper(), m.group(2).upper() or "N")),
        (RE_SUBST_G,     lambda m: (m.group(1).upper(), m.group(2).upper())),
        (RE_IVS_DEL,     lambda m: ("N" + m.group(1).upper(), "N")),
    ]:
        m = pat.search(name)
        if m:
            return fn(m)

    if RE_SUBST_NOALT.search(name): return None, None
    if RE_DOTDOT_DEL.search(name):  return "N", "."
    m = RE_DEL_SEQ.search(name)
    if m:                           return "N" + m.group(1).upper(), "N"
    if RE_DEL_NOSEQ.search(name):   return "N", "."
    m = RE_DUP_SEQ.search(name)
    if m:                           return "N", "N" + m.group(1).upper()
    if RE_DUP_NOSEQ.search(name):   return "N", "DUP"
    m = RE_INS.search(name)
    if m:                           return "N", "N" + m.group(1).upper()

    return None, None


def skip_reason(name):
    if re.match(r'^[A-Z0-9]+$', name): return "accession-only"
    if "Unknown" in name:               return "Unknown mutation"
    if RE_SUBST_NOALT.search(name):     return "truncated (missing ALT)"
    return "unparseable — review manually"


def parse_bed_files():
    files = sorted(glob.glob(os.path.join(BED_DIR, BED_GLOB)))
    print(f"Found {len(files)} BED files in {BED_DIR}\n")

    seen           = {}   # (chrom, pos, ref, alt, gene) → [acc, ...]
    variants       = []
    skipped        = []
    excluded_non_g = 0

    for fpath in files:
        fname = os.path.basename(fpath)
        m = re.match(r"(.+?)_(hg\d+)(?:_\d+)?\.BED$", fname, re.I)
        gene      = m.group(1) if m else fname
        hg_source = m.group(2) if m else "unknown"

        for line in open(fpath, encoding="utf-8"):
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue

            chrom  = parts[0]
            pos    = int(parts[1]) + 1
            name   = parts[3]
            strand = parts[5] if len(parts) > 5 else "."

            if not re.search(r"g\.", name, re.I):
                excluded_non_g += 1
                continue

            ref, alt = parse_ref_alt(name)
            acc = name.split("_")[0]

            if ref is None:
                skipped.append({"file": fname, "gene": gene, "hg_source": hg_source,
                                 "chrom": chrom, "pos": pos, "name": name,
                                 "reason": skip_reason(name)})
                continue

            key = (chrom, pos, ref, alt, gene)
            if key not in seen:
                seen[key] = [acc]
                variants.append({"chrom": chrom, "pos": pos, "name": name,
                                  "ref": ref, "alt": alt, "strand": strand,
                                  "source_file": fname, "gene": gene, "hg_source": hg_source})
            else:
                if acc not in seen[key]:
                    seen[key].append(acc)

    print(f"Non-g. rows excluded (c./r./p./bare): {excluded_non_g}")
    return variants, skipped, seen


def write_vcf(variants):
    with open(OUT_VCF, "w", encoding="utf-8") as f:
        f.write("\n".join([
            "##fileformat=VCFv4.2", "##reference=GRCh38/hg38",
            '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">',
            '##INFO=<ID=SRC,Number=1,Type=String,Description="Source BED file">',
            '##INFO=<ID=HG,Number=1,Type=String,Description="Original assembly">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ]) + "\n")
        for v in variants:
            f.write(f"{v['chrom']}\t{v['pos']}\t{v['name']}\t{v['ref']}\t{v['alt']}"
                    f"\t.\t.\tGENE={v['gene']};SRC={v['source_file']};HG={v['hg_source']}\n")
    print(f"VCF written      : {OUT_VCF}  ({len(variants)} variants)")


def write_accessions(seen):
    chrom_order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "M"])}
    rows = sorted(seen.items(), key=lambda x: (chrom_order.get(x[0][0], 99), x[0][1]))
    shared = sum(1 for accs in seen.values() if len(accs) > 1)
    with open(OUT_ACCESSIONS, "w", encoding="utf-8") as f:
        f.write("chrom\tpos\tref\talt\tgene\taccession_count\taccessions\n")
        for (chrom, pos, ref, alt, gene), accs in rows:
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{len(accs)}\t{','.join(accs)}\n")
    print(f"Accession map    : {OUT_ACCESSIONS}  ({len(seen)} variants, {shared} shared by >1 accession)")


def write_skipped(skipped):
    with open(OUT_SKIPPED, "w", encoding="utf-8") as f:
        f.write("file\tgene\thg_source\tchrom\tpos\tname\treason\n")
        for s in skipped:
            f.write(f"{s['file']}\t{s['gene']}\t{s['hg_source']}\t"
                    f"{s['chrom']}\t{s['pos']}\t{s['name']}\t{s['reason']}\n")
    print(f"Skipped log      : {OUT_SKIPPED}  ({len(skipped)} entries)")


def main():
    os.makedirs(_OUT_DIR, exist_ok=True)
    print("BED -> VCF  |  IDRefseq hg38\n")

    variants, skipped, seen = parse_bed_files()

    print(f"Variants parsed  : {len(variants)}")
    print(f"Skipped          : {len(skipped)}")
    for reason, count in Counter(s["reason"] for s in skipped).most_common():
        print(f"  {count:5d}  {reason}")
    print()

    chrom_order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "M"])}
    variants.sort(key=lambda v: (chrom_order.get(v["chrom"], 99), v["pos"]))

    write_vcf(variants)
    write_accessions(seen)
    write_skipped(skipped)


if __name__ == "__main__":
    main()
