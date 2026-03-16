"""
bed_to_vcf.py
=============
Converts hg38 LiftOver BED files into a VCF 4.2 file.

Input  : BED files produced by generate_bed.py → UCSC LiftOver.
         Mutation NAME fields contain sysname_per_allele values
         from Output/all_mutations/all_mutations.tsv.

         Source 1 (IDRefSeq) : 03_BED_Files/BED_source_1/hg38        (*_hg38_*.BED)
         Source 3 (LRG)      : 03_BED_Files/BED_LRG_source_3/bed_hg38/      	(*_hg38_*.BED)

Output : written to 04_Mutation_Processing/Output/
           Source 1 → all_variants_hg38_source_1.vcf        / bed_to_vcf_skipped_source_1.tsv
           Source 3 → all_variants_hg38_source_3_LRG.vcf    / bed_to_vcf_skipped_source_3_LRG.tsv

Usage  : python3 bed_to_vcf.py --source 1
         python3 bed_to_vcf.py --source 3
"""

import argparse
import os
import re
import glob
from collections import Counter

# ── Configuration ──────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_THESIS_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", ".."))
_OUT_DIR    = os.path.join(_SCRIPT_DIR, "..", "..", "Output", "bed_to_vcf")

BED_GLOB = "*_hg38_*.BED"

SOURCE_CONFIG = {
    1: {
        "bed_dir":     os.path.join(_THESIS_DIR, "03_BED_Files", "BED_source_1", "hg38"),
        "out_vcf":     os.path.join(_OUT_DIR, "all_variants_hg38_source_1.vcf"),
        "out_skipped": os.path.join(_OUT_DIR, "bed_to_vcf_skipped_source_1.tsv"),
        "label":       "Source 1 (IDRefSeq)",
    },
    3: {
        "bed_dir":     os.path.join(_THESIS_DIR, "03_BED_Files", "BED_LRG_source_3", "bed_hg38"),
        "out_vcf":     os.path.join(_OUT_DIR, "all_variants_hg38_source_3_LRG.vcf"),
        "out_skipped": os.path.join(_OUT_DIR, "bed_to_vcf_skipped_source_3_LRG.tsv"),
        "label":       "Source 3 (LRG)",
    },
}

# ── Regexes ────────────────────────────────────────────────────────────────────
RE_SUBST_G     = re.compile(r"g\.(?:IVS[\d_+\-]+|[-+*\d]+)([ACGT])>([ACGT])", re.I)
RE_SUBST_C     = re.compile(r"c\.(?:IVS[\d_+\-]+|[-+*\d]+)([ACGT])>([ACGT])", re.I)
RE_SUBST_NOALT = re.compile(r"[gc]\.[^\t\n]+[ACGT]>(?:_|\s*$)", re.I)
RE_DELINS      = re.compile(r"[gc]\.\d+(?:_\d+)?delins([ACGT]+)", re.I)
RE_INTRON_DEL  = re.compile(r"[gc]\.\d+[-+]\d+_\d+del([ACGTacgt]+)", re.I)
RE_INTRON_DINS = re.compile(r"[gc]\.\d+[-+]\d+([ACGT]{2,})>([ACGT]*)", re.I)
RE_DEL_SEQ     = re.compile(r"[gc]\.\d+(?:_\d+)?del([ACGT]+)", re.I)
RE_DEL_NOSEQ   = re.compile(r"[gc]\.\d+(?:_\d+)?del(?![ACGT])", re.I)
RE_DOTDOT_DEL  = re.compile(r"g\.(\d+)\.\.(\d+)del", re.I)
RE_DUP_SEQ     = re.compile(r"[gc]\.\d+(?:_\d+)?dup([ACGT]+)", re.I)
RE_DUP_NOSEQ   = re.compile(r"[gcr]\.\d+(?:_\d+)?dup(?![ACGT])", re.I)
RE_INS         = re.compile(r"[gc]\.\d+_\d+ins([ACGT]+)", re.I)
RE_PROT_NOPFX  = re.compile(r"^[A-Z0-9]+_[A-Z][a-z]{2}\d+")


# ── Parsers ────────────────────────────────────────────────────────────────────

def parse_ref_alt(name):
    """Return (REF, ALT) or (None, None) if unparseable.
    Sentinels: ("N", ".")   = range deletion, no explicit sequence
               ("N", "DUP") = duplication, no explicit sequence"""
    for pat, fn in [
        (RE_DELINS,      lambda m: ("N", m.group(1).upper())),
        (RE_INTRON_DEL,  lambda m: ("N" + m.group(1).upper(), "N")),
        (RE_INTRON_DINS, lambda m: ("N" + m.group(1).upper(), m.group(2).upper() or "N")),
        (RE_SUBST_G,     lambda m: (m.group(1).upper(), m.group(2).upper())),
        (RE_SUBST_C,     lambda m: (m.group(1).upper(), m.group(2).upper())),
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
    """Reason why a name field could not be parsed."""
    if re.match(r'^[A-Z0-9]+$', name):
        return "accession-only"
    if "Unknown" in name:
        return "Unknown mutation"
    if RE_SUBST_NOALT.search(name):
        return "truncated (missing ALT)"
    if re.search(r"_p\.", name) and not re.search(r"[gc]\.", name, re.I):
        return "protein-only"
    if RE_PROT_NOPFX.match(name) and not re.search(r"[gc]\.", name, re.I):
        return "protein-only (no p. prefix)"
    if re.search(r"\br\.[a-z]", name) and not re.search(r"[gc]\.\S+[ACGT]>[ACGT]", name, re.I):
        return "RNA-only"
    return "unparseable — review manually"


# ── I/O ───────────────────────────────────────────────────────────────────────

def parse_bed_files(bed_dir, bed_glob):
    files = sorted(glob.glob(os.path.join(bed_dir, bed_glob)))
    print(f"Found {len(files)} hg38 BED files in {bed_dir}\n")
    seen, variants, skipped = set(), [], []

    for fpath in files:
        fname = os.path.basename(fpath)
        m = re.match(r"(.+)_hg38_(\d+)\.BED$", fname, re.I)
        gene      = m.group(1)       if m else fname
        hg_source = f"hg{m.group(2)}" if m else "unknown"

        with open(fpath, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue

                chrom  = parts[0]
                pos    = int(parts[1]) + 1   # BED 0-based → VCF 1-based
                name   = parts[3]
                strand = parts[5] if len(parts) > 5 else "."
                ref, alt = parse_ref_alt(name)

                if ref is None:
                    skipped.append({
                        "file": fname, "gene": gene, "hg_source": hg_source,
                        "chrom": chrom, "pos": pos, "name": name,
                        "reason": skip_reason(name),
                    })
                    continue

                key = (chrom, pos, ref, alt)
                if key in seen:
                    continue
                seen.add(key)
                variants.append({
                    "chrom": chrom, "pos": pos, "name": name,
                    "ref": ref, "alt": alt, "strand": strand,
                    "source_file": fname, "gene": gene, "hg_source": hg_source,
                })

    return variants, skipped


def write_vcf(variants, out_path):
    VCF_HEADER = (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38/hg38\n"
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">\n'
        '##INFO=<ID=SRC,Number=1,Type=String,Description="Source BED file">\n'
        '##INFO=<ID=HG,Number=1,Type=String,Description="Original assembly">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(VCF_HEADER)
        for v in variants:
            fh.write(
                f"{v['chrom']}\t{v['pos']}\t{v['name']}\t{v['ref']}\t{v['alt']}"
                f"\t.\t.\tGENE={v['gene']};SRC={v['source_file']};HG={v['hg_source']}\n"
            )
    print(f"VCF written      : {out_path}  ({len(variants)} variants)")


def write_skipped(skipped, out_path):
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("file\tgene\thg_source\tchrom\tpos\tname\treason\n")
        for s in skipped:
            fh.write(
                f"{s['file']}\t{s['gene']}\t{s['hg_source']}\t"
                f"{s['chrom']}\t{s['pos']}\t{s['name']}\t{s['reason']}\n"
            )
    print(f"Skipped log      : {out_path}  ({len(skipped)} entries)")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Convert hg38 LiftOver BED files to VCF 4.2.")
    parser.add_argument(
        "--source", type=int, choices=[1, 3], required=True,
        help="Reference source: 1 = IDRefSeq (BED/hg38/), 3 = LRG (BED (LRG_Source 3)/bed_hg38/)",
    )
    args = parser.parse_args()

    cfg = SOURCE_CONFIG[args.source]
    print(f"BED → VCF  |  {cfg['label']}\n")

    variants, skipped = parse_bed_files(cfg["bed_dir"], BED_GLOB)

    print(f"Variants parsed  : {len(variants)}")
    print(f"Skipped          : {len(skipped)}")
    for reason, count in Counter(s["reason"] for s in skipped).most_common():
        print(f"  {count:5d}  {reason}")
    print()

    chrom_order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "M"])}
    variants.sort(key=lambda v: (chrom_order.get(v["chrom"], 99), v["pos"]))

    write_vcf(variants, cfg["out_vcf"])
    write_skipped(skipped, cfg["out_skipped"])


if __name__ == "__main__":
    main()
