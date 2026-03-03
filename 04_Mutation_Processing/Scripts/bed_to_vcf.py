"""
bed_to_vcf.py
=============
Converts hg38 LiftOver BED files into a VCF 4.2 file.

Usage:  python3 bed_to_vcf.py

Output (written to BED_DIR):
  all_variants_hg38.vcf   — VCF 4.2, deduplicated, sorted by chrom/pos
  bed_to_vcf_skipped.tsv  — variants that could not be parsed
"""

import os, re, glob
from collections import Counter

# ── Configuration ──────────────────────────────────────────────────────────────
BED_DIR  = "/home/khoi1996/Documents/BED/hg38"
BED_GLOB = "*_hg38_*.BED"

OUT_VCF     = os.path.join(BED_DIR, "all_variants_hg38.vcf")
OUT_SKIPPED = os.path.join(BED_DIR, "bed_to_vcf_skipped.tsv")

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


def parse_ref_alt(name):
    """Return (REF, ALT) or (None, None) if unparseable.
    Sentinels: ("N",".") = range del no seq; ("N","DUP") = dup no seq."""
    for pat, fn in [
        (RE_DELINS,      lambda m: ("N", m.group(1).upper())),
        (RE_INTRON_DEL,  lambda m: ("N" + m.group(1).upper(), "N")),
        (RE_INTRON_DINS, lambda m: ("N" + m.group(1).upper(), m.group(2).upper() or "N")),
        (RE_SUBST_G,     lambda m: (m.group(1).upper(), m.group(2).upper())),
        (RE_SUBST_C,     lambda m: (m.group(1).upper(), m.group(2).upper())),
    ]:
        m = pat.search(name)
        if m: return fn(m)

    if RE_SUBST_NOALT.search(name):  return None, None
    if RE_DOTDOT_DEL.search(name):   return "N", "."
    m = RE_DEL_SEQ.search(name)
    if m:                            return "N" + m.group(1).upper(), "N"
    if RE_DEL_NOSEQ.search(name):    return "N", "."
    m = RE_DUP_SEQ.search(name)
    if m:                            return "N", "N" + m.group(1).upper()
    if RE_DUP_NOSEQ.search(name):    return "N", "DUP"
    m = RE_INS.search(name)
    if m:                            return "N", "N" + m.group(1).upper()
    return None, None


def skip_reason(name):
    if re.match(r'^[A-Z0-9]+$', name):                                     return "accession-only"
    if "Unknown" in name:                                                   return "Unknown mutation"
    if RE_SUBST_NOALT.search(name):                                         return "truncated (missing ALT)"
    if re.search(r"_p\.", name) and not re.search(r"[gc]\.", name, re.I):  return "protein-only"
    if RE_PROT_NOPFX.match(name) and not re.search(r"[gc]\.", name, re.I): return "protein-only (no p. prefix)"
    if re.search(r"\br\.[a-z]", name) and not re.search(r"[gc]\.\S+[ACGT]>[ACGT]", name, re.I):
                                                                            return "RNA-only"
    return "unparseable — review manually"


def parse_bed_files(bed_dir, bed_glob):
    files = sorted(glob.glob(os.path.join(bed_dir, bed_glob)))
    print(f"Found {len(files)} hg38 BED files\n")
    seen, variants, skipped = set(), [], []

    for fpath in files:
        fname = os.path.basename(fpath)
        m = re.match(r"(.+)_hg38_(\d+)\.BED$", fname, re.I)
        gene, hg_source = (m.group(1), f"hg{m.group(2)}") if m else (fname, "unknown")

        for line in open(fpath, encoding="utf-8"):
            line = line.rstrip()
            if not line or line.startswith("#"): continue
            parts = line.split("\t")
            if len(parts) < 4: continue

            chrom  = parts[0]
            pos    = int(parts[1]) + 1       # BED 0-based → VCF 1-based
            name   = parts[3]
            strand = parts[5] if len(parts) > 5 else "."
            ref, alt = parse_ref_alt(name)

            if ref is None:
                skipped.append({"file": fname, "gene": gene, "hg_source": hg_source,
                                 "chrom": chrom, "pos": pos, "name": name,
                                 "reason": skip_reason(name)})
                continue

            key = (chrom, pos, ref, alt)
            if key in seen: continue
            seen.add(key)
            variants.append({"chrom": chrom, "pos": pos, "name": name,
                              "ref": ref, "alt": alt, "strand": strand,
                              "source_file": fname, "gene": gene, "hg_source": hg_source})
    return variants, skipped


def write_vcf(variants, out_path):
    with open(out_path, "w", encoding="utf-8") as f:
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
    print(f"VCF written      : {out_path}  ({len(variants)} variants)")


def write_skipped(skipped, out_path):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("file\tgene\thg_source\tchrom\tpos\tname\treason\n")
        for s in skipped:
            f.write(f"{s['file']}\t{s['gene']}\t{s['hg_source']}\t"
                    f"{s['chrom']}\t{s['pos']}\t{s['name']}\t{s['reason']}\n")
    print(f"Skipped log      : {out_path}  ({len(skipped)} entries)")


def main():
    print(f"BED → VCF  |  source: {BED_DIR}\n")
    variants, skipped = parse_bed_files(BED_DIR, BED_GLOB)

    print(f"Variants parsed  : {len(variants)}")
    print(f"Skipped          : {len(skipped)}")
    for reason, count in Counter(s["reason"] for s in skipped).most_common():
        print(f"  {count:5d}  {reason}")
    print()

    chrom_order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "M"])}
    variants.sort(key=lambda v: (chrom_order.get(v["chrom"], 99), v["pos"]))

    write_vcf(variants, OUT_VCF)
    write_skipped(skipped, OUT_SKIPPED)


if __name__ == "__main__":
    main()