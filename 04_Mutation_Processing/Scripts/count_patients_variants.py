"""
count_patients_variants.py
Count unique patients, variants, and zygosity categories per gene
from IDbases *pub.html files.
"""

import os
import re
import csv
from pathlib import Path

IDBASE_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\idbase_pub_html_only"
OUT_CSV    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\idbase_patient_variant_counts.csv"

X_LINKED = {
    "BTK", "CD40LG", "CYBB", "DKC1", "ELA2", "FOXP3", "G6PD", "GATA1",
    "HPRT1", "IKBKG", "IL2RG", "MAGT1", "MECP2", "MSN", "NEMO", "ORAI1",
    "PIGA", "PROP1", "RPL10", "SASH3", "SH2D1A", "TAZ", "TLR7", "TLR8",
    "WAS", "WASP", "XIAP", "ATP6AP1", "FLNA", "CFP", "SPTLC2", "DOCK11",
    "TFE3", "EBP",
}

EXCLUDE_GENES = {"BTK", "SH2"}

RE_STRIP  = re.compile(r"<[^>]+>")
RE_GENE   = re.compile(r"^Gene\s+(.+)$")
RE_ACC    = re.compile(r"^Accession\s+(\S+)")
RE_SYSNAME = re.compile(r"^Systematic name\s+(.+)$")
RE_SEX    = re.compile(r"^Sex\s+(\S+)")


def strip_html(text):
    return RE_STRIP.sub("", text)


def gene_from_filename(path):
    stem = Path(path).stem.lower()
    for suffix in ("_pub", "pub", "_base", "base"):
        while stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    return stem.upper()


def is_excluded(filename, gene):
    if gene_from_filename(filename) in EXCLUDE_GENES:
        return True
    if gene and gene.upper() in EXCLUDE_GENES:
        return True
    return False


def parse_file(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    text = strip_html(text)

    gene = None
    for line in text.splitlines():
        m = RE_GENE.match(line.strip())
        if m:
            gene = m.group(1).strip().split()[0]
            break
    if gene is None:
        gene = gene_from_filename(path)

    records = re.split(r"^//\s*$", text, flags=re.MULTILINE)

    patients = []
    for rec in records:
        acc = None
        sysnames = []
        sex = None
        for line in rec.splitlines():
            s = line.strip()
            if acc is None:
                m = RE_ACC.match(s)
                if m:
                    acc = m.group(1)
                    continue
            m = RE_SYSNAME.match(s)
            if m:
                sysnames.append(m.group(1).strip())
                continue
            m = RE_SEX.match(s)
            if m:
                sex = m.group(1).strip()
        if acc:
            patients.append({"accession": acc, "sysnames": sysnames, "sex": sex})
    return gene, patients


def classify_zygosity(patient, is_x_linked):
    sysnames = patient["sysnames"]
    sex = patient["sex"]
    n = len(sysnames)

    if n == 0:
        return "no_variant"

    if is_x_linked and sex and sex.upper() == "XY":
        return "hemizygous"

    if n == 1:
        return "heterozygous"

    unique_hgvs = set()
    for sn in sysnames:
        hgvs = re.sub(r"^Allele\s*\d+\s*:\s*", "", sn, flags=re.IGNORECASE).strip()
        unique_hgvs.add(hgvs)

    if len(unique_hgvs) == 1:
        return "homozygous"
    return "heterozygous"


def main():
    rows = []
    skipped = []
    files = sorted(Path(IDBASE_DIR).glob("*pub.html"))

    for fp in files:
        gene, patients = parse_file(fp)
        if is_excluded(fp.name, gene):
            skipped.append((fp.name, gene))
            continue

        is_x = gene.upper() in X_LINKED
        unique_patients = {p["accession"] for p in patients}

        all_sysnames    = [sn for p in patients for sn in p["sysnames"]]
        total_variants  = len(all_sysnames)
        unique_variants = len(set(all_sysnames))

        zyg_counts = {"homozygous": 0, "heterozygous": 0, "hemizygous": 0,
                      "no_variant": 0}
        sex_counts = {"XX": 0, "XY": 0, "other_or_unknown": 0}

        seen_acc = set()
        for p in patients:
            if p["accession"] in seen_acc:
                continue
            seen_acc.add(p["accession"])
            z = classify_zygosity(p, is_x)
            zyg_counts[z] = zyg_counts.get(z, 0) + 1
            s = (p["sex"] or "").upper()
            if s == "XX":
                sex_counts["XX"] += 1
            elif s == "XY":
                sex_counts["XY"] += 1
            else:
                sex_counts["other_or_unknown"] += 1

        rows.append({
            "gene": gene,
            "file": fp.name,
            "total_patients": len(unique_patients),
            "total_variants": total_variants,
            "unique_variants": unique_variants,
            "homozygous": zyg_counts["homozygous"],
            "heterozygous": zyg_counts["heterozygous"],
            "hemizygous": zyg_counts["hemizygous"],
            "no_variant": zyg_counts["no_variant"],
            "sex_XX": sex_counts["XX"],
            "sex_XY": sex_counts["XY"],
            "sex_other_or_unknown": sex_counts["other_or_unknown"],
            "is_x_linked": "yes" if is_x else "no",
        })

    rows.sort(key=lambda r: r["gene"].upper())

    total = {
        "gene": "TOTAL",
        "file": f"{len(rows)} files",
        "total_patients": sum(r["total_patients"] for r in rows),
        "total_variants": sum(r["total_variants"] for r in rows),
        "unique_variants": sum(r["unique_variants"] for r in rows),
        "homozygous":   sum(r["homozygous"]   for r in rows),
        "heterozygous": sum(r["heterozygous"] for r in rows),
        "hemizygous":   sum(r["hemizygous"]   for r in rows),
        "no_variant":   sum(r["no_variant"]   for r in rows),
        "sex_XX": sum(r["sex_XX"] for r in rows),
        "sex_XY": sum(r["sex_XY"] for r in rows),
        "sex_other_or_unknown": sum(r["sex_other_or_unknown"] for r in rows),
        "is_x_linked": "",
    }
    rows.append(total)

    fieldnames = ["gene", "file", "total_patients", "total_variants", "unique_variants",
                  "homozygous", "heterozygous", "hemizygous", "no_variant",
                  "sex_XX", "sex_XY", "sex_other_or_unknown", "is_x_linked"]
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {OUT_CSV}")
    print(f"Files processed: {len(rows)-1}")
    print(f"Files skipped (excluded): {len(skipped)}")
    for fn, g in skipped:
        print(f"  - {fn} (gene={g})")
    print(f"Total patients:  {total['total_patients']}")
    print(f"Total variants:  {total['total_variants']}")
    print(f"Unique variants: {total['unique_variants']}")
    print(f"  homozygous:    {total['homozygous']}")
    print(f"  heterozygous:  {total['heterozygous']}")
    print(f"  hemizygous:    {total['hemizygous']}")


if __name__ == "__main__":
    main()
