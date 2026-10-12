"""
build_nmd_review.py
===================
Curator table of the premature-stop (PTC) variants and their NMD prediction, one row per
distinct final variant (gene, record), from Step9_NMD/nmd_annotation.tsv and
transcript_structure.tsv. Offline.

Rule (Step9_NMD/annotate_nmd.py): a stop more than 50 nt upstream of the junction before the
last protein-coding exon is predicted to trigger NMD (missing RNA and protein: r.0?, p.0?); a
stop in the last coding exon, or within 50 nt of that junction, escapes. Variants whose stop
could not be assessed (intronic or non-coding descriptions, no stop before the normal one,
unparsed changes) are not listed; their count is printed.

Usage:
    python Step11_Export/build_nmd_review.py [--nmd-dir Output/Step9_NMD] [--out Output/Step11_Export/nmd_review_for_mauno.csv]
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import build_canonical_export as E

CALLS = {
    "NMD_predicted": "predicted NMD: missing RNA and protein",
    "NMD_escape_last_exon": "escapes NMD: stop in last coding exon",
    "NMD_escape_50nt_rule": "escapes NMD: within 50 nt of last junction",
    "NMD_escape_single_exon": "escapes NMD: single-exon transcript",
}
COLS = ["gene", "mane_transcript", "c_hgvs", "source_c", "source_p", "mutalyzer_p", "ptc_codon",
        "last_coding_exon_junction_c", "distance_to_junction_nt", "nmd_call", "proposed_r", "proposed_p",
        "flags", "n_patients", "accessions"]
E_THRESHOLD = 50
BORDERLINE = 5
START_PROXIMAL = 150


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--nmd-dir", type=Path, default=E.OUT_ROOT / "Step9_NMD")
    ap.add_argument("--out", type=Path, default=E.OUT_ROOT / "Step11_Export" / "nmd_review_for_mauno.csv")
    args = ap.parse_args()
    nmd = E.read_tsv(args.nmd_dir / "nmd_annotation.tsv")
    junction = {r["transcriptid"]: r.get("last_junction_c", "") for r in E.read_tsv(args.nmd_dir / "transcript_structure.tsv")}
    final = E.read_tsv(E.OUT_ROOT / "Step8_Rescue" / "lovd_flat_with_patients_final.tsv")
    pairs = {(r["individual_id"], r["dedup_key"]): r for r in E.read_tsv(E.OUT_ROOT / "Step8_Rescue" / "resolution_pairs.tsv")}
    units = E.Units(E.read_tsv(E.OUT_ROOT / "Step10_Audit" / "audit_variants.tsv"))
    alias = {r["gene"].strip().upper(): r["alias"].strip() for r in csv.DictReader(open(E.OUT_ROOT / "alias.csv", encoding="utf-8"))}
    members = defaultdict(list)
    for f in final:
        members[(f["_gene"], f["_dedup_key"])].append(f)
    rows, skipped = [], Counter()
    for r in nmd:
        call = r["nmd_call"]
        if call == "not_ptc":
            continue
        if call not in CALLS:
            skipped[call] += 1
            continue
        k = (r["_gene"], r["_dedup_key"])
        fs = members.get(k, [])
        src_c, src_p = [], []
        for f in fs:
            sid = pairs.get((f["Individual/ID"], f["_dedup_key"]), {}).get("source_change_ids", "").split(" ; ")[0]
            u = units.of(sid) if sid else None
            if u:
                src_c.append(u["src_c"].strip(" []"))
                src_p.append(u["src_p"].strip(" []"))
        dist = r["dist_to_last_junction_nt"]
        flags = [x for x in r["flags"].split(";") if x]
        if dist and abs(int(float(dist)) - E_THRESHOLD) <= BORDERLINE and "borderline" not in flags:
            flags.append("borderline")
        if r["ptc_c"] and int(float(r["ptc_c"])) <= START_PROXIMAL and "start_proximal" not in flags:
            flags.append("start_proximal")
        predicted = call == "NMD_predicted"
        rows.append({
            "gene": alias.get(r["_gene"].strip().upper(), r["_gene"]),
            "mane_transcript": r["transcriptid"], "c_hgvs": r["VariantOnTranscript/DNA"],
            "source_c": "; ".join(dict.fromkeys(src_c)), "source_p": "; ".join(dict.fromkeys(src_p)),
            "mutalyzer_p": r["protein_current"], "ptc_codon": r["ptc_aa"],
            "last_coding_exon_junction_c": junction.get(r["transcriptid"], ""),
            "distance_to_junction_nt": dist, "nmd_call": CALLS[call],
            "proposed_r": "r.0?" if predicted else r["rna_current"],
            "proposed_p": "p.0?" if predicted else r["protein_current"],
            "flags": ";".join(flags), "n_patients": len({f["Individual/ID"] for f in fs}),
            "accessions": ";".join(sorted({f["Individual/ID"].split("_", 1)[1] for f in fs})),
            "_order": 0 if predicted else 1,
        })
    rows.sort(key=lambda x: (x["_order"], x["gene"], x["c_hgvs"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"rows {len(rows)}:", dict(Counter(x["nmd_call"] for x in rows)))
    print("borderline:", sum("borderline" in x["flags"] for x in rows), "| start_proximal:", sum("start_proximal" in x["flags"] for x in rows))
    print("PTC variants not assessable (not listed):", dict(skipped))


if __name__ == "__main__":
    main()
