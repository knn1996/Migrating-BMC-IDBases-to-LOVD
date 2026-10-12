"""
build_nmd_sample.py
===================
Stratified, deterministic sample of 20 variants from nmd_review_for_mauno.csv for curator
review: 10 predicted NMD and 10 escapes, tier A only. Named cases are taken first, then a
seeded draw fills each stratum with distinct genes and a balance of frameshift and nonsense
changes.

    predicted   JAK3 c.2788del (curator's example), TNFRSF13B c.581_582delinsAA (borderline),
                one start-proximal stop, one stop in the top decile of distance to the junction,
                then the rest to about 5 frameshift / 5 nonsense
    escapes     about 6 "stop in last coding exon" and 4 "within 50 nt"; ITGB2 c.2200G>T
                (borderline), CXCR4 c.1000C>T, IFNGR1 c.819_822del, one escape in a transcript
                whose last exon is entirely 3' UTR when there is one, then the rest to about
                half frameshift / half nonsense

Usage:
    python Step11_Export/build_nmd_sample.py [--in <nmd_review_for_mauno.csv>] [--out <nmd_sample_for_mauno.csv>] [--seed 20261013]
"""

import argparse
import csv
import random
from collections import Counter

import build_canonical_export as E

PRED = "predicted NMD: missing RNA and protein"
LAST = "escapes NMD: stop in last coding exon"
NEAR = "escapes NMD: within 50 nt of last junction"
UTR_LAST_EXON = {"NM_000073.3", "NM_000246.4", "NM_001770.6", "NM_002187.3", "NM_002524.5", "NM_002834.5", "NM_002872.5"}
FIXED_PRED = [("JAK3", "c.2788del", "curator's example"),
              ("TNFRSF13B", "c.581_582delinsAA", "borderline: stop 52 nt before the last coding exon junction")]
FIXED_ESC = [("ITGB2", "c.2200G>T", "borderline: stop 48 nt before the last coding exon junction"),
             ("CXCR4", "c.1000C>T", "known truncation mechanism (stop in the last coding exon)"),
             ("IFNGR1", "c.819_822del", "known truncation mechanism (stop just before the last coding exon junction)")]


def kind(row):
    return "frameshift" if "fs" in row["mutalyzer_p"] else "nonsense"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--in", dest="inp", default=str(E.OUT_ROOT / "Step11_Export" / "nmd_review_for_mauno.csv"))
    ap.add_argument("--out", default=str(E.OUT_ROOT / "Step11_Export" / "nmd_sample_for_mauno.csv"))
    ap.add_argument("--seed", type=int, default=20261013)
    args = ap.parse_args()
    with open(args.inp, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    cols = list(rows[0].keys()) + ["why_selected"]
    pool = [r for r in rows if r["confidence_tier"] == "A"]
    rng = random.Random(args.seed)
    picked, genes = [], set()

    def take(r, why):
        picked.append(dict(r, why_selected=why))
        genes.add(r["gene"])

    def find(g, c):
        return next(r for r in pool if r["gene"] == g and r["c_hgvs"] == c)

    def draw(cands, why, n=1, want=None):
        cands = sorted(cands, key=lambda r: (r["gene"], r["c_hgvs"]))
        rng.shuffle(cands)
        got = 0
        for r in cands:
            if got == n:
                break
            if r["gene"] in genes or (want and kind(r) not in want):
                continue
            take(r, why if isinstance(why, str) else why(r))
            got += 1
        return got

    pred = [r for r in pool if r["nmd_call"] == PRED]
    for g, c, why in FIXED_PRED:
        take(find(g, c), why)
    draw([r for r in pred if "start_proximal" in r["flags"]], "start-proximal stop (within the first 150 nt of the coding sequence)")
    dists = sorted(int(r["distance_to_junction_nt"]) for r in pred)
    top = dists[int(0.9 * len(dists))]
    draw([r for r in pred if int(r["distance_to_junction_nt"]) >= top], f"stop far from the last coding exon junction (top decile, >= {top} nt)")
    n_pred = 10
    while len(picked) < n_pred:
        c = Counter(kind(r) for r in picked)
        want = {"frameshift"} if c["frameshift"] < 5 else ({"nonsense"} if c["nonsense"] < 5 else None)
        if not draw(pred, lambda r: f"random draw ({kind(r)})", want=want):
            draw(pred, lambda r: f"random draw ({kind(r)})")
    for g, c, why in FIXED_ESC:
        take(find(g, c), why)
    esc = [r for r in pool if r["nmd_call"] in (LAST, NEAR)]
    utr = [r for r in esc if r["mane_transcript"] in UTR_LAST_EXON]
    utr_note = "" if utr else "no escaping PTC variant lies in a transcript whose last exon is entirely 3' UTR"
    if utr:
        draw(utr, "escape in a transcript whose last exon is entirely 3' UTR")
    while len(picked) < 20:
        e = [r for r in picked[n_pred:]]
        c = Counter(r["nmd_call"] for r in e)
        k = Counter(kind(r) for r in e)
        call = LAST if c[LAST] < 6 else NEAR
        want = {"frameshift"} if k["frameshift"] < 5 else ({"nonsense"} if k["nonsense"] < 5 else None)
        cands = [r for r in esc if r["nmd_call"] == call]
        if not draw(cands, lambda r: f"random draw ({kind(r)})", want=want):
            draw(cands, lambda r: f"random draw ({kind(r)})")
    picked.sort(key=lambda r: (0 if r["nmd_call"] == PRED else 1, r["gene"], r["c_hgvs"]))
    with open(args.out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(picked)
    calls = Counter(r["nmd_call"] for r in rows)
    sc = Counter(r["nmd_call"] for r in picked)
    print(f"All PTC variants: {len(rows)} ({calls[PRED]} predicted NMD, {calls[LAST]} escape in last coding exon, {calls[NEAR]} escape within 50 nt).")
    print(f"Sample of {len(picked)}: {sc[PRED]} predicted NMD, {sc[LAST]} escape in last coding exon, {sc[NEAR]} escape within 50 nt; "
          f"{Counter(kind(r) for r in picked)['frameshift']} frameshift / {Counter(kind(r) for r in picked)['nonsense']} nonsense; {len({r['gene'] for r in picked})} genes.")
    print(f"Drawn from tier-A variants only: named cases first (curator example, borderline, known mechanisms), then a seeded draw (seed {args.seed}) within each stratum, one gene per variant"
          + (f"; {utr_note}." if utr_note else "."))


if __name__ == "__main__":
    main()
