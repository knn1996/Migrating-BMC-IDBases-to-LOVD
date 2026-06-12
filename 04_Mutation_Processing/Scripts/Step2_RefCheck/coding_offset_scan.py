import csv
import os
import re
from pathlib import Path
from collections import Counter

ALL_MUTATIONS = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step1_Extraction\all_mutations.tsv"
NM_DIR        = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\DNA sequences\Mane_Select_NM"
OUT_CSV       = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step2_RefCheck\coding_offset_scan_MANE.csv"

SCAN_MIN, SCAN_MAX = -1500, 1500
MIN_N        = 3
MIN_RELIABLE = 10
THRESH       = 0.9
NEAR         = 0.7
PEAK_GAP     = 0.3
STOPS        = {"TAA", "TAG", "TGA"}


def read_mrna_index(directory):
    idx = {}
    for p in Path(directory).iterdir():
        if p.suffix.lower() not in (".fasta", ".fa"):
            continue
        gene = p.stem.split("_", 1)[0].upper()
        seq = "".join(l.strip() for l in open(p, encoding="utf-8") if not l.startswith(">")).upper()
        idx[gene] = seq
    return idx


def find_cds(mrna):
    best = None
    n = len(mrna)
    for i in range(n - 2):
        if mrna[i:i + 3] != "ATG":
            continue
        j = i
        while j < n - 2:
            if mrna[j:j + 3] in STOPS:
                length = j + 3 - i
                if best is None or length > best[2]:
                    best = (i, j + 3, length)
                break
            j += 3
    if best is None:
        return None, None
    return mrna[best[0]:best[1]], best[0]


def load_subs():
    genes = {}
    with open(ALL_MUTATIONS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("variant_type", "").strip() != "coding":
                continue
            m = re.match(r"^c\.(\d+)([ACGT])>([ACGT])$", row["notation"].strip(), re.I)
            if not m:
                continue
            g = row["gene"].strip().upper()
            genes.setdefault(g, set()).add((int(m.group(1)), m.group(2).upper()))
    return genes


def scan(cds, subs):
    L = len(cds)

    def rate(d):
        ok = tot = 0
        for p, r in subs:
            i = (p - 1) + d
            if 0 <= i < L:
                tot += 1
                if cds[i] == r:
                    ok += 1
        return ok, tot

    results = []
    need = max(MIN_N, int(0.5 * len(subs)))
    for d in range(SCAN_MIN, SCAN_MAX + 1):
        ok, tot = rate(d)
        if tot >= need:
            results.append((ok / tot, ok, tot, d))
    if not results:
        return None
    results.sort(reverse=True)
    best = results[0]
    second = next((r for r in results[1:] if abs(r[3] - best[3]) > 2), (0.0, 0, 0, 0))
    ok0, tot0 = rate(0)
    m0 = ok0 / tot0 if tot0 else 0.0
    return best, second, m0


def classify(n_subs, m0, best_match, best_offset, second_match, n_tested):
    if m0 >= THRESH:
        return "aligned"
    if (best_match >= THRESH and best_offset != 0 and n_tested >= MIN_RELIABLE
            and (best_match - second_match) >= PEAK_GAP):
        return "recoverable_offset"
    if n_subs < MIN_RELIABLE:
        return "inconclusive_low_n"
    if m0 >= NEAR:
        return "near_aligned"
    return "isoform_or_mismatch"


def main():
    mrna_idx = read_mrna_index(NM_DIR)
    subs_by_gene = load_subs()

    fields = ["gene", "n_subs", "n_tested", "match_at_0", "best_offset", "best_match",
              "second_match", "cds_start", "cds_len", "orf_ok",
              "mane_c_equals_idbases_c_plus", "classification"]
    rows = []

    for gene in sorted(subs_by_gene):
        subs = subs_by_gene[gene]
        rec = {f: "" for f in fields}
        rec["gene"] = gene
        rec["n_subs"] = len(subs)

        if gene not in mrna_idx:
            rec["classification"] = "no_mane_fasta"
            rows.append(rec); continue

        cds, cds_start = find_cds(mrna_idx[gene])
        if not cds:
            rec["classification"] = "no_orf"
            rows.append(rec); continue

        sc = scan(cds, subs)
        if sc is None:
            rec["classification"] = "inconclusive_low_n"
            rows.append(rec); continue

        (bm, _, bt, bd), (sm, _, _, _), m0 = sc
        orf_ok = len(cds) % 3 == 0 and cds.startswith("ATG") and cds[-3:] in STOPS
        cls = classify(len(subs), m0, bm, bd, sm, bt)

        rec.update({"n_tested": bt, "match_at_0": round(m0, 3), "best_offset": bd,
                    "best_match": round(bm, 3), "second_match": round(sm, 3),
                    "cds_start": cds_start, "cds_len": len(cds), "orf_ok": orf_ok,
                    "mane_c_equals_idbases_c_plus": bd, "classification": cls})
        rows.append(rec)

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    counts = Counter(r["classification"] for r in rows)
    print(f"Genes scanned: {len(rows)}  -> {OUT_CSV}")
    for k, v in counts.most_common():
        print(f"  {k:22} {v}")

    rec_off = [r for r in rows if r["classification"] == "recoverable_offset"]
    if rec_off:
        print("\nRECOVERABLE OFFSETS (apply MANE_c = IDbases_c + offset):")
        for r in sorted(rec_off, key=lambda r: r["best_offset"]):
            print(f"  {r['gene']:10} offset={r['best_offset']:+5}  "
                  f"{r['match_at_0']:.0%} -> {r['best_match']:.0%}  (n={r['n_tested']})")

    near = [r for r in rows if r["classification"] == "near_aligned"]
    if near:
        print("\nNEAR-ALIGNED (no offset; ~match_at_0 resolve per-variant, rest logged):")
        for r in sorted(near, key=lambda r: -r["match_at_0"]):
            print(f"  {r['gene']:10} match_at_0={r['match_at_0']:.0%}  (n={r['n_tested']})")

    iso = [r for r in rows if r["classification"] == "isoform_or_mismatch"]
    if iso:
        print("\nISOFORM / MISMATCH (large n, no offset, <70% at 0):")
        for r in sorted(iso, key=lambda r: -r["n_tested"]):
            print(f"  {r['gene']:10} match_at_0={r['match_at_0']:.0%}  best={r['best_match']:.0%}@{r['best_offset']:+}  (n={r['n_tested']})")


if __name__ == "__main__":
    main()
