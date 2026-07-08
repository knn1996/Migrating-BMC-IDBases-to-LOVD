import os
import csv
from pathlib import Path
from collections import defaultdict

TRACK_A = os.environ["TRACK_A"]
TRACK_B = os.environ["TRACK_B"]
TRACK_C = os.environ["TRACK_C"]

OUT_MERGED     = os.environ["OUT_MERGED"]
OUT_UNRESOLVED = os.environ["OUT_UNRESOLVED"]
OUT_CROSSCHECK = os.environ["OUT_CROSSCHECK"]
LOG_PATH       = os.environ["LOG_PATH"]

for _d in (OUT_MERGED, LOG_PATH):
    os.makedirs(os.path.dirname(_d), exist_ok=True)

OUTPUT_COLS = [
    "gene", "accession", "allele_num", "sysname",
    "mut_type", "ref", "alt", "chrom", "pos_hg38", "strand",
    "status", "source_track",
    "nc_hgvs", "c_hgvs", "g_hgvs", "r_hgvs", "p_hgvs",
    "protein_pos_first", "mane_select", "mutalyzer_gene",
    "hgvs_input", "normalized", "errors",
]

CROSSCHECK_COLS = [
    "gene", "accession", "sysname", "field",
    "track_A_value", "track_B_value", "track_C_value",
    "AB_agree", "AC_agree", "BC_agree",
]

CHECK_FIELDS = ["c_hgvs", "g_hgvs", "r_hgvs", "p_hgvs"]

def normalise(row, track):
    out = {c: "" for c in OUTPUT_COLS}
    for col in OUTPUT_COLS:
        out[col] = row.get(col, "").strip()
    out["gene"]         = row.get("gene", "").strip().upper()
    out["source_track"] = track
    return out

def load_track(path, track_name):
    rows = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            n = normalise(row, track_name)
            key = (n["gene"], n["accession"], n["sysname"])
            if key not in rows:
                rows[key] = []
            rows[key].append(n)
    return rows

def best_ok(candidates):
    ok = [r for r in candidates if r["status"] == "ok"]
    return ok[0] if ok else None

def strip_version(hgvs):
    """Strip NM_/NC_ version numbers for loose comparison: NM_000022.4 -> NM_000022"""
    import re
    return re.sub(r'\.\d+', '', hgvs)

def main():
    log_lines = []
    def log(msg):
        print(msg)
        log_lines.append(msg)

    log("Loading tracks...")
    track_a = load_track(TRACK_A, "NG_IDRefseq")
    track_b = load_track(TRACK_B, "NM_MANE")
    track_c = load_track(TRACK_C, "NM_IDRefseq")

    log(f"  Track A (NG_IDRefseq): {sum(len(v) for v in track_a.values())} rows, {len(track_a)} unique keys")
    log(f"  Track B (NM_MANE):     {sum(len(v) for v in track_b.values())} rows, {len(track_b)} unique keys")
    log(f"  Track C (NM_IDRefseq): {sum(len(v) for v in track_c.values())} rows, {len(track_c)} unique keys")

    all_keys = sorted(set(list(track_a.keys()) + list(track_b.keys()) + list(track_c.keys())))
    log(f"  Total unique (gene, accession, sysname) keys: {len(all_keys)}")

    # ── Phase 1: Cross-check where multiple tracks have ok results ────────────
    log("\nPhase 1: Cross-checking agreement between tracks...")
    disagree_rows = []
    crosscheck_stats = defaultdict(lambda: {"compared": 0, "disagree": 0})

    for key in all_keys:
        ra = best_ok(track_a.get(key, []))
        rb = best_ok(track_b.get(key, []))
        rc = best_ok(track_c.get(key, []))

        ok_tracks = [(t, r) for t, r in [("A", ra), ("B", rb), ("C", rc)] if r is not None]
        if len(ok_tracks) < 2:
            continue

        for field in CHECK_FIELDS:
            vals = {t: r[field] for t, r in ok_tracks if r[field]}
            if len(vals) < 2:
                continue

            va = ra[field] if ra else ""
            vb = rb[field] if rb else ""
            vc = rc[field] if rc else ""

            va_s = strip_version(va)
            vb_s = strip_version(vb)
            vc_s = strip_version(vc)

            ab_agree = (va_s == vb_s) if va and vb else "N/A"
            ac_agree = (va_s == vc_s) if va and vc else "N/A"
            bc_agree = (vb_s == vc_s) if vb and vc else "N/A"

            track_pair_names = [t for t, _ in ok_tracks]
            all_agree = all(
                x is True for x in [ab_agree, ac_agree, bc_agree] if x != "N/A"
            )

            crosscheck_stats[field]["compared"] += 1
            if not all_agree:
                crosscheck_stats[field]["disagree"] += 1
                disagree_rows.append({
                    "gene":        key[0],
                    "accession":   key[1],
                    "sysname":     key[2],
                    "field":       field,
                    "track_A_value": va,
                    "track_B_value": vb,
                    "track_C_value": vc,
                    "AB_agree":    str(ab_agree),
                    "AC_agree":    str(ac_agree),
                    "BC_agree":    str(bc_agree),
                })

    with open(OUT_CROSSCHECK, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CROSSCHECK_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(disagree_rows)

    log(f"\n  Cross-check results (version-stripped comparison):")
    for field in CHECK_FIELDS:
        s = crosscheck_stats[field]
        log(f"    {field:<12}  compared: {s['compared']:>5}   disagreements: {s['disagree']:>4}")
    log(f"  Total disagreements written to: {OUT_CROSSCHECK}")

    # ── Phase 2: Priority merge ───────────────────────────────────────────────
    log("\nPhase 2: Priority merge (NM_MANE > NG_IDRefseq > NM_IDRefseq)...")
    merged     = []
    unresolved = []
    source_counts = {"NM_MANE": 0, "NG_IDRefseq": 0, "NM_IDRefseq": 0, "unresolved": 0}

    for key in all_keys:
        ok_b = best_ok(track_b.get(key, []))
        ok_a = best_ok(track_a.get(key, []))
        ok_c = best_ok(track_c.get(key, []))

        if ok_b:
            merged.append(ok_b)
            source_counts["NM_MANE"] += 1
        elif ok_a:
            merged.append(ok_a)
            source_counts["NG_IDRefseq"] += 1
        elif ok_c:
            merged.append(ok_c)
            source_counts["NM_IDRefseq"] += 1
        else:
            all_cands = track_b.get(key, []) + track_a.get(key, []) + track_c.get(key, [])
            if all_cands:
                row = all_cands[0]
            else:
                row = {c: "" for c in OUTPUT_COLS}
                row["gene"] = key[0]; row["accession"] = key[1]; row["sysname"] = key[2]
                row["status"] = "no_result"; row["source_track"] = "none"
            unresolved.append(row)
            source_counts["unresolved"] += 1

    with open(OUT_MERGED, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(merged)

    with open(OUT_UNRESOLVED, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(unresolved)

    log("")
    log("=== MERGE RESULTS ===")
    log(f"  Total merged (ok):                {len(merged)}")
    log(f"    from NM_MANE (priority 1):       {source_counts['NM_MANE']}")
    log(f"    from NG_IDRefseq (priority 2):   {source_counts['NG_IDRefseq']}")
    log(f"    from NM_IDRefseq (priority 3):   {source_counts['NM_IDRefseq']}")
    log(f"  Unresolved (all tracks failed):    {source_counts['unresolved']}")
    log(f"  Merged written to:    {OUT_MERGED}")
    log(f"  Unresolved written to: {OUT_UNRESOLVED}")
    log(f"  Crosscheck written to: {OUT_CROSSCHECK}")

    with open(LOG_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines))
    log(f"  Log: {LOG_PATH}")

if __name__ == "__main__":
    main()