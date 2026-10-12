"""
flag_review.py
==============
Last step on the resolved records before the final table. Records are never remapped or removed.

    review     a pair goes to confidence_tier Review, with the reason in review_reason, when
               (a) its record's genomic description projects outside the MANE Select transcript
                   or more than 20 nt into an intron (MANE v1.5 GFF) while a source c. of the pair
                   is exonic, or
               (b) its record keeps the literal source c. (shift 0) in a gene with a verified
                   legacy shift (refs/legacy_shift_evidence.tsv), the source c. at that shift
                   places the reference residue of a source p., and the pair is inconsistent;
                   the note gives the c. expected at the shift, or
               (c) its record's c. is "=" and its g. is "=" or absent: no change is described.

Env vars
--------
IN_DEDUP, IN_PAIRS        normalize_mane outputs
OUT_DEDUP, OUT_PAIRS      Step8_Rescue/dedup_merged_variants.resolved.tsv, resolution_pairs.tsv
OUT_LOG                   Step8_Rescue/review_log.tsv
MANE_GFF, MANE_NM_DIR, GB_CACHE, SHIFT_EVIDENCE, PUB_DIR, AUDIT_SCRIPT
"""

import csv
import importlib.util
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "Utilities"))
import normalize_mane as N
import legacy_remap as R

LOG_COLS = ["gene", "dedup_key", "individual_id", "rule", "tier_before", "reason"]
RE_SEG = re.compile(r"^[\[\s]*(?:([gcrp])\s*\.|(?=IVS))", re.I)
RE_EXONIC = re.compile(r"^c\.\d+(?:_\d+)?[A-Za-z>]")


def load_audit(path):
    spec = importlib.util.spec_from_file_location("audit_against_source", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def top_level_segments(txt):
    out, depth, cur = [], 0, []
    for ch in txt:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        if ch == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return [s.strip() for s in out if s.strip()]


def source_proteins(txt):
    out = []
    for seg in top_level_segments(txt):
        m = RE_SEG.match(seg)
        if not m or (m.group(1) or "g").lower() != "p":
            continue
        for comp in re.split(r";|\s\+\s", seg[m.end():]):
            comp = re.sub(r"^[gcrp]\.", "", re.sub(r"[\[\]\s]", "", comp), flags=re.I)
            if comp:
                out.append("p." + comp)
    return out


def load_source_proteins(audit, pub_dir):
    proteins = defaultdict(list)
    seen = set()
    for name in sorted(f for f in os.listdir(pub_dir) if f.lower().endswith("pub.html")):
        gene = audit.gene_from_filename(name)
        records, _ = audit.read_records(os.path.join(pub_dir, name))
        for lines in records:
            acc, fields, _ = audit.parse_record(lines)
            if not acc or (gene, acc) in seen:
                continue
            seen.add((gene, acc))
            for _, txt in audit.split_alleles(fields):
                proteins[f"{gene}_{acc}"].extend(source_proteins(txt))
    return proteins


def pair_inconsistent(final_p, src_ps, audit):
    if not final_p.strip() or not src_ps:
        return False
    res = [audit.compare_protein(s, final_p) for s in src_ps]
    return not any(st in audit.GOOD_P for st, _ in res) and any(st in audit.BAD_P for st, _ in res)


def genbank_cds(path):
    text = path.read_text(encoding="utf-8")
    m = re.search(r"^ {5}CDS {13}<?(\d+)\.\.>?(\d+)", text, re.M)
    if not m or "\nORIGIN" not in text:
        return None
    seq = re.sub(r"[^acgtnACGTN]", "", text.split("\nORIGIN", 1)[1].split("\n//")[0]).upper()
    return seq, int(m.group(1)) - 1


def mane_cds(nm, gb_dir, mane_dir):
    gb = Path(gb_dir) / f"{nm}.gb"
    if gb.exists():
        x = genbank_cds(gb)
        if x:
            return x
    seq = N.fasta_seq(mane_dir, nm)
    cds = N.longest_cds(seq) if seq else ""
    return (seq, seq.find(cds)) if cds else None


def genomic_of(r):
    g = r.get("g_hgvs", "") or ""
    if g:
        return g
    n = r.get("normalized", "") or ""
    return n if n.startswith("NC_") and ":g." in n else ""


def coding_of(r):
    c = r.get("c_hgvs", "") or r.get("normalized", "") or ""
    return c.rsplit(":", 1)[-1] if ":c." in c else ""


def is_no_change(r):
    c = r.get("c_hgvs", "") or r.get("normalized", "") or ""
    g = genomic_of(r)
    return c.endswith(":c.=") and (not g or g.endswith(":g.="))


def read_shifts(path):
    return {r["gene"].strip().upper(): int(r["shift"]) for r in N.read_tsv(path) if r.get("status") == "kept" and r.get("shift")}


def literal_placement(gene, rec_c, src_c, shift, prots, cds_of):
    rc, sc = R.parse_c(rec_c), R.parse_c(src_c)
    if not rc or not sc or (rc[1], rc[2]) != (sc[1], sc[2]):
        return False
    ref = cds_of(gene)
    if not ref:
        return False
    seq, cds0 = ref
    if R.mutate(seq, cds0, sc, shift) is None:
        return False
    wt = R.translate(seq[cds0:])
    for p in (R.parse_p(x) for x in prots):
        if p and (sc[1] + shift - 1) // 3 - 1 <= p[1] <= (sc[2] + shift) // 3 + 2 and 0 < p[1] <= len(wt) and wt[p[1] - 1] == p[2]:
            return True
    return False


def main():
    env = os.environ
    dedup = N.read_tsv(env["IN_DEDUP"])
    cols = list(dedup[0].keys()) if dedup else []
    pairs = N.read_tsv(env["IN_PAIRS"])
    pair_cols = (list(pairs[0].keys()) if pairs else []) + ["review_reason"]
    mane = N.mane_transcripts(env["MANE_GFF"])
    shifts = read_shifts(env["SHIFT_EVIDENCE"])
    audit = load_audit(env["AUDIT_SCRIPT"])
    src_p = load_source_proteins(audit, env["PUB_DIR"])
    cds_cache = {}

    def mane_of(r):
        return mane.get((r.get("mutalyzer_gene") or r["gene"]).strip()) or mane.get(r["gene"].strip().upper())

    def cds_of(gene_row):
        t = mane_of(gene_row)
        if not t:
            return None
        if t["nm"] not in cds_cache:
            cds_cache[t["nm"]] = mane_cds(t["nm"], env["GB_CACHE"], env["MANE_NM_DIR"])
        return cds_cache[t["nm"]]

    records = {(r["gene"], r["dedup_key"]): r for r in dedup}
    no_change = {k for k, r in records.items() if is_no_change(r)}
    log, counts = [], Counter()
    out_pairs = []
    for p in pairs:
        k = (p["gene"], p["dedup_key"])
        r = records.get(k, {})
        sids = [s.strip() for s in p.get("source_change_ids", "").split(" ; ") if s.strip()]
        src_cs = [s.split("|", 4)[4] for s in sids if s.count("|") >= 4]
        reasons = [("no_change", "no change described (c.=)")] if k in no_change else []
        t = mane_of(r) if r else None
        g = genomic_of(r) if r else ""
        if t and g and any(RE_EXONIC.match(c) for c in src_cs):
            _, why = N.project(t, g)
            if why in ("outside the MANE transcript", "deep intronic on MANE"):
                reasons.append(("genomic", f"genomic record {g} lies {why.replace(' on MANE', '')} on {t['nm']}, the source c. ({', '.join(src_cs)}) is exonic"))
        gene = p["gene"].strip().upper()
        if r and gene in shifts:
            rec_c = coding_of(r)
            prots = src_p.get(p["individual_id"], [])
            final_p = r.get("p_hgvs", "").split(":", 1)[1] if ":" in r.get("p_hgvs", "") else ""
            for c in src_cs:
                if literal_placement(gene, rec_c, c, shifts[gene], prots, lambda _g: cds_of(r)) and pair_inconsistent(final_p, prots, audit):
                    expected = R.shift_notation(c, shifts[gene]) or "?"
                    reasons.append(("literal", f"record keeps the literal source numbering {c}; at the {gene} legacy shift {shifts[gene]:+d} the expected c. is {expected}"))
                    break
        p = dict(p, review_reason="")
        if reasons:
            log.append({"gene": p["gene"], "dedup_key": p["dedup_key"], "individual_id": p["individual_id"],
                        "rule": "+".join(x for x, _ in reasons), "tier_before": p["confidence_tier"], "reason": " ; ".join(y for _, y in reasons)})
            p["confidence_tier"] = "Review"
            p["review_reason"] = " ; ".join(y for _, y in reasons)
            for x, _ in reasons:
                counts[f"review_{x}"] += 1
        out_pairs.append(p)
    N.write_tsv(env["OUT_DEDUP"], dedup, cols)
    N.write_tsv(env["OUT_PAIRS"], out_pairs, pair_cols)
    N.write_tsv(env["OUT_LOG"], log, LOG_COLS)
    print(f"pairs to Review: {len(log)} (" + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())) + ")")
    print("tier before:", dict(Counter(e["tier_before"] for e in log)))
    print("tiers now:", dict(Counter(p["confidence_tier"] for p in out_pairs)))


if __name__ == "__main__":
    main()
