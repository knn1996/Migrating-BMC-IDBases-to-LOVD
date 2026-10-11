"""
resolve_tracks.py
=================
Resolves conflicting track descriptions of one source change of one individual
(rule R2) and assigns a confidence tier to every final (individual, variant) pair.

A pair of the rescue-deduped table belongs to the source changes its augmented
rows carry (source_change_id). When one source change of one individual is
described by two or more final records, one record is kept:

    1. protein: a record whose protein is consistent with any source protein
       change of the individual (audit compare_protein, GOOD_P) beats the others;
       an inconsistent record (BAD_P) loses to an uncheckable one
    2. c_region: the record's c. region (exonic, intronic, 5' UTR, 3' UTR) matches
       the source c. notation of the change
    3. mane_select: the record is on the MANE Select transcript
    4. gene_prior: the record's track is protein-consistent more often in the other
       protein-decidable conflicts of the same gene (leave-one-out, Laplace)
    5. track_priority: NM_MANE > NG_IDRefseq > NM_IDRefseq > others, then key

A dropped record leaves the individual only if every source change it describes
is such a conflict. When all records of a conflict describe the same genomic
change the dropped descriptions are kept as alt_descriptions of the kept pair.

confidence_tier per pair: A (no conflict, or decided by protein), B (decided by
signals 2-5), Review (different genomic positions and every record inconsistent,
or every record uncheckable). A pair takes its worst tier over its conflicts.

Env vars
--------
DEDUP_TSV          Step8_Rescue/dedup_merged_variants.rescued.tsv
AUGMENT_TSV        Step8_Rescue/augmented_merged_variants.tsv (with source_change_id)
PUB_DIR            IDbases pub.html directory (source protein changes)
AUDIT_SCRIPT       Step10_Audit/audit_against_source.py (protein comparison, source parsing)
OUT_DEDUP          Step8_Rescue/dedup_merged_variants.resolved.tsv
OUT_PAIRS          Step8_Rescue/resolution_pairs.tsv
OUT_LOG            Step8_Rescue/resolution_log.tsv
OUT_STATS          Step8_Rescue/resolution_stats.txt
"""

import csv
import importlib.util
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

DEDUP_TSV    = os.environ["DEDUP_TSV"]
AUGMENT_TSV  = os.environ["AUGMENT_TSV"]
PUB_DIR      = os.environ["PUB_DIR"]
AUDIT_SCRIPT = os.environ.get("AUDIT_SCRIPT", str(Path(__file__).resolve().parents[1] / "Step10_Audit" / "audit_against_source.py"))
OUT_DEDUP    = os.environ["OUT_DEDUP"]
OUT_PAIRS    = os.environ["OUT_PAIRS"]
OUT_LOG      = os.environ["OUT_LOG"]
OUT_STATS    = os.environ["OUT_STATS"]

TRACK_PRIORITY = {"NM_MANE": 0, "NG_IDRefseq": 1, "NM_IDRefseq": 2}
RE_SEG = re.compile(r"^[\[\s]*(?:([gcrp])\s*\.|(?=IVS))", re.I)
TIER_RANK = {"A": 0, "B": 1, "Review": 2}
PAIR_COLS = ["gene", "accession", "individual_id", "dedup_key", "confidence_tier", "alt_descriptions", "source_change_ids"]
LOG_COLS = ["source_change_id", "individual_id", "kept", "dropped", "deciding_signal", "tier", "group_kind",
            "labels", "tracks"]


def load_audit():
    spec = importlib.util.spec_from_file_location("audit_against_source", AUDIT_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_tsv(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def variant_desc(s):
    return s.rsplit(":", 1)[1] if ":" in s else s


def record_key(row):
    c = row.get("c_hgvs", "").strip()
    return variant_desc(c) if c else variant_desc(row.get("normalized", "").strip())


def split_body(s):
    return s.split(":", 1)[1].strip() if ":" in s else ""


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


def load_source_proteins(audit):
    proteins = defaultdict(list)
    seen = set()
    for name in sorted(f for f in os.listdir(PUB_DIR) if f.lower().endswith("pub.html")):
        gene = audit.gene_from_filename(name)
        records, _ = audit.read_records(os.path.join(PUB_DIR, name))
        for lines in records:
            acc, fields, _ = audit.parse_record(lines)
            if not acc or (gene, acc) in seen:
                continue
            seen.add((gene, acc))
            for _, txt in audit.split_alleles(fields):
                proteins[f"{gene}_{acc}"].extend(source_proteins(txt))
    return proteins


def pair_label(final_p, src_ps, audit):
    if not final_p.strip():
        return "uncheckable"
    if not src_ps:
        return "uncheckable"
    res = [audit.compare_protein(s, final_p) for s in src_ps]
    if any(st in audit.GOOD_P for st, _ in res):
        return "consistent"
    if any(st in audit.BAD_P for st, _ in res):
        return "inconsistent"
    return "uncheckable"


def c_region(c):
    c = variant_desc(c)
    if not c.startswith("c."):
        return "none"
    body = c[2:]
    if body.startswith("*"):
        return "utr3"
    if body.startswith("-") and not re.match(r"^-\d+[+-]", body):
        return "utr5"
    if re.search(r"\d[+-]\d", body):
        return "intronic"
    return "exonic"


def genomic_of(row):
    g = row.get("g_hgvs", "").strip()
    n = row.get("normalized", "").strip()
    return g if g else (n if n.startswith("NC_") else "")


def coding_token(source_change_id):
    parts = source_change_id.split("|", 4)
    return parts[4] if len(parts) == 5 else "-"


def choose(group, tally):
    labels = [r["label"] for r in group]
    if "consistent" in labels:
        consistent = [r for r in group if r["label"] == "consistent"]
        if len(consistent) == 1:
            return consistent[0], "protein"
        pool = consistent
    else:
        pool = [r for r in group if r["label"] != "inconsistent"] or group
        if len(pool) == 1 and len(group) > 1:
            return pool[0], "protein"
    region = [r for r in pool if r["c_region_match"]]
    if len(region) == 1:
        return region[0], "c_region"
    mane = [r for r in pool if r["mane_select"]]
    if len(mane) == 1:
        return mane[0], "mane_select"
    prior = []
    for r in pool:
        ok, bad = tally[(r["gene"], r["track"], "consistent")], tally[(r["gene"], r["track"], "inconsistent")]
        if r["conflict"]:
            ok -= r["label"] == "consistent"
            bad -= r["label"] == "inconsistent"
        prior.append((ok + 1) / (ok + bad + 2))
    if len(set(prior)) > 1:
        return pool[prior.index(max(prior))], "gene_prior"
    return sorted(pool, key=lambda r: (r["priority"], r["key"]))[0], "track_priority"


def main():
    audit = load_audit()
    dedup_fields, dedup_rows = read_tsv(DEDUP_TSV)
    _, augment_rows = read_tsv(AUGMENT_TSV)
    src_proteins = load_source_proteins(audit)

    sids = defaultdict(set)
    for r in augment_rows:
        sid = r.get("source_change_id", "").strip()
        if sid:
            sids[(r["gene"].strip().upper(), r["accession"].strip(), record_key(r))].add(sid)

    records = {}
    pairs = {}
    for i, row in enumerate(dedup_rows):
        gene = row["gene"].strip().upper()
        key = row["dedup_key"]
        records[i] = row
        for acc in [a.strip() for a in row["accession_list"].split(";") if a.strip()]:
            iid = f"{gene}_{acc}"
            pairs[(iid, key)] = {
                "record": i, "gene": gene, "accession": acc, "iid": iid, "key": key,
                "track": row.get("source_track", ""), "mane_select": row.get("mane_select", "") == "MANE Select",
                "genomic": genomic_of(row), "c_hgvs": row.get("c_hgvs", "").strip(),
                "label": pair_label(split_body(row.get("p_hgvs", "")), src_proteins.get(iid, []), audit),
                "priority": TRACK_PRIORITY.get(row.get("source_track", ""), 99),
                "sids": sids.get((gene, acc, key), set()),
            }

    by_change = defaultdict(list)
    for p in pairs.values():
        for sid in p["sids"]:
            by_change[(p["iid"], sid)].append(p)
    groups = {k: v for k, v in by_change.items() if len({p["key"] for p in v}) > 1}

    members = []
    for (iid, sid), group in groups.items():
        conflict = {"consistent", "inconsistent"} <= {p["label"] for p in group}
        src_region = c_region(coding_token(sid))
        for p in group:
            members.append(dict(p, sid=sid, conflict=conflict, c_region_match=c_region(p["key"]) == src_region))
    tally = Counter((m["gene"], m["track"], m["label"]) for m in members if m["conflict"])
    by_group = defaultdict(list)
    for m in members:
        by_group[(m["iid"], m["sid"])].append(m)

    kept_in, dropped_in = defaultdict(set), defaultdict(set)
    tier = {}
    alt = defaultdict(list)
    log_rows = []
    decided = Counter()
    for (iid, sid), group in sorted(by_group.items()):
        kept, signal = choose(group, tally)
        same_genomic = len({m["genomic"] for m in group}) == 1
        labels = {m["label"] for m in group}
        if not same_genomic and (labels == {"inconsistent"} or labels == {"uncheckable"}):
            group_tier = "Review"
        elif signal == "protein":
            group_tier = "A"
        else:
            group_tier = "B"
        decided[(signal, group_tier)] += 1
        dropped = [m for m in group if m["key"] != kept["key"]]
        kept_in[iid].add(kept["key"])
        for m in dropped:
            dropped_in[iid].add(m["key"])
        pk = (iid, kept["key"])
        if TIER_RANK[group_tier] > TIER_RANK[tier.get(pk, "A")]:
            tier[pk] = group_tier
        tier.setdefault(pk, group_tier)
        if same_genomic:
            for m in dropped:
                desc = m["c_hgvs"] or m["key"]
                if desc not in alt[pk]:
                    alt[pk].append(desc)
        log_rows.append({
            "source_change_id": sid, "individual_id": iid, "kept": kept["key"],
            "dropped": " | ".join(m["key"] for m in dropped), "deciding_signal": signal, "tier": group_tier,
            "group_kind": "same_genomic_change" if same_genomic else "different_position",
            "labels": " | ".join(f"{m['key']}:{m['track']}:{m['label']}" for m in group),
            "tracks": " | ".join(sorted({m["track"] for m in group})),
        })

    in_groups = defaultdict(set)
    for (iid, sid), group in by_group.items():
        for m in group:
            in_groups[(iid, m["key"])].add(sid)
    removed = set()
    for iid, keys in dropped_in.items():
        for key in keys - kept_in[iid]:
            p = pairs[(iid, key)]
            if p["sids"] and p["sids"] <= in_groups[(iid, key)]:
                removed.add((iid, key))

    out_rows = []
    pair_rows = []
    for i, row in records.items():
        gene = row["gene"].strip().upper()
        accs = [a.strip() for a in row["accession_list"].split(";") if a.strip()]
        keep = [a for a in accs if (f"{gene}_{a}", row["dedup_key"]) not in removed]
        if not keep:
            continue
        new = dict(row)
        new["accession_list"] = ";".join(keep)
        new["patient_count"] = str(len(keep))
        out_rows.append(new)
        for a in keep:
            pk = (f"{gene}_{a}", row["dedup_key"])
            pair_rows.append({
                "gene": row["gene"], "accession": a, "individual_id": pk[0], "dedup_key": row["dedup_key"],
                "confidence_tier": tier.get(pk, "A"), "alt_descriptions": " | ".join(alt.get(pk, [])),
                "source_change_ids": " ; ".join(sorted(pairs[pk]["sids"])),
            })

    write_tsv(OUT_DEDUP, dedup_fields, out_rows)
    write_tsv(OUT_PAIRS, PAIR_COLS, pair_rows)
    write_tsv(OUT_LOG, LOG_COLS, log_rows)

    tiers = Counter(r["confidence_tier"] for r in pair_rows)
    kinds = Counter(r["group_kind"] for r in log_rows)
    lines = [
        "resolve_tracks summary",
        "======================",
        f"  input records / pairs                : {len(dedup_rows)} / {len(pairs)}",
        f"  pairs without a source_change_id     : {sum(1 for p in pairs.values() if not p['sids'])}",
        f"  conflicting source changes           : {len(log_rows)}",
        f"    same genomic change                : {kinds['same_genomic_change']}",
        f"    different position                 : {kinds['different_position']}",
    ]
    for (signal, t), n in sorted(decided.items()):
        lines.append(f"    decided by {signal:<16} tier {t:<6}: {n}")
    lines += [
        f"  pairs removed                        : {len(removed)}",
        f"  output records / pairs               : {len(out_rows)} / {len(pair_rows)}",
        f"  pairs with alt_descriptions          : {sum(1 for r in pair_rows if r['alt_descriptions'])}",
        f"  confidence_tier A / B / Review       : {tiers['A']} / {tiers['B']} / {tiers['Review']}",
        f"  output                               : {OUT_DEDUP}",
    ]
    Path(OUT_STATS).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
