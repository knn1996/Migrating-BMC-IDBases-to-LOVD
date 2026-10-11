import os
import re
import sys

import pandas as pd

INPUT_TSV  = os.environ["INPUT_TSV"]
OUTPUT_TSV = os.environ["OUTPUT_TSV"]
STATS_TXT  = os.environ["STATS_TXT"]
DEDUP_LOG  = os.environ.get("DEDUP_LOG", os.path.join(os.path.dirname(os.path.abspath(OUTPUT_TSV)), "dedup_log.tsv"))
TRACK_PRIORITY = {"NM_MANE": 0, "NG_IDRefseq": 1, "NM_IDRefseq": 2}
GENOMIC_HGVS = re.compile(r"^NC_\d+\.\d+:g\.\S+$")
TRANSCRIPT = re.compile(r"\(((?:NM|XM|NR|XR)_\d+\.\d+)\)")
LOG_COLS = ["gene", "rule", "kept_record", "kept_track", "dropped_record", "dropped_track",
            "source_units", "patients_moved"]

df = pd.read_csv(INPUT_TSV, sep="\t", low_memory=False)


def variant_desc(s):
    if pd.isna(s):
        return None
    s = str(s)
    if ":" in s:
        return s.rsplit(":", 1)[1]
    return s


def make_dedup_key(row):
    if pd.notna(row["c_hgvs"]):
        return variant_desc(row["c_hgvs"]), "c_hgvs"
    return variant_desc(row["normalized"]), "normalized"


def text(v):
    if pd.isna(v):
        return ""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v).strip()


def genomic_ids(row):
    ids = set()
    for col in ("g_hgvs", "nc_hgvs", "normalized"):
        v = text(row.get(col))
        if GENOMIC_HGVS.match(v):
            ids.add(v)
    return ids


def coordinate_ids(row):
    vals = [text(row.get(c)) for c in ("chrom", "pos_hg38", "ref", "alt")]
    return {tuple(vals)} if all(vals) else set()


def is_xm(row):
    return "XM_" in text(row.get("c_hgvs")) or "XM_" in text(row.get("nc_hgvs"))


def describe(row):
    m = TRANSCRIPT.search(text(row.get("c_hgvs")))
    tx = m.group(1) if m else ""
    return f"{text(row.get('dedup_key'))} [{tx}]" if tx else text(row.get("dedup_key"))


keys = df.apply(make_dedup_key, axis=1, result_type="expand")
df["_dedup_key"]      = keys[0]
df["_dedup_key_src"]  = keys[1]
df["_track_priority"] = df["source_track"].map(TRACK_PRIORITY).fillna(99)

input_units = {(text(r["gene"]), text(r["accession"]), text(r["allele_num"])) for _, r in df.iterrows()}

deduped_rows, units, gids, cids = [], [], [], []
for (gene, key), grp in df.groupby(["gene", "_dedup_key"], sort=False):
    best_idx = grp["_track_priority"].idxmin()
    rep = grp.loc[best_idx].copy()
    unique_accessions = grp["accession"].dropna().unique().tolist()
    rep["patient_count"]    = len(unique_accessions)
    rep["accession_list"]   = ";".join(sorted(unique_accessions))
    rep["dedup_key"]        = key
    rep["dedup_key_source"] = grp.loc[best_idx, "_dedup_key_src"]
    rep["tracks_merged"]    = ";".join(sorted(grp["source_track"].dropna().unique().tolist()))
    deduped_rows.append(rep)
    units.append({(text(r["accession"]), text(r["allele_num"])) for _, r in grp.iterrows()})
    g, c = set(), set()
    for _, r in grp.iterrows():
        g |= genomic_ids(r)
        c |= coordinate_ids(r)
    gids.append(g)
    cids.append(c)

result = pd.DataFrame(deduped_rows).reset_index(drop=True)
result = result.drop(columns=["_dedup_key", "_dedup_key_src", "_track_priority"], errors="ignore")
n_after_cnotation = len(result)

xm_flag = [is_xm(result.loc[i]) for i in range(len(result))]
priority = [TRACK_PRIORITY.get(text(result.loc[i, "source_track"]), 99) for i in range(len(result))]
dropped = set()
log_rows = []
nm_before = sum(1 for f in xm_flag if not f)


def accessions_of(i):
    return sorted({a for a, _ in units[i] if a})


def tracks_union(i, j):
    t = set(text(result.loc[i, "tracks_merged"]).split(";")) | set(text(result.loc[j, "tracks_merged"]).split(";"))
    t.discard("")
    return ";".join(sorted(t))


def refresh(i):
    accs = accessions_of(i)
    result.loc[i, "accession_list"] = ";".join(accs)
    result.loc[i, "patient_count"] = len(accs)


def unit_text(us):
    return ";".join(f"{a}:{al or '?'}" for a, al in sorted(us))


def log(gene, rule, kept, dropped_i, us, moved):
    log_rows.append({
        "gene": gene, "rule": rule,
        "kept_record": describe(result.loc[kept]) if kept is not None else "",
        "kept_track": text(result.loc[kept, "source_track"]) if kept is not None else "",
        "dropped_record": describe(result.loc[dropped_i]) if dropped_i is not None else "",
        "dropped_track": text(result.loc[dropped_i, "source_track"]) if dropped_i is not None else "",
        "source_units": unit_text(us), "patients_moved": moved,
    })


for gene, idx in result.groupby("gene", sort=False).groups.items():
    idx = list(idx)
    nm = [i for i in idx if not xm_flag[i]]
    xm = [i for i in idx if xm_flag[i]]
    if not nm or not xm:
        continue
    for x in xm:
        cands = [n for n in nm if (gids[x] & gids[n]) or (cids[x] & cids[n])]
        if not cands:
            continue
        n = min(cands, key=lambda k: (priority[k], -len(accessions_of(k)), k))
        before = set(accessions_of(n))
        units[n] |= units[x]
        result.loc[n, "tracks_merged"] = tracks_union(n, x)
        refresh(n)
        log(gene, "identical_change", n, x, units[x], len(set(accessions_of(n)) - before))
        dropped.add(x)
    for x in xm:
        if x in dropped:
            continue
        moves = {}
        for unit in sorted(units[x]):
            acc, allele = unit
            if not allele:
                continue
            hits = [n for n in nm if unit in units[n]]
            if len(hits) > 1:
                log(gene, "skip_in_cis_multi_record", None, None, {unit}, 0)
                continue
            if len(hits) != 1:
                continue
            n = hits[0]
            if gids[x] and gids[n] and not (gids[x] & gids[n]):
                log(gene, "skip_different_change_same_unit", n, None, {unit}, 0)
                continue
            moves.setdefault(n, set()).add(unit)
        for n, us in moves.items():
            before = set(accessions_of(n))
            units[n] |= us
            units[x] -= us
            result.loc[n, "tracks_merged"] = tracks_union(n, x)
            refresh(n)
            remaining = bool(units[x])
            log(gene, "same_source_unit_partial" if remaining else "same_source_unit", n,
                None if remaining else x, us, len(set(accessions_of(n)) - before))
        if moves:
            if units[x]:
                refresh(x)
            else:
                dropped.add(x)

kept = [i for i in range(len(result)) if i not in dropped]
nm_after = sum(1 for i in kept if not xm_flag[i])
output_units = set()
for i in kept:
    g = text(result.loc[i, "gene"])
    output_units |= {(g, a, al) for a, al in units[i]}

errors = []
if output_units != input_units:
    lost = sorted(input_units - output_units)[:5]
    added = sorted(output_units - input_units)[:5]
    errors.append(f"source units changed: lost {lost} added {added}")
if nm_after != nm_before:
    errors.append(f"NM records changed: {nm_before} before stage 2, {nm_after} after")
for i in kept:
    labels = {al for _, al in units[i]}
    if not labels <= {text(v) for v in df["allele_num"]}:
        errors.append(f"unknown allele label in record {describe(result.loc[i])}")
        break
if errors:
    sys.exit("DEDUP GUARD FAIL: " + " | ".join(errors))

n_xm_dropped = len(dropped)
result = result.loc[kept]

front_cols = ["gene", "dedup_key", "dedup_key_source", "patient_count", "accession_list", "tracks_merged"]
other_cols = [c for c in result.columns if c not in front_cols]
result = result[front_cols + other_cols]
result.to_csv(OUTPUT_TSV, sep="\t", index=False)
pd.DataFrame(log_rows, columns=LOG_COLS).to_csv(DEDUP_LOG, sep="\t", index=False)

rule_counts = pd.Series([r["rule"] for r in log_rows], dtype=object).value_counts()
patients_linked = df[["gene", "accession"]].dropna().drop_duplicates().shape[0]

lines = []
lines.append(f"Input rows                  : {len(df)}")
lines.append(f"After c.-notation dedup     : {n_after_cnotation}")
lines.append(f"XM/NM duplicate rows removed : {n_xm_dropped}")
lines.append(f"Unique variants (final)     : {len(result)}")
lines.append(f"Patients linked             : {patients_linked}")
lines.append(f"Patient-variant pairs       : {int(result['patient_count'].sum())}")
lines.append("")
lines.append("Stage-2 dedup log by rule:")
lines.append(rule_counts.to_string() if len(rule_counts) else "(none)")
lines.append("")
lines.append("Dedup key source breakdown:")
lines.append(result["dedup_key_source"].value_counts().to_string())
lines.append("")
lines.append("Representative track chosen:")
lines.append(result["source_track"].value_counts().to_string())
lines.append("")
lines.append("Track coverage (variants resolved by N distinct tracks):")
_ntracks = result["tracks_merged"].fillna("").apply(lambda s: len([t for t in str(s).split(";") if t]))
lines.append(_ntracks.value_counts().sort_index().to_string())
lines.append("")
lines.append("Patient-sharing distribution (patients per unique variant):")
lines.append(result["patient_count"].value_counts().sort_index().to_string())

report = "\n".join(lines)
print(report)
with open(STATS_TXT, "w") as f:
    f.write(report + "\n")
