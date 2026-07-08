import pandas as pd
import os
import sys
import re

INPUT_TSV  = os.environ["INPUT_TSV"]
OUTPUT_TSV = os.environ["OUTPUT_TSV"]
STATS_TXT  = os.environ["STATS_TXT"]
TRACK_PRIORITY = {"NM_MANE": 0, "NG_IDRefseq": 1, "NM_IDRefseq": 2}

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

keys = df.apply(make_dedup_key, axis=1, result_type="expand")
df["_dedup_key"]      = keys[0]
df["_dedup_key_src"]  = keys[1]
df["_track_priority"] = df["source_track"].map(TRACK_PRIORITY).fillna(99)

deduped_rows = []
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

result = pd.DataFrame(deduped_rows)
result = result.drop(columns=["_dedup_key", "_dedup_key_src", "_track_priority"], errors="ignore")
n_after_cnotation = len(result)

def genomic_pos(row):
    for col in ("g_hgvs", "nc_hgvs"):
        v = row.get(col)
        if pd.notna(v):
            m = re.search(r":g\.(\d+)", str(v))
            if m:
                return m.group(1)
    if pd.notna(row.get("pos_hg38")):
        return str(row["pos_hg38"]).split(".")[0]
    return None

result["_gpos"] = result.apply(genomic_pos, axis=1)
result["_xm"]   = (result["c_hgvs"].fillna("").str.contains("XM_") |
                   result["nc_hgvs"].fillna("").str.contains("XM_"))
result["_gid"]  = (result["gene"].astype(str) + ":" +
                   result["chrom"].astype(str) + ":" +
                   result["_gpos"].astype(str))

drop_idx = []
for gid, grp in result[result["_gpos"].notna()].groupby("_gid"):
    nm = grp[~grp["_xm"]]
    xm = grp[grp["_xm"]]
    if len(nm) > 0 and len(xm) > 0:
        best = nm["patient_count"].idxmax()
        nm_acc = set(a for a in str(result.loc[best, "accession_list"]).split(";") if a)
        xm_acc = set()
        for al in xm["accession_list"]:
            xm_acc.update(a for a in str(al).split(";") if a)
        merged = sorted(nm_acc | xm_acc)
        result.loc[best, "accession_list"] = ";".join(merged)
        result.loc[best, "patient_count"]  = len(merged)
        tracks = set(str(result.loc[best, "tracks_merged"]).split(";"))
        for tm in xm["tracks_merged"]:
            tracks.update(str(tm).split(";"))
        tracks.discard("")
        result.loc[best, "tracks_merged"] = ";".join(sorted(tracks))
        drop_idx.extend(xm.index.tolist())

n_xm_dropped = len(drop_idx)
result = result.drop(index=drop_idx).drop(columns=["_gpos", "_xm", "_gid"])

front_cols = ["gene", "dedup_key", "dedup_key_source", "patient_count", "accession_list", "tracks_merged"]
other_cols = [c for c in result.columns if c not in front_cols]
result = result[front_cols + other_cols]
result.to_csv(OUTPUT_TSV, sep="\t", index=False)

lines = []
lines.append(f"Input rows                  : {len(df)}")
lines.append(f"After c.-notation dedup     : {n_after_cnotation}")
lines.append(f"XM/NM duplicate rows removed : {n_xm_dropped}")
lines.append(f"Unique variants (final)     : {len(result)}")
lines.append(f"Patients linked             : {df['accession'].nunique()}")
lines.append(f"Patient-variant pairs       : {int(result['patient_count'].sum())}")
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
