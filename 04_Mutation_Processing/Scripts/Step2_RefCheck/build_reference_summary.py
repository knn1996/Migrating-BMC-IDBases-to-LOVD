import csv
import os

SOURCE1_CHECK = os.environ["SOURCE1_CHECK"]
LRG_OFFSET    = os.environ["LRG_OFFSET"]
LRG_NM_CSV    = os.environ["LRG_NM_CSV"]
BLAST_HITS    = os.environ["BLAST_HITS"]
OUT_CSV       = os.environ["OUT_CSV"]

THRESHOLD = 0.90


def read_csv(path, delimiter=","):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def load_blast_first_hits(hits_path):
    build_rank = {"hg16": 0, "hg17": 1, "hg18": 2}
    best = {}
    for row in read_csv(hits_path, delimiter="\t"):
        gene    = row["gene"].strip().upper()
        subject = row.get("subject_id", "").strip()
        if not subject or subject == "NO_HITS":
            continue
        try:
            pct = float(row["pct_identity"])
        except (ValueError, KeyError):
            continue
        rank = build_rank.get(row.get("build", "").strip().lower(), -1)
        cur  = best.get(gene)
        if cur is None or (pct, rank) > (cur["pct"], cur["rank"]):
            best[gene] = {"ref": subject, "sstart": row.get("s_start", "").strip(),
                          "send": row.get("s_end", "").strip(), "pct": pct, "rank": rank}
    return {g: {"ref": v["ref"], "sstart": v["sstart"], "send": v["send"]} for g, v in best.items()}


def main():
    s1       = {r["gene"]: r for r in read_csv(SOURCE1_CHECK, delimiter="\t")}
    lrg_off  = {r["gene"]: r for r in read_csv(LRG_OFFSET)}
    lrg_nm   = {r["name"]: r for r in read_csv(LRG_NM_CSV)}
    blasts   = load_blast_first_hits(BLAST_HITS)

    out_rows = []
    for gene in sorted(set(list(s1.keys()) + list(lrg_off.keys()))):
        r1  = s1.get(gene, {})
        off = lrg_off.get(gene, {})

        if float(r1.get("percentage_match", 0)) / 100 >= THRESHOLD:
            hit = blasts.get(gene) or blasts.get(gene + "_DNA")
            out_rows.append({
                "gene":         gene,
                "source":       "source_1",
                "ref":          hit["ref"] if hit else "not_available_yet",
                "sstart":       hit["sstart"] if hit else "",
                "send":         hit["send"] if hit else "",
                "match_pct":    r1.get("percentage_match", ""),
                "non_matching": r1.get("non_matching_accessions", ""),
            })

        if off.get("status") == "ok":
            pct = float(off["match_pct"]) * 100 if off.get("match_pct") else ""
            out_rows.append({
                "gene":         gene,
                "source":       "source_3",
                "ref":          lrg_nm.get(gene, {}).get("RSG", "") or "unknown",
                "sstart":       off.get("sstart", ""),
                "send":         off.get("send", ""),
                "match_pct":    str(round(pct, 2)) if pct != "" else "",
                "non_matching": off.get("non_matching_accessions", ""),
            })

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "source", "ref", "sstart", "send", "match_pct", "non_matching"])
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Written {len(out_rows)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
