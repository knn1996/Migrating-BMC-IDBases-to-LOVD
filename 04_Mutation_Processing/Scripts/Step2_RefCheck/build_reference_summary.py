import csv
import glob
import os
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
_THESIS_DIR = (_SCRIPT_DIR / ".." / ".." / "..").resolve()
_PROC_DIR   = _THESIS_DIR / "04_Mutation_Processing"

SOURCE1_CHECK = _PROC_DIR / "Output" / "Step2_RefCheck" / "reference_check_source_1.tsv"
LRG_OFFSET    = _PROC_DIR / "Output" / "Step2_RefCheck" / "lrg_offset_results.csv"
LRG_NM_CSV    = _PROC_DIR / "Output" / "Step1b_RSG_Mapping" / "LRG_with_NM.csv"
BLAST_DIR     = _PROC_DIR / "Output" / "Step3_BLAST"
OUT_CSV       = _PROC_DIR / "Output" / "Step2_RefCheck" / "reference_summary.csv"

THRESHOLD = 0.90


def read_csv(path, delimiter=","):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def load_blast_first_hits(blast_dir):
    hits = {}
    for fpath in glob.glob(os.path.join(blast_dir, "*_blast.tsv")):
        gene = os.path.basename(fpath).replace("_blast.tsv", "").upper()
        with open(fpath, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 12:
                    hits[gene] = {"ref": parts[1].strip(), "sstart": parts[10].strip(), "send": parts[11].strip()}
                break
    return hits


def main():
    s1       = {r["gene"]: r for r in read_csv(SOURCE1_CHECK, delimiter="\t")}
    lrg_off  = {r["gene"]: r for r in read_csv(LRG_OFFSET)}
    lrg_nm   = {r["name"]: r for r in read_csv(LRG_NM_CSV)}
    blasts   = load_blast_first_hits(BLAST_DIR)

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

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "source", "ref", "sstart", "send", "match_pct", "non_matching"])
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Written {len(out_rows)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
