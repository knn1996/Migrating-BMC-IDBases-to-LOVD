import csv
import re
from collections import defaultdict

BASE = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\04_Mutation_Processing\Output\Step7_Mutalyzer"

FILES = {
    "NG_IDRefseq": BASE + r"\mutalyzer_results.tsv",
    "NM_IDRefseq": BASE + r"\mutalyzer_results_NM.tsv",
    "NM_MANE":     BASE + r"\mutalyzer_results_NM_MANE.tsv",
}
OUT = BASE + r"\mutalyzer_per_gene_summary.tsv"

RESCUABLE = {"ESEQUENCEMISMATCH", "EINTRONIC", "ENOSELECTORFOUND", "ESYNTAXUEOF"}

def extract_codes(s):
    return re.findall(r"'code':\s*'([^']+)'", s)

summary = {}
for track, path in FILES.items():
    summary[track] = defaultdict(lambda: {"total": 0, "ok": 0, "errors": defaultdict(int)})
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene = row["gene"].strip().upper()
            d = summary[track][gene]
            d["total"] += 1
            if row["status"].strip() == "ok":
                d["ok"] += 1
            else:
                codes = extract_codes(row.get("errors", "")) or ["UNKNOWN"]
                for c in codes:
                    d["errors"][c] += 1

all_genes = sorted(set(g for t in summary.values() for g in t))

fields = ["gene", "track", "total", "ok", "pct_ok",
          "ESEQUENCEMISMATCH", "EINTRONIC", "ENOSELECTORFOUND",
          "ESYNTAXUEOF", "OTHER", "rescuable", "non_rescuable"]

with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(fields)
    for gene in all_genes:
        for track in FILES:
            d = summary[track].get(gene)
            if not d or d["total"] == 0:
                continue
            tot   = d["total"]
            ok    = d["ok"]
            errs  = d["errors"]
            esm   = errs.get("ESEQUENCEMISMATCH", 0)
            ein   = errs.get("EINTRONIC", 0)
            enos  = errs.get("ENOSELECTORFOUND", 0)
            esyn  = errs.get("ESYNTAXUEOF", 0)
            other = sum(v for k, v in errs.items() if k not in RESCUABLE)
            resc  = esm + ein + enos + esyn
            w.writerow([gene, track, tot, ok, f"{ok/tot*100:.1f}",
                        esm, ein, enos, esyn, other, resc, other])

print(f"Written to {OUT}")
print(f"Genes: {len(all_genes)}, Tracks: {list(FILES.keys())}")

# Print overall summary by track
print("\n=== Overall by track ===")
print(f"{'Track':<15} {'Total':>7} {'OK':>7} {'%OK':>6} {'Rescuable':>10} {'Non-resc':>9}")
for track in FILES:
    rows = [summary[track][g] for g in all_genes if summary[track].get(g) and summary[track][g]["total"] > 0]
    tot  = sum(r["total"] for r in rows)
    ok   = sum(r["ok"] for r in rows)
    resc = sum(sum(v for k,v in r["errors"].items() if k in RESCUABLE) for r in rows)
    non  = sum(sum(v for k,v in r["errors"].items() if k not in RESCUABLE) for r in rows)
    print(f"{track:<15} {tot:>7} {ok:>7} {ok/tot*100:>5.1f}% {resc:>10} {non:>9}")
