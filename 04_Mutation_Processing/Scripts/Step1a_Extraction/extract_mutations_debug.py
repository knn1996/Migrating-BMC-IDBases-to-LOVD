import re
from pathlib import Path

RE_STRIP = re.compile(r"<[^>]+>")
RE_LOC   = re.compile(r"/loc:", re.IGNORECASE)

IDBASE_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\02_Source_Database\idbase"
OUT_TXT    = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\debug_zero_rows.txt"

ZERO_ROW_GENES = {"CARD9", "DCLRE1C", "DNMT3B", "IGHM", "PIK3R1", "PTPN11",
                  "RASA1", "SH2D1A", "TAPBP", "TCIRG1", "WAS"}

CONTEXT = 5

with open(OUT_TXT, "w", encoding="utf-8") as out:
    for folder in sorted(Path(IDBASE_DIR).iterdir()):
        if not folder.is_dir() or not folder.name.lower().endswith("base"):
            continue
        gene = re.sub(r"base$", "", folder.name, flags=re.IGNORECASE).upper()
        if gene not in ZERO_ROW_GENES:
            continue

        pub = next((f for f in folder.iterdir() if f.name.lower().endswith("pub.html")), None)
        if not pub:
            out.write(f"\n{'='*60}\n{gene}: NO PUB.HTML\n")
            continue

        lines = RE_STRIP.sub("", pub.read_text(encoding="utf-8", errors="replace")).splitlines()
        out.write(f"\n{'='*60}\n{gene}: {len(lines)} lines\n{'='*60}\n")

        for i, line in enumerate(lines):
            s = line.strip()
            if (RE_LOC.search(s)
                    or "Systematic name" in s
                    or s.startswith("Accession")
                    or re.search(r"\bFeature\b", s)
                    or re.match(r"/name:", s)):
                start = max(0, i - CONTEXT)
                end   = min(len(lines), i + CONTEXT + 1)
                out.write(f"\n  --- context around line {i} ---\n")
                for j in range(start, end):
                    marker = ">>" if j == i else "  "
                    out.write(f"  {marker} {j:6d}  {lines[j].strip()[:120]}\n")

print(f"Written to {OUT_TXT}")