import os
import re
from collections import Counter, defaultdict
from pathlib import Path

THESIS_DIR  = r"C:\Users\BornLoser\Desktop\Assignment\Thesis"
IDBASE_DIR  = os.path.join(THESIS_DIR, "02_Source_Database", "idbase")
LOG_DIR     = os.path.join(THESIS_DIR, "04_Mutation_Processing", "Logs")
OUT_LABELS  = os.path.join(LOG_DIR, "pub_html_label_survey.tsv")
OUT_PUBMED  = os.path.join(LOG_DIR, "pub_html_pubmed_samples.tsv")

RE_HTML  = re.compile(r"<[^>]+>")
RE_FIELD = re.compile(r"^([A-Z][A-Za-z][A-Za-z /]*?)\s{2,}(.+?)\s*$", re.MULTILINE)


def html_to_text(raw):
    text = raw.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">").replace("&nbsp;", " ")
    return RE_HTML.sub("", text)


def main():
    os.makedirs(LOG_DIR, exist_ok=True)

    label_genes = defaultdict(set)
    pubmed_lines = []

    folders = sorted(
        d for d in os.listdir(IDBASE_DIR)
        if os.path.isdir(os.path.join(IDBASE_DIR, d))
        and d.lower().endswith("base")
        and d.lower() != "immunomebase"
    )

    for folder in folders:
        gene = re.sub(r"base$", "", folder, flags=re.IGNORECASE).upper()
        folder_path = os.path.join(IDBASE_DIR, folder)
        pub = next((f for f in os.listdir(folder_path)
                    if f.lower().endswith("pub.html")), None)
        if not pub:
            continue
        text = html_to_text(Path(os.path.join(folder_path, pub)).read_text(encoding="utf-8", errors="replace"))

        for m in RE_FIELD.finditer(text):
            label = m.group(1).strip()
            label_genes[label].add(gene)

        for line in text.split("\n"):
            if "PUBMED" in line.upper() or "PMID" in line.upper():
                pubmed_lines.append((gene, line.strip()[:200]))

    with open(OUT_LABELS, "w", encoding="utf-8") as f:
        f.write("label\tgene_count\tgenes\n")
        for label, genes in sorted(label_genes.items(), key=lambda x: -len(x[1])):
            f.write(f"{label}\t{len(genes)}\t{';'.join(sorted(genes))}\n")
    print(f"Wrote {len(label_genes)} unique labels -> {OUT_LABELS}")

    with open(OUT_PUBMED, "w", encoding="utf-8") as f:
        f.write("gene\tline\n")
        for gene, line in pubmed_lines[:500]:
            f.write(f"{gene}\t{line}\n")
    print(f"Wrote {len(pubmed_lines)} PubMed-containing lines (first 500 -> file)")

    print("\n=== Top 30 most common labels ===")
    for label, genes in sorted(label_genes.items(), key=lambda x: -len(x[1]))[:30]:
        print(f"  {label:<35} {len(genes):>4} genes")

    print("\n=== Diagnosis-like labels (containing 'iagnos', 'henotyp', 'isease', 'linical') ===")
    for label, genes in sorted(label_genes.items(), key=lambda x: -len(x[1])):
        if re.search(r"iagnos|henotyp|isease|linical", label, re.IGNORECASE):
            print(f"  {label:<35} {len(genes):>4} genes  e.g. {sorted(genes)[:5]}")

    print("\n=== Sex-like labels ===")
    for label, genes in sorted(label_genes.items(), key=lambda x: -len(x[1])):
        if re.search(r"\bsex\b|gender|karyotype", label, re.IGNORECASE):
            print(f"  {label:<35} {len(genes):>4} genes  e.g. {sorted(genes)[:5]}")

    print("\n=== Origin-like labels ===")
    for label, genes in sorted(label_genes.items(), key=lambda x: -len(x[1])):
        if re.search(r"origin|ethnic|race|ancestry|population|country|geograph", label, re.IGNORECASE):
            print(f"  {label:<35} {len(genes):>4} genes  e.g. {sorted(genes)[:5]}")

    print("\n=== First 10 PubMed-containing lines ===")
    for gene, line in pubmed_lines[:10]:
        print(f"  {gene:<10} {line}")


if __name__ == "__main__":
    main()
