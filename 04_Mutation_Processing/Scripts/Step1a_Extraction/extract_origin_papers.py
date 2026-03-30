import os
import re
from pathlib import Path
from collections import defaultdict

_SCRIPT_DIR = Path(__file__).parent
_THESIS_DIR = (_SCRIPT_DIR / ".." / ".." / "..").resolve()

IDBASE_DIR  = _THESIS_DIR / "02_Source_Database" / "idbase"
OUT_DIR     = _THESIS_DIR / "05_Testing" / "reference_check"
OUT_TSV     = OUT_DIR / "origin_paper.tsv"

RE_STRIP = re.compile(r"<[^>]+>")
RE_HREF  = re.compile(r'href="([^"]+)"', re.IGNORECASE)


def strip_tags(text):
    return RE_STRIP.sub("", text).strip()


def parse_pub_html(path):
    raw   = Path(path).read_text(encoding="utf-8", errors="replace")
    lines = raw.splitlines()

    entries = []
    acc = title_parts = link = None
    in_entry = False

    for line in lines:
        s = strip_tags(line).strip()

        if s.startswith("Accession"):
            m = re.match(r"Accession\s+(\S+)", s)
            if m:
                if acc and title_parts is not None:
                    entries.append((acc, " ".join(title_parts).strip(), link or ""))
                acc, title_parts, link, in_entry = m.group(1), [], "", True

        elif in_entry and s.startswith("RefCrossRef"):
            urls = RE_HREF.findall(line)
            for u in urls:
                if any(k in u.lower() for k in ("ncbi", "pubmed", "doi")):
                    link = u
                    break
            if not link and urls:
                link = urls[0]

        elif in_entry and s.startswith("RefTitle"):
            m = re.match(r"RefTitle\s+(.*)", s)
            if m:
                title_parts.append(m.group(1).strip())

        elif in_entry and s == "//":
            if acc and title_parts is not None:
                entries.append((acc, " ".join(title_parts).strip(), link or ""))
            acc = title_parts = link = None
            in_entry = False

    if acc and title_parts is not None:
        entries.append((acc, " ".join(title_parts).strip(), link or ""))

    return entries


def deduplicate(gene, entries):
    title_to_link = {title: link for _, title, link in entries if link}

    # Group accessions by (title, link), resolving missing links via fallback
    groups = defaultdict(list)
    for acc, title, link in entries:
        resolved_link = link or title_to_link.get(title, "None") or "None"
        groups[(title, resolved_link)].append(acc)

    return [
        {"idbase": gene, "accessions": ",".join(accs), "title": title, "link": link}
        for (title, link), accs in groups.items()
    ]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    folders = sorted(
        d for d in os.listdir(IDBASE_DIR)
        if os.path.isdir(os.path.join(IDBASE_DIR, d))
        and d.lower().endswith("base")
        and d.lower() != "immunomebase"
    )

    rows = []
    for folder in folders:
        gene        = re.sub(r"base$", "", folder, flags=re.IGNORECASE)
        folder_path = os.path.join(IDBASE_DIR, folder)
        pub         = next((f for f in os.listdir(folder_path) if f.lower().endswith("pub.html")), None)

        if not pub:
            print(f"  SKIP {gene}: no pub.html")
            continue

        entries = parse_pub_html(os.path.join(folder_path, pub))
        deduped = deduplicate(gene, entries)
        rows.extend(deduped)
        print(f"  {gene:<15}  {len(entries):4d} mutations  →  {len(deduped):3d} unique articles")

    with open(OUT_TSV, "w", encoding="utf-8") as fh:
        fh.write("idbase\taccessions\tarticle_title\tlink\n")
        for r in rows:
            fh.write(f"{r['idbase']}\t{r['accessions']}\t{r['title']}\t{r['link']}\n")

    print(f"\nTotal rows : {len(rows)}")
    print(f"Output     : {OUT_TSV}")


if __name__ == "__main__":
    main()
