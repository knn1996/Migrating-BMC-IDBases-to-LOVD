import os
import re
import pandas as pd
from pathlib import Path

IDBASE_DIR = os.environ["IDBASE_DIR"]
IN_PATH    = os.environ["IN_PATH"]
OUT_PATH   = os.environ["OUT_PATH"]
LOG_PATH   = os.environ["LOG_PATH"]

RE_HTML       = re.compile(r"<[^>]+>")
RE_FIELD      = re.compile(r"^([A-Z][A-Za-z][A-Za-z /]*?)(\s{2,})(.+?)\s*$", re.MULTILINE)
RE_PUBMED     = re.compile(r"PUBMED\s*[;:]\s*(\d+)", re.IGNORECASE)
RE_ACCESSION  = re.compile(r"^Accession\s+(\S+)", re.MULTILINE)

LABEL_PRIORITY = [
    ("Sex",            "Individual/Gender"),
    ("Ethnic origin",  "Individual/Origin/Geographic"),
    ("Age",            "Individual/Age"),
    ("Relative",       "Individual/Remarks"),
    ("Disease",        "Phenotype/Disease"),
    ("Diagnosis",      "Phenotype/Disease"),
    ("Phenotype",      "Phenotype/Disease"),
    ("Symptoms",       "Phenotype/Additional"),
    ("Description",    "VariantOnGenome/Remarks"),
    ("RefAuthors",     "Reference/Authors"),
    ("RefTitle",       "Reference/Title"),
    ("RefLoc",         "Reference/Location"),
]

LOVD_COLS = sorted({lovd for _, lovd in LABEL_PRIORITY} | {"Reference/PubMed_ID"})


def normalize_gender(s):
    if not s:
        return ""
    upper = s.strip().upper()
    if upper in {"XY", "MALE", "M"}:
        return "M"
    if upper in {"XX", "FEMALE", "F"}:
        return "F"
    return s.strip()


def html_to_text(raw):
    text = raw.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">").replace("&nbsp;", " ")
    return RE_HTML.sub("", text)


def parse_record(record_text):
    out = {col: "" for col in LOVD_COLS}
    column_lines = {col: [] for col in LOVD_COLS}
    label_owner  = {}

    for m in RE_FIELD.finditer(record_text):
        label  = m.group(1).strip()
        indent = len(m.group(2))
        value  = m.group(3).strip()

        for src_label, lovd_col in LABEL_PRIORITY:
            if label == src_label:
                if lovd_col not in label_owner:
                    label_owner[lovd_col] = label
                if label_owner[lovd_col] == label:
                    column_lines[lovd_col].append((indent, value))
                break

        if label.startswith("RefCrossRef"):
            pm = RE_PUBMED.search(value)
            if pm and not out["Reference/PubMed_ID"]:
                out["Reference/PubMed_ID"] = pm.group(1)

    for col, lines in column_lines.items():
        if not lines:
            continue
        max_indent = max(ind for ind, _ in lines)
        min_indent = min(ind for ind, _ in lines)
        if max_indent == min_indent:
            out[col] = " ".join(v for _, v in lines)
        else:
            out[col] = " ".join(v for ind, v in lines if ind == max_indent)

    out["Individual/Gender"] = normalize_gender(out["Individual/Gender"])
    return out


def parse_pub_html(path):
    raw = Path(path).read_text(encoding="utf-8", errors="replace")
    text = html_to_text(raw)
    records = re.split(r"\n//\s*\n", text)
    out = {}
    for rec in records:
        m = RE_ACCESSION.search(rec)
        if not m:
            continue
        out[m.group(1).strip()] = parse_record(rec)
    return out


def build_metadata_cache():
    cache = {}
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
            print(f"  SKIP {gene}: no pub.html in {folder}")
            continue
        records = parse_pub_html(os.path.join(folder_path, pub))
        cache[gene] = records
        print(f"  {gene:<15} {len(records):>4} accessions parsed")
    return cache


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    df = pd.read_csv(IN_PATH, sep="\t", dtype=str).fillna("")

    df["_accession"] = df["_accession_list"].apply(
        lambda s: [a.strip() for a in s.split(";") if a.strip()] if s else [""]
    )
    exploded = df.explode("_accession", ignore_index=True)
    exploded["_accession"] = exploded["_accession"].fillna("")
    print(f"Exploded {len(df)} variants -> {len(exploded)} patient-variant rows\n")

    print("Parsing pub.html files...")
    cache = build_metadata_cache()

    for col in LOVD_COLS:
        exploded[col] = ""

    missing_log = []
    for i, row in exploded.iterrows():
        gene = row["_gene"].upper()
        acc  = row["_accession"]
        if not acc:
            continue
        gene_cache = cache.get(gene)
        if gene_cache is None:
            missing_log.append((gene, acc, "gene_not_in_cache"))
            continue
        meta = gene_cache.get(acc)
        if meta is None:
            missing_log.append((gene, acc, "accession_not_in_pub"))
            continue
        for col, value in meta.items():
            if value:
                exploded.at[i, col] = value

    exploded["Individual/ID"] = exploded["_gene"] + "_" + exploded["_accession"]

    leading = ["Individual/ID", "_gene", "_accession", "_dedup_key",
               "_source_track", "_patient_count", "_status", "_mane_select"]
    rest    = [c for c in exploded.columns
               if c not in leading and c not in ("_accession_list",)]
    exploded = exploded[leading + rest]

    exploded.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"\nWrote {len(exploded)} rows -> {OUT_PATH}")

    if missing_log:
        log_df = pd.DataFrame(missing_log, columns=["gene", "accession", "reason"])
        log_df.to_csv(LOG_PATH, sep="\t", index=False)
        print(f"Logged {len(missing_log)} missing accessions -> {LOG_PATH}")
        print(f"  by reason:\n{log_df['reason'].value_counts().to_string()}")
    else:
        print("All accessions matched to pub.html records.")

    print("\n=== Fill rates (patient/phenotype/reference cols) ===")
    for col in sorted(LOVD_COLS):
        n = (exploded[col] != "").sum()
        print(f"  {col:<35} {n:>6}/{len(exploded)}  ({100*n/len(exploded):>5.1f}%)")

    print("\n=== Individual/Gender distribution ===")
    print(exploded["Individual/Gender"].value_counts().head(10).to_string())


if __name__ == "__main__":
    main()
