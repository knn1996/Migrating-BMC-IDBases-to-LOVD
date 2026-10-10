import os
import re
import html
from collections import defaultdict, Counter
from pathlib import Path
import pandas as pd

THESIS_DIR  = os.environ.get("THESIS_DIR", r"D:\Assignment\Thesis")
PUB_DIR     = os.environ.get("PUB_DIR", os.path.join(THESIS_DIR, "02_Source_Database", "idbase_pub_html_only"))
OUT_ROOT    = os.environ.get("OUT_ROOT", os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output"))
STEP8       = os.path.join(OUT_ROOT, "Step8_Rescue")
FINAL_TSV   = os.environ.get("FINAL_TSV", os.path.join(STEP8, "lovd_flat_with_patients_final.tsv"))
AUGMENT_TSV = os.environ.get("AUGMENT_TSV", os.path.join(STEP8, "augmented_merged_variants.tsv"))
UNRES_TSV   = os.environ.get("UNRES_TSV", os.path.join(STEP8, "unresolved_disposition.tsv"))
EXTRACT_TSV = os.environ.get("EXTRACT_TSV", os.path.join(OUT_ROOT, "Step1a_Extraction", "all_mutations.tsv"))
RESCUED_TSV = os.environ.get("RESCUED_TSV", os.path.join(STEP8, "resolved_unresolved.tsv"))
OUT_DIR     = os.environ.get("OUT_DIR", os.path.join(OUT_ROOT, "Step10_Audit"))
AUDIT_GENES = {g.strip().upper() for g in os.environ.get("AUDIT_GENES", "").split(",") if g.strip()}
OUT_OF_SCOPE = {"BTK", "SH2D1A"}
GOOD_P = {"MATCH", "MATCH_FS", "MATCH_INFRAME", "MATCH_START_LOSS"}
BAD_P = {"DIFFERENT", "POS_SHIFT", "FS_START_DIFF", "FS_STOP_DIFF", "TYPE_DIFF_SAME_POS"}

LABEL_W = 16

FIELD_MAP = [
    (("Sex",), "Individual/Gender"),
    (("Ethnic origin",), "Individual/Origin/Geographic"),
    (("Age",), "Individual/Age"),
    (("Relative",), "Individual/Remarks"),
    (("Disease", "Diagnosis", "Phenotype"), "Phenotype/Disease"),
    (("Symptoms",), "Phenotype/Additional"),
    (("Description",), "VariantOnGenome/Remarks"),
]
REF_MAP = [
    ("RefAuthors", "Reference/Authors"),
    ("RefTitle", "Reference/Title"),
    ("RefLoc", "Reference/Location"),
]

AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
       "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
       "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
       "X": "Ter", "*": "Ter", "U": "Sec"}
AA3_SET = set(AA3.values()) | {"Ter"}

RE_PUBMED = re.compile(r"PUBMED\s*[;:]\s*(\d+)", re.I)
RE_ALLELE = re.compile(r"Allele\s+(1\s+and\s+2|\d+)\s*:\s*", re.I)
RE_CPART  = re.compile(r"(c\.[^\s]+)$")


def norm(s):
    s = "" if s is None or (isinstance(s, float) and pd.isna(s)) else str(s)
    s = re.sub(r"\s+", " ", s).strip().strip(";,. ").casefold()
    return s


def norm_gender(s):
    u = norm(s).upper()
    if u in {"XY", "MALE", "M"}:
        return "m"
    if u in {"XX", "FEMALE", "F"}:
        return "f"
    return norm(s)


def gene_from_filename(name):
    g = re.sub(r"pub\.html$", "", name, flags=re.I)
    g = re.sub(r"base$", "", g, flags=re.I)
    return g.upper()


RE_GLUED_ID = re.compile(r"(?<=[^\n])(ID {14}\S)")


def read_records(path):
    raw = Path(path).read_text(encoding="utf-8", errors="replace")
    text = html.unescape(re.sub(r"<[^>]+>", "", raw))
    text, n_glued = RE_GLUED_ID.subn(r"\n//\n\1", text)
    records, cur, split_extra = [], [], 0
    for line in text.splitlines():
        if line.strip() == "//":
            if cur:
                records.append(cur)
            cur = []
            continue
        if not line.strip():
            continue
        label, value = line[:LABEL_W].strip(), line[LABEL_W:].rstrip()
        if label == "ID" and any(l == "Accession" for l, _ in cur):
            records.append(cur)
            cur = []
            split_extra += 1
        cur.append((label, value))
    if cur:
        records.append(cur)
    return records, n_glued + split_extra


def parse_record(lines):
    fields = defaultdict(list)
    refs = []
    for label, value in lines:
        if label == "RefNumber":
            refs.append(defaultdict(list))
            continue
        if label.startswith("Ref"):
            if not refs:
                refs.append(defaultdict(list))
            refs[-1][label].append(value.strip())
            continue
        fields[label].append(value)
    acc = " ".join(fields.get("Accession", [""])).strip().split(" ")[0]
    return acc, fields, refs


def joined(values):
    return re.sub(r"\s+", " ", " ".join(v.strip() for v in values)).strip()


def source_metadata(fields, refs):
    out = {}
    for labels, col in FIELD_MAP:
        val = ""
        for lab in labels:
            if fields.get(lab):
                val = joined(fields[lab])
                break
        out[col] = val
    pmids, primary = [], None
    for r in refs:
        m = RE_PUBMED.search(" ".join(r.get("RefCrossRef", [])))
        if m:
            pmids.append(m.group(1))
            if primary is None:
                primary = r
    if primary is None and refs:
        primary = refs[0]
    out["Reference/PubMed_ID"] = pmids[0] if pmids else ""
    out["_all_pmids"] = pmids
    for lab, col in REF_MAP:
        out[col] = joined(primary.get(lab, [])) if primary else ""
        out["_all_" + col] = joined(sum((r.get(lab, []) for r in refs), []))
    out["_n_refs"] = len(refs)
    return out


def split_alleles(fields):
    sysname = joined(fields.get("Systematic name", []))
    if not sysname:
        return []
    parts = RE_ALLELE.split(sysname)
    units = []
    if parts[0].strip():
        units.append(("0", parts[0].strip()))
    for i in range(1, len(parts) - 1, 2):
        lab = re.sub(r"\s+", " ", parts[i].lower())
        units.append(("1+2" if "and" in lab else lab, parts[i + 1].strip().rstrip(",;")))
    seen, dedup = {}, []
    for allele, txt in units:
        key = norm(txt)
        if key in seen:
            seen[key][0] = seen[key][0] + "+" + allele if allele not in seen[key][0] else seen[key][0]
            continue
        rec = [allele, txt]
        seen[key] = rec
        dedup.append(rec)
    return [(a, t) for a, t in dedup]


def levels(txt):
    out = defaultdict(list)
    for tok in re.split(r"[,;]\s*", txt):
        tok = tok.strip()
        if not tok:
            continue
        for p, lv in (("g.", "g"), ("c.", "c"), ("r.", "r"), ("p.", "p")):
            if tok.startswith(p):
                out[lv].append(tok)
                break
        else:
            if tok.upper().startswith("IVS"):
                out["ivs"].append(tok)
            else:
                out["other"].append(tok)
    return out


def aa3(tok):
    if tok in AA3:
        return AA3[tok]
    t = tok[:1].upper() + tok[1:].lower()
    if t in ("Stop", "Xaa"):
        return "Ter"
    return t if t in AA3_SET else None


def norm_protein(p):
    if not p or (isinstance(p, float) and pd.isna(p)):
        return None
    s = str(p)
    s = s.split(":p.")[-1] if ":p." in s else s
    s = re.sub(r"^p\.", "", s.strip()).strip("()")
    if s in ("", "?", "0", "=", "0?"):
        return {"kind": "special", "raw": s}
    m = re.match(r"^([A-Z][a-z]{2}|[A-Z*])(\d+)(.*)$", s)
    if not m or aa3(m.group(1)) is None:
        return {"kind": "unparsed", "raw": s}
    ref, pos, rest = aa3(m.group(1)), int(m.group(2)), m.group(3)
    fs = re.search(r"fs(?:\*|X|Ter)?(\d+)?", rest, re.I)
    if fs:
        return {"kind": "fs", "ref": ref, "pos": pos, "n": int(fs.group(1)) if fs.group(1) else None}
    m2 = re.match(r"^([A-Z][a-z]{2}|[A-Z*])$", rest)
    if m2 and aa3(m2.group(1)):
        return {"kind": "sub", "ref": ref, "pos": pos, "alt": aa3(m2.group(1))}
    if rest == "=":
        return {"kind": "sub", "ref": ref, "pos": pos, "alt": ref}
    op = re.search(r"delins|del|dup|ins|ext", rest)
    return {"kind": "other", "ref": ref, "pos": pos, "op": op.group(0) if op else rest.lower()}


def compare_protein(src, fin):
    a, b = norm_protein(src), norm_protein(fin)
    if a is None:
        return "SOURCE_NO_P", None
    if b is None:
        return "FINAL_NO_P", None
    if a["kind"] == "unparsed" or b["kind"] == "unparsed":
        return "NOT_COMPARABLE", None
    if b["kind"] == "special":
        if b["raw"] == "?" and a["kind"] == "sub" and a["pos"] == 1 and a["ref"] == "Met":
            return "MATCH_START_LOSS", 0
        if b["raw"] == "=" and a["kind"] == "sub":
            return ("MATCH", 0) if a["ref"] == a["alt"] else ("DIFFERENT", None)
        return "NOT_COMPARABLE", None
    if a["kind"] == "special":
        return "NOT_COMPARABLE", None
    if a["kind"] == "sub" and b["kind"] == "sub":
        if a["ref"] == b["ref"] and a["alt"] == b["alt"]:
            return ("MATCH", 0) if a["pos"] == b["pos"] else ("POS_SHIFT", b["pos"] - a["pos"])
        return "DIFFERENT", None
    if a["kind"] == "fs" and b["kind"] == "fs":
        stop_b = b["pos"] + b["n"] - 1 if b["n"] else None
        cands = set()
        if a["n"]:
            cands |= {a["n"], a["pos"] + a["n"] - 1, a["pos"] + a["n"]}
        stop_ok = stop_b is not None and any(abs(stop_b - c) <= 1 for c in cands)
        near = abs(a["pos"] - b["pos"]) <= 3
        if near and stop_ok:
            return "MATCH_FS", b["pos"] - a["pos"]
        if near:
            return "FS_STOP_DIFF", b["pos"] - a["pos"]
        if stop_ok:
            return "FS_START_DIFF", b["pos"] - a["pos"]
        return "DIFFERENT", b["pos"] - a["pos"]
    if a["kind"] == "other" and b["kind"] == "other" and a["op"] == b["op"]:
        if abs(a["pos"] - b["pos"]) <= 3:
            return "MATCH_INFRAME", b["pos"] - a["pos"]
        return "POS_SHIFT", b["pos"] - a["pos"]
    if {a["kind"], b["kind"]} <= {"fs", "other", "sub"} and abs(a["pos"] - b["pos"]) <= 1:
        return "TYPE_DIFF_SAME_POS", 0
    return "DIFFERENT", b["pos"] - a["pos"] if "pos" in a and "pos" in b else None


def cpart(s):
    if not isinstance(s, str):
        return ""
    m = RE_CPART.search(s.split(":")[-1] if ":" in s else s)
    return m.group(1) if m else s.split(":")[-1]


def load_tables():
    final = pd.read_csv(FINAL_TSV, sep="\t", dtype=str).fillna("")
    aug = pd.read_csv(AUGMENT_TSV, sep="\t", dtype=str).fillna("")
    unres = pd.read_csv(UNRES_TSV, sep="\t", dtype=str).fillna("")
    resc = pd.read_csv(RESCUED_TSV, sep="\t", dtype=str).fillna("") if os.path.exists(RESCUED_TSV) else pd.DataFrame(columns=["gene", "accession", "sysname"])
    ext = pd.read_csv(EXTRACT_TSV, sep="\t", dtype=str).fillna("")
    for df, g in ((final, "_gene"), (aug, "gene"), (unres, "gene"), (resc, "gene"), (ext, "gene")):
        df[g] = df[g].str.upper()
    return final, aug, unres, resc, ext


def index_by(df, gcol, acol):
    idx = defaultdict(list)
    for r in df.to_dict("records"):
        idx[(r[gcol], r[acol])].append(r)
    return idx


def audit_metadata(gene, acc, src, final_rows):
    out = []
    for col in [c for _, c in FIELD_MAP] + ["Reference/PubMed_ID"] + [c for _, c in REF_MAP]:
        fvals = sorted({r.get(col, "") for r in final_rows})
        if len(fvals) > 1:
            out.append((gene, acc, col, "INCONSISTENT_ACROSS_ROWS", src.get(col, ""), " || ".join(fvals)))
        fval = fvals[0] if fvals else ""
        sval = src.get(col, "")
        cmp = norm_gender if col == "Individual/Gender" else norm
        ns, nf = cmp(sval), cmp(fval)
        if not ns and not nf:
            continue
        if ns == nf:
            status = "MATCH"
        elif not nf:
            status = "MISSING_IN_FINAL"
        elif not ns:
            status = "EXTRA_IN_FINAL"
        elif "_all_" + col in src and src["_n_refs"] > 1 and nf == norm(src["_all_" + col]):
            status = "MERGED_ALL_REFS"
        elif all(it.strip() in ns for it in nf.split(";") if it.strip()):
            status = "SUBHEADERS_STRIPPED" if col == "Phenotype/Additional" else "PARTIAL_SUBSET"
        elif ns in nf:
            status = "PARTIAL_SUPERSET"
        else:
            status = "DIFFERENT"
        if status != "MATCH":
            out.append((gene, acc, col, status, sval, fval))
    extra = src["_all_pmids"][1:]
    if extra:
        out.append((gene, acc, "Reference/PubMed_ID(additional)", "EXTRA_REFS_NOT_CARRIED", ";".join(src["_all_pmids"]), src.get("Reference/PubMed_ID", "")))
    return out


def audit_unit(gene, acc, allele, txt, final_rows, aug_rows, unres_rows, resc_keys):
    lv = levels(txt)
    dna = set(lv["g"] + lv["c"] + lv["ivs"])
    src_p = lv["p"][0] if lv["p"] else ""
    row = {"gene": gene, "accession": acc, "allele": allele, "source_systematic_name": txt,
           "src_g": ";".join(lv["g"]), "src_c": ";".join(lv["c"]), "src_p": src_p,
           "fate": "", "fate_detail": "", "final_dedup_key": "", "final_p": "",
           "final_source_track": "", "protein_check": "", "protein_offset": ""}
    fin_keys = {r["_dedup_key"]: r for r in final_rows}
    fin_g = {r.get("VariantOnGenome/DNA", ""): r for r in final_rows if r.get("VariantOnGenome/DNA")}
    lineage = [r for r in aug_rows if r["sysname"] in dna]
    hit = None
    for r in lineage:
        k = cpart(r.get("c_hgvs", ""))
        if k in fin_keys:
            hit = fin_keys[k]
            break
        if r.get("g_hgvs") in fin_g:
            hit = fin_g[r["g_hgvs"]]
            break
    if hit is not None:
        row["fate"], row["fate_detail"] = "IN_FINAL", "traced_through_pipeline"
    else:
        best = None
        for r in final_rows:
            st, off = compare_protein(src_p, r.get("VariantOnTranscript/Protein", ""))
            if st in GOOD_P:
                best = r
                break
        if best is not None:
            hit = best
            row["fate"], row["fate_detail"] = "IN_FINAL", "protein_matched_untraced" if not lineage else "protein_matched_key_changed"
        elif lineage:
            row["fate"], row["fate_detail"] = "NOT_REGISTERED", "mutalyzer_ok_but_absent_from_final"
        else:
            unres = [r for r in unres_rows if r["sysname"] in dna]
            rescued = [s for s in dna if (gene, acc, s) in resc_keys]
            if unres and not rescued:
                row["fate"] = "KNOWN_UNRESOLVED"
                row["fate_detail"] = "|".join(sorted({f'{r.get("disposition","")}:{r.get("category","")}' for r in unres}))
            elif rescued:
                row["fate"], row["fate_detail"] = "NOT_REGISTERED", "rescued_but_absent_from_final"
            elif not dna:
                row["fate"], row["fate_detail"] = "NOT_REGISTERED", "source_has_no_dna_notation"
            elif any(re.search(r"ins$|del$|dup$", s) for s in dna) or any(re.search(r"_\d+ins$", s) for s in dna):
                row["fate"], row["fate_detail"] = "NOT_REGISTERED", "source_incomplete_notation"
            elif not final_rows:
                row["fate"], row["fate_detail"] = "NOT_REGISTERED", "accession_absent_from_final"
            else:
                row["fate"], row["fate_detail"] = "NOT_REGISTERED", "allele_absent_accession_present"
    if hit is not None:
        row["final_dedup_key"] = hit.get("_dedup_key", "")
        row["final_p"] = hit.get("VariantOnTranscript/Protein", "")
        row["final_source_track"] = hit.get("_source_track", "")
        st, off = compare_protein(src_p, row["final_p"])
        row["protein_check"] = st
        row["protein_offset"] = "" if off is None else off
    return row


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    final, aug, unres, resc, ext = load_tables()
    ext_tokens = defaultdict(set)
    for r in ext.to_dict("records"):
        ext_tokens[(r["gene"], r["accession"])].add(r["notation"])
    proc_acc = {(r["gene"], r["accession"]) for r in aug.to_dict("records")} | {(r["gene"], r["accession"]) for r in unres.to_dict("records")}
    f_idx = index_by(final, "_gene", "_accession")
    a_idx = index_by(aug, "gene", "accession")
    u_idx = index_by(unres, "gene", "accession")
    resc_keys = {(r["gene"], r["accession"], r["sysname"]) for r in resc.to_dict("records")}

    meta_rows, var_rows, gene_rows = [], [], []
    seen_pairs = set()
    repairs = {}
    files = sorted(f for f in os.listdir(PUB_DIR) if f.lower().endswith("pub.html"))
    for f in files:
        gene = gene_from_filename(f)
        if AUDIT_GENES and gene not in AUDIT_GENES:
            continue
        recs, n_repaired = read_records(os.path.join(PUB_DIR, f))
        repairs[gene] = n_repaired
        n_acc = 0
        for lines in recs:
            acc, fields, refs = parse_record(lines)
            if not acc:
                continue
            if (gene, acc) in seen_pairs:
                continue
            seen_pairs.add((gene, acc))
            n_acc += 1
            frows = f_idx.get((gene, acc), [])
            if frows:
                meta_rows += audit_metadata(gene, acc, source_metadata(fields, refs), frows)
            units = split_alleles(fields)
            n_units = len(units)
            for allele, txt in units:
                r = audit_unit(gene, acc, allele, txt, frows, a_idx.get((gene, acc), []), u_idx.get((gene, acc), []), resc_keys)
                r["n_source_alleles"] = n_units
                toks = set(filter(None, (r["src_g"] + ";" + r["src_c"] + ";" + r["src_p"]).split(";")))
                et = ext_tokens.get((gene, acc), set())
                r["stage_extracted"] = "all" if toks and toks <= et else ("partial" if toks & et else "none")
                r["stage_processed"] = "yes" if (gene, acc) in proc_acc else "no"
                if gene in OUT_OF_SCOPE:
                    r["fate"], r["fate_detail"] = "OUT_OF_SCOPE", "excluded_gene"
                var_rows.append(r)
        gene_rows.append((gene, f, n_acc))

    pub_genes = {g for g, _, _ in gene_rows}
    orphans = [(g, a, "", "ORPHAN_IN_FINAL", "", "accession not found in pub.html")
               for (g, a) in f_idx if (not AUDIT_GENES or g in AUDIT_GENES) and (g, a) not in seen_pairs]
    meta_rows += [(g, a, "Individual/ID", s, x, y) for g, a, _, s, x, y in orphans]

    meta = pd.DataFrame(meta_rows, columns=["gene", "accession", "field", "status", "source_value", "final_value"])
    var = pd.DataFrame(var_rows)
    meta.to_csv(os.path.join(OUT_DIR, "audit_metadata_flags.tsv"), sep="\t", index=False)
    var.to_csv(os.path.join(OUT_DIR, "audit_variants.tsv"), sep="\t", index=False)
    var[var["fate"] == "NOT_REGISTERED"].to_csv(os.path.join(OUT_DIR, "audit_not_registered.tsv"), sep="\t", index=False)
    var[var["protein_check"].isin(BAD_P)].to_csv(os.path.join(OUT_DIR, "audit_protein_flags.tsv"), sep="\t", index=False)

    bygene = var.groupby("gene").agg(
        source_units=("fate", "size"),
        in_final=("fate", lambda x: (x == "IN_FINAL").sum()),
        known_unresolved=("fate", lambda x: (x == "KNOWN_UNRESOLVED").sum()),
        not_registered=("fate", lambda x: (x == "NOT_REGISTERED").sum()),
        compound_het_units=("n_source_alleles", lambda x: (x > 1).sum()),
        protein_flags=("protein_check", lambda x: x.isin(BAD_P).sum()),
        protein_ok=("protein_check", lambda x: x.isin(GOOD_P).sum()),
    ).reset_index()
    if len(meta):
        mg = meta[~meta["status"].isin(["SUBHEADERS_STRIPPED", "EXTRA_REFS_NOT_CARRIED"])].groupby("gene").size().rename("metadata_flags")
        bygene = bygene.merge(mg, on="gene", how="left").fillna({"metadata_flags": 0})
    bygene["protein_flag_rate"] = (bygene["protein_flags"] / (bygene["protein_flags"] + bygene["protein_ok"]).where(lambda x: x > 0)).round(3)
    bygene.sort_values(["not_registered", "protein_flags"], ascending=False).to_csv(os.path.join(OUT_DIR, "audit_by_gene.tsv"), sep="\t", index=False)
    shifts = var[var["protein_check"] == "POS_SHIFT"].groupby("gene")["protein_offset"].agg(lambda s: Counter(s).most_common(3))
    summ = []
    summ.append(("pub_genes_parsed", len(pub_genes)))
    summ.append(("source_accessions", len(seen_pairs)))
    summ.append(("source_variant_units", len(var)))
    summ.append(("source_glued_records_repaired", sum(repairs.values())))
    for g, n in sorted(repairs.items()):
        if n:
            summ.append((f"glued_records_repaired:{g}", n))
    for k, v in var["fate"].value_counts().items():
        summ.append((f"fate:{k}", v))
    for k, v in var.loc[var["fate"] == "NOT_REGISTERED", "fate_detail"].value_counts().items():
        summ.append((f"not_registered:{k}", v))
    nrv = var[var["fate"] == "NOT_REGISTERED"]
    for k, v in nrv.groupby(["stage_extracted", "stage_processed"]).size().items():
        summ.append((f"not_registered_stage:extracted={k[0]}:processed={k[1]}", v))
    for k, v in var.groupby("stage_extracted").size().items():
        summ.append((f"all_units_extracted:{k}", v))
    for k, v in var.loc[var["fate"] == "IN_FINAL", "protein_check"].value_counts().items():
        summ.append((f"protein:{k}", v))
    for k, v in meta[meta["status"] != "EXTRA_REFS_NOT_CARRIED"].groupby(["field", "status"]).size().items():
        summ.append((f"meta:{k[0]}:{k[1]}", v))
    summ.append(("meta:accessions_with_extra_refs", int((meta["status"] == "EXTRA_REFS_NOT_CARRIED").sum())))
    for g, s in shifts.items():
        summ.append((f"pos_shift_offsets:{g}", str(s)))
    final_genes = set(final["_gene"]) if not AUDIT_GENES else set(final["_gene"]) & AUDIT_GENES
    for g in sorted(final_genes - pub_genes):
        summ.append(("final_gene_without_pub_html", g))
    sdf = pd.DataFrame(summ, columns=["metric", "value"])
    sdf.to_csv(os.path.join(OUT_DIR, "audit_summary.tsv"), sep="\t", index=False)
    print(sdf.to_string(index=False))
    print(f"\nOutputs in {OUT_DIR}")


if __name__ == "__main__":
    main()
