import os
import re
import sys
import html
import logging
from pathlib import Path
from datetime import datetime

IDBASE_DIR = os.environ["IDBASE_DIR"]
OUT_PATH   = os.environ["OUT_PATH"]
LOG_PATH   = os.environ["LOG_PATH"]

RE_STRIP     = re.compile(r"<[^>]+>")
RE_LOC       = re.compile(r"/loc:.*?:\s*(\d+)(?:\.\.(\d+))?", re.IGNORECASE)
RE_ALLELE    = re.compile(r"Allele\s+(1\s+and\s+2|\d+)\s*:\s*", re.IGNORECASE)
RE_FEAT_HDR  = re.compile(r"^FeatureHeader\s+allele;\s*(\d+)", re.IGNORECASE)
RE_FEAT_DNA  = re.compile(r"^Feature\s+dna;", re.IGNORECASE)
RE_FEAT_OTHR = re.compile(r"^Feature\s+(rna|aa);", re.IGNORECASE)
RE_FEAT_CHG  = re.compile(r"^Feature\s+/change:\s*(.*)$", re.IGNORECASE)
RE_FEAT_REG  = re.compile(r"^Feature\s+/genomic_region:\s*(.*)$", re.IGNORECASE)
RE_IVS_TOKEN = re.compile(r"^g\.IVS(\d+)\s*([+-])\s*(\d+)\s*(.*)$", re.IGNORECASE)

VARIANT_PREFIXES = [
    ("g.", "genomic"),
    ("c.", "coding"),
    ("r.", "rna"),
    ("p.", "protein"),
    ("m.", "mitochondrial"),
    ("n.", "non_coding"),
]

REQUIRED_LEVELS = {"genomic", "coding", "protein"}


def setup_logging():
    os.makedirs(os.path.dirname(LOG_PATH), exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(LOG_PATH, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger(__name__)


def detect_variant_type(notation):
    for prefix, vtype in VARIANT_PREFIXES:
        if notation.startswith(prefix):
            return vtype
    return "unknown"


def classify_mutation(notation):
    n = notation.lower()
    if "delins" in n:
        return "indel"
    if "dup" in n:
        return "duplication"
    if "ins" in n:
        return "insertion"
    if "del" in n:
        return "deletion"
    if ">" in n:
        return "substitution"
    return "unknown"


def extract_pos_from_notation(notation):
    if notation.startswith("p."):
        m = re.search(r"[A-Za-z*]+(\d+)", notation[2:])
        if m:
            pos = int(m.group(1))
            m2  = re.search(r"_[A-Za-z*]*(\d+)", notation[2:])
            return pos, int(m2.group(1)) if m2 else pos
        return None, None
    m = re.search(r"[gcrnm]\.[-*]?(\d+)", notation)
    if m:
        pos = int(m.group(1))
        m2  = re.search(r"_(\d+)", notation)
        return pos, int(m2.group(1)) if m2 else pos
    return None, None


def ivs_change(rest):
    rest = rest.upper()
    s = re.match(r"^([ACGT]+)>([ACGT]*)$", rest)
    if s and s.group(2):
        return f"{s.group(1)}>{s.group(2)}"
    if s:
        return f"-{s.group(1)}"
    d = re.match(r"^DEL([ACGT]*)", rest)
    if d:
        return f"-{d.group(1)}"
    i = re.match(r"^INS([ACGT]+)", rest)
    if i:
        return f"+{i.group(1)}"
    return rest


def feature_change(change):
    c = change.upper().replace(" ", "")
    s = re.match(r"^([ACGT]+)->([ACGT]+)$", c)
    return f"{s.group(1)}>{s.group(2)}" if s else c


def own_feature_starts(notation, allele_num, features):
    m = RE_IVS_TOKEN.match(notation)
    if not m:
        return None
    intron, ch = int(m.group(1)), ivs_change(m.group(4))
    located = [f for f in features if f["loc"] is not None]
    pool = [f for f in located if f["allele"] == allele_num] or located
    cand = [f for f in pool if re.search(rf"intron;\s*{intron}\b", f["region"])]
    if not cand:
        cand = [f for f in pool if feature_change(f["change"]) == ch]
    if len(cand) > 1:
        narrowed = [f for f in cand if feature_change(f["change"]) == ch
                    or (ch and feature_change(f["change"]).endswith(ch.lstrip("+-")))]
        cand = narrowed or cand
    return {f["loc"][0] for f in cand}


def position_is_own(notation, allele_num, loc, from_own_allele, loc_kind, features):
    if loc[0] is None:
        return True
    starts = own_feature_starts(notation, allele_num, features)
    if starts is None:
        return from_own_allele or loc_kind == "dna"
    return len(starts) == 1 and loc[0] in starts


RE_PREF_BRACKET = re.compile(r"([gcrpnm])\.\s*\[([^\]]*)\]")
RE_BARE_BRACKET = re.compile(r"\[([^\]]*)\]")


def join_sysname_lines(lines):
    out, depth = "", 0
    for ln in lines:
        ln = ln.strip()
        if not out:
            out = ln
        elif depth > 0:
            out = out + " " + ln
        else:
            out = out.rstrip(" ,;") + ", " + ln
        depth += ln.count("[") - ln.count("]")
    return out


def expand_brackets(text):
    def pref(m):
        p = m.group(1)
        items = [re.sub(r"^[gcrpnm]\.\s*", "", x.strip()) for x in re.split(r"[;+]", m.group(2)) if x.strip()]
        return "; ".join(f"{p}.{x}" for x in items)
    text = RE_PREF_BRACKET.sub(pref, text)
    return RE_BARE_BRACKET.sub(lambda m: "; ".join(x.strip() for x in re.split(r"[;+]", m.group(1)) if x.strip()), text)


def split_alleles(sysname):
    parts = RE_ALLELE.split(expand_brackets(sysname))
    units = []
    if parts[0].strip(" ,;"):
        units.append(("0", parts[0]))
    for i in range(1, len(parts) - 1, 2):
        label = re.sub(r"\s+", " ", parts[i].strip().lower())
        units.append(("both" if "and" in label else label, parts[i + 1]))
    return units


def split_by_level(text):
    levels = {}
    for seg in text.split(","):
        tokens = [t.strip() for t in seg.split(";") if t.strip()]
        for tok in tokens:
            vtype = detect_variant_type(tok)
            if vtype != "unknown":
                levels.setdefault(vtype, []).append(tok)
    return levels


def compute_missing_flag(levels):
    missing = REQUIRED_LEVELS - set(levels.keys())
    return ",".join(sorted(missing)) if missing else ""


def expand_rows(gene, accession, sysname, loc_by_allele, first_loc, first_loc_kind, features, stats):
    units = split_alleles(sysname)
    distinct = {re.sub(r"\s+", " ", t.strip(" ,;")) for a, t in units if a != "both"}
    if len(distinct) > 1:
        stats["compound_het"] += 1

    for allele_mode, text in units:
        clean  = re.sub(r"\s+", " ", text).strip(" ,;")
        levels = split_by_level(clean)
        allele_nums = [1, 2] if allele_mode == "both" else [int(allele_mode)]

        if not levels:
            for an in allele_nums:
                loc = loc_by_allele.get(an, first_loc)
                ps, pe = extract_pos_from_notation(clean)
                yield {
                    "gene": gene, "accession": accession, "allele_num": an,
                    "variant_type": "unknown", "notation": clean,
                    "pos_start": ps if ps is not None else loc[0],
                    "pos_end":   pe if pe is not None else loc[1],
                    "mut_class": classify_mutation(clean),
                    "missing_flag": "genomic,coding,protein",
                }
            continue

        missing_flag = compute_missing_flag(levels)
        genomic = levels.get("genomic", [])
        if any(len(v) > 1 for v in levels.values()):
            stats["cis_multi"] += 1

        for an in allele_nums:
            from_own_allele = an in loc_by_allele
            loc = loc_by_allele.get(an, first_loc)
            loc_kind = "dna" if from_own_allele else first_loc_kind
            for vtype, tokens in levels.items():
                for idx, notation in enumerate(tokens):
                    ps, pe = extract_pos_from_notation(notation)
                    if vtype == "genomic" and ps is None:
                        if position_is_own(notation, an, loc, from_own_allele, loc_kind, features):
                            ps, pe = loc
                        else:
                            stats["position_not_own"] += 1
                    if idx < len(genomic):
                        mut_class = classify_mutation(genomic[idx])
                    elif genomic:
                        mut_class = classify_mutation(genomic[0])
                    else:
                        mut_class = classify_mutation(notation)
                    yield {
                        "gene": gene, "accession": accession, "allele_num": an,
                        "variant_type": vtype, "notation": notation,
                        "pos_start": ps,
                        "pos_end":   pe,
                        "mut_class": mut_class,
                        "missing_flag": missing_flag,
                    }


def parse_pub_html(pub_html_path, gene, log, stats):
    text = html.unescape(RE_STRIP.sub("", Path(pub_html_path).read_text(encoding="utf-8", errors="replace")))
    lines = text.splitlines()

    state = {}
    emitted = set()
    total   = [0]

    def reset(acc=""):
        state.update({"acc": acc, "sys": [], "loc_by_allele": {}, "first_loc": (None, None),
                      "first_loc_kind": "none", "cur_allele": None, "in_dna": False, "kind": None,
                      "features": []})

    def flush():
        acc = state["acc"]
        if acc and state["sys"] and acc not in emitted:
            emitted.add(acc)
            total[0] += 1
            sysname = join_sysname_lines(state["sys"])
            yield from expand_rows(gene, acc, sysname, state["loc_by_allele"], state["first_loc"],
                                   state["first_loc_kind"], state["features"], stats)

    reset()
    for line in lines:
        s = line.strip()

        if s.startswith("Accession"):
            m = re.match(r"Accession\s+(\S+)", s)
            if m:
                if state["acc"] and state["acc"] not in emitted:
                    stats["glued"] += 1
                    log.warning(f"  {gene}: record separator missing before {m.group(1)}; emitting {state['acc']} separately")
                yield from flush()
                reset(m.group(1))

        elif s.startswith("Systematic name"):
            m = re.match(r"Systematic name\s+(.+)", s)
            if m:
                state["sys"].append(m.group(1).strip())

        elif RE_FEAT_HDR.match(s):
            state["cur_allele"] = int(RE_FEAT_HDR.match(s).group(1))
            state["in_dna"] = False
            state["kind"] = None

        elif RE_FEAT_DNA.match(s):
            state["in_dna"] = True
            state["kind"] = "dna"
            state["features"].append({"allele": state["cur_allele"], "loc": None, "change": "", "region": ""})

        elif RE_FEAT_OTHR.match(s):
            state["in_dna"] = False
            state["kind"] = RE_FEAT_OTHR.match(s).group(1).lower()

        elif state["in_dna"] and state["features"] and RE_FEAT_CHG.match(s):
            f = state["features"][-1]
            f["change"] = (f["change"] + " " + RE_FEAT_CHG.match(s).group(1).strip()).strip()

        elif state["in_dna"] and state["features"] and RE_FEAT_REG.match(s):
            state["features"][-1]["region"] = RE_FEAT_REG.match(s).group(1).strip()

        elif RE_LOC.search(s):
            m = RE_LOC.search(s)
            loc = (int(m.group(1)), int(m.group(2)) if m.group(2) else int(m.group(1)))
            if state["first_loc"] == (None, None):
                state["first_loc"] = loc
                state["first_loc_kind"] = state["kind"] or "none"
            if state["in_dna"] and state["features"] and state["features"][-1]["loc"] is None:
                state["features"][-1]["loc"] = loc
            if state["in_dna"] and state["cur_allele"] is not None:
                state["loc_by_allele"].setdefault(state["cur_allele"], loc)

        elif s == "//":
            yield from flush()
            reset()

    yield from flush()

    if total[0] == 0:
        log.warning(f"  {gene}: 0 entries emitted — check pub.html structure")


def main():
    log   = setup_logging()
    start = datetime.now()

    log.info("=" * 60)
    log.info("extract_mutations.py — starting")
    log.info(f"IDBASE_DIR : {IDBASE_DIR}")
    log.info(f"OUT_PATH   : {OUT_PATH}")
    log.info(f"Started at : {start.strftime('%Y-%m-%d %H:%M:%S')}")
    log.info("=" * 60)

    folders = sorted(
        d for d in os.listdir(IDBASE_DIR)
        if os.path.isdir(os.path.join(IDBASE_DIR, d))
        and d.lower().endswith("base")
        and d.lower() != "immunomebase"
    )

    log.info(f"Found {len(folders)} gene folders\n")

    rows           = []
    skipped_genes  = []
    zero_row_genes = []
    totals = {"compound_het": 0, "glued": 0, "cis_multi": 0, "position_not_own": 0}

    for folder in folders:
        gene        = re.sub(r"base$", "", folder, flags=re.IGNORECASE)
        folder_path = os.path.join(IDBASE_DIR, folder)
        pub         = next((f for f in os.listdir(folder_path) if f.lower().endswith("pub.html")), None)

        if not pub:
            log.warning(f"  SKIP {gene}: no pub.html found")
            skipped_genes.append(gene)
            continue

        stats = {"compound_het": 0, "glued": 0, "cis_multi": 0, "position_not_own": 0}
        gene_rows = list(parse_pub_html(os.path.join(folder_path, pub), gene, log, stats))
        rows.extend(gene_rows)
        for k in totals:
            totals[k] += stats[k]

        missing_flags = sum(1 for r in gene_rows if r["missing_flag"])
        unknown_types = sum(1 for r in gene_rows if r["variant_type"] == "unknown")

        log.info(
            f"  {gene:<15}  {len(gene_rows):5d} rows"
            + (f"  [{stats['compound_het']} compound_het]" if stats["compound_het"] else "")
            + (f"  [{stats['glued']} glued_repaired]" if stats["glued"] else "")
            + (f"  [{missing_flags} missing_flag]" if missing_flags else "")
            + (f"  [{unknown_types} unknown_type]"  if unknown_types  else "")
            + (f"  [{stats['position_not_own']} position_not_own]" if stats["position_not_own"] else "")
        )

        if not gene_rows:
            zero_row_genes.append(gene)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    header = "gene\taccession\tallele_num\tvariant_type\tnotation\tpos_start\tpos_end\tmut_class\tmissing_flag\n"
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        f.write(header)
        for r in rows:
            f.write(
                f"{r['gene']}\t{r['accession']}\t{r['allele_num']}\t{r['variant_type']}\t"
                f"{r['notation']}\t{r['pos_start']}\t{r['pos_end']}\t{r['mut_class']}\t{r['missing_flag']}\n"
            )

    elapsed = datetime.now() - start
    log.info("")
    log.info("=" * 60)
    log.info("SUMMARY")
    log.info("=" * 60)
    log.info(f"  Genes processed          : {len(folders) - len(skipped_genes)}")
    log.info(f"  Genes skipped            : {len(skipped_genes)}  {skipped_genes or ''}")
    log.info(f"  Genes with 0 rows        : {len(zero_row_genes)}  {zero_row_genes or ''}")
    log.info(f"  Compound-het accessions  : {totals['compound_het']}")
    log.info(f"  Glued records repaired   : {totals['glued']}")
    log.info(f"  Alleles with cis changes : {totals['cis_multi']}")
    log.info(f"  Positions not own (empty): {totals['position_not_own']}")
    log.info(f"  Total rows               : {len(rows)}")
    log.info(f"  Output                   : {OUT_PATH}")
    log.info(f"  Log                      : {LOG_PATH}")
    log.info(f"  Elapsed                  : {str(elapsed).split('.')[0]}")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
