"""
normalize_mane.py
=================
Moves resolved records that are not on the gene's MANE Select transcript (MANE v1.5 GFF)
onto it, keeping the original description in the pair's alt_descriptions. dedup_key is
the record id and stays unchanged.

    relabel   the record is on another version of the MANE NM, its c. is purely coding and
              the old version's CDS equals the MANE CDS (FASTA): transcript and protein
              accessions are replaced, c. and p. are kept. No query.
    query     otherwise, a record with a genomic description is projected onto the MANE
              exons; a projection in the CDS or within 20 nt of an exon in an intron is
              queried as NC(NM_MANE):c. A record on another version without a genomic
              description is queried with its own c. The response is accepted when it is
              ok, on the MANE NM, and describes the same change (same NC g. as the record,
              or the same protein description when the record has no g.); c., r. and p.
              then come from the response. Anything else keeps the record and is logged.

Stages
------
--stage input   IN_DEDUP -> OUT_INPUT (Mutalyzer input, run_mutalyzer_NM.py columns)
--stage apply   IN_DEDUP, IN_PAIRS, IN_RESULTS -> OUT_DEDUP, OUT_PAIRS, OUT_LOG

Env vars
--------
IN_DEDUP, IN_PAIRS, IN_RESULTS, OUT_INPUT, OUT_DEDUP, OUT_PAIRS, OUT_LOG
MANE_GFF          MANE v1.5 ensembl_genomic GFF (gz)
MANE_NM_DIR       MANE NM FASTA directory (GENE_NM_x.fasta)
IDREFSEQ_NM_DIR   IDRefseq NM FASTA directory
"""

import argparse
import csv
import gzip
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

MAX_INTRONIC = 20
INPUT_COLS = ["gene", "accession", "allele_num", "sysname", "hgvs_input", "mut_type", "ref", "alt", "chrom",
              "pos_hg38", "strand", "dedup_key"]
LOG_COLS = ["gene", "dedup_key", "action", "reason", "old_description", "new_description", "query"]
RE_TX = re.compile(r"\(((?:NM|XM|NR|XR)_\d+\.\d+)\)|^((?:NM|XM|NR|XR)_\d+\.\d+):")
RE_G = re.compile(r"^(NC_\d+\.\d+):g\.(\d+)(?:_(\d+))?(.*)$")
_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, cols):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def mane_transcripts(gff):
    tx = {}
    with gzip.open(gff, "rt", encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or "tag=MANE_Select" not in p[8]:
                continue
            a = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
            if p[2] == "transcript":
                nm = re.search(r"RefSeq:(NM_\d+\.\d+)", a.get("Dbxref", ""))
                if nm:
                    tx[a["ID"]] = {"gene": a["gene_name"], "nm": nm.group(1), "np": "", "chrom": p[0],
                                   "strand": p[6], "exons": [], "cds": []}
            elif p[2] in ("exon", "CDS") and a.get("Parent") in tx:
                t = tx[a["Parent"]]
                t["exons" if p[2] == "exon" else "cds"].append((int(p[3]), int(p[4])))
                np_ = re.search(r"RefSeq:(NP_\d+\.\d+)", a.get("Dbxref", ""))
                if p[2] == "CDS" and np_:
                    t["np"] = np_.group(1)
    return {t["gene"]: t for t in tx.values()}


def translate(seq):
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON.get(seq[i:i + 3], "X")
        out.append(aa)
        if aa == "*":
            break
    return "".join(out)


def longest_cds(seq):
    best = ""
    for s in range(len(seq) - 2):
        if seq[s:s + 3] == "ATG":
            prot = translate(seq[s:])
            if prot.endswith("*") and 3 * len(prot) > len(best):
                best = seq[s:s + 3 * len(prot)]
    return best


def fasta_seq(directory, acc):
    hits = sorted(Path(directory).glob(f"*_{acc}.fasta"))
    if not hits:
        return ""
    return "".join(l.strip() for l in open(hits[0], encoding="utf-8") if not l.startswith(">")).upper()


def transcript_of(row):
    for field in ("c_hgvs", "normalized"):
        m = RE_TX.search(row.get(field, "") or "")
        if m:
            return next(g for g in m.groups() if g), field
    return "", ""


def comp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def g_to_c(t):
    ex = sorted(t["exons"], reverse=t["strand"] == "-")
    cds = sorted(t["cds"])
    cs, ce = cds[0][0], cds[-1][1]

    def tpos(x):
        n = 0
        for a, b in ex:
            if a <= x <= b:
                return n + (x - a + 1 if t["strand"] == "+" else b - x + 1)
            n += b - a + 1
        return None

    s0 = tpos(cs if t["strand"] == "+" else ce)
    last = tpos(ce if t["strand"] == "+" else cs) - s0 + 1
    lo, hi = min(a for a, _ in ex), max(b for _, b in ex)

    def c_of(x):
        if x < lo or x > hi:
            return None
        tp = tpos(x)
        if tp is not None:
            v = tp - s0 + 1 if tp - s0 + 1 > 0 else tp - s0
            return (f"*{v - last}" if v > last else str(v)), 0
        e, d = min(((e, x - e) for a, b in ex for e in (a, b)), key=lambda z: abs(z[1]))
        base, _ = c_of(e)
        off = d if t["strand"] == "+" else -d
        return f"{base}{'+' if off > 0 else '-'}{abs(off)}", abs(off)
    return c_of


def project(t, g_desc):
    m = RE_G.match(g_desc)
    if not m:
        return None, "no genomic description"
    a, b = int(m.group(2)), int(m.group(3) or m.group(2))
    rest = m.group(4)
    if t["strand"] == "-":
        a, b = b, a
        ms = re.match(r"^([ACGT])>([ACGT])$", rest)
        mi = re.match(r"^(ins|delins)([ACGT]+)$", rest)
        if ms:
            rest = f"{comp(ms.group(1))}>{comp(ms.group(2))}"
        elif mi:
            rest = mi.group(1) + comp(mi.group(2))
        elif rest not in ("del", "dup"):
            return None, "unprojectable change"
    c_of = g_to_c(t)
    pa, pb = c_of(a), c_of(b)
    if pa is None or pb is None:
        return None, "outside the MANE transcript"
    if max(pa[1], pb[1]) > MAX_INTRONIC:
        return None, "deep intronic on MANE"
    if any(x[0].startswith("-") or x[0].startswith("*") for x in (pa, pb)):
        return None, "UTR on MANE"
    return f"c.{pa[0]}{'' if a == b else '_' + pb[0]}{rest}", ""


class Planner:
    def __init__(self, dedup, mane, mane_dir, idr_dir):
        self.mane, self.mane_dir, self.idr_dir = mane, mane_dir, idr_dir
        self.nc = {}
        for r in dedup:
            m = re.match(r"^(NC_\d+\.\d+)", r.get("normalized", "") or r.get("c_hgvs", "") or "")
            ch = (r.get("chrom", "") or "").strip()
            if m and ch:
                self.nc.setdefault(ch if ch.startswith("chr") else "chr" + ch, m.group(1))
        self.cds_same = {}

    def same_cds(self, old, new):
        if (old, new) not in self.cds_same:
            a, b = fasta_seq(self.idr_dir, old), fasta_seq(self.mane_dir, new)
            self.cds_same[(old, new)] = bool(a and b) and longest_cds(a) == longest_cds(b) != ""
        return self.cds_same[(old, new)]

    def plan(self, r):
        gene = (r.get("mutalyzer_gene") or r["gene"]).strip()
        t = self.mane.get(gene) or self.mane.get(r["gene"].strip().upper())
        if not t:
            return None
        tx, field = transcript_of(r)
        if not tx or tx == t["nm"]:
            return None
        desc = r.get(field, "")
        c = desc.rsplit(":", 1)[-1]
        nc = self.nc.get(t["chrom"], "")
        same_nm = tx.split(".")[0] == t["nm"].split(".")[0]
        coding = c.startswith("c.") and not re.search(r"[-*+]", c[2:])
        if same_nm and coding and self.same_cds(tx, t["nm"]):
            return {"action": "relabel", "t": t, "tx": tx, "field": field, "desc": desc}
        g = r.get("normalized", "") if (r.get("normalized", "") or "").startswith("NC_") and ":g." in r.get("normalized", "") else ""
        if not nc:
            return {"action": "skip", "reason": "no NC accession for the chromosome", "desc": desc, "t": t}
        if g:
            pc, why = project(t, g)
            if not pc:
                return {"action": "skip", "reason": why, "desc": desc, "t": t}
            return {"action": "query", "t": t, "desc": desc, "g": g, "query": f"{nc}({t['nm']}):{pc}"}
        if same_nm and c.startswith("c."):
            return {"action": "query", "t": t, "desc": desc, "g": "", "query": f"{nc}({t['nm']}):{c}"}
        return {"action": "skip", "reason": "no genomic description", "desc": desc, "t": t}


def p_body(s):
    return s.split(":", 1)[1] if ":" in (s or "") else (s or "")


def stage_input(args, env):
    dedup = read_tsv(env["IN_DEDUP"])
    planner = Planner(dedup, mane_transcripts(env["MANE_GFF"]), env["MANE_NM_DIR"], env["IDREFSEQ_NM_DIR"])
    rows = []
    for r in dedup:
        pl = planner.plan(r)
        if pl and pl["action"] == "query":
            rows.append({**{k: r.get(k, "") for k in INPUT_COLS}, "hgvs_input": pl["query"]})
    write_tsv(env["OUT_INPUT"], rows, INPUT_COLS)
    print(f"records to query on MANE: {len(rows)} ({len({x['hgvs_input'] for x in rows})} distinct queries)")


def stage_apply(args, env):
    dedup = read_tsv(env["IN_DEDUP"])
    cols = list(dedup[0].keys()) if dedup else []
    pairs = read_tsv(env["IN_PAIRS"])
    pair_cols = list(pairs[0].keys()) if pairs else []
    results = {(r["gene"], r["dedup_key"]): r for r in read_tsv(env["IN_RESULTS"])}
    planner = Planner(dedup, mane_transcripts(env["MANE_GFF"]), env["MANE_NM_DIR"], env["IDREFSEQ_NM_DIR"])
    log, alt, counts = [], {}, Counter()
    for r in dedup:
        pl = planner.plan(r)
        if not pl:
            continue
        key = (r["gene"], r["dedup_key"])
        t = pl["t"]
        entry = {"gene": r["gene"], "dedup_key": r["dedup_key"], "old_description": pl["desc"], "query": pl.get("query", "")}
        if pl["action"] == "relabel":
            old_np = re.search(r"\((NP_\d+\.\d+)\)", r.get("p_hgvs", "") or "")
            for f in ("c_hgvs", "normalized", "r_hgvs"):
                if r.get(f) and pl["tx"] in r[f]:
                    r[f] = r[f].replace(pl["tx"], t["nm"])
            if r.get("p_hgvs"):
                r["p_hgvs"] = r["p_hgvs"].replace(pl["tx"], t["nm"])
                if old_np and t["np"]:
                    r["p_hgvs"] = r["p_hgvs"].replace(old_np.group(1), t["np"])
            r["mane_select"] = "MANE Select"
            alt[key] = pl["desc"]
            entry.update(action="relabel", reason="same CDS, coding c.", new_description=r[pl["field"]])
        elif pl["action"] == "query":
            res = results.get(key, {})
            ok = res.get("status") == "ok"
            res_tx, _ = transcript_of({"c_hgvs": res.get("c_hgvs", "")})
            same = (res.get("nc_hgvs", "") == pl["g"]) if pl["g"] else (p_body(res.get("p_hgvs", "")) == p_body(r.get("p_hgvs", "")))
            if ok and res_tx == t["nm"] and same:
                for f in ("c_hgvs", "r_hgvs", "p_hgvs", "protein_pos_first"):
                    r[f] = res.get(f, "")
                r["mane_select"] = "MANE Select"
                r["dedup_key_source"] = "c_hgvs"
                alt[key] = pl["desc"]
                entry.update(action="normalized", reason="MANE response, same change", new_description=r["c_hgvs"])
            else:
                why = ("query status " + (res.get("status") or "missing")) if not ok else (
                    "response not on MANE" if res_tx != t["nm"] else "response describes a different change")
                entry.update(action="kept", reason=why, new_description="")
        else:
            entry.update(action="kept", reason=pl["reason"], new_description="")
        counts[entry["action"]] += 1
        log.append(entry)
    for p in pairs:
        k = (p["gene"], p["dedup_key"])
        if k in alt:
            p["alt_descriptions"] = " | ".join(dict.fromkeys([x for x in p["alt_descriptions"].split(" | ") if x] + [alt[k]]))
    write_tsv(env["OUT_DEDUP"], dedup, cols)
    write_tsv(env["OUT_PAIRS"], pairs, pair_cols)
    write_tsv(env["OUT_LOG"], log, LOG_COLS)
    print(f"records off MANE Select: {len(log)}; " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    for reason, n in Counter(e["reason"] for e in log if e["action"] == "kept").most_common():
        print(f"  kept ({reason}): {n}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--stage", choices=["input", "apply"], required=True)
    args = ap.parse_args()
    env = defaultdict(str, os.environ)
    if args.stage == "input":
        stage_input(args, env)
    else:
        stage_apply(args, env)


if __name__ == "__main__":
    main()
