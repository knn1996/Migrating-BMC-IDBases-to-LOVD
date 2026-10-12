import os
import re
import time
import urllib.request
import pandas as pd

THESIS_DIR  = os.environ.get("THESIS_DIR", r"D:\Assignment\Thesis")
OUT_ROOT    = os.environ.get("OUT_ROOT", os.path.join(THESIS_DIR, "04_Mutation_Processing", "Output"))
IN_TSV      = os.environ.get("IN_TSV", os.path.join(OUT_ROOT, "Step8_Rescue", "lovd_flat_final_cols_vario.tsv"))
OUT_DIR     = os.environ.get("OUT_DIR", os.path.join(OUT_ROOT, "Step9_NMD"))
GB_CACHE    = os.environ.get("GB_CACHE", os.path.join(OUT_ROOT, "cache", "genbank"))
NCBI_KEY    = os.environ.get("NCBI_API_KEY", "")
AUDIT_GENES = {g.strip().upper() for g in os.environ.get("AUDIT_GENES", "").split(",") if g.strip()}

NMD_THRESHOLD   = int(os.environ.get("NMD_THRESHOLD", "50"))
BORDERLINE_NT   = int(os.environ.get("BORDERLINE_NT", "5"))
START_PROXIMAL  = int(os.environ.get("START_PROXIMAL_NT", "150"))

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=gb&retmode=text"

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for x in _b:
    for y in _b:
        for z in _b:
            CODON[x + y + z] = _aa[_i]
            _i += 1


def fetch_genbank(acc):
    os.makedirs(GB_CACHE, exist_ok=True)
    path = os.path.join(GB_CACHE, acc + ".gb")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return open(path, encoding="utf-8").read()
    url = EUTILS.format(acc) + (f"&api_key={NCBI_KEY}" if NCBI_KEY else "")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                txt = r.read().decode("utf-8", errors="replace")
            if txt.startswith("LOCUS"):
                open(path, "w", encoding="utf-8").write(txt)
                time.sleep(0.12 if NCBI_KEY else 0.4)
                return txt
        except Exception:
            pass
        time.sleep(1.5 * (attempt + 1))
    return None


def parse_genbank(txt):
    exons, cds = [], None
    for m in re.finditer(r"^ {5}exon {12}<?(\d+)\.\.>?(\d+)", txt, re.M):
        exons.append((int(m.group(1)), int(m.group(2))))
    m = re.search(r"^ {5}CDS {13}<?(\d+)\.\.>?(\d+)", txt, re.M)
    if m:
        cds = (int(m.group(1)), int(m.group(2)))
    seq = ""
    if "\nORIGIN" in txt:
        seq = re.sub(r"[^acgtnACGTN]", "", txt.split("\nORIGIN", 1)[1].split("\n//")[0]).upper()
    return sorted(exons), cds, seq


def last_coding_exon(exons, cds):
    return max(i for i, (a, b) in enumerate(exons) if a <= cds[1])


def translate_stop(seq, start):
    i = start
    while i + 3 <= len(seq):
        if CODON.get(seq[i:i + 3], "X") == "*":
            return i
        i += 3
    return None


RE_POS = r"(\d+)"
RE_C = re.compile(r"^c\.(\d+)(?:_(\d+))?(.*)$")


def apply_variant(seq, cds_start0, cdna):
    m = RE_C.match(cdna)
    if not m:
        return None, "not_assessable_position"
    a = int(m.group(1))
    b = int(m.group(2)) if m.group(2) else a
    rest = m.group(3)
    i0 = cds_start0 + a - 1
    i1 = cds_start0 + b - 1
    if i0 < 0 or i1 >= len(seq):
        return None, "out_of_range"
    ref_check = ""
    ms = re.match(r"^([ACGT])>([ACGT])$", rest)
    if ms:
        ref_check = "ok" if seq[i0] == ms.group(1) else "REF_MISMATCH"
        return (seq[:i0] + ms.group(2) + seq[i0 + 1:], i0, 0, ref_check), "ok"
    md = re.match(r"^delins([ACGT]+)$", rest)
    if md:
        ins = md.group(1)
        return (seq[:i0] + ins + seq[i1 + 1:], i0, len(ins) - (i1 - i0 + 1), ref_check), "ok"
    md = re.match(r"^del([ACGT]*)$", rest)
    if md:
        if md.group(1):
            ref_check = "ok" if seq[i0:i1 + 1] == md.group(1) else "REF_MISMATCH"
        return (seq[:i0] + seq[i1 + 1:], i0, -(i1 - i0 + 1), ref_check), "ok"
    md = re.match(r"^dup([ACGT]*)$", rest)
    if md:
        dup = seq[i0:i1 + 1]
        if md.group(1):
            ref_check = "ok" if dup == md.group(1) else "REF_MISMATCH"
        return (seq[:i1 + 1] + dup + seq[i1 + 1:], i1 + 1, len(dup), ref_check), "ok"
    md = re.match(r"^ins([ACGT]+)$", rest)
    if md:
        ins = md.group(1)
        return (seq[:i0 + 1] + ins + seq[i0 + 1:], i0 + 1, len(ins), ref_check), "ok"
    md = re.match(r"^([ACGT]+)\[(\d+)\]$", rest)
    if md:
        unit, k = md.group(1), int(md.group(2))
        j, r = i0, 0
        while seq[j:j + len(unit)] == unit:
            j += len(unit)
            r += 1
        if r == 0:
            return None, "REF_MISMATCH_repeat"
        return (seq[:i0] + unit * k + seq[j:], i0, (k - r) * len(unit), "ok"), "ok"
    return None, "unparsed_change"


def ptc_kind(protein, cdna):
    p = protein or ""
    if "fs" in p:
        return "frameshift"
    if re.search(r"(Ter|\*|X)\)?$", p) and "ext" not in p:
        return "nonsense"
    if re.search(r"(del|ins|dup)", cdna or "") and not re.search(r"(del|ins|dup)", p):
        return "frameshift_candidate"
    return None


def assess(row, tx):
    out = {}
    cdna, prot = row["VariantOnTranscript/DNA"], row["VariantOnTranscript/Protein"]
    kind = ptc_kind(prot, cdna)
    out["ptc_kind"] = kind or ""
    if not kind:
        out["nmd_call"] = "not_ptc"
        return out
    if re.search(r"c\.[-*]|[+-]\d", cdna):
        out["nmd_call"] = "not_assessable_noncoding_or_intronic"
        return out
    if tx is None:
        out["nmd_call"] = "not_assessable_no_transcript"
        return out
    exons, cds, seq = tx["exons"], tx["cds"], tx["seq"]
    if not exons or not cds or not seq:
        out["nmd_call"] = "not_assessable_no_structure"
        return out
    cds0 = cds[0] - 1
    res, status = apply_variant(seq, cds0, cdna)
    if res is None:
        out["nmd_call"] = "not_assessable_" + status
        return out
    mseq, vpos0, delta, ref_check = res
    out["ref_check"] = ref_check
    ref_stop0 = cds[1] - 3
    mstop0 = translate_stop(mseq, cds0)
    if mstop0 is None:
        out["nmd_call"] = "no_stop_nonstop_decay"
        return out
    ref_stop_mut0 = ref_stop0 + delta if vpos0 <= ref_stop0 else ref_stop0
    if mstop0 >= ref_stop_mut0:
        out["nmd_call"] = "no_PTC_stop_at_or_beyond_normal"
        out["ptc_c"] = mstop0 - cds0 + 1 - (delta if vpos0 <= mstop0 else 0)
        return out
    ptc_c_ref = mstop0 - cds0 + 1 - (delta if vpos0 <= mstop0 else 0)
    out["ptc_c"] = ptc_c_ref
    out["ptc_aa"] = (mstop0 - cds0) // 3 + 1
    m_ter = re.search(r"(\d+)\D*fs\D*(\d+)", prot or "")
    if m_ter:
        out["mutalyzer_stop_aa"] = int(m_ter.group(1)) + int(m_ter.group(2)) - 1
        out["stop_agrees_mutalyzer"] = "yes" if out["mutalyzer_stop_aa"] == out["ptc_aa"] else "NO"
    elif kind == "nonsense":
        m_n = re.search(r"(\d+)(?:Ter|\*|X)\)?$", prot or "")
        if m_n:
            out["mutalyzer_stop_aa"] = int(m_n.group(1))
            out["stop_agrees_mutalyzer"] = "yes" if out["mutalyzer_stop_aa"] == out["ptc_aa"] else "NO"
    if len(exons) == 1:
        out["nmd_call"] = "NMD_escape_single_exon"
        return out
    last_start0 = exons[last_coding_exon(exons, cds)][0] - 1
    last_start_mut0 = last_start0 + delta if vpos0 < last_start0 else last_start0
    junction_mut0 = last_start_mut0 - 1
    if mstop0 >= last_start_mut0:
        out["nmd_call"] = "NMD_escape_last_exon"
        out["dist_to_last_junction_nt"] = junction_mut0 - mstop0 + 1
        return out
    dist = junction_mut0 - mstop0 + 1
    out["dist_to_last_junction_nt"] = dist
    out["nmd_call"] = "NMD_predicted" if dist > NMD_THRESHOLD else "NMD_escape_50nt_rule"
    flags = []
    if abs(dist - NMD_THRESHOLD) <= BORDERLINE_NT:
        flags.append("borderline")
    if ptc_c_ref <= START_PROXIMAL:
        flags.append("start_proximal")
    if tx["ref_stop_triggers"]:
        flags.append("normal_stop_also_upstream_of_rule")
    out["flags"] = ";".join(flags)
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    df = pd.read_csv(IN_TSV, sep="\t", dtype=str).fillna("")
    if AUDIT_GENES:
        df = df[df["_gene"].str.upper().isin(AUDIT_GENES)]
    txs = sorted({t for t in df["transcriptid"] if t})
    structs, srows = {}, []
    for n, acc in enumerate(txs, 1):
        txt = fetch_genbank(acc)
        if not txt:
            structs[acc] = None
            srows.append({"transcriptid": acc, "status": "fetch_failed"})
            continue
        exons, cds, seq = parse_genbank(txt)
        ok = bool(exons and cds and seq)
        trig = False
        info = {"transcriptid": acc, "status": "ok" if ok else "no_structure", "n_exons": len(exons)}
        if ok:
            last_start = exons[last_coding_exon(exons, cds)][0]
            junction = last_start - 1
            ref_stop_first = cds[1] - 2
            trig = len(exons) > 1 and ref_stop_first < last_start and (junction - ref_stop_first + 1) > NMD_THRESHOLD
            info.update({
                "cds_start_tx": cds[0], "cds_end_tx": cds[1],
                "last_exon_start_c": last_start - cds[0] + 1,
                "last_junction_c": junction - cds[0] + 1,
                "nmd_boundary_c": junction - cds[0] + 1 - NMD_THRESHOLD,
                "normal_stop_c": ref_stop_first - cds[0] + 1,
                "normal_stop_triggers_rule": trig,
            })
        structs[acc] = {"exons": exons, "cds": cds, "seq": seq, "ref_stop_triggers": trig} if ok else None
        srows.append(info)
        if n % 25 == 0:
            print(f"  fetched {n}/{len(txs)} transcripts")

    rows = []
    for r in df.to_dict("records"):
        a = assess(r, structs.get(r["transcriptid"]))
        call = a.get("nmd_call", "")
        rows.append({
            "_gene": r["_gene"], "_dedup_key": r["_dedup_key"], "transcriptid": r["transcriptid"],
            "VariantOnTranscript/DNA": r["VariantOnTranscript/DNA"],
            "rna_current": r.get("VariantOnTranscript/RNA", ""),
            "protein_current": r["VariantOnTranscript/Protein"],
            "ptc_kind": a.get("ptc_kind", ""), "ref_check": a.get("ref_check", ""),
            "ptc_c": a.get("ptc_c", ""), "ptc_aa": a.get("ptc_aa", ""),
            "mutalyzer_stop_aa": a.get("mutalyzer_stop_aa", ""),
            "stop_agrees_mutalyzer": a.get("stop_agrees_mutalyzer", ""),
            "dist_to_last_junction_nt": a.get("dist_to_last_junction_nt", ""),
            "nmd_call": call, "flags": a.get("flags", ""),
            "rna_proposed": "r.0?" if call == "NMD_predicted" else r.get("VariantOnTranscript/RNA", ""),
            "protein_proposed": "p.0?" if call == "NMD_predicted" else r["VariantOnTranscript/Protein"],
        })
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT_DIR, "nmd_annotation.tsv"), sep="\t", index=False)
    pd.DataFrame(srows).to_csv(os.path.join(OUT_DIR, "transcript_structure.tsv"), sep="\t", index=False)
    review = out[(out["nmd_call"].str.startswith("NMD") & (out["flags"] != "")) |
                 (out["stop_agrees_mutalyzer"] == "NO") | (out["ref_check"].str.contains("MISMATCH")) |
                 (out["ptc_kind"] == "frameshift_candidate")]
    review.to_csv(os.path.join(OUT_DIR, "nmd_review.tsv"), sep="\t", index=False)

    ptc = out[out["ptc_kind"] != ""]
    print(f"\nVariants: {len(out)}  with PTC kind: {len(ptc)}  transcripts: {len(txs)}")
    print(ptc.groupby(["ptc_kind", "nmd_call"]).size().to_string())
    print("\nflags:", ptc["flags"].str.split(";").explode().replace("", pd.NA).dropna().value_counts().to_dict())
    print("stop agrees with Mutalyzer:", ptc["stop_agrees_mutalyzer"].value_counts().to_dict())
    print("ref_check:", ptc["ref_check"].value_counts().to_dict())
    print(f"\nOutputs in {OUT_DIR}: nmd_annotation.tsv, transcript_structure.tsv, nmd_review.tsv")


if __name__ == "__main__":
    main()
