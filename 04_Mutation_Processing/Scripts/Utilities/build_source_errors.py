"""
build_source_errors.py
======================
Classifies the source changes behind the inconsistent final pairs (Step10_Audit/audit_pairs.tsv),
plus the named residuals, as suspected IDbases source errors and writes private/source_errors.tsv
(git-ignored: record-level candidates are not published until the curator has reviewed them).

Each distinct source change (gene, g., c., p.) is tested on the gene's MANE Select CDS at shift 0
and at the gene's verified legacy shift (refs/legacy_shift_evidence.tsv), and on every cached
non-MANE GenBank transcript, with legacy_remap.mutate / reconcile. Categories, in order:

    pipeline           the source c. gives the source p. at some shift or isoform our record is not on
    not_error          the pair is consistent (notation equivalence)
    exon_boundary      our genomic record crosses an exon/intron junction the source c. writes as exonic
    g_c_mismatch       the source g. and c. are different change types, or g. maps elsewhere than c.
    notation           our record is the source c.; the source p. has the wrong type, extent or stop
    c_p_mismatch       our record is the source c.; it cannot give the source p.
    ref_base_mismatch  the bases stated in the source c. are absent at every tested shift
    untestable         no parseable change, protein or MANE CDS

pipeline_fault_possible is yes when an isoform, an unverified shift, a non-GenBank CDS, an intronic
record or a literal record in a legacy-shift gene could explain the discrepancy; an intrinsic c./p.
or g./c. type conflict is a source error regardless. our_record_wrong marks a literal kept record in
a verified-shift gene whose source c. at that shift places the source p. position. severity is low
when only the frameshift stop count or only the source g. differs, high otherwise.

private/source_errors.tsv holds the candidates for review (status candidate): every source change of
the rows with a genuine source issue (pipeline_fault_possible no, or an intrinsic conflict).

Usage:
    python Utilities/build_source_errors.py --scan <scan_changes.tsv> [--out private/source_errors.tsv] [--report-dir <dir>]
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import yaml

SC = Path(__file__).resolve().parents[1]
for p in ("Utilities", "Step10_Audit", "Step7_Merging"):
    sys.path.insert(0, str(SC / p))
import legacy_remap as R
import build_legacy_remap as B
import audit_ledger as L

OUT_COLS = ["source_change_id", "gene", "accession", "allele", "source_c", "source_p", "category", "severity",
            "our_record_wrong", "status", "evidence"]


def genbank_cds(text):
    m = re.search(r"^ {5}CDS {13}<?(\d+)\.\.>?(\d+)(.*?)(?=^ {5}\S|^ORIGIN)", text, re.M | re.S)
    if not m or "\nORIGIN" not in text:
        return None, None
    g = re.search(r'/gene="([^"]+)"', m.group(3))
    seq = re.sub(r"[^acgtnACGTN]", "", text.split("\nORIGIN", 1)[1].split("\n//")[0]).upper()
    return (g.group(1).strip().upper() if g else None), (seq[int(m.group(1)) - 1:int(m.group(2))] or None)


def transcript_store(genbank_dir):
    by_gene = {}
    for path in sorted(Path(genbank_dir).glob("*.gb")):
        gene, cds = genbank_cds(path.read_text(encoding="utf-8"))
        if gene and cds:
            by_gene.setdefault(gene, []).append((path.stem, cds))
    return by_gene


ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("--scan", type=Path, required=True)
ap.add_argument("--out", type=Path, default=SC.parents[1] / "private" / "source_errors.tsv")
ap.add_argument("--report-dir", type=Path, default=None)
args = ap.parse_args()
with open(os.environ.get("CONFIG", SC / "config" / "config.yaml"), encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
O = Path(cfg["output_root"])

KMAX = 1500
rd = lambda p: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
pairs = rd(O / "Step10_Audit/audit_pairs.tsv")
var = L.disambiguate_alleles(rd(O / "Step10_Audit/audit_variants.tsv")); sp = L.SourceProteins(var)
final = rd(O / "Step8_Rescue/lovd_flat_with_patients_final.tsv")
aug = rd(O / "Step8_Rescue/augmented_merged_variants.tsv")
dedup = rd(O / "Step8_Rescue/dedup_merged_variants.resolved.tsv")
ev = rd(SC / "refs/legacy_shift_evidence.tsv")
scan = rd(args.scan)
symbol = {g.upper(): (x.mutalyzer_gene[x.mutalyzer_gene != ""].mode().iloc[0] if (x.mutalyzer_gene != "").any() else g.upper()) for g, x in dedup.groupby("gene")}
mane = B.mane_select(Path(cfg["thesis_dir"]) / B.MANE_GFF)
offs = {r.gene.upper(): r.cds_start for r in pd.read_csv(O / "Step2_RefCheck/coding_offset_scan_MANE.csv", dtype=str, keep_default_na=False).itertuples()}
store_by_gene = transcript_store(O / "cache/genbank")
verified = {r.gene: int(r.shift) for r in ev.itertuples() if r.status == "kept"}
scan_shift = {g: int(x.gene_shift.iloc[0]) for g, x in scan.groupby("gene") if x.gene_shift.iloc[0]}
scan_cls = dict(zip(scan.sid, scan["class"])); scan_ks = dict(zip(scan.sid, scan.change_shifts))

src_rows = {}
for r in var.itertuples():
    src_rows[(r.gene.upper(), r.accession, r.allele)] = r
fin = {(i, d): r for i, d, r in zip(final["Individual/ID"], final["VariantOnTranscript/DNA"], final.itertuples())}
trackA = defaultdict(set)
for s, t, c in zip(aug.source_change_id, aug.source_track, aug.c_hgvs):
    if t == "NG_IDRefseq" and c: trackA[s].add(c.rsplit(":", 1)[-1])

targets = []
for r in pairs.itertuples():
    if r.label.startswith("inconsistent"):
        for s in r.source_change_ids.split(" ; "):
            targets.append((s, r.individual_id, r.variant_key, r.track, r.label))
extra = [("CD40L", "C0047"), ("WAS", "W0009")]
for r in pairs.itertuples():
    for s in r.source_change_ids.split(" ; "):
        if tuple(s.split("|")[:2]) in extra and not r.label.startswith("inconsistent"):
            targets.append((s, r.individual_id, r.variant_key, r.track, r.label))

refs = {}
def mane_ref(gene):
    if gene not in refs:
        nm = mane.get(symbol.get(gene, gene), "")
        refs[gene] = (nm, B.reference(gene, nm, O / "cache/genbank", cfg["dna_seq"]["mane_nm"], offs) if nm else None)
    return refs[gene]

def pos_of(c):
    m = re.match(r"^c\.(-?\*?\d+)([+-]\d+)?(?:_(-?\*?\d+)([+-]\d+)?)?", c)
    return m.groups() if m else None

def isoform_hit(gene, ch, prots, nm):
    hits = []
    for name, cds in store_by_gene.get(gene, []) + store_by_gene.get(symbol.get(gene, gene).upper(), []):
        if nm and name.startswith(nm.split(".")[0]):
            continue
        for k in range(-KMAX, KMAX + 1):
            for p in prots:
                if R.reconcile(cds, 0, ch, p, k) or (k % 3 == 0 and R.reconcile(cds, 0, ch, (p[0], p[1] + k // 3, p[2], p[3]), k)):
                    hits.append((name, k)); break
            if hits and hits[-1][0] == name:
                break
    return hits

THREE_OF = {v: k for k, v in R.THREE.items()}
def aa3(a): return THREE_OF.get(a, "Xaa")
def describe_p(seq, cds0, mut, frameshift):
    wt, mt = R.translate(seq[cds0:]), R.translate(mut[cds0:])
    n = next((i for i in range(min(len(wt), len(mt))) if wt[i] != mt[i]), None)
    if n is None:
        return "p.(=)" if len(wt) == len(mt) else "p.?"
    if mt[n] == "*":
        return f"p.({aa3(wt[n])}{n + 1}Ter)"
    if frameshift:
        m = mt.find("*", n)
        return f"p.({aa3(wt[n])}{n + 1}{aa3(mt[n])}fsTer{m - n + 1 if m >= 0 else '?'})"
    s = 0
    while s < min(len(wt), len(mt)) - n and wt[len(wt) - 1 - s] == mt[len(mt) - 1 - s]:
        s += 1
    d, ins = wt[n:len(wt) - s], mt[n:len(mt) - s]
    a = f"{aa3(wt[n])}{n + 1}"
    b = f"_{aa3(wt[n + len(d) - 1])}{n + len(d)}" if len(d) > 1 else ""
    if len(d) == 1 and len(ins) == 1:
        return f"p.({a}{aa3(ins)})"
    if d and not ins:
        return f"p.({a}{b}del)"
    if ins and not d:
        return f"p.({aa3(wt[n - 1])}{n}_{aa3(wt[n])}{n + 1}ins{''.join(aa3(x) for x in ins)})"
    return f"p.({a}{b}delins{''.join(aa3(x) for x in ins)})"

KIND_NAME = {"sub": "a substitution", "del": "a deletion", "dup": "a duplication", "ins": "an insertion", "delins": "a deletion-insertion"}
groups = defaultdict(list)
for s, iid, key, track, label in targets:
    gene, acc, al, g, c = s.split("|", 4)
    groups[(gene, g, c, tuple(sp.of(s)))].append((s, iid, key, track, label))

rows = []
for (gene, g, c, prots_s), members in sorted(groups.items()):
    nm, ref = mane_ref(gene)
    ch = R.parse_c(c) if c != "-" else None
    prots = [p for p in (R.parse_p(x) for x in prots_s) if p]
    gk = verified.get(gene, scan_shift.get(gene))
    ks = [0] + ([gk] if gk else [])
    labels = sorted({m[4] for m in members})
    keys = sorted({m[2] for m in members}); tracks = sorted({m[3] for m in members})
    fr = [fin.get((m[1], m[2])) for m in members]
    our_p = sorted({f._26 for f in fr if f is not None and f._26})
    tier = sorted({f.confidence_tier for f in fr if f is not None})
    tA = sorted(set().union(*[trackA[m[0]] for m in members]))
    sids = [m[0] for m in members]
    cls_scan = sorted({scan_cls.get(s, "") for s in sids})
    info = {}
    kind_of = lambda s: "" if s in ("", "-") else ("sub" if ">" in s else next((x for x in ("delins", "del", "dup", "ins") if x in s), ""))
    gk0, ck0 = kind_of(g), kind_of(c)
    if (not ch or not prots or not ref) and gk0 and ck0 and gk0 != ck0 and {gk0, ck0} != {"dup", "ins"}:
        cat = "g_c_mismatch"; fault = "no"
        evid = f"The source itself is contradictory: {g} is {KIND_NAME[gk0]} but {c} is {KIND_NAME[ck0]}; our record follows the g. ({', '.join(keys)}, {', '.join(our_p) or 'no p.'})."
        info = {"sev": "high", "intrinsic": True}
    elif not ch or not prots or not ref:
        cat = "untestable"
        why = "source c. not a simple parseable change" if not ch else ("no parseable source p." if not prots else "no MANE CDS offline")
        evid = f"Cannot be tested offline: {why}."
        fault = "yes" if ch and prots else "no"
    else:
        seq, cds0, src = ref
        self_k = [k for k in ks if any(R.reconcile(seq, cds0, ch, p, k) for p in prots)]
        ref_ok = [k for k in ks if R.mutate(seq, cds0, ch, k) is not None]
        iso = [] if self_k else isoform_hit(gene, ch, prots, nm)
        unver = any(x in ("unverified_shift",) for x in cls_scan)
        our_pos = [pos_of(k) for k in keys]
        exon_b = (ch[1] != ch[2] and any(pp and ((pp[1] is None) != (pp[3] is None)) and pp[2] for pp in our_pos))
        def shifted_pos(k): return (ch[1] + k, ch[2] + k)
        def plain(c2):
            pp = pos_of(c2)
            if not pp or pp[1] or pp[3] or "*" in (pp[0] or "") or (pp[0] or "").startswith("-"): return None
            return int(pp[0])
        tA_disagree = bool(tA) and all(plain(x) is None or all(abs(plain(x) - shifted_pos(k)[0]) > 3 for k in ks) for x in tA)
        same_k = [k for k in ks if (lambda m: m is not None and any(R.parse_c(kk) and R.mutate(seq, cds0, R.parse_c(kk), 0) == m for kk in keys))(R.mutate(seq, cds0, ch, k))]
        start_k = [k for k in ks if any((pp := pos_of(kk)) and not pp[1] and (pp[0] or "").isdigit() and int(pp[0]) == ch[1] + k for kk in keys)]
        notation_label = any(l.split(":")[-1] in ("FS_STOP_DIFF", "FS_START_DIFF", "TYPE_DIFF_SAME_POS") for l in labels)
        is_inconsistent = any(l.startswith("inconsistent") for l in labels)
        gkind = "" if g in ("", "-") else ("sub" if ">" in g else next((t for t in ("delins", "del", "dup", "ins") if t in g), ""))
        kind_eq = lambda a, b: a == b or {a, b} == {"dup", "ins"}
        intronic_ours = any(pp and (pp[1] or pp[3] or (pp[0] or "").startswith("-") or "*" in (pp[0] or "")) for pp in our_pos) and not exon_b
        weak = src != "genbank" or intronic_ours
        what = 'stop' if any('STOP' in l for l in labels) else ('frameshift start' if any('START' in l for l in labels) else 'change type or extent')
        shift_txt = f"shift 0{'' if not gk else f' or {gk:+d}'}"
        def indel_len():
            if ch[0] == "sub": return 0
            if ch[0] in ("ins",): return len(ch[3])
            if ch[0] == "delins": return None
            return len(ch[3]) if ch[3] and ch[1] == ch[2] else ch[2] - ch[1] + 1
        n_ind = indel_len()
        pk = {pp[0] for pp in prots}
        type_conflict = (ch[0] == "sub" and pk <= {"fs", "inframe"}) or (n_ind not in (None, 0) and n_ind % 3 != 0 and pk <= {"sub", "inframe"}) or (n_ind not in (None, 0) and n_ind % 3 == 0 and pk <= {"fs"})
        pos_agree = any((ch[1] + k - 1) // 3 - 1 <= pp[1] <= (ch[2] + k) // 3 + 2 for k in ks for pp in prots)
        one_letter = any(re.match(r"^p\.\(?[A-Z*]\d", s) for s in prots_s)
        used_k = (same_k or start_k or [None])[0]
        suspect = bool(gk) and used_k == 0 and one_letter
        gv = verified.get(gene)
        pos_ok_gv = bool(gv) and R.mutate(seq, cds0, ch, gv) is not None and any(
            (ch[1] + gv - 1) // 3 - 1 <= pp[1] <= (ch[2] + gv) // 3 + 2 and 0 < pp[1] <= len(R.translate(seq[cds0:])) and R.translate(seq[cds0:])[pp[1] - 1] == pp[2] for pp in prots)
        our_wrong = bool(gv) and used_k == 0 and pos_ok_gv and is_inconsistent
        def fault_of(base):
            if our_wrong: return "yes"
            if type_conflict: return "no"
            if not pos_agree or suspect: return "yes"
            return base
        if self_k and not (set(same_k) & set(self_k)) and not (set(start_k) & set(self_k)) and is_inconsistent:
            cat = "pipeline"; fault = "yes"
            evid = f"Source c. {c} gives source {'/'.join(prots_s)} on {nm} at shift {self_k[0]:+d}, but our record {', '.join(keys)} is not on that shift: a pipeline issue, not a source error."
        elif not is_inconsistent and not exon_b:
            cat = "not_error"; fault = "no"
            evid = f"Pair is consistent ({', '.join(labels)}); notation equivalence, not a source error."
        elif exon_b:
            cat = "exon_boundary"; fault = "no"
            evid = f"Source {g} removes bases across an exon/intron junction ({', '.join(keys)} on {nm}); the source {c} writes it as a contiguous exonic range, which on the transcript describes a different change."
        elif gkind and not kind_eq(gkind, ch[0]):
            cat = "g_c_mismatch"; fault = "no"
            evid = f"The source itself is contradictory: {g} is {KIND_NAME[gkind]} but {c} is {KIND_NAME[ch[0]]}."
        elif (same_k or (self_k and start_k)) and notation_label:
            k0 = (same_k or start_k)[0]
            cat = "notation"; fault = fault_of("yes" if unver or gene == "CYBB" or weak else "no")
            evid = f"Our record {', '.join(keys)} on {nm} is the source {c} (shift {k0:+d}) and translates to {', '.join(our_p) or '?'}; source {'/'.join(prots_s)} has the right position but the wrong {what}."
        elif same_k:
            cat = "c_p_mismatch"; fault = fault_of("yes" if unver or weak else "no")
            evid = f"Our record {', '.join(keys)} on {nm} is the source {c} (shift {same_k[0]:+d}) and translates to {', '.join(our_p) or '?'}; it cannot give source {'/'.join(prots_s)} on {nm}{' or on any cached isoform' if not iso else ''}."
        elif self_k:
            cat = "pipeline"; fault = "yes"
            evid = f"Source c. {c} gives source {'/'.join(prots_s)} on {nm} at shift {self_k[0]:+d}, but our record {', '.join(keys)} is a different change: a pipeline issue, not a source error."
        elif iso:
            cat = "pipeline"; fault = "yes"
            evid = f"Source c. {c} gives source {'/'.join(prots_s)} on non-MANE transcript {iso[0][0]} at shift {iso[0][1]:+d}; the source follows another isoform and our record {', '.join(keys)} is mis-numbered."
        elif ch[3] and not ref_ok:
            cat = "ref_base_mismatch"; fault = "yes" if weak or unver else "no"
            evid = f"The base(s) stated in source {c} are not at c.{ch[1]} of {nm} ({shift_txt})."
        elif tA_disagree:
            cat = "g_c_mismatch"; fault = "yes"
            evid = f"Source {g} maps to {', '.join(tA)} on {nm}, which is not source {c} at {shift_txt}; g. and c. point to different positions (our genomic mapping may also be at fault)."
        elif notation_label:
            cat = "notation"; fault = fault_of("yes" if unver or gene == "CYBB" or weak else "no")
            evid = f"Source c. {c} gives {', '.join(our_p) or '?'} on {nm}; source {'/'.join(prots_s)} has the right position but the wrong {what}."
        else:
            cat = "c_p_mismatch"; fault = fault_of("yes" if unver or weak else "no")
            evid = f"Source c. {c} gives {', '.join(our_p) or '?'} on {nm}; it cannot give source {'/'.join(prots_s)} at {shift_txt}{' or on any cached isoform' if not iso else ''}."
        intrinsic = type_conflict or (gkind and not kind_eq(gkind, ch[0]))
        exp_c = exp_p = ""
        if our_wrong:
            fault = "yes"
            exp_c = R.shift_notation(c, gv) or ""
            exp_p = describe_p(seq, cds0, R.mutate(seq, cds0, ch, gv), n_ind is not None and n_ind % 3 != 0)
            if not intrinsic:
                cat = "pipeline"
                evid = f"Our record {', '.join(keys)} keeps the literal numbering in a {gene} {gv:+d} legacy-shift gene; at {gv:+d} the source {c} is {exp_c} on {nm}, giving {exp_p} at the source position. Re-check the source after the pipeline fix."
            elif gkind and not kind_eq(gkind, ch[0]):
                evid = f"Only the source g. is off: {g} is {KIND_NAME[gkind]}, but the source {c} and {'/'.join(prots_s)} agree at the {gene} legacy shift {gv:+d} ({exp_c} on {nm}, {exp_p})."
            else:
                evid = f"At the {gene} legacy shift {gv:+d} the source {c} is {exp_c} on {nm} and gives {exp_p}; source {'/'.join(prots_s)} names that position but a {'/'.join(sorted(pk))} change, which a {n_ind}-nt {ch[0]} cannot produce."
        sev = ""
        if cat not in ("pipeline", "not_error", "untestable"):
            stop_only = labels and all(l.endswith("FS_STOP_DIFF") for l in labels)
            if cat == "g_c_mismatch" and our_wrong:
                sev = "low"
            elif cat == "notation" and stop_only:
                sev = "low"
                L0 = len(R.translate(seq[cds0:])) - 1
                for q in our_p:
                    mm = re.search(r"\(\D{3}(\d+)\D{3}fsTer(\d+)\)", q)
                    if mm:
                        end = int(mm.group(1)) + int(mm.group(2)) - 1
                        if end > L0:
                            evid += f" The new reading frame runs past the normal stop codon ({gene} is {L0} aa; the new stop is codon {end}), so the stop count depends on the 3' UTR reference sequence."
            else:
                sev = "high"
        info = dict(self_k=self_k, ref_ok=ref_ok, suspect=suspect, type_conflict=type_conflict, our_wrong=our_wrong, intrinsic=intrinsic, exp_c=exp_c, exp_p=exp_p, sev=sev)
        if suspect and cat in ("notation", "c_p_mismatch"):
            evid += f" Our record keeps the literal numbering although {gene} uses a {gk:+d} legacy shift for one-letter records, so it may itself be mis-numbered."
    sugg = ""
    if cat in ("notation", "c_p_mismatch") and info.get("our_wrong") and info.get("exp_p"):
        sugg = info["exp_p"] + " (expected after fix)"
    elif cat in ("notation", "c_p_mismatch") and fault == "no" and not info.get("suspect") and len(our_p) == 1 and set(tier) <= {"A", "B"}:
        sugg = our_p[0]
    if cat == "exon_boundary":
        sugg = f"{nm}:{keys[0]}"
    info = info if info else {}
    accs = sorted({s.split("|")[1] for s in sids}); als = sorted({s.split("|")[2] for s in sids})
    iids = sorted({m[1] for m in members})
    rows.append(dict(gene=gene, accessions=";".join(accs), allele=";".join(als), individual_ids=";".join(iids), n_patients=len(iids),
                     source_g=g, source_c=c, source_p="/".join(prots_s),
                     our_c_mane=(f"{nm}:{info['exp_c']} (expected after fix)" if info.get("our_wrong") else ";".join(f"{nm}:{k}" for k in keys)),
                     our_p=(f"{info['exp_p']} (expected after fix)" if info.get("our_wrong") else ";".join(our_p)), our_track=";".join(tracks),
                     kept_record=";".join(f"{nm}:{k}" for k in keys), kept_p=";".join(our_p),
                     confidence_tier=";".join(tier), pair_labels=";".join(labels), scan_class=";".join(cls_scan),
                     category=cat, severity=info.get("sev", ""), evidence=evid, suggested_correction=sugg, pipeline_fault_possible=fault,
                     our_record_wrong="yes" if info.get("our_wrong") else "no", literal_in_shift_gene="yes" if info.get("suspect") else "no",
                     intrinsic="yes" if info.get("intrinsic") else "no", sids=" ; ".join(sorted(set(sids)))))

x = pd.DataFrame(rows).sort_values(["gene", "accessions"])
cols = ["gene", "accessions", "allele", "individual_ids", "n_patients", "source_g", "source_c", "source_p", "our_c_mane", "our_p", "our_track",
        "confidence_tier", "category", "severity", "evidence", "suggested_correction", "pipeline_fault_possible", "our_record_wrong"]
cand = x[~x.category.isin(["pipeline", "not_error"])]
keep = cand[((cand.pipeline_fault_possible == "no") | (cand.intrinsic == "yes")) & (cand.category != "untestable")]
out = []
for r in keep.itertuples():
    for s in r.sids.split(" ; "):
        gene, acc, al, _, ctok = s.split("|", 4)
        out.append({"source_change_id": s, "gene": gene, "accession": acc, "allele": al, "source_c": ctok, "source_p": r.source_p,
                    "category": r.category, "severity": r.severity, "our_record_wrong": r.our_record_wrong, "status": "candidate",
                    "evidence": r.evidence})
table = pd.DataFrame(out, columns=OUT_COLS).sort_values(["gene", "accession", "allele", "source_c"])
args.out.parent.mkdir(parents=True, exist_ok=True)
table.to_csv(args.out, sep="\t", index=False)
if args.report_dir:
    args.report_dir.mkdir(parents=True, exist_ok=True)
    x.drop(columns=["sids"]).to_csv(args.report_dir / "source_error_all.tsv", sep="\t", index=False)
    cand[cols].to_csv(args.report_dir / "source_error_candidates.csv", index=False)
    keep[cols].to_csv(args.report_dir / "source_error_for_mauno.csv", index=False)
    tick = x[x.category.isin(["pipeline", "not_error"]) | (x.our_record_wrong == "yes") | (x.literal_in_shift_gene == "yes")]
    tick[cols + ["kept_record", "kept_p", "literal_in_shift_gene"]].to_csv(args.report_dir / "pipeline_or_not_error.csv", index=False)
print(f"source changes classified: {len(x)} ({sum(x.n_patients)} pairs); candidates written: {len(keep)} changes, {len(table)} source change ids")
print(dict(Counter(keep.category)), dict(Counter(keep.severity)))
