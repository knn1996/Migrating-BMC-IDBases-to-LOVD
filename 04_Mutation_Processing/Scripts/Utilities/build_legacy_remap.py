"""
build_legacy_remap.py
=====================
Derives the legacy c. numbering remap table from the current pipeline outputs.

For every inconsistent final pair (Step10_Audit/audit_pairs.tsv) the pair's own source
changes are tested with legacy_remap.reconcile on the gene's MANE Select CDS at every
shift k in [-K, K], k != 0. A gene's shift is chosen from the changes that reconcile:

    candidates     shifts supported by at least 2 independent source changes
    choice         the candidate with the most substitution supports, then the most
                   supports; a tie on that order drops the gene
    evidence       a gene whose chosen shift has exactly 2 supports keeps it only if at
                   least one support is a substitution

Every source change that reconciles at its gene's chosen shift is written to
legacy_remap.tsv for each accession whose inconsistent pair carries it. The per-gene
evidence (every supporting change and its evidence type, the decision and its reason)
goes to legacy_shift_evidence.tsv.

--all-changes [GENE ...] tests every source change of the named genes (all genes when
none are named) from all_mutations.tsv, whether or not its pair is inconsistent; a change
that reconciles at k = 0 is neither a support nor remapped. --base-dir keeps the rows of
an existing table for the genes that are not rebuilt. --manual adds reviewed decisions
(gene, accession, source_c, action shift|exclude, shift, reason); an exclude row has an
empty remapped_c and drops the change from the c. track inputs.

MANE CDS: GenBank cache record of the MANE Select NM (MANE v1.5 GFF); otherwise the
MANE NM FASTA with the CDS start from coding_offset_scan_MANE.csv, else the longest ORF.

Usage:
    python Utilities/build_legacy_remap.py --out-dir <dir> [--exclude JAK3 C1QB CYBB ADA] [--max-shift 1500]
    python Utilities/build_legacy_remap.py --out-dir <dir> --all-changes ADA --base-dir refs --manual <tsv>
"""

import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS / "Utilities"))
sys.path.insert(0, str(SCRIPTS / "Step10_Audit"))
sys.path.insert(0, str(SCRIPTS / "Step7_Merging"))
import legacy_remap as R
import audit_ledger as L
from source_change import SourceChanges

MANE_GFF = Path("01_Reference_Genomes") / "GRCh38_MANE" / "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz"
EVIDENCE_COLS = ["gene", "mane_nm", "cds_source", "status", "reason", "shift", "n_support", "n_substitution",
                 "supports", "other_candidates"]


def mane_select(gff):
    out = {}
    with gzip.open(gff, "rt", encoding="utf-8") as f:
        for line in f:
            p = line.split("\t")
            if len(p) > 8 and p[2] == "transcript" and "tag=MANE_Select" in p[8]:
                a = dict(x.split("=", 1) for x in p[8].strip().split(";") if "=" in x)
                nm = re.search(r"RefSeq:(NM_\d+\.\d+)", a.get("Dbxref", ""))
                if nm:
                    out[a["gene_name"]] = nm.group(1)
    return out


def genbank_cds(path):
    text = path.read_text(encoding="utf-8")
    m = re.search(r"^ {5}CDS {13}<?(\d+)\.\.>?(\d+)", text, re.M)
    if not m or "\nORIGIN" not in text:
        return None
    seq = re.sub(r"[^acgtnACGTN]", "", text.split("\nORIGIN", 1)[1].split("\n//")[0]).upper()
    return seq, int(m.group(1)) - 1


def longest_orf(seq):
    best = (0, 0)
    for s in range(len(seq) - 2):
        if seq[s:s + 3] == "ATG":
            i = s
            while i + 3 <= len(seq) and R.CODON.get(seq[i:i + 3]) != "*":
                i += 3
            if i - s > best[0]:
                best = (i - s, s)
    return best[1]


def reference(gene, nm, genbank_dir, mane_nm_dir, scan):
    gb = Path(genbank_dir) / f"{nm}.gb"
    if gb.exists():
        x = genbank_cds(gb)
        if x:
            return x[0], x[1], "genbank"
    fa = Path(mane_nm_dir) / f"{gene}_{nm}.fasta"
    if not fa.exists():
        return None
    seq = "".join(l.strip() for l in open(fa, encoding="utf-8") if not l.startswith(">")).upper()
    start = scan.get(gene, "")
    if start:
        return seq, int(float(start)) - 1, "fasta+scan"
    return seq, longest_orf(seq), "fasta+longest_orf"


def choose(support):
    cands = {k: v for k, v in support.items() if len(v) >= 2}
    if not cands:
        return None, "no shift supported by 2 or more changes", []
    order = sorted(cands, key=lambda k: (-sum(e == "substitution" for _, e in cands[k]), -len(cands[k]), abs(k)))
    top = order[0]
    key = lambda k: (sum(e == "substitution" for _, e in cands[k]), len(cands[k]))
    ties = [k for k in order[1:] if key(k) == key(top)]
    if ties:
        return None, f"tie between {', '.join(f'{k:+d}' for k in [top] + ties)} not pinned by substitutions", order
    n_sub = sum(e == "substitution" for _, e in cands[top])
    if len(cands[top]) == 2 and n_sub == 0:
        return None, f"{top:+d} has 2 supports and no substitution", order
    return top, "", order


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--max-shift", type=int, default=1500)
    ap.add_argument("--all-changes", nargs="*", default=None, metavar="GENE")
    ap.add_argument("--base-dir", type=Path, default=None)
    ap.add_argument("--manual", type=Path, default=None)
    args = ap.parse_args()
    with open(os.environ.get("CONFIG", SCRIPTS / "config" / "config.yaml"), encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    out_root, thesis = Path(cfg["output_root"]), Path(cfg["thesis_dir"])
    mane_nm_dir = cfg["dna_seq"]["mane_nm"]
    excluded = {g.upper() for g in args.exclude}

    pairs = L.read_tsv(out_root / "Step10_Audit" / "audit_pairs.tsv")
    var = L.disambiguate_alleles(L.read_tsv(out_root / "Step10_Audit" / "audit_variants.tsv"))
    proteins = L.SourceProteins(var)
    dedup = L.read_tsv(out_root / "Step8_Rescue" / "dedup_merged_variants.resolved.tsv")
    symbol = {}
    for g, x in dedup.groupby("gene"):
        named = x.mutalyzer_gene[x.mutalyzer_gene != ""]
        symbol[g.upper()] = named.mode().iloc[0] if len(named) else g.upper()
    mane = mane_select(thesis / MANE_GFF)
    scan_csv = out_root / "Step2_RefCheck" / "coding_offset_scan_MANE.csv"
    scan = {}
    if scan_csv.exists():
        with open(scan_csv, encoding="utf-8") as f:
            scan = {r["gene"].upper(): r.get("cds_start", "") for r in csv.DictReader(f)}

    changes = defaultdict(set)
    all_mode = args.all_changes is not None
    scope = {g.upper() for g in args.all_changes} if args.all_changes else None
    if all_mode:
        source = SourceChanges(str(out_root / "Step1a_Extraction" / "all_mutations.tsv"))
        for key, (sid, _) in source.by_coding.items():
            gene, acc, _, _, ctok = sid.split("|", 4)
            if ctok != "-" and (scope is None or gene in scope):
                changes[(gene, ctok, tuple(proteins.of(sid)))].add(acc)
    else:
        for iid, sids, label in zip(pairs.individual_id, pairs.source_change_ids, pairs.label):
            if not label.startswith("inconsistent:"):
                continue
            for sid in (s.strip() for s in sids.split(" ; ") if s.strip()):
                gene, acc, _, _, ctok = sid.split("|", 4)
                if ctok != "-":
                    changes[(gene, ctok, tuple(proteins.of(sid)))].add(acc)

    by_gene = defaultdict(list)
    for key in changes:
        by_gene[key[0]].append(key)

    evidence_rows, remap_rows = [], []
    for gene in sorted(by_gene):
        nm = mane.get(symbol.get(gene, gene), "")
        row = {"gene": gene, "mane_nm": nm, "cds_source": "", "status": "", "reason": "", "shift": "",
               "n_support": "", "n_substitution": "", "supports": "", "other_candidates": ""}
        if gene in excluded:
            row.update(status="excluded", reason="excluded by decision")
            evidence_rows.append(row)
            continue
        ref = reference(gene, nm, out_root / "cache" / "genbank", mane_nm_dir, scan) if nm else None
        if ref is None:
            row.update(status="dropped", reason="no MANE CDS available offline")
            evidence_rows.append(row)
            continue
        seq, cds0, source = ref
        row["cds_source"] = source
        hits = {}
        support = defaultdict(list)
        for key in by_gene[gene]:
            change = R.parse_c(key[1])
            prots = [p for p in (R.parse_p(x) for x in key[2]) if p]
            if not change or not prots:
                continue
            found = {}
            if all_mode and any(R.reconcile(seq, cds0, change, p, 0) for p in prots):
                hits[key] = found
                continue
            for k in range(-args.max_shift, args.max_shift + 1):
                if k == 0:
                    continue
                for p in prots:
                    ev = R.reconcile(seq, cds0, change, p, k)
                    if ev:
                        found[k] = ev
                        break
            hits[key] = found
            for k, ev in found.items():
                support[k].append((key, ev))
        shift, reason, order = choose(support)
        others = [f"{k:+d}({len(support[k])})" for k in order if k != shift][:6]
        row["other_candidates"] = " ".join(others)
        if shift is None:
            row.update(status="dropped", reason=reason)
            evidence_rows.append(row)
            continue
        sup = support[shift]
        row.update(status="kept", shift=f"{shift:+d}", n_support=len(sup),
                   n_substitution=sum(e == "substitution" for _, e in sup),
                   supports=" ; ".join(f"{k[1]} {'/'.join(k[2])} [{e}]" for k, e in sorted(sup)))
        evidence_rows.append(row)
        for key, ev in sorted(sup):
            new = R.shift_notation(key[1], shift)
            if not new:
                continue
            for acc in sorted(changes[key]):
                remap_rows.append({"gene": gene, "accession": acc, "source_c": key[1], "remapped_c": new,
                                   "shift": f"{shift:+d}", "evidence": ev, "mane_nm": nm})

    if args.manual:
        rebuilt = {r["gene"] for r in evidence_rows}
        with open(args.manual, encoding="utf-8") as f:
            for m in csv.DictReader(f, delimiter="\t"):
                gene = m["gene"].strip().upper()
                nm = mane.get(symbol.get(gene, gene), "")
                remap_rows = [r for r in remap_rows if (r["gene"], r["accession"], r["source_c"]) != (gene, m["accession"], m["source_c"])]
                if m["action"].strip() == "exclude":
                    new, shift = "", ""
                else:
                    shift = f"{int(m['shift']):+d}"
                    new = R.shift_notation(m["source_c"], int(m["shift"])) or ""
                remap_rows.append({"gene": gene, "accession": m["accession"].strip(), "source_c": m["source_c"].strip(),
                                   "remapped_c": new, "shift": shift, "evidence": f"{m['action'].strip()}: {m['reason'].strip()}",
                                   "mane_nm": nm})
                rebuilt.add(gene)
    if args.base_dir:
        rebuilt = {r["gene"] for r in evidence_rows} | {r["gene"] for r in remap_rows}
        base_remap = pd.read_csv(args.base_dir / "legacy_remap.tsv", sep="\t", dtype=str, keep_default_na=False)
        base_ev = pd.read_csv(args.base_dir / "legacy_shift_evidence.tsv", sep="\t", dtype=str, keep_default_na=False)
        remap_rows = base_remap[~base_remap.gene.isin(rebuilt)].to_dict("records") + remap_rows
        evidence_rows = base_ev[~base_ev.gene.isin(rebuilt)].to_dict("records") + evidence_rows
        evidence_rows.sort(key=lambda r: r["gene"])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(evidence_rows, columns=EVIDENCE_COLS).to_csv(args.out_dir / "legacy_shift_evidence.tsv", sep="\t", index=False)
    pd.DataFrame(remap_rows, columns=R.REMAP_COLS).to_csv(args.out_dir / "legacy_remap.tsv", sep="\t", index=False)
    kept = [r for r in evidence_rows if r["status"] == "kept"]
    print(f"genes with inconsistent pairs: {len(evidence_rows)}; kept {len(kept)}; "
          f"dropped {sum(r['status'] == 'dropped' for r in evidence_rows)}; excluded {sum(r['status'] == 'excluded' for r in evidence_rows)}")
    print(f"remap rows: {len(remap_rows)}")


if __name__ == "__main__":
    main()
