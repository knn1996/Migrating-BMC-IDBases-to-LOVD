"""
audit_ledger.py
===============
Append one audit run to the Step10 audit history and compare it with the previous run,
or re-score every recorded run with the current definitions (--rescore-all).

Reads the outputs of audit_against_source.py (audit_variants.tsv,
audit_metadata_flags.tsv) plus the final LOVD tables of the same pipeline state,
and maintains, under the history directory:

    audit_history.tsv                 one wide row per run; re-running a run_id replaces its row in place
    snapshots/<run_id>.tsv.gz         per source unit: gene, accession, allele, fate, fate_detail, protein_check
    snapshots/<run_id>.pairs.tsv.gz   per final (individual, variant) pair: individual_id, variant_key, track,
                                      label, equivalence, source_change_ids
    changes/<run_id>_vs_<prev>.tsv    units whose fate or protein_check differ from the previous run, and
                                      pairs that became inconsistent, became consistent, or were removed
                                      while inconsistent
    audit_trend.png                   rendered by Utilities/plot_audit_history.py

audit_pairs.tsv (the pair table of the run) is also written next to audit_variants.tsv
in --audit-dir.

Source units are keyed on (gene, accession, allele) from the IDbases source, never
on remapped HGVS, so runs stay comparable when normalisation changes. A few source
records repeat an allele label (two "Allele 1" entries); the second and later
occurrences, in source order, get an ordinal suffix ("1#2"). If
source_units differs from the previous run the ruler has changed and the run is
refused: re-score the history with the current audit script first.

Final counts follow Step8_Rescue/build_final_lovd.py (its definitions are
mirrored here because that module needs its pipeline environment to import):
final_unique_variants = rows of the final dedup table, final_pairs = sum of
patient_count, final_individuals = distinct individual_id (GENE_accession) over
accession_list. The final dedup table is dedup_merged_variants.resolved.tsv
(after resolve_tracks) where the pipeline state has it, and
dedup_merged_variants.rescued.tsv for states before track resolution. The wb_*
columns are the lovd_review workbook definitions on
lovd_flat_with_patients_final.tsv: distinct (gene, dedup_key), distinct
Individual/ID (gene_accession), and row count.

Precision metrics work on the final pairs, i.e. distinct (individual_id, _dedup_key)
in lovd_flat_with_patients_final.tsv with individual_id = GENE_accession:

    excess        expected = distinct source DNA changes over all alleles of the
                  individual (each in-cis component counts, the same change on both
                  alleles counts once; genomic notation is used when the allele has it,
                  coding otherwise). individuals_with_excess counts individuals with
                  more final pairs than expected, excess_pairs sums the surplus.
    link          every pair is linked to its own source changes. The pipeline state's
                  own dedup_merged_variants.py (git show at the run commit) is replayed
                  on its augmented_merged_variants.tsv; the replay must reproduce
                  dedup_merged_variants.rescued.tsv byte for byte. Each augmented row is
                  followed from its stage-1 record to the record that holds its unit
                  after the XM/NM merge (position merge before dd183da, logged stage-2
                  rules after), and gets its source_change_id from
                  Step7_Merging/source_change.py on the run's extraction. Where the
                  state stores ids (augmented source_change_id, resolution_pairs.tsv)
                  the rebuilt ids must equal them; a pair may differ only when its
                  stored link is empty and every rebuilt id came through a stage-2 merge.
    consistency   the final protein of a pair is compared with the source proteins of
                  its own source changes: the p. component at the index of the change's
                  c. token (g. token if the c. token is absent) in the allele notation,
                  where the p. list has one entry per token. A slash-joined or repeated
                  p. component counts as alternatives. consistent = some comparison is in
                  GOOD_P (compare_protein) or is one of the notation equivalences;
                  inconsistent:<label> = none good and at least one in BAD_P (the label
                  of the closest one); uncheckable:no_final_p, :unlinked, :source_error
                  (every source change of the pair is in SOURCE_ERRORS), :no_source_p or
                  :not_comparable otherwise. pair_consistency_pct = consistent /
                  (consistent + inconsistent).
    equivalences  E1 immediate stop: source XnfsX1, XnfsX (no distance) or XnfsX<n> with
                  n = start (absolute stop numbering) equals final p.(XnTer).
                  E2 insertion as dup: an ins whose right flank is the dup start, with the
                  dup length and the dup end residues as first and last inserted residues.

--rescore-all re-runs audit_against_source.py and the pair linkage for every recorded
run on inputs pulled from the run's pipeline_commit with git show into a scratch
directory, and recomputes the row. Every pre-existing column outside the precision
columns (pairs_*, pair_consistency_pct) must reproduce exactly; if any value moves
nothing is written. label, recorded_utc, pipeline_commit and outputs_dirty are kept
from the recorded row. all_mutations.tsv is not tracked, so the current file is used
unless --extract-tsv RUN_ID=PATH names another; it feeds stage_extracted and the
source change ids of the linkage.

Usage:
    python Step10_Audit/audit_ledger.py --label "after rescue rerun"
    python Step10_Audit/audit_ledger.py --rescore-all \
        --extract-tsv 20261011_7c7f6fc=99_Archive/04_Mutation_Processing/Output/Step1_Extraction/all_mutations_pre_fix.tsv
"""

import argparse
import contextlib
import gzip
import importlib.util
import io
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

SCRIPTS = Path(__file__).resolve().parents[1]
AUDIT_SCRIPT = SCRIPTS / "Step10_Audit" / "audit_against_source.py"
DEDUP_SCRIPT = SCRIPTS / "Step7_Merging" / "dedup_merged_variants.py"
SOURCE_CHANGE_MODULE = SCRIPTS / "Step7_Merging" / "source_change.py"
PLOT_SCRIPT = SCRIPTS / "Utilities" / "plot_audit_history.py"
KEY = ["gene", "accession", "allele"]
SNAP_COLS = KEY + ["fate", "fate_detail", "protein_check"]
PAIR_COLS = ["individual_id", "variant_key", "track", "label", "equivalence", "source_change_ids"]
META_EXCLUDED = {"SUBHEADERS_STRIPPED", "EXTRA_REFS_NOT_CARRIED"}
HEAD_COLS = ["run_id", "label", "recorded_utc", "rescored_utc", "pipeline_commit", "outputs_dirty",
             "audit_script_sha", "prev_run_id"]
KEPT_ON_RESCORE = ["label", "recorded_utc", "pipeline_commit", "outputs_dirty"]
NOT_VALIDATED = {"rescored_utc"}
FINAL_DEDUP_CANDIDATES = ("Step8_Rescue/dedup_merged_variants.resolved.tsv",
                          "Step8_Rescue/dedup_merged_variants.rescued.tsv")
RUN_INPUTS = {
    "FINAL_TSV": ("Step8_Rescue/lovd_flat_with_patients_final.tsv",),
    "AUGMENT_TSV": ("Step8_Rescue/augmented_merged_variants.tsv",),
    "UNRES_TSV": ("Step8_Rescue/unresolved_disposition.tsv",),
    "RESCUED_TSV": ("Step8_Rescue/resolved_unresolved.tsv",),
    "DEDUP_TSV": FINAL_DEDUP_CANDIDATES,
    "RESCUED_DEDUP_TSV": ("Step8_Rescue/dedup_merged_variants.rescued.tsv",),
    "RESOLUTION_PAIRS_TSV": ("Step8_Rescue/resolution_pairs.tsv",),
}
REQUIRED_INPUTS = {"FINAL_TSV", "AUGMENT_TSV", "UNRES_TSV", "DEDUP_TSV", "RESCUED_DEDUP_TSV"}
RE_SEG = re.compile(r"^[\[\s]*(?:([gcrp])\s*\.|(?=IVS))", re.I)
COUNT_KEYS = ["gains", "regressions", "protein_good_to_bad", "protein_bad_to_good",
              "pairs_became_inconsistent", "pairs_became_consistent", "pairs_inconsistent_removed"]
UNCHECKABLE = ["no_final_p", "unlinked", "source_error", "no_source_p", "not_comparable"]
EQUIVALENCES = ["E1", "E2"]
SOURCE_ERRORS = frozenset()
DEDUP_CAPTURE = (
    ("result = result.loc[kept]", "result_all = result.copy()\nresult = result.loc[kept]"),
    ("result = result.drop(index=drop_idx)", "result_all = result.copy()\nresult = result.drop(index=drop_idx)"),
)
AA1 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
       "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
       "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "X": "Ter", "*": "Ter"}
RES = r"(?:[A-Z][a-z]{2}|[A-Z*])"
RE_SRC_FS = re.compile(rf"^p\.\(?({RES})(\d+){RES}?fs(?:X|\*|Ter)?(\d*)\)?$")
RE_FIN_STOP = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)(?:Ter|\*|X)\)?$")
RE_INS = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)ins((?:[A-Z][a-z]{2})+)\)?$")
RE_DUP = re.compile(r"^p\.\(?([A-Z][a-z]{2})(\d+)(?:_([A-Z][a-z]{2})(\d+))?dup\)?$")


class LedgerError(Exception):
    pass


class Audit:
    def __init__(self):
        spec = importlib.util.spec_from_file_location("audit_against_source", AUDIT_SCRIPT)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        self.good_p = frozenset(mod.GOOD_P)
        self.bad_p = tuple(sorted(mod.BAD_P))
        self.compare = mod.compare_protein
        spec = importlib.util.spec_from_file_location("source_change", SOURCE_CHANGE_MODULE)
        sc = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(sc)
        self.source_changes = sc.SourceChanges


def git(repo, *args):
    r = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if r.returncode != 0:
        raise LedgerError(f"git {' '.join(args)} failed: {r.stderr.strip()}")
    return r.stdout.strip()


def final_dedup_path(output_root):
    for rel in FINAL_DEDUP_CANDIDATES:
        if (output_root / rel).exists():
            return output_root / rel
    return output_root / FINAL_DEDUP_CANDIDATES[-1]


def config_paths(config_path):
    with open(config_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    out = Path(cfg["output_root"])
    return {
        "audit_dir": out / "Step10_Audit",
        "history_dir": out / "Step10_Audit" / "history",
        "final_tsv": out / "Step8_Rescue" / "lovd_flat_with_patients_final.tsv",
        "final_dedup_tsv": final_dedup_path(out),
        "augment_tsv": out / "Step8_Rescue" / "augmented_merged_variants.tsv",
        "rescued_dedup_tsv": out / "Step8_Rescue" / "dedup_merged_variants.rescued.tsv",
        "resolution_pairs_tsv": out / "Step8_Rescue" / "resolution_pairs.tsv",
        "extract_tsv": out / "Step1a_Extraction" / "all_mutations.tsv",
        "thesis_dir": Path(cfg["thesis_dir"]),
        "output_root": out,
    }


def parse_args(defaults):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--label", default=None)
    ap.add_argument("--audit-dir", type=Path, default=defaults["audit_dir"])
    ap.add_argument("--final-tsv", type=Path, default=defaults["final_tsv"])
    ap.add_argument("--final-dedup-tsv", type=Path, default=defaults["final_dedup_tsv"])
    ap.add_argument("--augment-tsv", type=Path, default=defaults["augment_tsv"])
    ap.add_argument("--rescued-dedup-tsv", type=Path, default=defaults["rescued_dedup_tsv"])
    ap.add_argument("--resolution-pairs-tsv", type=Path, default=defaults["resolution_pairs_tsv"])
    ap.add_argument("--history-dir", type=Path, default=defaults["history_dir"])
    ap.add_argument("--pipeline-commit", default="HEAD")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--plot-python", default=sys.executable)
    ap.add_argument("--rescore-all", action="store_true")
    ap.add_argument("--scratch-dir", type=Path, default=None)
    ap.add_argument("--extract-tsv", action="append", default=[], metavar="RUN_ID=PATH")
    args = ap.parse_args()
    if not args.rescore_all and not args.label:
        ap.error("--label is required unless --rescore-all is given")
    return args


def resolve_commit(repo, ref, output_root, history_dir, allow_dirty):
    sha = git(repo, "rev-parse", "--verify", f"{ref}^{{commit}}")
    head = git(repo, "rev-parse", "HEAD")
    if sha != head:
        return sha, "n/a"
    rel_out = os.path.relpath(output_root, repo)
    rel_hist = os.path.relpath(Path(history_dir).resolve(), repo)
    spec = ["--", rel_out]
    if not rel_hist.startswith(".."):
        spec.append(f":(exclude){rel_hist}")
    dirty = git(repo, "status", "--porcelain", *spec)
    if dirty and not allow_dirty:
        lines = dirty.splitlines()
        shown = "\n  ".join(lines[:15]) + ("\n  ..." if len(lines) > 15 else "")
        raise LedgerError(
            f"pipeline outputs under {rel_out} are not committed ({len(lines)} paths):\n  {shown}\n"
            f"Commit them first so the run maps to a commit, or pass --allow-dirty.")
    return sha, "yes" if dirty else "no"


def read_tsv(path):
    if not Path(path).exists():
        raise LedgerError(f"missing input: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def text(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v).strip()


def disambiguate_alleles(var):
    var = var.copy()
    nth = var.groupby(KEY, sort=False).cumcount() + 1
    var["allele"] = var["allele"].where(nth == 1, var["allele"] + "#" + nth.astype(str))
    return var


def build_snapshot(var):
    missing = [c for c in SNAP_COLS if c not in var.columns]
    if missing:
        raise LedgerError(f"audit_variants.tsv lacks columns {missing}")
    snap = var[SNAP_COLS].copy()
    dup = snap.duplicated(KEY, keep=False)
    if dup.any():
        sample = snap.loc[dup, KEY].head(5).to_dict("records")
        raise LedgerError(f"(gene, accession, allele) is not unique in audit_variants.tsv: {int(dup.sum())} rows, e.g. {sample}")
    return snap.sort_values(KEY).reset_index(drop=True)


def top_level_segments(txt):
    out, depth, cur = [], 0, []
    for ch in txt:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        if ch == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return [s.strip() for s in out if s.strip()]


def source_components(txt):
    lv = defaultdict(list)
    for seg in top_level_segments(txt):
        m = RE_SEG.match(seg)
        if not m:
            continue
        level = (m.group(1) or "g").lower()
        for comp in re.split(r";|\s\+\s", seg[m.end():]):
            comp = re.sub(r"^[gcrp]\.", "", re.sub(r"[\[\]\s]", "", comp), flags=re.I)
            if comp:
                lv[level].append(comp)
    return lv


def source_dna_changes(txt):
    lv = source_components(txt)
    for level in ("g", "c"):
        if lv[level]:
            return {(level, comp) for comp in lv[level]}
    return set()


def protein_alternatives(comp):
    parts = [x for x in re.split(r"/|p\.", comp) if x and x != "?"]
    return list(dict.fromkeys("p." + x for x in parts))


def individual_ids(genes, accessions):
    return genes.str.strip().str.upper() + "_" + accessions.str.strip()


def notation_token(tok):
    return re.sub(r"^[gc]\.", "", re.sub(r"[\[\]\s]", "", tok), flags=re.I)


class SourceProteins:
    def __init__(self, var):
        self.alleles = defaultdict(list)
        for iid, txt in zip(individual_ids(var["gene"], var["accession"]), var["source_systematic_name"]):
            self.alleles[iid].append(source_components(txt))

    def of(self, sid):
        parts = sid.split("|", 4)
        if len(parts) != 5:
            return []
        gene, acc, _, g_tok, c_tok = parts
        out = []
        for lv in self.alleles.get(f"{gene}_{acc}", []):
            for level, tok in (("c", c_tok), ("g", g_tok)):
                if tok == "-":
                    continue
                toks, n = lv[level], notation_token(tok)
                if n not in toks:
                    continue
                if len(lv["p"]) == len(toks):
                    for i, x in enumerate(toks):
                        if x == n:
                            out.extend(protein_alternatives(lv["p"][i]))
                break
        return list(dict.fromkeys(out))


def residue(x):
    return AA1.get(x, x) if len(x) == 1 else x


def immediate_stop(src, fin):
    s, f = RE_SRC_FS.match(src), RE_FIN_STOP.match(fin)
    if not s or not f:
        return False
    pos = int(s.group(2))
    if residue(s.group(1)) != f.group(1) or pos != int(f.group(2)):
        return False
    return s.group(3) in ("", "1") or int(s.group(3)) == pos


def ins_as_dup(ins, dup):
    a, d = RE_INS.match(ins), RE_DUP.match(dup)
    if not a or not d:
        return False
    if int(a.group(4)) != int(a.group(2)) + 1:
        return False
    seq = re.findall(r"[A-Z][a-z]{2}", a.group(5))
    start_res, start = d.group(1), int(d.group(2))
    end_res, end = (d.group(3), int(d.group(4))) if d.group(3) else (start_res, start)
    return (int(a.group(4)) == start and a.group(3) == start_res and len(seq) == end - start + 1
            and seq[0] == start_res and seq[-1] == end_res)


def equivalence(src, fin):
    if immediate_stop(src, fin):
        return "E1"
    if ins_as_dup(src, fin) or ins_as_dup(fin, src):
        return "E2"
    return ""


def pair_label(final_p, sids, src_proteins, audit):
    if not final_p.strip():
        return "uncheckable:no_final_p", ""
    if not sids:
        return "uncheckable:unlinked", ""
    if all(s in SOURCE_ERRORS for s in sids):
        return "uncheckable:source_error", ""
    src_ps = list(dict.fromkeys(p for s in sids if s not in SOURCE_ERRORS for p in src_proteins.of(s)))
    if not src_ps:
        return "uncheckable:no_source_p", ""
    res = [audit.compare(s, final_p) for s in src_ps]
    if any(st in audit.good_p for st, _ in res):
        return "consistent", ""
    for s in src_ps:
        eq = equivalence(s, final_p.strip())
        if eq:
            return "consistent", eq
    bad = sorted((abs(off) if off is not None else float("inf"), st) for st, off in res if st in audit.bad_p)
    if bad:
        return f"inconsistent:{bad[0][1]}", ""
    return "uncheckable:not_comparable", ""


def replay_dedup(code, augment_tsv, rescued_dedup_tsv, work):
    hits = [c for c in DEDUP_CAPTURE if code.count(c[0]) == 1]
    if len(hits) != 1:
        raise LedgerError("dedup replay: capture anchor not found exactly once in dedup_merged_variants.py")
    work.mkdir(parents=True, exist_ok=True)
    env = {"INPUT_TSV": str(augment_tsv), "OUTPUT_TSV": str(work / "dedup_replay.tsv"),
           "STATS_TXT": str(work / "dedup_replay_stats.txt"), "DEDUP_LOG": str(work / "dedup_replay_log.tsv")}
    saved = {k: os.environ.get(k) for k in env}
    os.environ.update(env)
    ns = {"__name__": "dedup_replay"}
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(code.replace(*hits[0]), "dedup_merged_variants.py", "exec"), ns)
    except SystemExit as e:
        raise LedgerError(f"dedup replay exited: {e}")
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    if Path(env["OUTPUT_TSV"]).read_bytes() != Path(rescued_dedup_tsv).read_bytes():
        raise LedgerError(f"dedup replay does not reproduce {rescued_dedup_tsv}; see {env['OUTPUT_TSV']}")
    return ns


def row_targets(ns):
    df, ra, res = ns["df"], ns["result_all"], ns["result"]
    kept = set(res.index)
    keyidx = {(text(ra.at[i, "gene"]), text(ra.at[i, "dedup_key"])): i for i in ra.index}
    out = {}
    if "units" in ns:
        units, describe = ns["units"], ns["describe"]
        desc = {(text(ra.at[i, "gene"]), describe(ra.loc[i])): i for i in ra.index}
        ident, moved = {}, defaultdict(list)
        for lr in ns["log_rows"]:
            g = lr["gene"]
            if lr["rule"] == "identical_change":
                ident[desc[(g, lr["dropped_record"])]] = desc[(g, lr["kept_record"])]
            elif lr["rule"].startswith("same_source_unit"):
                n = desc[(g, lr["kept_record"])]
                for u in lr["source_units"].split(";"):
                    a, al = u.split(":", 1)
                    moved[(g, a, "" if al == "?" else al)].append(n)
        for r in df.index:
            g = text(df.at[r, "gene"])
            u = (text(df.at[r, "accession"]), text(df.at[r, "allele_num"]))
            i = keyidx.get((g, text(df.at[r, "_dedup_key"])))
            if i is None:
                out[r] = []
            elif i in kept and u in units[i]:
                out[r] = [(i, "stage1")]
            elif i in ident:
                out[r] = [(ident[i], "stage2")]
            else:
                out[r] = [(n, "stage2") for n in dict.fromkeys(moved.get((g,) + u, [])) if n in kept and u in units[n]]
    else:
        best = {}
        located = ra[ra["_gpos"].notna()]
        for _, grp in located.groupby("_gid"):
            nm, xm = grp[~grp["_xm"]], grp[grp["_xm"]]
            if len(nm) and len(xm):
                b = nm["patient_count"].idxmax()
                for x in xm.index:
                    best[x] = b
        accs = {i: set(text(res.at[i, "accession_list"]).split(";")) for i in res.index}
        for r in df.index:
            a = text(df.at[r, "accession"])
            i = keyidx.get((text(df.at[r, "gene"]), text(df.at[r, "_dedup_key"])))
            if i is None:
                out[r] = []
            elif i in kept:
                out[r] = [(i, "stage1")] if a in accs[i] else []
            elif i in best and a in accs[best[i]]:
                out[r] = [(best[i], "stage2")]
            else:
                out[r] = []
    return out


def build_links(link_inputs, audit):
    ns = replay_dedup(link_inputs["dedup_code"], link_inputs["augment_tsv"], link_inputs["rescued_dedup_tsv"],
                      link_inputs["work"])
    df, ra = ns["df"], ns["result_all"]
    changes = audit.source_changes(str(link_inputs["extract_tsv"]))
    has_stored = "source_change_id" in df.columns
    links = defaultdict(dict)
    stats = {"rows": len(df), "rows_unmapped": 0, "rows_without_source_change": 0}
    row_mismatch = []
    for r, targets in row_targets(ns).items():
        gene, acc = text(df.at[r, "gene"]).upper(), text(df.at[r, "accession"])
        sid, _ = changes.lookup(gene, acc, text(df.at[r, "allele_num"]), text(df.at[r, "sysname"]),
                                text(df.at[r, "source_track"]))
        if has_stored and text(df.at[r, "source_change_id"]) != sid:
            row_mismatch.append((gene, acc, text(df.at[r, "sysname"]), text(df.at[r, "source_change_id"]), sid))
        if not targets:
            stats["rows_unmapped"] += 1
        if not sid:
            stats["rows_without_source_change"] += 1
        for i, path in targets:
            d = links[(f"{gene}_{acc}", text(ra.at[i, "dedup_key"]))]
            if sid and (sid not in d or path == "stage1"):
                d[sid] = path
    if row_mismatch:
        shown = "\n  ".join(map(str, row_mismatch[:10]))
        raise LedgerError(f"rebuilt source_change_id differs from the stored one on {len(row_mismatch)} "
                          f"augmented rows:\n  {shown}")
    return dict(links), stats


def link_gate(links, pair_keys, resolution_pairs_tsv):
    if not resolution_pairs_tsv or not Path(resolution_pairs_tsv).exists():
        return ""
    rp = read_tsv(resolution_pairs_tsv)
    stored = {(i, k): {s.strip() for s in v.split(" ; ") if s.strip()}
              for i, k, v in zip(rp["individual_id"], rp["dedup_key"], rp["source_change_ids"])}
    bad, exceptions = [], 0
    for pk in pair_keys:
        rebuilt = links.get(pk, {})
        if pk not in stored:
            bad.append((pk, "pair absent from resolution_pairs", sorted(rebuilt), None))
            continue
        if set(rebuilt) == stored[pk]:
            continue
        if not stored[pk] and rebuilt and all(p == "stage2" for p in rebuilt.values()):
            exceptions += 1
            continue
        bad.append((pk, "link differs", sorted(rebuilt.items()), sorted(stored[pk])))
    if bad:
        shown = "\n  ".join(map(str, bad[:15]))
        raise LedgerError(f"link gate failed on {len(bad)} pairs (rebuilt vs resolution_pairs):\n  {shown}")
    return exceptions


def build_pairs(var, final, audit, links):
    need = ["_gene", "_accession", "_dedup_key", "_source_track", "VariantOnTranscript/Protein"]
    missing = [c for c in need if c not in final.columns]
    if missing:
        raise LedgerError(f"final table lacks columns {missing}")
    src_proteins = SourceProteins(var)
    f = final.assign(individual_id=individual_ids(final["_gene"], final["_accession"]))
    f = f.drop_duplicates(["individual_id", "_dedup_key"])
    rows = []
    for iid, key, track, p in zip(f["individual_id"], f["_dedup_key"], f["_source_track"],
                                  f["VariantOnTranscript/Protein"]):
        sids = sorted(links.get((iid, key), {}))
        label, eq = pair_label(p, sids, src_proteins, audit)
        rows.append({"individual_id": iid, "variant_key": key, "track": track, "label": label,
                     "equivalence": eq, "source_change_ids": " ; ".join(sids)})
    pairs = pd.DataFrame(rows, columns=PAIR_COLS)
    return pairs.sort_values(["individual_id", "variant_key"]).reset_index(drop=True)


def excess_metrics(var, pairs):
    expected = defaultdict(set)
    for iid, txt in zip(individual_ids(var["gene"], var["accession"]), var["source_systematic_name"]):
        expected[iid] |= source_dna_changes(txt)
    observed = pairs.groupby("individual_id").size()
    surplus = [n - len(expected.get(i, ())) for i, n in observed.items() if n > len(expected.get(i, ()))]
    return {"individuals_with_excess": len(surplus), "excess_pairs": int(sum(surplus))}


def pair_metrics(pairs, bad_p, link_stats):
    m = {}
    lab = pairs["label"]
    cons = int((lab == "consistent").sum())
    inc = lab.str.startswith("inconsistent:")
    m["pairs_checkable"] = cons + int(inc.sum())
    m["pairs_consistent"] = cons
    m["pair_consistency_pct"] = round(100 * cons / m["pairs_checkable"], 2) if m["pairs_checkable"] else 0.0
    m["pairs_inconsistent"] = int(inc.sum())
    for label in bad_p:
        m[f"pairs_inconsistent__{label}"] = int((lab == f"inconsistent:{label}").sum())
    m["pairs_uncheckable"] = int(lab.str.startswith("uncheckable:").sum())
    for reason in UNCHECKABLE:
        m[f"pairs_uncheckable__{reason}"] = int((lab == f"uncheckable:{reason}").sum())
    for eq in EQUIVALENCES:
        m[f"pairs_equivalent__{eq}"] = int((pairs["equivalence"] == eq).sum())
    linked = int((pairs["source_change_ids"] != "").sum())
    m["pairs_linked"] = linked
    m["pairs_unlinked"] = len(pairs) - linked
    m["pairs_link_gate_stage2_exceptions"] = link_stats["gate_exceptions"]
    return m


def compute_metrics(var, meta, final, dedup, pairs, audit, link_stats):
    good_p, bad_p = audit.good_p, audit.bad_p
    m = {}
    fate = var["fate"]
    m["source_units"] = len(var)
    m["out_of_scope"] = int((fate == "OUT_OF_SCOPE").sum())
    m["in_scope"] = m["source_units"] - m["out_of_scope"]
    m["in_final"] = int((fate == "IN_FINAL").sum())
    m["in_final_pct"] = round(100 * m["in_final"] / m["in_scope"], 2) if m["in_scope"] else 0.0
    m["known_unresolved"] = int((fate == "KNOWN_UNRESOLVED").sum())
    nr = var.loc[fate == "NOT_REGISTERED", "fate_detail"]
    m["not_registered"] = len(nr)
    for detail, n in sorted(nr.value_counts().items()):
        m[f"not_registered__{detail}"] = int(n)
    pc = var["protein_check"]
    good = int(pc.isin(good_p).sum())
    bad = int(pc.isin(bad_p).sum())
    m["protein_checkable"] = good + bad
    m["protein_good"] = good
    m["protein_good_pct"] = round(100 * good / (good + bad), 2) if good + bad else 0.0
    m["protein_bad"] = bad
    for label in bad_p:
        m[f"protein_bad__{label}"] = int((pc == label).sum())
    m["metadata_flags"] = int((~meta["status"].isin(META_EXCLUDED)).sum()) if len(meta) else 0
    m["final_unique_variants"] = len(dedup)
    m["final_individuals"] = len({f"{g.strip().upper()}_{a.strip()}" for g, s in zip(dedup["gene"], dedup["accession_list"])
                                  for a in s.split(";") if a.strip()})
    m["final_pairs"] = int(pd.to_numeric(dedup["patient_count"], errors="coerce").fillna(0).sum())
    m["wb_distinct_gene_key"] = len(final.drop_duplicates(["_gene", "_dedup_key"]))
    m["wb_individuals_gene_accession"] = int(final["Individual/ID"].nunique())
    m["wb_rows"] = len(final)
    m.update(excess_metrics(var, pairs))
    m.update(pair_metrics(pairs, bad_p, link_stats))
    if len(pairs) != m["final_pairs"]:
        print(f"WARNING: {len(pairs)} distinct (individual, variant) pairs in the final table "
              f"vs final_pairs {m['final_pairs']} from the dedup table")
    return m


def pair_class(label):
    if label == "consistent":
        return "consistent"
    if label.startswith("inconsistent:"):
        return "inconsistent"
    if label.startswith("uncheckable:"):
        return "uncheckable"
    return "absent"


def diff_pairs(prev, cur):
    keep = ["individual_id", "variant_key", "track", "label"]
    j = prev[keep].merge(cur[keep], on=["individual_id", "variant_key"], how="outer", suffixes=("_prev", "_now"))
    j["label_prev"] = j["label_prev"].fillna("")
    j["label_now"] = j["label_now"].fillna("")
    cp = j["label_prev"].map(pair_class)
    cn = j["label_now"].map(pair_class)
    change = pd.Series("", index=j.index)
    change[(cn == "inconsistent") & (cp != "inconsistent")] = "pair_became_inconsistent"
    change[(cn == "consistent") & cp.isin(["inconsistent", "uncheckable"])] = "pair_became_consistent"
    change[(cp == "inconsistent") & (cn == "absent")] = "pair_inconsistent_removed"
    j["change"] = change
    j = j[j["change"] != ""].copy()
    parts = j["individual_id"].str.split("_", n=1)
    out = pd.DataFrame({
        "gene": parts.str[0], "accession": parts.str[1], "allele": "", "change": j["change"],
        "individual_id": j["individual_id"], "variant_key": j["variant_key"],
        "track_prev": j["track_prev"].fillna(""), "track_now": j["track_now"].fillna(""),
        "pair_label_prev": j["label_prev"].replace("", "absent"),
        "pair_label_now": j["label_now"].replace("", "absent"),
    })
    counts = {k: int((out["change"] == f"pair_{k[6:]}").sum())
              for k in ("pairs_became_inconsistent", "pairs_became_consistent", "pairs_inconsistent_removed")}
    return out.sort_values(["individual_id", "variant_key"]).reset_index(drop=True), counts


def diff_snapshots(prev, cur, good_p, bad_p):
    j = prev.merge(cur, on=KEY, how="outer", suffixes=("_prev", "_now"), indicator=True)
    for c in ("fate", "fate_detail", "protein_check"):
        j[f"{c}_prev"] = j[f"{c}_prev"].fillna("")
        j[f"{c}_now"] = j[f"{c}_now"].fillna("")
    in_prev = j["fate_prev"] == "IN_FINAL"
    in_now = j["fate_now"] == "IN_FINAL"
    gp, bp = j["protein_check_prev"].isin(good_p), j["protein_check_prev"].isin(bad_p)
    gn, bn = j["protein_check_now"].isin(good_p), j["protein_check_now"].isin(bad_p)
    kinds = []
    for i in range(len(j)):
        k = []
        if j["_merge"].iat[i] == "left_only":
            k.append("key_dropped")
        elif j["_merge"].iat[i] == "right_only":
            k.append("key_new")
        else:
            if in_now.iat[i] and not in_prev.iat[i]:
                k.append("gain")
            elif in_prev.iat[i] and not in_now.iat[i]:
                k.append("regression")
            elif j["fate_prev"].iat[i] != j["fate_now"].iat[i] or j["fate_detail_prev"].iat[i] != j["fate_detail_now"].iat[i]:
                k.append("fate_detail_change")
            if gp.iat[i] and bn.iat[i]:
                k.append("protein_good_to_bad")
            elif bp.iat[i] and gn.iat[i]:
                k.append("protein_bad_to_good")
            elif j["protein_check_prev"].iat[i] != j["protein_check_now"].iat[i]:
                k.append("protein_check_change")
        kinds.append(";".join(k))
    j["change"] = kinds
    changes = j[j["change"] != ""].drop(columns="_merge")
    cols = KEY + ["change", "fate_prev", "fate_now", "fate_detail_prev", "fate_detail_now",
                  "protein_check_prev", "protein_check_now"]
    changes = changes[cols].sort_values(KEY).reset_index(drop=True)
    counts = {
        "gains": int(changes["change"].str.contains(r"\bgain\b").sum()),
        "regressions": int(changes["change"].str.contains(r"\bregression\b").sum()),
        "protein_good_to_bad": int(changes["change"].str.contains("protein_good_to_bad").sum()),
        "protein_bad_to_good": int(changes["change"].str.contains("protein_bad_to_good").sum()),
        "keys_new": int(changes["change"].str.contains("key_new").sum()),
        "keys_dropped": int(changes["change"].str.contains("key_dropped").sum()),
    }
    return changes, counts


def compare_runs(prev_id, prev_snap, prev_pairs, snap, pairs, audit):
    unit_changes, c = diff_snapshots(prev_snap, snap, audit.good_p, audit.bad_p)
    if c["keys_new"] or c["keys_dropped"]:
        print(f"WARNING: key set differs from {prev_id}: {c['keys_new']} new, {c['keys_dropped']} dropped "
              f"(same source_units); see the changes file")
    counts = {k: c[k] for k in COUNT_KEYS[:4]}
    if prev_pairs is None:
        counts.update({k: "" for k in COUNT_KEYS[4:]})
        return unit_changes, counts
    pair_changes, pc = diff_pairs(prev_pairs, pairs)
    counts.update(pc)
    return pd.concat([unit_changes, pair_changes], ignore_index=True).fillna(""), counts


def evaluate(audit_dir, final_tsv, dedup_tsv, audit, link_inputs):
    var = disambiguate_alleles(read_tsv(Path(audit_dir) / "audit_variants.tsv"))
    meta = read_tsv(Path(audit_dir) / "audit_metadata_flags.tsv")
    final = read_tsv(final_tsv)
    dedup = read_tsv(dedup_tsv)
    snap = build_snapshot(var)
    links, link_stats = build_links(link_inputs, audit)
    pair_keys = set(zip(individual_ids(final["_gene"], final["_accession"]), final["_dedup_key"]))
    link_stats["gate_exceptions"] = link_gate(links, pair_keys, link_inputs.get("resolution_pairs_tsv"))
    pairs = build_pairs(var, final, audit, links)
    metrics = compute_metrics(var, meta, final, dedup, pairs, audit, link_stats)
    return snap, pairs, metrics, link_stats


def load_history(path):
    if not path.exists():
        return pd.DataFrame(columns=HEAD_COLS)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def previous_row(history, run_id):
    ids = list(history["run_id"])
    if run_id in ids:
        pos = ids.index(run_id)
        return history.iloc[pos - 1] if pos > 0 else None
    return history.iloc[-1] if len(history) else None


def read_gz_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_atomic_tsv(df, path):
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def write_gz_tsv(df, path):
    text_out = df.to_csv(sep="\t", index=False)
    if path.exists():
        with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
            if f.read() == text_out:
                return
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(text_out.encode("utf-8"))
    os.replace(tmp, path)


def make_row(run_id, head, metrics, counts):
    return {"run_id": run_id, **head, **metrics, **counts}


def merge_history(history, rows):
    for row in rows:
        new_row = pd.DataFrame([{k: str(v) for k, v in row.items()}])
        ids = list(history["run_id"])
        if row["run_id"] in ids:
            pos = ids.index(row["run_id"])
            history = pd.concat([history.iloc[:pos], new_row, history.iloc[pos + 1:]], ignore_index=True)
        else:
            history = pd.concat([history, new_row], ignore_index=True)
    ordered = HEAD_COLS + [c for c in rows[-1] if c not in HEAD_COLS]
    history = history.reindex(columns=ordered + [c for c in history.columns if c not in ordered]).fillna("")
    keep = set(HEAD_COLS) | set(rows[-1])
    return history[[c for c in history.columns if c in keep or (history[c] != "").any()]]


def write_run(history_dir, run_id, prev_id, snap, pairs, changes):
    (history_dir / "snapshots").mkdir(parents=True, exist_ok=True)
    (history_dir / "changes").mkdir(parents=True, exist_ok=True)
    write_gz_tsv(snap, history_dir / "snapshots" / f"{run_id}.tsv.gz")
    write_gz_tsv(pairs, history_dir / "snapshots" / f"{run_id}.pairs.tsv.gz")
    if changes is not None:
        write_atomic_tsv(changes, history_dir / "changes" / f"{run_id}_vs_{prev_id}.tsv")


def print_row(row, label):
    print(f"run_id {row['run_id']}  ({label})  commit {row['pipeline_commit'][:7]}  prev {row['prev_run_id'] or '-'}")
    for k in [c for c in row if c not in HEAD_COLS]:
        print(f"  {k:<48} {row.get(k, '')}")


def print_coverage(entries):
    print("pair linkage coverage")
    print(f"  {'run_id':<20} {'pairs':>6} {'linked':>7} {'unlinked':>8} {'linked%':>8} {'no_source_p':>11} "
          f"{'stage2_exc':>10} {'rows_unmapped':>13} {'rows_no_sid':>11}")
    pct = []
    for run_id, m, ls in entries:
        n = m["pairs_linked"] + m["pairs_unlinked"]
        p = 100 * m["pairs_linked"] / n if n else 0.0
        pct.append(p)
        print(f"  {run_id:<20} {n:>6} {m['pairs_linked']:>7} {m['pairs_unlinked']:>8} {p:>7.2f}% "
              f"{m['pairs_uncheckable__no_source_p']:>11} {str(ls['gate_exceptions']):>10} "
              f"{ls['rows_unmapped']:>13} {ls['rows_without_source_change']:>11}")
    med = sorted(pct)[len(pct) // 2]
    for (run_id, _, _), p in zip(entries, pct):
        if abs(p - med) > 1.0:
            print(f"  FLAG {run_id}: linked {p:.2f}% vs median {med:.2f}%")


def render_plot(plot_python, history_dir):
    r = subprocess.run([plot_python, str(PLOT_SCRIPT), "--history-dir", str(history_dir)])
    if r.returncode != 0:
        print("WARNING: audit_trend.png not updated (plot script failed; is matplotlib installed?)")


def record_one(args, defaults, repo, audit):
    sha, dirty = resolve_commit(repo, args.pipeline_commit, defaults["output_root"], args.history_dir, args.allow_dirty)
    run_id = args.run_id or f"{datetime.now():%Y%m%d}_{sha[:7]}"
    audit_sha = git(repo, "hash-object", str(AUDIT_SCRIPT))
    print(f"final dedup table: {args.final_dedup_tsv}")
    overrides = extract_overrides(args.extract_tsv, repo)
    link_inputs = {
        "dedup_code": DEDUP_SCRIPT.read_text(encoding="utf-8"),
        "augment_tsv": args.augment_tsv, "rescued_dedup_tsv": args.rescued_dedup_tsv,
        "resolution_pairs_tsv": args.resolution_pairs_tsv,
        "extract_tsv": overrides.get(run_id, defaults["extract_tsv"]),
        "work": Path(tempfile.mkdtemp(prefix="audit_link_")),
    }
    snap, pairs, metrics, link_stats = evaluate(args.audit_dir, args.final_tsv, args.final_dedup_tsv, audit, link_inputs)

    hist_path = args.history_dir / "audit_history.tsv"
    history = load_history(hist_path)
    prev = previous_row(history, run_id)
    prev_id = prev["run_id"] if prev is not None else ""
    if prev is not None and str(prev.get("source_units", "")) != str(metrics["source_units"]):
        raise LedgerError(
            f"ruler changed — re-score history (source_units {prev.get('source_units')} in {prev_id} "
            f"vs {metrics['source_units']} now)")

    counts = {k: "" for k in COUNT_KEYS}
    changes = None
    if prev is not None:
        prev_snap_path = args.history_dir / "snapshots" / f"{prev_id}.tsv.gz"
        if not prev_snap_path.exists():
            raise LedgerError(f"previous snapshot missing: {prev_snap_path}")
        prev_pairs_path = args.history_dir / "snapshots" / f"{prev_id}.pairs.tsv.gz"
        prev_pairs = read_gz_tsv(prev_pairs_path) if prev_pairs_path.exists() else None
        if prev_pairs is None:
            print(f"WARNING: {prev_pairs_path.name} missing; pair changes not compared (run --rescore-all)")
        changes, counts = compare_runs(prev_id, read_gz_tsv(prev_snap_path), prev_pairs, snap, pairs, audit)

    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    head = {"label": args.label, "recorded_utc": now, "rescored_utc": "", "pipeline_commit": sha,
            "outputs_dirty": dirty, "audit_script_sha": audit_sha, "prev_run_id": prev_id}
    row = make_row(run_id, head, metrics, counts)

    write_atomic_tsv(pairs, Path(args.audit_dir) / "audit_pairs.tsv")
    write_run(args.history_dir, run_id, prev_id, snap, pairs, changes)
    write_atomic_tsv(merge_history(history, [row]), hist_path)

    print_row(row, args.label)
    print_coverage([(run_id, metrics, link_stats)])
    print(f"history : {hist_path}")
    print(f"snapshot: {args.history_dir / 'snapshots' / f'{run_id}.tsv.gz'}")
    render_plot(args.plot_python, args.history_dir)


def stage_inputs(repo, rel_out, sha, dest):
    dest.mkdir(parents=True, exist_ok=True)
    paths = {}
    for var, candidates in RUN_INPUTS.items():
        found = None
        for rel in candidates:
            r = subprocess.run(["git", "-C", str(repo), "show", f"{sha}:{rel_out}/{rel}"], capture_output=True)
            if r.returncode == 0:
                found = (rel, r.stdout)
                break
        if found is None:
            if var in REQUIRED_INPUTS:
                raise LedgerError(f"none of {', '.join(candidates)} found at {sha[:7]}")
            paths[var] = dest / f"absent_{Path(candidates[0]).name}"
            continue
        target = dest / Path(found[0]).name
        target.write_bytes(found[1])
        paths[var] = target
    return paths


def run_audit(defaults, inputs, extract_tsv, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if k != "AUDIT_GENES"}
    env.update({
        "THESIS_DIR": str(defaults["thesis_dir"]),
        "PUB_DIR": str(defaults["thesis_dir"] / "02_Source_Database" / "idbase_pub_html_only"),
        "OUT_ROOT": str(defaults["output_root"]),
        "FINAL_TSV": str(inputs["FINAL_TSV"]),
        "AUGMENT_TSV": str(inputs["AUGMENT_TSV"]),
        "UNRES_TSV": str(inputs["UNRES_TSV"]),
        "RESCUED_TSV": str(inputs["RESCUED_TSV"]),
        "EXTRACT_TSV": str(extract_tsv),
        "OUT_DIR": str(out_dir),
        "PYTHONDONTWRITEBYTECODE": "1",
    })
    log = out_dir / "audit.log"
    with open(log, "w", encoding="utf-8") as fh:
        r = subprocess.run([sys.executable, str(AUDIT_SCRIPT)], env=env, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise LedgerError(f"audit_against_source.py failed (exit {r.returncode}); see {log}")


def extract_overrides(items, repo):
    out = {}
    for item in items:
        if "=" not in item:
            raise LedgerError(f"--extract-tsv expects RUN_ID=PATH, got {item!r}")
        run_id, path = item.split("=", 1)
        p = Path(path)
        p = p if p.is_absolute() else repo / p
        if not p.exists():
            raise LedgerError(f"--extract-tsv {run_id}: {p} not found")
        out[run_id] = p
    return out


def is_precision(col):
    return col.startswith("pairs_") or col == "pair_consistency_pct"


def validate(old, new):
    moved, precision = [], []
    for col, old_val in old.items():
        if col in NOT_VALIDATED:
            continue
        new_val = str(new.get(col, ""))
        if str(old_val) != new_val:
            (precision if is_precision(col) else moved).append((col, old_val, new_val))
    return moved, precision


def rescore_all(args, defaults, repo, audit):
    hist_path = args.history_dir / "audit_history.tsv"
    history = load_history(hist_path)
    if history.empty:
        raise LedgerError(f"no recorded runs in {hist_path}")
    rel_out = os.path.relpath(defaults["output_root"], repo)
    rel_dedup = os.path.relpath(DEDUP_SCRIPT, repo)
    overrides = extract_overrides(args.extract_tsv, repo)
    unknown = set(overrides) - set(history["run_id"])
    if unknown:
        raise LedgerError(f"--extract-tsv names runs not in the history: {sorted(unknown)}")
    scratch = args.scratch_dir or Path(tempfile.mkdtemp(prefix="audit_rescore_"))
    audit_sha = git(repo, "hash-object", str(AUDIT_SCRIPT))
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    results, failures = [], []
    prev = None
    for old in history.to_dict("records"):
        run_id = old["run_id"]
        sha = git(repo, "rev-parse", "--verify", f"{old['pipeline_commit']}^{{commit}}")
        work = scratch / run_id
        inputs = stage_inputs(repo, rel_out, sha, work / "in")
        extract = overrides.get(run_id, defaults["extract_tsv"])
        run_audit(defaults, inputs, extract, work / "audit")
        link_inputs = {
            "dedup_code": subprocess.run(["git", "-C", str(repo), "show", f"{sha}:{rel_dedup}"],
                                         capture_output=True, check=True).stdout.decode("utf-8"),
            "augment_tsv": inputs["AUGMENT_TSV"], "rescued_dedup_tsv": inputs["RESCUED_DEDUP_TSV"],
            "resolution_pairs_tsv": inputs["RESOLUTION_PAIRS_TSV"], "extract_tsv": extract,
            "work": work / "link",
        }
        snap, pairs, metrics, link_stats = evaluate(work / "audit", inputs["FINAL_TSV"], inputs["DEDUP_TSV"],
                                                    audit, link_inputs)
        write_atomic_tsv(pairs, work / "audit" / "audit_pairs.tsv")

        prev_id = prev["run_id"] if prev is not None else ""
        if prev is not None and str(prev["metrics"]["source_units"]) != str(metrics["source_units"]):
            raise LedgerError(f"ruler changed between {prev_id} and {run_id}")
        changes, counts = None, {k: "" for k in COUNT_KEYS}
        if prev is not None:
            changes, counts = compare_runs(prev_id, prev["snap"], prev["pairs"], snap, pairs, audit)

        head = {k: old.get(k, "") for k in KEPT_ON_RESCORE}
        head.update({"rescored_utc": now, "audit_script_sha": audit_sha, "prev_run_id": prev_id})
        head = {k: head[k] for k in HEAD_COLS[1:]}
        row = make_row(run_id, head, metrics, counts)
        moved, precision = validate(old, row)
        if moved:
            failures.append((run_id, moved))
        cur = {"run_id": run_id, "row": row, "snap": snap, "pairs": pairs, "changes": changes,
               "prev_id": prev_id, "metrics": metrics, "extract": extract, "link_stats": link_stats,
               "precision": precision}
        results.append(cur)
        prev = cur
        print(f"rescored {run_id} on {sha[:7]} (extract {extract.name}, dedup {inputs['DEDUP_TSV'].name}); "
              f"{len(moved)} non-precision values moved, {len(precision)} precision values moved")

    if failures:
        lines = []
        for run_id, moved in failures:
            for col, a, b in moved:
                lines.append(f"  {run_id}  {col}: recorded {a!r} -> rescored {b!r}")
        raise LedgerError("rescore changed pre-existing non-precision values; nothing was written:\n"
                          + "\n".join(lines) + f"\nscratch kept at {scratch}")

    for r in results:
        write_run(args.history_dir, r["run_id"], r["prev_id"], r["snap"], r["pairs"], r["changes"])
    write_atomic_tsv(merge_history(history, [r["row"] for r in results]), hist_path)
    print(f"validation: every pre-existing non-precision column reproduced exactly for {len(results)} runs")
    for r in results:
        print(f"precision moves {r['run_id']}:")
        for col, a, b in r["precision"]:
            print(f"  {col:<48} {a} -> {b}")
    new_cols = [c for c in results[-1]["row"] if c not in history.columns]
    print(f"new columns: {', '.join(new_cols) if new_cols else 'none'}")
    print_coverage([(r["run_id"], r["metrics"], r["link_stats"]) for r in results])
    print(f"history : {hist_path}")
    print(f"scratch : {scratch}")
    render_plot(args.plot_python, args.history_dir)


def main():
    defaults = config_paths(Path(os.environ.get("CONFIG", SCRIPTS / "config" / "config.yaml")))
    args = parse_args(defaults)
    repo = Path(git(SCRIPTS, "rev-parse", "--show-toplevel"))
    audit = Audit()
    if args.rescore_all:
        rescore_all(args, defaults, repo, audit)
    else:
        record_one(args, defaults, repo, audit)


if __name__ == "__main__":
    try:
        main()
    except LedgerError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)
