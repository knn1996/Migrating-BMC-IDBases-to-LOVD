"""
audit_ledger.py
===============
Append one audit run to the Step10 audit history and compare it with the previous run.

Reads the outputs of audit_against_source.py (audit_variants.tsv,
audit_metadata_flags.tsv) plus the final LOVD tables of the same pipeline state,
and maintains, under the history directory:

    audit_history.tsv                 one wide row per run; re-running a run_id replaces its row in place
    snapshots/<run_id>.tsv.gz         per source unit: gene, accession, allele, fate, fate_detail, protein_check
    changes/<run_id>_vs_<prev>.tsv    units whose fate or protein_check differ from the previous run
    audit_trend.png                   rendered by Utilities/plot_audit_history.py

Source units are keyed on (gene, accession, allele) from the IDbases source, never
on remapped HGVS, so runs stay comparable when normalisation changes. A few source
records repeat an allele label (two "Allele 1" entries); the second and later
occurrences, in source order, get an ordinal suffix ("1#2"). If
source_units differs from the previous run the ruler has changed and the run is
refused: re-score the history with the current audit script first.

Final counts follow Step8_Rescue/build_final_lovd.py (its definitions are
mirrored here because that module needs its pipeline environment to import):
final_unique_variants = rows of dedup_merged_variants.rescued.tsv, final_pairs =
sum of patient_count, final_individuals = distinct individual_id (GENE_accession)
over accession_list. The wb_* columns are the lovd_review workbook definitions on
lovd_flat_with_patients_final.tsv: distinct (gene, dedup_key), distinct
Individual/ID (gene_accession), and row count.

Usage:
    python Step10_Audit/audit_ledger.py --label "after rescue rerun"
    python Step10_Audit/audit_ledger.py --label "pre-fix baseline" --audit-dir /tmp/audit0 \
        --final-tsv /tmp/audit0/lovd_flat_with_patients_final.tsv \
        --final-dedup-tsv /tmp/audit0/dedup_merged_variants.rescued.tsv --pipeline-commit 7c7f6fc
"""

import argparse
import gzip
import importlib.util
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

SCRIPTS = Path(__file__).resolve().parents[1]
AUDIT_SCRIPT = SCRIPTS / "Step10_Audit" / "audit_against_source.py"
PLOT_SCRIPT = SCRIPTS / "Utilities" / "plot_audit_history.py"
KEY = ["gene", "accession", "allele"]
SNAP_COLS = KEY + ["fate", "fate_detail", "protein_check"]
META_EXCLUDED = {"SUBHEADERS_STRIPPED", "EXTRA_REFS_NOT_CARRIED"}
HEAD_COLS = ["run_id", "label", "recorded_utc", "pipeline_commit", "outputs_dirty", "audit_script_sha", "prev_run_id"]


class LedgerError(Exception):
    pass


def load_protein_classes():
    spec = importlib.util.spec_from_file_location("audit_against_source", AUDIT_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return frozenset(mod.GOOD_P), tuple(sorted(mod.BAD_P))


def git(repo, *args):
    r = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if r.returncode != 0:
        raise LedgerError(f"git {' '.join(args)} failed: {r.stderr.strip()}")
    return r.stdout.strip()


def config_paths(config_path):
    with open(config_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    out = Path(cfg["output_root"])
    return {
        "audit_dir": out / "Step10_Audit",
        "history_dir": out / "Step10_Audit" / "history",
        "final_tsv": out / "Step8_Rescue" / "lovd_flat_with_patients_final.tsv",
        "final_dedup_tsv": out / "Step8_Rescue" / "dedup_merged_variants.rescued.tsv",
        "output_root": out,
    }


def parse_args(defaults):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--label", required=True)
    ap.add_argument("--audit-dir", type=Path, default=defaults["audit_dir"])
    ap.add_argument("--final-tsv", type=Path, default=defaults["final_tsv"])
    ap.add_argument("--final-dedup-tsv", type=Path, default=defaults["final_dedup_tsv"])
    ap.add_argument("--history-dir", type=Path, default=defaults["history_dir"])
    ap.add_argument("--pipeline-commit", default="HEAD")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--plot-python", default=sys.executable)
    return ap.parse_args()


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
    if not path.exists():
        raise LedgerError(f"missing input: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


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


def compute_metrics(var, meta, final, dedup, good_p, bad_p):
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
    return m


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


def write_atomic_tsv(df, path):
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def main():
    defaults = config_paths(Path(os.environ.get("CONFIG", SCRIPTS / "config" / "config.yaml")))
    args = parse_args(defaults)
    repo = Path(git(SCRIPTS, "rev-parse", "--show-toplevel"))
    good_p, bad_p = load_protein_classes()

    sha, dirty = resolve_commit(repo, args.pipeline_commit, defaults["output_root"], args.history_dir, args.allow_dirty)
    run_id = args.run_id or f"{datetime.now():%Y%m%d}_{sha[:7]}"
    audit_sha = git(repo, "hash-object", str(AUDIT_SCRIPT))

    var = disambiguate_alleles(read_tsv(args.audit_dir / "audit_variants.tsv"))
    meta = read_tsv(args.audit_dir / "audit_metadata_flags.tsv")
    final = read_tsv(args.final_tsv)
    dedup = read_tsv(args.final_dedup_tsv)
    snap = build_snapshot(var)
    metrics = compute_metrics(var, meta, final, dedup, good_p, bad_p)

    hist_path = args.history_dir / "audit_history.tsv"
    history = load_history(hist_path)
    prev = previous_row(history, run_id)
    prev_id = prev["run_id"] if prev is not None else ""
    if prev is not None and str(prev.get("source_units", "")) != str(metrics["source_units"]):
        raise LedgerError(
            f"ruler changed — re-score history (source_units {prev.get('source_units')} in {prev_id} "
            f"vs {metrics['source_units']} now)")

    counts = {"gains": "", "regressions": "", "protein_good_to_bad": "", "protein_bad_to_good": ""}
    changes = None
    if prev is not None:
        prev_snap_path = args.history_dir / "snapshots" / f"{prev_id}.tsv.gz"
        if not prev_snap_path.exists():
            raise LedgerError(f"previous snapshot missing: {prev_snap_path}")
        prev_snap = pd.read_csv(prev_snap_path, sep="\t", dtype=str, keep_default_na=False)
        changes, c = diff_snapshots(prev_snap, snap, good_p, bad_p)
        counts = {k: c[k] for k in counts}
        if c["keys_new"] or c["keys_dropped"]:
            print(f"WARNING: key set differs from {prev_id}: {c['keys_new']} new, {c['keys_dropped']} dropped "
                  f"(same source_units); see the changes file")

    row = {"run_id": run_id, "label": args.label,
           "recorded_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "pipeline_commit": sha, "outputs_dirty": dirty, "audit_script_sha": audit_sha,
           "prev_run_id": prev_id, **metrics, **counts}

    (args.history_dir / "snapshots").mkdir(parents=True, exist_ok=True)
    (args.history_dir / "changes").mkdir(parents=True, exist_ok=True)
    snap_path = args.history_dir / "snapshots" / f"{run_id}.tsv.gz"
    with gzip.open(snap_path, "wt", encoding="utf-8", newline="") as f:
        snap.to_csv(f, sep="\t", index=False)
    if changes is not None:
        write_atomic_tsv(changes, args.history_dir / "changes" / f"{run_id}_vs_{prev_id}.tsv")

    new_row = pd.DataFrame([{k: str(v) for k, v in row.items()}])
    if run_id in set(history["run_id"]):
        pos = list(history["run_id"]).index(run_id)
        history = pd.concat([history.iloc[:pos], new_row, history.iloc[pos + 1:]], ignore_index=True)
    else:
        history = pd.concat([history, new_row], ignore_index=True)
    ordered = HEAD_COLS + [c for c in row if c not in HEAD_COLS]
    history = history.reindex(columns=ordered + [c for c in history.columns if c not in ordered]).fillna("")
    keep = set(HEAD_COLS) | set(row)
    history = history[[c for c in history.columns if c in keep or (history[c] != "").any()]]
    write_atomic_tsv(history, hist_path)

    print(f"run_id {run_id}  ({args.label})  commit {sha[:7]}  prev {prev_id or '-'}")
    for k in ordered[len(HEAD_COLS):]:
        print(f"  {k:<48} {row.get(k, '')}")
    print(f"history : {hist_path}")
    print(f"snapshot: {snap_path}")

    r = subprocess.run([args.plot_python, str(PLOT_SCRIPT), "--history-dir", str(args.history_dir)])
    if r.returncode != 0:
        print("WARNING: audit_trend.png not updated (plot script failed; is matplotlib installed?)")


if __name__ == "__main__":
    try:
        main()
    except LedgerError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)
