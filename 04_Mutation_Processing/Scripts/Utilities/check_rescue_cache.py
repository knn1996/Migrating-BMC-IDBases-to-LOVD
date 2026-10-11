"""
check_rescue_cache.py
=====================
Read-only estimate of how many API queries the next rescue.smk run would send.

Rebuilds the rescue inputs from the current Step7 unresolved_variants.tsv in a
scratch directory, using the same scripts and environment variables as the
rescue.smk rules (classify_unresolved, generate_rescue_input). Nothing under
Output/ or the cache directories is written, and no network calls are made.

Forecasts:
    mutalyzer_rescue    exact: unique hgvs_input not in rescue_MANE.jsonl
    vv_rescue           exact: normalised hgvs_input of RESCUABLE_AUTO rows not in
                        vv_rescue_probe.jsonl
    resolve_unresolved  upper bound: uncached direct probes and +/-1 offset shifts,
                        plus ENOSELECTORFOUND anchor lookups that are either known
                        and uncached in vv_resolve.jsonl or depend on a Mutalyzer
                        rescue query that is not cached yet
Per-query timings are measured from the existing rescue logs.

Usage:
    python Utilities/check_rescue_cache.py
Environment:
    SCRIPTS_DIR  pipeline Scripts directory  (default: parent of this file's folder)
    CONFIG       path to config.yaml         (default: SCRIPTS_DIR/config/config.yaml)
    WORK_DIR     scratch directory to use    (default: a new temporary directory)
    KEEP_WORK    set to 1 to keep WORK_DIR after the run
"""

import csv
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import yaml

SCRIPTS = Path(os.environ.get("SCRIPTS_DIR", Path(__file__).resolve().parents[1]))
CONFIG = Path(os.environ.get("CONFIG", SCRIPTS / "config" / "config.yaml"))
sys.path.insert(0, str(SCRIPTS / "Step8_Rescue"))
from hgvs_norm import normalize_hgvs_query

SIMPLE_SUB = re.compile(r"^c\.(-?\d+)([ACGT]+>[ACGT]+)$", re.IGNORECASE)
VV_PROBE_SUBS = {"single_base_wrong", "single_base_complement", "multibase_rc_match",
                 "completely_different", "length_mismatch"}
SKIP_DISPOSITIONS = {"RESCUABLE_MANUAL", "UNKNOWN"}
UNRESCUABLE = {"UNRESCUABLE", "UNRESCUABLE_DATA_ERROR"}
DEFAULT_MUTALYZER_SPQ = 3.0
DEFAULT_VV_SPQ = 2.0


def run(label, script, env_extra, work):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", **{k: str(v) for k, v in env_extra.items()})
    log_path = work / "logs" / f"{label}.out"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as fh:
        rc = subprocess.run([sys.executable, str(SCRIPTS / script)], env=env, cwd=work,
                            stdout=fh, stderr=subprocess.STDOUT).returncode
    print(f"  {label:<24} {'ok' if rc == 0 else 'FAILED rc=' + str(rc)}", flush=True)
    if rc != 0:
        print(f"\nStep {label} failed; see {log_path}")
        sys.exit(1)


def read_rows(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_cache(path, key_field, value_field):
    cache = {}
    if not path.exists():
        return cache
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
                cache[obj[key_field]] = obj[value_field]
            except (json.JSONDecodeError, KeyError):
                continue
    return cache


def seconds_per_query(log_path, pattern, fmt):
    if not log_path.exists():
        return None
    stamps = []
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = pattern.search(line)
            if m:
                stamps.append((datetime.strptime(m.group(1), fmt), int(m.group(2))))
    gaps = [(b[0] - a[0]).total_seconds() for a, b in zip(stamps, stamps[1:])
            if b[1] == a[1] + 1 and 0 <= (b[0] - a[0]).total_seconds() < 120]
    return statistics.mean(gaps) if len(gaps) >= 20 else None


def fmt_dur(sec):
    sec = int(round(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    return f"{h}h{m:02d}m" if h else f"{m}m{s:02d}s"


def resolve_candidates(disp):
    cand = []
    n_anchor_rows = 0
    for r in disp:
        cat, sub, d = r.get("category", ""), r.get("subcategory", ""), r.get("disposition", "")
        if d in UNRESCUABLE and not (cat == "ESEQUENCEMISMATCH" and sub in VV_PROBE_SUBS):
            continue
        if d in SKIP_DISPOSITIONS or cat == "NO_REF":
            continue
        h = normalize_hgvs_query((r.get("hgvs_input") or "").strip())
        if cat == "ENOSELECTORFOUND":
            n_anchor_rows += 1
        elif cat == "ESEQUENCEMISMATCH" and sub == "positional_offset" and ":" in h:
            prefix, c_part = h.rsplit(":", 1)
            m = SIMPLE_SUB.match(c_part.strip())
            if m:
                for delta in (-1, 1):
                    p = int(m.group(1)) + delta
                    if p:
                        cand.append(f"{prefix}:c.{p}{m.group(2)}")
        elif cat == "ESEQUENCEMISMATCH" and sub in VV_PROBE_SUBS and h:
            cand.append(h)
    return list(dict.fromkeys(normalize_hgvs_query(c) for c in cand)), n_anchor_rows


def main():
    with open(CONFIG, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    out = Path(cfg["output_root"])
    logs = Path(cfg["logs_root"])
    mcache = Path(cfg["mutalyzer"]["cache_dir"])
    vcache = Path(cfg["variantvalidator"]["cache_dir"])
    unresolved = out / "Step7_Merging" / "unresolved_variants.tsv"
    if not unresolved.exists():
        sys.exit(f"Missing required input: {unresolved}")

    work = Path(os.environ["WORK_DIR"]) if os.environ.get("WORK_DIR") else Path(tempfile.mkdtemp(prefix="rescue_check_"))
    work.mkdir(parents=True, exist_ok=True)
    try:
        disp_path = work / "unresolved_disposition.tsv"
        resc_path = work / "rescue_mutalyzer_input.tsv"
        print(f"Rebuilding rescue inputs in {work}", flush=True)
        run("classify_unresolved", "Step8_Rescue/classify_unresolved.py",
            {"IN_PATH": unresolved, "OUT_PATH": disp_path,
             "SUMMARY": work / "unresolved_disposition_summary.tsv"}, work)
        run("generate_rescue_input", "Step8_Rescue/generate_rescue_input.py",
            {"IN_DISPOSITION": disp_path, "LRG_NM_CSV": out / "Step1b_RSG_Mapping" / "LRG_with_NM.csv",
             "MANE_NM_DIR": cfg["dna_seq"]["mane_nm"], "OUT_TSV": resc_path,
             "LOG_TSV": work / "rescue_mutalyzer_skipped.tsv"}, work)

        disp = read_rows(disp_path)
        old_disp_path = out / "Step8_Rescue" / "unresolved_disposition.tsv"
        old_disp = read_rows(old_disp_path) if old_disp_path.exists() else []
        mut_cache = load_cache(mcache / "rescue_MANE.jsonl", "key", "value")
        probe_cache = load_cache(vcache / "vv_rescue_probe.jsonl", "query", "response")
        resolve_cache = load_cache(vcache / "vv_resolve.jsonl", "query", "response")

        resc_q = list(dict.fromkeys(r["hgvs_input"] for r in read_rows(resc_path) if r.get("hgvs_input")))
        old_resc_path = out / "Step8_Rescue" / "rescue_mutalyzer_input.tsv"
        old_resc_q = list(dict.fromkeys(r["hgvs_input"] for r in read_rows(old_resc_path)
                                        if r.get("hgvs_input"))) if old_resc_path.exists() else []
        mut_miss = [q for q in resc_q if q not in mut_cache]

        vv_q = list(dict.fromkeys(normalize_hgvs_query(r["hgvs_input"].strip()) for r in disp
                                  if r.get("disposition") == "RESCUABLE_AUTO" and (r.get("hgvs_input") or "").strip()))
        vv_miss = [q for q in vv_q if q not in probe_cache]

        cand, n_anchor_rows = resolve_candidates(disp)
        res_miss = [q for q in cand if q not in resolve_cache]
        anchors_known_uncached = 0
        for q in resc_q:
            v = mut_cache.get(q)
            if not v or v.get("status") != "ok":
                continue
            a = v.get("normalized") or v.get("c_hgvs")
            if a and normalize_hgvs_query(a) not in resolve_cache:
                anchors_known_uncached += 1
        anchor_hi = min(n_anchor_rows, anchors_known_uncached + len(mut_miss))

        m_spq = seconds_per_query(logs / "mutalyzer_rescue.log",
                                  re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\s+\[(\d+)/\d+\]"),
                                  "%Y-%m-%d %H:%M:%S") or DEFAULT_MUTALYZER_SPQ
        v_spq = seconds_per_query(logs / "vv_rescue.log",
                                  re.compile(r"^\[(\d\d:\d\d:\d\d)\]\s+\[(\d+)/\d+\]"),
                                  "%H:%M:%S") or DEFAULT_VV_SPQ

        breakdown = {}
        for r in disp:
            breakdown[r.get("disposition", "")] = breakdown.get(r.get("disposition", ""), 0) + 1
        print(f"\ndisposition rows      : before {len(old_disp)}  after {len(disp)}")
        print(f"disposition breakdown : {dict(sorted(breakdown.items()))}")

        hi = len(res_miss) + anchor_hi
        hdr = f"{'call':<42}{'unique':>8}{'prev':>8}{'cached':>8}{'to send':>10}{'s/q':>7}{'est.':>15}"
        print("\n" + hdr)
        print("-" * len(hdr))
        print(f"{'mutalyzer_rescue (rescue_MANE)':<42}{len(resc_q):>8}{len(old_resc_q):>8}"
              f"{len(resc_q) - len(mut_miss):>8}{len(mut_miss):>10}{m_spq:>7.2f}{fmt_dur(len(mut_miss) * m_spq):>15}")
        print(f"{'vv_rescue (vv_rescue_probe)':<42}{len(vv_q):>8}{'-':>8}"
              f"{len(vv_q) - len(vv_miss):>8}{len(vv_miss):>10}{v_spq:>7.2f}{fmt_dur(len(vv_miss) * v_spq):>15}")
        print(f"{'resolve_unresolved VV (vv_resolve), max':<42}{len(cand) + n_anchor_rows:>8}{'-':>8}"
              f"{'-':>8}{f'0-{hi}':>10}{v_spq:>7.2f}{'0m00s-' + fmt_dur(hi * v_spq):>15}")
        print("-" * len(hdr))
        lo_t = len(mut_miss) * m_spq + len(vv_miss) * v_spq
        print(f"{'Estimated API time (sequential)':<83}{fmt_dur(lo_t) + '-' + fmt_dur(lo_t + hi * v_spq):>15}")
        print(f"\nresolve max = {len(res_miss)} uncached probe/offset-shift queries "
              f"+ {anchor_hi} anchor lookups ({anchors_known_uncached} known uncached, "
              f"up to {len(mut_miss)} from uncached Mutalyzer rescue results)")
        print("NMD GenBank fetches (annotate_nmd) are not forecast.")
    finally:
        if os.environ.get("KEEP_WORK") == "1" or os.environ.get("WORK_DIR"):
            print(f"\nScratch kept: {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
