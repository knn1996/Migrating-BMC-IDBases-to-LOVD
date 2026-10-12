"""
check_mutalyzer_cache.py
========================
Read-only estimate of how many Mutalyzer queries the next primary run would send.

Rebuilds the Step6 Mutalyzer input tables from the current Step1a extraction in a
scratch directory, using the same scripts and environment variables as the
Snakefile rules (refcheck_offset, reference_confirm_source1,
build_reference_summary, generate_bed, liftover, filter_bed_hg38,
bed_to_mutalyzer, mutalyzer_input_NM). Nothing under Output/, 03_BED_Files/ or
the cache directory is written, and no network calls are made.

Each track's unique hgvs_input set is then compared against its JSONL cache.
Track A pass 2 (c. level) depends on pass 1 responses, so for g. queries not yet
cached it can only report an upper bound. Per-query timings are measured from
the existing Mutalyzer run logs.

Usage:
    python Utilities/check_mutalyzer_cache.py
Environment:
    SCRIPTS_DIR  pipeline Scripts directory  (default: parent of this file's folder)
    CONFIG       path to config.yaml         (default: SCRIPTS_DIR/config/config.yaml)
    LEGACY_REMAP_TSV  legacy c. remap table  (default: SCRIPTS_DIR/refs/legacy_remap.tsv)
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
LEGACY_REMAP_TSV = Path(os.environ.get("LEGACY_REMAP_TSV", SCRIPTS / "refs" / "legacy_remap.tsv"))
DEFAULT_SEC_PER_QUERY = 3.0
SLEEP_FLOOR = 0.5


def load_config():
    with open(CONFIG, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    out = Path(cfg["output_root"])
    return {
        "out": out,
        "logs": Path(cfg["logs_root"]),
        "cache": Path(cfg["mutalyzer"]["cache_dir"]),
        "mut": out / "Step1a_Extraction" / "all_mutations.tsv",
        "lrg_nm": out / "Step1b_RSG_Mapping" / "LRG_with_NM.csv",
        "blast": out / "Step3_BLAST" / "blast_first_hits_source_IDRefseq.tsv",
        "s6": out / "Step6_Mutalyzer",
        "idrefseq_ng": cfg["ref_seq"]["idrefseq_ng"],
        "source1": cfg["ref_seq"]["source1_processed"],
        "mane_nm": cfg["dna_seq"]["mane_nm"],
        "idr_nm": cfg["dna_seq"]["idrefseq_nm"],
        "mane_ng": cfg["dna_seq"]["mane_ng"],
    }


def run(label, cmd, env_extra, work):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", **{k: str(v) for k, v in env_extra.items()})
    log_path = work / "logs" / f"{label}.out"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = datetime.now()
    with open(log_path, "w", encoding="utf-8") as fh:
        rc = subprocess.run(cmd, env=env, cwd=work, stdout=fh, stderr=subprocess.STDOUT).returncode
    secs = (datetime.now() - t0).total_seconds()
    print(f"  {label:<28} {'ok' if rc == 0 else 'FAILED rc=' + str(rc):<12} {secs:7.1f}s", flush=True)
    if rc != 0:
        print(f"\nStep {label} failed; see {log_path}")
        sys.exit(1)


def rebuild_inputs(p, work):
    py = sys.executable
    s2 = work / "Step2_RefCheck"
    bed_raw = work / "bed_raw"
    s5 = work / "Step5_LiftOver"
    s6 = work / "Step6_Mutalyzer"
    print(f"Rebuilding Step6 inputs in {work}", flush=True)
    run("refcheck_offset", [py, str(SCRIPTS / "Step2_RefCheck/find_lrg_offset.py"), "--source", "IDRefseq"],
        {"MUTATIONS_TSV": p["mut"], "SEQ_DIR": p["idrefseq_ng"], "OUT_CSV": s2 / "lrg_offset_results.csv",
         "LOG_PATH": work / "logs" / "find_lrg_offset_IDRefseq.log"}, work)
    run("reference_confirm_source1", [py, str(SCRIPTS / "Step2_RefCheck/reference_confirm.py"), "--source", "1"],
        {"SEQ_DIR": p["source1"], "MUTATIONS_TSV": p["mut"], "OUT_PATH": s2 / "reference_check_source_1.tsv"}, work)
    run("build_reference_summary", [py, str(SCRIPTS / "Step2_RefCheck/build_reference_summary.py")],
        {"SOURCE1_CHECK": s2 / "reference_check_source_1.tsv", "LRG_OFFSET": s2 / "lrg_offset_results.csv",
         "LRG_NM_CSV": p["lrg_nm"], "BLAST_HITS": p["blast"], "OUT_CSV": s2 / "reference_summary.csv"}, work)
    run("generate_bed", [py, str(SCRIPTS / "Step4_BED/generate_bed.py"), "--source", "IDRefseq"],
        {"MUTATIONS_TSV": p["mut"], "BLAST_HITS": p["blast"], "BED_DIR": bed_raw, "LOG_DIR": work / "logs"}, work)
    run("liftover", ["bash", str(SCRIPTS / "Step5_LiftOver/liftover.sh")],
        {"IN_DIR": bed_raw, "OUT_DIR": bed_raw / "hg38", "CHAIN_DIR": SCRIPTS / "Step5_LiftOver"}, work)
    run("filter_bed_hg38", [py, str(SCRIPTS / "Step4_BED/filter_bed.py")],
        {"REF_SUMMARY": s2 / "reference_summary.csv", "IN_DIR": bed_raw / "hg38", "OUT_DIR": s5,
         "OUT_EXCLUDED": work / "excluded_accessions_hg38.tsv"}, work)
    run("bed_to_mutalyzer", [py, str(SCRIPTS / "Step6_Mutalyzer/bed_to_mutalyzer.py")],
        {"BED_DIR": s5, "OUT_TSV": s6 / "mutalyzer_input.tsv", "LOG_TSV": s6 / "mutalyzer_skipped.tsv"}, work)
    for source in ("MANE", "IDRefseq"):
        run(f"mutalyzer_input_NM_{source}",
            [py, str(SCRIPTS / "Step6_Mutalyzer/generate_mutalyzer_input_NM.py"), "--source", source],
            {"ALL_MUTATIONS_TSV": p["mut"], "BED_DIR": s5, "MANE_NM_DIR": p["mane_nm"],
             "IDREFSEQ_NM_DIR": p["idr_nm"], "MANE_NG_DIR": p["mane_ng"],
             "OUT_TSV": s6 / f"mutalyzer_input_NM_{source}.tsv",
             "LOG_TSV": s6 / f"mutalyzer_input_NM_skipped_{source}.tsv", "LEGACY_REMAP_TSV": LEGACY_REMAP_TSV}, work)
    return s6


def unique_queries(tsv):
    if not tsv.exists():
        return []
    with open(tsv, encoding="utf-8") as f:
        return list(dict.fromkeys(r["hgvs_input"] for r in csv.DictReader(f, delimiter="\t") if r.get("hgvs_input")))


def load_cache(path):
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
            except json.JSONDecodeError:
                continue
            cache[obj["key"]] = obj["value"]
    return cache


def seconds_per_query(log_path, marker):
    if not log_path.exists():
        return None
    pat = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\s+" + re.escape(marker) + r"\[(\d+)/(\d+)\]")
    stamps = []
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = pat.match(line.strip())
            if m:
                stamps.append((datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S"), int(m.group(2))))
    gaps = [(b[0] - a[0]).total_seconds() for a, b in zip(stamps, stamps[1:])
            if b[1] == a[1] + 1 and 0 <= (b[0] - a[0]).total_seconds() < 120]
    return statistics.mean(gaps) if len(gaps) >= 20 else None


def fmt_dur(sec):
    sec = int(round(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    return f"{h}h{m:02d}m" if h else f"{m}m{s:02d}s"


def main():
    p = load_config()
    for key in ("mut", "lrg_nm", "blast"):
        if not p[key].exists():
            sys.exit(f"Missing required input: {p[key]}")

    work = Path(os.environ["WORK_DIR"]) if os.environ.get("WORK_DIR") else Path(tempfile.mkdtemp(prefix="mutalyzer_check_"))
    work.mkdir(parents=True, exist_ok=True)
    try:
        s6_new = rebuild_inputs(p, work)

        cache_ag = load_cache(p["cache"] / "track_a_g.jsonl")
        cache_ac = load_cache(p["cache"] / "track_a_c.jsonl")
        cache_b = load_cache(p["cache"] / "track_b_MANE.jsonl")
        cache_c = load_cache(p["cache"] / "track_c_IDRefseq.jsonl")

        new_a = unique_queries(s6_new / "mutalyzer_input.tsv")
        old_a = unique_queries(p["s6"] / "mutalyzer_input.tsv")
        miss_ag = [h for h in new_a if h not in cache_ag]
        c_from_cached = list(dict.fromkeys(cache_ag[h]["c_hgvs"] for h in new_a
                                           if h in cache_ag and cache_ag[h].get("c_hgvs")))
        miss_ac_known = [c for c in c_from_cached if c not in cache_ac]

        tracks = [
            ("A g.  (track_a_g)", new_a, old_a, cache_ag, len(miss_ag), len(miss_ag), "mutalyzer_run.log", "P1 "),
            ("A c.  (track_a_c)", c_from_cached, None, cache_ac, len(miss_ac_known),
             len(miss_ac_known) + len(miss_ag), "mutalyzer_run.log", "P2 "),
        ]
        for source, cache, stem in (("MANE", cache_b, "track_b_MANE"), ("IDRefseq", cache_c, "track_c_IDRefseq")):
            new_q = unique_queries(s6_new / f"mutalyzer_input_NM_{source}.tsv")
            old_q = unique_queries(p["s6"] / f"mutalyzer_input_NM_{source}.tsv")
            miss = sum(1 for h in new_q if h not in cache)
            label = f"{'B' if source == 'MANE' else 'C'} {source:<9}({stem})"
            tracks.append((label, new_q, old_q, cache, miss, miss, f"mutalyzer_run_NM_{source}.log", ""))

        print("\nPer-track query forecast (cache is read-only)")
        hdr = f"{'track':<32}{'unique':>8}{'prev':>8}{'cached':>8}{'to send':>10}{'s/query':>9}{'est.':>12}"
        print(hdr)
        print("-" * len(hdr))
        total_lo = total_hi = 0.0
        notes = []
        for label, qs, old, cache, lo, hi, log_name, marker in tracks:
            spq = seconds_per_query(p["logs"] / log_name, marker)
            if spq is None and log_name.startswith("mutalyzer_run_NM_"):
                spq = seconds_per_query(p["logs"] / "mutalyzer_run_NM.log", marker)
            src = "log"
            if spq is None:
                spq, src = DEFAULT_SEC_PER_QUERY, "default"
            spq = max(spq, SLEEP_FLOOR)
            to_send = str(lo) if lo == hi else f"{lo}-{hi}"
            est = fmt_dur(lo * spq) if lo == hi else f"{fmt_dur(lo * spq)}-{fmt_dur(hi * spq)}"
            prev = "-" if old is None else str(len(old))
            flag = "*" if src == "default" else " "
            print(f"{label:<32}{len(qs):>8}{prev:>8}{len(qs) - lo:>8}{to_send:>10}{spq:>8.2f}{flag}{est:>12}")
            total_lo += lo * spq
            total_hi += hi * spq
            bad = sum(1 for h in qs if h in cache and str(cache[h].get("status", "")) == "unreachable")
            if bad:
                notes.append(f"{label.strip()}: {bad} cached entries have status 'unreachable' and will NOT be retried")
        print("-" * len(hdr))
        total = fmt_dur(total_lo) if total_lo == total_hi else f"{fmt_dur(total_lo)}-{fmt_dur(total_hi)}"
        print(f"{'Estimated Mutalyzer time (sequential)':<75}{total:>12}")
        print("\nunique = distinct hgvs_input the next run would submit; prev = same count in the current Output inputs")
        print("A c. lower bound counts c. strings derivable from cached g. responses; the upper bound adds one per uncached g. query")
        print("* = no usable timings in logs, default used")
        for n in notes:
            print("NOTE:", n)
        print("Rescue track (rescue_MANE) is not forecast: its inputs depend on Step7 results of this run.")
    finally:
        if os.environ.get("KEEP_WORK") == "1" or os.environ.get("WORK_DIR"):
            print(f"\nScratch kept: {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
