"""
plot_audit_history.py
=====================
Render the Step10 audit trend from history/audit_history.tsv.

Top panel: in-scope recovery (in_final_pct) and protein-good rate
(protein_good_pct) per run. Bottom panel: gains and regressions versus the
previous run (regressions drawn below zero). Runs are plotted in history order
and labelled with run_id and label.

Usage:
    python Utilities/plot_audit_history.py
    python Utilities/plot_audit_history.py --history-dir <dir> --out <png>
"""

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml

SCRIPTS = Path(__file__).resolve().parents[1]
LINE_COLORS = {"in_final_pct": "#2a6fb0", "protein_good_pct": "#c0662b"}
BAR_COLORS = {"gains": "#3a8f5c", "regressions": "#b8433a"}


def default_history_dir():
    with open(os.environ.get("CONFIG", SCRIPTS / "config" / "config.yaml"), encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return Path(cfg["output_root"]) / "Step10_Audit" / "history"


def num(series):
    return pd.to_numeric(series, errors="coerce")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--history-dir", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    hist_dir = args.history_dir or default_history_dir()
    out = args.out or hist_dir / "audit_trend.png"

    h = pd.read_csv(hist_dir / "audit_history.tsv", sep="\t", dtype=str, keep_default_na=False)
    if h.empty:
        raise SystemExit("audit_history.tsv has no runs")
    x = list(range(len(h)))
    ticks = [f"{r}\n{l}" for r, l in zip(h["run_id"], h["label"])]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(6.5, 1.6 * len(h) + 3), 6.5), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2]})
    for col, name in (("in_final_pct", "In-scope recovery (%)"), ("protein_good_pct", "Protein check good (%)")):
        y = num(h[col])
        ax1.plot(x, y, marker="o", color=LINE_COLORS[col], label=name, linewidth=2)
        for xi, yi in zip(x, y):
            if pd.notna(yi):
                ax1.annotate(f"{yi:.1f}", (xi, yi), textcoords="offset points", xytext=(0, 7),
                             ha="center", fontsize=8, color=LINE_COLORS[col])
    lo = min(num(h["in_final_pct"]).min(), num(h["protein_good_pct"]).min())
    ax1.set_ylim(max(0, lo - 10), 100)
    ax1.set_ylabel("%")
    ax1.grid(axis="y", alpha=0.3)
    ax1.legend(loc="lower right", frameon=False, fontsize=9)
    ax1.set_title("IDbases → LOVD audit history", fontsize=11)

    gains = num(h["gains"]).fillna(0)
    regs = num(h["regressions"]).fillna(0)
    w = 0.35
    ax2.bar([i - w / 2 for i in x], gains, width=w, color=BAR_COLORS["gains"], label="Gains vs previous")
    ax2.bar([i + w / 2 for i in x], -regs, width=w, color=BAR_COLORS["regressions"], label="Regressions vs previous")
    for i, (g, r) in enumerate(zip(gains, regs)):
        if g:
            ax2.annotate(f"+{int(g)}", (i - w / 2, g), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8)
        if r:
            ax2.annotate(f"-{int(r)}", (i + w / 2, -r), textcoords="offset points", xytext=(0, -10), ha="center", fontsize=8)
    ax2.axhline(0, color="#555555", linewidth=0.8)
    ax2.set_ylabel("Source units")
    ax2.grid(axis="y", alpha=0.3)
    ax2.legend(loc="upper left", frameon=False, fontsize=9)
    if not (gains.any() or regs.any()):
        ax2.set_ylim(-1, 1)
        ax2.text(0.5, 0.78, "no previous run to compare", transform=ax2.transAxes,
                 ha="center", va="center", fontsize=9, color="#777777")
    ax2.set_xticks(x)
    ax2.set_xticklabels(ticks, fontsize=8)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
