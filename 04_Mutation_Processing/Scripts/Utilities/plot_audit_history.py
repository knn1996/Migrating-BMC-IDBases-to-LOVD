"""
plot_audit_history.py
=====================
Render the Step10 audit trend from history/audit_history.tsv.

Top panel: in-scope recovery (in_final_pct), protein-good rate
(protein_good_pct) and final pair consistency (pair_consistency_pct) per run.
Middle panel: gains and regressions versus the previous run (regressions drawn
below zero). Bottom panel: excess pairs (final pairs beyond the source DNA
changes of each individual). Runs are plotted in history order and labelled
with run_id and label; a metric a run did not record is left out.

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
LINES = [("in_final_pct", "In-scope recovery (%)", "#2a6fb0", 7),
         ("protein_good_pct", "Protein check good (%)", "#c0662b", 7),
         ("pair_consistency_pct", "Final pair consistency (%)", "#7a4fa3", -13)]
BAR_COLORS = {"gains": "#3a8f5c", "regressions": "#b8433a", "excess_pairs": "#8a8a3a"}


def default_history_dir():
    with open(os.environ.get("CONFIG", SCRIPTS / "config" / "config.yaml"), encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return Path(cfg["output_root"]) / "Step10_Audit" / "history"


def num(h, col):
    if col not in h.columns:
        return pd.Series([float("nan")] * len(h), index=h.index)
    return pd.to_numeric(h[col], errors="coerce")


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

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(max(6.5, 1.6 * len(h) + 3), 8.5), sharex=True,
                                        gridspec_kw={"height_ratios": [3, 2, 2]})
    lows = []
    for col, name, color, dy in LINES:
        y = num(h, col)
        if y.isna().all():
            continue
        lows.append(y.min())
        ax1.plot(x, y, marker="o", color=color, label=name, linewidth=2)
        for xi, yi in zip(x, y):
            if pd.notna(yi):
                ax1.annotate(f"{yi:.1f}", (xi, yi), textcoords="offset points", xytext=(0, dy),
                             ha="center", fontsize=8, color=color)
    ax1.set_ylim(max(0, min(lows) - 10), 100)
    ax1.set_ylabel("%")
    ax1.grid(axis="y", alpha=0.3)
    ax1.legend(loc="lower right", frameon=False, fontsize=9)
    ax1.set_title("IDbases → LOVD audit history", fontsize=11)

    gains = num(h, "gains").fillna(0)
    regs = num(h, "regressions").fillna(0)
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

    excess = num(h, "excess_pairs")
    ax3.bar(x, excess.fillna(0), width=0.5, color=BAR_COLORS["excess_pairs"], label="Excess final pairs")
    for xi, e in zip(x, excess):
        if pd.notna(e):
            ax3.annotate(f"{int(e)}", (xi, e), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8)
    ax3.set_ylabel("Pairs")
    ax3.grid(axis="y", alpha=0.3)
    ax3.legend(loc="upper right", frameon=False, fontsize=9)
    if excess.isna().all():
        ax3.text(0.5, 0.5, "not recorded", transform=ax3.transAxes, ha="center", va="center",
                 fontsize=9, color="#777777")
    else:
        ax3.set_ylim(0, excess.max() * 1.18)
    ax3.set_xticks(x)
    ax3.set_xticklabels(ticks, fontsize=8)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
