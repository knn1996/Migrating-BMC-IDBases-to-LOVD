"""
visualize_skipped_LRG.py
========================
Quick-check dashboard for bed_to_vcf_skipped_LRG.tsv.

Run this script after bed_to_vcf.py to instantly open an interactive
summary in your browser. Or import and call check() from bed_to_vcf.py.

Usage (standalone):
    python visualize_skipped_LRG.py

Usage (from bed_to_vcf.py — add at the bottom of that script):
    from visualize_skipped_LRG import check
    check()          # generates HTML and opens browser
"""

import os
import sys
import webbrowser
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR  = os.path.join(SCRIPTS_DIR, "..", "Output")   # ../Output/
TSV_PATH    = os.path.join(OUTPUT_DIR, "bed_to_vcf_skipped_LRG.tsv")
OUTPUT_HTML = os.path.join(OUTPUT_DIR, "bed_to_vcf_skipped_LRG_dashboard.html")

REASON_COLORS = {
    "truncated (missing ALT)":       "#e07b54",
    "unparseable — review manually": "#c0392b",
    "accession-only":                "#5b8db8",
}
GENE_PALETTE = px.colors.qualitative.Safe


def build_dashboard(tsv_path: str = TSV_PATH, output_html: str = OUTPUT_HTML) -> str:
    """Read the TSV and write a self-contained HTML dashboard. Returns the output path."""

    # ── Load ──────────────────────────────────────────────────────────────────
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"Cannot find TSV: {tsv_path}")

    df       = pd.read_csv(tsv_path, sep="\t")
    df_dedup = df.drop_duplicates()

    total_raw   = len(df)
    total_dedup = len(df_dedup)
    n_dup       = total_raw - total_dedup
    n_genes     = df["gene"].nunique()
    n_files     = df["file"].nunique()
    n_reasons   = df["reason"].nunique()

    genes   = sorted(df["gene"].unique())
    reasons = sorted(df["reason"].unique())
    gene_colors = {g: GENE_PALETTE[i % len(GENE_PALETTE)] for i, g in enumerate(genes)}

    # ── Subplots: row1=KPIs (indicators), row2=pie+bar, row3=stacked+chrom ────
    fig = make_subplots(
        rows=3, cols=6,
        row_heights=[0.14, 0.40, 0.46],
        specs=[
            [{"type": "indicator"}] * 6,
            [{"type": "pie", "colspan": 3}, None, None,
             {"type": "xy",  "colspan": 3}, None, None],
            [{"type": "xy",  "colspan": 3}, None, None,
             {"type": "xy",  "colspan": 3}, None, None],
        ],
        subplot_titles=(
            "Total rows", "Unique rows", "Duplicate rows",
            "Genes", "Source files", "Skip reasons",
            "Skip Reason Breakdown", "", "", "Skipped Rows per Gene", "", "",
            "Skip Reasons per Gene (stacked)", "", "", "Chromosome Distribution (top 15)", "", "",
        ),
        vertical_spacing=0.06,
        horizontal_spacing=0.04,
    )

    # ── Row 1 — KPI Indicators ────────────────────────────────────────────────
    kpi_values = [total_raw, total_dedup, n_dup, n_genes, n_files, n_reasons]
    kpi_colors = ["#2c3e50", "#27ae60", "#e07b54", "#2980b9", "#8e44ad", "#16a085"]
    for col_idx, (val, color) in enumerate(zip(kpi_values, kpi_colors), start=1):
        fig.add_trace(
            go.Indicator(
                mode="number",
                value=val,
                number=dict(font=dict(size=32, color=color, family="Segoe UI, Arial")),
            ),
            row=1, col=col_idx,
        )

    # ── Row 2 col 1-3 — Donut: reason breakdown ───────────────────────────────
    reason_counts = df["reason"].value_counts()
    fig.add_trace(
        go.Pie(
            labels=reason_counts.index.tolist(),
            values=reason_counts.values.tolist(),
            marker_colors=[REASON_COLORS.get(r, "#aaa") for r in reason_counts.index],
            hole=0.42,
            textinfo="percent+label",
            hovertemplate="<b>%{label}</b><br>Count: %{value:,}<br>Share: %{percent}<extra></extra>",
            showlegend=True,
        ),
        row=2, col=1,
    )

    # ── Row 2 col 4-6 — Bar: skipped per gene ─────────────────────────────────
    gene_counts = df["gene"].value_counts().reindex(genes)
    fig.add_trace(
        go.Bar(
            x=genes,
            y=gene_counts.values.tolist(),
            marker_color=[gene_colors[g] for g in genes],
            text=gene_counts.values.tolist(),
            textposition="outside",
            hovertemplate="<b>%{x}</b><br>Skipped rows: %{y:,}<extra></extra>",
            showlegend=False,
        ),
        row=2, col=4,
    )

    # ── Row 3 col 1-3 — Stacked bar: reasons per gene ─────────────────────────
    gene_reason = df.groupby(["gene", "reason"]).size().unstack(fill_value=0)
    for reason in reasons:
        if reason not in gene_reason.columns:
            gene_reason[reason] = 0
        fig.add_trace(
            go.Bar(
                name=reason,
                x=gene_reason.index.tolist(),
                y=gene_reason[reason].tolist(),
                marker_color=REASON_COLORS.get(reason, "#aaa"),
                hovertemplate=f"<b>%{{x}}</b><br>{reason}: %{{y:,}}<extra></extra>",
            ),
            row=3, col=1,
        )

    # ── Row 3 col 4-6 — Bar: chrom distribution (top 15) ─────────────────────
    chrom_counts = df["chrom"].value_counts().head(15)

    def chrom_sort_key(c):
        c = str(c).replace("chr", "")
        if c.isdigit(): return (0, int(c))
        if c == "X":    return (1, 0)
        if c == "Y":    return (1, 1)
        return (2, 0)

    sorted_chroms = sorted(chrom_counts.index, key=chrom_sort_key)
    chrom_counts  = chrom_counts.reindex(sorted_chroms)

    fig.add_trace(
        go.Bar(
            x=chrom_counts.index.tolist(),
            y=chrom_counts.values.tolist(),
            marker_color="#5b8db8",
            hovertemplate="<b>%{x}</b><br>Rows: %{y:,}<extra></extra>",
            showlegend=False,
        ),
        row=3, col=4,
    )

    # ── Axis labels ───────────────────────────────────────────────────────────
    fig.update_xaxes(title_text="Gene",       row=2, col=4)
    fig.update_yaxes(title_text="Row count",  row=2, col=4)
    fig.update_xaxes(title_text="Gene",       row=3, col=1)
    fig.update_yaxes(title_text="Row count",  row=3, col=1)
    fig.update_xaxes(title_text="Chromosome", row=3, col=4, tickangle=30)
    fig.update_yaxes(title_text="Row count",  row=3, col=4)

    # ── Layout ────────────────────────────────────────────────────────────────
    fig.update_layout(
        title=dict(
            text="<b>bed_to_vcf — Skipped LRG Entries Dashboard</b>",
            font=dict(family="Segoe UI, Arial", size=20, color="#2c3e50"),
            x=0.5,
        ),
        barmode="stack",
        template="plotly_white",
        font=dict(family="Segoe UI, Arial", size=12),
        legend=dict(
            title="<b>Skip Reason</b>",
            orientation="h",
            x=0.5, xanchor="center",
            y=-0.04,
            font=dict(size=11),
        ),
        margin=dict(t=100, b=80, l=60, r=40),
        height=950,
        paper_bgcolor="#f8f9fa",
        plot_bgcolor="#ffffff",
    )

    # ── Write HTML (Plotly bundled — works fully offline) ─────────────────────
    fig.write_html(
        output_html,
        full_html=True,
        include_plotlyjs=True,    # embeds ~3 MB of Plotly.js — no internet needed
        config={"scrollZoom": True, "displayModeBar": True},
    )

    return os.path.abspath(output_html)


def check(auto_open: bool = True) -> None:
    """Build the dashboard and (optionally) open it in the default browser."""
    print("[visualize_skipped_LRG] Reading:", TSV_PATH)
    try:
        html_path = build_dashboard()
        print(f"[visualize_skipped_LRG] Dashboard written → {html_path}")
        if auto_open:
            webbrowser.open(f"file:///{html_path.replace(os.sep, '/')}")
            print("[visualize_skipped_LRG] Opened in browser ✓")
    except FileNotFoundError as e:
        print(f"[visualize_skipped_LRG] ERROR: {e}", file=sys.stderr)
    except Exception as e:
        print(f"[visualize_skipped_LRG] ERROR building dashboard: {e}", file=sys.stderr)
        raise


# ── Standalone entry point ────────────────────────────────────────────────────
if __name__ == "__main__":
    check(auto_open=True)
