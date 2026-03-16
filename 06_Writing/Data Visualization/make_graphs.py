import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

OUT_DIR = r"C:\Users\BornLoser\Desktop\Assignment\Thesis\06_Writing\Data Visualization"
os.makedirs(OUT_DIR, exist_ok=True)

DATA = [
    ("ADA","seed_not_found",None,None,0,63),
    ("AICDA","below_threshold",4556,13535,23,28),
    ("AIRE","ok",422,16847,34,36),
    ("AK2","ok",4001,28547,10,10),
    ("AP3B1","ok",4034,183477,6,6),
    ("BIRC4","no_fasta",None,None,None,None),
    ("BLM","ok",4034,96861,51,51),
    ("BLNK","ok",7,10205,1,1),
    ("C1QA","ok",160,2872,4,4),
    ("C1QB","ok",118,9225,3,3),
    ("C1QC","ok",167,4805,4,4),
    ("C1S","ok",4073,14798,5,5),
    ("C2","ok",10,16309,3,3),
    ("C3","below_threshold",1874,45235,7,8),
    ("C5","ok",1371,89880,5,5),
    ("C6","ok",7032,78745,5,5),
    ("C7","below_threshold",11618,82920,6,22),
    ("C8B","ok",3750,29761,6,6),
    ("C9","ok",834,59512,5,5),
    ("CASP10","ok",8,27214,2,2),
    ("CASP8","ok",2,19792,1,1),
    ("CD19","ok",28,6707,2,2),
    ("CD247","ok",2,80236,1,1),
    ("CD3D","ok",112,3260,3,3),
    ("CD3E","ok",38,8944,1,1),
    ("CD3G","ok",4,6486,1,1),
    ("CD40","seed_not_found",None,None,0,3),
    ("CD40L","no_fasta",None,None,None,None),
    ("CD55","ok",942,7191,4,4),
    ("CD59","ok",2,20064,1,1),
    ("CD79A","ok",38,2961,2,2),
    ("CD79B","ok",2,2551,1,1),
    ("CD8A","ok",17,1060,1,1),
    ("CEBPE","below_threshold",381,3605,1,2),
    ("CFD","seed_not_found",None,None,0,3),
    ("CFH","seed_not_found",None,None,0,50),
    ("CFI","ok",4015,64500,13,13),
    ("CFP","ok",3256,8882,12,12),
    ("CIITA","below_threshold",643,47996,5,6),
    ("CTSC","ok",4004,48716,65,65),
    ("CXCR4","ok",4060,8255,6,6),
    ("CYBA","ok",4009,12593,44,45),
    ("CYBB","ok",4048,35890,475,479),
    ("DKC1","ok",4003,16951,6,6),
    ("ELA2","no_fasta",None,None,None,None),
    ("FASLG","seed_not_found",None,None,0,1),
    ("FCGR1A","ok",8,2500,1,1),
    ("FCGR3A","ok",3,2219,1,1),
    ("FERMT3","ok",4055,21415,7,7),
    ("FOXN1","seed_not_found",None,None,0,1),
    ("FOXP3","seed_not_found",None,None,0,28),
    ("G6PC3","ok",4021,10102,12,12),
    ("GFI1","ok",40,11832,3,3),
    ("HAX1","ok",3949,7603,10,10),
    ("ICOS","seed_not_found",None,None,0,1),
    ("IFNGR1","seed_not_found",None,None,0,24),
    ("IFNGR2","ok",21904,52792,5,5),
    ("IGHG2","no_fasta",None,None,None,None),
    ("IGLL1","ok",9,6912,2,2),
    ("IKBKG","below_threshold",5024,27219,8,34),
    ("IL12B","ok",2560,12360,1,1),
    ("IL12RB1","ok",16063,42830,21,21),
    ("IL2RA","ok",601,39616,2,2),
    ("IL7R","ok",6225,23863,6,6),
    ("IRAK4","ok",3121,32471,10,10),
    ("ITGB2","ok",11950,47056,10,10),
    ("JAK3","ok",3519,22773,39,39),
    ("LIG1","ok",5,50064,2,2),
    ("LIG4","ok",4753,11706,11,11),
    ("LYST","ok",20716,211514,28,28),
    ("MAPBPIP","no_fasta",None,None,None,None),
    ("MASP2","ok",3,1632,1,1),
    ("MLPH","ok",4,7246,1,1),
    ("MPO","below_threshold",2501,15064,5,12),
    ("MRE11A","no_fasta",None,None,None,None),
    ("MYO5A","ok",21,178704,2,2),
    ("NCF1","below_threshold",2294,17265,8,11),
    ("NCF2","ok",4001,35307,57,57),
    ("NFKBIA","ok",4,1192,2,2),
    ("NHEJ1","below_threshold",456,25547,5,6),
    ("NP","no_fasta",None,None,None,None),
    ("NRAS","seed_not_found",None,None,0,1),
    ("ORAI1","ok",24,15472,3,3),
    ("PRF1","ok",3601,9691,71,71),
    ("PRKDC","ok",145,102473,1,1),
    ("PTPRC","ok",4002,82747,4,4),
    ("RAB27A","ok",23200,89229,19,19),
    ("RAC2","ok",2,12414,1,1),
    ("RAD50","ok",9,86435,2,2),
    ("RAG1","ok",4001,13267,59,60),
    ("RAG2","ok",4001,10572,22,22),
    ("RASGRP2","ok",2,16937,1,1),
    ("RFX5","seed_not_found",None,None,0,6),
    ("RFXANK","ok",4089,11958,9,9),
    ("RFXAP","seed_not_found",None,None,0,2),
    ("RNF168","ok",237,32739,1,1),
    ("SBDS","ok",4001,12130,32,32),
    ("SERPING1","below_threshold",3949,22023,50,123),
    ("SLC35C1","ok",5145,11588,5,5),
    ("SMARCAL1","ok",6292,70803,33,33),
    ("SP110","ok",4,8570,2,2),
    ("SPINK5","ok",4074,64902,43,44),
    ("STAT1","ok",4563,43421,6,6),
    ("STAT2","ok",4,5456,1,1),
    ("STAT3","ok",6489,71211,31,31),
    ("STAT5B","ok",9,67283,2,2),
    ("STIM1","ok",3,200846,1,1),
    ("STX11","ok",883,38780,3,3),
    ("STXBP2","ok",4020,15332,10,10),
    ("TAP1","seed_not_found",None,None,0,4),
    ("TAP2","seed_not_found",None,None,0,3),
    ("TAZ","no_fasta",None,None,None,None),
    ("TCN2","seed_not_found",None,None,0,11),
    ("TLR3","ok",39,15230,2,2),
    ("TMC6","ok",2925,16950,6,6),
    ("TMC8","below_threshold",4011,9189,5,6),
    ("TNFRSF13B","ok",4001,36737,8,8),
    ("TYK2","ok",653,11682,1,1),
    ("UNC13D","ok",4001,21673,50,52),
    ("UNC93B1","ok",467,8935,2,2),
    ("UNG","ok",4001,10953,5,5),
    ("ZAP70","below_threshold",4001,29518,8,9),
]

df = pd.DataFrame(DATA, columns=["gene","status","sstart","send","matched","total"])

PALETTE = {
    "ok":              "#2E8B57",
    "below_threshold": "#D48B1A",
    "seed_not_found":  "#C0392B",
    "no_fasta":        "#7F8C8D",
}
LABELS = {
    "ok":              "OK (≥90%)",
    "below_threshold": "Below threshold",
    "seed_not_found":  "Seed not found",
    "no_fasta":        "No FASTA",
}

plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.color":         "#e8e8e8",
    "grid.linewidth":     0.6,
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
})

# ── FIGURE 1: Status breakdown donut + bar ────────────────────────────────────
counts = df["status"].value_counts()
ordered = ["ok","below_threshold","seed_not_found","no_fasta"]
counts = counts.reindex(ordered).fillna(0).astype(int)
colors = [PALETTE[s] for s in ordered]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("LRG Offset Search — Status Overview", fontsize=15, fontweight="bold", y=1.01)

ax1 = axes[0]
wedges, texts, autotexts = ax1.pie(
    counts, labels=None, colors=colors,
    autopct=lambda p: f"{p:.1f}%\n({int(round(p*sum(counts)/100))})",
    pctdistance=0.72, startangle=90,
    wedgeprops={"linewidth": 1.5, "edgecolor": "white"},
)
for at in autotexts:
    at.set_fontsize(9)
    at.set_fontweight("bold")
    at.set_color("white")
centre = plt.Circle((0, 0), 0.48, color="white")
ax1.add_patch(centre)
ax1.text(0, 0, f"{sum(counts)}\ngenes", ha="center", va="center",
         fontsize=12, fontweight="bold", color="#2c2c2c")
legend_handles = [mpatches.Patch(color=PALETTE[s], label=f"{LABELS[s]}  ({counts[s]})") for s in ordered]
ax1.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.12),
           ncol=2, frameon=False, fontsize=9)
ax1.set_title("Gene count by status", fontsize=11, pad=10)

ax2 = axes[1]
ok_df = df[df["status"].isin(["ok","below_threshold"])].copy()
ok_df["match_pct"] = ok_df["matched"] / ok_df["total"] * 100
bins = [0,10,20,30,40,50,60,70,80,90,100]
ok_df["bin"] = pd.cut(ok_df["match_pct"], bins=bins, right=True, include_lowest=True)
bin_counts = ok_df.groupby("bin", observed=False).size()
bin_colors = ["#C0392B" if b.right <= 50 else "#D48B1A" if b.right <= 90 else "#2E8B57"
              for b in bin_counts.index]
ax2.bar(range(len(bin_counts)), bin_counts.values, color=bin_colors,
        edgecolor="white", linewidth=1.2, width=0.75)
ax2.set_xticks(range(len(bin_counts)))
ax2.set_xticklabels([f"{int(b.left)}–{int(b.right)}%" for b in bin_counts.index],
                    rotation=40, ha="right", fontsize=8.5)
ax2.set_ylabel("Number of genes", fontsize=10)
ax2.set_title("Match rate distribution\n(ok + below threshold genes)", fontsize=11)
for i, v in enumerate(bin_counts.values):
    if v > 0:
        ax2.text(i, v + 0.3, str(v), ha="center", va="bottom", fontsize=9, fontweight="bold")

plt.tight_layout()
p1 = os.path.join(OUT_DIR, "fig1_status_overview.png")
plt.savefig(p1, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {p1}")

# ── FIGURE 2: Below-threshold genes — matched vs unmatched stacked bar ────────
bt = df[df["status"] == "below_threshold"].copy()
bt["match_pct"] = bt["matched"] / bt["total"] * 100
bt = bt.sort_values("match_pct")
unmatched = bt["total"] - bt["matched"]

fig, ax = plt.subplots(figsize=(10, 5))
fig.suptitle("Below-threshold Genes — Matched vs Unmatched Mutations",
             fontsize=14, fontweight="bold")

x = np.arange(len(bt))
b1 = ax.bar(x, bt["matched"].values, color="#2E8B57", label="Matched", width=0.6)
b2 = ax.bar(x, unmatched.values, bottom=bt["matched"].values, color="#C0392B",
            label="Unmatched", width=0.6)

for i, (_, row) in enumerate(bt.iterrows()):
    pct = row["matched"] / row["total"] * 100
    ax.text(i, row["total"] + 0.8, f"{pct:.0f}%",
            ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#2c2c2c")

ax.set_xticks(x)
ax.set_xticklabels(bt["gene"].values, rotation=35, ha="right", fontsize=9)
ax.set_ylabel("Mutation count", fontsize=10)
ax.legend(frameon=False, fontsize=10)
ax.set_xlim(-0.6, len(bt) - 0.4)

plt.tight_layout()
p2 = os.path.join(OUT_DIR, "fig2_below_threshold_detail.png")
plt.savefig(p2, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {p2}")

# ── FIGURE 3: Fragment size (send - sstart) for ok genes, sorted ──────────────
ok = df[df["status"] == "ok"].copy()
ok["frag_len"] = ok["send"] - ok["sstart"]
ok = ok.sort_values("frag_len", ascending=False)

fig, ax = plt.subplots(figsize=(14, 5))
fig.suptitle("IDRefSeq Fragment Length within LRG (ok genes)",
             fontsize=14, fontweight="bold")

cmap_vals = (ok["frag_len"] - ok["frag_len"].min()) / (ok["frag_len"].max() - ok["frag_len"].min())
bar_colors = plt.cm.YlGn(0.3 + 0.65 * cmap_vals.values)

ax.bar(range(len(ok)), ok["frag_len"].values / 1000, color=bar_colors,
       edgecolor="white", linewidth=0.8, width=0.8)
ax.set_xticks(range(len(ok)))
ax.set_xticklabels(ok["gene"].values, rotation=90, fontsize=7)
ax.set_ylabel("Fragment length (kb)", fontsize=10)
ax.axhline(ok["frag_len"].median() / 1000, color="#C0392B", linewidth=1.4,
           linestyle="--", label=f"Median: {ok['frag_len'].median()/1000:.1f} kb")
ax.legend(frameon=False, fontsize=10)

plt.tight_layout()
p3 = os.path.join(OUT_DIR, "fig3_fragment_lengths.png")
plt.savefig(p3, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {p3}")

print("\nAll figures saved.")
