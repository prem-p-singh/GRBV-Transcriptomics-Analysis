#!/usr/bin/env python
"""Regenerate main-paper Figures 49 (mediation diagrams) and 50 (direct vs
indirect effects) from the fix03 tv=2000 mediation_results.csv.

The original causal_mediation_analysis.py wrote both figures using a column
schema that fix03 replaced (gene/total_c/path_b_padj/sig_joint instead of
gene/total_effect_c/total_padj/category). It also used anonymous LOC IDs as
panel titles. This script reads the new CSV plus the fix11 NCBI annotations
and writes drop-in replacements alongside the rest of the main-paper PNGs.

Run from any directory; paths are absolute.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

ap = argparse.ArgumentParser()
ap.add_argument("--mediation", default="/Users/prempratapsingh/Data_Analysis/mediation_results.csv")
ap.add_argument("--annot",     default="/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/fix11_annotations/gene_annotations.csv")
ap.add_argument("--fig49",     default="/Users/prempratapsingh/Data_Analysis/49_mediation_diagrams.png")
ap.add_argument("--fig50",     default="/Users/prempratapsingh/Data_Analysis/50_direct_vs_indirect_effects.png")
ap.add_argument("--top_var",   type=int, default=2000)
args = ap.parse_args()

plt.rcParams.update({
    "figure.dpi": 300, "font.size": 14, "axes.titlesize": 15,
    "axes.labelsize": 14, "figure.facecolor": "white",
})

med = pd.read_csv(args.mediation)
annot = pd.read_csv(args.annot)
name_map = dict(zip(annot["gene_id"], annot["product"].fillna("")))

n_total = len(med)
n_sig = int(med["sig_joint"].sum())

sig = med[med["sig_joint"]].copy()
sig = sig.sort_values("path_b_padj")

cat_temp_only = (med["total_padj"] < 0.05) & (~med["sig_joint"])
cat_titer_med = med["sig_joint"]
cat_neither   = (med["total_padj"] >= 0.05) & (~med["sig_joint"])

med["category"] = np.select(
    [cat_titer_med, cat_temp_only],
    ["Titer-mediated (joint sig.)", "Temperature-only"],
    default="Not significant",
)

# =====================================================================
# FIGURE 49 — six-panel mediation diagrams, MBF1c front-panel
# =====================================================================
top6 = sig.head(6)
print(f"[fig49] top 6 mediators by joint p:")
for _, r in top6.iterrows():
    print(f"  {r['gene']:<16} {name_map.get(r['gene'], '')[:60]}  path_b_padj={r['path_b_padj']:.2e}")

fig, axes = plt.subplots(2, 3, figsize=(24, 16))
fig.suptitle(
    "Causal mediation: Temperature  →  Viral titer  →  Gene expression\n"
    f"(top {len(top6)} of {n_sig} joint-significant mediators at tv = {args.top_var}; "
    f"BH FDR < 0.05 on both Temp→Gene and Titer→Gene paths)",
    fontsize=18, fontweight="bold", y=0.995,
)

for idx, (_, row) in enumerate(top6.iterrows()):
    ax = axes[idx // 3, idx % 3]
    gene = row["gene"]
    gname = name_map.get(gene, "")
    if len(gname) > 50:
        gname = gname[:47] + "..."

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect("equal")
    ax.axis("off")

    temp_box = FancyBboxPatch((0.5, 4), 2.5, 1.5, boxstyle="round,pad=0.2",
                              facecolor="#f39c12", alpha=0.85,
                              edgecolor="black", linewidth=1.5)
    ax.add_patch(temp_box)
    ax.text(1.75, 4.75, "Temperature", ha="center", va="center",
            fontsize=12, fontweight="bold")

    gene_box = FancyBboxPatch((7, 4), 2.5, 1.5, boxstyle="round,pad=0.2",
                              facecolor="#2ecc71", alpha=0.85,
                              edgecolor="black", linewidth=1.5)
    ax.add_patch(gene_box)
    ax.text(8.25, 4.75, gene, ha="center", va="center",
            fontsize=11, fontweight="bold")

    titer_box = FancyBboxPatch((3.5, 7.5), 2.5, 1.5, boxstyle="round,pad=0.2",
                               facecolor="#e74c3c", alpha=0.85,
                               edgecolor="black", linewidth=1.5)
    ax.add_patch(titer_box)
    ax.text(4.75, 8.25, "Viral Titer", ha="center", va="center",
            fontsize=12, fontweight="bold")

    a_val = row["indirect_ab"] / row["path_b"] if abs(row["path_b"]) > 1e-20 else 0.0
    ax.annotate("", xy=(3.5, 8), xytext=(3, 5.5),
                arrowprops=dict(arrowstyle="->", color="#2c3e50", lw=2))
    ax.text(2.4, 7.0, f"a = {a_val:,.0f}", fontsize=10, rotation=50,
            ha="center", color="#2c3e50")

    ax.annotate("", xy=(7, 5.5), xytext=(6, 8),
                arrowprops=dict(arrowstyle="->", color="#2c3e50", lw=2))
    ax.text(7.1, 7.0, f"b = {row['path_b']:.2e}", fontsize=10, rotation=-50,
            ha="center", color="#2c3e50")

    ax.annotate("", xy=(7, 4.75), xytext=(3, 4.75),
                arrowprops=dict(arrowstyle="->", color="#e74c3c",
                                lw=2.5, linestyle="--"))
    ax.text(5, 4, f"c' = {row['direct_cprime']:.3f}",
            fontsize=11, ha="center", color="#e74c3c", fontweight="bold")

    ax.text(5, 2.6,
            f"Total c = {row['total_c']:.3f}  (p = {row['total_p']:.1e})",
            ha="center", fontsize=10)
    ax.text(5, 1.95,
            f"Indirect a×b = {row['indirect_ab']:.3f}  "
            f"(joint padj = {row['path_b_padj']:.1e})",
            ha="center", fontsize=10)
    pct = row["prop_mediated"] * 100
    ax.text(5, 1.30, f"% mediated = {pct:+.1f}%",
            ha="center", fontsize=11, fontweight="bold", color="#8e44ad")

    if gname:
        ax.text(5, 0.45, gname, ha="center", fontsize=10,
                style="italic", color="#555")

    ax.set_title(gene, fontsize=14, fontweight="bold", pad=4)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(args.fig49, dpi=300, bbox_inches="tight")
plt.close()
print(f"[fig49] wrote {args.fig49}")

# =====================================================================
# FIGURE 50 — direct vs indirect scatter + category bar
# =====================================================================
fig, axes = plt.subplots(1, 2, figsize=(20, 9))
fig.suptitle(
    f"Temperature effects on gene expression: direct vs titer-mediated "
    f"(tv = {args.top_var}, n = {n_total} variable genes)",
    fontsize=18, fontweight="bold", y=0.99,
)

cat_colors = {
    "Titer-mediated (joint sig.)": "#e74c3c",
    "Temperature-only":            "#f39c12",
    "Not significant":             "#bdc3c7",
}

ax = axes[0]
for cat, color in cat_colors.items():
    sub = med[med["category"] == cat]
    if not sub.empty:
        ax.scatter(sub["direct_cprime"], sub["indirect_ab"],
                   c=color, alpha=0.55, s=30,
                   label=f"{cat}  ({len(sub)})", edgecolors="none")

ax.axhline(0, color="grey", linewidth=0.5)
ax.axvline(0, color="grey", linewidth=0.5)
ax.set_xlabel("Direct effect  c'  (Temp → Gene | Titer)")
ax.set_ylabel("Indirect effect  a × b  (Temp → Titer → Gene)")
ax.set_title("A. Direct vs indirect temperature effects",
             fontsize=15, fontweight="bold")
ax.legend(fontsize=10, loc="best")

label_top = sig.head(8)
for _, r in label_top.iterrows():
    annot_label = r["gene"]
    short = name_map.get(r["gene"], "")
    if short:
        annot_label = f"{annot_label}\n{short.split(',')[0][:24]}"
    ax.annotate(annot_label,
                (r["direct_cprime"], r["indirect_ab"]),
                fontsize=8, alpha=0.85,
                xytext=(5, 5), textcoords="offset points")

ax = axes[1]
cat_counts = med["category"].value_counts()
order = ["Titer-mediated (joint sig.)", "Temperature-only", "Not significant"]
cat_counts = cat_counts.reindex([c for c in order if c in cat_counts.index])

colors_pie = [cat_colors[c] for c in cat_counts.index]
ax.barh(range(len(cat_counts)), cat_counts.values, color=colors_pie,
        alpha=0.9, edgecolor="white")
ax.set_yticks(range(len(cat_counts)))
ax.set_yticklabels(cat_counts.index, fontsize=12)
ax.set_xlabel("Number of genes")
ax.set_title("B. Gene classification by mediation",
             fontsize=15, fontweight="bold")
ax.invert_yaxis()
xmax = float(cat_counts.values.max())
for i, v in enumerate(cat_counts.values):
    ax.text(v + xmax * 0.02, i, f"{v}", va="center", fontsize=12)

ax.text(0.97, 0.07,
        f"{n_sig} / {n_total}  joint-significant\n(BH FDR<0.05 on both paths,\n |a×b|>0.05)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.85))

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(args.fig50, dpi=300, bbox_inches="tight")
plt.close()
print(f"[fig50] wrote {args.fig50}")
