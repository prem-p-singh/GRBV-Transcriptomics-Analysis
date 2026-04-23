#!/usr/bin/env python
"""12. Comparative meta-analysis: GRBV DEG signature vs harmonized literature.

Operates on harmonized DEG tables (rule 11 outputs) plus a GRBV reference
signature. Produces:
  - Per-pair overlap statistics (Fisher exact, signed concordance)
  - Per-pair log2FC Spearman correlation on shared DEGs
  - Per-gene specificity classification
  - UpSet plot of DEG set overlaps
  - Fold-change correlation heatmap

Generates the GRBV signature on the fly from the counts matrix + metadata
(pooled RB+ vs RB- Welch t-test with BH FDR) to guarantee it uses the same
gene-ID space (LOC_) as the harmonized literature tables.
"""
from __future__ import annotations
import argparse, os, sys, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import from_contents, UpSet

ap = argparse.ArgumentParser()
ap.add_argument("--counts",       required=True)
ap.add_argument("--meta",         required=True)
ap.add_argument("--harmonized_dir", required=True)
ap.add_argument("--coverage",     required=True)
ap.add_argument("--out_dir",      required=True)
ap.add_argument("--padj_thresh",  type=float, default=0.05)
ap.add_argument("--lfc_thresh",   type=float, default=1.0)
args = ap.parse_args()

os.makedirs(args.out_dir, exist_ok=True)


# -----------------------------------------------------------------------------
# 1. Build GRBV DEG signature (pooled RB+ vs RB-, Welch t-test, BH FDR)
# -----------------------------------------------------------------------------
print("[12] computing GRBV DEG signature (pooled RB+ vs RB-)", file=sys.stderr)
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta   = pd.read_csv(args.meta,   sep="\t")

# Gene-level filter (>=10 counts in >=10% samples) mirrors the main paper
keep = (counts >= 10).sum(axis=1) >= 0.10 * counts.shape[1]
counts = counts.loc[keep]
lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)

rb_pos = meta.loc[meta["Infection_Group"] == "RB(+)", "Sample"].tolist()
rb_neg = meta.loc[meta["Infection_Group"] == "RB(-)", "Sample"].tolist()
rb_pos = [s for s in rb_pos if s in lcpm.columns]
rb_neg = [s for s in rb_neg if s in lcpm.columns]

a = lcpm[rb_pos].values
b = lcpm[rb_neg].values
t, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
lfc = a.mean(axis=1) - b.mean(axis=1)
padj = multipletests(np.nan_to_num(p, nan=1.0), method="fdr_bh")[1]

grbv = pd.DataFrame({
    "LOC_id": lcpm.index,
    "log2FC": lfc,
    "pvalue": p,
    "padj":   padj,
    "tested": 1,
})
grbv = grbv.dropna(subset=["pvalue"])
grbv.to_csv(Path(args.out_dir) / "GRBV_signature.csv", index=False)

def sig_set(df, padj_thr=args.padj_thresh, lfc_thr=args.lfc_thresh):
    return set(df.loc[(df["padj"] < padj_thr) & (df["log2FC"].abs() >= lfc_thr),
                      "LOC_id"].dropna())

grbv_sig = sig_set(grbv)
print(f"[12] GRBV signature: {len(grbv_sig)} DEGs (padj<{args.padj_thresh}, |logFC|>={args.lfc_thresh})",
      file=sys.stderr)


# -----------------------------------------------------------------------------
# 2. Load harmonized literature datasets
# -----------------------------------------------------------------------------
cov = pd.read_csv(args.coverage)
datasets: dict[str, pd.DataFrame] = {"GRBV": grbv}
labels: dict[str, str] = {"GRBV": "GRBV (this study)"}
stress_of = {"GRBV": "GRBV"}

for _, row in cov.iterrows():
    name = row["dataset"]
    f = Path(args.harmonized_dir) / f"{name}.harmonized.csv"
    if not f.exists():
        continue
    df = pd.read_csv(f)
    df = df.dropna(subset=["LOC_id"])
    # Aggregate duplicate LOC_ids by taking most-significant hit
    df = df.sort_values("padj").drop_duplicates("LOC_id", keep="first")
    datasets[name] = df
    labels[name] = f"{row['stress']} ({name})"
    stress_of[name] = row["stress"]

print(f"[12] loaded {len(datasets)} datasets (including GRBV)", file=sys.stderr)


# -----------------------------------------------------------------------------
# 3. Pairwise overlap (Fisher) and FC correlation (Spearman)
# -----------------------------------------------------------------------------
names = list(datasets)
overlap_rows = []

for a_name, b_name in itertools.combinations(names, 2):
    A = datasets[a_name]; B = datasets[b_name]
    # Universe = intersection of tested genes
    uA = set(A["LOC_id"]); uB = set(B["LOC_id"])
    universe = uA & uB
    if len(universe) < 100:
        continue
    sigA = sig_set(A) & universe
    sigB = sig_set(B) & universe
    both = sigA & sigB
    onlyA = sigA - sigB; onlyB = sigB - sigA
    neither = universe - sigA - sigB

    table = [[len(both), len(onlyA)], [len(onlyB), len(neither)]]
    odds, p = stats.fisher_exact(table, alternative="greater")

    # FC correlation on shared tested genes
    shared = A.set_index("LOC_id").loc[list(universe), "log2FC"].to_frame("A")
    shared["B"] = B.set_index("LOC_id").loc[list(universe), "log2FC"].values
    shared = shared.dropna()
    rho = stats.spearmanr(shared["A"], shared["B"]).statistic if len(shared) > 10 else np.nan

    # Signed concordance on overlapping DEGs
    if both:
        fcA = A.set_index("LOC_id").loc[list(both), "log2FC"]
        fcB = B.set_index("LOC_id").loc[list(both), "log2FC"]
        mask = fcA.notna() & fcB.notna()
        if mask.any():
            same_sign = float((np.sign(fcA[mask].values) == np.sign(fcB[mask].values)).mean())
        else:
            same_sign = np.nan
    else:
        same_sign = np.nan

    overlap_rows.append({
        "dataset_A": a_name, "dataset_B": b_name,
        "universe_size": len(universe),
        "sig_A": len(sigA), "sig_B": len(sigB),
        "overlap_sig": len(both),
        "odds_ratio": round(float(odds), 3) if np.isfinite(odds) else np.nan,
        "fisher_p": float(p),
        "signed_concordance": round(float(same_sign), 3) if same_sign == same_sign else np.nan,
        "spearman_rho_logFC": round(float(rho), 3) if rho == rho else np.nan,
    })

overlap_df = pd.DataFrame(overlap_rows)
# BH-adjust Fisher p-values
if len(overlap_df):
    overlap_df["fisher_q"] = multipletests(overlap_df["fisher_p"], method="fdr_bh")[1]
    overlap_df = overlap_df.sort_values("fisher_q")
overlap_df.to_csv(Path(args.out_dir) / "overlap_matrix.csv", index=False)
print(f"[12] wrote {len(overlap_df)} pairwise overlap rows", file=sys.stderr)


# -----------------------------------------------------------------------------
# 4. Per-GRBV-DEG specificity classification
# -----------------------------------------------------------------------------
def classify(loc_id: str) -> dict:
    row = {"LOC_id": loc_id}
    for name, df in datasets.items():
        sub = df[df["LOC_id"] == loc_id]
        if len(sub) == 0 or sub["tested"].iloc[0] == 0:
            row[f"{name}_logFC"] = np.nan; row[f"{name}_padj"] = np.nan
            row[f"{name}_sig"] = np.nan
        else:
            row[f"{name}_logFC"] = sub["log2FC"].iloc[0]
            row[f"{name}_padj"] = sub["padj"].iloc[0]
            row[f"{name}_sig"] = int(
                (sub["padj"].iloc[0] < args.padj_thresh)
                and (abs(sub["log2FC"].iloc[0]) >= args.lfc_thresh)
            )
    return row

class_rows = [classify(g) for g in grbv_sig]
cls = pd.DataFrame(class_rows)

# Stress groupings — use user-provided labels to partition datasets
viral = [n for n in datasets if stress_of.get(n, "").upper() in {"GLRAV-3", "GPGV", "GRBV"} and n != "GRBV"]
mildew = [n for n in datasets if "mildew" in stress_of.get(n, "").lower()]
abiotic = [n for n in datasets if stress_of.get(n, "").lower() in {"heat", "drought", "cold"}]

def count_sig(row, cols): return sum(int(row[f"{n}_sig"]) == 1 for n in cols
                                     if f"{n}_sig" in row and not pd.isna(row[f"{n}_sig"]))

if len(cls):
    cls["n_viral_sig"]   = cls.apply(lambda r: count_sig(r, viral), axis=1)
    cls["n_mildew_sig"]  = cls.apply(lambda r: count_sig(r, mildew), axis=1)
    cls["n_abiotic_sig"] = cls.apply(lambda r: count_sig(r, abiotic), axis=1)
    cls["n_other_sig"]   = cls["n_viral_sig"] + cls["n_mildew_sig"] + cls["n_abiotic_sig"]

    def classify_specificity(r):
        if r["n_other_sig"] == 0:
            return "GRBV-unique"
        if r["n_mildew_sig"] == 0 and r["n_abiotic_sig"] == 0 and r["n_viral_sig"] >= 1:
            return "Viral-shared"
        if r["n_mildew_sig"] >= 1 and r["n_abiotic_sig"] == 0:
            return "Pathogen-generic"
        if r["n_abiotic_sig"] >= 1:
            return "Stress-generic"
        return "Mixed"

    cls["specificity_class"] = cls.apply(classify_specificity, axis=1)

cls.to_csv(Path(args.out_dir) / "specificity_classes.csv", index=False)
if len(cls):
    summary = cls["specificity_class"].value_counts().rename_axis("class").reset_index(name="n_GRBV_DEGs")
    summary.to_csv(Path(args.out_dir) / "specificity_summary.csv", index=False)
    print(f"[12] specificity breakdown:\n{summary}", file=sys.stderr)


# -----------------------------------------------------------------------------
# 5. UpSet plot of DEG set overlaps
# -----------------------------------------------------------------------------
try:
    contents = {labels[n]: sig_set(datasets[n]) for n in names if len(sig_set(datasets[n]))}
    up_data = from_contents({k: v for k, v in contents.items() if len(v)})
    fig = plt.figure(figsize=(11, 5))
    UpSet(up_data, subset_size="count", sort_by="cardinality", show_counts=True).plot(fig=fig)
    fig.suptitle("DEG overlaps across Vitis stress datasets", y=1.02, fontsize=11)
    fig.savefig(Path(args.out_dir) / "upset_plot.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
except Exception as e:
    print(f"[12] UpSet plot failed ({e})", file=sys.stderr)


# -----------------------------------------------------------------------------
# 6. Fold-change correlation heatmap
# -----------------------------------------------------------------------------
mat = pd.DataFrame(np.nan, index=names, columns=names)
for _, r in overlap_df.iterrows():
    mat.loc[r["dataset_A"], r["dataset_B"]] = r["spearman_rho_logFC"]
    mat.loc[r["dataset_B"], r["dataset_A"]] = r["spearman_rho_logFC"]
np.fill_diagonal(mat.values, 1.0)
mat.to_csv(Path(args.out_dir) / "fc_correlation_matrix.csv")

fig, ax = plt.subplots(figsize=(7, 5.5))
sns.heatmap(mat.astype(float), annot=True, fmt=".2f", center=0, vmin=-1, vmax=1,
            cmap="RdBu_r", cbar_kws={"label": "Spearman rho (log2FC)"}, ax=ax,
            square=True, linewidths=0.5)
ax.set_title("Fold-change concordance across Vitis stress datasets", fontsize=11)
fig.tight_layout()
fig.savefig(Path(args.out_dir) / "fc_correlation_heatmap.png", dpi=200)
plt.close(fig)


# -----------------------------------------------------------------------------
# 7. Per-stress signature-enrichment barplot (top summary)
# -----------------------------------------------------------------------------
if len(overlap_df):
    grbv_rows = overlap_df[(overlap_df["dataset_A"] == "GRBV") | (overlap_df["dataset_B"] == "GRBV")].copy()
    grbv_rows["other"] = np.where(grbv_rows["dataset_A"] == "GRBV",
                                  grbv_rows["dataset_B"], grbv_rows["dataset_A"])
    grbv_rows = grbv_rows.sort_values("odds_ratio", ascending=True)

    fig, ax = plt.subplots(figsize=(7, 3.2 + 0.3 * max(1, len(grbv_rows))))
    ax.barh(grbv_rows["other"], np.log2(grbv_rows["odds_ratio"].replace([np.inf, 0], np.nan)))
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("log2 odds ratio (GRBV vs dataset overlap)")
    ax.set_title("DEG overlap enrichment of GRBV signature with each comparator")
    for y, row in enumerate(grbv_rows.itertuples()):
        ax.text(0.02, y, f"q={row.fisher_q:.2e}", transform=ax.get_yaxis_transform(),
                va="center", fontsize=8, color="dimgray")
    fig.tight_layout()
    fig.savefig(Path(args.out_dir) / "grbv_vs_others_barplot.png", dpi=200)
    plt.close(fig)

print("[12] done.", file=sys.stderr)
