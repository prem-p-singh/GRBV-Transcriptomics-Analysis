#!/usr/bin/env python
"""13. Viral ORF dissection extensions.

Two analyses requested by the agenda audit (rec #1 extensions):

  A. Viral ORF expression as phenotype predictor.  Per-ORF log2(CPM+1) regressed
     against Harvest juice Brix (and TA, malic) on RB(+) Harvest samples,
     compared head-to-head with qPCR titer as predictor.  Reports adjusted R^2,
     AIC, and likelihood-ratio test of "ORF expression adds beyond titer".

  B. Per-ORF x WGCNA module correlation.  Spearman rho between each viral ORF's
     log2(CPM+1) and each of 13 module eigengenes (MEs.csv), restricted to
     RB(+) samples to avoid the trivial RB+ vs RB- artefact.  Identifies which
     host module each viral gene tracks.

Outputs:
  - orf_phenotype_prediction.csv   : per-ORF regression vs Brix/TA/malic
  - orf_vs_titer_LRT.csv           : LRT for "ORF beyond titer" per ORF, phenotype
  - orf_module_correlation.csv     : per-ORF x module Spearman matrix
  - orf_module_heatmap.png         : visual of B
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from openpyxl import load_workbook

ap = argparse.ArgumentParser()
ap.add_argument("--counts",   required=True)
ap.add_argument("--meta",     required=True)
ap.add_argument("--titer",    required=True)
ap.add_argument("--mes",      required=True)
ap.add_argument("--biochem",  required=True)
ap.add_argument("--out_dir",  required=True)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)

# ---- load -----------------------------------------------------------------
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta   = pd.read_csv(args.meta,   sep="\t")
mes    = pd.read_csv(args.mes)

# Viral ORFs: gene IDs starting with N761_gp
viral = [g for g in counts.index if g.startswith("N761_gp")]
print(f"Viral ORFs: {viral}")

# log2(CPM+1) per sample
lib = counts.sum(axis=0)
cpm = counts.div(lib, axis=1) * 1e6
logcpm = np.log2(cpm + 1)
viral_expr = logcpm.loc[viral].T   # samples x ORFs
viral_expr.index.name = "Sample"

# ---- load titer (wide -> long, then merge per Sample) ---------------------
ti = pd.read_csv(args.titer, sep="\t")
# Columns: Season, Year, Biorep, Pre_véraison, Véraison, Post_véraison, Harvest
stage_cols = ["Pre_véraison", "Véraison", "Post_véraison", "Harvest"]
ti_long = ti.melt(id_vars=["Season","Year","Biorep"], value_vars=stage_cols,
                  var_name="Stage_raw", value_name="titer")
stage_norm = {"Pre_véraison":"Pre_veraison","Véraison":"Veraison",
              "Post_véraison":"Post_veraison","Harvest":"Harvest"}
ti_long["Maturity_stage"] = ti_long["Stage_raw"].map(stage_norm)
ti_long["Sampling_Season"] = ti_long["Season"].astype(str)
ti_long["Years_Since_Infection"] = ti_long["Year"]
ti_long["biorep_id"] = ti_long["Biorep"]

# Coerce shared keys to consistent dtypes (numeric-first, then to int-string)
for df_ in (meta, ti_long):
    df_["Sampling_Season"] = pd.to_numeric(df_["Sampling_Season"], errors="coerce").astype("Int64").astype(str)
    df_["Years_Since_Infection"] = pd.to_numeric(df_["Years_Since_Infection"], errors="coerce").astype("Int64")
    df_["biorep_id"] = pd.to_numeric(df_["biorep_id"], errors="coerce").astype("Int64")
merged = meta.merge(
    ti_long[["Sampling_Season","Years_Since_Infection","Maturity_stage","biorep_id","titer"]],
    on=["Sampling_Season","Years_Since_Infection","Maturity_stage","biorep_id"],
    how="left",
)
merged["log10_titer"] = np.log10(merged["titer"].clip(lower=1))

# ---- load tidied biochem (Harvest only) -----------------------------------
# Use already-prepped biochemistry_harvest.csv (sample_id parser lives in the
# repo-root biochemistry_harvest_prep.py).  Schema:
#   Sample, Sampling_Season, Infection_Group, Cohort, Infection_Year,
#   YSI_years, Biorep, Brix, pH, TA, Malic, NOPA, NH3, YAN
bc = pd.read_csv(args.biochem)
print(f"Biochem rows: {len(bc)}; columns: {bc.columns.tolist()[:10]}")

bc_clean = bc.rename(columns={"Biorep":"biorep_id"})
bc_clean["Sampling_Season"] = bc_clean["Sampling_Season"].astype(str)
# Match transcriptome's Years_Since_Infection (= cohort onset year) to biochem's Infection_Year
bc_clean["Years_Since_Infection"] = pd.to_numeric(bc_clean["Infection_Year"], errors="coerce").astype("Int64")
bc_clean["biorep_id"] = pd.to_numeric(bc_clean["biorep_id"], errors="coerce").astype("Int64")
keep = ["Sampling_Season","Years_Since_Infection","biorep_id","Infection_Group","Brix"]
if "TA"    in bc_clean.columns: keep.append("TA")
if "Malic" in bc_clean.columns: keep.append("Malic")

merged_h = merged[merged["Maturity_stage"] == "Harvest"].copy()
merged_h["biorep_id"] = pd.to_numeric(merged_h["biorep_id"], errors="coerce").astype("Int64")
merged_h["Years_Since_Infection"] = pd.to_numeric(merged_h["Years_Since_Infection"], errors="coerce").astype("Int64")

mh = merged_h.merge(bc_clean[keep], on=["Sampling_Season","Years_Since_Infection","biorep_id","Infection_Group"], how="inner")
mh = mh.merge(viral_expr, left_on="Sample", right_index=True, how="left")
print(f"Harvest+biochem+viral merged rows: {len(mh)}")

# ---- A. ORF as phenotype predictor (vs titer) -----------------------------
phenotypes = [c for c in ["Brix","TA","Malic"] if c in mh.columns]
results_A, lrt_rows = [], []
for orf in viral:
    for pheno in phenotypes:
        sub = mh[[orf, "log10_titer", pheno]].dropna()
        if len(sub) < 10: continue
        # Model 1: pheno ~ log10_titer
        X1 = sm.add_constant(sub[["log10_titer"]])
        m1 = sm.OLS(sub[pheno], X1).fit()
        # Model 2: pheno ~ log10_titer + ORF
        X2 = sm.add_constant(sub[["log10_titer", orf]])
        m2 = sm.OLS(sub[pheno], X2).fit()
        # Model 3: pheno ~ ORF (alone)
        X3 = sm.add_constant(sub[[orf]])
        m3 = sm.OLS(sub[pheno], X3).fit()
        # LRT m1 vs m2
        lrt = -2 * (m1.llf - m2.llf)
        p_lrt = stats.chi2.sf(lrt, df=1)
        results_A.append({
            "ORF": orf, "phenotype": pheno, "n": len(sub),
            "R2_titer_only": m1.rsquared, "AIC_titer_only": m1.aic,
            "R2_orf_only":   m3.rsquared, "AIC_orf_only":   m3.aic,
            "R2_titer_plus_orf": m2.rsquared, "AIC_titer_plus_orf": m2.aic,
            "beta_orf_in_full": m2.params.get(orf, np.nan),
            "p_orf_in_full":    m2.pvalues.get(orf, np.nan),
            "LRT_p_orf_added":  p_lrt,
        })

resA = pd.DataFrame(results_A)
resA.to_csv(os.path.join(args.out_dir, "orf_phenotype_prediction.csv"), index=False)
print(f"A. {len(resA)} (ORF, phenotype) rows -> orf_phenotype_prediction.csv")

# ---- B. Per-ORF x module correlation (RB+ only) --------------------------
mes_cols = [c for c in mes.columns if c.startswith("ME") and c != "MEgrey"]
mes_idx = mes.set_index("sample")
viral_idx = viral_expr.copy()
viral_idx.index = viral_idx.index.astype(str)
common = mes_idx.index.intersection(viral_idx.index)
print(f"Samples shared between MEs and viral expr: {len(common)}")

# Restrict to RB(+).  MEs.csv has 'group' field.
mes_rb = mes_idx.loc[common]
rb_pos_samples = mes_rb[mes_rb["group"].astype(str).str.startswith("RB.")].index
ve_rb = viral_idx.loc[rb_pos_samples]
me_rb = mes_idx.loc[rb_pos_samples, mes_cols]

corr_rows = []
for orf in viral:
    for me in mes_cols:
        s = pd.concat([ve_rb[orf], me_rb[me]], axis=1).dropna()
        if len(s) < 10: continue
        rho, p = stats.spearmanr(s.iloc[:, 0], s.iloc[:, 1])
        corr_rows.append({"ORF": orf, "module": me, "n": len(s), "rho": rho, "p": p})
corr = pd.DataFrame(corr_rows)
# BH-adjust per ORF
corr["padj"] = np.nan
for orf in viral:
    sel = corr["ORF"] == orf
    if sel.sum() == 0: continue
    pv = corr.loc[sel, "p"].values
    order = np.argsort(pv)
    ranks = np.empty_like(order); ranks[order] = np.arange(1, len(pv)+1)
    padj = pv * len(pv) / ranks
    padj = np.minimum.accumulate(padj[order[::-1]])[::-1]
    out = np.empty_like(padj); out[order] = padj
    corr.loc[sel, "padj"] = np.clip(out, 0, 1)

corr.to_csv(os.path.join(args.out_dir, "orf_module_correlation.csv"), index=False)
print(f"B. {len(corr)} (ORF, module) rows -> orf_module_correlation.csv")

# Heatmap
piv = corr.pivot(index="ORF", columns="module", values="rho")
fig, ax = plt.subplots(figsize=(8, 4))
sns.heatmap(piv, cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6, ax=ax,
            annot=True, fmt=".2f", cbar_kws={"label":"Spearman rho"})
ax.set_title("Per-ORF x module correlation (RB+ only)")
plt.tight_layout()
plt.savefig(os.path.join(args.out_dir, "orf_module_heatmap.png"), dpi=150)
plt.close(fig)
print("Done.")
