#!/usr/bin/env python
"""10. Characterise subtypes from rule 03.

For each NMF subtype within RB+ samples:
  - DE vs all other subtypes (Welch t-test, BH FDR)
  - Clinical associations (stage, season, YSI, titer)
Outputs subtype × metadata summary plus per-subtype marker gene list.
"""
import argparse, os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True); ap.add_argument("--meta", required=True)
ap.add_argument("--titer",  required=True); ap.add_argument("--subtypes", required=True)
ap.add_argument("--out_de", required=True); ap.add_argument("--out_clinical", required=True)
ap.add_argument("--out_plot", required=True)
args = ap.parse_args()
os.makedirs(os.path.dirname(args.out_de), exist_ok=True)

counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta   = pd.read_csv(args.meta,   sep="\t")
subs   = pd.read_csv(args.subtypes)

lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)

merged = subs.merge(meta, on="Sample", how="left")
samples = merged["Sample"].tolist()
subt    = merged["Subtype"].tolist()

# --------- DE per subtype vs rest ---------
de_all = []
genes = lcpm.index
for s in sorted(set(subt)):
    in_s  = [x for i, x in enumerate(samples) if subt[i] == s and x in lcpm.columns]
    out_s = [x for i, x in enumerate(samples) if subt[i] != s and x in lcpm.columns]
    if len(in_s) < 3 or len(out_s) < 3:
        continue
    a = lcpm[in_s].values
    b = lcpm[out_s].values
    t, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    fc = a.mean(axis=1) - b.mean(axis=1)
    q = multipletests(np.nan_to_num(p, nan=1), method="fdr_bh")[1]
    df = pd.DataFrame({"gene": genes, "subtype": s, "log2FC": fc,
                       "pvalue": p, "qvalue": q})
    df = df.dropna().sort_values("qvalue")
    de_all.append(df.head(1000))
de_out = pd.concat(de_all, ignore_index=True) if de_all else pd.DataFrame()
de_out.to_csv(args.out_de, index=False)

# --------- Clinical associations ---------
# Stage / season composition
tab_stage  = pd.crosstab(merged["Subtype"], merged["Maturity_stage"], normalize="index")
tab_season = pd.crosstab(merged["Subtype"], merged["Sampling_Season"], normalize="index")
# YSI mean per subtype
ysi_stats = merged.groupby("Subtype")["Years_Since_Infection"].agg(["mean","std","count"])

# Titer mean per subtype (use per-sample titer via metadata-titer join, same as 07/08)
titer_raw = pd.read_csv(args.titer, sep="\t")
rename = {"Pre_véraison":"Pre_veraison","Véraison":"Veraison",
          "Post_véraison":"Post_veraison","Harvest":"Harvest"}
long = titer_raw.melt(id_vars=["Season","Year","Biorep"],
                      value_vars=list(rename), var_name="stage_raw", value_name="titer")
long["Maturity_stage"] = long["stage_raw"].map(rename)
avg = long.groupby(["Season","Year","Maturity_stage"], as_index=False)["titer"].mean()
m2  = merged.merge(avg, left_on=["Sampling_Season","Years_Since_Infection","Maturity_stage"],
                         right_on=["Season","Year","Maturity_stage"], how="left")
titer_stats = m2.groupby("Subtype")["titer"].agg(["mean","std","count"]).rename(columns=lambda c: f"titer_{c}")

clinical = pd.concat([ysi_stats.rename(columns=lambda c: f"YSI_{c}"),
                      titer_stats,
                      tab_stage.add_prefix("stage_"),
                      tab_season.add_prefix("season_")], axis=1)
clinical.to_csv(args.out_clinical)

# --------- Overview plot ---------
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
sns.heatmap(tab_stage, annot=True, fmt=".2f", cmap="Blues", ax=axes[0]); axes[0].set_title("Stage composition")
sns.heatmap(tab_season, annot=True, fmt=".2f", cmap="Purples", ax=axes[1]); axes[1].set_title("Season composition")
merged.boxplot(column="Years_Since_Infection", by="Subtype", ax=axes[2], grid=False)
axes[2].set_title("YSI by subtype"); axes[2].set_xlabel("Subtype"); plt.suptitle("")
fig.tight_layout(); fig.savefig(args.out_plot, dpi=200); plt.close(fig)
print("Done.")
