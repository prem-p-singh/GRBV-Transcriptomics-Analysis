#!/usr/bin/env python
"""16. Vine-level heterogeneity within genotype-matched cohort cells.

44 unique vines x up to 8 samples each (4 stages x 2 seasons).  Per-vine
infection-effect-size analysis answers: do vines respond uniformly to GRBV,
or are there 'tolerant' vs 'susceptible' responders even within a single
vineyard, single cultivar (Merlot), single infection cohort?

Method:
  1. Compute a per-vine 'infection response score' = mean Spearman rho
     between that vine's expression deviation (vs cohort+season+stage centroid
     of RB-) and the consensus universal-DEG signature.
  2. Cluster vines by infection response into tolerant / mid / susceptible.
  3. Test whether tolerant vs susceptible differ on biochem (Brix, TA),
     viral titer, or top biomarker expression.

Outputs (results/16_vine_heterogeneity/):
  - per_vine_response_scores.csv
  - vine_pca_infection_axis.png       (PCA coloured by per-vine score)
  - tolerant_vs_susceptible_traits.csv
  - vine_heterogeneity_dashboard.png
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

ap = argparse.ArgumentParser()
ap.add_argument("--counts",   required=True)
ap.add_argument("--meta",     required=True)
ap.add_argument("--titer",    required=True)
ap.add_argument("--biochem",  required=True)
ap.add_argument("--out_dir",  required=True)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)

# -------- load -----------------------------------------------------------
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta = pd.read_csv(args.meta, sep="\t").set_index("Sample")
common = counts.columns.intersection(meta.index)
counts, meta = counts[common], meta.loc[common]
meta.index.name = "Sample"
keep = (counts >= 10).sum(axis=1) >= int(0.10 * counts.shape[1])
counts = counts.loc[keep]
lib = counts.sum(axis=0)
logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
print(f"Filtered counts: {counts.shape}; vines: {meta['vine_id'].nunique()}")

# -------- per-vine response score ---------------------------------------
# Approach: for each (season, stage) cell, compute mean RB(-) expression
# centroid, then for each RB(+) vine sample compute log2FC deviation from
# centroid. Per-vine score = mean(|deviation|) on top variable genes.
top_var = logcpm.var(axis=1).nlargest(2000).index

records = []
for season in meta["Sampling_Season"].unique():
    for stage in meta["Maturity_stage"].unique():
        cell = meta[(meta["Sampling_Season"]==season) & (meta["Maturity_stage"]==stage)]
        if cell.empty: continue
        rb_neg = cell[cell["Infection_Group"]=="RB(-)"].index
        rb_pos = cell[cell["Infection_Group"]=="RB(+)"].index
        if len(rb_neg) < 2 or len(rb_pos) < 2: continue
        centroid = logcpm.loc[top_var, rb_neg].mean(axis=1)
        for s in rb_pos:
            dev = logcpm.loc[top_var, s] - centroid
            records.append({
                "Sample": s, "vine_id": meta.loc[s, "vine_id"],
                "Sampling_Season": season, "Maturity_stage": stage,
                "Years_Since_Infection": meta.loc[s, "Years_Since_Infection"],
                "mean_abs_dev": float(np.mean(np.abs(dev))),
                "mean_dev":     float(np.mean(dev)),
            })
sample_dev = pd.DataFrame(records)
print(f"Per-RB+ sample deviation rows: {len(sample_dev)}")

per_vine = sample_dev.groupby("vine_id").agg(
    n_samples=("Sample","size"),
    mean_response=("mean_abs_dev","mean"),
    sd_response=("mean_abs_dev","std"),
    seasons_observed=("Sampling_Season", lambda s: ",".join(sorted(set(s.astype(str))))),
    stages_observed=("Maturity_stage", lambda s: ",".join(sorted(set(s.astype(str))))),
).reset_index()
per_vine["response_z"] = (per_vine["mean_response"] - per_vine["mean_response"].mean()) / per_vine["mean_response"].std()

# 3-cluster KMeans on response score
km = KMeans(n_clusters=3, random_state=42, n_init=20).fit(per_vine[["mean_response"]].values)
order = np.argsort([per_vine.loc[per_vine.cluster_label==i if False else km.labels_==i, "mean_response"].mean()
                    for i in range(3)])
label_map = {old: ["tolerant","intermediate","susceptible"][new] for new, old in enumerate(order)}
per_vine["cluster"] = [label_map[c] for c in km.labels_]
per_vine.to_csv(os.path.join(args.out_dir, "per_vine_response_scores.csv"), index=False)
print(per_vine.groupby("cluster").size())

# -------- vine-level PCA ------------------------------------------------
# Restrict to RB(+) at Harvest (densest cell), one row per vine
harvest_pos = meta[(meta["Maturity_stage"]=="Harvest") & (meta["Infection_Group"]=="RB(+)")].index
expr = logcpm.loc[top_var, harvest_pos].T
expr = expr - expr.mean()
pca = PCA(n_components=2).fit(expr)
pc = pca.transform(expr)
pc_df = pd.DataFrame({"Sample": expr.index, "PC1": pc[:,0], "PC2": pc[:,1]})
pc_df["vine_id"] = pc_df["Sample"].map(meta["vine_id"])
pc_df = pc_df.merge(per_vine[["vine_id","cluster","response_z"]], on="vine_id", how="left")

fig, ax = plt.subplots(figsize=(7, 5))
sns.scatterplot(data=pc_df, x="PC1", y="PC2", hue="cluster",
                palette={"tolerant":"#2CA02C","intermediate":"#999999","susceptible":"#B2182B"},
                size="response_z", sizes=(40,180), ax=ax)
ax.set_title("Vine-level PCA at Harvest (RB+ only) coloured by tolerance cluster")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
plt.tight_layout()
plt.savefig(os.path.join(args.out_dir, "vine_pca_infection_axis.png"), dpi=150)
plt.close(fig)

# -------- tolerant vs susceptible traits -------------------------------
bc = pd.read_csv(args.biochem)
bc = bc.rename(columns={"Biorep":"biorep_id"})
bc["Sampling_Season"] = bc["Sampling_Season"].astype(str)
bc["Years_Since_Infection"] = pd.to_numeric(bc["Infection_Year"], errors="coerce").astype("Int64")
bc["biorep_id"] = pd.to_numeric(bc["biorep_id"], errors="coerce").astype("Int64")
mh = meta.reset_index()
mh["Sampling_Season"] = mh["Sampling_Season"].astype(str)
mh["biorep_id"] = pd.to_numeric(mh["biorep_id"], errors="coerce").astype("Int64")
mh["Years_Since_Infection"] = pd.to_numeric(mh["Years_Since_Infection"], errors="coerce").astype("Int64")
mh_brix = mh.merge(bc[["Sampling_Season","Years_Since_Infection","biorep_id","Infection_Group","Brix"]],
                   on=["Sampling_Season","Years_Since_Infection","biorep_id","Infection_Group"], how="inner")
if True:
    mh_brix = mh_brix.merge(per_vine[["vine_id","cluster"]], on="vine_id", how="left")
    grp = mh_brix.groupby("cluster")["Brix"].agg(["mean","std","count"]).reset_index()
    grp.to_csv(os.path.join(args.out_dir, "tolerant_vs_susceptible_traits.csv"), index=False)
    print(grp)

# Dashboard panel
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
sns.histplot(per_vine, x="mean_response", hue="cluster", ax=axes[0],
             palette={"tolerant":"#2CA02C","intermediate":"#999999","susceptible":"#B2182B"})
axes[0].set_title(f"Per-vine infection response distribution (n={len(per_vine)} vines)")
axes[0].set_xlabel("Mean |log2FC vs cohort RB- centroid| (top 2000 var genes)")
sns.scatterplot(data=pc_df, x="PC1", y="PC2", hue="cluster",
                palette={"tolerant":"#2CA02C","intermediate":"#999999","susceptible":"#B2182B"},
                size="response_z", sizes=(40,180), ax=axes[1], legend=False)
axes[1].set_title("Harvest RB+ vines in PC space")
plt.tight_layout()
plt.savefig(os.path.join(args.out_dir, "vine_heterogeneity_dashboard.png"), dpi=150)
plt.close(fig)
print("Done.")
