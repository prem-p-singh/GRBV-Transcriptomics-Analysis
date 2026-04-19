#!/usr/bin/env python
"""07. MOFA+ multi-omics integration.

Views (all per-sample):
  - expression (top-variable log-CPM)
  - titer      (log10 viral load)
  - weather    (window aggregates around each sample's sampling date)

Groups: Sampling_Season (2021 / 2022).
"""
import argparse, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mofapy2.run.entry_point import entry_point

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True); ap.add_argument("--meta", required=True)
ap.add_argument("--titer",  required=True); ap.add_argument("--weather", required=True)
ap.add_argument("--out_factors", required=True); ap.add_argument("--out_varexp", required=True)
ap.add_argument("--out_model", required=True)
ap.add_argument("--n_factors", type=int, default=10); ap.add_argument("--n_iter", type=int, default=1000)
ap.add_argument("--top_var", type=int, default=3000)
args = ap.parse_args()

os.makedirs(os.path.dirname(args.out_factors), exist_ok=True)

# ----- Expression view -----
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta   = pd.read_csv(args.meta,   sep="\t").set_index("Sample")
lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
top = lcpm.var(axis=1).nlargest(args.top_var).index
expr = lcpm.loc[top, :].T                     # samples × genes

samples = expr.index.intersection(meta.index)
expr = expr.loc[samples]; meta = meta.loc[samples]

# ----- Titer view (per-sample via season×year×stage average of bioreps) -----
titer = pd.read_csv(args.titer, sep="\t")
stages = ["Pre_véraison", "Véraison", "Post_véraison", "Harvest"]
rename = {"Pre_véraison": "Pre_veraison", "Véraison": "Veraison",
          "Post_véraison": "Post_veraison", "Harvest": "Harvest"}
melt = titer.melt(id_vars=["Season", "Year", "Biorep"],
                  value_vars=stages, var_name="stage_raw", value_name="titer")
melt["Maturity_stage"] = melt["stage_raw"].map(rename)
avg = melt.groupby(["Season", "Year", "Maturity_stage"], as_index=False)["titer"].mean()
key = meta.reset_index()[["Sample", "Sampling_Season", "Years_Since_Infection", "Maturity_stage"]]
key = key.rename(columns={"Sampling_Season": "Season", "Years_Since_Infection": "Year"})
per_sample_titer = key.merge(avg, on=["Season", "Year", "Maturity_stage"], how="left")\
                      .set_index("Sample")["titer"].reindex(samples)
titer_view = pd.DataFrame({"log_titer": np.log10(per_sample_titer + 1)}, index=samples)

# ----- Weather view (from rule 00) -----
weather_view = pd.read_csv(args.weather, index_col=0).reindex(samples)

# ----- Build MOFA long-format input -----
def to_long(df, view_name):
    recs = []
    for sample in df.index:
        grp = str(meta.loc[sample, "Sampling_Season"])
        for feat in df.columns:
            val = df.loc[sample, feat]
            if pd.isna(val):
                continue
            recs.append((sample, grp, feat, view_name, float(val)))
    return pd.DataFrame(recs, columns=["sample", "group", "feature", "view", "value"])

long_df = pd.concat([to_long(expr,         "expression"),
                     to_long(titer_view,   "titer"),
                     to_long(weather_view, "weather")], ignore_index=True)

ent = entry_point()
ent.set_data_options(scale_views=True)
ent.set_data_df(long_df)
ent.set_model_options(factors=args.n_factors, spikeslab_weights=True, ard_weights=True)
ent.set_train_options(iter=args.n_iter, convergence_mode="medium", seed=42, verbose=True)
ent.build()
ent.run()
ent.save(args.out_model)

# ----- Factors + variance explained -----
Z = ent.model.getExpectations("Z")["E"]
factors_df = pd.concat(
    [pd.DataFrame(Z[g], index=long_df.loc[long_df["group"] == g, "sample"].unique(),
                  columns=[f"F{i+1}" for i in range(args.n_factors)]) for g in Z],
    axis=0)
factors_df.to_csv(args.out_factors)

var_exp = ent.model.calculate_variance_explained()
fig, ax = plt.subplots(figsize=(8, 5))
try:
    r2_pf = var_exp["r2_per_factor"]
    mat = np.stack([np.array(r2_pf[g]) for g in r2_pf], axis=0).mean(axis=0)
    views_ = ["expression", "titer", "weather"]
    im = ax.imshow(mat, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(len(views_))); ax.set_xticklabels(views_, rotation=30)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels([f"F{i+1}" for i in range(mat.shape[0])])
    fig.colorbar(im, ax=ax, label="R²")
except Exception as e:
    ax.text(0.5, 0.5, f"var-exp plot failed: {e}", ha="center", va="center")
fig.tight_layout(); fig.savefig(args.out_varexp, dpi=200); plt.close(fig)
print("Done.")
