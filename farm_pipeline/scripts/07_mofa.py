#!/usr/bin/env python
"""07. MOFA+ multi-omics integration.

Views:
  - expression (top variable log-CPM genes)
  - viral titer (log10)
  - weather aggregates (mean/max temp + GDD in N days prior to sampling stage)

Groups: Sampling_Season (2021 / 2022).
"""
import argparse, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mofapy2.run.entry_point import entry_point

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True); ap.add_argument("--meta", required=True)
ap.add_argument("--titer",  required=True)
ap.add_argument("--weather_2021", required=True); ap.add_argument("--weather_2022", required=True)
ap.add_argument("--out_factors", required=True); ap.add_argument("--out_varexp", required=True)
ap.add_argument("--out_model", required=True)
ap.add_argument("--n_factors", type=int, default=10); ap.add_argument("--n_iter", type=int, default=1000)
ap.add_argument("--top_var", type=int, default=3000)
args = ap.parse_args()

os.makedirs(os.path.dirname(args.out_factors), exist_ok=True)

# ----- Expression view -----
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta   = pd.read_csv(args.meta,   sep="\t")
lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
top = lcpm.var(axis=1).nlargest(args.top_var).index
expr = lcpm.loc[top, :].T                # samples × genes
expr.index = expr.index.astype(str)

# Align metadata to samples with expression
meta = meta.set_index("Sample").loc[expr.index.intersection(meta["Sample"]) if "Sample" in meta.columns else expr.index]
# fallback: meta already indexed by Sample if Sample column missing
meta = pd.read_csv(args.meta, sep="\t").set_index("Sample")
meta = meta.loc[meta.index.intersection(expr.index)]
expr = expr.loc[meta.index]

# ----- Titer view -----
titer_raw = pd.read_csv(args.titer, sep="\t")
# Long form: one value per (Season, Year, Biorep, stage)
stages = ["Pre_véraison", "Véraison", "Post_véraison", "Harvest"]
# Normalise stage strings to match metadata (Maturity_stage uses Pre_veraison etc.)
rename = {"Pre_véraison": "Pre_veraison", "Véraison": "Veraison",
          "Post_véraison": "Post_veraison", "Harvest": "Harvest"}
melt = titer_raw.melt(id_vars=["Season", "Year", "Biorep"],
                      value_vars=stages, var_name="Maturity_stage", value_name="titer")
melt["Maturity_stage"] = melt["Maturity_stage"].map(rename)
melt = melt.rename(columns={"Season": "Sampling_Season", "Year": "Years_Since_Infection"})
key = ["Sampling_Season", "Years_Since_Infection", "Maturity_stage"]
# Average replicates per condition (closest we can get without sample-level titer)
titer_avg = melt.groupby(key, as_index=False)["titer"].mean()
meta_key = meta.reset_index()[["Sample"] + key]
titer_per_sample = meta_key.merge(titer_avg, on=key, how="left").set_index("Sample")["titer"]
titer_per_sample = titer_per_sample.reindex(meta.index)
titer_view = pd.DataFrame({"log_titer": np.log10(titer_per_sample + 1)}, index=meta.index)

# ----- Weather view -----
def agg_weather(path, season):
    w = pd.read_csv(path)
    w.columns = [c.strip() for c in w.columns]
    # heuristic: find columns by keyword
    tmax = next((c for c in w.columns if "tmax" in c.lower() or "max" in c.lower()), None)
    tmin = next((c for c in w.columns if "tmin" in c.lower() or "min" in c.lower()), None)
    tavg = next((c for c in w.columns if "tavg" in c.lower() or "mean" in c.lower()), None)
    if tavg is None and tmax and tmin:
        w["_tavg"] = (w[tmax] + w[tmin]) / 2
        tavg = "_tavg"
    w["_gdd"] = np.clip(w[tavg] - 10, 0, None)
    return {"season_mean_tavg": float(w[tavg].mean()),
            "season_max_tmax": float(w[tmax].max()) if tmax else np.nan,
            "season_GDD": float(w["_gdd"].sum())}

w21 = agg_weather(args.weather_2021, 2021)
w22 = agg_weather(args.weather_2022, 2022)
wdf = pd.DataFrame([w21, w22], index=[2021, 2022])
weather_view = pd.DataFrame(
    [wdf.loc[int(s)] for s in meta["Sampling_Season"].astype(str).str.extract(r"(\d+)")[0].astype(int)],
    index=meta.index)

# ----- Build MOFA input (long format) -----
def to_long(df, view_name, group_col):
    recs = []
    for i, sample in enumerate(df.index):
        grp = str(meta.loc[sample, group_col])
        for feat in df.columns:
            val = df.loc[sample, feat]
            if pd.isna(val):
                continue
            recs.append((sample, grp, feat, view_name, float(val)))
    return pd.DataFrame(recs, columns=["sample", "group", "feature", "view", "value"])

long_df = pd.concat([
    to_long(expr,         "expression",  "Sampling_Season"),
    to_long(titer_view,   "titer",       "Sampling_Season"),
    to_long(weather_view, "weather",     "Sampling_Season"),
], ignore_index=True)

ent = entry_point()
ent.set_data_options(scale_views=True)
ent.set_data_df(long_df)
ent.set_model_options(factors=args.n_factors, spikeslab_weights=True, ard_weights=True)
ent.set_train_options(iter=args.n_iter, convergence_mode="medium", seed=42, verbose=True)
ent.build()
ent.run()
ent.save(args.out_model)

# ----- Extract factors & variance explained -----
Z = ent.model.getExpectations("Z")["E"]    # dict of group->matrix
factors_df = pd.concat(
    [pd.DataFrame(Z[g], index=long_df.loc[long_df["group"] == g, "sample"].unique(),
                  columns=[f"F{i+1}" for i in range(args.n_factors)]) for g in Z],
    axis=0)
factors_df.to_csv(args.out_factors)

var_exp = ent.model.calculate_variance_explained()
# var_exp: dict with keys 'r2_total' and 'r2_per_factor'; shapes vary by version. Plot per-view totals.
fig, ax = plt.subplots(figsize=(8, 5))
try:
    r2_pf = var_exp["r2_per_factor"]     # list of 2D arrays (factor × view) per group
    mat = np.stack([np.array(r2_pf[g]) for g in r2_pf], axis=0).mean(axis=0)
    views_ = list(ent.model.data.keys()) if hasattr(ent.model, "data") else ["expression","titer","weather"]
    im = ax.imshow(mat, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(len(views_))); ax.set_xticklabels(views_, rotation=30)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels([f"F{i+1}" for i in range(mat.shape[0])])
    fig.colorbar(im, ax=ax, label="R²")
except Exception as e:
    ax.text(0.5, 0.5, f"Var-exp plot failed: {e}", ha="center", va="center")
fig.tight_layout(); fig.savefig(args.out_varexp, dpi=200); plt.close(fig)
print("Done.")
