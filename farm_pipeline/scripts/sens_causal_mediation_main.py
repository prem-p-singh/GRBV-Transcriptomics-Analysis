#!/usr/bin/env python
"""Sensitivity wrapper for the main-paper causal_mediation_analysis.py.

Replicates only the mediation-regression core, parametrised by --top_var.
Runs locally on the Mac, reads the same inputs the main script uses.

Mediation direction tested (as in the main paper): Temperature -> Titer -> Gene.
  path a: Temperature -> Titer
  path b: Titer -> Gene (controlling for Temperature)
  path c: Temperature -> Gene (total effect)
  indirect = a * b
"""
from __future__ import annotations
import argparse, sys, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata
from numpy.linalg import lstsq

ap = argparse.ArgumentParser()
ap.add_argument("--counts",      default="/Users/prempratapsingh/Data_Analysis/counts_matrix.txt")
ap.add_argument("--meta",        default="/Users/prempratapsingh/Data_Analysis/sample_metadata.txt")
ap.add_argument("--titer",       default="/Users/prempratapsingh/Data_Analysis/Virus titer.txt")
ap.add_argument("--w2021",       default="/Users/prempratapsingh/Data_Analysis/daily_2021.csv")
ap.add_argument("--w2022",       default="/Users/prempratapsingh/Data_Analysis/daily_2022.csv")
ap.add_argument("--top_var",     type=int, required=True)
ap.add_argument("--out",         required=True)
args = ap.parse_args()

# ---- Load expression ----
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
keep = (counts >= 10).sum(axis=1) >= 0.10 * counts.shape[1]
counts = counts.loc[keep]
lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)

# ---- Metadata ----
meta = pd.read_csv(args.meta, sep="\t").set_index("Sample")
meta = meta.loc[meta["Infection_Group"] == "RB(+)"]

# ---- Weather: 7-day-mean temp per (season, stage) ----
sampling_dates = {
    (2021, "Pre_veraison"): "2021-08-02",
    (2021, "Veraison"):     "2021-08-16",
    (2021, "Post_veraison"):"2021-09-08",
    (2021, "Harvest"):      "2021-09-21",
    (2022, "Pre_veraison"): "2022-07-26",
    (2022, "Veraison"):     "2022-08-09",
    (2022, "Post_veraison"):"2022-08-24",
    (2022, "Harvest"):      "2022-09-13",
}

def load_cimis(p):
    w = pd.read_csv(p)
    w.columns = [c.strip() for c in w.columns]
    w["Date"] = pd.to_datetime(w["Date"], format="%m/%d/%Y", errors="coerce")
    w["Tavg"] = pd.to_numeric(w["Avg Air Temp (C)"], errors="coerce")
    return w.dropna(subset=["Date","Tavg"])

w_by_season = {2021: load_cimis(args.w2021), 2022: load_cimis(args.w2022)}
temp_map = {}
for (season, stage), date_str in sampling_dates.items():
    d = pd.Timestamp(date_str)
    w = w_by_season[season]
    mask = (w["Date"] >= d - pd.Timedelta(days=7)) & (w["Date"] <= d)
    if mask.any():
        temp_map[(season, stage)] = w.loc[mask, "Tavg"].mean()

# ---- Titer: average across bioreps for each (Season, Year, Stage) ----
titer_raw = pd.read_csv(args.titer, sep="\t")
rename = {"Pre_véraison":"Pre_veraison", "Véraison":"Veraison",
          "Post_véraison":"Post_veraison", "Harvest":"Harvest"}
tl = titer_raw.melt(id_vars=["Season","Year","Biorep"],
                    value_vars=list(rename), var_name="stage_raw", value_name="titer")
tl["Maturity_stage"] = tl["stage_raw"].map(rename)
titer_avg = tl.groupby(["Season","Year","Maturity_stage"], as_index=False)["titer"].mean()

# ---- Merge metadata with per-sample temperature + titer ----
m = meta.reset_index()
m["Temperature"] = m.apply(
    lambda r: temp_map.get((int(r["Sampling_Season"]), r["Maturity_stage"]), np.nan), axis=1)
m = m.merge(titer_avg,
            left_on=["Sampling_Season","Years_Since_Infection","Maturity_stage"],
            right_on=["Season","Year","Maturity_stage"], how="left")
m = m.set_index("Sample")

# Samples for full analysis: temp AND titer present AND in counts matrix
samples_full = [s for s in m.index if s in lcpm.columns and
                not pd.isna(m.loc[s, "Temperature"]) and
                not pd.isna(m.loc[s, "titer"])]
print(f"[sens] n={len(samples_full)} complete-data samples", file=sys.stderr)

# ---- Top-variable gene filter ----
top = lcpm[samples_full].var(axis=1).nlargest(args.top_var).index
print(f"[sens] tested {len(top)} top-variable genes", file=sys.stderr)

T = m.loc[samples_full, "Temperature"].values
M = m.loc[samples_full, "titer"].values

def reg1(X, y):
    X1 = np.column_stack([np.ones(len(X)), X])
    b, *_ = lstsq(X1, y, rcond=None)
    pred = X1 @ b
    n, p = X1.shape
    ssr = np.sum((y - pred) ** 2)
    mse = ssr / max(n - p, 1e-10)
    se = np.sqrt(np.diag(mse * np.linalg.pinv(X1.T @ X1)))
    t = b[1] / se[1] if se[1] > 0 else 0.0
    pv = 2 * (1 - stats.t.cdf(abs(t), n - p))
    return float(b[1]), float(pv)

def reg_multi(X, y):
    X1 = np.column_stack([np.ones(len(y)), X])
    b, *_ = lstsq(X1, y, rcond=None)
    pred = X1 @ b
    n, p = X1.shape
    ssr = np.sum((y - pred) ** 2)
    mse = ssr / max(n - p, 1e-10)
    try:
        se = np.sqrt(np.diag(mse * np.linalg.pinv(X1.T @ X1)))
    except Exception:
        return b[1:], np.ones(p - 1)
    t = b[1:] / np.where(se[1:] > 0, se[1:], 1e-10)
    pv = 2 * (1 - stats.t.cdf(np.abs(t), n - p))
    return b[1:], pv

a_coef, a_pval = reg1(T, M)
print(f"[sens] path a (T -> titer): beta={a_coef:.4f}, p={a_pval:.4g}", file=sys.stderr)

rows = []
for g in top:
    Y = lcpm.loc[g, samples_full].values
    c_coef, c_pval = reg1(T, Y)
    betas, pvals = reg_multi(np.column_stack([T, M]), Y)
    c_prime = betas[0]; b_coef = betas[1]; b_pval = pvals[1]
    indirect = a_coef * b_coef
    rows.append({"gene": g, "total_c": c_coef, "total_p": c_pval,
                 "direct_cprime": c_prime, "path_b": b_coef, "path_b_p": b_pval,
                 "indirect_ab": indirect,
                 "prop_mediated": indirect / c_coef if abs(c_coef) > 1e-10 else 0})

df = pd.DataFrame(rows)
df["total_padj"] = np.minimum(df["total_p"] * len(df) / rankdata(df["total_p"]), 1.0)
df["path_b_padj"] = np.minimum(df["path_b_p"] * len(df) / rankdata(df["path_b_p"]), 1.0)
df["sig_mediation"] = (df["path_b_padj"] < 0.05) & (df["total_padj"] < 0.05) & (df["indirect_ab"].abs() > 0.05)
df = df.sort_values("path_b_padj")
df.to_csv(args.out, index=False)
print(f"[sens] wrote {args.out}: {len(df)} genes tested; {df['sig_mediation'].sum()} pass joint significance",
      file=sys.stderr)
