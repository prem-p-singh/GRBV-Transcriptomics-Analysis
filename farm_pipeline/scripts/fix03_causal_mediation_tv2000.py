#!/usr/bin/env python
"""Fix #3: rerun the main-paper causal mediation analysis at tv=2000 (the
first threshold where the signal stabilises; the paper's default tv=500
finds ZERO significant mediators). Writes a drop-in replacement for
mediation_results.csv and a paper-ready summary.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import rankdata
from numpy.linalg import lstsq

ap = argparse.ArgumentParser()
ap.add_argument("--counts", default="/Users/prempratapsingh/Data_Analysis/counts_matrix.txt")
ap.add_argument("--meta",   default="/Users/prempratapsingh/Data_Analysis/sample_metadata.txt")
ap.add_argument("--titer",  default="/Users/prempratapsingh/Data_Analysis/Virus titer.txt")
ap.add_argument("--w2021",  default="/Users/prempratapsingh/Data_Analysis/daily_2021.csv")
ap.add_argument("--w2022",  default="/Users/prempratapsingh/Data_Analysis/daily_2022.csv")
ap.add_argument("--top_var", type=int, default=2000)
ap.add_argument("--out_full", default="/Users/prempratapsingh/Data_Analysis/mediation_results.csv")
ap.add_argument("--out_summary", default="/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/fix03_mediation/mediation_summary_tv2000.csv")
args = ap.parse_args()

Path(args.out_full).parent.mkdir(parents=True, exist_ok=True)
Path(args.out_summary).parent.mkdir(parents=True, exist_ok=True)

counts = pd.read_csv(args.counts, sep="\t", index_col=0)
keep = (counts >= 10).sum(axis=1) >= 0.10 * counts.shape[1]
counts = counts.loc[keep]
lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)

meta = pd.read_csv(args.meta, sep="\t").set_index("Sample")
meta = meta.loc[meta["Infection_Group"] == "RB(+)"]

sampling_dates = {
    (2021,"Pre_veraison"):"2021-08-02", (2021,"Veraison"):"2021-08-16",
    (2021,"Post_veraison"):"2021-09-08", (2021,"Harvest"):"2021-09-21",
    (2022,"Pre_veraison"):"2022-07-26", (2022,"Veraison"):"2022-08-09",
    (2022,"Post_veraison"):"2022-08-24", (2022,"Harvest"):"2022-09-13",
}
def cimis(p):
    w = pd.read_csv(p); w.columns=[c.strip() for c in w.columns]
    w["Date"] = pd.to_datetime(w["Date"], format="%m/%d/%Y", errors="coerce")
    w["Tavg"] = pd.to_numeric(w["Avg Air Temp (C)"], errors="coerce")
    return w.dropna(subset=["Date","Tavg"])
w_by = {2021: cimis(args.w2021), 2022: cimis(args.w2022)}
temp_map = {}
for (s, st), d in sampling_dates.items():
    wd = pd.Timestamp(d); w = w_by[s]
    mask = (w["Date"] >= wd - pd.Timedelta(days=7)) & (w["Date"] <= wd)
    if mask.any(): temp_map[(s, st)] = w.loc[mask, "Tavg"].mean()

titer = pd.read_csv(args.titer, sep="\t")
rn = {"Pre_véraison":"Pre_veraison","Véraison":"Veraison","Post_véraison":"Post_veraison","Harvest":"Harvest"}
lg = titer.melt(id_vars=["Season","Year","Biorep"], value_vars=list(rn),
                var_name="s_raw", value_name="titer")
lg["Maturity_stage"] = lg["s_raw"].map(rn)
titer_avg = lg.groupby(["Season","Year","Maturity_stage"], as_index=False)["titer"].mean()

m = meta.reset_index()
m["Temperature"] = m.apply(lambda r: temp_map.get((int(r["Sampling_Season"]), r["Maturity_stage"]), np.nan), axis=1)
m = m.merge(titer_avg, left_on=["Sampling_Season","Years_Since_Infection","Maturity_stage"],
            right_on=["Season","Year","Maturity_stage"], how="left").set_index("Sample")

samples = [s for s in m.index if s in lcpm.columns
           and not pd.isna(m.loc[s,"Temperature"]) and not pd.isna(m.loc[s,"titer"])]
print(f"[fix03] n={len(samples)} complete-data RB+ samples", file=sys.stderr)

top = lcpm[samples].var(axis=1).nlargest(args.top_var).index
T = m.loc[samples,"Temperature"].values; M = m.loc[samples,"titer"].values

def r1(X,y):
    X1=np.column_stack([np.ones(len(X)),X]); b,*_=lstsq(X1,y,rcond=None); p=X1@b
    n,k=X1.shape; s=np.sum((y-p)**2); mse=s/max(n-k,1e-10)
    se=np.sqrt(np.diag(mse*np.linalg.pinv(X1.T@X1)))
    t=b[1]/se[1] if se[1]>0 else 0; return float(b[1]), float(2*(1-stats.t.cdf(abs(t),n-k)))

def rm(X,y):
    X1=np.column_stack([np.ones(len(y)),X]); b,*_=lstsq(X1,y,rcond=None); p=X1@b
    n,k=X1.shape; s=np.sum((y-p)**2); mse=s/max(n-k,1e-10)
    try: se=np.sqrt(np.diag(mse*np.linalg.pinv(X1.T@X1)))
    except Exception: return b[1:], np.ones(k-1)
    t=b[1:]/np.where(se[1:]>0,se[1:],1e-10); return b[1:], 2*(1-stats.t.cdf(np.abs(t),n-k))

a, ap_ = r1(T, M)
print(f"[fix03] path a (T->titer): beta={a:.4f}, p={ap_:.4g}", file=sys.stderr)

rows = []
for g in top:
    Y = lcpm.loc[g, samples].values
    c, cp = r1(T, Y)
    bs, ps = rm(np.column_stack([T, M]), Y)
    rows.append({"gene": g, "total_c": c, "total_p": cp,
                 "direct_cprime": bs[0], "path_b": bs[1], "path_b_p": ps[1],
                 "indirect_ab": a*bs[1],
                 "prop_mediated": (a*bs[1]) / c if abs(c) > 1e-10 else 0})

df = pd.DataFrame(rows)
df["total_padj"]  = np.minimum(df["total_p"] * len(df) / rankdata(df["total_p"]), 1.0)
df["path_b_padj"] = np.minimum(df["path_b_p"] * len(df) / rankdata(df["path_b_p"]), 1.0)
df["sig_joint"]   = (df["total_padj"]<0.05) & (df["path_b_padj"]<0.05) & (df["indirect_ab"].abs()>0.05)
df = df.sort_values("path_b_padj")
df.to_csv(args.out_full, index=False)

sig = df[df["sig_joint"]]
print(f"[fix03] {len(sig)} significant mediators at BH FDR<0.05 + |ab|>0.05", file=sys.stderr)

summary = {
    "n_samples": len(samples), "top_var": args.top_var,
    "path_a_beta": a, "path_a_p": ap_,
    "n_significant_mediators": int(len(sig)),
    "top5_by_joint_p": ",".join(sig.nsmallest(5,"path_b_padj")["gene"].tolist()),
    "LOC100249249_in_sig": bool("LOC100249249" in sig["gene"].values),
    "LOC100853778_in_sig": bool("LOC100853778" in sig["gene"].values),
    "LOC100261367_in_sig": bool("LOC100261367" in sig["gene"].values),
    "prior_tv500_n_sig": 0,
}
pd.DataFrame([summary]).to_csv(args.out_summary, index=False)
print(f"[fix03] wrote {args.out_full}  and  {args.out_summary}", file=sys.stderr)
