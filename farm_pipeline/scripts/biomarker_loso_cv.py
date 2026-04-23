#!/usr/bin/env python
"""Fix #4 + #13: rerun ML biomarker analysis with (a) host-only as the
primary reporting panel, (b) standard random 5-fold CV AND (c) leave-one-
infection-cohort-out (LOCO) CV as the proper generalisation test.

Writes two CSVs that should replace the paper's advanced_ml_comparison.csv
and stability_selection_results.csv equivalents.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold, LeaveOneGroupOut
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import warnings; warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--counts", default="/Users/prempratapsingh/Data_Analysis/counts_matrix.txt")
ap.add_argument("--meta",   default="/Users/prempratapsingh/Data_Analysis/sample_metadata.txt")
ap.add_argument("--top_var", type=int, default=500)  # smaller pool than 200 to avoid overfit
ap.add_argument("--out_auc", default="/Users/prempratapsingh/Data_Analysis/advanced_ml_comparison_FIX04.csv")
ap.add_argument("--out_stab", default="/Users/prempratapsingh/Data_Analysis/stability_selection_FIX04.csv")
args = ap.parse_args()

counts = pd.read_csv(args.counts, sep="\t", index_col=0)
keep = (counts >= 10).sum(axis=1) >= 0.10 * counts.shape[1]
counts = counts.loc[keep]
lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
meta = pd.read_csv(args.meta, sep="\t").set_index("Sample")
samples = [s for s in meta.index if s in lcpm.columns]
meta = meta.loc[samples]; y = (meta["Infection_Group"] == "RB(+)").astype(int).values

var_genes = lcpm.var(axis=1).nlargest(args.top_var).index.tolist()
viral = [g for g in var_genes if g.startswith("N761_") or g.startswith("ViviM_")]
host  = [g for g in var_genes if g not in viral]
print(f"[biomarker] tv={args.top_var}: {len(viral)} viral + {len(host)} host features", file=sys.stderr)

def scale(X): return StandardScaler().fit_transform(X)
def feats(panel): return scale(lcpm.loc[panel, samples].T.values)

def random_cv(X, y, clf, seed=42):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    aucs, probs = [], np.zeros(len(y))
    for tr, te in cv.split(X, y):
        clf.fit(X[tr], y[tr]); probs[te] = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], probs[te]))
    return float(np.mean(aucs)), float(np.std(aucs)), float(roc_auc_score(y, probs))

def loco_cv(X, y, clf, groups):
    """Leave-one-infection-cohort-out (LOCO). Groups are the Years_Since_Infection column."""
    cv = LeaveOneGroupOut()
    aucs, probs = [], np.full(len(y), np.nan)
    for tr, te in cv.split(X, y, groups):
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            continue
        clf.fit(X[tr], y[tr]); probs[te] = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], probs[te]))
    overall = roc_auc_score(y[~np.isnan(probs)], probs[~np.isnan(probs)]) if np.any(~np.isnan(probs)) else np.nan
    return (float(np.mean(aucs)) if aucs else np.nan,
            float(np.std(aucs))  if aucs else np.nan,
            float(overall))

classifiers = {
    "LR":      LogisticRegression(max_iter=2000, class_weight="balanced"),
    "RF":      RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42, n_jobs=-1),
    "LASSO":   LogisticRegression(penalty="l1", solver="saga", max_iter=3000, class_weight="balanced", C=0.1),
    "SVM-RBF": SVC(kernel="rbf", probability=True, class_weight="balanced"),
    "GB":      GradientBoostingClassifier(n_estimators=100, max_depth=3, random_state=42),
}

# Infection cohort = the *year the vine got infected* (sampling_season - YSI + 1).
# Controls (RB-) share a single bucket "control" so each fold contains both classes.
year_infected = (meta["Sampling_Season"].astype(int) - meta["Years_Since_Infection"].astype(int)).astype(str)
year_infected.loc[meta["Infection_Group"] == "RB(-)"] = "control"
groups = year_infected.values
print(f"[biomarker] LOCO groups: {np.unique(groups, return_counts=True)}", file=sys.stderr)

# Alternate grouping: leave-one-sampling-season-out (2 folds: 2021, 2022)
season_groups = meta["Sampling_Season"].astype(str).values

rows = []
for panel_name, panel in [("all_vars", var_genes), ("host_only", host)]:
    X = feats(panel)
    for name, clf in classifiers.items():
        r_mean, r_sd, r_overall = random_cv(X, y, clf)
        l_mean, l_sd, l_overall = loco_cv(X, y, clf, groups)
        s_mean, s_sd, s_overall = loco_cv(X, y, clf, season_groups)
        rows.append({
            "panel": panel_name, "method": name, "n_features": len(panel),
            "random_5fold_mean_auc": r_mean, "random_5fold_std": r_sd,
            "random_5fold_overall_auc": r_overall,
            "loco_cohort_mean_auc": l_mean, "loco_cohort_std": l_sd,
            "loco_cohort_overall_auc": l_overall,
            "loco_season_mean_auc": s_mean, "loco_season_std": s_sd,
            "loco_season_overall_auc": s_overall,
        })
auc_df = pd.DataFrame(rows)
auc_df.to_csv(args.out_auc, index=False)
print("[biomarker] AUC table:")
print(auc_df.to_string(index=False))

# Stability selection on HOST-only
X_host = feats(host)
n_sub = 100
freq = np.zeros(len(host))
for i in range(n_sub):
    rng = np.random.RandomState(i); idx = rng.choice(len(y), int(0.8*len(y)), replace=False)
    l = LogisticRegression(penalty="l1", solver="saga", max_iter=2000,
                            class_weight="balanced", C=0.05, random_state=i)
    l.fit(X_host[idx], y[idx]); freq += (np.abs(l.coef_[0]) > 0).astype(float)
freq /= n_sub
stab = pd.DataFrame({"panel":"host_only", "gene": host, "lasso_selection_freq": freq})\
          .sort_values("lasso_selection_freq", ascending=False)
stab.to_csv(args.out_stab, index=False)
print(f"[biomarker] stable host genes (freq>=0.6): {(freq>=0.6).sum()}")
