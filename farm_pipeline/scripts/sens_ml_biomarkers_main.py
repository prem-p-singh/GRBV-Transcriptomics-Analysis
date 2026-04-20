#!/usr/bin/env python
"""Sensitivity wrapper for main-paper advanced_ml_biomarkers.py.

Replicates only the 5-method CV AUC and LASSO stability-selection, parameterised
by --top_var. Runs locally on the Mac on counts_matrix.txt + sample_metadata.txt.

Outputs one CSV per TV level with:
  - per-method 5-fold CV AUC (Logistic, Random Forest, LASSO, SVM-RBF, GB)
  - top-20 stable host-only genes (viral genes stripped) with LASSO selection freq
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import warnings; warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--counts",    default="/Users/prempratapsingh/Data_Analysis/counts_matrix.txt")
ap.add_argument("--meta",      default="/Users/prempratapsingh/Data_Analysis/sample_metadata.txt")
ap.add_argument("--top_var",   type=int, required=True)
ap.add_argument("--out_auc",   required=True)
ap.add_argument("--out_stab",  required=True)
args = ap.parse_args()

# Load
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
keep = (counts >= 10).sum(axis=1) >= 0.10 * counts.shape[1]
counts = counts.loc[keep]
lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)

meta = pd.read_csv(args.meta, sep="\t").set_index("Sample")
samples = [s for s in meta.index if s in lcpm.columns]
meta = meta.loc[samples]
y = (meta["Infection_Group"] == "RB(+)").astype(int).values

# Top variable genes
var_series = lcpm.var(axis=1).nlargest(args.top_var)
var_genes = var_series.index.tolist()
# Identify viral genes (N761_* from GRBV genome) so we can check host-only performance
viral_genes = [g for g in var_genes if g.startswith("N761_") or g.startswith("ViviM_")]
host_genes = [g for g in var_genes if g not in viral_genes]
print(f"[sens] tv={args.top_var}: {len(viral_genes)} viral + {len(host_genes)} host genes in pool",
      file=sys.stderr)

def build(panel):
    X = lcpm.loc[panel, samples].T.values
    return StandardScaler().fit_transform(X)

def cv_auc(X, y, clf, seed=42):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    aucs = []; probs = np.zeros(len(y))
    for tr, te in cv.split(X, y):
        clf.fit(X[tr], y[tr])
        probs[te] = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], probs[te]))
    return float(np.mean(aucs)), float(np.std(aucs)), float(roc_auc_score(y, probs))

classifiers = {
    "LR":      LogisticRegression(max_iter=2000, class_weight="balanced"),
    "RF":      RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42, n_jobs=-1),
    "LASSO":   LogisticRegression(penalty="l1", solver="saga", max_iter=3000, class_weight="balanced", C=0.1),
    "SVM-RBF": SVC(kernel="rbf", probability=True, class_weight="balanced"),
    "GB":      GradientBoostingClassifier(n_estimators=100, max_depth=3, random_state=42),
}

rows = []
for panel_name, panel in [("all_vars", var_genes), ("host_only", host_genes)]:
    X = build(panel)
    for name, clf in classifiers.items():
        mean_auc, sd_auc, overall = cv_auc(X, y, clf)
        rows.append({
            "top_var": args.top_var, "panel": panel_name, "method": name,
            "n_features": len(panel), "cv_mean_auc": mean_auc,
            "cv_std_auc": sd_auc, "overall_auc": overall,
        })
pd.DataFrame(rows).to_csv(args.out_auc, index=False)
print(f"[sens] wrote {args.out_auc}", file=sys.stderr)

# Stability selection on HOST-only panel at this tv
X_host = build(host_genes)
n_sub = 100
freq = np.zeros(len(host_genes))
for i in range(n_sub):
    rng = np.random.RandomState(i)
    idx = rng.choice(len(y), int(0.8 * len(y)), replace=False)
    lasso = LogisticRegression(penalty="l1", solver="saga", max_iter=2000,
                                class_weight="balanced", C=0.05, random_state=i)
    lasso.fit(X_host[idx], y[idx])
    freq += (np.abs(lasso.coef_[0]) > 0).astype(float)
freq /= n_sub

stab = pd.DataFrame({"gene": host_genes, "selection_freq": freq}).sort_values("selection_freq", ascending=False)
stab.insert(0, "top_var", args.top_var)
stab.to_csv(args.out_stab, index=False)
print(f"[sens] wrote {args.out_stab}: {(freq > 0.6).sum()} stable host genes (freq>0.6)",
      file=sys.stderr)
