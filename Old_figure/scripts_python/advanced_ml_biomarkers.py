#!/usr/bin/env python3
"""
Advanced Multi-method ML Biomarker Comparison for GRBV Transcriptomics.

Extends biomarker_discovery.py with LASSO, SVM-RBF, Gradient Boosting,
stability selection, and class-imbalance handling to validate the AUC=1.0
finding and identify a consensus biomarker panel.

Outputs:
  - Fig 46: Multi-method AUC comparison (bar + CI)
  - Fig 47: Stability selection frequency plot
  - Fig 48: Consensus panel Venn + ROC
  - advanced_ml_comparison.csv
  - stability_selection_results.csv
  - consensus_biomarker_panel.csv
"""

import pandas as pd
import numpy as np
import sqlite3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import rankdata
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'figure.dpi': 300, 'font.size': 15, 'axes.titlesize': 16,
    'axes.labelsize': 15, 'figure.facecolor': 'white'
})

# ============================================================
# 1. LOAD DATA
# ============================================================
print("=" * 70)
print("ADVANCED ML BIOMARKERS: LOADING DATA")
print("=" * 70)

counts = pd.read_csv('counts_matrix.txt', sep='\t', index_col=0)
meta = pd.read_csv('sample_metadata.txt', sep='\t')

common_samples = list(set(counts.columns) & set(meta['Sample']))
counts = counts[common_samples]
meta = meta[meta['Sample'].isin(common_samples)].set_index('Sample').loc[common_samples]

min_samples = int(0.1 * counts.shape[1])
gene_mask = (counts >= 10).sum(axis=1) >= min_samples
counts_filtered = counts.loc[gene_mask]

lib_sizes = counts_filtered.sum(axis=0)
cpm = counts_filtered.div(lib_sizes, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

# Labels: RB(+) = 1, RB(-) = 0
y = (meta['Infection_Group'] == 'RB(+)').astype(int).values
sample_names = meta.index.tolist()

print(f"Genes: {counts_filtered.shape[0]}, Samples: {counts_filtered.shape[1]}")
print(f"Class distribution: RB(+)={y.sum()}, RB(-)={(y == 0).sum()}")

# Load gene names
conn = sqlite3.connect('org.Vvinifera.eg.db/inst/extdata/org.Vvinifera.eg.sqlite')
gene_info = pd.read_sql("SELECT SYMBOL, GENENAME FROM gene_info", conn)
gene_name_map = dict(zip(gene_info['SYMBOL'], gene_info['GENENAME']))
conn.close()


# ============================================================
# 2. FEATURE SELECTION — Top 200 most variable genes
# ============================================================
print("\n" + "=" * 70)
print("FEATURE SELECTION")
print("=" * 70)

# Use top 200 most variable genes as candidate features
var_genes = log_cpm.var(axis=1).nlargest(200).index.tolist()
X = log_cpm.loc[var_genes].T.values  # (327, 200)
gene_names_arr = np.array(var_genes)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"Feature matrix: {X_scaled.shape}")


# ============================================================
# 3. MULTI-METHOD COMPARISON
# ============================================================
print("\n" + "=" * 70)
print("MULTI-METHOD AUC COMPARISON")
print("=" * 70)

# Define classifiers
classifiers = {
    'Logistic Regression': LogisticRegression(
        max_iter=2000, class_weight='balanced', C=1.0, solver='lbfgs'),
    'Random Forest': RandomForestClassifier(
        n_estimators=200, class_weight='balanced', random_state=42, n_jobs=-1),
    'LASSO (L1)': LogisticRegression(
        penalty='l1', solver='saga', max_iter=3000, class_weight='balanced', C=0.1),
    'SVM-RBF': SVC(
        kernel='rbf', probability=True, class_weight='balanced', C=1.0),
    'Gradient Boosting': GradientBoostingClassifier(
        n_estimators=100, max_depth=3, learning_rate=0.1, random_state=42),
}

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
method_results = {}

for name, clf in classifiers.items():
    print(f"\n--- {name} ---")

    # 5-fold CV AUC
    fold_aucs = []
    fold_preds = np.zeros(len(y))
    fold_probs = np.zeros(len(y))
    feature_importances = np.zeros(len(var_genes))

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X_scaled, y)):
        X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        clf.fit(X_train, y_train)
        y_prob = clf.predict_proba(X_test)[:, 1]
        fold_auc = roc_auc_score(y_test, y_prob)
        fold_aucs.append(fold_auc)
        fold_probs[test_idx] = y_prob

        # Feature importance
        if hasattr(clf, 'feature_importances_'):
            feature_importances += clf.feature_importances_
        elif hasattr(clf, 'coef_'):
            feature_importances += np.abs(clf.coef_[0])

    feature_importances /= 5

    mean_auc = np.mean(fold_aucs)
    std_auc = np.std(fold_aucs)
    overall_auc = roc_auc_score(y, fold_probs)

    method_results[name] = {
        'mean_auc': mean_auc,
        'std_auc': std_auc,
        'overall_auc': overall_auc,
        'fold_aucs': fold_aucs,
        'probs': fold_probs,
        'importances': feature_importances,
    }

    print(f"  5-fold CV AUC: {mean_auc:.4f} ± {std_auc:.4f}")
    print(f"  Overall AUC: {overall_auc:.4f}")

    # Top 10 features
    top_idx = np.argsort(feature_importances)[::-1][:10]
    print(f"  Top features:")
    for i in top_idx:
        gname = gene_name_map.get(var_genes[i], '')[:40]
        print(f"    {var_genes[i]:20s} ({gname:40s}) imp={feature_importances[i]:.4f}")


# ============================================================
# 4. STABILITY SELECTION
# ============================================================
print("\n" + "=" * 70)
print("STABILITY SELECTION (100 subsamples)")
print("=" * 70)

n_subsamples = 100
n_features_select = 20
gene_selection_freq = np.zeros(len(var_genes))

for i in range(n_subsamples):
    # Random subsample: 80% of samples
    np.random.seed(i)
    idx = np.random.choice(len(y), size=int(0.8 * len(y)), replace=False)
    X_sub = X_scaled[idx]
    y_sub = y[idx]

    # LASSO feature selection
    lasso = LogisticRegression(penalty='l1', solver='saga', max_iter=2000,
                                class_weight='balanced', C=0.05, random_state=i)
    lasso.fit(X_sub, y_sub)
    selected = np.abs(lasso.coef_[0]) > 0
    gene_selection_freq += selected

gene_selection_freq /= n_subsamples

# Stable genes: selected in >60% of subsamples
stable_mask = gene_selection_freq > 0.6
stable_genes = gene_names_arr[stable_mask]
stable_freqs = gene_selection_freq[stable_mask]

print(f"Stable genes (>60% selection): {len(stable_genes)}")
for g, f in sorted(zip(stable_genes, stable_freqs), key=lambda x: -x[1]):
    gname = gene_name_map.get(g, '')[:40]
    print(f"  {g:20s} ({gname:40s}) freq={f:.2f}")

# Save stability results
stab_df = pd.DataFrame({
    'gene': var_genes,
    'selection_frequency': gene_selection_freq,
    'gene_name': [gene_name_map.get(g, '') for g in var_genes],
    'stable': gene_selection_freq > 0.6,
})
stab_df = stab_df.sort_values('selection_frequency', ascending=False)
stab_df.to_csv('stability_selection_results.csv', index=False)
print(f"Saved stability_selection_results.csv")


# ============================================================
# 5. CONSENSUS PANEL
# ============================================================
print("\n" + "=" * 70)
print("CONSENSUS BIOMARKER PANEL")
print("=" * 70)

# Genes selected by ≥3 methods in top 20
method_top20 = {}
for name, res in method_results.items():
    top20 = gene_names_arr[np.argsort(res['importances'])[::-1][:20]]
    method_top20[name] = set(top20)

gene_method_count = defaultdict(int)
gene_method_list = defaultdict(list)
for name, genes in method_top20.items():
    for g in genes:
        gene_method_count[g] += 1
        gene_method_list[g].append(name)

# Consensus: in ≥3 methods AND stable
consensus_genes = []
for g, count in sorted(gene_method_count.items(), key=lambda x: -x[1]):
    is_stable = g in stable_genes
    consensus_genes.append({
        'gene': g,
        'gene_name': gene_name_map.get(g, ''),
        'n_methods': count,
        'methods': ';'.join(gene_method_list[g]),
        'stability_freq': gene_selection_freq[list(var_genes).index(g)] if g in var_genes else 0,
        'is_stable': is_stable,
        'consensus': count >= 3 or (count >= 2 and is_stable),
    })

consensus_df = pd.DataFrame(consensus_genes)
consensus_df.to_csv('consensus_biomarker_panel.csv', index=False)
n_consensus = consensus_df['consensus'].sum()
print(f"Consensus panel: {n_consensus} genes (≥3 methods OR ≥2 methods + stable)")
print(f"Saved consensus_biomarker_panel.csv")

# Evaluate consensus panel AUC
consensus_gene_list = consensus_df[consensus_df['consensus']]['gene'].tolist()
if len(consensus_gene_list) >= 2:
    cons_idx = [list(var_genes).index(g) for g in consensus_gene_list if g in var_genes]
    X_cons = X_scaled[:, cons_idx]

    lr = LogisticRegression(max_iter=2000, class_weight='balanced')
    cons_aucs = cross_val_score(lr, X_cons, y, cv=cv, scoring='roc_auc')
    print(f"Consensus panel AUC (5-fold CV): {cons_aucs.mean():.4f} ± {cons_aucs.std():.4f}")


# ============================================================
# 6. SAVE METHOD COMPARISON
# ============================================================
comp_rows = []
for name, res in method_results.items():
    comp_rows.append({
        'method': name,
        'mean_auc': res['mean_auc'],
        'std_auc': res['std_auc'],
        'overall_auc': res['overall_auc'],
        'fold_1': res['fold_aucs'][0],
        'fold_2': res['fold_aucs'][1],
        'fold_3': res['fold_aucs'][2],
        'fold_4': res['fold_aucs'][3],
        'fold_5': res['fold_aucs'][4],
    })
comp_df = pd.DataFrame(comp_rows)
comp_df.to_csv('advanced_ml_comparison.csv', index=False)
print(f"Saved advanced_ml_comparison.csv")


# ============================================================
# 7. FIGURE 46 — Multi-method AUC Comparison
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 46: MULTI-METHOD AUC COMPARISON")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(20, 8))
fig.suptitle('Multi-Method ML Biomarker Comparison',
             fontsize=20, fontweight='bold', y=0.98)

# Panel A: AUC bar chart with CI
ax = axes[0]
methods = list(method_results.keys())
mean_aucs = [method_results[m]['mean_auc'] for m in methods]
std_aucs = [method_results[m]['std_auc'] for m in methods]

colors = ['#e74c3c', '#2ecc71', '#3498db', '#9b59b6', '#f39c12']
bars = ax.bar(range(len(methods)), mean_aucs, yerr=std_aucs,
              color=colors[:len(methods)], alpha=0.85, edgecolor='white',
              capsize=5, linewidth=1.5)

ax.set_xticks(range(len(methods)))
ax.set_xticklabels(methods, fontsize=11, rotation=25, ha='right')
ax.set_ylabel('AUC-ROC (5-fold CV)', fontsize=14)
ax.set_title('A. Classification Performance', fontsize=16, fontweight='bold')
ax.set_ylim(0.5, 1.05)
ax.axhline(0.5, color='grey', linestyle='--', linewidth=0.5, label='Random')

for i, (m, s) in enumerate(zip(mean_aucs, std_aucs)):
    ax.text(i, m + s + 0.01, f'{m:.3f}', ha='center', fontsize=11, fontweight='bold')

# Panel B: ROC curves
ax = axes[1]
for i, (name, res) in enumerate(method_results.items()):
    fpr, tpr, _ = roc_curve(y, res['probs'])
    ax.plot(fpr, tpr, color=colors[i], linewidth=2,
            label=f"{name} (AUC={res['overall_auc']:.3f})")

ax.plot([0, 1], [0, 1], 'k--', linewidth=0.5, alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=14)
ax.set_ylabel('True Positive Rate', fontsize=14)
ax.set_title('B. ROC Curves (Overall)', fontsize=16, fontweight='bold')
ax.legend(fontsize=10, loc='lower right')

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('46_ml_method_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 46_ml_method_comparison.png")


# ============================================================
# 8. FIGURE 47 — Stability Selection Frequency Plot
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 47: STABILITY SELECTION")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(20, 9))
fig.suptitle('Stability Selection Analysis',
             fontsize=20, fontweight='bold', y=0.98)

# Panel A: Frequency histogram
ax = axes[0]
sorted_freq = np.sort(gene_selection_freq)[::-1]
ax.bar(range(len(sorted_freq)), sorted_freq, color='#3498db', alpha=0.7, width=1.0)
ax.axhline(0.6, color='#e74c3c', linestyle='--', linewidth=2, label='Stability threshold (60%)')
ax.set_xlabel('Gene rank', fontsize=14)
ax.set_ylabel('Selection frequency', fontsize=14)
ax.set_title('A. LASSO Selection Frequency (100 subsamples)', fontsize=16, fontweight='bold')
ax.legend(fontsize=12)

n_stable = (gene_selection_freq > 0.6).sum()
ax.text(0.95, 0.95, f'{n_stable} stable genes',
        transform=ax.transAxes, ha='right', va='top', fontsize=14,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Panel B: Top 20 genes with labels
ax = axes[1]
top20_idx = np.argsort(gene_selection_freq)[::-1][:20]
top20_genes = gene_names_arr[top20_idx]
top20_freqs = gene_selection_freq[top20_idx]

colors_bar = ['#e74c3c' if f > 0.6 else '#95a5a6' for f in top20_freqs]
bars = ax.barh(range(len(top20_genes)), top20_freqs, color=colors_bar, alpha=0.85)
ax.set_yticks(range(len(top20_genes)))
labels = [f"{g} ({gene_name_map.get(g, '')[:30]})" for g in top20_genes]
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel('Selection frequency', fontsize=14)
ax.set_title('B. Top 20 Most Stable Genes', fontsize=16, fontweight='bold')
ax.invert_yaxis()
ax.axvline(0.6, color='#e74c3c', linestyle='--', linewidth=1.5)

for i, v in enumerate(top20_freqs):
    ax.text(v + 0.01, i, f'{v:.2f}', va='center', fontsize=10)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('47_stability_selection.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 47_stability_selection.png")


# ============================================================
# 9. FIGURE 48 — Consensus Panel
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 48: CONSENSUS PANEL")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(20, 9))
fig.suptitle('Consensus Biomarker Panel',
             fontsize=20, fontweight='bold', y=0.98)

# Panel A: Method overlap matrix
ax = axes[0]
method_names = list(method_top20.keys())
overlap_mat = np.zeros((len(method_names), len(method_names)), dtype=int)
for i, m1 in enumerate(method_names):
    for j, m2 in enumerate(method_names):
        overlap_mat[i, j] = len(method_top20[m1] & method_top20[m2])

short_names = [m.split(' ')[0] if len(m.split(' ')) <= 2 else m[:12] for m in method_names]
sns.heatmap(overlap_mat, ax=ax, annot=True, fmt='d', cmap='Blues',
            xticklabels=short_names, yticklabels=short_names,
            linewidths=0.5, linecolor='white')
ax.set_title('A. Top-20 Gene Overlap Between Methods', fontsize=16, fontweight='bold')

# Panel B: Consensus genes ROC
ax = axes[1]
if len(consensus_gene_list) >= 2:
    cons_idx = [list(var_genes).index(g) for g in consensus_gene_list if g in var_genes]
    X_cons = X_scaled[:, cons_idx]

    # Full model ROC
    lr = LogisticRegression(max_iter=2000, class_weight='balanced')
    lr.fit(X_cons, y)
    y_prob_full = lr.predict_proba(X_cons)[:, 1]
    fpr_full, tpr_full, _ = roc_curve(y, y_prob_full)
    auc_full = roc_auc_score(y, y_prob_full)

    # CV ROC
    cv_probs = np.zeros(len(y))
    for train_idx, test_idx in cv.split(X_cons, y):
        lr_cv = LogisticRegression(max_iter=2000, class_weight='balanced')
        lr_cv.fit(X_cons[train_idx], y[train_idx])
        cv_probs[test_idx] = lr_cv.predict_proba(X_cons[test_idx])[:, 1]

    fpr_cv, tpr_cv, _ = roc_curve(y, cv_probs)
    auc_cv = roc_auc_score(y, cv_probs)

    ax.plot(fpr_full, tpr_full, 'b-', linewidth=2,
            label=f'Training (AUC={auc_full:.3f})')
    ax.plot(fpr_cv, tpr_cv, 'r-', linewidth=2,
            label=f'5-fold CV (AUC={auc_cv:.3f})')
    ax.plot([0, 1], [0, 1], 'k--', linewidth=0.5)
    ax.set_xlabel('False Positive Rate', fontsize=14)
    ax.set_ylabel('True Positive Rate', fontsize=14)
    ax.set_title(f'B. Consensus Panel ROC ({len(consensus_gene_list)} genes)',
                fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)

    # List consensus genes
    gene_text = '\n'.join([f"  {g}" for g in consensus_gene_list[:8]])
    if len(consensus_gene_list) > 8:
        gene_text += f"\n  +{len(consensus_gene_list)-8} more"
    ax.text(0.55, 0.25, gene_text, transform=ax.transAxes,
            fontsize=9, fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('48_consensus_panel.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 48_consensus_panel.png")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("ADVANCED ML BIOMARKER ANALYSIS COMPLETE")
print("=" * 70)

print("\nMethod AUC summary:")
for name, res in method_results.items():
    print(f"  {name:25s}: AUC={res['mean_auc']:.4f} ± {res['std_auc']:.4f}")

print(f"\nStable genes (>60% LASSO selection): {n_stable}")
print(f"Consensus panel: {n_consensus} genes")

print(f"\nOutput files:")
print(f"  - 46_ml_method_comparison.png")
print(f"  - 47_stability_selection.png")
print(f"  - 48_consensus_panel.png")
print(f"  - advanced_ml_comparison.csv")
print(f"  - stability_selection_results.csv")
print(f"  - consensus_biomarker_panel.csv")
