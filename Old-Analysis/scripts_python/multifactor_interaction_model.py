#!/usr/bin/env python3
"""
Multi-factor Interaction Model for GRBV Transcriptomics.

Tests Infection × Season × Stage interactions using per-gene ANOVA to formally
quantify whether the inverted seasonal trajectory is a statistically significant
interaction effect.

Outputs:
  - Fig 40: Variance decomposition (stacked bar, eta-squared per factor)
  - Fig 41: Interaction term gene sets (UpSet-style / Venn)
  - Fig 42: Infection × Season interaction expression profiles
  - multifactor_anova_results.csv
  - interaction_gene_categories.csv
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import rankdata
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
print("MULTI-FACTOR INTERACTION MODEL: LOADING DATA")
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

print(f"Genes: {counts_filtered.shape[0]}, Samples: {counts_filtered.shape[1]}")

# Encode factors
meta['Infection'] = (meta['Infection_Group'] == 'RB(+)').astype(int)
meta['Season'] = meta['Sampling_Season'].astype(int)

# Encode maturity stage as ordered numeric
stage_order = {'Pre_veraison': 0, 'Veraison': 1, 'Post_veraison': 2, 'Harvest': 3}
meta['Stage_num'] = meta['Maturity_stage'].map(stage_order)

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
stage_labels = {'Pre_veraison': 'Pre-vér', 'Veraison': 'Véraison',
                'Post_veraison': 'Post-vér', 'Harvest': 'Harvest'}

print(f"\nInfection: {meta['Infection'].value_counts().to_dict()}")
print(f"Seasons: {meta['Season'].value_counts().to_dict()}")
print(f"Stages: {meta['Maturity_stage'].value_counts().to_dict()}")

# ============================================================
# 2. PER-GENE ANOVA (Type II)
# ============================================================
print("\n" + "=" * 70)
print("RUNNING PER-GENE ANOVA (Infection × Season × Stage)")
print("=" * 70)

# Use a simplified but robust ANOVA approach:
# For each gene, fit OLS: expression ~ Infection + Season + Stage + I×S + I×St + S×St + I×S×St
# Using dummy-coded factors and computing sum of squares

try:
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    from statsmodels.stats.anova import anova_lm
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("WARNING: statsmodels not installed. Using scipy-based ANOVA.")

# Build design DataFrame
design = meta[['Infection', 'Season', 'Maturity_stage']].copy()
design['Season'] = design['Season'].astype(str)
design['Infection'] = design['Infection'].map({1: 'Infected', 0: 'Control'})

anova_results = []
n_genes = len(log_cpm.index)

print(f"Testing {n_genes} genes...")

for idx, gene in enumerate(log_cpm.index):
    if idx % 2000 == 0:
        print(f"  Processing gene {idx}/{n_genes}...")

    expr = log_cpm.loc[gene].values
    df = design.copy()
    df['expr'] = expr

    if HAS_STATSMODELS:
        try:
            model = ols('expr ~ C(Infection) * C(Season) * C(Maturity_stage)', data=df).fit()
            aov = anova_lm(model, typ=2)

            # Extract sum of squares for eta-squared calculation
            ss_total = aov['sum_sq'].sum()
            result = {'gene': gene}

            for term in aov.index:
                if term == 'Residual':
                    continue
                clean_term = term.replace('C(', '').replace(')', '')
                clean_term = clean_term.replace('Infection', 'Inf').replace('Season', 'Sea').replace('Maturity_stage', 'Stg')
                result[f'F_{clean_term}'] = aov.loc[term, 'F']
                result[f'p_{clean_term}'] = aov.loc[term, 'PR(>F)']
                result[f'eta2_{clean_term}'] = aov.loc[term, 'sum_sq'] / ss_total

            anova_results.append(result)
        except Exception:
            pass
    else:
        # Fallback: simple factorial ANOVA using scipy
        # Test main effects and key interactions
        result = {'gene': gene}

        # Main effect: Infection
        groups = [expr[design['Infection'] == g] for g in design['Infection'].unique()]
        if all(len(g) > 1 for g in groups):
            f, p = stats.f_oneway(*groups)
            result['F_Inf'] = f
            result['p_Inf'] = p

        # Main effect: Season
        groups = [expr[design['Season'] == g] for g in design['Season'].unique()]
        if all(len(g) > 1 for g in groups):
            f, p = stats.f_oneway(*groups)
            result['F_Sea'] = f
            result['p_Sea'] = p

        # Main effect: Stage
        groups = [expr[design['Maturity_stage'] == g] for g in design['Maturity_stage'].unique()]
        if all(len(g) > 1 for g in groups):
            f, p = stats.f_oneway(*groups)
            result['F_Stg'] = f
            result['p_Stg'] = p

        anova_results.append(result)

anova_df = pd.DataFrame(anova_results)
print(f"ANOVA completed for {len(anova_df)} genes")

# BH correction for each term
p_cols = [c for c in anova_df.columns if c.startswith('p_')]
for pc in p_cols:
    padj_col = pc.replace('p_', 'padj_')
    valid = anova_df[pc].dropna()
    if len(valid) > 0:
        anova_df.loc[valid.index, padj_col] = np.minimum(
            valid.values * len(valid) / rankdata(valid.values), 1.0)

# Save
anova_df.to_csv('multifactor_anova_results.csv', index=False)
print(f"Saved multifactor_anova_results.csv ({len(anova_df)} genes)")

# ============================================================
# 3. CATEGORIZE GENES BY SIGNIFICANT TERMS
# ============================================================
print("\n" + "=" * 70)
print("CATEGORIZING GENES BY SIGNIFICANT INTERACTION TERMS")
print("=" * 70)

padj_thresh = 0.05
categories = {}

# Identify genes significant for each term
padj_cols = [c for c in anova_df.columns if c.startswith('padj_')]
term_names = [c.replace('padj_', '') for c in padj_cols]

for term, padj_col in zip(term_names, padj_cols):
    if padj_col in anova_df.columns:
        sig = anova_df[anova_df[padj_col] < padj_thresh]
        categories[term] = set(sig['gene'].tolist())
        print(f"  {term}: {len(sig)} significant genes (padj < {padj_thresh})")

# Key interaction categories
interaction_terms = [t for t in term_names if ':' in t]
main_terms = [t for t in term_names if ':' not in t]

# Genes with significant interactions
interaction_genes = []
for gene_idx, row in anova_df.iterrows():
    gene = row['gene']
    sig_terms = []
    for term, padj_col in zip(term_names, padj_cols):
        if padj_col in anova_df.columns and pd.notna(row.get(padj_col)) and row[padj_col] < padj_thresh:
            sig_terms.append(term)
    if sig_terms:
        interaction_genes.append({
            'gene': gene,
            'significant_terms': ';'.join(sig_terms),
            'n_sig_terms': len(sig_terms),
            'has_interaction': any(':' in t for t in sig_terms),
        })

interaction_df = pd.DataFrame(interaction_genes)
interaction_df.to_csv('interaction_gene_categories.csv', index=False)
print(f"\nSaved interaction_gene_categories.csv ({len(interaction_df)} genes with ≥1 significant term)")

if len(interaction_df) > 0:
    n_any_interaction = interaction_df['has_interaction'].sum()
    print(f"Genes with significant interaction effect: {n_any_interaction}")


# ============================================================
# 4. FIGURE 40 — Variance Decomposition
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 40: VARIANCE DECOMPOSITION")
print("=" * 70)

eta2_cols = [c for c in anova_df.columns if c.startswith('eta2_')]
if len(eta2_cols) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    fig.suptitle('Variance Decomposition: Infection × Season × Stage',
                 fontsize=20, fontweight='bold', y=0.98)

    # Panel A: Mean eta-squared across all genes
    ax = axes[0]
    mean_eta2 = anova_df[eta2_cols].mean()
    clean_names = [c.replace('eta2_', '').replace(':', '×') for c in eta2_cols]

    # Sort by effect size
    order = mean_eta2.values.argsort()[::-1]
    sorted_names = [clean_names[i] for i in order]
    sorted_vals = mean_eta2.values[order]

    colors_bar = ['#e74c3c' if '×' not in n else '#3498db' if n.count('×') == 1 else '#2ecc71'
                  for n in sorted_names]

    bars = ax.barh(range(len(sorted_names)), sorted_vals, color=colors_bar, alpha=0.85, edgecolor='white')
    ax.set_yticks(range(len(sorted_names)))
    ax.set_yticklabels(sorted_names, fontsize=12)
    ax.set_xlabel('Mean η² (proportion of variance)', fontsize=14)
    ax.set_title('A. Mean Variance Explained per Factor', fontsize=16, fontweight='bold')
    ax.invert_yaxis()

    for i, v in enumerate(sorted_vals):
        ax.text(v + 0.001, i, f'{v:.3f}', va='center', fontsize=11)

    # Legend
    from matplotlib.patches import Patch
    ax.legend([Patch(color='#e74c3c'), Patch(color='#3498db'), Patch(color='#2ecc71')],
              ['Main effects', '2-way interactions', '3-way interaction'],
              loc='lower right', fontsize=11)

    # Panel B: Stacked bar — top genes with largest interaction effects
    ax = axes[1]

    # Find interaction eta2 column
    int_cols = [c for c in eta2_cols if ':' in c]
    if int_cols:
        # Sum all interaction eta2 for each gene
        anova_df['total_interaction_eta2'] = anova_df[int_cols].sum(axis=1)
        top_int = anova_df.nlargest(20, 'total_interaction_eta2')

        # Stacked bar of eta2 components for top interaction genes
        bottom = np.zeros(len(top_int))
        color_map = {}
        palette = plt.cm.Set2(np.linspace(0, 1, len(eta2_cols)))
        for i, col in enumerate(eta2_cols):
            clean = col.replace('eta2_', '').replace(':', '×')
            color_map[col] = palette[i]
            vals = top_int[col].values
            ax.barh(range(len(top_int)), vals, left=bottom, color=palette[i],
                    label=clean, alpha=0.85, edgecolor='white', linewidth=0.3)
            bottom += vals

        ax.set_yticks(range(len(top_int)))
        ax.set_yticklabels([g[:20] for g in top_int['gene']], fontsize=9)
        ax.set_xlabel('η² (proportion of variance)', fontsize=14)
        ax.set_title('B. Top 20 Genes by Interaction Effect', fontsize=16, fontweight='bold')
        ax.invert_yaxis()
        ax.legend(fontsize=9, loc='lower right', ncol=2)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig('40_variance_decomposition.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved 40_variance_decomposition.png")


# ============================================================
# 5. FIGURE 41 — Interaction Term Gene Sets
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 41: INTERACTION TERM GENE SETS")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(20, 9))
fig.suptitle('Genes Significant for Each ANOVA Term',
             fontsize=20, fontweight='bold', y=0.98)

# Panel A: Bar chart of significant gene counts per term
ax = axes[0]
term_counts = {}
for term in term_names:
    padj_col = f'padj_{term}'
    if padj_col in anova_df.columns:
        n = (anova_df[padj_col] < padj_thresh).sum()
        term_counts[term.replace(':', '×')] = n

if term_counts:
    names = list(term_counts.keys())
    vals = list(term_counts.values())
    colors_bar2 = ['#e74c3c' if '×' not in n else '#3498db' if n.count('×') == 1 else '#2ecc71'
                   for n in names]

    ax.barh(range(len(names)), vals, color=colors_bar2, alpha=0.85, edgecolor='white')
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=12)
    ax.set_xlabel('Number of Significant Genes (padj < 0.05)', fontsize=14)
    ax.set_title('A. Significant Genes per ANOVA Term', fontsize=16, fontweight='bold')
    ax.invert_yaxis()

    for i, v in enumerate(vals):
        ax.text(v + max(vals) * 0.01, i, str(v), va='center', fontsize=11)

# Panel B: Overlap heatmap between interaction terms
ax = axes[1]

if len(categories) >= 2:
    cat_names = sorted(categories.keys(), key=lambda x: len(categories[x]), reverse=True)[:8]
    overlap_matrix = pd.DataFrame(0, index=cat_names, columns=cat_names)

    for i, t1 in enumerate(cat_names):
        for j, t2 in enumerate(cat_names):
            s1 = categories.get(t1, set())
            s2 = categories.get(t2, set())
            if i == j:
                overlap_matrix.iloc[i, j] = len(s1)
            else:
                overlap_matrix.iloc[i, j] = len(s1 & s2)

    display_names = [n.replace(':', '×') for n in cat_names]
    sns.heatmap(overlap_matrix.values.astype(int), ax=ax,
                xticklabels=display_names, yticklabels=display_names,
                annot=True, fmt='d', cmap='YlOrRd',
                linewidths=0.5, linecolor='white')
    ax.set_title('B. Gene Overlap Between Terms', fontsize=16, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    ax.tick_params(axis='y', rotation=0)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('41_interaction_gene_sets.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 41_interaction_gene_sets.png")


# ============================================================
# 6. FIGURE 42 — Infection × Season Interaction Profiles
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 42: INFECTION × SEASON INTERACTION PROFILES")
print("=" * 70)

# Find genes with the strongest Infection × Season interaction
inf_sea_col = None
for c in anova_df.columns:
    if 'padj' in c and 'Inf' in c and 'Sea' in c and 'Stg' not in c:
        inf_sea_col = c
        break

if inf_sea_col is None:
    # Fallback: look for any interaction column
    for c in anova_df.columns:
        if 'padj' in c and ':' in c:
            inf_sea_col = c
            break

if inf_sea_col is not None:
    fig, axes = plt.subplots(3, 4, figsize=(24, 18))
    fig.suptitle('Top Genes with Infection × Season Interaction Effects',
                 fontsize=20, fontweight='bold', y=0.98)

    # Get top 12 genes by interaction p-value
    top_genes = anova_df.nsmallest(12, inf_sea_col)

    for idx, (_, row) in enumerate(top_genes.iterrows()):
        ax = axes[idx // 4, idx % 4]
        gene = row['gene']

        # Get expression data for this gene across conditions
        plot_data = []
        for season in [2021, 2022]:
            for stage in stages:
                for inf_status in ['RB(+)', 'RB(-)']:
                    mask = (meta['Sampling_Season'] == season) & \
                           (meta['Maturity_stage'] == stage) & \
                           (meta['Infection_Group'] == inf_status)
                    samples = meta[mask].index
                    expr = log_cpm.loc[gene, samples].values
                    for v in expr:
                        plot_data.append({
                            'Season': str(season),
                            'Stage': stage_labels.get(stage, stage),
                            'Infection': inf_status,
                            'Expression': v
                        })

        plot_df = pd.DataFrame(plot_data)

        # Line plot: mean expression across stages, by season × infection
        for season, ls in [('2021', '-'), ('2022', '--')]:
            for inf, color in [('RB(+)', '#e74c3c'), ('RB(-)', '#2980b9')]:
                sub = plot_df[(plot_df['Season'] == season) & (plot_df['Infection'] == inf)]
                if len(sub) == 0:
                    continue
                means = sub.groupby('Stage')['Expression'].mean()
                sems = sub.groupby('Stage')['Expression'].sem()
                stage_order_labels = [stage_labels[s] for s in stages]
                ordered_means = [means.get(s, np.nan) for s in stage_order_labels]
                ordered_sems = [sems.get(s, np.nan) for s in stage_order_labels]

                ax.errorbar(range(4), ordered_means, yerr=ordered_sems,
                           marker='o', linestyle=ls, color=color,
                           label=f"{season} {inf}", markersize=5, linewidth=1.5,
                           capsize=3)

        padj_val = row[inf_sea_col]
        ax.set_title(f"{gene[:18]}\npadj={padj_val:.2e}", fontsize=12)
        ax.set_xticks(range(4))
        ax.set_xticklabels([stage_labels[s] for s in stages], fontsize=9, rotation=30)
        ax.set_ylabel('log2(CPM+1)', fontsize=11)

        if idx == 0:
            ax.legend(fontsize=7, loc='best')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig('42_interaction_profiles.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved 42_interaction_profiles.png")
else:
    print("No interaction terms found in ANOVA results")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("MULTI-FACTOR ANOVA COMPLETE")
print("=" * 70)

if len(eta2_cols) > 0:
    print("\nMean η² across all genes:")
    for col in eta2_cols:
        clean = col.replace('eta2_', '').replace(':', '×')
        print(f"  {clean}: {anova_df[col].mean():.4f}")

print(f"\nSignificant genes per term (padj < {padj_thresh}):")
for term in term_names:
    padj_col = f'padj_{term}'
    if padj_col in anova_df.columns:
        n = (anova_df[padj_col] < padj_thresh).sum()
        clean = term.replace(':', '×')
        print(f"  {clean}: {n}")

print(f"\nOutput files:")
print(f"  - 40_variance_decomposition.png")
print(f"  - 41_interaction_gene_sets.png")
print(f"  - 42_interaction_profiles.png")
print(f"  - multifactor_anova_results.csv")
print(f"  - interaction_gene_categories.csv")
