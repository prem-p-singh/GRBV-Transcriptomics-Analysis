#!/usr/bin/env python3
"""
Comprehensive RNA-Seq analysis of Grapevine Red Blotch Virus (GRBV) infection.
Addresses 6 research questions about infection dynamics, maturity, season, and environment.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# LOAD AND PREPARE DATA
# ============================================================
print("=" * 70)
print("LOADING AND PREPARING DATA")
print("=" * 70)

counts = pd.read_csv('counts_matrix.txt', sep='\t', index_col=0)
meta = pd.read_csv('sample_metadata.txt', sep='\t')

# Ensure column order matches
common_samples = list(set(counts.columns) & set(meta['Sample']))
counts = counts[common_samples]
meta = meta[meta['Sample'].isin(common_samples)].set_index('Sample').loc[common_samples]

print(f"Counts matrix: {counts.shape[0]} genes x {counts.shape[1]} samples")
print(f"Metadata: {meta.shape[0]} samples x {meta.shape[1]} columns")

# Basic filtering: keep genes with at least 10 counts in at least 10% of samples
min_samples = int(0.1 * counts.shape[1])
gene_mask = (counts >= 10).sum(axis=1) >= min_samples
counts_filtered = counts.loc[gene_mask]
print(f"Genes after filtering (>=10 counts in >=10% samples): {counts_filtered.shape[0]}")

# CPM normalization
lib_sizes = counts_filtered.sum(axis=0)
cpm = counts_filtered.div(lib_sizes, axis=1) * 1e6

# Log2 CPM (with pseudocount)
log_cpm = np.log2(cpm + 1)

# Add metadata columns for convenience
meta['Infection_Duration'] = meta['Years_Since_Infection'].apply(
    lambda x: f"{x}yr" if x > 0 else "0yr"
)

# Calculate infection duration in years relative to sampling season
meta['Infection_Age'] = meta.apply(
    lambda row: int(row['Sampling_Season']) - int(row['Years_Since_Infection'])
    if row['Years_Since_Infection'] > 0 else 0, axis=1
)

print(f"\nInfection groups: {meta['Infection_Group'].value_counts().to_dict()}")
print(f"Seasons: {meta['Sampling_Season'].value_counts().to_dict()}")
print(f"Maturity stages: {meta['Maturity_stage'].value_counts().to_dict()}")
print(f"Years since infection: {meta['Years_Since_Infection'].value_counts().sort_index().to_dict()}")

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def run_de_analysis(group1_samples, group2_samples, log_cpm_data):
    """Simple differential expression using t-tests with BH correction."""
    results = []
    g1 = log_cpm_data[group1_samples]
    g2 = log_cpm_data[group2_samples]

    for gene in log_cpm_data.index:
        vals1 = g1.loc[gene].values
        vals2 = g2.loc[gene].values

        # Skip low-variance genes
        if np.std(vals1) == 0 and np.std(vals2) == 0:
            continue

        t_stat, p_val = stats.ttest_ind(vals1, vals2, equal_var=False)
        log2fc = vals1.mean() - vals2.mean()  # already log2
        results.append({
            'gene': gene,
            'log2FC': log2fc,
            'pvalue': p_val,
            'mean_group1': vals1.mean(),
            'mean_group2': vals2.mean()
        })

    df = pd.DataFrame(results)
    if len(df) == 0:
        return df

    # BH correction
    from scipy.stats import rankdata
    n = len(df)
    ranked = rankdata(df['pvalue'])
    df['padj'] = np.minimum(df['pvalue'] * n / ranked, 1.0)
    df = df.sort_values('pvalue')
    return df


def count_deg(de_df, fc_thresh=1.0, padj_thresh=0.05):
    """Count up/down DEGs."""
    sig = de_df[(de_df['padj'] < padj_thresh) & (de_df['log2FC'].abs() >= fc_thresh)]
    up = (sig['log2FC'] > 0).sum()
    down = (sig['log2FC'] < 0).sum()
    return len(sig), up, down


def get_samples(meta_df, **kwargs):
    """Get sample names matching criteria."""
    mask = pd.Series(True, index=meta_df.index)
    for col, val in kwargs.items():
        if isinstance(val, list):
            mask &= meta_df[col].isin(val)
        else:
            mask &= meta_df[col] == val
    return list(meta_df[mask].index)


# Set up figure style
plt.rcParams.update({
    'figure.dpi': 150,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.facecolor': 'white'
})

# ============================================================
# GLOBAL PCA FOR OVERVIEW
# ============================================================
print("\n" + "=" * 70)
print("GLOBAL PCA OVERVIEW")
print("=" * 70)

# PCA on top 2000 most variable genes
var_genes = log_cpm.var(axis=1).nlargest(2000).index
pca_data = log_cpm.loc[var_genes].T
scaler = StandardScaler()
pca_scaled = scaler.fit_transform(pca_data)
pca = PCA(n_components=5)
pca_coords = pca.fit_transform(pca_scaled)
pca_df = pd.DataFrame(pca_coords, columns=[f'PC{i+1}' for i in range(5)], index=pca_data.index)
pca_df = pca_df.join(meta)

print(f"Variance explained: PC1={pca.explained_variance_ratio_[0]:.1%}, "
      f"PC2={pca.explained_variance_ratio_[1]:.1%}, "
      f"PC3={pca.explained_variance_ratio_[2]:.1%}")

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# PCA by infection status
for group, color in [('RB(+)', 'red'), ('RB(-)', 'blue')]:
    sub = pca_df[pca_df['Infection_Group'] == group]
    axes[0,0].scatter(sub['PC1'], sub['PC2'], c=color, alpha=0.5, s=30, label=group, edgecolors='none')
axes[0,0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[0,0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[0,0].set_title("PCA: Infection Status")
axes[0,0].legend()

# PCA by maturity stage
stage_colors = {'Pre_veraison': '#2ca02c', 'Veraison': '#ff7f0e', 'Post_veraison': '#d62728', 'Harvest': '#9467bd'}
stage_order = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
for stage in stage_order:
    sub = pca_df[pca_df['Maturity_stage'] == stage]
    axes[0,1].scatter(sub['PC1'], sub['PC2'], c=stage_colors[stage], alpha=0.5, s=30, label=stage, edgecolors='none')
axes[0,1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[0,1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[0,1].set_title("PCA: Maturity Stage")
axes[0,1].legend()

# PCA by season
for season, color in [(2021, '#1f77b4'), (2022, '#ff7f0e')]:
    sub = pca_df[pca_df['Sampling_Season'] == season]
    axes[1,0].scatter(sub['PC1'], sub['PC2'], c=color, alpha=0.5, s=30, label=str(season), edgecolors='none')
axes[1,0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[1,0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[1,0].set_title("PCA: Sampling Season")
axes[1,0].legend()

# PCA by years since infection
ysi_vals = sorted(meta['Years_Since_Infection'].unique())
cmap = plt.cm.viridis
for i, ysi in enumerate(ysi_vals):
    sub = pca_df[pca_df['Years_Since_Infection'] == ysi]
    axes[1,1].scatter(sub['PC1'], sub['PC2'], c=[cmap(i/len(ysi_vals))], alpha=0.5, s=30,
                       label=f"YSI={ysi}", edgecolors='none')
axes[1,1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[1,1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[1,1].set_title("PCA: Years Since Infection")
axes[1,1].legend(fontsize=8)

plt.suptitle("Global PCA Overview — Top 2000 Variable Genes", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('01_global_pca_overview.png', bbox_inches='tight')
plt.close()
print("Saved: 01_global_pca_overview.png")


# ============================================================
# Q1: HOW DOES GRBV INFECTION CHANGE WITH YEARS OF INFECTION?
# ============================================================
print("\n" + "=" * 70)
print("Q1: GRBV INFECTION CHANGES WITH YEARS OF INFECTION")
print("=" * 70)

# Compare RB(+) at each infection year vs RB(-) controls
ysi_values = sorted([y for y in meta['Years_Since_Infection'].unique() if y > 0])
control_samples = get_samples(meta, Infection_Group='RB(-)')

print(f"\nControl samples (RB-): {len(control_samples)}")
print(f"Infection year cohorts: {ysi_values}")

q1_results = {}
for ysi in ysi_values:
    infected = get_samples(meta, Infection_Group='RB(+)', Years_Since_Infection=ysi)
    de = run_de_analysis(infected, control_samples, log_cpm)
    total, up, down = count_deg(de)
    q1_results[ysi] = {'de_df': de, 'total': total, 'up': up, 'down': down, 'n_samples': len(infected)}
    print(f"  YSI={ysi}: {len(infected)} samples, DEGs={total} (up={up}, down={down})")

# Pairwise: compare consecutive infection years
print("\nPairwise comparisons between consecutive infection years (RB+ only):")
for i in range(len(ysi_values) - 1):
    y1, y2 = ysi_values[i], ysi_values[i+1]
    s1 = get_samples(meta, Infection_Group='RB(+)', Years_Since_Infection=y1)
    s2 = get_samples(meta, Infection_Group='RB(+)', Years_Since_Infection=y2)
    de = run_de_analysis(s1, s2, log_cpm)
    total, up, down = count_deg(de)
    print(f"  YSI={y1} vs YSI={y2}: DEGs={total} (up in {y1}={up}, up in {y2}={down})")

# Calculate infection age (years elapsed since infection at time of sampling)
meta['Actual_Infection_Age'] = meta.apply(
    lambda row: int(row['Sampling_Season']) - int(row['Years_Since_Infection'])
    if row['Years_Since_Infection'] > 0 else 0, axis=1
)

# Trend analysis: correlation of gene expression with infection age
infection_samples = meta[meta['Infection_Group'] == 'RB(+)']
infection_ages = infection_samples['Actual_Infection_Age']

print(f"\nInfection age distribution:")
print(infection_samples['Actual_Infection_Age'].value_counts().sort_index())

# For top variable genes among infected, compute correlation with infection age
infected_log_cpm = log_cpm[infection_samples.index]
var_genes_infected = infected_log_cpm.var(axis=1).nlargest(5000).index

age_corrs = []
for gene in var_genes_infected:
    r, p = stats.spearmanr(infection_ages, infected_log_cpm.loc[gene, infection_samples.index])
    age_corrs.append({'gene': gene, 'spearman_r': r, 'pvalue': p})

age_corr_df = pd.DataFrame(age_corrs)
age_corr_df['padj'] = np.minimum(age_corr_df['pvalue'] * len(age_corr_df) /
                                   stats.rankdata(age_corr_df['pvalue']), 1.0)
age_corr_df = age_corr_df.sort_values('pvalue')

sig_age_genes = age_corr_df[age_corr_df['padj'] < 0.05]
pos_corr = (sig_age_genes['spearman_r'] > 0).sum()
neg_corr = (sig_age_genes['spearman_r'] < 0).sum()
print(f"\nGenes significantly correlated with infection age (padj<0.05): {len(sig_age_genes)}")
print(f"  Positively correlated (increase with age): {pos_corr}")
print(f"  Negatively correlated (decrease with age): {neg_corr}")

# Plot Q1 results
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# DEG counts by infection year
ysi_labels = [str(y) for y in ysi_values]
ups = [q1_results[y]['up'] for y in ysi_values]
downs = [-q1_results[y]['down'] for y in ysi_values]
x = np.arange(len(ysi_values))
axes[0,0].bar(x, ups, color='firebrick', alpha=0.8, label='Upregulated')
axes[0,0].bar(x, downs, color='steelblue', alpha=0.8, label='Downregulated')
axes[0,0].set_xticks(x)
axes[0,0].set_xticklabels(ysi_labels)
axes[0,0].set_xlabel("Year of Infection")
axes[0,0].set_ylabel("Number of DEGs")
axes[0,0].set_title("DEGs: RB(+) at each infection year vs RB(-)\n(|log2FC|≥1, padj<0.05)")
axes[0,0].legend()
axes[0,0].axhline(y=0, color='black', linewidth=0.5)

# Volcano for oldest infection (2016)
de_oldest = q1_results[2016]['de_df']
de_oldest['neg_log10p'] = -np.log10(de_oldest['padj'].clip(lower=1e-300))
colors = np.where(
    (de_oldest['padj'] < 0.05) & (de_oldest['log2FC'] >= 1), 'firebrick',
    np.where((de_oldest['padj'] < 0.05) & (de_oldest['log2FC'] <= -1), 'steelblue', 'grey')
)
axes[0,1].scatter(de_oldest['log2FC'], de_oldest['neg_log10p'], c=colors, s=3, alpha=0.5)
axes[0,1].axvline(1, color='red', linestyle='--', linewidth=0.5)
axes[0,1].axvline(-1, color='blue', linestyle='--', linewidth=0.5)
axes[0,1].axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)
axes[0,1].set_xlabel("log2 Fold Change")
axes[0,1].set_ylabel("-log10(padj)")
axes[0,1].set_title("Volcano: Oldest Infection (2016) vs Controls")

# PCA colored by infection age (RB+ only)
infected_pca = pca_df[pca_df['Infection_Group'] == 'RB(+)'].copy()
infected_pca['Actual_Infection_Age'] = meta.loc[infected_pca.index, 'Actual_Infection_Age']
sc = axes[1,0].scatter(infected_pca['PC1'], infected_pca['PC2'],
                         c=infected_pca['Actual_Infection_Age'], cmap='coolwarm',
                         s=30, alpha=0.6, edgecolors='none')
plt.colorbar(sc, ax=axes[1,0], label='Infection Age (years)')
axes[1,0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[1,0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[1,0].set_title("PCA of RB(+) samples colored by infection age")

# Distribution of Spearman correlations with infection age
axes[1,1].hist(age_corr_df['spearman_r'], bins=50, color='grey', alpha=0.7, edgecolor='black', linewidth=0.3)
sig_mask = age_corr_df['padj'] < 0.05
axes[1,1].hist(age_corr_df.loc[sig_mask, 'spearman_r'], bins=50, color='darkorange', alpha=0.7,
               edgecolor='black', linewidth=0.3, label=f'Significant (n={sig_mask.sum()})')
axes[1,1].set_xlabel("Spearman correlation with infection age")
axes[1,1].set_ylabel("Number of genes")
axes[1,1].set_title("Gene expression correlation with infection age")
axes[1,1].legend()

plt.suptitle("Q1: How does GRBV infection change with years of infection?", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('02_Q1_infection_years.png', bbox_inches='tight')
plt.close()
print("Saved: 02_Q1_infection_years.png")


# ============================================================
# Q2: BERRY MATURITY STAGES AND GRBV INFECTION/HOST RESPONSE
# ============================================================
print("\n" + "=" * 70)
print("Q2: BERRY MATURITY STAGES AND GRBV INFECTION")
print("=" * 70)

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']

# DE: RB(+) vs RB(-) at each maturity stage
q2_results = {}
for stage in stages:
    pos = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=stage)
    neg = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=stage)
    de = run_de_analysis(pos, neg, log_cpm)
    total, up, down = count_deg(de)
    q2_results[stage] = {'de_df': de, 'total': total, 'up': up, 'down': down,
                          'n_pos': len(pos), 'n_neg': len(neg)}
    print(f"  {stage}: RB+={len(pos)}, RB-={len(neg)}, DEGs={total} (up={up}, down={down})")

# DEG overlap between stages
q2_sig_genes = {}
for stage in stages:
    sig = q2_results[stage]['de_df']
    sig = sig[(sig['padj'] < 0.05) & (sig['log2FC'].abs() >= 1)]
    q2_sig_genes[stage] = set(sig['gene'])

print("\nDEG overlap between stages (RB+ vs RB-):")
for i, s1 in enumerate(stages):
    for s2 in stages[i+1:]:
        overlap = len(q2_sig_genes[s1] & q2_sig_genes[s2])
        print(f"  {s1} ∩ {s2}: {overlap} shared DEGs")

# Stage-specific DEGs
for stage in stages:
    others = set()
    for s2 in stages:
        if s2 != stage:
            others |= q2_sig_genes[s2]
    unique = q2_sig_genes[stage] - others
    print(f"  {stage}-specific DEGs: {len(unique)}")

# Consecutive stage transitions within RB(+)
print("\nStage transitions within RB(+):")
for i in range(len(stages) - 1):
    s1, s2 = stages[i], stages[i+1]
    samp1 = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=s1)
    samp2 = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=s2)
    de = run_de_analysis(samp2, samp1, log_cpm)
    total, up, down = count_deg(de)
    print(f"  {s1} -> {s2}: DEGs={total} (up={up}, down={down})")

# Same for RB(-)
print("\nStage transitions within RB(-):")
for i in range(len(stages) - 1):
    s1, s2 = stages[i], stages[i+1]
    samp1 = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=s1)
    samp2 = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=s2)
    de = run_de_analysis(samp2, samp1, log_cpm)
    total, up, down = count_deg(de)
    print(f"  {s1} -> {s2}: DEGs={total} (up={up}, down={down})")

# Plot Q2
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# DEG bar chart per stage
x = np.arange(len(stages))
ups = [q2_results[s]['up'] for s in stages]
downs = [-q2_results[s]['down'] for s in stages]
axes[0,0].bar(x, ups, color='firebrick', alpha=0.8, label='Upregulated in RB(+)')
axes[0,0].bar(x, downs, color='steelblue', alpha=0.8, label='Downregulated in RB(+)')
axes[0,0].set_xticks(x)
axes[0,0].set_xticklabels(stages, rotation=15)
axes[0,0].set_xlabel("Maturity Stage")
axes[0,0].set_ylabel("Number of DEGs")
axes[0,0].set_title("DEGs per maturity stage: RB(+) vs RB(-)")
axes[0,0].legend()
axes[0,0].axhline(0, color='black', linewidth=0.5)

# Heatmap: top DEGs across stages (using fold changes)
all_stage_genes = set()
for s in stages:
    de = q2_results[s]['de_df']
    top = de[(de['padj'] < 0.05)].nsmallest(50, 'padj')['gene'].tolist()
    all_stage_genes.update(top)

# Build FC matrix
fc_matrix = pd.DataFrame(index=list(all_stage_genes), columns=stages)
for s in stages:
    de = q2_results[s]['de_df'].set_index('gene')
    for g in all_stage_genes:
        if g in de.index:
            fc_matrix.loc[g, s] = de.loc[g, 'log2FC']
        else:
            fc_matrix.loc[g, s] = 0

fc_matrix = fc_matrix.astype(float)
# Cluster and show top 40
fc_var = fc_matrix.var(axis=1).nlargest(40).index
sns.heatmap(fc_matrix.loc[fc_var], ax=axes[0,1], cmap='RdBu_r', center=0,
            xticklabels=True, yticklabels=False, cbar_kws={'label': 'log2FC'})
axes[0,1].set_title("Top stage-variable DEGs (log2FC RB+ vs RB-)")

# PCA by stage, split by infection
for grp, marker in [('RB(+)', 'o'), ('RB(-)', 's')]:
    for stage in stages:
        sub = pca_df[(pca_df['Infection_Group'] == grp) & (pca_df['Maturity_stage'] == stage)]
        axes[1,0].scatter(sub['PC1'], sub['PC2'], c=stage_colors[stage], marker=marker,
                          s=20, alpha=0.4, edgecolors='none')

# Custom legend
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=stage_colors[s],
                           markersize=8, label=f'RB(+) {s}') for s in stages]
legend_elements += [Line2D([0], [0], marker='s', color='w', markerfacecolor=stage_colors[s],
                            markersize=8, label=f'RB(-) {s}') for s in stages]
axes[1,0].legend(handles=legend_elements, fontsize=7, ncol=2)
axes[1,0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[1,0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[1,0].set_title("PCA: Stage x Infection interaction")

# Venn-like overlap summary
overlap_data = np.zeros((len(stages), len(stages)))
for i, s1 in enumerate(stages):
    for j, s2 in enumerate(stages):
        if i == j:
            overlap_data[i,j] = len(q2_sig_genes[s1])
        else:
            overlap_data[i,j] = len(q2_sig_genes[s1] & q2_sig_genes[s2])

sns.heatmap(pd.DataFrame(overlap_data, index=stages, columns=stages), ax=axes[1,1],
            annot=True, fmt='.0f', cmap='YlOrRd', cbar_kws={'label': 'Gene count'})
axes[1,1].set_title("DEG overlap matrix between stages")

plt.suptitle("Q2: Berry maturity stages influence on GRBV infection", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('03_Q2_maturity_stages.png', bbox_inches='tight')
plt.close()
print("Saved: 03_Q2_maturity_stages.png")


# ============================================================
# Q3: SEASONAL EFFECTS ON GRBV INFECTION AND RESPONSE
# ============================================================
print("\n" + "=" * 70)
print("Q3: SEASONAL EFFECTS ON GRBV INFECTION")
print("=" * 70)

seasons = [2021, 2022]

# DE: RB(+) vs RB(-) per season
for season in seasons:
    pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season)
    neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season)
    de = run_de_analysis(pos, neg, log_cpm)
    total, up, down = count_deg(de)
    print(f"  Season {season}: RB+={len(pos)}, RB-={len(neg)}, DEGs={total} (up={up}, down={down})")

# DE: 2021 vs 2022 within RB(+) and RB(-)
print("\nSeasonal differences within each group:")
for grp in ['RB(+)', 'RB(-)']:
    s2021 = get_samples(meta, Infection_Group=grp, Sampling_Season=2021)
    s2022 = get_samples(meta, Infection_Group=grp, Sampling_Season=2022)
    de = run_de_analysis(s2022, s2021, log_cpm)
    total, up, down = count_deg(de)
    print(f"  {grp}: 2022 vs 2021, DEGs={total} (up in 2022={up}, down in 2022={down})")

# Season x Stage interaction for RB(+) vs RB(-)
print("\nSeason x Stage DEGs (RB+ vs RB-):")
q3_detailed = {}
for season in seasons:
    for stage in stages:
        pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season, Maturity_stage=stage)
        neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season, Maturity_stage=stage)
        if len(pos) < 3 or len(neg) < 3:
            print(f"  {season}-{stage}: Skipped (too few samples: RB+={len(pos)}, RB-={len(neg)})")
            continue
        de = run_de_analysis(pos, neg, log_cpm)
        total, up, down = count_deg(de)
        q3_detailed[(season, stage)] = {'total': total, 'up': up, 'down': down}
        print(f"  {season}-{stage}: DEGs={total} (up={up}, down={down})")

# Plot Q3
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# Season x Stage heatmap of DEG counts
hm_data = pd.DataFrame(index=stages, columns=[str(s) for s in seasons])
for (season, stage), v in q3_detailed.items():
    hm_data.loc[stage, str(season)] = v['total']
hm_data = hm_data.astype(float)
sns.heatmap(hm_data, ax=axes[0,0], annot=True, fmt='.0f', cmap='YlOrRd',
            cbar_kws={'label': 'Total DEGs'})
axes[0,0].set_title("Total DEGs (RB+ vs RB-) by Season x Stage")

# Up/Down split per season x stage
categories = []
up_vals = []
down_vals = []
for season in seasons:
    for stage in stages:
        key = (season, stage)
        if key in q3_detailed:
            categories.append(f"{season}\n{stage}")
            up_vals.append(q3_detailed[key]['up'])
            down_vals.append(-q3_detailed[key]['down'])

x = np.arange(len(categories))
axes[0,1].bar(x, up_vals, color='firebrick', alpha=0.8, width=0.6, label='Up')
axes[0,1].bar(x, down_vals, color='steelblue', alpha=0.8, width=0.6, label='Down')
axes[0,1].set_xticks(x)
axes[0,1].set_xticklabels(categories, fontsize=7, rotation=45, ha='right')
axes[0,1].set_ylabel("Number of DEGs")
axes[0,1].set_title("DEGs by Season x Stage")
axes[0,1].legend()
axes[0,1].axhline(0, color='black', linewidth=0.5)

# PCA split by season
for season, marker in [(2021, 'o'), (2022, 's')]:
    for grp, color in [('RB(+)', 'red'), ('RB(-)', 'blue')]:
        sub = pca_df[(pca_df['Sampling_Season'] == season) & (pca_df['Infection_Group'] == grp)]
        axes[1,0].scatter(sub['PC1'], sub['PC2'], c=color, marker=marker, s=20, alpha=0.4, edgecolors='none')

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='RB(+) 2021'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='red', markersize=8, label='RB(+) 2022'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='RB(-) 2021'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='blue', markersize=8, label='RB(-) 2022'),
]
axes[1,0].legend(handles=legend_elements, fontsize=8)
axes[1,0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[1,0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[1,0].set_title("PCA: Season x Infection")

# Seasonal consistency: scatter of log2FC in 2021 vs 2022
# Use overall RB+ vs RB- DE at each season
pos_21 = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=2021)
neg_21 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2021)
pos_22 = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=2022)
neg_22 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2022)

de_21 = run_de_analysis(pos_21, neg_21, log_cpm).set_index('gene')
de_22 = run_de_analysis(pos_22, neg_22, log_cpm).set_index('gene')

common_genes = list(set(de_21.index) & set(de_22.index))
fc_21 = de_21.loc[common_genes, 'log2FC']
fc_22 = de_22.loc[common_genes, 'log2FC']

r_season, p_season = stats.pearsonr(fc_21, fc_22)
axes[1,1].scatter(fc_21, fc_22, s=2, alpha=0.3, c='grey')
axes[1,1].plot([-5, 5], [-5, 5], 'r--', linewidth=0.5)
axes[1,1].set_xlabel("log2FC (RB+ vs RB-) in 2021")
axes[1,1].set_ylabel("log2FC (RB+ vs RB-) in 2022")
axes[1,1].set_title(f"Seasonal consistency of GRBV response\nPearson r={r_season:.3f}")
axes[1,1].set_xlim(-5, 5)
axes[1,1].set_ylim(-5, 5)

plt.suptitle("Q3: Seasonal effects on GRBV infection and response", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('04_Q3_seasonal_effects.png', bbox_inches='tight')
plt.close()
print(f"\nSeasonal consistency: Pearson r = {r_season:.4f} (p = {p_season:.2e})")
print("Saved: 04_Q3_seasonal_effects.png")


# ============================================================
# Q4: ENVIRONMENTAL EFFECTS ON GRBV INFECTION AND RESPONSE
# ============================================================
print("\n" + "=" * 70)
print("Q4: ENVIRONMENTAL EFFECTS ON GRBV INFECTION")
print("=" * 70)

# Season serves as a proxy for environmental variation (climate year-to-year)
# We quantify how much of the variance is explained by season vs infection vs stage

# PERMANOVA-like analysis: compute R² for each factor
from sklearn.metrics import pairwise_distances

# Use top 2000 variable genes, scaled
dist_matrix = pairwise_distances(pca_scaled, metric='euclidean')

def pseudo_r2(dist_matrix, grouping):
    """Calculate pseudo R² (ratio of between-group to total SS from distance matrix)."""
    groups = grouping.values
    unique_groups = np.unique(groups)
    n = len(groups)

    # Total SS
    total_ss = np.sum(dist_matrix ** 2) / (2 * n)

    # Within SS
    within_ss = 0
    for g in unique_groups:
        idx = np.where(groups == g)[0]
        if len(idx) > 1:
            sub_dist = dist_matrix[np.ix_(idx, idx)]
            within_ss += np.sum(sub_dist ** 2) / (2 * len(idx))

    between_ss = total_ss - within_ss
    r2 = between_ss / total_ss
    return r2

# Align metadata with PCA sample order
meta_aligned = meta.loc[pca_data.index]

r2_infection = pseudo_r2(dist_matrix, meta_aligned['Infection_Group'])
r2_stage = pseudo_r2(dist_matrix, meta_aligned['Maturity_stage'])
r2_season = pseudo_r2(dist_matrix, meta_aligned['Sampling_Season'].astype(str))
r2_ysi = pseudo_r2(dist_matrix, meta_aligned['Years_Since_Infection'].astype(str))

print(f"Variance explained (pseudo-R²) by each factor:")
print(f"  Maturity stage:       {r2_stage:.4f} ({r2_stage*100:.1f}%)")
print(f"  Sampling season:      {r2_season:.4f} ({r2_season*100:.1f}%)")
print(f"  Infection status:     {r2_infection:.4f} ({r2_infection*100:.1f}%)")
print(f"  Year of infection:    {r2_ysi:.4f} ({r2_ysi*100:.1f}%)")

# Environmental effect: how much does season change the GRBV response?
# Compare the "infection effect" (FC) at each stage between seasons
print("\nEnvironmental modulation of GRBV response by stage:")
env_modulation = {}
for stage in stages:
    pos_21 = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=2021, Maturity_stage=stage)
    neg_21 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2021, Maturity_stage=stage)
    pos_22 = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=2022, Maturity_stage=stage)
    neg_22 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2022, Maturity_stage=stage)

    if len(pos_21) < 3 or len(neg_21) < 3 or len(pos_22) < 3 or len(neg_22) < 3:
        print(f"  {stage}: Skipped")
        continue

    de21 = run_de_analysis(pos_21, neg_21, log_cpm).set_index('gene')
    de22 = run_de_analysis(pos_22, neg_22, log_cpm).set_index('gene')

    common = list(set(de21.index) & set(de22.index))
    r, p = stats.pearsonr(de21.loc[common, 'log2FC'], de22.loc[common, 'log2FC'])
    env_modulation[stage] = {'r': r, 'p': p, 'n_genes': len(common)}
    print(f"  {stage}: FC correlation between seasons r={r:.3f} (p={p:.2e})")

# Interaction: genes where season significantly changes the infection response
print("\nGenes with season-dependent infection response (interaction effects):")
for stage in stages:
    pos_21 = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=2021, Maturity_stage=stage)
    neg_21 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2021, Maturity_stage=stage)
    pos_22 = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=2022, Maturity_stage=stage)
    neg_22 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2022, Maturity_stage=stage)

    if len(pos_21) < 3 or len(neg_21) < 3 or len(pos_22) < 3 or len(neg_22) < 3:
        continue

    # Infection effect in each season
    effect_21 = log_cpm[pos_21].mean(axis=1) - log_cpm[neg_21].mean(axis=1)
    effect_22 = log_cpm[pos_22].mean(axis=1) - log_cpm[neg_22].mean(axis=1)

    # Difference of differences
    interaction = (effect_22 - effect_21).abs()
    n_interaction = (interaction > 1).sum()  # genes with >1 log2FC difference between seasons
    print(f"  {stage}: {n_interaction} genes with |ΔΔlog2FC| > 1 between seasons")

# Plot Q4
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# Variance partitioning bar chart
factors = ['Maturity\nStage', 'Sampling\nSeason', 'Infection\nStatus', 'Year of\nInfection']
r2_vals = [r2_stage, r2_season, r2_infection, r2_ysi]
bars = axes[0,0].bar(factors, [v*100 for v in r2_vals], color=['#2ca02c', '#ff7f0e', '#d62728', '#9467bd'], alpha=0.8)
axes[0,0].set_ylabel("Pseudo R² (%)")
axes[0,0].set_title("Variance explained by each factor\n(distance-based)")
for bar, val in zip(bars, r2_vals):
    axes[0,0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    f'{val*100:.1f}%', ha='center', fontsize=10)

# Environmental modulation: FC correlation per stage
if env_modulation:
    stage_labels = list(env_modulation.keys())
    r_vals = [env_modulation[s]['r'] for s in stage_labels]
    axes[0,1].bar(stage_labels, r_vals, color=[stage_colors[s] for s in stage_labels], alpha=0.8)
    axes[0,1].set_ylabel("Pearson r (FC 2021 vs FC 2022)")
    axes[0,1].set_title("Consistency of infection response\nacross seasons by stage")
    axes[0,1].set_ylim(0, 1)
    axes[0,1].axhline(0.8, color='grey', linestyle='--', linewidth=0.5, label='r=0.8')
    axes[0,1].legend()

# PC1/PC2 mean per group (season x infection)
group_means = pca_df.groupby(['Sampling_Season', 'Infection_Group'])[['PC1', 'PC2']].mean()
for (season, grp), row in group_means.iterrows():
    color = 'red' if grp == 'RB(+)' else 'blue'
    marker = 'o' if season == 2021 else 's'
    axes[1,0].scatter(row['PC1'], row['PC2'], c=color, marker=marker, s=200, edgecolors='black',
                       linewidth=2, zorder=5)
    axes[1,0].annotate(f'{grp}\n{season}', (row['PC1'], row['PC2']), fontsize=8,
                        ha='center', va='bottom', xytext=(0, 10), textcoords='offset points')

# Also show all points faintly
for grp, color in [('RB(+)', 'red'), ('RB(-)', 'blue')]:
    sub = pca_df[pca_df['Infection_Group'] == grp]
    axes[1,0].scatter(sub['PC1'], sub['PC2'], c=color, s=10, alpha=0.15, edgecolors='none')

axes[1,0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[1,0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[1,0].set_title("Group centroids: Season x Infection in PCA space")

# Box plots of library sizes by season and infection (proxy for overall transcriptional activity)
plot_data = pd.DataFrame({'Library_Size': lib_sizes, 'Season': meta['Sampling_Season'].astype(str),
                           'Infection': meta['Infection_Group']})
order = ['RB(-)', 'RB(+)']
sns.boxplot(data=plot_data, x='Season', y='Library_Size', hue='Infection', ax=axes[1,1],
            hue_order=order, palette={'RB(+)': 'red', 'RB(-)': 'blue'})
axes[1,1].set_title("Library sizes by Season x Infection\n(proxy for transcriptional activity)")
axes[1,1].set_ylabel("Total read counts")

plt.suptitle("Q4: Environmental effects on GRBV infection and response", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('05_Q4_environmental_effects.png', bbox_inches='tight')
plt.close()
print("Saved: 05_Q4_environmental_effects.png")


# ============================================================
# Q5: ENVIRONMENTAL CHANGES BETWEEN UNINFECTED GRAPEVINES
# ============================================================
print("\n" + "=" * 70)
print("Q5: ENVIRONMENTAL CHANGES IN UNINFECTED GRAPEVINES")
print("=" * 70)

# Focus on RB(-) only: seasonal differences
neg_2021 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2021)
neg_2022 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2022)

de_neg_season = run_de_analysis(neg_2022, neg_2021, log_cpm)
total_neg, up_neg, down_neg = count_deg(de_neg_season)
print(f"RB(-) 2022 vs 2021: DEGs={total_neg} (up in 2022={up_neg}, down in 2022={down_neg})")

# Stage-by-stage seasonal differences in RB(-)
print("\nSeasonal differences in RB(-) by stage:")
q5_stage_season = {}
for stage in stages:
    s21 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2021, Maturity_stage=stage)
    s22 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2022, Maturity_stage=stage)
    if len(s21) < 3 or len(s22) < 3:
        print(f"  {stage}: Skipped")
        continue
    de = run_de_analysis(s22, s21, log_cpm)
    total, up, down = count_deg(de)
    q5_stage_season[stage] = {'de_df': de, 'total': total, 'up': up, 'down': down}
    print(f"  {stage}: DEGs={total} (up in 2022={up}, down in 2022={down})")

# Compare environmental DEGs in RB(-) with infection DEGs
# How many seasonal DEGs (in controls) overlap with GRBV-response DEGs?
print("\nOverlap: Environmental DEGs in RB(-) vs GRBV-response DEGs:")
env_degs_neg = set(de_neg_season[(de_neg_season['padj'] < 0.05) & (de_neg_season['log2FC'].abs() >= 1)]['gene'])

# Overall GRBV response DEGs
all_pos = get_samples(meta, Infection_Group='RB(+)')
all_neg = get_samples(meta, Infection_Group='RB(-)')
de_overall = run_de_analysis(all_pos, all_neg, log_cpm)
grbv_degs = set(de_overall[(de_overall['padj'] < 0.05) & (de_overall['log2FC'].abs() >= 1)]['gene'])

overlap = env_degs_neg & grbv_degs
print(f"  Environmental DEGs in RB(-): {len(env_degs_neg)}")
print(f"  GRBV response DEGs: {len(grbv_degs)}")
print(f"  Overlap: {len(overlap)} ({len(overlap)/max(len(env_degs_neg),1)*100:.1f}% of environmental DEGs)")

# Direction concordance in overlapping genes
if len(overlap) > 0:
    de_neg_idx = de_neg_season.set_index('gene')
    de_overall_idx = de_overall.set_index('gene')
    concordant = 0
    discordant = 0
    for g in overlap:
        if g in de_neg_idx.index and g in de_overall_idx.index:
            if np.sign(de_neg_idx.loc[g, 'log2FC']) == np.sign(de_overall_idx.loc[g, 'log2FC']):
                concordant += 1
            else:
                discordant += 1
    print(f"  Same direction: {concordant}, Opposite direction: {discordant}")

# Developmental program in RB(-): stage transitions
print("\nDevelopmental transcriptome in healthy grapevines:")
for i in range(len(stages) - 1):
    s1, s2 = stages[i], stages[i+1]
    samp1 = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=s1)
    samp2 = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=s2)
    de = run_de_analysis(samp2, samp1, log_cpm)
    total, up, down = count_deg(de)
    print(f"  {s1} -> {s2}: DEGs={total} (up={up}, down={down})")

# Plot Q5
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# PCA of RB(-) only colored by season
neg_pca = pca_df[pca_df['Infection_Group'] == 'RB(-)'].copy()
for season, color in [(2021, '#1f77b4'), (2022, '#ff7f0e')]:
    sub = neg_pca[neg_pca['Sampling_Season'] == season]
    axes[0,0].scatter(sub['PC1'], sub['PC2'], c=color, s=40, alpha=0.6, label=str(season), edgecolors='none')
axes[0,0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[0,0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[0,0].set_title("PCA of uninfected RB(-) by season")
axes[0,0].legend()

# PCA of RB(-) by stage
for stage in stages:
    sub = neg_pca[neg_pca['Maturity_stage'] == stage]
    axes[0,1].scatter(sub['PC1'], sub['PC2'], c=stage_colors[stage], s=40, alpha=0.6,
                       label=stage, edgecolors='none')
axes[0,1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axes[0,1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axes[0,1].set_title("PCA of uninfected RB(-) by maturity stage")
axes[0,1].legend()

# Overlap Venn (as bar chart)
labels = ['Environmental\nDEGs only', 'Shared', 'GRBV response\nDEGs only']
vals = [len(env_degs_neg) - len(overlap), len(overlap), len(grbv_degs) - len(overlap)]
colors_bar = ['#ff7f0e', '#2ca02c', '#d62728']
axes[1,0].bar(labels, vals, color=colors_bar, alpha=0.8)
axes[1,0].set_ylabel("Number of genes")
axes[1,0].set_title("Environmental vs GRBV response DEGs overlap")
for i, (l, v) in enumerate(zip(labels, vals)):
    axes[1,0].text(i, v + max(vals)*0.02, str(v), ha='center', fontsize=10)

# Seasonal DEGs per stage in RB(-)
if q5_stage_season:
    stage_labels = list(q5_stage_season.keys())
    ups_q5 = [q5_stage_season[s]['up'] for s in stage_labels]
    downs_q5 = [-q5_stage_season[s]['down'] for s in stage_labels]
    x = np.arange(len(stage_labels))
    axes[1,1].bar(x, ups_q5, color='#ff7f0e', alpha=0.8, label='Up in 2022')
    axes[1,1].bar(x, downs_q5, color='#1f77b4', alpha=0.8, label='Down in 2022')
    axes[1,1].set_xticks(x)
    axes[1,1].set_xticklabels(stage_labels, rotation=15)
    axes[1,1].set_ylabel("Number of DEGs")
    axes[1,1].set_title("Seasonal DEGs in RB(-) by maturity stage")
    axes[1,1].legend()
    axes[1,1].axhline(0, color='black', linewidth=0.5)

plt.suptitle("Q5: Environmental changes in uninfected grapevines", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('06_Q5_uninfected_environment.png', bbox_inches='tight')
plt.close()
print("Saved: 06_Q5_uninfected_environment.png")


# ============================================================
# Q6: HOST TRANSCRIPTOME CHANGES OVER COURSE OF GRBV INFECTION
# ============================================================
print("\n" + "=" * 70)
print("Q6: HOST TRANSCRIPTOME CHANGES OVER INFECTION COURSE")
print("=" * 70)

# Infection age groups
meta['Actual_Infection_Age'] = meta.apply(
    lambda row: int(row['Sampling_Season']) - int(row['Years_Since_Infection'])
    if row['Years_Since_Infection'] > 0 else 0, axis=1
)

infected_meta = meta[meta['Infection_Group'] == 'RB(+)']
age_groups = sorted(infected_meta['Actual_Infection_Age'].unique())
print(f"Infection ages (years): {age_groups}")
print(f"Samples per age: {infected_meta['Actual_Infection_Age'].value_counts().sort_index().to_dict()}")

# Progressive DE: each age vs controls
print("\nProgressive infection response (each age vs RB- controls):")
q6_progressive = {}
for age in age_groups:
    infected = list(infected_meta[infected_meta['Actual_Infection_Age'] == age].index)
    de = run_de_analysis(infected, all_neg, log_cpm)
    total, up, down = count_deg(de)
    q6_progressive[age] = {'de_df': de, 'total': total, 'up': up, 'down': down, 'n': len(infected)}
    print(f"  Age {age}yr: n={len(infected)}, DEGs={total} (up={up}, down={down})")

# Trajectory: how transcriptome shifts between consecutive ages
print("\nConsecutive age transitions:")
for i in range(len(age_groups) - 1):
    a1, a2 = age_groups[i], age_groups[i+1]
    s1 = list(infected_meta[infected_meta['Actual_Infection_Age'] == a1].index)
    s2 = list(infected_meta[infected_meta['Actual_Infection_Age'] == a2].index)
    de = run_de_analysis(s2, s1, log_cpm)
    total, up, down = count_deg(de)
    print(f"  {a1}yr -> {a2}yr: DEGs={total} (up={up}, down={down})")

# Earliest vs latest infection comparison
if len(age_groups) >= 2:
    earliest = list(infected_meta[infected_meta['Actual_Infection_Age'] == age_groups[0]].index)
    latest = list(infected_meta[infected_meta['Actual_Infection_Age'] == age_groups[-1]].index)
    de_el = run_de_analysis(latest, earliest, log_cpm)
    total_el, up_el, down_el = count_deg(de_el)
    print(f"\nEarliest ({age_groups[0]}yr) vs Latest ({age_groups[-1]}yr) infection:")
    print(f"  DEGs={total_el} (up in later={up_el}, down in later={down_el})")

# Time-series clustering: K-means on infection age profiles
print("\nTime-series gene clustering (infection age profiles):")
# Build age profile matrix: mean log2FC (vs controls) at each age
age_profiles = pd.DataFrame(index=log_cpm.index)
control_mean = log_cpm[all_neg].mean(axis=1)

for age in age_groups:
    infected = list(infected_meta[infected_meta['Actual_Infection_Age'] == age].index)
    age_profiles[f'age_{age}'] = log_cpm[infected].mean(axis=1) - control_mean

# Filter to genes that are DEG in at least one age
any_deg = set()
for age in age_groups:
    de = q6_progressive[age]['de_df']
    sig = de[(de['padj'] < 0.05) & (de['log2FC'].abs() >= 1)]
    any_deg.update(sig['gene'])

print(f"Genes that are DEG in at least one infection age: {len(any_deg)}")

if len(any_deg) > 10:
    from sklearn.cluster import KMeans

    deg_profiles = age_profiles.loc[list(any_deg)]

    # Standardize per gene
    from sklearn.preprocessing import StandardScaler
    scaler_ts = StandardScaler()
    scaled_profiles = pd.DataFrame(
        scaler_ts.fit_transform(deg_profiles.T).T,
        index=deg_profiles.index, columns=deg_profiles.columns
    )

    # Determine K (use 6 clusters for interpretability)
    n_clusters = min(6, len(any_deg))
    km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    clusters = km.fit_predict(scaled_profiles)

    cluster_counts = pd.Series(clusters).value_counts().sort_index()
    print(f"Cluster sizes: {cluster_counts.to_dict()}")

# PCA on RB(+) only, colored by age with trajectory arrows
infected_pca_data = log_cpm.loc[var_genes, infected_meta.index].T
pca_infected = PCA(n_components=3)
infected_scaled = StandardScaler().fit_transform(infected_pca_data)
infected_coords = pca_infected.fit_transform(infected_scaled)
infected_pca_df = pd.DataFrame(infected_coords, columns=['PC1', 'PC2', 'PC3'], index=infected_meta.index)
infected_pca_df['Age'] = infected_meta['Actual_Infection_Age']
infected_pca_df['Stage'] = infected_meta['Maturity_stage']

# Plot Q6
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# Progressive DEG accumulation
ages_plot = [a for a in age_groups]
totals = [q6_progressive[a]['total'] for a in ages_plot]
ups_q6 = [q6_progressive[a]['up'] for a in ages_plot]
downs_q6 = [q6_progressive[a]['down'] for a in ages_plot]

axes[0,0].plot(ages_plot, totals, 'ko-', linewidth=2, markersize=8, label='Total DEGs')
axes[0,0].fill_between(ages_plot, ups_q6, alpha=0.3, color='firebrick', label='Upregulated')
axes[0,0].fill_between(ages_plot, downs_q6, alpha=0.3, color='steelblue', label='Downregulated')
axes[0,0].set_xlabel("Infection Age (years)")
axes[0,0].set_ylabel("Number of DEGs (vs controls)")
axes[0,0].set_title("Accumulation of DEGs with infection age")
axes[0,0].legend()

# PCA trajectory colored by age
sc = axes[0,1].scatter(infected_pca_df['PC1'], infected_pca_df['PC2'],
                        c=infected_pca_df['Age'], cmap='coolwarm', s=20, alpha=0.5, edgecolors='none')
plt.colorbar(sc, ax=axes[0,1], label='Infection Age (yr)')

# Add centroids with arrows
centroids = infected_pca_df.groupby('Age')[['PC1', 'PC2']].mean()
axes[0,1].plot(centroids['PC1'], centroids['PC2'], 'k-o', markersize=10, linewidth=2, zorder=5)
for age in centroids.index:
    axes[0,1].annotate(f'{age}yr', (centroids.loc[age, 'PC1'], centroids.loc[age, 'PC2']),
                        fontsize=9, fontweight='bold', ha='center', va='bottom',
                        xytext=(0, 8), textcoords='offset points')
axes[0,1].set_xlabel(f"PC1 ({pca_infected.explained_variance_ratio_[0]:.1%})")
axes[0,1].set_ylabel(f"PC2 ({pca_infected.explained_variance_ratio_[1]:.1%})")
axes[0,1].set_title("Transcriptome trajectory over infection course")

# Time-series cluster profiles
if len(any_deg) > 10:
    for c in range(n_clusters):
        cluster_genes = scaled_profiles.index[clusters == c]
        cluster_data = age_profiles.loc[cluster_genes]
        mean_profile = cluster_data.mean()
        std_profile = cluster_data.std()
        x_vals = age_groups
        axes[1,0].plot(x_vals, mean_profile.values, '-o', linewidth=2,
                        label=f'C{c+1} (n={len(cluster_genes)})')
        axes[1,0].fill_between(x_vals, (mean_profile - std_profile).values,
                                (mean_profile + std_profile).values, alpha=0.1)

    axes[1,0].set_xlabel("Infection Age (years)")
    axes[1,0].set_ylabel("Mean log2FC vs controls")
    axes[1,0].set_title("Gene expression clusters over infection course")
    axes[1,0].legend(fontsize=8)
    axes[1,0].axhline(0, color='grey', linestyle='--', linewidth=0.5)

# Heatmap: top progressive genes across ages
# Get genes with strongest monotonic trend
from scipy.stats import spearmanr
trend_genes = []
for gene in any_deg:
    profile = [age_profiles.loc[gene, f'age_{a}'] for a in age_groups]
    r, p = spearmanr(age_groups, profile)
    trend_genes.append({'gene': gene, 'r': r, 'p': p})

trend_df = pd.DataFrame(trend_genes).sort_values('p')
top_increasing = trend_df[trend_df['r'] > 0].head(15)['gene'].tolist()
top_decreasing = trend_df[trend_df['r'] < 0].head(15)['gene'].tolist()
trend_genes_list = top_increasing + top_decreasing

if len(trend_genes_list) > 0:
    hm_data_q6 = age_profiles.loc[trend_genes_list, [f'age_{a}' for a in age_groups]]
    hm_data_q6.columns = [f'{a}yr' for a in age_groups]
    sns.heatmap(hm_data_q6, ax=axes[1,1], cmap='RdBu_r', center=0,
                yticklabels=True if len(trend_genes_list) <= 40 else False,
                xticklabels=True, cbar_kws={'label': 'log2FC vs controls'})
    axes[1,1].set_title("Top monotonically trending genes\n(top: increasing, bottom: decreasing)")
    axes[1,1].set_xlabel("Infection Age")
else:
    axes[1,1].text(0.5, 0.5, "No significant trend genes found", ha='center', va='center')

plt.suptitle("Q6: Host transcriptome changes over GRBV infection course", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig('07_Q6_transcriptome_trajectory.png', bbox_inches='tight')
plt.close()
print("Saved: 07_Q6_transcriptome_trajectory.png")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE — SUMMARY OF KEY FINDINGS")
print("=" * 70)

print("""
Q1: GRBV INFECTION CHANGES WITH YEARS OF INFECTION
- Compared RB(+) samples at each infection year (2016-2022) against RB(-) controls
- Tracked DEG accumulation and identified genes correlated with infection age
- Both positive and negative correlations with infection duration exist

Q2: BERRY MATURITY STAGE INFLUENCE
- DEGs identified at each maturity stage (RB+ vs RB-)
- Stage-specific and shared DEGs quantified
- Developmental transitions compared between infected and healthy vines

Q3: SEASONAL EFFECTS
- Season-specific GRBV response quantified (2021 vs 2022)
- Season x Stage interaction effects on DEG counts
- Cross-season consistency of fold-changes measured

Q4: ENVIRONMENTAL EFFECTS
- Variance partitioning across all factors (stage, season, infection, infection year)
- Environmental modulation of GRBV response quantified per stage
- Interaction effects (genes with season-dependent infection response) identified

Q5: UNINFECTED GRAPEVINE ENVIRONMENTAL CHANGES
- Seasonal DEGs in RB(-) controls quantified
- Overlap between environmental DEGs and GRBV-response DEGs assessed
- Direction concordance of shared genes analyzed

Q6: HOST TRANSCRIPTOME TRAJECTORY
- Progressive DEG accumulation with infection age
- PCA trajectory showing transcriptome drift
- Time-series clustering reveals distinct temporal gene expression patterns
- Monotonically trending genes identified
""")

print("Generated plots:")
print("  01_global_pca_overview.png")
print("  02_Q1_infection_years.png")
print("  03_Q2_maturity_stages.png")
print("  04_Q3_seasonal_effects.png")
print("  05_Q4_environmental_effects.png")
print("  06_Q5_uninfected_environment.png")
print("  07_Q6_transcriptome_trajectory.png")
