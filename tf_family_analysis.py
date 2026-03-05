#!/usr/bin/env python3
"""
Transcription Factor Family Profiling for GRBV Transcriptomics.

Systematically profiles 6 major TF families (WRKY, MYB, NAC, bHLH, ERF/AP2, HSF)
across all season × stage conditions.  Identifies TF families most disrupted by GRBV
and correlates TF expression with downstream target genes.

Outputs:
  - Fig 43: TF family heatmap across conditions
  - Fig 44: TF family DEG enrichment balloon plot
  - Fig 45: Top TF gene expression boxplots
  - tf_family_DE_summary.csv
  - tf_family_gene_details.csv
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
from scipy.stats import rankdata, fisher_exact
from collections import defaultdict
import textwrap
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
print("TF FAMILY ANALYSIS: LOADING DATA")
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

expressed = set(log_cpm.index)
print(f"Genes: {counts_filtered.shape[0]}, Samples: {counts_filtered.shape[1]}")

# ============================================================
# 2. QUERY TF FAMILIES FROM ANNOTATION DATABASE
# ============================================================
print("\n" + "=" * 70)
print("QUERYING TF FAMILIES FROM ANNOTATION DB")
print("=" * 70)

conn = sqlite3.connect('org.Vvinifera.eg.db/inst/extdata/org.Vvinifera.eg.sqlite')
cursor = conn.cursor()

# Define TF family queries based on gene name patterns
tf_queries = {
    'WRKY': "LOWER(gi.GENENAME) LIKE '%wrky%'",
    'MYB': "(LOWER(gi.GENENAME) LIKE '%myb%' AND LOWER(gi.GENENAME) NOT LIKE '%myb-like%dna%')",
    'NAC': "(LOWER(gi.GENENAME) LIKE '%nac domain%' OR LOWER(gi.GENENAME) LIKE '%nac transcription%' OR LOWER(gi.GENENAME) LIKE 'nac%')",
    'bHLH': "LOWER(gi.GENENAME) LIKE '%bhlh%'",
    'ERF/AP2': "(LOWER(gi.GENENAME) LIKE '%ethylene-responsive%' OR LOWER(gi.GENENAME) LIKE '%erf%' OR LOWER(gi.GENENAME) LIKE '%ap2%')",
    'HSF': "(LOWER(gi.GENENAME) LIKE '%heat stress transcription%' OR LOWER(gi.GENENAME) LIKE '%heat shock factor%')",
}

tf_families = {}
gene_names = {}

for family, where_clause in tf_queries.items():
    cursor.execute(f"SELECT gi.SYMBOL, gi.GENENAME FROM gene_info gi WHERE {where_clause}")
    results = {r[0]: r[1] for r in cursor.fetchall()}
    # Filter to expressed genes
    expressed_tf = {g: n for g, n in results.items() if g in expressed}
    tf_families[family] = expressed_tf
    gene_names.update(expressed_tf)
    print(f"  {family}: {len(results)} total, {len(expressed_tf)} expressed")

# Also get all gene names for context
cursor.execute("SELECT SYMBOL, GENENAME FROM gene_info")
all_gene_names = {r[0]: r[1] for r in cursor.fetchall()}
conn.close()

# Total TF genes
all_tf_genes = set()
for fam, genes in tf_families.items():
    all_tf_genes.update(genes.keys())
print(f"\nTotal unique TF genes: {len(all_tf_genes)}")


# ============================================================
# 3. DE ANALYSIS ACROSS CONDITIONS
# ============================================================
print("\n" + "=" * 70)
print("DE ANALYSIS FOR TF FAMILIES ACROSS CONDITIONS")
print("=" * 70)

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
seasons = [2021, 2022]
stage_labels = {'Pre_veraison': 'Pre-vér', 'Veraison': 'Véraison',
                'Post_veraison': 'Post-vér', 'Harvest': 'Harvest'}

def get_samples(**kwargs):
    mask = pd.Series(True, index=meta.index)
    for col, val in kwargs.items():
        if isinstance(val, list):
            mask &= meta[col].isin(val)
        else:
            mask &= meta[col] == val
    return list(meta[mask].index)


def run_de(group1, group2):
    """Welch's t-test DE with BH correction."""
    g1 = log_cpm[group1]
    g2 = log_cpm[group2]
    results = []
    for gene in log_cpm.index:
        v1 = g1.loc[gene].values
        v2 = g2.loc[gene].values
        if np.std(v1) == 0 and np.std(v2) == 0:
            continue
        t, p = stats.ttest_ind(v1, v2, equal_var=False)
        results.append({'gene': gene, 'log2FC': v1.mean() - v2.mean(), 'pvalue': p})
    df = pd.DataFrame(results)
    if len(df) > 0:
        df['padj'] = np.minimum(df['pvalue'] * len(df) / rankdata(df['pvalue']), 1.0)
    return df.sort_values('pvalue')


# Run DE for each season × stage condition
condition_de = {}
condition_degs = {}

for season in seasons:
    for stage in stages:
        label = f"{season}_{stage}"
        infected = get_samples(Infection_Group='RB(+)', Sampling_Season=season, Maturity_stage=stage)
        control = get_samples(Infection_Group='RB(-)', Sampling_Season=season, Maturity_stage=stage)

        if len(infected) < 5 or len(control) < 3:
            continue

        de = run_de(infected, control)
        condition_de[label] = de

        degs = de[(de['padj'] < 0.05) & (de['log2FC'].abs() >= 1)]
        condition_degs[label] = set(degs['gene'])
        print(f"  {label}: {len(degs)} DEGs total")


# ============================================================
# 4. TF FAMILY DE ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("TF FAMILY DE ENRICHMENT")
print("=" * 70)

summary_rows = []
detail_rows = []

for condition, degs in condition_degs.items():
    de_df = condition_de[condition]
    n_total_degs = len(degs)

    for family, tf_genes in tf_families.items():
        tf_set = set(tf_genes.keys())
        tf_degs = tf_set & degs
        n_tf = len(tf_set)
        n_tf_deg = len(tf_degs)

        # Fisher's exact test for enrichment
        a = n_tf_deg
        b = n_total_degs - a
        c = n_tf - a
        d = len(expressed) - n_total_degs - c

        if a + c > 0:
            _, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        else:
            pval = 1.0

        # Direction bias
        if n_tf_deg > 0:
            tf_de = de_df[de_df['gene'].isin(tf_degs)]
            n_up = (tf_de['log2FC'] > 0).sum()
            n_down = (tf_de['log2FC'] < 0).sum()
            mean_fc = tf_de['log2FC'].mean()

            # Module score: mean z-score across TF family
            tf_expr_infected = log_cpm.loc[list(tf_set)].mean(axis=0)
            infected_samples = get_samples(
                Infection_Group='RB(+)',
                Sampling_Season=int(condition.split('_')[0]),
                Maturity_stage='_'.join(condition.split('_')[1:]))
            control_samples = get_samples(
                Infection_Group='RB(-)',
                Sampling_Season=int(condition.split('_')[0]),
                Maturity_stage='_'.join(condition.split('_')[1:]))

            if infected_samples and control_samples:
                inf_score = tf_expr_infected[infected_samples].mean()
                ctrl_score = tf_expr_infected[control_samples].mean()
                module_diff = inf_score - ctrl_score
            else:
                module_diff = 0
        else:
            n_up, n_down, mean_fc, module_diff = 0, 0, 0, 0

        summary_rows.append({
            'condition': condition,
            'tf_family': family,
            'n_tf_total': n_tf,
            'n_tf_deg': n_tf_deg,
            'n_up': n_up,
            'n_down': n_down,
            'mean_log2FC': mean_fc,
            'module_diff': module_diff,
            'enrichment_pval': pval,
            'percent_de': n_tf_deg / n_tf * 100 if n_tf > 0 else 0,
        })

        # Detail: each TF gene's DE results
        for tf_gene in tf_set:
            gene_de = de_df[de_df['gene'] == tf_gene]
            if len(gene_de) > 0:
                r = gene_de.iloc[0]
                detail_rows.append({
                    'condition': condition,
                    'tf_family': family,
                    'gene': tf_gene,
                    'gene_name': tf_genes.get(tf_gene, ''),
                    'log2FC': r['log2FC'],
                    'pvalue': r['pvalue'],
                    'padj': r['padj'],
                    'is_deg': abs(r['log2FC']) >= 1 and r['padj'] < 0.05,
                })

summary_df = pd.DataFrame(summary_rows)
detail_df = pd.DataFrame(detail_rows)

summary_df.to_csv('tf_family_DE_summary.csv', index=False)
detail_df.to_csv('tf_family_gene_details.csv', index=False)
print(f"Saved tf_family_DE_summary.csv ({len(summary_df)} rows)")
print(f"Saved tf_family_gene_details.csv ({len(detail_df)} rows)")


# ============================================================
# 5. FIGURE 43 — TF Family Heatmap Across Conditions
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 43: TF FAMILY HEATMAP")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(22, 8))
fig.suptitle('Transcription Factor Family Disruption by GRBV',
             fontsize=20, fontweight='bold', y=0.98)

# Panel A: Percent DE heatmap
ax = axes[0]
families = list(tf_families.keys())
conditions = sorted(condition_degs.keys())

pct_matrix = pd.DataFrame(0.0, index=families, columns=conditions)
for _, row in summary_df.iterrows():
    if row['condition'] in conditions:
        pct_matrix.loc[row['tf_family'], row['condition']] = row['percent_de']

col_labels = [f"{c.split('_')[0]}\n{'_'.join(c.split('_')[1:])}" for c in conditions]
sns.heatmap(pct_matrix.values, ax=ax,
            xticklabels=col_labels, yticklabels=families,
            cmap='YlOrRd', annot=True, fmt='.0f',
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': '% DE genes'})
ax.set_title('A. Percent DE Genes per TF Family', fontsize=16, fontweight='bold')
ax.tick_params(axis='x', rotation=45, labelsize=10)
ax.tick_params(axis='y', labelsize=12)

# Panel B: Mean log2FC heatmap (direction)
ax = axes[1]
fc_matrix = pd.DataFrame(0.0, index=families, columns=conditions)
for _, row in summary_df.iterrows():
    if row['condition'] in conditions:
        fc_matrix.loc[row['tf_family'], row['condition']] = row['mean_log2FC']

sns.heatmap(fc_matrix.values, ax=ax,
            xticklabels=col_labels, yticklabels=families,
            cmap='RdBu_r', center=0, vmin=-2, vmax=2,
            annot=True, fmt='.2f',
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': 'Mean log2FC'})
ax.set_title('B. Mean log2FC per TF Family', fontsize=16, fontweight='bold')
ax.tick_params(axis='x', rotation=45, labelsize=10)
ax.tick_params(axis='y', labelsize=12)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('43_tf_family_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 43_tf_family_heatmap.png")


# ============================================================
# 6. FIGURE 44 — TF Family DEG Enrichment Balloon Plot
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 44: TF FAMILY ENRICHMENT BALLOON PLOT")
print("=" * 70)

fig, ax = plt.subplots(figsize=(16, 10))

# Balloon plot: x = condition, y = family, size = n_tf_deg, color = direction
for i, family in enumerate(families):
    for j, cond in enumerate(conditions):
        row = summary_df[(summary_df['tf_family'] == family) & (summary_df['condition'] == cond)]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        n_deg = row['n_tf_deg']
        if n_deg == 0:
            continue

        # Size proportional to number of DE TFs
        size = n_deg * 40

        # Color by direction bias
        if row['n_up'] > row['n_down']:
            color = '#e74c3c'  # red for upregulated
        elif row['n_down'] > row['n_up']:
            color = '#2980b9'  # blue for downregulated
        else:
            color = '#95a5a6'  # grey for mixed

        # Alpha by significance
        alpha = 0.9 if row['enrichment_pval'] < 0.05 else 0.4

        ax.scatter(j, i, s=size, c=color, alpha=alpha, edgecolors='black',
                   linewidth=0.5, zorder=3)
        ax.text(j, i, str(n_deg), ha='center', va='center', fontsize=8, fontweight='bold')

ax.set_xticks(range(len(conditions)))
ax.set_xticklabels([f"{c.split('_')[0]}\n{'_'.join(c.split('_')[1:])}" for c in conditions],
                    fontsize=10, rotation=45)
ax.set_yticks(range(len(families)))
ax.set_yticklabels(families, fontsize=13)
ax.set_title('TF Family DEG Enrichment\n(size = count, red = up, blue = down)',
             fontsize=18, fontweight='bold')
ax.set_xlabel('Condition', fontsize=14)
ax.set_ylabel('TF Family', fontsize=14)
ax.grid(True, alpha=0.2)
ax.set_xlim(-0.5, len(conditions) - 0.5)
ax.set_ylim(-0.5, len(families) - 0.5)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c',
           markersize=12, label='Predominantly UP'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2980b9',
           markersize=12, label='Predominantly DOWN'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#95a5a6',
           markersize=12, label='Mixed'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=11)

plt.tight_layout()
plt.savefig('44_tf_enrichment_balloon.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 44_tf_enrichment_balloon.png")


# ============================================================
# 7. FIGURE 45 — Top TF Gene Expression Boxplots
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 45: TOP TF EXPRESSION BOXPLOTS")
print("=" * 70)

# Select top 2 most disrupted genes per family
fig, axes = plt.subplots(2, 3, figsize=(22, 14))
fig.suptitle('Most Disrupted Transcription Factor Genes by GRBV',
             fontsize=20, fontweight='bold', y=0.98)

for fam_idx, family in enumerate(families):
    ax = axes[fam_idx // 3, fam_idx % 3]
    tf_genes_list = list(tf_families[family].keys())

    # Find most consistently significant TF genes across conditions
    gene_scores = defaultdict(float)
    for cond, de in condition_de.items():
        for gene in tf_genes_list:
            gene_de = de[de['gene'] == gene]
            if len(gene_de) > 0:
                r = gene_de.iloc[0]
                if r['padj'] < 0.05 and abs(r['log2FC']) >= 0.5:
                    gene_scores[gene] += abs(r['log2FC'])

    if not gene_scores:
        ax.text(0.5, 0.5, f'{family}\nNo significant genes', ha='center', va='center',
                transform=ax.transAxes, fontsize=12)
        ax.set_title(family, fontsize=14, fontweight='bold')
        continue

    top_tf = sorted(gene_scores, key=gene_scores.get, reverse=True)[:4]

    # Build boxplot data
    plot_data = []
    for gene in top_tf:
        for season in seasons:
            for inf in ['RB(+)', 'RB(-)']:
                mask = (meta['Sampling_Season'] == season) & (meta['Infection_Group'] == inf)
                samples = meta[mask].index.tolist()
                if samples:
                    vals = log_cpm.loc[gene, samples].values
                    for v in vals:
                        short_name = gene_names.get(gene, gene)[:25]
                        plot_data.append({
                            'Gene': f"{gene}\n({short_name})" if len(short_name) < 20 else gene,
                            'Group': f"{season} {inf}",
                            'Expression': v
                        })

    if plot_data:
        pdf = pd.DataFrame(plot_data)
        palette = {'2021 RB(+)': '#e74c3c', '2021 RB(-)': '#2980b9',
                   '2022 RB(+)': '#c0392b', '2022 RB(-)': '#3498db'}
        sns.boxplot(data=pdf, x='Gene', y='Expression', hue='Group',
                    palette=palette, ax=ax, linewidth=0.8, fliersize=2)
        ax.set_title(f'{family} Family', fontsize=14, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('log2(CPM+1)', fontsize=11)
        ax.tick_params(axis='x', labelsize=8, rotation=15)
        ax.legend(fontsize=7, loc='best')

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('45_tf_expression_boxplots.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 45_tf_expression_boxplots.png")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("TF FAMILY ANALYSIS COMPLETE")
print("=" * 70)

print("\nTF Family Summary:")
for family in families:
    tf_set = set(tf_families[family].keys())
    n_ever_de = 0
    for cond, degs in condition_degs.items():
        n_ever_de += len(tf_set & degs)
    print(f"  {family}: {len(tf_set)} expressed, "
          f"{n_ever_de} total DE events across conditions")

print(f"\nOutput files:")
print(f"  - 43_tf_family_heatmap.png")
print(f"  - 44_tf_enrichment_balloon.png")
print(f"  - 45_tf_expression_boxplots.png")
print(f"  - tf_family_DE_summary.csv")
print(f"  - tf_family_gene_details.csv")
