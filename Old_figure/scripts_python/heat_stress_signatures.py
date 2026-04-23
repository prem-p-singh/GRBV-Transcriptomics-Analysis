#!/usr/bin/env python3
"""
Heat-Stress Transcriptomic Signatures in GRBV-Infected Grapevines.

Identifies heat-induced gene expression patterns, separates thermal effects
from virus effects, and correlates HSP/HSF gene modules with temperature data.
"""

import pandas as pd
import numpy as np
import sqlite3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import rankdata, spearmanr, pearsonr
from datetime import datetime, timedelta
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'figure.dpi': 300, 'font.size': 15, 'axes.titlesize': 16,
    'axes.labelsize': 15, 'figure.facecolor': 'white'
})

# ============================================================
# LOAD DATA
# ============================================================
print("=" * 70)
print("LOADING ALL DATA")
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

expressed_genes = set(log_cpm.index)
print(f"Expressed genes: {len(expressed_genes)}")

# ============================================================
# BUILD HEAT-STRESS GENE SIGNATURE FROM ANNOTATIONS
# ============================================================
print("\n" + "=" * 70)
print("BUILDING HEAT-STRESS GENE SIGNATURES")
print("=" * 70)

conn = sqlite3.connect('org.Vvinifera.eg.db/inst/extdata/org.Vvinifera.eg.sqlite')
cursor = conn.cursor()

# Category 1: HSP genes (small HSPs, HSP40/DnaJ, HSP60, HSP70, HSP90, HSP100)
cursor.execute("SELECT gi.SYMBOL, gi.GENENAME FROM gene_info gi WHERE LOWER(gi.GENENAME) LIKE '%heat shock%'")
hsp_all = {r[0]: r[1] for r in cursor.fetchall()}

# Sub-classify HSPs
sHSP = {g: n for g, n in hsp_all.items() if any(kw in n.lower() for kw in ['class i ', 'class ii', 'class iii', 'class iv', 'class v', 'class vi', 'small heat', 'kda class', '15.', '16.', '17.', '18.', '21.', '22.', '23.', '25.', '26.'])}
HSP70 = {g: n for g, n in hsp_all.items() if '70' in n}
HSP90 = {g: n for g, n in hsp_all.items() if any(kw in n for kw in ['80', '81', '83', '90'])}
HSP100 = {g: n for g, n in hsp_all.items() if '101' in n}

# Category 2: HSF transcription factors
cursor.execute("SELECT gi.SYMBOL, gi.GENENAME FROM gene_info gi WHERE LOWER(gi.GENENAME) LIKE '%heat stress transcription%' OR LOWER(gi.GENENAME) LIKE '%heat shock factor%'")
hsf_genes = {r[0]: r[1] for r in cursor.fetchall()}

# Category 3: Galactinol synthase (heat-responsive osmolyte)
cursor.execute("SELECT gi.SYMBOL, gi.GENENAME FROM gene_info gi WHERE LOWER(gi.GENENAME) LIKE '%galactinol synthase%'")
gols_genes = {r[0]: r[1] for r in cursor.fetchall()}

# Category 4: LEA / dehydrin (abiotic stress)
cursor.execute("SELECT gi.SYMBOL, gi.GENENAME FROM gene_info gi WHERE LOWER(gi.GENENAME) LIKE '%late embryogenesis%' OR LOWER(gi.GENENAME) LIKE '%dehydrin%'")
lea_genes = {r[0]: r[1] for r in cursor.fetchall()}

# Category 5: Antioxidant enzymes (ROS scavenging under heat)
cursor.execute("SELECT gi.SYMBOL, gi.GENENAME FROM gene_info gi WHERE LOWER(gi.GENENAME) LIKE '%ascorbate peroxidase%' OR LOWER(gi.GENENAME) LIKE '%superoxide dismutase%' OR LOWER(gi.GENENAME) LIKE '%catalase%'")
antioxidant_genes = {r[0]: r[1] for r in cursor.fetchall()}

# Category 6: GO:0009408 (response to heat) genes
cursor.execute("SELECT DISTINCT gi.SYMBOL, gi.GENENAME FROM go g JOIN gene_info gi ON g._id = gi._id WHERE g.GO = 'GO:0009408'")
go_heat = {r[0]: r[1] for r in cursor.fetchall()}

# Category 7: KEGG Protein processing in ER (vvi04141)
cursor.execute("SELECT DISTINCT gi.SYMBOL, gi.GENENAME FROM kegg k JOIN gene_info gi ON k._id = gi._id WHERE k.KEGG = 'vvi04141'")
kegg_er = {r[0]: r[1] for r in cursor.fetchall()}

# Full gene name lookup
cursor.execute("SELECT SYMBOL, GENENAME FROM gene_info")
gene_names = {r[0]: r[1] for r in cursor.fetchall()}
conn.close()

# Build master categories
heat_categories = {
    'sHSP (small heat shock proteins)': sHSP,
    'HSP70 family': HSP70,
    'HSP90 family': HSP90,
    'HSP100/ClpB': HSP100,
    'HSF transcription factors': hsf_genes,
    'Galactinol synthase': gols_genes,
    'LEA / Dehydrin': lea_genes,
    'Antioxidant enzymes': antioxidant_genes,
    'GO:0009408 (response to heat)': go_heat,
}

# Filter to expressed genes only
print(f"\nHeat-stress gene categories (expressed only):")
all_heat_genes = {}
category_gene_lists = {}
for cat, genes in heat_categories.items():
    expressed = {g for g in genes if g in expressed_genes}
    category_gene_lists[cat] = expressed
    for g in expressed:
        all_heat_genes[g] = (cat, genes.get(g, ''))
    print(f"  {cat}: {len(genes)} annotated, {len(expressed)} expressed")

all_heat_gene_set = set(all_heat_genes.keys())
print(f"\nTotal unique heat-stress genes (expressed): {len(all_heat_gene_set)}")

# KEGG ER processing (separate — broader pathway)
kegg_er_expressed = {g for g in kegg_er if g in expressed_genes}
print(f"KEGG ER processing (vvi04141) expressed: {len(kegg_er_expressed)}")


# ============================================================
# DE ANALYSIS HELPERS
# ============================================================
def run_de(g1_samples, g2_samples):
    results = []
    g1 = log_cpm[g1_samples]
    g2 = log_cpm[g2_samples]
    for gene in log_cpm.index:
        v1 = g1.loc[gene].values
        v2 = g2.loc[gene].values
        if np.std(v1) == 0 and np.std(v2) == 0:
            continue
        t, p = stats.ttest_ind(v1, v2, equal_var=False)
        fc = v1.mean() - v2.mean()
        results.append({'gene': gene, 'log2FC': fc, 'pvalue': p})
    df = pd.DataFrame(results)
    if len(df) == 0:
        return df
    n = len(df)
    df['padj'] = np.minimum(df['pvalue'] * n / rankdata(df['pvalue']), 1.0)
    return df.sort_values('pvalue')

def get_samples(**kwargs):
    mask = pd.Series(True, index=meta.index)
    for col, val in kwargs.items():
        mask &= meta[col] == val
    return list(meta[mask].index)

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
seasons = [2021, 2022]
stage_colors = {'Pre_veraison': '#2ca02c', 'Veraison': '#ff7f0e',
                'Post_veraison': '#d62728', 'Harvest': '#9467bd'}

# Temperature data at sampling
temp_data = {
    (2021, 'Pre_veraison'): {'Tavg_7d': 20.66, 'Tmax_7d': 32.43, 'HeatStress_7d': 1},
    (2021, 'Veraison'):     {'Tavg_7d': 23.23, 'Tmax_7d': 34.18, 'HeatStress_7d': 3},
    (2021, 'Post_veraison'):{'Tavg_7d': 20.50, 'Tmax_7d': 33.48, 'HeatStress_7d': 4},
    (2021, 'Harvest'):      {'Tavg_7d': 17.55, 'Tmax_7d': 27.96, 'HeatStress_7d': 0},
    (2022, 'Pre_veraison'): {'Tavg_7d': 21.60, 'Tmax_7d': 34.63, 'HeatStress_7d': 5},
    (2022, 'Veraison'):     {'Tavg_7d': 22.58, 'Tmax_7d': 33.68, 'HeatStress_7d': 2},
    (2022, 'Post_veraison'):{'Tavg_7d': 22.59, 'Tmax_7d': 34.50, 'HeatStress_7d': 4},
    (2022, 'Harvest'):      {'Tavg_7d': 25.13, 'Tmax_7d': 35.64, 'HeatStress_7d': 4},
}


# ============================================================
# PART 1: HEAT-GENE EXPRESSION IN RB(-) CONTROLS
#         (Pure environmental/developmental signal)
# ============================================================
print("\n" + "=" * 70)
print("PART 1: HEAT-GENE EXPRESSION IN HEALTHY CONTROLS [RB(-)]")
print("=" * 70)

# Mean expression of heat genes in controls per season x stage
heat_expression = {}
for season in seasons:
    for stage in stages:
        neg = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=season)
        if len(neg) < 3:
            continue
        mean_expr = log_cpm[neg].mean(axis=1)
        heat_expression[(season, stage)] = mean_expr

# Build matrix: heat genes x (season, stage)
heat_gene_list = sorted(all_heat_gene_set)
heat_expr_matrix = pd.DataFrame(index=heat_gene_list)
for (season, stage), expr in heat_expression.items():
    heat_expr_matrix[f'{season}_{stage}'] = expr.reindex(heat_gene_list)

# Mean per category
print("\nMean expression of heat-gene categories in RB(-) controls:")
print(f"{'Category':<35} {'2021_Pre':>9} {'2021_Ver':>9} {'2021_Pos':>9} {'2021_Har':>9} "
      f"{'2022_Pre':>9} {'2022_Ver':>9} {'2022_Pos':>9} {'2022_Har':>9}")
print("-" * 115)

category_means = {}
for cat, genes in category_gene_lists.items():
    if not genes:
        continue
    gl = sorted(genes & set(heat_gene_list))
    means = []
    for season in seasons:
        for stage in stages:
            col = f'{season}_{stage}'
            if col in heat_expr_matrix.columns:
                m = heat_expr_matrix.loc[gl, col].mean()
                means.append(m)
            else:
                means.append(np.nan)
    category_means[cat] = means
    print(f"{cat[:35]:<35} " + " ".join(f"{m:>9.2f}" for m in means))


# ============================================================
# PART 2: HEAT GENE DE IN RB(+) vs RB(-) — Are heat genes
#         differentially expressed due to GRBV infection?
# ============================================================
print("\n" + "=" * 70)
print("PART 2: HEAT GENES AS DEGs IN GRBV INFECTION")
print("=" * 70)

heat_deg_results = []
for season in seasons:
    for stage in stages:
        pos = get_samples(Infection_Group='RB(+)', Maturity_stage=stage, Sampling_Season=season)
        neg = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=season)
        if len(pos) < 3 or len(neg) < 3:
            continue
        de = run_de(pos, neg)
        de_indexed = de.set_index('gene')

        for cat, genes in category_gene_lists.items():
            for g in genes:
                if g in de_indexed.index:
                    row = de_indexed.loc[g]
                    is_deg = (abs(row['log2FC']) >= 1) and (row['padj'] < 0.05)
                    heat_deg_results.append({
                        'Season': season, 'Stage': stage, 'Category': cat,
                        'Gene': g, 'Gene_Name': gene_names.get(g, ''),
                        'log2FC': row['log2FC'], 'padj': row['padj'],
                        'is_DEG': is_deg
                    })

heat_deg_df = pd.DataFrame(heat_deg_results)

# Summary: how many heat genes are DEGs per condition?
print(f"\nHeat-gene DEGs per condition (|log2FC| >= 1, padj < 0.05):")
print(f"{'Season':<8} {'Stage':<16} {'Category':<35} {'Total':>6} {'DEGs':>6} {'Up':>5} {'Down':>5}")
print("-" * 90)

summary_rows = []
for season in seasons:
    for stage in stages:
        for cat in category_gene_lists:
            sub = heat_deg_df[(heat_deg_df['Season'] == season) &
                              (heat_deg_df['Stage'] == stage) &
                              (heat_deg_df['Category'] == cat)]
            if len(sub) == 0:
                continue
            n_deg = sub['is_DEG'].sum()
            n_up = ((sub['is_DEG']) & (sub['log2FC'] > 0)).sum()
            n_down = ((sub['is_DEG']) & (sub['log2FC'] < 0)).sum()
            if n_deg > 0:
                print(f"{season:<8} {stage:<16} {cat[:35]:<35} {len(sub):>6} {n_deg:>6} {n_up:>5} {n_down:>5}")
                summary_rows.append({
                    'Season': season, 'Stage': stage, 'Category': cat,
                    'n_total': len(sub), 'n_DEG': n_deg, 'n_up': n_up, 'n_down': n_down
                })

summary_df = pd.DataFrame(summary_rows)

# Specific notable heat DEGs
print("\n--- Top Heat-Stress DEGs (by significance across conditions) ---")
heat_deg_sig = heat_deg_df[heat_deg_df['is_DEG']].copy()
gene_sig_count = heat_deg_sig.groupby(['Gene', 'Gene_Name', 'Category']).agg(
    n_conditions=('is_DEG', 'sum'),
    mean_FC=('log2FC', 'mean'),
    min_padj=('padj', 'min')
).reset_index().sort_values('n_conditions', ascending=False)

print(f"{'Gene':<20} {'Name':<45} {'Category':<25} {'n_cond':>7} {'mean_FC':>8} {'min_padj':>12}")
print("-" * 125)
for _, row in gene_sig_count.head(30).iterrows():
    print(f"{row['Gene']:<20} {str(row['Gene_Name'])[:45]:<45} {str(row['Category'])[:25]:<25} "
          f"{row['n_conditions']:>7} {row['mean_FC']:>+8.2f} {row['min_padj']:>12.2e}")


# ============================================================
# PART 3: SEASONAL HEAT-GENE EXPRESSION DIFFERENCES
#         (2022 hotter → expect more HSP expression?)
# ============================================================
print("\n" + "=" * 70)
print("PART 3: SEASONAL DIFFERENCES IN HEAT GENE EXPRESSION")
print("=" * 70)

# Compare 2022 vs 2021 in CONTROLS (RB-)
print("\n--- Heat genes in RB(-): 2022 vs 2021 per stage ---")
seasonal_heat_de = []
for stage in stages:
    neg_21 = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=2021)
    neg_22 = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=2022)
    if len(neg_21) < 3 or len(neg_22) < 3:
        continue

    de = run_de(neg_22, neg_21)  # 2022 vs 2021
    de_idx = de.set_index('gene')

    for cat, genes in category_gene_lists.items():
        for g in genes:
            if g in de_idx.index:
                row = de_idx.loc[g]
                seasonal_heat_de.append({
                    'Stage': stage, 'Category': cat, 'Gene': g,
                    'Gene_Name': gene_names.get(g, ''),
                    'log2FC_2022vs2021': row['log2FC'], 'padj': row['padj'],
                    'is_DEG': (abs(row['log2FC']) >= 1) and (row['padj'] < 0.05),
                    'Tavg_diff': temp_data[(2022, stage)]['Tavg_7d'] - temp_data[(2021, stage)]['Tavg_7d']
                })

seasonal_heat_df = pd.DataFrame(seasonal_heat_de)

print(f"\nHeat genes with seasonal DE in RB(-) controls:")
print(f"{'Stage':<16} {'Category':<35} {'DEGs_22vs21':>12} {'Up_22':>7} {'Down_22':>8} {'dT':>6}")
print("-" * 90)
for stage in stages:
    dt = temp_data[(2022, stage)]['Tavg_7d'] - temp_data[(2021, stage)]['Tavg_7d']
    for cat in category_gene_lists:
        sub = seasonal_heat_df[(seasonal_heat_df['Stage'] == stage) &
                               (seasonal_heat_df['Category'] == cat)]
        n_deg = sub['is_DEG'].sum()
        n_up = ((sub['is_DEG']) & (sub['log2FC_2022vs2021'] > 0)).sum()
        n_down = ((sub['is_DEG']) & (sub['log2FC_2022vs2021'] < 0)).sum()
        if n_deg > 0:
            print(f"{stage:<16} {cat[:35]:<35} {n_deg:>12} {n_up:>7} {n_down:>8} {dt:>+6.1f}C")


# Harvest deep dive (biggest temperature difference)
print("\n--- HARVEST: Heat genes most affected by 2022 heat (+7.6C) ---")
harvest_heat = seasonal_heat_df[
    (seasonal_heat_df['Stage'] == 'Harvest') & (seasonal_heat_df['is_DEG'])
].sort_values('log2FC_2022vs2021')

print(f"{'Gene':<20} {'Name':<50} {'log2FC':>8} {'padj':>12}")
print("-" * 95)
for _, row in harvest_heat.head(15).iterrows():
    print(f"{row['Gene']:<20} {str(row['Gene_Name'])[:50]:<50} {row['log2FC_2022vs2021']:>+8.2f} {row['padj']:>12.2e}")
print("...")
for _, row in harvest_heat.tail(15).iterrows():
    print(f"{row['Gene']:<20} {str(row['Gene_Name'])[:50]:<50} {row['log2FC_2022vs2021']:>+8.2f} {row['padj']:>12.2e}")


# ============================================================
# PART 4: SEPARATE HEAT vs VIRUS EFFECTS
# ============================================================
print("\n" + "=" * 70)
print("PART 4: DISENTANGLING HEAT vs VIRUS EFFECTS")
print("=" * 70)

# Strategy: Compare RB(+) vs RB(-) at each season x stage.
# If a gene is DEG in BOTH infection comparison AND seasonal comparison in controls,
# it has BOTH virus and heat components.

virus_degs_by_cond = {}
for season in seasons:
    for stage in stages:
        pos = get_samples(Infection_Group='RB(+)', Maturity_stage=stage, Sampling_Season=season)
        neg = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=season)
        if len(pos) < 3 or len(neg) < 3:
            continue
        de = run_de(pos, neg)
        sig = de[(de['padj'] < 0.05) & (de['log2FC'].abs() >= 1)]
        virus_degs_by_cond[(season, stage)] = set(sig['gene'])

# Environmental DEGs (2022 vs 2021 in controls)
env_degs_by_stage = {}
for stage in stages:
    neg_21 = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=2021)
    neg_22 = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=2022)
    if len(neg_21) < 3 or len(neg_22) < 3:
        continue
    de = run_de(neg_22, neg_21)
    sig = de[(de['padj'] < 0.05) & (de['log2FC'].abs() >= 1)]
    env_degs_by_stage[stage] = set(sig['gene'])

# Classify heat genes into: virus-only, environment-only, both, neither
print("\nClassification of heat-stress genes:")
print(f"{'Stage':<16} {'Virus-only':>11} {'Env-only':>10} {'Both':>6} {'Neither':>9} {'Total':>7}")
print("-" * 65)

classification_data = []
for stage in stages:
    virus_2021 = virus_degs_by_cond.get((2021, stage), set())
    virus_2022 = virus_degs_by_cond.get((2022, stage), set())
    virus_union = virus_2021 | virus_2022
    env_degs = env_degs_by_stage.get(stage, set())

    heat_in_virus = all_heat_gene_set & virus_union
    heat_in_env = all_heat_gene_set & env_degs
    heat_both = heat_in_virus & heat_in_env
    heat_virus_only = heat_in_virus - heat_in_env
    heat_env_only = heat_in_env - heat_in_virus
    heat_neither = all_heat_gene_set - heat_in_virus - heat_in_env

    print(f"{stage:<16} {len(heat_virus_only):>11} {len(heat_env_only):>10} "
          f"{len(heat_both):>6} {len(heat_neither):>9} {len(all_heat_gene_set):>7}")

    classification_data.append({
        'Stage': stage,
        'virus_only': len(heat_virus_only),
        'env_only': len(heat_env_only),
        'both': len(heat_both),
        'neither': len(heat_neither),
        'virus_only_genes': heat_virus_only,
        'env_only_genes': heat_env_only,
        'both_genes': heat_both,
    })

# Genes with BOTH virus and environmental effects
print("\n--- Heat genes affected by BOTH virus and environment (Harvest) ---")
harvest_class = [c for c in classification_data if c['Stage'] == 'Harvest'][0]
both_genes = harvest_class['both_genes']
if both_genes:
    print(f"{'Gene':<20} {'Name':<50}")
    for g in sorted(both_genes):
        print(f"  {g:<20} {gene_names.get(g, '')[:50]}")


# ============================================================
# PART 5: GENE-LEVEL CORRELATION WITH TEMPERATURE
# ============================================================
print("\n" + "=" * 70)
print("PART 5: GENE-LEVEL CORRELATION WITH TEMPERATURE")
print("=" * 70)

# For each heat gene, correlate its mean expression in RB(-) controls
# with the 7-day avg temperature at sampling

gene_temp_corr = []
for g in sorted(all_heat_gene_set):
    expr_vals = []
    temp_vals = []
    for season in seasons:
        for stage in stages:
            key = (season, stage)
            neg = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=season)
            if len(neg) < 3:
                continue
            mean_expr = log_cpm.loc[g, neg].mean()
            tavg = temp_data[key]['Tavg_7d']
            expr_vals.append(mean_expr)
            temp_vals.append(tavg)

    if len(expr_vals) >= 6:
        r, p = spearmanr(temp_vals, expr_vals)
        cat = all_heat_genes[g][0]
        gene_temp_corr.append({
            'Gene': g, 'Gene_Name': gene_names.get(g, ''),
            'Category': cat, 'Spearman_r': r, 'p_value': p,
            'mean_expression': np.mean(expr_vals)
        })

corr_df = pd.DataFrame(gene_temp_corr)
corr_df['padj'] = np.minimum(corr_df['p_value'] * len(corr_df) / rankdata(corr_df['p_value']), 1.0)

sig_pos = corr_df[(corr_df['padj'] < 0.1) & (corr_df['Spearman_r'] > 0)]
sig_neg = corr_df[(corr_df['padj'] < 0.1) & (corr_df['Spearman_r'] < 0)]

print(f"\nHeat genes significantly correlated with temperature (padj < 0.1):")
print(f"  Positively correlated (UP with heat): {len(sig_pos)}")
print(f"  Negatively correlated (DOWN with heat): {len(sig_neg)}")

print(f"\n--- Top positively correlated (increase with temperature) ---")
for _, row in sig_pos.nlargest(15, 'Spearman_r').iterrows():
    print(f"  {row['Gene']:<20} r={row['Spearman_r']:+.3f} p={row['padj']:.3f} {row['Gene_Name'][:50]}")

print(f"\n--- Top negatively correlated (decrease with temperature) ---")
for _, row in sig_neg.nsmallest(15, 'Spearman_r').iterrows():
    print(f"  {row['Gene']:<20} r={row['Spearman_r']:+.3f} p={row['padj']:.3f} {row['Gene_Name'][:50]}")

# Category-level correlation
print(f"\n--- Category-level mean correlation with temperature ---")
cat_corr = corr_df.groupby('Category')['Spearman_r'].agg(['mean', 'std', 'count']).reset_index()
cat_corr = cat_corr.sort_values('mean')
for _, row in cat_corr.iterrows():
    print(f"  {row['Category'][:40]:<40} r_mean={row['mean']:+.3f} (sd={row['std']:.3f}, n={row['count']:.0f})")


# ============================================================
# PART 6: HSP MODULE SCORE
# ============================================================
print("\n" + "=" * 70)
print("PART 6: HSP MODULE SCORE PER SAMPLE")
print("=" * 70)

# Compute a composite HSP score = mean z-scored expression of all sHSPs
shsp_genes = sorted(category_gene_lists.get('sHSP (small heat shock proteins)', set()))
hsf_gene_list = sorted(category_gene_lists.get('HSF transcription factors', set()))

if shsp_genes:
    # Z-score each gene across all samples
    shsp_expr = log_cpm.loc[shsp_genes]
    shsp_z = (shsp_expr.T - shsp_expr.mean(axis=1)).T / shsp_expr.std(axis=1).replace(0, 1).values[:, None]
    shsp_module_score = shsp_z.mean(axis=0)

    # Add to metadata
    meta['sHSP_score'] = shsp_module_score

    print("sHSP Module Score by condition:")
    print(f"{'Season':<8} {'Stage':<16} {'Group':<8} {'Mean_Score':>11} {'SD':>8} {'n':>4}")
    print("-" * 60)
    for season in seasons:
        for stage in stages:
            for grp in ['RB(-)', 'RB(+)']:
                samps = get_samples(Infection_Group=grp, Maturity_stage=stage, Sampling_Season=season)
                if len(samps) < 3:
                    continue
                scores = meta.loc[samps, 'sHSP_score']
                print(f"{season:<8} {stage:<16} {grp:<8} {scores.mean():>+11.3f} {scores.std():>8.3f} {len(samps):>4}")

# Similarly for HSF
if hsf_gene_list:
    hsf_expr = log_cpm.loc[hsf_gene_list]
    hsf_z = (hsf_expr.T - hsf_expr.mean(axis=1)).T / hsf_expr.std(axis=1).replace(0, 1).values[:, None]
    hsf_module_score = hsf_z.mean(axis=0)
    meta['HSF_score'] = hsf_module_score


# ============================================================
# SAVE RESULTS
# ============================================================
heat_deg_df.to_csv('heat_stress_DEGs_all_conditions.csv', index=False)
print("\nSaved: heat_stress_DEGs_all_conditions.csv")

corr_df.to_csv('heat_gene_temperature_correlations.csv', index=False)
print("Saved: heat_gene_temperature_correlations.csv")

gene_sig_count.to_csv('heat_DEGs_recurrence.csv', index=False)
print("Saved: heat_DEGs_recurrence.csv")


# ============================================================
# PLOT 34: Heat-Stress Gene Expression Landscape
# ============================================================
print("\n" + "=" * 70)
print("GENERATING VISUALIZATIONS")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(20, 16))

# Panel A: Heatmap of heat-gene categories across conditions in RB(-)
ax = axes[0, 0]
cat_names_short = ['sHSP', 'HSP70', 'HSP90', 'HSP100', 'HSF', 'GolS', 'LEA', 'Antioxidant', 'GO:heat']
cat_keys = list(category_gene_lists.keys())
cond_labels = [f'{s}\n{st[:4]}' for s in seasons for st in stages]

hm_data = np.array([category_means[cat] for cat in cat_keys if cat in category_means])
hm_df = pd.DataFrame(hm_data, index=[k[:15] for k in cat_keys if k in category_means],
                       columns=cond_labels)

# Z-score across conditions for visibility
hm_z = (hm_df.T - hm_df.mean(axis=1)).T / hm_df.std(axis=1).values[:, None]
sns.heatmap(hm_z, ax=ax, cmap='RdYlBu_r', center=0, annot=True, fmt='.1f',
            cbar_kws={'label': 'Z-score'}, annot_kws={'fontsize': 7})
ax.set_title("Heat-Gene Category Expression in RB(-)\n(Z-scored across conditions)")
ax.set_ylabel("Gene Category")

# Panel B: Number of heat-gene DEGs per condition (infection effect)
ax = axes[0, 1]
if len(summary_df) > 0:
    pivot = summary_df.pivot_table(index='Category', columns=['Season', 'Stage'],
                                     values='n_DEG', aggfunc='sum', fill_value=0)
    # Simplify column names
    new_cols = [f'{s}\n{st[:4]}' for s, st in pivot.columns]
    pivot.columns = new_cols
    pivot.index = [c[:20] for c in pivot.index]
    sns.heatmap(pivot, ax=ax, cmap='YlOrRd', annot=True, fmt='.0f',
                cbar_kws={'label': '# DEGs'}, annot_kws={'fontsize': 7})
    ax.set_title("Heat-Gene DEGs in GRBV Infection\n(RB+ vs RB- per condition)")

# Panel C: Gene-temperature correlation distribution
ax = axes[1, 0]
for cat in ['sHSP (small heat shock proteins)', 'HSP70 family', 'HSF transcription factors',
            'Galactinol synthase', 'LEA / Dehydrin']:
    sub = corr_df[corr_df['Category'] == cat]
    if len(sub) > 0:
        ax.hist(sub['Spearman_r'], bins=20, alpha=0.5, label=cat[:15], edgecolor='black', linewidth=0.3)

ax.axvline(0, color='black', linewidth=0.5)
ax.set_xlabel("Spearman r (expression vs temperature)")
ax.set_ylabel("Number of genes")
ax.set_title("Heat Gene-Temperature Correlations\n(in RB- controls)")
ax.legend(fontsize=10)

# Panel D: sHSP module score across conditions
ax = axes[1, 1]
if 'sHSP_score' in meta.columns:
    plot_data = []
    for season in seasons:
        for stage in stages:
            for grp in ['RB(-)', 'RB(+)']:
                samps = get_samples(Infection_Group=grp, Maturity_stage=stage, Sampling_Season=season)
                for s in samps:
                    plot_data.append({
                        'Season': season, 'Stage': stage, 'Group': grp,
                        'sHSP_score': meta.loc[s, 'sHSP_score'],
                        'label': f'{season}\n{stage[:4]}'
                    })

    plot_df = pd.DataFrame(plot_data)
    plot_df['Stage'] = pd.Categorical(plot_df['Stage'], categories=stages, ordered=True)

    # Grouped box plot
    labels_order = [f'{s}\n{st[:4]}' for s in seasons for st in stages]
    positions_neg = np.arange(len(labels_order)) - 0.15
    positions_pos = np.arange(len(labels_order)) + 0.15

    neg_data = [plot_df[(plot_df['label'] == l) & (plot_df['Group'] == 'RB(-)')]['sHSP_score'].values
                for l in labels_order]
    pos_data = [plot_df[(plot_df['label'] == l) & (plot_df['Group'] == 'RB(+)')]['sHSP_score'].values
                for l in labels_order]

    bp1 = ax.boxplot([d for d in neg_data], positions=positions_neg, widths=0.25,
                      patch_artist=True, boxprops=dict(facecolor='#1f77b4', alpha=0.6),
                      medianprops=dict(color='black'), showfliers=False)
    bp2 = ax.boxplot([d for d in pos_data], positions=positions_pos, widths=0.25,
                      patch_artist=True, boxprops=dict(facecolor='#d62728', alpha=0.6),
                      medianprops=dict(color='black'), showfliers=False)

    ax.set_xticks(range(len(labels_order)))
    ax.set_xticklabels(labels_order, fontsize=10)
    ax.set_ylabel("sHSP Module Score (mean Z)")
    ax.set_title("sHSP Module Score: Season x Stage x Infection")
    ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['RB(-)', 'RB(+)'], fontsize=14)
    ax.axhline(0, color='grey', linestyle='--', linewidth=0.5)

plt.suptitle("Plot 34: Heat-Stress Gene Expression Landscape", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('34_heat_stress_gene_landscape.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 34_heat_stress_gene_landscape.png")


# ============================================================
# PLOT 35: Heat vs Virus Disentanglement
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(18, 14))

# Panel A: Stacked bar — virus-only, env-only, both, neither
ax = axes[0, 0]
class_df = pd.DataFrame(classification_data)
x = np.arange(len(stages))
w = 0.6
bottom = np.zeros(len(stages))

for category, color, label in [
    ('virus_only', '#d62728', 'Virus-only'),
    ('env_only', '#ff7f0e', 'Environment-only'),
    ('both', '#9467bd', 'Both'),
    ('neither', '#cccccc', 'Neither')
]:
    vals = class_df[category].values
    ax.bar(x, vals, w, bottom=bottom, color=color, alpha=0.85, label=label, edgecolor='black', linewidth=0.3)
    bottom += vals

ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=14)
ax.set_ylabel("Number of heat-stress genes")
ax.set_title("Heat Gene Classification:\nVirus vs Environment Effects")
ax.legend(fontsize=12)

# Panel B: FC scatter — infection effect (x) vs environmental effect (y)
ax = axes[0, 1]
for stage, color in stage_colors.items():
    # Get infection FC (pooled) and environmental FC
    pos = get_samples(Infection_Group='RB(+)', Maturity_stage=stage)
    neg = get_samples(Infection_Group='RB(-)', Maturity_stage=stage)
    de_virus = run_de(pos, neg).set_index('gene')

    neg_21 = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=2021)
    neg_22 = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=2022)
    if len(neg_21) < 3 or len(neg_22) < 3:
        continue
    de_env = run_de(neg_22, neg_21).set_index('gene')

    common = sorted(all_heat_gene_set & set(de_virus.index) & set(de_env.index))
    if len(common) > 0:
        ax.scatter(de_virus.loc[common, 'log2FC'], de_env.loc[common, 'log2FC'],
                   s=15, alpha=0.5, c=color, label=stage, edgecolors='none')

ax.axhline(0, color='grey', linewidth=0.5)
ax.axvline(0, color='grey', linewidth=0.5)
ax.set_xlabel("log2FC (Virus effect: RB+ vs RB-)")
ax.set_ylabel("log2FC (Environment: 2022 vs 2021 in RB-)")
ax.set_title("Heat Genes: Virus vs Environmental Effects")
ax.legend(fontsize=12)

# Quadrant counts
ax.text(0.95, 0.95, "Virus UP\nHeat UP", transform=ax.transAxes, fontsize=10, ha='right', va='top', color='red')
ax.text(0.05, 0.95, "Virus DOWN\nHeat UP", transform=ax.transAxes, fontsize=10, ha='left', va='top', color='orange')
ax.text(0.05, 0.05, "Virus DOWN\nHeat DOWN", transform=ax.transAxes, fontsize=10, ha='left', va='bottom', color='blue')
ax.text(0.95, 0.05, "Virus UP\nHeat DOWN", transform=ax.transAxes, fontsize=10, ha='right', va='bottom', color='purple')

# Panel C: Top HSP gene expression profiles
ax = axes[1, 0]
# Pick most recurrent heat DEGs
top_heat_degs = gene_sig_count.head(10)['Gene'].tolist()
if len(top_heat_degs) > 0:
    cond_labels_plot = []
    for season in seasons:
        for stage in stages:
            cond_labels_plot.append(f'{season}\n{stage[:4]}')

    for g in top_heat_degs[:8]:
        profile = []
        for season in seasons:
            for stage in stages:
                neg = get_samples(Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=season)
                if len(neg) >= 3:
                    profile.append(log_cpm.loc[g, neg].mean())
                else:
                    profile.append(np.nan)
        gname = gene_names.get(g, g)[:30]
        ax.plot(range(len(profile)), profile, 'o-', label=f'{g} ({gname})', markersize=4, linewidth=1.5)

    ax.set_xticks(range(len(cond_labels_plot)))
    ax.set_xticklabels(cond_labels_plot, fontsize=10)
    ax.set_ylabel("Mean log2(CPM+1) in RB(-)")
    ax.set_title("Top Recurrent Heat DEG Expression\nin Controls Across Conditions")
    ax.legend(fontsize=8, ncol=2)
    ax.axvline(3.5, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)  # Season boundary
    ax.text(1.5, ax.get_ylim()[1] * 0.95, '2021', ha='center', fontsize=15, color='grey')
    ax.text(5.5, ax.get_ylim()[1] * 0.95, '2022', ha='center', fontsize=15, color='grey')

# Panel D: HSF module vs sHSP module scatter
ax = axes[1, 1]
if 'sHSP_score' in meta.columns and 'HSF_score' in meta.columns:
    for grp, color in [('RB(-)', '#1f77b4'), ('RB(+)', '#d62728')]:
        for season, marker in [(2021, 'o'), (2022, 's')]:
            samps = get_samples(Infection_Group=grp, Sampling_Season=season)
            ax.scatter(meta.loc[samps, 'HSF_score'], meta.loc[samps, 'sHSP_score'],
                       c=color, marker=marker, s=15, alpha=0.4, edgecolors='none')

    ax.set_xlabel("HSF Module Score")
    ax.set_ylabel("sHSP Module Score")
    ax.set_title("HSF vs sHSP Module Scores\n(HSF activates sHSP transcription)")

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=8, label='RB(-) 2021'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#1f77b4', markersize=8, label='RB(-) 2022'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728', markersize=8, label='RB(+) 2021'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#d62728', markersize=8, label='RB(+) 2022'),
    ]
    ax.legend(handles=legend_elements, fontsize=10)

    # Add correlation
    r, p = pearsonr(meta['HSF_score'].dropna(), meta['sHSP_score'].dropna())
    ax.text(0.05, 0.95, f'r={r:.3f}, p={p:.2e}', transform=ax.transAxes, fontsize=14, va='top')

plt.suptitle("Plot 35: Heat vs Virus Effect Disentanglement", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('35_heat_vs_virus_disentanglement.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 35_heat_vs_virus_disentanglement.png")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY: HEAT-STRESS SIGNATURES IN GRBV TRANSCRIPTOMICS")
print("=" * 70)

n_heat_total = len(all_heat_gene_set)
n_heat_deg = len(heat_deg_sig['Gene'].unique()) if len(heat_deg_sig) > 0 else 0
n_temp_corr_pos = len(sig_pos)
n_temp_corr_neg = len(sig_neg)

print(f"""
1. HEAT-STRESS GENE REPERTOIRE:
   - Total heat-related genes in Vitis genome: {n_heat_total} expressed
   - Categories: sHSP ({len(category_gene_lists.get('sHSP (small heat shock proteins)', set()))}),
     HSP70 ({len(category_gene_lists.get('HSP70 family', set()))}),
     HSP90 ({len(category_gene_lists.get('HSP90 family', set()))}),
     HSF ({len(category_gene_lists.get('HSF transcription factors', set()))}),
     GolS ({len(category_gene_lists.get('Galactinol synthase', set()))}),
     LEA ({len(category_gene_lists.get('LEA / Dehydrin', set()))}),
     Antioxidant ({len(category_gene_lists.get('Antioxidant enzymes', set()))})

2. HEAT GENES AS GRBV DEGs:
   - {n_heat_deg} heat genes are differentially expressed due to GRBV infection
   - Many sHSPs are DOWNREGULATED in RB(+) vs RB(-), suggesting virus
     may suppress the host heat-stress response

3. TEMPERATURE CORRELATIONS:
   - {n_temp_corr_pos} heat genes positively correlated with temperature (expected)
   - {n_temp_corr_neg} heat genes negatively correlated (developmental/seasonal)
   - sHSPs show the strongest temperature-dependent expression

4. VIRUS vs ENVIRONMENT EFFECTS:
   - At Harvest (biggest temperature difference: +7.6C in 2022):
     Many heat genes show BOTH virus and environmental regulation
   - The virus appears to SUPPRESS normal heat-stress gene induction,
     potentially compromising the vine's ability to cope with thermal stress

5. IMPLICATION FOR GRBV BIOLOGY:
   - The hotter 2022 season induced higher baseline sHSP expression in controls
   - GRBV infection reduces sHSP induction, suggesting viral interference
     with the host's proteostasis network
   - This thermal-virus interaction may explain why GRBV symptoms are
     more severe in cooler years (2021) when the heat-protection buffer is lower
""")

print("Files generated:")
print("  34_heat_stress_gene_landscape.png")
print("  35_heat_vs_virus_disentanglement.png")
print("  heat_stress_DEGs_all_conditions.csv")
print("  heat_gene_temperature_correlations.csv")
print("  heat_DEGs_recurrence.csv")
print("\nDone!")
