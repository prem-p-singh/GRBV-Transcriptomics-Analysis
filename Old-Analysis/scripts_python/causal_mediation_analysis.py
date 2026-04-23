#!/usr/bin/env python3
"""
Causal Mediation Analysis: Temperature → Titer → Transcriptome.

Tests whether the effect of temperature on gene expression is mediated through
viral titer, providing causal pathway hypotheses for the GRBV-temperature
interaction.

Outputs:
  - Fig 49: Mediation path diagrams (top genes)
  - Fig 50: Direct vs indirect effect scatter
  - mediation_results.csv
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import rankdata, spearmanr
from datetime import datetime
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'figure.dpi': 300, 'font.size': 15, 'axes.titlesize': 16,
    'axes.labelsize': 15, 'figure.facecolor': 'white'
})

# ============================================================
# 1. LOAD EXPRESSION DATA
# ============================================================
print("=" * 70)
print("CAUSAL MEDIATION ANALYSIS: LOADING DATA")
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

# ============================================================
# 2. LOAD WEATHER DATA
# ============================================================
print("\nLoading weather data...")

sampling_dates = {
    2021: {
        'Pre_veraison': datetime(2021, 8, 2),
        'Veraison':     datetime(2021, 8, 16),
        'Post_veraison': datetime(2021, 9, 9),
        'Harvest':      datetime(2021, 9, 21),
    },
    2022: {
        'Pre_veraison': datetime(2022, 7, 28),
        'Veraison':     datetime(2022, 8, 9),
        'Post_veraison': datetime(2022, 8, 24),
        'Harvest':      datetime(2022, 9, 13),
    }
}

def load_cimis(filepath, year):
    df = pd.read_csv(filepath)
    df.columns = [c.strip() for c in df.columns]
    df['Date'] = pd.to_datetime(df['Date'].astype(str).str.strip(), format='%m/%d/%Y')
    result = pd.DataFrame()
    result['Date'] = df['Date']
    for orig, new in {'Max Air Temp (C)': 'Tmax', 'Min Air Temp (C)': 'Tmin',
                       'Avg Air Temp (C)': 'Tavg'}.items():
        if orig in df.columns:
            result[new] = pd.to_numeric(df[orig], errors='coerce')
    result['Season'] = year
    return result

weather_2021 = load_cimis('daily_2021.csv', 2021)
weather_2022 = load_cimis('daily_2022.csv', 2022)

# Compute 7-day mean temperature at each sampling date
temp_at_sampling = {}
for season, dates in sampling_dates.items():
    weather = weather_2021 if season == 2021 else weather_2022
    for stage, date in dates.items():
        mask = (weather['Date'] >= date - pd.Timedelta(days=7)) & (weather['Date'] <= date)
        if mask.sum() > 0:
            temp_at_sampling[(season, stage)] = weather.loc[mask, 'Tavg'].mean()
        else:
            temp_at_sampling[(season, stage)] = np.nan

print("Temperature at sampling (7-day mean):")
for (s, st), t in sorted(temp_at_sampling.items()):
    print(f"  {s} {st}: {t:.1f}°C")


# ============================================================
# 3. LOAD TITER DATA
# ============================================================
print("\nLoading virus titer data...")

try:
    titer_raw = pd.read_csv('Virus titer.txt', sep='\t')
    print(f"Titer data: {len(titer_raw)} rows, columns: {list(titer_raw.columns)}")

    # Parse titer data — adapt to actual column names
    titer_cols = titer_raw.columns.tolist()

    # Try to find sample and titer columns
    titer_data = {}
    if 'Sample' in titer_cols:
        for _, row in titer_raw.iterrows():
            sample = str(row['Sample']).strip()
            # Find numeric titer column
            for col in titer_cols:
                if col != 'Sample' and pd.api.types.is_numeric_dtype(titer_raw[col]):
                    titer_data[sample] = float(row[col])
                    break
    else:
        # Try first column as sample, second as titer
        for _, row in titer_raw.iterrows():
            sample = str(row.iloc[0]).strip()
            try:
                titer_data[sample] = float(row.iloc[1])
            except (ValueError, IndexError):
                pass

    print(f"Parsed titer for {len(titer_data)} samples")
except Exception as e:
    print(f"Titer loading failed: {e}")
    titer_data = {}


# ============================================================
# 4. MAP TEMPERATURE + TITER TO INFECTED SAMPLES
# ============================================================
print("\n" + "=" * 70)
print("MAPPING VARIABLES TO SAMPLES")
print("=" * 70)

# Only RB(+) samples
infected_meta = meta[meta['Infection_Group'] == 'RB(+)'].copy()

# Map temperature
infected_meta['Temperature'] = infected_meta.apply(
    lambda row: temp_at_sampling.get(
        (int(row['Sampling_Season']), row['Maturity_stage']), np.nan),
    axis=1)

# Map titer (try direct sample match)
infected_meta['Titer'] = infected_meta.index.map(
    lambda s: titer_data.get(s, np.nan))

# If no direct match, try partial matching
if infected_meta['Titer'].isna().all() and titer_data:
    print("  No direct titer match, trying fuzzy matching...")
    for sample in infected_meta.index:
        for titer_key, titer_val in titer_data.items():
            if sample in titer_key or titer_key in sample:
                infected_meta.loc[sample, 'Titer'] = titer_val
                break

n_with_temp = infected_meta['Temperature'].notna().sum()
n_with_titer = infected_meta['Titer'].notna().sum()
n_complete = (infected_meta['Temperature'].notna() & infected_meta['Titer'].notna()).sum()

print(f"Infected samples: {len(infected_meta)}")
print(f"  With temperature: {n_with_temp}")
print(f"  With titer: {n_with_titer}")
print(f"  Complete cases: {n_complete}")


# ============================================================
# 5. MEDIATION ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("RUNNING MEDIATION ANALYSIS")
print("=" * 70)

# Use samples with temperature data (all infected should have it)
analysis_samples = infected_meta[infected_meta['Temperature'].notna()].index.tolist()
has_titer = n_with_titer > 10

if has_titer:
    analysis_samples_full = infected_meta[
        infected_meta['Temperature'].notna() & infected_meta['Titer'].notna()
    ].index.tolist()
    print(f"Full mediation analysis: {len(analysis_samples_full)} samples")
else:
    analysis_samples_full = analysis_samples
    print(f"Temperature-only analysis (insufficient titer): {len(analysis_samples)} samples")

# For each gene, compute:
# Path c:  Temperature → Gene (total effect)
# Path a:  Temperature → Titer
# Path b:  Titer → Gene (controlling for Temperature)
# Indirect: a × b
# Direct:   c - indirect (c')

from numpy.linalg import lstsq

def simple_regression(X, y):
    """OLS regression, returns (coef, pvalue, r_squared)."""
    X = np.column_stack([np.ones(len(X)), X])
    beta, residuals, _, _ = lstsq(X, y, rcond=None)
    y_pred = X @ beta
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    n = len(y)
    p = X.shape[1]
    mse = ss_res / (n - p) if n > p else 1e-10
    se = np.sqrt(np.diag(mse * np.linalg.pinv(X.T @ X)))
    t_stat = beta[1] / se[1] if se[1] > 0 else 0
    p_val = 2 * (1 - stats.t.cdf(abs(t_stat), n - p))
    return beta[1], p_val, r2


def multiple_regression(X_multi, y):
    """Multiple OLS regression. X_multi: (n, k). Returns coefs, pvals."""
    X = np.column_stack([np.ones(len(y)), X_multi])
    beta, _, _, _ = lstsq(X, y, rcond=None)
    y_pred = X @ beta
    ss_res = np.sum((y - y_pred) ** 2)
    n, p = X.shape
    mse = ss_res / (n - p) if n > p else 1e-10
    try:
        se = np.sqrt(np.diag(mse * np.linalg.pinv(X.T @ X)))
    except Exception:
        return beta[1:], np.ones(len(beta) - 1)
    t_stats = beta[1:] / se[1:]
    p_vals = 2 * (1 - stats.t.cdf(np.abs(t_stats), n - p))
    return beta[1:], p_vals


# Temperature variable
T = infected_meta.loc[analysis_samples, 'Temperature'].values

# Select top variable genes for mediation analysis
var_genes = log_cpm[analysis_samples].var(axis=1).nlargest(500).index.tolist()

mediation_results = []

if has_titer and len(analysis_samples_full) > 20:
    T_full = infected_meta.loc[analysis_samples_full, 'Temperature'].values
    M_full = infected_meta.loc[analysis_samples_full, 'Titer'].values

    # Path a: Temperature → Titer
    a_coef, a_pval, a_r2 = simple_regression(T_full, M_full)
    print(f"\nPath a (Temp → Titer): β={a_coef:.4f}, p={a_pval:.4f}, R²={a_r2:.4f}")

    for gene in var_genes:
        Y = log_cpm.loc[gene, analysis_samples_full].values

        # Path c: Temperature → Gene (total)
        c_coef, c_pval, c_r2 = simple_regression(T_full, Y)

        # Path b + c': Titer → Gene, controlling for Temperature
        betas, pvals = multiple_regression(np.column_stack([T_full, M_full]), Y)
        c_prime = betas[0]  # direct effect of temp
        b_coef = betas[1]   # effect of titer controlling for temp
        b_pval = pvals[1]

        # Indirect effect
        indirect = a_coef * b_coef
        # Proportion mediated
        prop_mediated = indirect / c_coef if abs(c_coef) > 1e-10 else 0

        mediation_results.append({
            'gene': gene,
            'total_effect_c': c_coef,
            'total_pval': c_pval,
            'direct_effect_c_prime': c_prime,
            'indirect_effect_ab': indirect,
            'path_a': a_coef,
            'path_b': b_coef,
            'path_b_pval': b_pval,
            'prop_mediated': np.clip(prop_mediated, -1, 1),
            'temp_r2': c_r2,
        })
else:
    print("\nInsufficient titer data for full mediation. Running temperature correlation only.")
    for gene in var_genes:
        Y = log_cpm.loc[gene, analysis_samples].values
        c_coef, c_pval, c_r2 = simple_regression(T, Y)

        mediation_results.append({
            'gene': gene,
            'total_effect_c': c_coef,
            'total_pval': c_pval,
            'direct_effect_c_prime': c_coef,
            'indirect_effect_ab': 0,
            'path_a': 0,
            'path_b': 0,
            'path_b_pval': 1.0,
            'prop_mediated': 0,
            'temp_r2': c_r2,
        })

med_df = pd.DataFrame(mediation_results)

# BH correction
if len(med_df) > 0:
    med_df['total_padj'] = np.minimum(
        med_df['total_pval'] * len(med_df) / rankdata(med_df['total_pval']), 1.0)

# Classify genes
med_df['category'] = 'Not significant'
med_df.loc[med_df['total_padj'] < 0.05, 'category'] = 'Temperature-responsive'
if has_titer:
    med_df.loc[
        (med_df['total_padj'] < 0.05) & (med_df['prop_mediated'].abs() > 0.3),
        'category'] = 'Titer-mediated'
    med_df.loc[
        (med_df['total_padj'] < 0.05) & (med_df['prop_mediated'].abs() <= 0.3) &
        (med_df['path_b_pval'] > 0.05),
        'category'] = 'Direct temperature'

med_df = med_df.sort_values('total_pval')
med_df.to_csv('mediation_results.csv', index=False)
print(f"\nSaved mediation_results.csv ({len(med_df)} genes)")

# Summarize
for cat in med_df['category'].unique():
    n = (med_df['category'] == cat).sum()
    print(f"  {cat}: {n} genes")


# ============================================================
# 6. LOAD GENE NAMES FOR PLOTS
# ============================================================
import sqlite3
conn = sqlite3.connect('org.Vvinifera.eg.db/inst/extdata/org.Vvinifera.eg.sqlite')
gene_info = pd.read_sql("SELECT SYMBOL, GENENAME FROM gene_info", conn)
gene_name_map = dict(zip(gene_info['SYMBOL'], gene_info['GENENAME']))
conn.close()


# ============================================================
# 7. FIGURE 49 — Mediation Path Diagrams
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 49: MEDIATION PATH DIAGRAMS")
print("=" * 70)

fig, axes = plt.subplots(2, 3, figsize=(24, 16))
fig.suptitle('Causal Mediation: Temperature → Titer → Gene Expression',
             fontsize=20, fontweight='bold', y=0.98)

# Top 6 genes by total effect significance
top_genes = med_df.head(6)

for idx, (_, row) in enumerate(top_genes.iterrows()):
    ax = axes[idx // 3, idx % 3]
    gene = row['gene']
    gname = gene_name_map.get(gene, '')[:35]

    # Draw mediation diagram
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')

    # Boxes
    from matplotlib.patches import FancyBboxPatch
    # Temperature box
    temp_box = FancyBboxPatch((0.5, 4), 2.5, 1.5, boxstyle="round,pad=0.2",
                               facecolor='#f39c12', alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.add_patch(temp_box)
    ax.text(1.75, 4.75, 'Temperature', ha='center', va='center', fontsize=12, fontweight='bold')

    # Gene box
    gene_box = FancyBboxPatch((7, 4), 2.5, 1.5, boxstyle="round,pad=0.2",
                               facecolor='#2ecc71', alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.add_patch(gene_box)
    ax.text(8.25, 4.75, gene[:15], ha='center', va='center', fontsize=11, fontweight='bold')

    # Titer box (mediator)
    titer_box = FancyBboxPatch((3.5, 7.5), 2.5, 1.5, boxstyle="round,pad=0.2",
                                facecolor='#e74c3c', alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.add_patch(titer_box)
    ax.text(4.75, 8.25, 'Viral Titer', ha='center', va='center', fontsize=12, fontweight='bold')

    # Arrows with coefficients
    # Path a: Temp → Titer
    ax.annotate('', xy=(3.5, 8), xytext=(3, 5.5),
                arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=2))
    ax.text(2.5, 7, f'a={row["path_a"]:.3f}', fontsize=10, rotation=50,
            ha='center', color='#2c3e50')

    # Path b: Titer → Gene
    ax.annotate('', xy=(7, 5.5), xytext=(6, 8),
                arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=2))
    ax.text(7, 7, f'b={row["path_b"]:.3f}', fontsize=10, rotation=-50,
            ha='center', color='#2c3e50')

    # Path c' (direct): Temp → Gene
    ax.annotate('', xy=(7, 4.75), xytext=(3, 4.75),
                arrowprops=dict(arrowstyle='->', color='#e74c3c', lw=2.5, linestyle='--'))
    ax.text(5, 4, f"c'={row['direct_effect_c_prime']:.3f}", fontsize=11,
            ha='center', color='#e74c3c', fontweight='bold')

    # Total effect and mediation info
    ax.text(5, 2.5, f"Total (c) = {row['total_effect_c']:.3f}  (p={row['total_pval']:.2e})",
            ha='center', fontsize=10)
    ax.text(5, 1.8, f"Indirect (a×b) = {row['indirect_effect_ab']:.3f}",
            ha='center', fontsize=10)
    prop = row['prop_mediated']
    ax.text(5, 1.1, f"% Mediated = {prop*100:.1f}%", ha='center', fontsize=11,
            fontweight='bold', color='#8e44ad')

    # Gene name
    ax.text(5, 0.3, gname, ha='center', fontsize=9, style='italic', color='grey')

    ax.set_title(f'{gene}', fontsize=13, fontweight='bold', pad=5)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('49_mediation_diagrams.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 49_mediation_diagrams.png")


# ============================================================
# 8. FIGURE 50 — Direct vs Indirect Effect Scatter
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 50: DIRECT vs INDIRECT EFFECTS")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(20, 9))
fig.suptitle('Temperature Effects on Gene Expression: Direct vs Titer-Mediated',
             fontsize=20, fontweight='bold', y=0.98)

# Panel A: Scatter of direct vs indirect effects
ax = axes[0]
sig = med_df[med_df['total_padj'] < 0.05]

if len(sig) > 0:
    cat_colors = {
        'Temperature-responsive': '#f39c12',
        'Titer-mediated': '#e74c3c',
        'Direct temperature': '#2980b9',
        'Not significant': '#bdc3c7',
    }

    for cat, color in cat_colors.items():
        sub = sig[sig['category'] == cat]
        if len(sub) > 0:
            ax.scatter(sub['direct_effect_c_prime'], sub['indirect_effect_ab'],
                      c=color, alpha=0.5, s=30, label=f"{cat} ({len(sub)})",
                      edgecolors='none')

    ax.axhline(0, color='grey', linewidth=0.5)
    ax.axvline(0, color='grey', linewidth=0.5)
    ax.set_xlabel("Direct effect (c')", fontsize=14)
    ax.set_ylabel('Indirect effect (a×b)', fontsize=14)
    ax.set_title("A. Direct vs Indirect Temperature Effects", fontsize=16, fontweight='bold')
    ax.legend(fontsize=10, loc='best')

    # Annotate top genes
    top5 = sig.nsmallest(5, 'total_pval')
    for _, r in top5.iterrows():
        ax.annotate(r['gene'][:12], (r['direct_effect_c_prime'], r['indirect_effect_ab']),
                   fontsize=8, alpha=0.7)
else:
    ax.text(0.5, 0.5, 'No significant genes', ha='center', va='center',
            transform=ax.transAxes, fontsize=14)
    ax.set_title("A. Direct vs Indirect Temperature Effects", fontsize=16, fontweight='bold')

# Panel B: Category pie/bar chart
ax = axes[1]
cat_counts = med_df['category'].value_counts()

if len(cat_counts) > 0:
    colors_pie = [cat_colors.get(c, '#bdc3c7') for c in cat_counts.index]
    bars = ax.barh(range(len(cat_counts)), cat_counts.values, color=colors_pie, alpha=0.85,
                    edgecolor='white')
    ax.set_yticks(range(len(cat_counts)))
    ax.set_yticklabels(cat_counts.index, fontsize=12)
    ax.set_xlabel('Number of genes', fontsize=14)
    ax.set_title('B. Gene Classification by Mediation', fontsize=16, fontweight='bold')
    ax.invert_yaxis()

    for i, v in enumerate(cat_counts.values):
        ax.text(v + max(cat_counts.values) * 0.02, i, str(v), va='center', fontsize=12)

    # Summary stats
    n_temp = (med_df['total_padj'] < 0.05).sum()
    n_total = len(med_df)
    ax.text(0.95, 0.95, f'{n_temp}/{n_total} genes\ntemperature-responsive\n(padj < 0.05)',
            transform=ax.transAxes, ha='right', va='top', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('50_direct_vs_indirect_effects.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 50_direct_vs_indirect_effects.png")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("CAUSAL MEDIATION ANALYSIS COMPLETE")
print("=" * 70)

n_sig = (med_df['total_padj'] < 0.05).sum()
print(f"\nTemperature-responsive genes: {n_sig}/{len(med_df)}")

if has_titer:
    n_mediated = (med_df['category'] == 'Titer-mediated').sum()
    n_direct = (med_df['category'] == 'Direct temperature').sum()
    print(f"Titer-mediated: {n_mediated}")
    print(f"Direct temperature: {n_direct}")

print(f"\nOutput files:")
print(f"  - 49_mediation_diagrams.png")
print(f"  - 50_direct_vs_indirect_effects.png")
print(f"  - mediation_results.csv")
