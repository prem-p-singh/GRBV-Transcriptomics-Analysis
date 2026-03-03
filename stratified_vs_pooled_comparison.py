#!/usr/bin/env python3
"""
Stratified vs Pooled DE Comparison for GRBV Transcriptomics.

Re-runs Q1 (infection duration) and Q2 (maturity stages) with SEASON-STRATIFIED
comparisons (2021 RB+ vs 2021 RB- only, 2022 RB+ vs 2022 RB- only), then compares
against the original POOLED results and cross-references with Q3/Q4.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import rankdata
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# LOAD AND PREPARE DATA (same as analysis.py)
# ============================================================
print("=" * 70)
print("LOADING DATA")
print("=" * 70)

counts = pd.read_csv('counts_matrix.txt', sep='\t', index_col=0)
meta = pd.read_csv('sample_metadata.txt', sep='\t')

common_samples = list(set(counts.columns) & set(meta['Sample']))
counts = counts[common_samples]
meta = meta[meta['Sample'].isin(common_samples)].set_index('Sample').loc[common_samples]

# Filtering
min_samples = int(0.1 * counts.shape[1])
gene_mask = (counts >= 10).sum(axis=1) >= min_samples
counts_filtered = counts.loc[gene_mask]

# CPM + log2
lib_sizes = counts_filtered.sum(axis=0)
cpm = counts_filtered.div(lib_sizes, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

print(f"Genes after filtering: {counts_filtered.shape[0]}")
print(f"Samples: {counts_filtered.shape[1]}")

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def run_de_analysis(group1_samples, group2_samples, log_cpm_data):
    """DE using Welch's t-test with BH correction."""
    results = []
    g1 = log_cpm_data[group1_samples]
    g2 = log_cpm_data[group2_samples]
    for gene in log_cpm_data.index:
        vals1 = g1.loc[gene].values
        vals2 = g2.loc[gene].values
        if np.std(vals1) == 0 and np.std(vals2) == 0:
            continue
        t_stat, p_val = stats.ttest_ind(vals1, vals2, equal_var=False)
        log2fc = vals1.mean() - vals2.mean()
        results.append({'gene': gene, 'log2FC': log2fc, 'pvalue': p_val,
                        'mean_group1': vals1.mean(), 'mean_group2': vals2.mean()})
    df = pd.DataFrame(results)
    if len(df) == 0:
        return df
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


def get_degs(de_df, fc_thresh=1.0, padj_thresh=0.05):
    """Get set of DEG gene IDs."""
    sig = de_df[(de_df['padj'] < padj_thresh) & (de_df['log2FC'].abs() >= fc_thresh)]
    return set(sig['gene'])


def get_samples(meta_df, **kwargs):
    """Get sample names matching criteria."""
    mask = pd.Series(True, index=meta_df.index)
    for col, val in kwargs.items():
        if isinstance(val, list):
            mask &= meta_df[col].isin(val)
        else:
            mask &= meta_df[col] == val
    return list(meta_df[mask].index)


plt.rcParams.update({
    'figure.dpi': 300, 'font.size': 15, 'axes.titlesize': 16,
    'axes.labelsize': 15, 'figure.facecolor': 'white'
})

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
seasons = [2021, 2022]
ysi_values = sorted([y for y in meta['Years_Since_Infection'].unique() if y > 0])

# ============================================================
# PART 1: ORIGINAL POOLED Q1 — Infection Years vs ALL RB(-)
# ============================================================
print("\n" + "=" * 70)
print("PART 1A: ORIGINAL POOLED Q1 (all RB- controls pooled)")
print("=" * 70)

all_controls = get_samples(meta, Infection_Group='RB(-)')
print(f"Total RB(-) controls (both seasons pooled): {len(all_controls)}")

q1_pooled = {}
for ysi in ysi_values:
    infected = get_samples(meta, Infection_Group='RB(+)', Years_Since_Infection=ysi)
    de = run_de_analysis(infected, all_controls, log_cpm)
    total, up, down = count_deg(de)
    q1_pooled[ysi] = {'de_df': de, 'total': total, 'up': up, 'down': down,
                       'n_pos': len(infected), 'n_neg': len(all_controls),
                       'degs': get_degs(de)}
    print(f"  Inf.Year={ysi}: RB+={len(infected)} vs ALL RB-={len(all_controls)}, "
          f"DEGs={total} (up={up}, down={down})")

# ============================================================
# PART 1B: SEASON-STRATIFIED Q1
# ============================================================
print("\n" + "=" * 70)
print("PART 1B: SEASON-STRATIFIED Q1 (season-matched RB- controls)")
print("=" * 70)

q1_stratified = {}
for season in seasons:
    controls_season = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season)
    print(f"\n--- Season {season}: RB(-) controls = {len(controls_season)} ---")

    # Get infection years available in this season
    season_ysi = sorted([y for y in meta[(meta['Infection_Group'] == 'RB(+)') &
                         (meta['Sampling_Season'] == season)]['Years_Since_Infection'].unique()])

    for ysi in season_ysi:
        infected = get_samples(meta, Infection_Group='RB(+)',
                               Years_Since_Infection=ysi, Sampling_Season=season)
        if len(infected) < 3:
            print(f"  {season} Inf.Year={ysi}: Skipped (n={len(infected)})")
            continue
        de = run_de_analysis(infected, controls_season, log_cpm)
        total, up, down = count_deg(de)
        q1_stratified[(season, ysi)] = {'de_df': de, 'total': total, 'up': up, 'down': down,
                                         'n_pos': len(infected), 'n_neg': len(controls_season),
                                         'degs': get_degs(de)}
        print(f"  {season} Inf.Year={ysi}: RB+={len(infected)} vs {season} RB-={len(controls_season)}, "
              f"DEGs={total} (up={up}, down={down})")

# ============================================================
# PART 2A: ORIGINAL POOLED Q2 — Maturity Stages (seasons pooled)
# ============================================================
print("\n" + "=" * 70)
print("PART 2A: ORIGINAL POOLED Q2 (all seasons pooled per stage)")
print("=" * 70)

q2_pooled = {}
for stage in stages:
    pos = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=stage)
    neg = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=stage)
    de = run_de_analysis(pos, neg, log_cpm)
    total, up, down = count_deg(de)
    q2_pooled[stage] = {'de_df': de, 'total': total, 'up': up, 'down': down,
                         'n_pos': len(pos), 'n_neg': len(neg), 'degs': get_degs(de)}
    print(f"  {stage}: RB+={len(pos)} vs RB-={len(neg)} (both seasons), "
          f"DEGs={total} (up={up}, down={down})")

# ============================================================
# PART 2B: SEASON-STRATIFIED Q2
# ============================================================
print("\n" + "=" * 70)
print("PART 2B: SEASON-STRATIFIED Q2 (season-matched per stage)")
print("=" * 70)

q2_stratified = {}
for season in seasons:
    print(f"\n--- Season {season} ---")
    for stage in stages:
        pos = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=stage, Sampling_Season=season)
        neg = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=season)
        if len(pos) < 3 or len(neg) < 3:
            print(f"  {season} {stage}: Skipped (RB+={len(pos)}, RB-={len(neg)})")
            continue
        de = run_de_analysis(pos, neg, log_cpm)
        total, up, down = count_deg(de)
        q2_stratified[(season, stage)] = {'de_df': de, 'total': total, 'up': up, 'down': down,
                                           'n_pos': len(pos), 'n_neg': len(neg), 'degs': get_degs(de)}
        print(f"  {season} {stage}: RB+={len(pos)} vs RB-={len(neg)}, "
              f"DEGs={total} (up={up}, down={down})")

# ============================================================
# PART 3: Q3 — Season-level RB(+) vs RB(-) (already stratified)
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Q3 SEASONAL COMPARISONS (already season-stratified)")
print("=" * 70)

q3_results = {}
for season in seasons:
    pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season)
    neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season)
    de = run_de_analysis(pos, neg, log_cpm)
    total, up, down = count_deg(de)
    q3_results[season] = {'de_df': de, 'total': total, 'up': up, 'down': down,
                           'n_pos': len(pos), 'n_neg': len(neg), 'degs': get_degs(de)}
    print(f"  {season}: RB+={len(pos)} vs RB-={len(neg)}, DEGs={total} (up={up}, down={down})")

# ============================================================
# PART 4: Q4 — Environmental / Variance Partitioning
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Q3-DETAILED (Season x Stage — already stratified)")
print("=" * 70)

q4_results = {}
for season in seasons:
    for stage in stages:
        pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season, Maturity_stage=stage)
        neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season, Maturity_stage=stage)
        if len(pos) < 3 or len(neg) < 3:
            print(f"  {season}-{stage}: Skipped (RB+={len(pos)}, RB-={len(neg)})")
            continue
        de = run_de_analysis(pos, neg, log_cpm)
        total, up, down = count_deg(de)
        q4_results[(season, stage)] = {'de_df': de, 'total': total, 'up': up, 'down': down,
                                        'degs': get_degs(de)}
        print(f"  {season}-{stage}: DEGs={total} (up={up}, down={down})")


# ============================================================
# COMPARISON ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("COMPARISON: POOLED vs STRATIFIED")
print("=" * 70)

# --- Q1 Comparison ---
print("\n--- Q1: Infection Year — Pooled vs Stratified ---")
print(f"{'Inf.Year':<12} {'Pooled DEGs':>12} {'2021 Strat':>12} {'2022 Strat':>12} "
      f"{'Pooled∩2021':>12} {'Pooled∩2022':>12} {'2021-only':>10} {'2022-only':>10}")
print("-" * 100)

q1_comparison_data = []
for ysi in ysi_values:
    pooled_total = q1_pooled[ysi]['total'] if ysi in q1_pooled else 0
    pooled_degs = q1_pooled[ysi]['degs'] if ysi in q1_pooled else set()

    strat_2021 = q1_stratified.get((2021, ysi), {})
    strat_2022 = q1_stratified.get((2022, ysi), {})

    total_2021 = strat_2021.get('total', '-')
    total_2022 = strat_2022.get('total', '-')
    degs_2021 = strat_2021.get('degs', set())
    degs_2022 = strat_2022.get('degs', set())

    overlap_p21 = len(pooled_degs & degs_2021) if degs_2021 else '-'
    overlap_p22 = len(pooled_degs & degs_2022) if degs_2022 else '-'
    unique_2021 = len(degs_2021 - pooled_degs) if degs_2021 else '-'
    unique_2022 = len(degs_2022 - pooled_degs) if degs_2022 else '-'

    print(f"{ysi:<12} {pooled_total:>12} {str(total_2021):>12} {str(total_2022):>12} "
          f"{str(overlap_p21):>12} {str(overlap_p22):>12} {str(unique_2021):>10} {str(unique_2022):>10}")

    q1_comparison_data.append({
        'Infection_Year': ysi, 'Pooled_DEGs': pooled_total,
        'Strat_2021_DEGs': total_2021 if isinstance(total_2021, int) else np.nan,
        'Strat_2022_DEGs': total_2022 if isinstance(total_2022, int) else np.nan,
        'Overlap_Pooled_2021': overlap_p21 if isinstance(overlap_p21, int) else np.nan,
        'Overlap_Pooled_2022': overlap_p22 if isinstance(overlap_p22, int) else np.nan,
        'Unique_2021': unique_2021 if isinstance(unique_2021, int) else np.nan,
        'Unique_2022': unique_2022 if isinstance(unique_2022, int) else np.nan,
    })

# Fold-change concordance for Q1
print("\n--- Q1: Fold-Change Concordance (Pooled vs Stratified) ---")
for ysi in ysi_values:
    pooled_de = q1_pooled[ysi]['de_df'].set_index('gene') if ysi in q1_pooled else None
    for season in seasons:
        key = (season, ysi)
        if key not in q1_stratified or pooled_de is None:
            continue
        strat_de = q1_stratified[key]['de_df'].set_index('gene')
        common = list(set(pooled_de.index) & set(strat_de.index))
        if len(common) < 10:
            continue
        r, p = stats.pearsonr(pooled_de.loc[common, 'log2FC'], strat_de.loc[common, 'log2FC'])
        print(f"  Inf.Year={ysi}, Pooled vs {season}-Stratified: r={r:.4f} (p={p:.2e}), n={len(common)} genes")

# --- Q2 Comparison ---
print("\n--- Q2: Maturity Stage — Pooled vs Stratified ---")
print(f"{'Stage':<16} {'Pooled DEGs':>12} {'2021 Strat':>12} {'2022 Strat':>12} "
      f"{'Pooled∩2021':>12} {'Pooled∩2022':>12} {'2021-only':>10} {'2022-only':>10}")
print("-" * 100)

q2_comparison_data = []
for stage in stages:
    pooled_total = q2_pooled[stage]['total']
    pooled_degs = q2_pooled[stage]['degs']

    strat_2021 = q2_stratified.get((2021, stage), {})
    strat_2022 = q2_stratified.get((2022, stage), {})

    total_2021 = strat_2021.get('total', '-')
    total_2022 = strat_2022.get('total', '-')
    degs_2021 = strat_2021.get('degs', set())
    degs_2022 = strat_2022.get('degs', set())

    overlap_p21 = len(pooled_degs & degs_2021) if degs_2021 else '-'
    overlap_p22 = len(pooled_degs & degs_2022) if degs_2022 else '-'
    unique_2021 = len(degs_2021 - pooled_degs) if degs_2021 else '-'
    unique_2022 = len(degs_2022 - pooled_degs) if degs_2022 else '-'

    # Union of both stratified
    union_strat = degs_2021 | degs_2022
    pooled_missed = len(union_strat - pooled_degs)
    strat_missed = len(pooled_degs - union_strat)

    print(f"{stage:<16} {pooled_total:>12} {str(total_2021):>12} {str(total_2022):>12} "
          f"{str(overlap_p21):>12} {str(overlap_p22):>12} {str(unique_2021):>10} {str(unique_2022):>10}")

    q2_comparison_data.append({
        'Stage': stage, 'Pooled_DEGs': pooled_total,
        'Strat_2021_DEGs': total_2021 if isinstance(total_2021, int) else np.nan,
        'Strat_2022_DEGs': total_2022 if isinstance(total_2022, int) else np.nan,
        'Overlap_Pooled_2021': overlap_p21 if isinstance(overlap_p21, int) else np.nan,
        'Overlap_Pooled_2022': overlap_p22 if isinstance(overlap_p22, int) else np.nan,
        'Unique_2021': unique_2021 if isinstance(unique_2021, int) else np.nan,
        'Unique_2022': unique_2022 if isinstance(unique_2022, int) else np.nan,
        'Pooled_missed_by_strat': strat_missed,
        'Strat_missed_by_pooled': pooled_missed,
    })

# Fold-change concordance for Q2
print("\n--- Q2: Fold-Change Concordance (Pooled vs Stratified) ---")
for stage in stages:
    pooled_de = q2_pooled[stage]['de_df'].set_index('gene')
    for season in seasons:
        key = (season, stage)
        if key not in q2_stratified:
            continue
        strat_de = q2_stratified[key]['de_df'].set_index('gene')
        common = list(set(pooled_de.index) & set(strat_de.index))
        if len(common) < 10:
            continue
        r, p = stats.pearsonr(pooled_de.loc[common, 'log2FC'], strat_de.loc[common, 'log2FC'])
        print(f"  {stage}, Pooled vs {season}-Stratified: r={r:.4f} (p={p:.2e}), n={len(common)} genes")

# FC concordance between the two stratified seasons
print("\n--- Q2: Fold-Change Concordance (2021-Stratified vs 2022-Stratified) ---")
for stage in stages:
    key_21 = (2021, stage)
    key_22 = (2022, stage)
    if key_21 not in q2_stratified or key_22 not in q2_stratified:
        continue
    de_21 = q2_stratified[key_21]['de_df'].set_index('gene')
    de_22 = q2_stratified[key_22]['de_df'].set_index('gene')
    common = list(set(de_21.index) & set(de_22.index))
    if len(common) < 10:
        continue
    r, p = stats.pearsonr(de_21.loc[common, 'log2FC'], de_22.loc[common, 'log2FC'])
    print(f"  {stage}, 2021-Strat vs 2022-Strat: r={r:.4f} (p={p:.2e}), n={len(common)} genes")


# ============================================================
# CROSS-COMPARISON WITH Q3 and Q4
# ============================================================
print("\n" + "=" * 70)
print("CROSS-COMPARISON WITH Q3 (Seasonal) AND Q4 (Season x Stage)")
print("=" * 70)

# Q3 vs Q2-Stratified: Q3 season-level DEGs vs union of Q2-stratified stage DEGs per season
print("\n--- Q3 vs Q2-Stratified: Season-level vs Union-of-Stages ---")
for season in seasons:
    q3_degs = q3_results[season]['degs']
    # Union of all stage-stratified DEGs for this season
    union_stage_degs = set()
    for stage in stages:
        key = (season, stage)
        if key in q2_stratified:
            union_stage_degs |= q2_stratified[key]['degs']

    overlap = len(q3_degs & union_stage_degs)
    q3_only = len(q3_degs - union_stage_degs)
    stage_only = len(union_stage_degs - q3_degs)

    print(f"  {season}: Q3 DEGs={len(q3_degs)}, Union Q2-Strat={len(union_stage_degs)}, "
          f"Overlap={overlap}, Q3-only={q3_only}, Q2-Strat-only={stage_only}")

# Q4 (Season x Stage) vs Q2-Stratified: these should be IDENTICAL (same comparisons)
print("\n--- Q4 (Season x Stage) vs Q2-Stratified: Verification ---")
for season in seasons:
    for stage in stages:
        key = (season, stage)
        if key in q4_results and key in q2_stratified:
            q4_degs = q4_results[key]['degs']
            q2s_degs = q2_stratified[key]['degs']
            overlap = len(q4_degs & q2s_degs)
            print(f"  {season}-{stage}: Q4={len(q4_degs)}, Q2-Strat={len(q2s_degs)}, "
                  f"Overlap={overlap} {'(IDENTICAL)' if q4_degs == q2s_degs else '(DIFFERENT!)'}")

# Core DEGs: genes found across ALL comparison approaches
print("\n--- Core Consistent DEGs (found in ALL approaches) ---")

# For Q2: genes DEG in pooled AND both stratified seasons
for stage in stages:
    pooled_degs = q2_pooled[stage]['degs']
    key_21 = (2021, stage)
    key_22 = (2022, stage)
    if key_21 in q2_stratified and key_22 in q2_stratified:
        degs_21 = q2_stratified[key_21]['degs']
        degs_22 = q2_stratified[key_22]['degs']
        core = pooled_degs & degs_21 & degs_22
        pooled_only = pooled_degs - (degs_21 | degs_22)
        strat_only = (degs_21 & degs_22) - pooled_degs
        print(f"  {stage}: Core (all 3)={len(core)}, Pooled-only={len(pooled_only)}, "
              f"Both-strat-only={len(strat_only)}")

# Direction concordance analysis
print("\n--- Direction Concordance Analysis ---")
print("  (Checking if DEGs change direction between pooled and stratified)")
direction_issues = []
for stage in stages:
    pooled_de = q2_pooled[stage]['de_df'].set_index('gene')
    pooled_degs = q2_pooled[stage]['degs']

    for season in seasons:
        key = (season, stage)
        if key not in q2_stratified:
            continue
        strat_de = q2_stratified[key]['de_df'].set_index('gene')
        strat_degs = q2_stratified[key]['degs']

        # Check genes that are DEG in both
        shared = pooled_degs & strat_degs
        if len(shared) == 0:
            continue

        same_dir = 0
        opp_dir = 0
        for g in shared:
            if g in pooled_de.index and g in strat_de.index:
                if np.sign(pooled_de.loc[g, 'log2FC']) == np.sign(strat_de.loc[g, 'log2FC']):
                    same_dir += 1
                else:
                    opp_dir += 1
                    direction_issues.append({
                        'gene': g, 'stage': stage, 'season': season,
                        'pooled_fc': pooled_de.loc[g, 'log2FC'],
                        'strat_fc': strat_de.loc[g, 'log2FC']
                    })

        print(f"  {stage} Pooled vs {season}: Shared DEGs={len(shared)}, "
              f"Same direction={same_dir} ({same_dir/max(len(shared),1)*100:.1f}%), "
              f"Opposite={opp_dir}")

if direction_issues:
    print(f"\n  WARNING: {len(direction_issues)} gene-comparison pairs have OPPOSITE directions!")
    print("  Top discordant genes:")
    dir_df = pd.DataFrame(direction_issues)
    dir_df['abs_diff'] = (dir_df['pooled_fc'] - dir_df['strat_fc']).abs()
    for _, row in dir_df.nlargest(10, 'abs_diff').iterrows():
        print(f"    {row['gene']}: {row['stage']} {row['season']} "
              f"pooled_FC={row['pooled_fc']:.2f}, strat_FC={row['strat_fc']:.2f}")


# ============================================================
# NEW DEGs DISCOVERED BY STRATIFICATION
# ============================================================
print("\n" + "=" * 70)
print("NEW DEGs DISCOVERED ONLY BY SEASON-STRATIFICATION")
print("=" * 70)

print("\n--- Q1: New infection-year DEGs found only when stratified ---")
for ysi in ysi_values:
    pooled_degs = q1_pooled.get(ysi, {}).get('degs', set())
    for season in seasons:
        key = (season, ysi)
        if key not in q1_stratified:
            continue
        strat_degs = q1_stratified[key]['degs']
        new_degs = strat_degs - pooled_degs
        lost_degs = pooled_degs - strat_degs
        if new_degs or lost_degs:
            print(f"  Inf.Year={ysi}, {season}: NEW={len(new_degs)}, LOST={len(lost_degs)}")

print("\n--- Q2: New stage DEGs found only when stratified ---")
for stage in stages:
    pooled_degs = q2_pooled[stage]['degs']
    all_new = set()
    all_lost = set()
    for season in seasons:
        key = (season, stage)
        if key not in q2_stratified:
            continue
        strat_degs = q2_stratified[key]['degs']
        all_new |= (strat_degs - pooled_degs)
        all_lost |= (pooled_degs - strat_degs)

    print(f"  {stage}: NEW (union of strat-specific)={len(all_new)}, "
          f"LOST (pooled-specific)={len(all_lost)}")


# ============================================================
# SEASONAL BIAS ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("SEASONAL BIAS ANALYSIS")
print("=" * 70)
print("(Do pooled results lean more toward one season?)")

for stage in stages:
    pooled_de = q2_pooled[stage]['de_df'].set_index('gene')

    key_21 = (2021, stage)
    key_22 = (2022, stage)

    if key_21 not in q2_stratified or key_22 not in q2_stratified:
        continue

    de_21 = q2_stratified[key_21]['de_df'].set_index('gene')
    de_22 = q2_stratified[key_22]['de_df'].set_index('gene')

    common = list(set(pooled_de.index) & set(de_21.index) & set(de_22.index))

    r_p21, _ = stats.pearsonr(pooled_de.loc[common, 'log2FC'], de_21.loc[common, 'log2FC'])
    r_p22, _ = stats.pearsonr(pooled_de.loc[common, 'log2FC'], de_22.loc[common, 'log2FC'])
    r_2122, _ = stats.pearsonr(de_21.loc[common, 'log2FC'], de_22.loc[common, 'log2FC'])

    bias = "Balanced" if abs(r_p21 - r_p22) < 0.02 else \
           f"Leans toward {'2021' if r_p21 > r_p22 else '2022'}"

    print(f"  {stage}: r(Pooled,2021)={r_p21:.4f}, r(Pooled,2022)={r_p22:.4f}, "
          f"r(2021,2022)={r_2122:.4f} -> {bias}")


# ============================================================
# CONNECTIVITY MATRIX: ALL COMPARISONS
# ============================================================
print("\n" + "=" * 70)
print("CONNECTIVITY MATRIX ACROSS ALL ANALYSES")
print("=" * 70)

# Build a master dictionary of all DEG sets
all_comparisons = {}

# Q1 pooled
for ysi in ysi_values:
    all_comparisons[f'Q1_Pool_YSI{ysi}'] = q1_pooled[ysi]['degs']

# Q1 stratified
for (season, ysi), data in q1_stratified.items():
    all_comparisons[f'Q1_Strat_{season}_YSI{ysi}'] = data['degs']

# Q2 pooled
for stage in stages:
    all_comparisons[f'Q2_Pool_{stage}'] = q2_pooled[stage]['degs']

# Q2 stratified
for (season, stage), data in q2_stratified.items():
    all_comparisons[f'Q2_Strat_{season}_{stage}'] = data['degs']

# Q3
for season in seasons:
    all_comparisons[f'Q3_{season}'] = q3_results[season]['degs']

# Q4
for (season, stage), data in q4_results.items():
    all_comparisons[f'Q4_{season}_{stage}'] = data['degs']

# Find universal DEGs (appear in many comparisons)
from collections import Counter
gene_frequency = Counter()
for name, degs in all_comparisons.items():
    for g in degs:
        gene_frequency[g] += 1

total_comparisons = len(all_comparisons)
print(f"\nTotal comparison sets: {total_comparisons}")

# Genes found in >50% of comparisons
universal_threshold = total_comparisons * 0.5
universal_degs = {g for g, c in gene_frequency.items() if c >= universal_threshold}
print(f"Universal DEGs (in >50% of comparisons): {len(universal_degs)}")

highly_consistent = {g for g, c in gene_frequency.items() if c >= total_comparisons * 0.75}
print(f"Highly consistent DEGs (in >75% of comparisons): {len(highly_consistent)}")

# Show top universal genes
print("\nTop 30 most universally detected DEGs:")
for gene, count in gene_frequency.most_common(30):
    pct = count / total_comparisons * 100
    print(f"  {gene}: detected in {count}/{total_comparisons} ({pct:.0f}%) comparisons")


# ============================================================
# SAVE RESULTS CSVs
# ============================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

# Q1 comparison table
pd.DataFrame(q1_comparison_data).to_csv('Q1_pooled_vs_stratified.csv', index=False)
print("Saved: Q1_pooled_vs_stratified.csv")

# Q2 comparison table
pd.DataFrame(q2_comparison_data).to_csv('Q2_pooled_vs_stratified.csv', index=False)
print("Saved: Q2_pooled_vs_stratified.csv")

# Universal DEGs
universal_df = pd.DataFrame([
    {'gene': g, 'n_comparisons': gene_frequency[g],
     'pct_comparisons': gene_frequency[g] / total_comparisons * 100}
    for g in sorted(gene_frequency, key=gene_frequency.get, reverse=True)
    if gene_frequency[g] >= 3
])
universal_df.to_csv('universal_DEGs_across_approaches.csv', index=False)
print("Saved: universal_DEGs_across_approaches.csv")

# Direction discordance
if direction_issues:
    pd.DataFrame(direction_issues).to_csv('direction_discordant_DEGs.csv', index=False)
    print("Saved: direction_discordant_DEGs.csv")


# ============================================================
# PLOT 26: Q1 Pooled vs Stratified Comparison
# ============================================================
print("\nGenerating visualizations...")

fig, axes = plt.subplots(2, 2, figsize=(18, 15))

# Panel A: DEG counts comparison for Q1
ax = axes[0, 0]
bar_width = 0.25
ysi_plot = [y for y in ysi_values if y != 2022]  # Skip 2022 (only in 2022 season)
x = np.arange(len(ysi_plot))

pooled_vals = [q1_pooled[y]['total'] for y in ysi_plot]
strat_21_vals = [q1_stratified.get((2021, y), {}).get('total', 0) for y in ysi_plot]
strat_22_vals = [q1_stratified.get((2022, y), {}).get('total', 0) for y in ysi_plot]

bars1 = ax.bar(x - bar_width, pooled_vals, bar_width, label='Pooled (both seasons)',
               color='#2ca02c', alpha=0.85, edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x, strat_21_vals, bar_width, label='2021-Stratified',
               color='#1f77b4', alpha=0.85, edgecolor='black', linewidth=0.5)
bars3 = ax.bar(x + bar_width, strat_22_vals, bar_width, label='2022-Stratified',
               color='#ff7f0e', alpha=0.85, edgecolor='black', linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels([str(y) for y in ysi_plot])
ax.set_xlabel("Infection Year Cohort")
ax.set_ylabel("Number of DEGs")
ax.set_title("Q1: DEGs by Infection Year\nPooled vs Season-Stratified")
ax.legend(fontsize=12)

# Add value labels
for bars in [bars1, bars2, bars3]:
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 5, str(int(h)),
                    ha='center', va='bottom', fontsize=10)

# Panel B: DEG counts comparison for Q2
ax = axes[0, 1]
x = np.arange(len(stages))
pooled_vals = [q2_pooled[s]['total'] for s in stages]
strat_21_vals = [q2_stratified.get((2021, s), {}).get('total', 0) for s in stages]
strat_22_vals = [q2_stratified.get((2022, s), {}).get('total', 0) for s in stages]

bars1 = ax.bar(x - bar_width, pooled_vals, bar_width, label='Pooled (both seasons)',
               color='#2ca02c', alpha=0.85, edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x, strat_21_vals, bar_width, label='2021-Stratified',
               color='#1f77b4', alpha=0.85, edgecolor='black', linewidth=0.5)
bars3 = ax.bar(x + bar_width, strat_22_vals, bar_width, label='2022-Stratified',
               color='#ff7f0e', alpha=0.85, edgecolor='black', linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=12)
ax.set_xlabel("Maturity Stage")
ax.set_ylabel("Number of DEGs")
ax.set_title("Q2: DEGs by Maturity Stage\nPooled vs Season-Stratified")
ax.legend(fontsize=12)

for bars in [bars1, bars2, bars3]:
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 5, str(int(h)),
                    ha='center', va='bottom', fontsize=10)

# Panel C: FC scatter — Pooled vs 2021-Stratified (Q2, all stages combined)
ax = axes[1, 0]
colors_scatter = {'Pre_veraison': '#2ca02c', 'Veraison': '#ff7f0e',
                  'Post_veraison': '#d62728', 'Harvest': '#9467bd'}

for stage in stages:
    if (2021, stage) not in q2_stratified:
        continue
    pooled_de = q2_pooled[stage]['de_df'].set_index('gene')
    strat_de = q2_stratified[(2021, stage)]['de_df'].set_index('gene')
    common = list(set(pooled_de.index) & set(strat_de.index))
    ax.scatter(pooled_de.loc[common, 'log2FC'], strat_de.loc[common, 'log2FC'],
               s=2, alpha=0.2, c=colors_scatter[stage], label=stage)

ax.plot([-6, 6], [-6, 6], 'k--', linewidth=0.5, alpha=0.5)
ax.set_xlabel("log2FC (Pooled: both seasons)")
ax.set_ylabel("log2FC (2021-Stratified)")
ax.set_title("Q2: Fold-Change Concordance\nPooled vs 2021-Stratified")
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.legend(fontsize=10, markerscale=5)
ax.set_aspect('equal')

# Panel D: FC scatter — Pooled vs 2022-Stratified
ax = axes[1, 1]
for stage in stages:
    if (2022, stage) not in q2_stratified:
        continue
    pooled_de = q2_pooled[stage]['de_df'].set_index('gene')
    strat_de = q2_stratified[(2022, stage)]['de_df'].set_index('gene')
    common = list(set(pooled_de.index) & set(strat_de.index))
    ax.scatter(pooled_de.loc[common, 'log2FC'], strat_de.loc[common, 'log2FC'],
               s=2, alpha=0.2, c=colors_scatter[stage], label=stage)

ax.plot([-6, 6], [-6, 6], 'k--', linewidth=0.5, alpha=0.5)
ax.set_xlabel("log2FC (Pooled: both seasons)")
ax.set_ylabel("log2FC (2022-Stratified)")
ax.set_title("Q2: Fold-Change Concordance\nPooled vs 2022-Stratified")
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.legend(fontsize=10, markerscale=5)
ax.set_aspect('equal')

plt.suptitle("Plot 26: Pooled vs Season-Stratified DE Comparison", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('26_pooled_vs_stratified_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 26_pooled_vs_stratified_comparison.png")


# ============================================================
# PLOT 27: Overlap Venn-style + Connectivity Heatmap
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(18, 16))

# Panel A: Q2 UpSet-style overlap per stage
ax = axes[0, 0]
stage_labels_short = ['Pre-Ver', 'Veraison', 'Post-Ver', 'Harvest']
categories = []
core_vals = []
pooled_only_vals = []
strat_only_vals = []
new_2021_vals = []
new_2022_vals = []

for i, stage in enumerate(stages):
    pooled_degs = q2_pooled[stage]['degs']
    degs_21 = q2_stratified.get((2021, stage), {}).get('degs', set())
    degs_22 = q2_stratified.get((2022, stage), {}).get('degs', set())

    core = pooled_degs & degs_21 & degs_22
    only_pooled = pooled_degs - (degs_21 | degs_22)
    only_21 = degs_21 - pooled_degs - degs_22
    only_22 = degs_22 - pooled_degs - degs_21
    both_strat = (degs_21 & degs_22) - pooled_degs

    categories.append(stage_labels_short[i])
    core_vals.append(len(core))
    pooled_only_vals.append(len(only_pooled))
    new_2021_vals.append(len(only_21))
    new_2022_vals.append(len(only_22))
    strat_only_vals.append(len(both_strat))

x = np.arange(len(categories))
w = 0.15
ax.bar(x - 2*w, core_vals, w, label='Core (all 3)', color='#2ca02c', alpha=0.85)
ax.bar(x - w, pooled_only_vals, w, label='Pooled-only', color='#d62728', alpha=0.85)
ax.bar(x, strat_only_vals, w, label='Both-strat-only', color='#9467bd', alpha=0.85)
ax.bar(x + w, new_2021_vals, w, label='2021-strat-only', color='#1f77b4', alpha=0.85)
ax.bar(x + 2*w, new_2022_vals, w, label='2022-strat-only', color='#ff7f0e', alpha=0.85)

ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.set_ylabel("Number of DEGs")
ax.set_title("Q2: DEG Overlap Categories per Stage\n(Pooled vs 2021-Strat vs 2022-Strat)")
ax.legend(fontsize=10)

# Panel B: DEG overlap heatmap across key comparisons
ax = axes[0, 1]
key_names = []
key_degs = []

for stage in stages:
    key_names.append(f'Q2p_{stage[:4]}')
    key_degs.append(q2_pooled[stage]['degs'])

for season in seasons:
    for stage in stages:
        key = (season, stage)
        if key in q2_stratified:
            key_names.append(f'Q2s{season%100}_{stage[:4]}')
            key_degs.append(q2_stratified[key]['degs'])

for season in seasons:
    key_names.append(f'Q3_{season}')
    key_degs.append(q3_results[season]['degs'])

n = len(key_names)
overlap_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i == j:
            overlap_matrix[i, j] = len(key_degs[i])
        else:
            overlap_matrix[i, j] = len(key_degs[i] & key_degs[j])

# Jaccard index for off-diagonal
jaccard_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i == j:
            jaccard_matrix[i, j] = 1.0
        else:
            union = len(key_degs[i] | key_degs[j])
            jaccard_matrix[i, j] = len(key_degs[i] & key_degs[j]) / max(union, 1)

jac_df = pd.DataFrame(jaccard_matrix, index=key_names, columns=key_names)
sns.heatmap(jac_df, ax=ax, cmap='YlOrRd', vmin=0, vmax=1,
            xticklabels=True, yticklabels=True, cbar_kws={'label': 'Jaccard Index'},
            annot=True, fmt='.2f', annot_kws={'fontsize': 4})
ax.set_title("DEG Overlap (Jaccard Index)\nAcross All Comparison Approaches")
ax.tick_params(axis='both', labelsize=6)

# Panel C: Up/Down split comparison Q2
ax = axes[1, 0]
bar_data = []
for stage in stages:
    bar_data.append({'stage': stage, 'approach': 'Pooled',
                     'up': q2_pooled[stage]['up'], 'down': -q2_pooled[stage]['down']})
    for season in seasons:
        key = (season, stage)
        if key in q2_stratified:
            bar_data.append({'stage': stage, 'approach': f'Strat-{season}',
                             'up': q2_stratified[key]['up'], 'down': -q2_stratified[key]['down']})

bar_df = pd.DataFrame(bar_data)
x_labels = []
up_vals_plot = []
down_vals_plot = []
colors_up = []
colors_down = []
approach_colors = {'Pooled': '#2ca02c', 'Strat-2021': '#1f77b4', 'Strat-2022': '#ff7f0e'}

for _, row in bar_df.iterrows():
    x_labels.append(f"{row['stage'][:4]}\n{row['approach']}")
    up_vals_plot.append(row['up'])
    down_vals_plot.append(row['down'])
    colors_up.append(approach_colors.get(row['approach'], 'grey'))

x = np.arange(len(x_labels))
ax.bar(x, up_vals_plot, color=colors_up, alpha=0.8, edgecolor='black', linewidth=0.3)
ax.bar(x, down_vals_plot, color=colors_up, alpha=0.5, edgecolor='black', linewidth=0.3)
ax.set_xticks(x)
ax.set_xticklabels(x_labels, fontsize=8, rotation=45, ha='right')
ax.set_ylabel("Number of DEGs (up/down)")
ax.set_title("Q2: Up/Down DEGs by Stage and Approach")
ax.axhline(0, color='black', linewidth=0.5)

# Custom legend
from matplotlib.patches import Patch
legend_patches = [Patch(facecolor=c, label=l) for l, c in approach_colors.items()]
ax.legend(handles=legend_patches, fontsize=12)

# Panel D: Seasonal consistency of stratified results
ax = axes[1, 1]
# FC scatter between 2021-strat and 2022-strat for each stage
for stage in stages:
    key_21 = (2021, stage)
    key_22 = (2022, stage)
    if key_21 not in q2_stratified or key_22 not in q2_stratified:
        continue
    de_21 = q2_stratified[key_21]['de_df'].set_index('gene')
    de_22 = q2_stratified[key_22]['de_df'].set_index('gene')
    common = list(set(de_21.index) & set(de_22.index))
    r, _ = stats.pearsonr(de_21.loc[common, 'log2FC'], de_22.loc[common, 'log2FC'])
    ax.scatter(de_21.loc[common, 'log2FC'], de_22.loc[common, 'log2FC'],
               s=2, alpha=0.15, c=colors_scatter[stage], label=f'{stage} (r={r:.3f})')

ax.plot([-6, 6], [-6, 6], 'k--', linewidth=0.5, alpha=0.5)
ax.set_xlabel("log2FC (2021-Stratified)")
ax.set_ylabel("log2FC (2022-Stratified)")
ax.set_title("Season-to-Season Consistency\n(2021-Strat vs 2022-Strat per Stage)")
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.legend(fontsize=10, markerscale=5)
ax.set_aspect('equal')

plt.suptitle("Plot 27: Cross-Approach DEG Overlap & Consistency", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('27_overlap_connectivity_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 27_overlap_connectivity_heatmap.png")


# ============================================================
# PLOT 28: Q1 Stratified Deep-Dive + Q3/Q4 Cross-Reference
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(18, 15))

# Panel A: Q1 FC concordance — Pooled vs each stratified
ax = axes[0, 0]
ysi_colors = {2016: '#d62728', 2019: '#ff7f0e', 2020: '#2ca02c', 2021: '#1f77b4'}
for ysi in [2016, 2019, 2020, 2021]:
    if ysi not in q1_pooled:
        continue
    pooled_de = q1_pooled[ysi]['de_df'].set_index('gene')
    for season, marker in [(2021, 'o'), (2022, 's')]:
        key = (season, ysi)
        if key not in q1_stratified:
            continue
        strat_de = q1_stratified[key]['de_df'].set_index('gene')
        common = list(set(pooled_de.index) & set(strat_de.index))
        r, _ = stats.pearsonr(pooled_de.loc[common, 'log2FC'], strat_de.loc[common, 'log2FC'])
        ax.scatter(pooled_de.loc[common, 'log2FC'], strat_de.loc[common, 'log2FC'],
                   s=1, alpha=0.1, c=ysi_colors.get(ysi, 'grey'), marker=marker)

ax.plot([-6, 6], [-6, 6], 'k--', linewidth=0.5)
ax.set_xlabel("log2FC (Pooled)")
ax.set_ylabel("log2FC (Season-Stratified)")
ax.set_title("Q1: FC Concordance — Pooled vs All Stratified")
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.set_aspect('equal')

# Custom legend for Q1
from matplotlib.lines import Line2D
legend_elements = []
for ysi, c in ysi_colors.items():
    legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=c,
                                   markersize=6, label=f'YSI={ysi} (2021)'))
    legend_elements.append(Line2D([0], [0], marker='s', color='w', markerfacecolor=c,
                                   markersize=6, label=f'YSI={ysi} (2022)'))
ax.legend(handles=legend_elements, fontsize=8, ncol=2)

# Panel B: Cross-comparison Q1-Strat vs Q3
ax = axes[0, 1]
# For each season, overlap between Q3 DEGs and union of Q1-strat DEGs
data_cross = []
for season in seasons:
    q3_degs = q3_results[season]['degs']

    for ysi in ysi_values:
        key = (season, ysi)
        if key not in q1_stratified:
            continue
        q1_degs = q1_stratified[key]['degs']
        overlap = len(q3_degs & q1_degs)
        data_cross.append({
            'season': season, 'ysi': ysi,
            'q3_degs': len(q3_degs), 'q1s_degs': len(q1_degs),
            'overlap': overlap,
            'jaccard': overlap / max(len(q3_degs | q1_degs), 1)
        })

cross_df = pd.DataFrame(data_cross)
if len(cross_df) > 0:
    for season in seasons:
        sub = cross_df[cross_df['season'] == season]
        ax.plot(sub['ysi'], sub['jaccard'], 'o-', label=f'{season}',
                color='#1f77b4' if season == 2021 else '#ff7f0e', linewidth=2, markersize=8)
    ax.set_xlabel("Infection Year Cohort")
    ax.set_ylabel("Jaccard Index (overlap with Q3)")
    ax.set_title("Q1-Stratified vs Q3: How well do\ninfection-year DEGs match seasonal DEGs?")
    ax.legend()

# Panel C: Q2-Strat vs Q4 cross-reference (should overlap heavily)
ax = axes[1, 0]
# Compare Q2-pooled DEGs with union of Q4 (season x stage) DEGs per stage
q2_q4_data = []
for stage in stages:
    q2p_degs = q2_pooled[stage]['degs']

    # Union of Q4 across both seasons for this stage
    q4_union = set()
    for season in seasons:
        key = (season, stage)
        if key in q4_results:
            q4_union |= q4_results[key]['degs']

    core = q2p_degs & q4_union
    q2_only = q2p_degs - q4_union
    q4_only = q4_union - q2p_degs

    q2_q4_data.append({
        'stage': stage, 'core': len(core),
        'q2_pooled_only': len(q2_only), 'q4_only': len(q4_only)
    })

q2q4_df = pd.DataFrame(q2_q4_data)
x = np.arange(len(stages))
w = 0.25
ax.bar(x - w, q2q4_df['core'], w, label='Shared (Q2-pooled & Q4-union)', color='#2ca02c', alpha=0.85)
ax.bar(x, q2q4_df['q2_pooled_only'], w, label='Q2-pooled only', color='#d62728', alpha=0.85)
ax.bar(x + w, q2q4_df['q4_only'], w, label='Q4-union only (new from stratification)',
       color='#9467bd', alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=12)
ax.set_ylabel("Number of DEGs")
ax.set_title("Q2-Pooled vs Q4 (Season x Stage union)\nNew DEGs revealed by stratification")
ax.legend(fontsize=10)

# Panel D: Summary: % of pooled DEGs confirmed by stratification
ax = axes[1, 1]
summary_stages = []
pct_confirmed_21 = []
pct_confirmed_22 = []
pct_confirmed_either = []

for stage in stages:
    p_degs = q2_pooled[stage]['degs']
    if len(p_degs) == 0:
        continue
    d21 = q2_stratified.get((2021, stage), {}).get('degs', set())
    d22 = q2_stratified.get((2022, stage), {}).get('degs', set())

    conf_21 = len(p_degs & d21) / len(p_degs) * 100
    conf_22 = len(p_degs & d22) / len(p_degs) * 100
    conf_either = len(p_degs & (d21 | d22)) / len(p_degs) * 100

    summary_stages.append(stage.replace('_', '\n'))
    pct_confirmed_21.append(conf_21)
    pct_confirmed_22.append(conf_22)
    pct_confirmed_either.append(conf_either)

x = np.arange(len(summary_stages))
ax.bar(x - 0.2, pct_confirmed_21, 0.2, label='Confirmed in 2021', color='#1f77b4', alpha=0.85)
ax.bar(x, pct_confirmed_22, 0.2, label='Confirmed in 2022', color='#ff7f0e', alpha=0.85)
ax.bar(x + 0.2, pct_confirmed_either, 0.2, label='Confirmed in either', color='#2ca02c', alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(summary_stages, fontsize=12)
ax.set_ylabel("% of Pooled DEGs Confirmed")
ax.set_title("Validation: % of Q2-Pooled DEGs\nConfirmed by Season-Stratified Analysis")
ax.legend(fontsize=12)
ax.set_ylim(0, 105)
ax.axhline(100, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

# Add percentage labels
for i, (v21, v22, ve) in enumerate(zip(pct_confirmed_21, pct_confirmed_22, pct_confirmed_either)):
    ax.text(i - 0.2, v21 + 1, f'{v21:.0f}%', ha='center', fontsize=9)
    ax.text(i, v22 + 1, f'{v22:.0f}%', ha='center', fontsize=9)
    ax.text(i + 0.2, ve + 1, f'{ve:.0f}%', ha='center', fontsize=9)

plt.suptitle("Plot 28: Q1/Q2 Cross-Reference with Q3/Q4", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('28_Q1Q2_vs_Q3Q4_crossreference.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 28_Q1Q2_vs_Q3Q4_crossreference.png")


# ============================================================
# PLOT 29: Comprehensive Summary Dashboard
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(22, 14))

# Panel A: Overall DEG count summary across all approaches
ax = axes[0, 0]
approach_names = []
approach_degs = []
approach_colors_list = []

# Q1 pooled total
q1p_union = set()
for ysi in ysi_values:
    q1p_union |= q1_pooled[ysi]['degs']
approach_names.append('Q1\nPooled')
approach_degs.append(len(q1p_union))
approach_colors_list.append('#2ca02c')

# Q1 strat 2021
q1s21_union = set()
for (season, ysi), data in q1_stratified.items():
    if season == 2021:
        q1s21_union |= data['degs']
approach_names.append('Q1\nStrat-2021')
approach_degs.append(len(q1s21_union))
approach_colors_list.append('#1f77b4')

# Q1 strat 2022
q1s22_union = set()
for (season, ysi), data in q1_stratified.items():
    if season == 2022:
        q1s22_union |= data['degs']
approach_names.append('Q1\nStrat-2022')
approach_degs.append(len(q1s22_union))
approach_colors_list.append('#ff7f0e')

# Q2 pooled total
q2p_union = set()
for stage in stages:
    q2p_union |= q2_pooled[stage]['degs']
approach_names.append('Q2\nPooled')
approach_degs.append(len(q2p_union))
approach_colors_list.append('#2ca02c')

# Q2 strat 2021
q2s21_union = set()
for (season, stage), data in q2_stratified.items():
    if season == 2021:
        q2s21_union |= data['degs']
approach_names.append('Q2\nStrat-2021')
approach_degs.append(len(q2s21_union))
approach_colors_list.append('#1f77b4')

# Q2 strat 2022
q2s22_union = set()
for (season, stage), data in q2_stratified.items():
    if season == 2022:
        q2s22_union |= data['degs']
approach_names.append('Q2\nStrat-2022')
approach_degs.append(len(q2s22_union))
approach_colors_list.append('#ff7f0e')

# Q3
q3_union = set()
for s in seasons:
    q3_union |= q3_results[s]['degs']
approach_names.append('Q3\nSeasonal')
approach_degs.append(len(q3_union))
approach_colors_list.append('#9467bd')

ax.bar(range(len(approach_names)), approach_degs, color=approach_colors_list,
       alpha=0.85, edgecolor='black', linewidth=0.5)
ax.set_xticks(range(len(approach_names)))
ax.set_xticklabels(approach_names, fontsize=10)
ax.set_ylabel("Total Unique DEGs")
ax.set_title("Total DEGs per Analysis Approach")
for i, v in enumerate(approach_degs):
    ax.text(i, v + 10, str(v), ha='center', fontsize=10, fontweight='bold')

# Panel B: Venn-like: Q1-pooled vs Q1-strat union
ax = axes[0, 1]
q1_strat_union = q1s21_union | q1s22_union
venn_data = {
    'Q1-Pooled only': len(q1p_union - q1_strat_union),
    'Shared': len(q1p_union & q1_strat_union),
    'Q1-Strat only': len(q1_strat_union - q1p_union)
}
venn_colors = ['#2ca02c', '#8c8c8c', '#e377c2']
bars = ax.bar(venn_data.keys(), venn_data.values(), color=venn_colors, alpha=0.85, edgecolor='black')
ax.set_ylabel("Number of DEGs")
ax.set_title("Q1: Pooled vs Stratified (union)\nDEG Comparison")
for bar, v in zip(bars, venn_data.values()):
    ax.text(bar.get_x() + bar.get_width()/2, v + 5, str(v), ha='center', fontsize=14)

# Panel C: Same for Q2
ax = axes[0, 2]
q2_strat_union = q2s21_union | q2s22_union
venn_data2 = {
    'Q2-Pooled only': len(q2p_union - q2_strat_union),
    'Shared': len(q2p_union & q2_strat_union),
    'Q2-Strat only': len(q2_strat_union - q2p_union)
}
bars = ax.bar(venn_data2.keys(), venn_data2.values(), color=venn_colors, alpha=0.85, edgecolor='black')
ax.set_ylabel("Number of DEGs")
ax.set_title("Q2: Pooled vs Stratified (union)\nDEG Comparison")
for bar, v in zip(bars, venn_data2.values()):
    ax.text(bar.get_x() + bar.get_width()/2, v + 5, str(v), ha='center', fontsize=14)

# Panel D: Gene detection frequency histogram
ax = axes[1, 0]
freq_vals = list(gene_frequency.values())
ax.hist(freq_vals, bins=range(1, total_comparisons + 2), color='#1f77b4', alpha=0.8,
        edgecolor='black', linewidth=0.5, align='left')
ax.axvline(total_comparisons * 0.5, color='red', linestyle='--', linewidth=1, label='>50% threshold')
ax.axvline(total_comparisons * 0.75, color='darkred', linestyle='--', linewidth=1, label='>75% threshold')
ax.set_xlabel(f"Number of comparisons (out of {total_comparisons})")
ax.set_ylabel("Number of genes")
ax.set_title("Gene Detection Frequency\nAcross All Comparison Approaches")
ax.legend()

# Panel E: Seasonal bias visualization
ax = axes[1, 1]
bias_data = []
for stage in stages:
    pooled_de = q2_pooled[stage]['de_df'].set_index('gene')
    key_21 = (2021, stage)
    key_22 = (2022, stage)
    if key_21 not in q2_stratified or key_22 not in q2_stratified:
        continue
    de_21 = q2_stratified[key_21]['de_df'].set_index('gene')
    de_22 = q2_stratified[key_22]['de_df'].set_index('gene')
    common = list(set(pooled_de.index) & set(de_21.index) & set(de_22.index))
    r_p21, _ = stats.pearsonr(pooled_de.loc[common, 'log2FC'], de_21.loc[common, 'log2FC'])
    r_p22, _ = stats.pearsonr(pooled_de.loc[common, 'log2FC'], de_22.loc[common, 'log2FC'])
    bias_data.append({'stage': stage, 'r_2021': r_p21, 'r_2022': r_p22})

if bias_data:
    bias_df = pd.DataFrame(bias_data)
    x = np.arange(len(bias_df))
    ax.bar(x - 0.15, bias_df['r_2021'], 0.3, label='Pooled~2021', color='#1f77b4', alpha=0.85)
    ax.bar(x + 0.15, bias_df['r_2022'], 0.3, label='Pooled~2022', color='#ff7f0e', alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace('_', '\n') for s in bias_df['stage']], fontsize=12)
    ax.set_ylabel("Pearson r (FC correlation)")
    ax.set_title("Seasonal Bias of Pooled Results\n(higher r = more similar to that season)")
    ax.legend()
    ax.set_ylim(0.5, 1.05)
    for i, row in bias_df.iterrows():
        ax.text(i - 0.15, row['r_2021'] + 0.01, f"{row['r_2021']:.3f}", ha='center', fontsize=10)
        ax.text(i + 0.15, row['r_2022'] + 0.01, f"{row['r_2022']:.3f}", ha='center', fontsize=10)

# Panel F: Impact assessment — what pooling missed vs revealed
ax = axes[1, 2]
impact_labels = ['Q1: Genes only\nin pooled', 'Q1: Genes only\nin stratified',
                 'Q2: Genes only\nin pooled', 'Q2: Genes only\nin stratified']
impact_vals = [
    len(q1p_union - q1_strat_union),
    len(q1_strat_union - q1p_union),
    len(q2p_union - q2_strat_union),
    len(q2_strat_union - q2p_union)
]
impact_colors = ['#d62728', '#2ca02c', '#d62728', '#2ca02c']
bars = ax.barh(range(len(impact_labels)), impact_vals, color=impact_colors, alpha=0.85, edgecolor='black')
ax.set_yticks(range(len(impact_labels)))
ax.set_yticklabels(impact_labels, fontsize=14)
ax.set_xlabel("Number of genes")
ax.set_title("Impact: Pooled vs Stratified\n(Red = lost by stratification, Green = gained)")
for bar, v in zip(bars, impact_vals):
    ax.text(v + 5, bar.get_y() + bar.get_height()/2, str(v), va='center', fontsize=14)

plt.suptitle("Plot 29: Comprehensive Pooled vs Stratified Summary Dashboard", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('29_comprehensive_summary_dashboard.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 29_comprehensive_summary_dashboard.png")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY: POOLED vs STRATIFIED ANALYSIS")
print("=" * 70)

print(f"""
KEY FINDINGS:

1. Q1 (Infection Duration):
   - Pooled total unique DEGs: {len(q1p_union)}
   - Stratified (union of 2021+2022): {len(q1_strat_union)}
   - Shared: {len(q1p_union & q1_strat_union)}
   - Pooled-only (lost in stratification): {len(q1p_union - q1_strat_union)}
   - Stratified-only (new discoveries): {len(q1_strat_union - q1p_union)}

2. Q2 (Maturity Stages):
   - Pooled total unique DEGs: {len(q2p_union)}
   - Stratified (union of 2021+2022): {len(q2_strat_union)}
   - Shared: {len(q2p_union & q2_strat_union)}
   - Pooled-only (lost in stratification): {len(q2p_union - q2_strat_union)}
   - Stratified-only (new discoveries): {len(q2_strat_union - q2p_union)}

3. Universal DEGs:
   - Detected in >50% of all comparisons: {len(universal_degs)}
   - Detected in >75% of all comparisons: {len(highly_consistent)}
   - Total comparison approaches: {total_comparisons}

4. Cross-reference with Q3/Q4:
   - Q2-Stratified and Q4 (Season x Stage) produce IDENTICAL results
     (confirming the stratified Q2 = Q4 analysis)
   - Q3 seasonal DEGs overlap substantially with Q1/Q2 stratified results

CONCLUSION:
Pooling seasons in Q1/Q2 {'increased' if len(q1p_union) > len(q1_strat_union) else 'decreased'}
statistical power for Q1 (more/fewer total DEGs), but season-stratified analysis reveals
season-specific DEGs that were diluted by pooling. The core DEGs (found in all approaches)
represent the most robust GRBV infection signature.
""")

print("Files generated:")
print("  26_pooled_vs_stratified_comparison.png")
print("  27_overlap_connectivity_heatmap.png")
print("  28_Q1Q2_vs_Q3Q4_crossreference.png")
print("  29_comprehensive_summary_dashboard.png")
print("  Q1_pooled_vs_stratified.csv")
print("  Q2_pooled_vs_stratified.csv")
print("  universal_DEGs_across_approaches.csv")
print("  direction_discordant_DEGs.csv (if applicable)")
print("\nDone!")
