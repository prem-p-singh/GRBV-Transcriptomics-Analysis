#!/usr/bin/env python3
"""
Gene Set Enrichment Analysis (GSEA) for GRBV Transcriptomics.

Complements the Fisher's exact ORA by using full ranked gene lists rather
than hard DE cutoffs.  Runs permutation-based GSEA across 8 season × stage
conditions and compares with existing Fisher ORA results.

OPTIMIZED VERSION: NumPy-vectorized for HPC execution.
  - Vectorized t-tests (no per-gene Python loop)
  - Batch ES computation across all gene sets per permutation
  - Pre-computed hit/miss masks as boolean matrices
  - ~50-100x faster than pure Python loops

Outputs:
  - Fig 37: GSEA enrichment plots (running ES + top pathways)
  - Fig 38: Fisher ORA vs GSEA concordance scatter
  - Fig 39: GSEA NES heatmap across season × stage conditions
  - gsea_results_all_conditions.csv
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
from collections import defaultdict
import textwrap
import time
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
print("GSEA ANALYSIS: LOADING DATA (VECTORIZED)")
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
# 2. LOAD ANNOTATIONS — GO + KEGG gene sets
# ============================================================
print("\nLoading annotations...")

db_path = 'org.Vvinifera.eg.db/inst/extdata/org.Vvinifera.eg.sqlite'
conn = sqlite3.connect(db_path)

go_df = pd.read_sql(
    "SELECT g.SYMBOL, go.GO, go.ONTOLOGY "
    "FROM go JOIN gene_info g ON go._id = g._id", conn)
kegg_df = pd.read_sql(
    "SELECT g.SYMBOL, k.KEGG "
    "FROM kegg k JOIN gene_info g ON k._id = g._id", conn)
conn.close()

# GO term names
from goatools.obo_parser import GODag
godag = GODag('go-basic.obo')

# KEGG pathway names
import urllib.request
try:
    url = "https://rest.kegg.jp/list/pathway/vvi"
    response = urllib.request.urlopen(url, timeout=30)
    kegg_names = {}
    for line in response.read().decode().strip().split('\n'):
        parts = line.split('\t')
        if len(parts) == 2:
            pid = parts[0].replace('path:', '')
            kegg_names[pid] = parts[1].replace(' - Vitis vinifera (wine grape)', '').strip()
    print(f"KEGG pathway names: {len(kegg_names)}")
except Exception:
    kegg_names = {}
    print("KEGG names fetch failed, using IDs")

# Build gene set dictionaries  { term_id : set(genes) }
# Only keep terms with 15-500 genes (standard GSEA filter)
expressed = set(counts_filtered.index)

go2genes = defaultdict(set)
for _, row in go_df.iterrows():
    gene, term, ont = row['SYMBOL'], row['GO'], row['ONTOLOGY']
    if gene in expressed and ont in ('BP', 'MF'):
        go2genes[term].add(gene)

kegg2genes = defaultdict(set)
for _, row in kegg_df.iterrows():
    gene, pathway = row['SYMBOL'], row['KEGG']
    if gene in expressed:
        kegg2genes[pathway].add(gene)

# Filter by size
MIN_SET = 15
MAX_SET = 500

gene_sets = {}
for term, genes in go2genes.items():
    if MIN_SET <= len(genes) <= MAX_SET:
        name = godag[term].name if term in godag else term
        gene_sets[f"GO:{name}"] = genes

for pathway, genes in kegg2genes.items():
    if MIN_SET <= len(genes) <= MAX_SET:
        name = kegg_names.get(pathway, pathway)
        gene_sets[f"KEGG:{name}"] = genes

print(f"Gene sets for GSEA: {len(gene_sets)} (GO-BP/MF + KEGG, {MIN_SET}-{MAX_SET} genes)")

# ============================================================
# 3. VECTORIZED GSEA ALGORITHM
# ============================================================

def rank_genes_vectorized(group1_samples, group2_samples, log_cpm_data):
    """Rank all genes by sign(log2FC) * -log10(p) — VECTORIZED.

    Uses scipy.stats.ttest_ind with axis parameter to compute all
    t-tests in a single call (~100x faster than per-gene loop).
    """
    g1 = log_cpm_data[group1_samples].values  # (n_genes, n1)
    g2 = log_cpm_data[group2_samples].values  # (n_genes, n2)
    genes = log_cpm_data.index.values

    # Vectorized t-test across all genes at once
    t_stats, p_vals = stats.ttest_ind(g1, g2, axis=1, equal_var=False)

    # Handle NaN p-values (constant expression)
    p_vals = np.where(np.isnan(p_vals), 1.0, p_vals)
    t_stats = np.where(np.isnan(t_stats), 0.0, t_stats)

    # log2FC = mean(group1) - mean(group2)
    log2fc = g1.mean(axis=1) - g2.mean(axis=1)

    # Rank metric: sign(log2FC) * -log10(p)
    rank_metric = np.sign(log2fc) * -np.log10(np.maximum(p_vals, 1e-300))

    df = pd.DataFrame({
        'gene': genes,
        'log2FC': log2fc,
        'rank_metric': rank_metric,
        'pvalue': p_vals
    })
    df = df.sort_values('rank_metric', ascending=False).reset_index(drop=True)
    return df


def compute_enrichment_score(ranked_genes, gene_set_members, weighted=True):
    """Compute running enrichment score for a single gene set."""
    N = len(ranked_genes)
    gene_list = ranked_genes['gene'].values
    metrics = ranked_genes['rank_metric'].values.astype(float)

    hits = np.isin(gene_list, list(gene_set_members))
    n_hits = hits.sum()
    if n_hits == 0:
        return 0.0, np.zeros(N)

    if weighted:
        hit_weights = np.abs(metrics) * hits
        nr = hit_weights.sum()
        if nr == 0:
            return 0.0, np.zeros(N)
        p_hit = np.cumsum(hit_weights) / nr
    else:
        p_hit = np.cumsum(hits) / n_hits

    p_miss = np.cumsum(~hits) / (N - n_hits)
    running_es = p_hit - p_miss

    # ES is the maximum deviation from 0
    es_max = running_es.max()
    es_min = running_es.min()
    es = es_max if abs(es_max) >= abs(es_min) else es_min

    return es, running_es


def compute_es_batch_vectorized(rank_metrics_sorted, hit_matrix, n_genes):
    """Compute ES for ALL gene sets simultaneously given sorted rank metrics.

    Parameters:
      rank_metrics_sorted: (n_genes,) — absolute rank metrics in sorted order
      hit_matrix: (n_gene_sets, n_genes) — boolean, True where gene is in set
      n_genes: int

    Returns:
      es_values: (n_gene_sets,) — enrichment score for each gene set
    """
    n_sets = hit_matrix.shape[0]

    # Number of hits per gene set
    n_hits = hit_matrix.sum(axis=1)  # (n_sets,)
    n_miss = n_genes - n_hits  # (n_sets,)

    # Weighted hit scores: |rank_metric| * hit indicator
    abs_metrics = np.abs(rank_metrics_sorted)  # (n_genes,)
    # (n_sets, n_genes) element-wise multiply
    hit_weights = hit_matrix * abs_metrics[np.newaxis, :]

    # Normalization factor per gene set
    nr = hit_weights.sum(axis=1)  # (n_sets,)
    # Avoid division by zero
    nr_safe = np.where(nr == 0, 1.0, nr)

    # Cumulative hit score
    p_hit = np.cumsum(hit_weights, axis=1) / nr_safe[:, np.newaxis]  # (n_sets, n_genes)

    # Cumulative miss score
    miss_matrix = ~hit_matrix
    n_miss_safe = np.where(n_miss == 0, 1, n_miss)
    p_miss = np.cumsum(miss_matrix.astype(float), axis=1) / n_miss_safe[:, np.newaxis]

    # Running ES
    running_es = p_hit - p_miss  # (n_sets, n_genes)

    # ES = max deviation from zero
    es_max = running_es.max(axis=1)  # (n_sets,)
    es_min = running_es.min(axis=1)  # (n_sets,)

    # Choose the one with larger absolute value
    es = np.where(np.abs(es_max) >= np.abs(es_min), es_max, es_min)

    # Gene sets with 0 hits get ES = 0
    es = np.where(n_hits == 0, 0.0, es)

    return es


def run_gsea_vectorized(ranked_df, gene_sets_dict, n_perm=500):
    """Run GSEA with gene-label permutation — FULLY VECTORIZED.

    Key optimizations:
      1. Pre-compute boolean hit matrix (n_sets × n_genes) once
      2. All gene sets computed simultaneously per permutation
      3. Permutation loop: only shuffle + argsort + matrix index (all NumPy)
      4. ~50-100x faster than nested Python loops
    """
    n_genes = len(ranked_df)
    gene_list = ranked_df['gene'].values
    rank_metrics = ranked_df['rank_metric'].values.astype(np.float64)

    gs_names = list(gene_sets_dict.keys())
    n_sets = len(gs_names)

    # Build hit matrix: (n_sets, n_genes) — True where gene is in gene set
    # Gene order matches ranked_df sorted order
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    hit_matrix = np.zeros((n_sets, n_genes), dtype=bool)
    set_sizes = np.zeros(n_sets, dtype=int)

    for s_idx, gs_name in enumerate(gs_names):
        gs_genes = gene_sets_dict[gs_name]
        for g in gs_genes:
            if g in gene_to_idx:
                hit_matrix[s_idx, gene_to_idx[g]] = True
        set_sizes[s_idx] = hit_matrix[s_idx].sum()

    # Compute observed ES (all gene sets at once)
    print(f"  Computing observed ES for {n_sets} gene sets...")
    observed_es = compute_es_batch_vectorized(rank_metrics, hit_matrix, n_genes)

    # Permutations: shuffle rank metrics, re-sort, recompute hit matrix positions, compute ES
    print(f"  Running {n_perm} permutations for {n_sets} gene sets (vectorized)...")
    perm_es_all = np.zeros((n_perm, n_sets))

    t_start = time.time()
    for i in range(n_perm):
        # Shuffle rank metrics
        perm_metrics = np.random.permutation(rank_metrics)

        # Sort indices (descending)
        sort_idx = np.argsort(-perm_metrics)

        # Re-order metrics and hit matrix by new sort order
        perm_metrics_sorted = perm_metrics[sort_idx]

        # Re-order hit matrix columns: gene at original position j moves to new position
        perm_hit_matrix = hit_matrix[:, sort_idx]

        # Compute ES for all gene sets simultaneously
        perm_es_all[i, :] = compute_es_batch_vectorized(
            perm_metrics_sorted, perm_hit_matrix, n_genes)

        if (i + 1) % 50 == 0:
            elapsed = time.time() - t_start
            rate = (i + 1) / elapsed
            remaining = (n_perm - i - 1) / rate
            print(f"    Perm {i+1}/{n_perm} "
                  f"({elapsed:.1f}s elapsed, ~{remaining:.0f}s remaining)")

    elapsed_total = time.time() - t_start
    print(f"  Permutations complete in {elapsed_total:.1f}s "
          f"({n_perm * n_sets:,} ES computations)")

    # Compute NES, p-value, FDR for each gene set
    results = []
    for s_idx, gs_name in enumerate(gs_names):
        obs = observed_es[s_idx]
        null = perm_es_all[:, s_idx]

        # NES: normalize by mean of same-sign permutation scores
        if obs >= 0:
            pos_null = null[null >= 0]
            mean_null = pos_null.mean() if len(pos_null) > 0 else 1.0
        else:
            neg_null = null[null < 0]
            mean_null = abs(neg_null.mean()) if len(neg_null) > 0 else 1.0

        nes = obs / mean_null if mean_null != 0 else 0.0

        # p-value: fraction of permutations with |ES| >= |observed| (same sign)
        if obs >= 0:
            pval = (null >= obs).sum() / max(1, (null >= 0).sum())
        else:
            pval = (null <= obs).sum() / max(1, (null < 0).sum())

        results.append({
            'gene_set': gs_name,
            'es': obs,
            'nes': nes,
            'pvalue': pval,
            'set_size': int(set_sizes[s_idx]),
        })

    res_df = pd.DataFrame(results)
    if len(res_df) > 0:
        # BH FDR correction
        res_df['fdr'] = np.minimum(
            res_df['pvalue'] * len(res_df) / rankdata(res_df['pvalue'].clip(lower=1e-300)),
            1.0)
        res_df = res_df.sort_values('pvalue')
    return res_df


# ============================================================
# 4. HELPER
# ============================================================
def get_samples(**kwargs):
    mask = pd.Series(True, index=meta.index)
    for col, val in kwargs.items():
        if isinstance(val, list):
            mask &= meta[col].isin(val)
        else:
            mask &= meta[col] == val
    return list(meta[mask].index)


# ============================================================
# 5. RUN GSEA ACROSS 8 SEASON × STAGE CONDITIONS
# ============================================================
print("\n" + "=" * 70)
print("RUNNING GSEA: 8 SEASON × STAGE CONDITIONS (VECTORIZED)")
print("=" * 70)

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
seasons = [2021, 2022]
stage_labels = {'Pre_veraison': 'Pre-vér', 'Veraison': 'Véraison',
                'Post_veraison': 'Post-vér', 'Harvest': 'Harvest'}

all_gsea_results = {}
all_ranked = {}

overall_start = time.time()
for season in seasons:
    for stage in stages:
        label = f"{season}_{stage}"
        print(f"\n--- {label} ---")

        infected = get_samples(Infection_Group='RB(+)',
                               Sampling_Season=season,
                               Maturity_stage=stage)
        control = get_samples(Infection_Group='RB(-)',
                              Sampling_Season=season,
                              Maturity_stage=stage)

        if len(infected) < 5 or len(control) < 3:
            print(f"  Skipping: infected={len(infected)}, control={len(control)}")
            continue

        print(f"  Samples: {len(infected)} RB(+) vs {len(control)} RB(-)")

        t0 = time.time()
        ranked = rank_genes_vectorized(infected, control, log_cpm)
        print(f"  Gene ranking: {time.time() - t0:.1f}s")
        all_ranked[label] = ranked

        t0 = time.time()
        gsea_res = run_gsea_vectorized(ranked, gene_sets, n_perm=500)
        print(f"  GSEA total: {time.time() - t0:.1f}s")
        all_gsea_results[label] = gsea_res

        n_sig = (gsea_res['fdr'] < 0.25).sum()
        n_sig05 = (gsea_res['fdr'] < 0.05).sum()
        print(f"  Significant gene sets: {n_sig} (FDR<0.25), {n_sig05} (FDR<0.05)")

        if n_sig05 > 0:
            top = gsea_res[gsea_res['fdr'] < 0.05].head(5)
            for _, r in top.iterrows():
                print(f"    {r['gene_set'][:55]:55s}  NES={r['nes']:+.2f}  FDR={r['fdr']:.3f}")

overall_elapsed = time.time() - overall_start
print(f"\n  Total GSEA time for all conditions: {overall_elapsed:.1f}s "
      f"({overall_elapsed/60:.1f} min)")


# ============================================================
# 6. ALSO RUN FISHER ORA FOR CONCORDANCE COMPARISON
# ============================================================
print("\n" + "=" * 70)
print("RUNNING FISHER ORA FOR CONCORDANCE COMPARISON")
print("=" * 70)

def fisher_enrichment_for_set(query_genes, term2genes_dict, background):
    """Fisher's exact test enrichment."""
    query = set(query_genes) & background
    n_q = len(query)
    n_bg = len(background)
    results = []
    for term, t_genes in term2genes_dict.items():
        t_in_bg = t_genes & background
        if len(t_in_bg) < MIN_SET:
            continue
        a = len(query & t_in_bg)
        if a < 2:
            continue
        b = n_q - a
        c = len(t_in_bg) - a
        d = n_bg - a - b - c
        _, pval = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        results.append({
            'term': term,
            'overlap': a,
            'term_size': len(t_in_bg),
            'pvalue': pval,
            'fold_enrichment': (a / n_q) / (len(t_in_bg) / n_bg) if n_q > 0 else 0
        })
    if not results:
        return pd.DataFrame()
    df = pd.DataFrame(results)
    df['padj'] = np.minimum(df['pvalue'] * len(df) / rankdata(df['pvalue']), 1.0)
    return df.sort_values('pvalue')

# Build combined gene set dict matching GSEA keys
combined_t2g = {}
for term, genes in go2genes.items():
    if MIN_SET <= len(genes) <= MAX_SET:
        name = godag[term].name if term in godag else term
        combined_t2g[f"GO:{name}"] = genes
for pathway, genes in kegg2genes.items():
    if MIN_SET <= len(genes) <= MAX_SET:
        name = kegg_names.get(pathway, pathway)
        combined_t2g[f"KEGG:{name}"] = genes

all_fisher_results = {}
for label, ranked in all_ranked.items():
    # BH correction on full ranked list
    ranked_bh = ranked.copy()
    ranked_bh['padj'] = np.minimum(
        ranked_bh['pvalue'] * len(ranked_bh) / rankdata(ranked_bh['pvalue']), 1.0)
    degs = ranked_bh[(ranked_bh['padj'] < 0.05) & (ranked_bh['log2FC'].abs() >= 1)]
    deg_genes = set(degs['gene'])
    if len(deg_genes) >= 10:
        fisher_res = fisher_enrichment_for_set(deg_genes, combined_t2g, expressed)
        all_fisher_results[label] = fisher_res
        print(f"  {label}: {len(deg_genes)} DEGs, {(fisher_res['padj'] < 0.05).sum() if len(fisher_res) > 0 else 0} sig terms (Fisher)")


# ============================================================
# 7. SAVE RESULTS
# ============================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

all_rows = []
for label, res in all_gsea_results.items():
    res_copy = res.copy()
    res_copy['condition'] = label
    all_rows.append(res_copy)

if all_rows:
    gsea_csv = pd.concat(all_rows, ignore_index=True)
    gsea_csv.to_csv('gsea_results_all_conditions.csv', index=False)
    print(f"Saved gsea_results_all_conditions.csv ({len(gsea_csv)} rows)")
else:
    gsea_csv = pd.DataFrame()
    print("No GSEA results to save")


# ============================================================
# 8. FIGURE 37 — GSEA Enrichment Plots
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 37: GSEA ENRICHMENT PLOTS")
print("=" * 70)

fig, axes = plt.subplots(2, 4, figsize=(28, 14))
fig.suptitle('GSEA Enrichment Scores Across Season × Stage',
             fontsize=20, fontweight='bold', y=0.98)

for col_idx, stage in enumerate(stages):
    for row_idx, season in enumerate(seasons):
        ax = axes[row_idx, col_idx]
        label = f"{season}_{stage}"

        if label not in all_gsea_results or label not in all_ranked:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                    transform=ax.transAxes, fontsize=14)
            ax.set_title(f"{season} {stage_labels.get(stage, stage)}", fontsize=14)
            continue

        gsea_res = all_gsea_results[label]
        ranked = all_ranked[label]

        # Get top 3 positive and top 3 negative by NES
        sig = gsea_res[gsea_res['fdr'] < 0.25]
        top_pos = sig[sig['nes'] > 0].nsmallest(3, 'pvalue') if len(sig[sig['nes'] > 0]) > 0 else pd.DataFrame()
        top_neg = sig[sig['nes'] < 0].nsmallest(3, 'pvalue') if len(sig[sig['nes'] < 0]) > 0 else pd.DataFrame()
        top_sets = pd.concat([top_pos, top_neg])

        if len(top_sets) == 0:
            # Show top by p-value regardless
            top_sets = gsea_res.head(4)

        colors_pos = ['#e74c3c', '#c0392b', '#a93226']
        colors_neg = ['#2980b9', '#2471a3', '#1a5276']

        pos_idx = 0
        neg_idx = 0
        for _, row in top_sets.iterrows():
            gs_name = row['gene_set']
            gs_genes = gene_sets.get(gs_name, set())
            if len(gs_genes) == 0:
                continue
            _, running_es = compute_enrichment_score(ranked, gs_genes)

            if row['nes'] >= 0:
                c = colors_pos[pos_idx % 3]
                pos_idx += 1
            else:
                c = colors_neg[neg_idx % 3]
                neg_idx += 1

            short_name = gs_name.split(':')[-1][:30]
            ax.plot(running_es, color=c, linewidth=1.2, alpha=0.8,
                    label=f"{short_name} ({row['nes']:+.1f})")

        ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
        ax.set_title(f"{season} {stage_labels.get(stage, stage)}", fontsize=14)
        if row_idx == 1:
            ax.set_xlabel('Gene rank', fontsize=12)
        if col_idx == 0:
            ax.set_ylabel('Running ES', fontsize=12)
        ax.legend(fontsize=7, loc='best', framealpha=0.7)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('37_gsea_enrichment_plots.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 37_gsea_enrichment_plots.png")


# ============================================================
# 9. FIGURE 38 — Fisher ORA vs GSEA Concordance
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 38: FISHER vs GSEA CONCORDANCE")
print("=" * 70)

fig, axes = plt.subplots(2, 4, figsize=(28, 14))
fig.suptitle('Fisher ORA vs GSEA Concordance Across Conditions',
             fontsize=20, fontweight='bold', y=0.98)

concordance_stats = []

for col_idx, stage in enumerate(stages):
    for row_idx, season in enumerate(seasons):
        ax = axes[row_idx, col_idx]
        label = f"{season}_{stage}"

        if label not in all_gsea_results or label not in all_fisher_results:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center',
                    transform=ax.transAxes, fontsize=12)
            ax.set_title(f"{season} {stage_labels.get(stage, stage)}", fontsize=14)
            continue

        gsea_res = all_gsea_results[label].set_index('gene_set')
        fisher_res = all_fisher_results[label].set_index('term')

        # Find common terms
        common = list(set(gsea_res.index) & set(fisher_res.index))
        if len(common) < 5:
            ax.text(0.5, 0.5, f'Only {len(common)} common\nterms', ha='center',
                    va='center', transform=ax.transAxes, fontsize=12)
            ax.set_title(f"{season} {stage_labels.get(stage, stage)}", fontsize=14)
            continue

        gsea_vals = -np.log10(gsea_res.loc[common, 'pvalue'].clip(lower=1e-30).values)
        fisher_vals = -np.log10(fisher_res.loc[common, 'padj'].clip(lower=1e-30).values)

        # Color by GSEA direction
        nes_vals = gsea_res.loc[common, 'nes'].values
        colors = ['#e74c3c' if n > 0 else '#2980b9' for n in nes_vals]

        ax.scatter(fisher_vals, gsea_vals, c=colors, alpha=0.4, s=20, edgecolors='none')

        # Correlation
        r, p = stats.spearmanr(fisher_vals, gsea_vals)
        concordance_stats.append({'condition': label, 'spearman_r': r, 'p': p, 'n_common': len(common)})

        ax.set_title(f"{season} {stage_labels.get(stage, stage)}\nr={r:.2f}", fontsize=14)
        ax.set_xlabel('-log10(Fisher padj)', fontsize=11)
        ax.set_ylabel('-log10(GSEA pval)', fontsize=11)

        # Diagonal
        lims = [0, max(fisher_vals.max(), gsea_vals.max()) * 1.1]
        ax.plot(lims, lims, 'k--', alpha=0.3, linewidth=0.5)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('38_gsea_vs_fisher_concordance.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved 38_gsea_vs_fisher_concordance.png")

if concordance_stats:
    print("\nConcordance summary (Spearman r):")
    for s in concordance_stats:
        print(f"  {s['condition']}: r={s['spearman_r']:.3f} (n={s['n_common']})")


# ============================================================
# 10. FIGURE 39 — NES Heatmap Across Conditions
# ============================================================
print("\n" + "=" * 70)
print("PLOTTING FIGURE 39: GSEA NES HEATMAP")
print("=" * 70)

# Build NES matrix: rows = gene sets, columns = conditions
conditions = [f"{s}_{st}" for s in seasons for st in stages]
conditions = [c for c in conditions if c in all_gsea_results]

if len(conditions) >= 2:
    # Collect all significant gene sets across conditions
    sig_sets = set()
    for c in conditions:
        res = all_gsea_results[c]
        sig = res[res['fdr'] < 0.25]
        sig_sets.update(sig['gene_set'].tolist())

    if len(sig_sets) > 60:
        # Take top terms from each condition
        sig_sets = set()
        for c in conditions:
            res = all_gsea_results[c]
            top = res.nsmallest(8, 'pvalue')
            sig_sets.update(top['gene_set'].tolist())

    sig_sets = sorted(sig_sets)

    if len(sig_sets) >= 3:
        nes_matrix = pd.DataFrame(0.0, index=sig_sets, columns=conditions)
        fdr_matrix = pd.DataFrame(1.0, index=sig_sets, columns=conditions)

        for c in conditions:
            res = all_gsea_results[c].set_index('gene_set')
            for gs in sig_sets:
                if gs in res.index:
                    nes_matrix.loc[gs, c] = res.loc[gs, 'nes']
                    fdr_matrix.loc[gs, c] = res.loc[gs, 'fdr']

        # Shorten labels for display
        short_labels = [s.split(':')[-1][:45] for s in nes_matrix.index]
        col_labels = [f"{c.split('_')[0]}\n{'_'.join(c.split('_')[1:])}" for c in conditions]

        fig, ax = plt.subplots(figsize=(14, max(8, len(sig_sets) * 0.35)))

        # Mask non-significant (FDR > 0.25) as lighter
        mask_arr = fdr_matrix.values > 0.25

        sns.heatmap(nes_matrix.values, ax=ax,
                    xticklabels=col_labels,
                    yticklabels=short_labels,
                    cmap='RdBu_r', center=0, vmin=-3, vmax=3,
                    linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'NES', 'shrink': 0.7})

        # Add significance markers
        for i in range(len(sig_sets)):
            for j in range(len(conditions)):
                fdr_val = fdr_matrix.iloc[i, j]
                if fdr_val < 0.05:
                    ax.text(j + 0.5, i + 0.5, '***', ha='center', va='center',
                            fontsize=8, fontweight='bold', color='black')
                elif fdr_val < 0.10:
                    ax.text(j + 0.5, i + 0.5, '**', ha='center', va='center',
                            fontsize=8, color='black')
                elif fdr_val < 0.25:
                    ax.text(j + 0.5, i + 0.5, '*', ha='center', va='center',
                            fontsize=8, color='black')

        ax.set_title('GSEA Normalized Enrichment Scores (NES)\nAcross Season × Stage Conditions',
                     fontsize=18, fontweight='bold', pad=15)
        ax.set_xlabel('Condition', fontsize=14)
        ax.tick_params(axis='y', labelsize=10)
        ax.tick_params(axis='x', labelsize=11)

        plt.tight_layout()
        plt.savefig('39_gsea_nes_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved 39_gsea_nes_heatmap.png ({len(sig_sets)} gene sets × {len(conditions)} conditions)")
    else:
        print("Too few significant gene sets for heatmap")
else:
    print("Too few conditions for heatmap")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("GSEA ANALYSIS COMPLETE (VECTORIZED)")
print("=" * 70)

total_sig = 0
for c, res in all_gsea_results.items():
    n = (res['fdr'] < 0.05).sum()
    total_sig += n
    print(f"  {c}: {n} gene sets at FDR<0.05")

print(f"\nTotal significant gene set-condition pairs: {total_sig}")
print(f"Total runtime: {(time.time() - overall_start)/60:.1f} min")
print(f"Output files:")
print(f"  - 37_gsea_enrichment_plots.png")
print(f"  - 38_gsea_vs_fisher_concordance.png")
print(f"  - 39_gsea_nes_heatmap.png")
print(f"  - gsea_results_all_conditions.csv")
