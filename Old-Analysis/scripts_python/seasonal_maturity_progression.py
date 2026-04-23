#!/usr/bin/env python3
"""
Seasonal Maturity Progression Analysis for GRBV Study.
Examines how GRBV-induced transcriptomic changes progress through maturity stages
in 2021 vs 2022, with DE analysis + GO/KEGG pathway enrichment.
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
from scipy.stats import fisher_exact, rankdata
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from collections import defaultdict
import textwrap
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'figure.dpi': 300, 'font.size': 14, 'axes.titlesize': 16,
    'axes.labelsize': 15, 'figure.facecolor': 'white'
})

# ============================================================
# 1. LOAD DATA & ANNOTATIONS
# ============================================================
print("=" * 70)
print("LOADING DATA AND ANNOTATIONS")
print("=" * 70)

counts = pd.read_csv('counts_matrix.txt', sep='\t', index_col=0)
meta = pd.read_csv('sample_metadata.txt', sep='\t')

common_samples = list(set(counts.columns) & set(meta['Sample']))
counts = counts[common_samples]
meta = meta[meta['Sample'].isin(common_samples)].set_index('Sample').loc[common_samples]

# Filter genes
min_samples = int(0.1 * counts.shape[1])
gene_mask = (counts >= 10).sum(axis=1) >= min_samples
counts_filtered = counts.loc[gene_mask]

# Normalize
lib_sizes = counts_filtered.sum(axis=0)
cpm = counts_filtered.div(lib_sizes, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

print(f"Filtered genes: {counts_filtered.shape[0]}")
print(f"Total samples: {counts_filtered.shape[1]}")

# Load annotations
db_path = 'org.Vvinifera.eg.db/inst/extdata/org.Vvinifera.eg.sqlite'
conn = sqlite3.connect(db_path)

go_df = pd.read_sql("SELECT g.SYMBOL, go.GO, go.ONTOLOGY FROM go JOIN gene_info g ON go._id = g._id", conn)
kegg_df = pd.read_sql("SELECT g.SYMBOL, k.KEGG FROM kegg k JOIN gene_info g ON k._id = g._id", conn)

import urllib.request
try:
    url = "https://rest.kegg.jp/list/pathway/vvi"
    response = urllib.request.urlopen(url, timeout=30)
    kegg_names = {}
    for line in response.read().decode().strip().split('\n'):
        parts = line.split('\t')
        if len(parts) == 2:
            path_id = parts[0].replace('path:', '')
            kegg_names[path_id] = parts[1].replace(' - Vitis vinifera (wine grape)', '').strip()
    print(f"KEGG pathway names loaded: {len(kegg_names)}")
except Exception as e:
    print(f"KEGG names fetch failed: {e}")
    kegg_names = {}

from goatools.obo_parser import GODag
godag = GODag('go-basic.obo')
print(f"GO terms loaded: {len(godag)}")

conn.close()

# Build annotation dictionaries
gene2go_bp, gene2go_mf, gene2go_cc, gene2kegg = defaultdict(set), defaultdict(set), defaultdict(set), defaultdict(set)
go2genes_bp, go2genes_mf, go2genes_cc, kegg2genes = defaultdict(set), defaultdict(set), defaultdict(set), defaultdict(set)

for _, row in go_df.iterrows():
    gene, term, ont = row['SYMBOL'], row['GO'], row['ONTOLOGY']
    if ont == 'BP':
        gene2go_bp[gene].add(term); go2genes_bp[term].add(gene)
    elif ont == 'MF':
        gene2go_mf[gene].add(term); go2genes_mf[term].add(gene)
    elif ont == 'CC':
        gene2go_cc[gene].add(term); go2genes_cc[term].add(gene)

for _, row in kegg_df.iterrows():
    gene2kegg[row['SYMBOL']].add(row['KEGG']); kegg2genes[row['KEGG']].add(row['SYMBOL'])

all_genes = set(counts_filtered.index)
bg_go_bp = all_genes & set(gene2go_bp.keys())
bg_go_mf = all_genes & set(gene2go_mf.keys())
bg_go_cc = all_genes & set(gene2go_cc.keys())
bg_kegg = all_genes & set(gene2kegg.keys())

print(f"Background: GO-BP={len(bg_go_bp)}, GO-MF={len(bg_go_mf)}, GO-CC={len(bg_go_cc)}, KEGG={len(bg_kegg)}")


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def run_de(group1_samples, group2_samples):
    """T-test DE with BH correction."""
    g1, g2 = log_cpm[group1_samples], log_cpm[group2_samples]
    results = []
    for gene in log_cpm.index:
        v1, v2 = g1.loc[gene].values, g2.loc[gene].values
        if np.std(v1) == 0 and np.std(v2) == 0:
            continue
        t, p = stats.ttest_ind(v1, v2, equal_var=False)
        results.append({'gene': gene, 'log2FC': v1.mean() - v2.mean(), 'pvalue': p})
    df = pd.DataFrame(results)
    if len(df) == 0:
        return df
    df['padj'] = np.minimum(df['pvalue'] * len(df) / rankdata(df['pvalue']), 1.0)
    return df.sort_values('pvalue')


def get_degs(de_df, fc=1.0, padj=0.05):
    return de_df[(de_df['padj'] < padj) & (de_df['log2FC'].abs() >= fc)]


def get_up_down(de_df, fc=1.0, padj=0.05):
    sig = get_degs(de_df, fc, padj)
    return sig[sig['log2FC'] > 0]['gene'].tolist(), sig[sig['log2FC'] < 0]['gene'].tolist()


def get_samples(meta_df, **kwargs):
    mask = pd.Series(True, index=meta_df.index)
    for col, val in kwargs.items():
        mask &= meta_df[col].isin(val) if isinstance(val, list) else meta_df[col] == val
    return list(meta_df[mask].index)


def fisher_enrichment(query_genes, term2genes, background_genes, min_genes=3):
    query = set(query_genes) & background_genes
    n_query, n_bg = len(query), len(background_genes)
    results = []
    for term, term_genes in term2genes.items():
        term_in_bg = term_genes & background_genes
        if len(term_in_bg) < min_genes:
            continue
        a = len(query & term_in_bg)
        if a < 2:
            continue
        b = len(query - term_in_bg)
        c = len(term_in_bg - query)
        d = n_bg - len(query | term_in_bg)
        odds, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        results.append({
            'term': term, 'overlap': a, 'term_size': len(term_in_bg),
            'query_size': n_query, 'bg_size': n_bg,
            'fold_enrichment': (a/n_query)/(len(term_in_bg)/n_bg) if n_query > 0 and len(term_in_bg) > 0 else 0,
            'pvalue': pval
        })
    if not results:
        return pd.DataFrame()
    df = pd.DataFrame(results)
    df['padj'] = np.minimum(df['pvalue'] * len(df) / rankdata(df['pvalue']), 1.0)
    return df.sort_values('pvalue')


def run_full_enrichment(query_genes):
    query = set(query_genes) & all_genes
    res = {}
    for name, t2g, bg in [('GO_BP', go2genes_bp, bg_go_bp), ('GO_MF', go2genes_mf, bg_go_mf),
                            ('GO_CC', go2genes_cc, bg_go_cc), ('KEGG', kegg2genes, bg_kegg)]:
        df = fisher_enrichment(query, t2g, bg)
        if len(df) > 0:
            if name.startswith('GO'):
                df['term_name'] = df['term'].apply(lambda x: godag[x].name if x in godag else x)
            else:
                df['term_name'] = df['term'].apply(lambda x: kegg_names.get(x, x))
            df['category'] = name
        res[name] = df
    return res


def get_combined_top(enrichment_dict, n_per_cat=5, padj_thresh=0.05):
    combined = []
    for cat in ['GO_BP', 'GO_MF', 'GO_CC', 'KEGG']:
        df = enrichment_dict.get(cat, pd.DataFrame())
        if len(df) > 0:
            combined.append(df[df['padj'] < padj_thresh].head(n_per_cat))
    return pd.concat(combined, ignore_index=True) if combined else pd.DataFrame()


stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
stage_colors = {'Pre_veraison': '#2ca02c', 'Veraison': '#ff7f0e', 'Post_veraison': '#d62728', 'Harvest': '#9467bd'}
seasons = [2021, 2022]


# ============================================================
# 3. CORE ANALYSIS: SEASON × STAGE DE + ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("SEASONAL MATURITY PROGRESSION: DE ANALYSIS")
print("=" * 70)

# Store all results
de_results = {}      # (season, stage) -> de_df
enrichments = {}     # (season, stage, direction) -> enrichment_dict
deg_counts = {}      # (season, stage) -> {total, up, down}
deg_gene_sets = {}   # (season, stage, direction) -> set of genes

for season in seasons:
    print(f"\n{'='*50}")
    print(f"SEASON {season}")
    print(f"{'='*50}")
    for stage in stages:
        pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season, Maturity_stage=stage)
        neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season, Maturity_stage=stage)

        print(f"\n  {stage}: RB(+)={len(pos)}, RB(-)={len(neg)}")

        if len(pos) < 3 or len(neg) < 3:
            print(f"    Skipped (insufficient samples)")
            continue

        de = run_de(pos, neg)
        de_results[(season, stage)] = de
        up, down = get_up_down(de)

        deg_counts[(season, stage)] = {'total': len(up)+len(down), 'up': len(up), 'down': len(down)}
        deg_gene_sets[(season, stage, 'up')] = set(up)
        deg_gene_sets[(season, stage, 'down')] = set(down)

        print(f"    DEGs: {len(up)+len(down)} total (UP={len(up)}, DOWN={len(down)})")

        for direction, genes in [('up', up), ('down', down)]:
            if len(genes) >= 5:
                enr = run_full_enrichment(genes)
                enrichments[(season, stage, direction)] = enr
                # Print top terms
                combined = get_combined_top(enr, n_per_cat=3)
                if len(combined) > 0:
                    sig = combined[combined['padj'] < 0.05]
                    if len(sig) > 0:
                        terms = '; '.join(sig['term_name'].head(5).tolist())
                        print(f"    {direction.upper()} pathways: {terms}")
            else:
                enrichments[(season, stage, direction)] = {}


# ============================================================
# 4. CROSS-SEASON COMPARISON AT EACH STAGE
# ============================================================
print("\n" + "=" * 70)
print("CROSS-SEASON COMPARISON AT EACH STAGE")
print("=" * 70)

# FC correlation between seasons at each stage
fc_correlations = {}
for stage in stages:
    de1 = de_results.get((2021, stage))
    de2 = de_results.get((2022, stage))
    if de1 is None or de2 is None:
        continue

    de1_idx = de1.set_index('gene')
    de2_idx = de2.set_index('gene')
    common = list(set(de1_idx.index) & set(de2_idx.index))

    r, p = stats.pearsonr(de1_idx.loc[common, 'log2FC'], de2_idx.loc[common, 'log2FC'])
    fc_correlations[stage] = {'r': r, 'p': p, 'n': len(common)}
    print(f"\n  {stage}: FC correlation (2021 vs 2022) r={r:.4f} (p={p:.2e})")

    # DEG overlap between seasons
    up_21 = deg_gene_sets.get((2021, stage, 'up'), set())
    up_22 = deg_gene_sets.get((2022, stage, 'up'), set())
    down_21 = deg_gene_sets.get((2021, stage, 'down'), set())
    down_22 = deg_gene_sets.get((2022, stage, 'down'), set())

    up_shared = up_21 & up_22
    down_shared = down_21 & down_22
    up_only21 = up_21 - up_22
    up_only22 = up_22 - up_21
    down_only21 = down_21 - down_22
    down_only22 = down_22 - down_21

    print(f"    UP genes:   2021={len(up_21)}, 2022={len(up_22)}, Shared={len(up_shared)}, "
          f"2021-only={len(up_only21)}, 2022-only={len(up_only22)}")
    print(f"    DOWN genes: 2021={len(down_21)}, 2022={len(down_22)}, Shared={len(down_shared)}, "
          f"2021-only={len(down_only21)}, 2022-only={len(down_only22)}")


# ============================================================
# 5. STAGE TRANSITIONS WITHIN EACH SEASON (RB+ vs RB-)
# ============================================================
print("\n" + "=" * 70)
print("STAGE TRANSITION ENRICHMENT (How GRBV effect changes with maturity)")
print("=" * 70)

# Track how infection effect changes from one stage to the next within each season
transition_enrichments = {}
transition_deg_counts = {}

for season in seasons:
    print(f"\n--- Season {season} ---")
    for i in range(len(stages) - 1):
        s1, s2 = stages[i], stages[i+1]

        # Infection effect at stage 1
        pos1 = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season, Maturity_stage=s1)
        neg1 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season, Maturity_stage=s1)
        pos2 = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season, Maturity_stage=s2)
        neg2 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season, Maturity_stage=s2)

        if len(pos1) < 3 or len(neg1) < 3 or len(pos2) < 3 or len(neg2) < 3:
            continue

        # Infection effect = mean(RB+) - mean(RB-) at each stage
        effect1 = log_cpm[pos1].mean(axis=1) - log_cpm[neg1].mean(axis=1)
        effect2 = log_cpm[pos2].mean(axis=1) - log_cpm[neg2].mean(axis=1)

        # Change in infection effect from s1 to s2
        delta_effect = effect2 - effect1

        # Use permutation-like approach: genes with large delta
        # Identify genes whose GRBV effect intensifies or weakens
        intensifying = delta_effect[delta_effect > 1.0].index.tolist()  # Effect becomes more positive
        weakening = delta_effect[delta_effect < -1.0].index.tolist()     # Effect becomes more negative

        trans_key = f"{s1} → {s2}"
        transition_deg_counts[(season, trans_key)] = {
            'intensifying': len(intensifying), 'weakening': len(weakening)
        }

        print(f"  {trans_key}: Effect intensifying={len(intensifying)}, weakening={len(weakening)}")

        for direction, genes in [('intensifying', intensifying), ('weakening', weakening)]:
            if len(genes) >= 10:
                enr = run_full_enrichment(genes)
                transition_enrichments[(season, trans_key, direction)] = enr
                combined = get_combined_top(enr, n_per_cat=2)
                if len(combined) > 0:
                    sig = combined[combined['padj'] < 0.05]
                    if len(sig) > 0:
                        terms = '; '.join(sig['term_name'].head(3).tolist())
                        print(f"    {direction}: {terms}")


# ============================================================
# 6. PATHWAY LANDSCAPE COMPARISON: 2021 vs 2022
# ============================================================
print("\n" + "=" * 70)
print("BUILDING PATHWAY LANDSCAPE MATRICES")
print("=" * 70)

# Collect all significant pathway terms across both seasons and all stages
pathway_data_by_season = {2021: {}, 2022: {}}

for season in seasons:
    for stage in stages:
        for direction in ['up', 'down']:
            enr = enrichments.get((season, stage, direction), {})
            for cat in ['GO_BP', 'GO_MF', 'KEGG']:
                df = enr.get(cat, pd.DataFrame())
                if len(df) > 0:
                    for _, row in df[df['padj'] < 0.05].head(8).iterrows():
                        key = f"[{cat}] {row['term_name']}"
                        if key not in pathway_data_by_season[season]:
                            pathway_data_by_season[season][key] = {}
                        sign = 1 if direction == 'up' else -1
                        pathway_data_by_season[season][key][stage] = sign * (-np.log10(max(row['padj'], 1e-30)))


# ============================================================
# PLOT 1: DEG COUNTS COMPARISON (2021 vs 2022)
# ============================================================
print("\n" + "=" * 70)
print("GENERATING PLOTS")
print("=" * 70)

fig = plt.figure(figsize=(24, 20))
gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.45)

# Row 0, Col 0-1: Side-by-side DEG barplots
for s_idx, season in enumerate(seasons):
    ax = fig.add_subplot(gs[0, s_idx])
    x = np.arange(len(stages))
    ups = [deg_counts.get((season, st), {}).get('up', 0) for st in stages]
    downs = [-deg_counts.get((season, st), {}).get('down', 0) for st in stages]

    ax.bar(x, ups, color='firebrick', alpha=0.85, label='Upregulated', width=0.6)
    ax.bar(x, downs, color='steelblue', alpha=0.85, label='Downregulated', width=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=12)
    ax.set_ylabel("Number of DEGs")
    ax.set_title(f"Season {season}\nDEGs: RB(+) vs RB(-)", fontsize=16)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.legend(fontsize=12)

    # Annotate counts
    for i, (u, d) in enumerate(zip(ups, downs)):
        if u > 0:
            ax.text(i, u + 5, str(u), ha='center', fontsize=10, color='firebrick')
        if d < 0:
            ax.text(i, d - 15, str(abs(d)), ha='center', fontsize=10, color='steelblue')

# Row 0, Col 2: DEG comparison grouped bars
ax = fig.add_subplot(gs[0, 2])
x = np.arange(len(stages))
width = 0.35
totals_21 = [deg_counts.get((2021, st), {}).get('total', 0) for st in stages]
totals_22 = [deg_counts.get((2022, st), {}).get('total', 0) for st in stages]
ax.bar(x - width/2, totals_21, width, color='#3498db', alpha=0.85, label='2021')
ax.bar(x + width/2, totals_22, width, color='#e67e22', alpha=0.85, label='2022')
ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=12)
ax.set_ylabel("Total DEGs")
ax.set_title("Total DEGs: 2021 vs 2022", fontsize=16)
ax.legend()

# Row 0, Col 3: FC correlation between seasons
ax = fig.add_subplot(gs[0, 3])
stage_labels = [st for st in stages if st in fc_correlations]
r_vals = [fc_correlations[st]['r'] for st in stage_labels]
bars = ax.bar(range(len(stage_labels)), r_vals,
              color=[stage_colors[st] for st in stage_labels], alpha=0.85)
ax.set_xticks(range(len(stage_labels)))
ax.set_xticklabels([s.replace('_', '\n') for s in stage_labels], fontsize=12)
ax.set_ylabel("Pearson r")
ax.set_title("FC Consistency\n(2021 vs 2022 per stage)", fontsize=16)
ax.set_ylim(0, 1)
ax.axhline(0.8, color='grey', linestyle='--', linewidth=0.5)
for i, r in enumerate(r_vals):
    ax.text(i, r + 0.02, f'{r:.3f}', ha='center', fontsize=12)

# Row 1: DEG Overlap Venn-like plots per stage (2021 vs 2022)
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[1, idx])

    up_21 = deg_gene_sets.get((2021, stage, 'up'), set())
    up_22 = deg_gene_sets.get((2022, stage, 'up'), set())
    down_21 = deg_gene_sets.get((2021, stage, 'down'), set())
    down_22 = deg_gene_sets.get((2022, stage, 'down'), set())

    categories = ['UP\n2021 only', 'UP\nShared', 'UP\n2022 only',
                   'DOWN\n2021 only', 'DOWN\nShared', 'DOWN\n2022 only']
    values = [len(up_21 - up_22), len(up_21 & up_22), len(up_22 - up_21),
              len(down_21 - down_22), len(down_21 & down_22), len(down_22 - down_21)]
    colors = ['#aed6f1', '#2ecc71', '#f5cba7',
              '#aed6f1', '#2ecc71', '#f5cba7']

    bars = ax.bar(range(len(categories)), values, color=colors, edgecolor='black', linewidth=0.5, alpha=0.85)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, fontsize=10)
    ax.set_ylabel("Gene count")
    ax.set_title(f"{stage}\nDEG Overlap 2021 vs 2022", fontsize=15)

    for i, v in enumerate(values):
        if v > 0:
            ax.text(i, v + max(values)*0.03, str(v), ha='center', fontsize=10)

# Row 2: FC scatter plots (2021 vs 2022) for each stage
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[2, idx])

    de1 = de_results.get((2021, stage))
    de2 = de_results.get((2022, stage))

    if de1 is None or de2 is None:
        ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
        continue

    de1_idx = de1.set_index('gene')
    de2_idx = de2.set_index('gene')
    common = list(set(de1_idx.index) & set(de2_idx.index))

    fc1 = de1_idx.loc[common, 'log2FC']
    fc2 = de2_idx.loc[common, 'log2FC']

    # Color: sig in both = red, sig in one = orange, neither = grey
    sig1 = set(get_degs(de1)['gene'])
    sig2 = set(get_degs(de2)['gene'])

    colors = []
    for g in common:
        if g in sig1 and g in sig2:
            colors.append('#e74c3c')  # sig in both
        elif g in sig1 or g in sig2:
            colors.append('#f39c12')  # sig in one
        else:
            colors.append('#bdc3c7')  # neither

    ax.scatter(fc1.values, fc2.values, c=colors, s=2, alpha=0.4)
    ax.plot([-6, 6], [-6, 6], 'k--', linewidth=0.5)
    ax.axhline(0, color='grey', linewidth=0.3)
    ax.axvline(0, color='grey', linewidth=0.3)

    r = fc_correlations.get(stage, {}).get('r', 0)
    ax.set_xlabel("log2FC (2021)")
    ax.set_ylabel("log2FC (2022)")
    ax.set_title(f"{stage}\nr = {r:.3f}", fontsize=15)
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)

fig.suptitle("GRBV Infection Response Across Maturity Stages: 2021 vs 2022\n"
             "(DEGs = RB+ vs RB- at each Season × Stage combination)",
             fontsize=21, y=1.02)
plt.savefig('15_seasonal_maturity_DEG_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 15_seasonal_maturity_DEG_comparison.png")


# ============================================================
# PLOT 2: PATHWAY ENRICHMENT PER SEASON × STAGE
# ============================================================

fig = plt.figure(figsize=(26, 28))
gs = gridspec.GridSpec(4, 4, figure=fig, hspace=0.55, wspace=0.6)

cat_colors = {'GO_BP': '#e74c3c', 'GO_MF': '#3498db', 'GO_CC': '#2ecc71', 'KEGG': '#9b59b6'}

def plot_enrichment_bars(enrichment_dict, ax, title="", n=10, padj_thresh=0.05):
    combined = get_combined_top(enrichment_dict, n_per_cat=4, padj_thresh=padj_thresh)
    if len(combined) == 0:
        ax.text(0.5, 0.5, 'No significant\nenrichment', ha='center', va='center',
                transform=ax.transAxes, fontsize=14, color='grey')
        ax.set_title(title, fontsize=14)
        return

    combined = combined.head(n)
    labels = [textwrap.fill(str(row['term_name'])[:50], 35) for _, row in combined.iterrows()]
    colors = [cat_colors.get(row['category'], 'grey') for _, row in combined.iterrows()]
    y = np.arange(len(combined))

    ax.barh(y, -np.log10(combined['padj'].clip(lower=1e-30).values), color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('-log10(padj)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)

    for i, (_, row) in enumerate(combined.iterrows()):
        ax.text(0.2, i, f" {int(row['overlap'])}/{int(row['term_size'])}",
                va='center', fontsize=8, color='white', fontweight='bold')

# 2021 UP (Row 0)
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[0, idx])
    enr = enrichments.get((2021, stage, 'up'), {})
    plot_enrichment_bars(enr, ax, f"2021 — {stage}\nUP in RB(+)", n=10)

# 2021 DOWN (Row 1)
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[1, idx])
    enr = enrichments.get((2021, stage, 'down'), {})
    plot_enrichment_bars(enr, ax, f"2021 — {stage}\nDOWN in RB(+)", n=10)

# 2022 UP (Row 2)
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[2, idx])
    enr = enrichments.get((2022, stage, 'up'), {})
    plot_enrichment_bars(enr, ax, f"2022 — {stage}\nUP in RB(+)", n=10)

# 2022 DOWN (Row 3)
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[3, idx])
    enr = enrichments.get((2022, stage, 'down'), {})
    plot_enrichment_bars(enr, ax, f"2022 — {stage}\nDOWN in RB(+)", n=10)

fig.suptitle("Pathway Enrichment: GRBV Response Progression Through Maturity Stages\n"
             "Top: 2021 Season | Bottom: 2022 Season\n"
             "(Red=GO:BP, Blue=GO:MF, Green=GO:CC, Purple=KEGG)",
             fontsize=21, y=1.02)
plt.savefig('16_seasonal_maturity_pathway_enrichment.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 16_seasonal_maturity_pathway_enrichment.png")


# ============================================================
# PLOT 3: PATHWAY LANDSCAPE HEATMAPS (2021 vs 2022 side by side)
# ============================================================

fig = plt.figure(figsize=(24, 20))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.4)

for s_idx, season in enumerate(seasons):
    ax = fig.add_subplot(gs[0, s_idx])

    data = pathway_data_by_season[season]
    if data:
        hm = pd.DataFrame(data).T.fillna(0)
        hm = hm.reindex(columns=[s for s in stages if s in hm.columns])
        # Select most variable terms
        if len(hm) > 30:
            hm_var = hm.var(axis=1).nlargest(30).index
            hm = hm.loc[hm_var]

        sns.heatmap(hm, ax=ax, cmap='RdBu_r', center=0,
                    yticklabels=[textwrap.fill(str(t)[:55], 45) for t in hm.index],
                    xticklabels=True,
                    cbar_kws={'label': 'Signed -log10(padj)\n(+Up in RB+, -Down)', 'shrink': 0.8},
                    linewidths=0.5)
        ax.tick_params(axis='y', labelsize=6.5)
    ax.set_title(f"Season {season}: Pathway Landscape\n(RB+ vs RB- across maturity stages)", fontsize=16)

# Bottom: Combined heatmap with both seasons side by side
ax = fig.add_subplot(gs[1, :])

# Build combined matrix: columns = (Season, Stage)
all_terms_combined = set()
for season in seasons:
    all_terms_combined.update(pathway_data_by_season[season].keys())

combined_matrix = pd.DataFrame(0.0, index=sorted(all_terms_combined),
                                columns=[f"{s}\n{st}" for s in seasons for st in stages])

for season in seasons:
    for term, stage_vals in pathway_data_by_season[season].items():
        for stage, val in stage_vals.items():
            combined_matrix.loc[term, f"{season}\n{stage}"] = val

# Remove rows with all zeros
combined_matrix = combined_matrix.loc[(combined_matrix != 0).any(axis=1)]

# Select most variable
if len(combined_matrix) > 35:
    var_idx = combined_matrix.var(axis=1).nlargest(35).index
    combined_matrix = combined_matrix.loc[var_idx]

if len(combined_matrix) > 0:
    sns.heatmap(combined_matrix, ax=ax, cmap='RdBu_r', center=0,
                yticklabels=[textwrap.fill(str(t)[:55], 45) for t in combined_matrix.index],
                xticklabels=True,
                cbar_kws={'label': 'Signed -log10(padj)', 'shrink': 0.6},
                linewidths=0.5)
    ax.tick_params(axis='y', labelsize=6)
    # Add separator line between seasons
    ax.axvline(x=4, color='black', linewidth=2)
    ax.set_title("Combined Pathway Landscape: 2021 vs 2022 Across Maturity Stages", fontsize=18)

fig.suptitle("Pathway Enrichment Landscape: Season-Specific Maturity Progression",
             fontsize=21, y=1.01)
plt.savefig('17_seasonal_maturity_pathway_landscape.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 17_seasonal_maturity_pathway_landscape.png")


# ============================================================
# PLOT 4: STAGE TRANSITION ANALYSIS + PCA
# ============================================================

fig = plt.figure(figsize=(24, 18))
gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.5)

# Row 0: PCA by season, split into subplots with stage coloring
for s_idx, season in enumerate(seasons):
    ax = fig.add_subplot(gs[0, s_idx*2])

    # Get samples for this season (infected only)
    season_mask = (meta['Sampling_Season'] == season) & (meta['Infection_Group'] == 'RB(+)')
    season_samples = meta[season_mask].index.tolist()

    if len(season_samples) > 5:
        pca_data_s = log_cpm.loc[log_cpm.var(axis=1).nlargest(2000).index, season_samples].T
        scaler_s = StandardScaler()
        pca_s = PCA(n_components=3)
        coords_s = pca_s.fit_transform(scaler_s.fit_transform(pca_data_s))

        for stage in stages:
            mask = meta.loc[season_samples, 'Maturity_stage'] == stage
            idx = np.where(mask)[0]
            if len(idx) > 0:
                ax.scatter(coords_s[idx, 0], coords_s[idx, 1], c=stage_colors[stage],
                          s=30, alpha=0.6, label=stage, edgecolors='none')

        ax.set_xlabel(f"PC1 ({pca_s.explained_variance_ratio_[0]:.1%})")
        ax.set_ylabel(f"PC2 ({pca_s.explained_variance_ratio_[1]:.1%})")
        ax.set_title(f"{season} RB(+)\nby Maturity Stage", fontsize=15)
        ax.legend(fontsize=10)

    # Same for controls
    ax2 = fig.add_subplot(gs[0, s_idx*2 + 1])
    season_neg_mask = (meta['Sampling_Season'] == season) & (meta['Infection_Group'] == 'RB(-)')
    season_neg = meta[season_neg_mask].index.tolist()

    if len(season_neg) > 5:
        pca_data_n = log_cpm.loc[log_cpm.var(axis=1).nlargest(2000).index, season_neg].T
        scaler_n = StandardScaler()
        pca_n = PCA(n_components=3)
        coords_n = pca_n.fit_transform(scaler_n.fit_transform(pca_data_n))

        for stage in stages:
            mask = meta.loc[season_neg, 'Maturity_stage'] == stage
            idx = np.where(mask)[0]
            if len(idx) > 0:
                ax2.scatter(coords_n[idx, 0], coords_n[idx, 1], c=stage_colors[stage],
                           s=30, alpha=0.6, label=stage, edgecolors='none')

        ax2.set_xlabel(f"PC1 ({pca_n.explained_variance_ratio_[0]:.1%})")
        ax2.set_ylabel(f"PC2 ({pca_n.explained_variance_ratio_[1]:.1%})")
        ax2.set_title(f"{season} RB(-)\nby Maturity Stage", fontsize=15)
        ax2.legend(fontsize=10)

# Row 1: Stage transition — change in GRBV effect
transitions = [f'{stages[i]} → {stages[i+1]}' for i in range(len(stages)-1)]

for s_idx, season in enumerate(seasons):
    ax = fig.add_subplot(gs[1, s_idx])

    intens = [transition_deg_counts.get((season, t), {}).get('intensifying', 0) for t in transitions]
    weak = [-transition_deg_counts.get((season, t), {}).get('weakening', 0) for t in transitions]

    x = np.arange(len(transitions))
    ax.bar(x, intens, color='#e74c3c', alpha=0.85, label='Effect intensifying')
    ax.bar(x, weak, color='#3498db', alpha=0.85, label='Effect weakening')
    ax.set_xticks(x)
    ax.set_xticklabels([t.replace(' → ', '\n→\n') for t in transitions], fontsize=12)
    ax.set_ylabel("Number of genes")
    ax.set_title(f"Season {season}\nChange in GRBV Effect Across Stages", fontsize=15)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.legend(fontsize=10)

# Row 1, Col 2-3: Transition enrichment comparison
for s_idx, season in enumerate(seasons):
    ax = fig.add_subplot(gs[1, s_idx + 2])

    # Collect transition pathway data
    trans_terms = {}
    for trans in transitions:
        for direction in ['intensifying', 'weakening']:
            enr = transition_enrichments.get((season, trans, direction), {})
            for cat in ['GO_BP', 'KEGG']:
                df = enr.get(cat, pd.DataFrame())
                if len(df) > 0:
                    for _, row in df[df['padj'] < 0.05].head(5).iterrows():
                        key = row['term_name']
                        if key not in trans_terms:
                            trans_terms[key] = {}
                        sign = 1 if direction == 'intensifying' else -1
                        trans_terms[key][trans] = sign * (-np.log10(max(row['padj'], 1e-30)))

    if trans_terms:
        hm = pd.DataFrame(trans_terms).T.fillna(0)
        hm = hm.reindex(columns=[t for t in transitions if t in hm.columns])
        if len(hm) > 15:
            hm = hm.loc[hm.var(axis=1).nlargest(15).index]

        sns.heatmap(hm, ax=ax, cmap='RdBu_r', center=0,
                    yticklabels=[textwrap.fill(str(t)[:40], 30) for t in hm.index],
                    xticklabels=[t.replace(' → ', '→\n') for t in hm.columns],
                    cbar_kws={'label': 'Signed -log10(padj)', 'shrink': 0.7},
                    linewidths=0.5)
        ax.tick_params(axis='y', labelsize=6)
        ax.set_title(f"{season}: Transition Pathways\n(+Red=intensifying, -Blue=weakening)", fontsize=14)
    else:
        ax.text(0.5, 0.5, 'No significant\ntransition pathways', ha='center', va='center',
                transform=ax.transAxes)
        ax.set_title(f"{season}: Transition Pathways", fontsize=14)

# Row 2: Season-specific unique DEGs at each stage — what's different between 2021 and 2022?
# Build a summary: per stage, which pathways are unique to each season
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[2, idx])

    # Get season-specific DEGs
    shared_up = deg_gene_sets.get((2021, stage, 'up'), set()) & deg_gene_sets.get((2022, stage, 'up'), set())
    only21_up = deg_gene_sets.get((2021, stage, 'up'), set()) - deg_gene_sets.get((2022, stage, 'up'), set())
    only22_up = deg_gene_sets.get((2022, stage, 'up'), set()) - deg_gene_sets.get((2021, stage, 'up'), set())

    shared_down = deg_gene_sets.get((2021, stage, 'down'), set()) & deg_gene_sets.get((2022, stage, 'down'), set())
    only21_down = deg_gene_sets.get((2021, stage, 'down'), set()) - deg_gene_sets.get((2022, stage, 'down'), set())
    only22_down = deg_gene_sets.get((2022, stage, 'down'), set()) - deg_gene_sets.get((2021, stage, 'down'), set())

    # Enrichment for 2021-only vs 2022-only DEGs
    terms_21 = set()
    terms_22 = set()

    for gene_set, term_set, label in [(only21_up | only21_down, terms_21, '2021'),
                                       (only22_up | only22_down, terms_22, '2022')]:
        if len(gene_set) >= 10:
            enr = run_full_enrichment(list(gene_set))
            for cat in ['GO_BP', 'KEGG']:
                df = enr.get(cat, pd.DataFrame())
                if len(df) > 0:
                    term_set.update(df[df['padj'] < 0.05]['term_name'].head(10).tolist())

    # Display as bar plot
    shared_total = len(shared_up) + len(shared_down)
    only21_total = len(only21_up) + len(only21_down)
    only22_total = len(only22_up) + len(only22_down)

    bars = ax.bar(['Shared', '2021\nonly', '2022\nonly'],
                  [shared_total, only21_total, only22_total],
                  color=['#2ecc71', '#3498db', '#e67e22'], alpha=0.85, edgecolor='black', linewidth=0.5)

    ax.set_ylabel("Number of DEGs")
    ax.set_title(f"{stage}\nSeason-specific DEGs", fontsize=15)

    for b, v in zip(bars, [shared_total, only21_total, only22_total]):
        if v > 0:
            ax.text(b.get_x() + b.get_width()/2, v + max(shared_total, only21_total, only22_total)*0.03,
                    str(v), ha='center', fontsize=12)

fig.suptitle("GRBV Response: Stage Transitions & Season-Specific Patterns",
             fontsize=21, y=1.02)
plt.savefig('18_seasonal_maturity_transitions_PCA.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 18_seasonal_maturity_transitions_PCA.png")


# ============================================================
# PLOT 5: DETAILED PATHWAY COMPARISON (Focused dotplots)
# ============================================================

fig = plt.figure(figsize=(24, 22))
gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.5, wspace=0.4)

# For each stage, show a dotplot comparing 2021 vs 2022 enrichment
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[idx, 0])

    # Collect all terms from both seasons
    all_terms = {}
    for season in seasons:
        for direction in ['up', 'down']:
            enr = enrichments.get((season, stage, direction), {})
            for cat in ['GO_BP', 'GO_MF', 'KEGG']:
                df = enr.get(cat, pd.DataFrame())
                if len(df) > 0:
                    for _, row in df[df['padj'] < 0.1].head(6).iterrows():
                        key = f"[{cat[:2]}] {row['term_name']}"
                        if key not in all_terms:
                            all_terms[key] = {}
                        sign = 1 if direction == 'up' else -1
                        all_terms[key][f"{season}"] = sign * (-np.log10(max(row['padj'], 1e-30)))

    if all_terms:
        hm = pd.DataFrame(all_terms).T.fillna(0)
        hm = hm.reindex(columns=['2021', '2022'])
        hm = hm.dropna(axis=1, how='all')
        if len(hm) > 20:
            hm = hm.loc[hm.var(axis=1).nlargest(20).index]

        if len(hm) > 0:
            sns.heatmap(hm, ax=ax, cmap='RdBu_r', center=0,
                        yticklabels=[textwrap.fill(str(t)[:55], 40) for t in hm.index],
                        xticklabels=True,
                        cbar_kws={'label': 'Signed -log10(padj)', 'shrink': 0.7},
                        linewidths=0.5, annot=True, fmt='.1f', annot_kws={'fontsize': 7})
            ax.tick_params(axis='y', labelsize=6.5)

    ax.set_title(f"{stage}: Pathway Enrichment 2021 vs 2022\n(+Red=Up in RB+, -Blue=Down in RB+)", fontsize=15)

    # Right column: season-specific pathway enrichment
    ax2 = fig.add_subplot(gs[idx, 1])

    # Find pathways that differ between seasons
    season_diff_terms = {}
    for season in seasons:
        for direction in ['up', 'down']:
            enr = enrichments.get((season, stage, direction), {})
            for cat in ['GO_BP', 'KEGG']:
                df = enr.get(cat, pd.DataFrame())
                if len(df) > 0:
                    for _, row in df[df['padj'] < 0.05].head(8).iterrows():
                        key = row['term_name']
                        if key not in season_diff_terms:
                            season_diff_terms[key] = {'2021': 0, '2022': 0}
                        sign = 1 if direction == 'up' else -1
                        season_diff_terms[key][str(season)] = sign * row['fold_enrichment']

    if season_diff_terms:
        diff_df = pd.DataFrame(season_diff_terms).T
        # Sort by difference magnitude
        diff_df['delta'] = (diff_df['2022'] - diff_df['2021']).abs()
        diff_df = diff_df.nlargest(15, 'delta')

        y = np.arange(len(diff_df))
        ax2.barh(y - 0.15, diff_df['2021'].values, height=0.3, color='#3498db', alpha=0.85, label='2021')
        ax2.barh(y + 0.15, diff_df['2022'].values, height=0.3, color='#e67e22', alpha=0.85, label='2022')
        ax2.set_yticks(y)
        ax2.set_yticklabels([textwrap.fill(str(t)[:40], 30) for t in diff_df.index], fontsize=10)
        ax2.set_xlabel("Signed Fold Enrichment\n(+Up in RB+, -Down in RB+)")
        ax2.set_title(f"{stage}: Season Comparison\n(Fold Enrichment)", fontsize=15)
        ax2.legend(fontsize=12)
        ax2.axvline(0, color='black', linewidth=0.5)
        ax2.invert_yaxis()
    else:
        ax2.text(0.5, 0.5, 'No terms to compare', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title(f"{stage}: Season Comparison", fontsize=15)

fig.suptitle("Stage-by-Stage Pathway Comparison: 2021 vs 2022\n"
             "How does the pathway signature of GRBV change across seasons?",
             fontsize=21, y=1.01)
plt.savefig('19_seasonal_maturity_pathway_detail.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 19_seasonal_maturity_pathway_detail.png")


# ============================================================
# COMPREHENSIVE PRINTOUT
# ============================================================
print("\n" + "=" * 70)
print("COMPREHENSIVE SUMMARY: SEASONAL MATURITY PROGRESSION OF GRBV")
print("=" * 70)

print("""
╔══════════════════════════════════════════════════════════════════════╗
║  HOW GRBV-INDUCED CHANGES PROGRESS THROUGH MATURITY STAGES        ║
║  IN SEASONAL SETTINGS (2021 vs 2022)                               ║
╚══════════════════════════════════════════════════════════════════════╝
""")

for stage in stages:
    print(f"\n{'─'*60}")
    print(f"  ■ {stage.upper()}")
    print(f"{'─'*60}")

    for season in seasons:
        dc = deg_counts.get((season, stage), {})
        print(f"\n    Season {season}:")
        print(f"      Total DEGs: {dc.get('total', 0)} (UP={dc.get('up', 0)}, DOWN={dc.get('down', 0)})")

        for direction in ['up', 'down']:
            enr = enrichments.get((season, stage, direction), {})
            combined = get_combined_top(enr, n_per_cat=3)
            if len(combined) > 0:
                sig = combined[combined['padj'] < 0.05]
                if len(sig) > 0:
                    terms = '; '.join(sig['term_name'].head(4).tolist())
                    print(f"      {direction.upper()} pathways: {terms}")

    # Cross-season comparison
    fc_info = fc_correlations.get(stage, {})
    if fc_info:
        print(f"\n    Cross-season FC correlation: r = {fc_info['r']:.4f}")

    up_21 = deg_gene_sets.get((2021, stage, 'up'), set())
    up_22 = deg_gene_sets.get((2022, stage, 'up'), set())
    down_21 = deg_gene_sets.get((2021, stage, 'down'), set())
    down_22 = deg_gene_sets.get((2022, stage, 'down'), set())

    shared_total = len(up_21 & up_22) + len(down_21 & down_22)
    total_21 = len(up_21) + len(down_21)
    total_22 = len(up_22) + len(down_22)

    if total_21 > 0 and total_22 > 0:
        jaccard = shared_total / len((up_21 | down_21) | (up_22 | down_22)) if len((up_21 | down_21) | (up_22 | down_22)) > 0 else 0
        print(f"    DEG overlap: {shared_total} shared / 2021={total_21} / 2022={total_22} (Jaccard={jaccard:.3f})")

print("\n\n" + "=" * 70)
print("STAGE TRANSITION SUMMARY (How GRBV effect changes with ripening)")
print("=" * 70)

for season in seasons:
    print(f"\n  Season {season}:")
    for trans in transitions:
        tc = transition_deg_counts.get((season, trans), {})
        print(f"    {trans}: Intensifying={tc.get('intensifying', 0)}, Weakening={tc.get('weakening', 0)}")


print("\n\n" + "=" * 70)
print("GENERATED PLOTS:")
print("  15_seasonal_maturity_DEG_comparison.png     — DEG counts, overlaps, FC correlations")
print("  16_seasonal_maturity_pathway_enrichment.png  — Full pathway enrichment (4×4 grid)")
print("  17_seasonal_maturity_pathway_landscape.png   — Pathway heatmaps (2021 vs 2022)")
print("  18_seasonal_maturity_transitions_PCA.png     — PCA, stage transitions, unique DEGs")
print("  19_seasonal_maturity_pathway_detail.png      — Stage-by-stage pathway comparison")
print("=" * 70)
