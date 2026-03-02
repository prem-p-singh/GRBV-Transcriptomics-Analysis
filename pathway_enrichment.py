#!/usr/bin/env python3
"""
Comprehensive GO + KEGG pathway enrichment analysis for GRBV study.
Uses the org.Vvinifera.eg.db SQLite database for annotations.
Answers all 6 research questions at the pathway level.
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
    'figure.dpi': 150, 'font.size': 9, 'axes.titlesize': 11,
    'axes.labelsize': 10, 'figure.facecolor': 'white'
})

# ============================================================
# 1. LOAD DATA
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

# Infection age
meta['Actual_Infection_Age'] = meta.apply(
    lambda row: int(row['Sampling_Season']) - int(row['Years_Since_Infection'])
    if row['Years_Since_Infection'] > 0 else 0, axis=1
)

print(f"Filtered genes: {counts_filtered.shape[0]}")

# ============================================================
# 2. LOAD ANNOTATIONS FROM orgdb
# ============================================================
db_path = 'org.Vvinifera.eg.db/inst/extdata/org.Vvinifera.eg.sqlite'
conn = sqlite3.connect(db_path)

# Build symbol -> _id mapping
symbol_map = pd.read_sql("SELECT gi._id, gi.SYMBOL, g.GID FROM gene_info gi JOIN genes g ON gi._id = g._id", conn)
symbol_to_id = dict(zip(symbol_map['SYMBOL'], symbol_map['_id']))

# GO annotations
go_df = pd.read_sql("SELECT g.SYMBOL, go.GO, go.ONTOLOGY FROM go JOIN gene_info g ON go._id = g._id", conn)
print(f"GO annotation rows: {len(go_df)}")

# KEGG annotations
kegg_df = pd.read_sql("SELECT g.SYMBOL, k.KEGG FROM kegg k JOIN gene_info g ON k._id = g._id", conn)
print(f"KEGG annotation rows: {len(kegg_df)}")

# Get KEGG pathway names via REST API (or build from IDs)
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
    print(f"KEGG names fetch failed: {e}, using IDs as names")
    kegg_names = {}

# Get GO term names
from goatools.obo_parser import GODag
import os

obo_path = 'go-basic.obo'
if not os.path.exists(obo_path):
    print("Downloading GO ontology...")
    urllib.request.urlretrieve("http://purl.obolibrary.org/obo/go-basic.obo", obo_path)

godag = GODag(obo_path)
print(f"GO terms loaded: {len(godag)}")

# Build annotation dictionaries
# gene -> set of GO terms (by ontology)
gene2go_bp = defaultdict(set)
gene2go_mf = defaultdict(set)
gene2go_cc = defaultdict(set)
gene2kegg = defaultdict(set)

# GO term -> set of genes
go2genes_bp = defaultdict(set)
go2genes_mf = defaultdict(set)
go2genes_cc = defaultdict(set)
kegg2genes = defaultdict(set)

for _, row in go_df.iterrows():
    gene, term, ont = row['SYMBOL'], row['GO'], row['ONTOLOGY']
    if ont == 'BP':
        gene2go_bp[gene].add(term)
        go2genes_bp[term].add(gene)
    elif ont == 'MF':
        gene2go_mf[gene].add(term)
        go2genes_mf[term].add(gene)
    elif ont == 'CC':
        gene2go_cc[gene].add(term)
        go2genes_cc[term].add(gene)

for _, row in kegg_df.iterrows():
    gene, pathway = row['SYMBOL'], row['KEGG']
    gene2kegg[gene].add(pathway)
    kegg2genes[pathway].add(gene)

# Background: all genes in filtered counts that have annotations
all_genes = set(counts_filtered.index)
bg_go_bp = all_genes & set(gene2go_bp.keys())
bg_go_mf = all_genes & set(gene2go_mf.keys())
bg_go_cc = all_genes & set(gene2go_cc.keys())
bg_kegg = all_genes & set(gene2kegg.keys())

print(f"\nBackground genes with GO-BP: {len(bg_go_bp)}")
print(f"Background genes with GO-MF: {len(bg_go_mf)}")
print(f"Background genes with GO-CC: {len(bg_go_cc)}")
print(f"Background genes with KEGG: {len(bg_kegg)}")

conn.close()

# ============================================================
# 3. ENRICHMENT FUNCTIONS
# ============================================================

def fisher_enrichment(query_genes, term2genes, background_genes, min_genes=3):
    """Fisher's exact test for enrichment."""
    query = set(query_genes) & background_genes
    n_query = len(query)
    n_bg = len(background_genes)

    results = []
    for term, term_genes in term2genes.items():
        term_in_bg = term_genes & background_genes
        if len(term_in_bg) < min_genes:
            continue

        a = len(query & term_in_bg)       # query AND term
        b = len(query - term_in_bg)        # query NOT term
        c = len(term_in_bg - query)        # term NOT query
        d = n_bg - len(query | term_in_bg) # neither

        if a < 2:
            continue

        odds, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        results.append({
            'term': term,
            'overlap': a,
            'term_size': len(term_in_bg),
            'query_size': n_query,
            'bg_size': n_bg,
            'fold_enrichment': (a / n_query) / (len(term_in_bg) / n_bg) if n_query > 0 and len(term_in_bg) > 0 else 0,
            'pvalue': pval
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    df['padj'] = np.minimum(df['pvalue'] * len(df) / rankdata(df['pvalue']), 1.0)
    return df.sort_values('pvalue')


def run_full_enrichment(query_genes, label=""):
    """Run GO (BP, MF, CC) + KEGG enrichment for a gene set."""
    query = set(query_genes) & all_genes

    res = {}
    for name, t2g, bg in [('GO_BP', go2genes_bp, bg_go_bp),
                            ('GO_MF', go2genes_mf, bg_go_mf),
                            ('GO_CC', go2genes_cc, bg_go_cc),
                            ('KEGG', kegg2genes, bg_kegg)]:
        df = fisher_enrichment(query, t2g, bg)
        if len(df) > 0:
            # Add term names
            if name.startswith('GO'):
                df['term_name'] = df['term'].apply(
                    lambda x: godag[x].name if x in godag else x)
            else:
                df['term_name'] = df['term'].apply(
                    lambda x: kegg_names.get(x, x))
            df['category'] = name
        res[name] = df

    return res


def format_top_terms(enrichment_dict, n=8, padj_thresh=0.05):
    """Format top enriched terms for printing."""
    lines = []
    for cat in ['GO_BP', 'GO_MF', 'GO_CC', 'KEGG']:
        df = enrichment_dict.get(cat, pd.DataFrame())
        if len(df) == 0:
            lines.append(f"  {cat}: No significant enrichment")
            continue
        sig = df[df['padj'] < padj_thresh].head(n)
        if len(sig) == 0:
            lines.append(f"  {cat}: No terms at padj < {padj_thresh}")
            continue
        lines.append(f"  {cat} ({len(df[df['padj'] < padj_thresh])} significant terms):")
        for _, row in sig.iterrows():
            name = row['term_name'][:60]
            lines.append(f"    {name:60s} FE={row['fold_enrichment']:.1f}  padj={row['padj']:.2e}  [{row['overlap']}/{row['term_size']}]")
    return '\n'.join(lines)


def get_combined_top(enrichment_dict, n_per_cat=5, padj_thresh=0.05):
    """Get combined top terms across categories."""
    combined = []
    for cat in ['GO_BP', 'GO_MF', 'GO_CC', 'KEGG']:
        df = enrichment_dict.get(cat, pd.DataFrame())
        if len(df) > 0:
            sig = df[df['padj'] < padj_thresh].head(n_per_cat)
            combined.append(sig)
    if combined:
        return pd.concat(combined, ignore_index=True)
    return pd.DataFrame()


# ============================================================
# 4. DE HELPER FUNCTIONS
# ============================================================

def run_de(group1_samples, group2_samples):
    """T-test DE with BH correction."""
    g1 = log_cpm[group1_samples]
    g2 = log_cpm[group2_samples]
    results = []
    for gene in log_cpm.index:
        v1 = g1.loc[gene].values
        v2 = g2.loc[gene].values
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
    """Extract DEGs."""
    sig = de_df[(de_df['padj'] < padj) & (de_df['log2FC'].abs() >= fc)]
    return sig


def get_up_down(de_df, fc=1.0, padj=0.05):
    """Get up and down gene lists."""
    sig = get_degs(de_df, fc, padj)
    up = sig[sig['log2FC'] > 0]['gene'].tolist()
    down = sig[sig['log2FC'] < 0]['gene'].tolist()
    return up, down


def get_samples(meta_df, **kwargs):
    mask = pd.Series(True, index=meta_df.index)
    for col, val in kwargs.items():
        if isinstance(val, list):
            mask &= meta_df[col].isin(val)
        else:
            mask &= meta_df[col] == val
    return list(meta_df[mask].index)


# ============================================================
# PLOTTING HELPER
# ============================================================

def plot_enrichment_bars(enrichment_dict, ax, title="", n=10, padj_thresh=0.05):
    """Plot enrichment bar chart on given axis."""
    combined = get_combined_top(enrichment_dict, n_per_cat=4, padj_thresh=padj_thresh)
    if len(combined) == 0:
        ax.text(0.5, 0.5, 'No significant\nenrichment', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return

    combined = combined.head(n)
    cat_colors = {'GO_BP': '#e74c3c', 'GO_MF': '#3498db', 'GO_CC': '#2ecc71', 'KEGG': '#9b59b6'}

    labels = [textwrap.fill(str(row['term_name'])[:50], 35) for _, row in combined.iterrows()]
    values = [-np.log10(combined['padj'].clip(lower=1e-30))]
    colors = [cat_colors.get(row['category'], 'grey') for _, row in combined.iterrows()]

    y = np.arange(len(combined))
    bars = ax.barh(y, -np.log10(combined['padj'].clip(lower=1e-30).values), color=colors, alpha=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('-log10(padj)')
    ax.set_title(title, fontsize=10)
    ax.invert_yaxis()
    ax.axvline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)

    # Add overlap info
    for i, (_, row) in enumerate(combined.iterrows()):
        ax.text(0.2, i, f" {int(row['overlap'])}/{int(row['term_size'])}", va='center', fontsize=6, color='white', fontweight='bold')


def plot_enrichment_dotplot(enrichment_dict, ax, title="", n=12, padj_thresh=0.05):
    """Dot plot for enrichment results."""
    combined = get_combined_top(enrichment_dict, n_per_cat=5, padj_thresh=padj_thresh)
    if len(combined) == 0:
        ax.text(0.5, 0.5, 'No significant\nenrichment', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return

    combined = combined.head(n)
    cat_colors = {'GO_BP': '#e74c3c', 'GO_MF': '#3498db', 'GO_CC': '#2ecc71', 'KEGG': '#9b59b6'}

    y = np.arange(len(combined))
    colors = [cat_colors.get(row['category'], 'grey') for _, row in combined.iterrows()]
    sizes = combined['overlap'].values * 8

    sc = ax.scatter(combined['fold_enrichment'].values, y, s=sizes,
                     c=-np.log10(combined['padj'].clip(lower=1e-30).values),
                     cmap='YlOrRd', edgecolors='black', linewidth=0.5, zorder=3)

    labels = [textwrap.fill(str(row['term_name'])[:45], 30) for _, row in combined.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('Fold Enrichment')
    ax.set_title(title, fontsize=10)
    ax.invert_yaxis()
    plt.colorbar(sc, ax=ax, label='-log10(padj)', shrink=0.7)


# ============================================================
# Q1: INFECTION YEARS — PATHWAY ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("Q1: PATHWAY ENRICHMENT — INFECTION YEARS")
print("=" * 70)

control_samples = get_samples(meta, Infection_Group='RB(-)')
ysi_values = sorted([y for y in meta['Years_Since_Infection'].unique() if y > 0])

q1_enrichments = {}
for ysi in ysi_values:
    infected = get_samples(meta, Infection_Group='RB(+)', Years_Since_Infection=ysi)
    de = run_de(infected, control_samples)
    up, down = get_up_down(de)

    print(f"\n--- YSI={ysi} (up={len(up)}, down={len(down)}) ---")

    if len(up) >= 5:
        enr_up = run_full_enrichment(up, f"YSI={ysi} UP")
        print(f"UP-regulated:\n{format_top_terms(enr_up)}")
        q1_enrichments[(ysi, 'up')] = enr_up
    else:
        print(f"UP: Too few DEGs ({len(up)})")
        q1_enrichments[(ysi, 'up')] = {}

    if len(down) >= 5:
        enr_down = run_full_enrichment(down, f"YSI={ysi} DOWN")
        print(f"DOWN-regulated:\n{format_top_terms(enr_down)}")
        q1_enrichments[(ysi, 'down')] = enr_down
    else:
        print(f"DOWN: Too few DEGs ({len(down)})")
        q1_enrichments[(ysi, 'down')] = {}

# Pathway trend across infection years
print("\n\nPathway trends across infection cohorts:")
# Collect all significant GO-BP terms across years
all_bp_terms = set()
for ysi in ysi_values:
    for direction in ['up', 'down']:
        enr = q1_enrichments.get((ysi, direction), {})
        bp = enr.get('GO_BP', pd.DataFrame())
        if len(bp) > 0:
            sig = bp[bp['padj'] < 0.1].head(15)
            all_bp_terms.update(sig['term'].tolist())

# Build matrix: YSI x term -> -log10(padj) * sign
bp_matrix = pd.DataFrame(0.0, index=list(all_bp_terms), columns=[f"{y}_up" for y in ysi_values] + [f"{y}_down" for y in ysi_values])
for ysi in ysi_values:
    for direction in ['up', 'down']:
        enr = q1_enrichments.get((ysi, direction), {})
        bp = enr.get('GO_BP', pd.DataFrame())
        if len(bp) > 0:
            for _, row in bp.iterrows():
                if row['term'] in all_bp_terms:
                    val = -np.log10(max(row['padj'], 1e-30))
                    if direction == 'down':
                        val = -val
                    bp_matrix.loc[row['term'], f"{ysi}_{direction}"] = val

# Plot Q1
fig = plt.figure(figsize=(20, 16))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.5)

# Enrichment bars for key infection years
for idx, ysi in enumerate([2016, 2019, 2021]):
    ax = fig.add_subplot(gs[0, idx])
    enr = q1_enrichments.get((ysi, 'up'), {})
    plot_enrichment_bars(enr, ax, f"YSI={ysi} Upregulated", n=8)

# Enrichment for down across years
for idx, ysi in enumerate([2016, 2019, 2021]):
    ax = fig.add_subplot(gs[1, idx])
    enr = q1_enrichments.get((ysi, 'down'), {})
    plot_enrichment_bars(enr, ax, f"YSI={ysi} Downregulated", n=8)

fig.suptitle("Q1: Pathway Enrichment by Infection Year (vs Controls)", fontsize=14, y=1.01)
plt.savefig('08_Q1_pathway_enrichment.png', bbox_inches='tight')
plt.close()
print("\nSaved: 08_Q1_pathway_enrichment.png")


# ============================================================
# Q2: MATURITY STAGES — PATHWAY ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("Q2: PATHWAY ENRICHMENT — MATURITY STAGES")
print("=" * 70)

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
q2_enrichments = {}

for stage in stages:
    pos = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=stage)
    neg = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=stage)
    de = run_de(pos, neg)
    up, down = get_up_down(de)

    print(f"\n--- {stage} (up={len(up)}, down={len(down)}) ---")

    for direction, genes in [('up', up), ('down', down)]:
        if len(genes) >= 5:
            enr = run_full_enrichment(genes)
            print(f"{direction.upper()}:\n{format_top_terms(enr)}")
            q2_enrichments[(stage, direction)] = enr
        else:
            print(f"{direction.upper()}: Too few DEGs")
            q2_enrichments[(stage, direction)] = {}

# Stage transition enrichment within RB(+)
print("\n\nStage transition pathways in RB(+):")
q2_transitions = {}
for i in range(len(stages) - 1):
    s1, s2 = stages[i], stages[i+1]
    samp1 = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=s1)
    samp2 = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=s2)
    de = run_de(samp2, samp1)
    up, down = get_up_down(de)
    print(f"\n--- {s1} -> {s2} (up={len(up)}, down={len(down)}) ---")
    for direction, genes in [('up', up), ('down', down)]:
        if len(genes) >= 10:
            enr = run_full_enrichment(genes)
            print(f"{direction.upper()}:\n{format_top_terms(enr)}")
            q2_transitions[(f"{s1}->{s2}", direction)] = enr

# Plot Q2
fig = plt.figure(figsize=(22, 20))
gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.5, wspace=0.6)

# Row 1: Upregulated per stage
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[0, idx])
    enr = q2_enrichments.get((stage, 'up'), {})
    plot_enrichment_bars(enr, ax, f"{stage}\nUP in RB(+)", n=8)

# Row 2: Downregulated per stage
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[1, idx])
    enr = q2_enrichments.get((stage, 'down'), {})
    plot_enrichment_bars(enr, ax, f"{stage}\nDOWN in RB(+)", n=8)

# Row 3: Build pathway heatmap across stages
# Collect significant terms
all_terms_q2 = {}
for stage in stages:
    for direction in ['up', 'down']:
        enr = q2_enrichments.get((stage, direction), {})
        for cat in ['GO_BP', 'KEGG']:
            df = enr.get(cat, pd.DataFrame())
            if len(df) > 0:
                for _, row in df[df['padj'] < 0.05].head(8).iterrows():
                    key = row['term_name']
                    if key not in all_terms_q2:
                        all_terms_q2[key] = {}
                    sign = 1 if direction == 'up' else -1
                    all_terms_q2[key][stage] = sign * (-np.log10(max(row['padj'], 1e-30)))

if all_terms_q2:
    hm_df = pd.DataFrame(all_terms_q2).T.fillna(0)
    hm_df = hm_df.reindex(columns=[s for s in stages if s in hm_df.columns])
    # Select most variable terms
    hm_var = hm_df.var(axis=1).nlargest(25).index
    hm_df = hm_df.loc[hm_var]

    ax = fig.add_subplot(gs[2, :])
    sns.heatmap(hm_df, ax=ax, cmap='RdBu_r', center=0,
                yticklabels=[textwrap.fill(str(t)[:55], 50) for t in hm_df.index],
                xticklabels=True, cbar_kws={'label': 'Signed -log10(padj)\n(+Up in RB+, -Down in RB+)'},
                linewidths=0.5)
    ax.set_title("Pathway enrichment landscape across maturity stages")
    ax.tick_params(axis='y', labelsize=7)

fig.suptitle("Q2: Pathway Enrichment by Maturity Stage (RB+ vs RB-)", fontsize=14, y=1.01)
plt.savefig('09_Q2_pathway_enrichment.png', bbox_inches='tight')
plt.close()
print("\nSaved: 09_Q2_pathway_enrichment.png")


# ============================================================
# Q3: SEASONAL EFFECTS — PATHWAY ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("Q3: PATHWAY ENRICHMENT — SEASONAL EFFECTS")
print("=" * 70)

q3_enrichments = {}

# Season-specific infection response
for season in [2021, 2022]:
    pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season)
    neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season)
    de = run_de(pos, neg)
    up, down = get_up_down(de)
    print(f"\n--- Season {season}: Infection response (up={len(up)}, down={len(down)}) ---")
    for direction, genes in [('up', up), ('down', down)]:
        if len(genes) >= 5:
            enr = run_full_enrichment(genes)
            print(f"{direction.upper()}:\n{format_top_terms(enr)}")
            q3_enrichments[(season, 'infection', direction)] = enr
        else:
            q3_enrichments[(season, 'infection', direction)] = {}

# Seasonal change within RB(+) and RB(-)
for grp in ['RB(+)', 'RB(-)']:
    s21 = get_samples(meta, Infection_Group=grp, Sampling_Season=2021)
    s22 = get_samples(meta, Infection_Group=grp, Sampling_Season=2022)
    de = run_de(s22, s21)
    up, down = get_up_down(de)
    grp_label = 'pos' if grp == 'RB(+)' else 'neg'
    print(f"\n--- {grp}: 2022 vs 2021 (up={len(up)}, down={len(down)}) ---")
    for direction, genes in [('up', up), ('down', down)]:
        if len(genes) >= 5:
            enr = run_full_enrichment(genes)
            print(f"{direction.upper()}:\n{format_top_terms(enr)}")
            q3_enrichments[('seasonal', grp_label, direction)] = enr
        else:
            q3_enrichments[('seasonal', grp_label, direction)] = {}

# Season x Stage interaction
print("\n\nSeason x Stage pathway enrichment:")
q3_sxs = {}
for season in [2021, 2022]:
    for stage in stages:
        pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season, Maturity_stage=stage)
        neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season, Maturity_stage=stage)
        if len(pos) < 3 or len(neg) < 3:
            continue
        de = run_de(pos, neg)
        up, down = get_up_down(de)
        print(f"\n--- {season}-{stage} (up={len(up)}, down={len(down)}) ---")
        for direction, genes in [('up', up), ('down', down)]:
            if len(genes) >= 10:
                enr = run_full_enrichment(genes)
                print(f"{direction.upper()}:\n{format_top_terms(enr, n=5)}")
                q3_sxs[(season, stage, direction)] = enr

# Plot Q3
fig = plt.figure(figsize=(22, 16))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.55)

# Infection response per season
ax1 = fig.add_subplot(gs[0, 0])
plot_enrichment_bars(q3_enrichments.get((2021, 'infection', 'up'), {}), ax1, "2021 UP (RB+ vs RB-)", n=10)

ax2 = fig.add_subplot(gs[0, 1])
plot_enrichment_bars(q3_enrichments.get((2022, 'infection', 'up'), {}), ax2, "2022 UP (RB+ vs RB-)", n=10)

ax3 = fig.add_subplot(gs[0, 2])
plot_enrichment_bars(q3_enrichments.get((2021, 'infection', 'down'), {}), ax3, "2021 DOWN (RB+ vs RB-)", n=10)

# Season x Stage heatmap
all_terms_q3 = {}
for season in [2021, 2022]:
    for stage in stages:
        for direction in ['up', 'down']:
            enr = q3_sxs.get((season, stage, direction), {})
            for cat in ['GO_BP', 'KEGG']:
                df = enr.get(cat, pd.DataFrame())
                if len(df) > 0:
                    for _, row in df[df['padj'] < 0.05].head(5).iterrows():
                        key = row['term_name']
                        if key not in all_terms_q3:
                            all_terms_q3[key] = {}
                        sign = 1 if direction == 'up' else -1
                        col = f"{season}\n{stage}"
                        all_terms_q3[key][col] = sign * (-np.log10(max(row['padj'], 1e-30)))

if all_terms_q3:
    hm_q3 = pd.DataFrame(all_terms_q3).T.fillna(0)
    col_order = [f"{s}\n{st}" for s in [2021, 2022] for st in stages]
    hm_q3 = hm_q3.reindex(columns=[c for c in col_order if c in hm_q3.columns])
    hm_var = hm_q3.var(axis=1).nlargest(20).index
    hm_q3 = hm_q3.loc[hm_var]

    ax = fig.add_subplot(gs[1, :])
    sns.heatmap(hm_q3, ax=ax, cmap='RdBu_r', center=0,
                yticklabels=[textwrap.fill(str(t)[:50], 40) for t in hm_q3.index],
                xticklabels=True, cbar_kws={'label': 'Signed -log10(padj)'},
                linewidths=0.5)
    ax.set_title("Season × Stage pathway enrichment landscape")
    ax.tick_params(axis='y', labelsize=7)

fig.suptitle("Q3: Seasonal Effects — Pathway Enrichment", fontsize=14, y=1.01)
plt.savefig('10_Q3_pathway_enrichment.png', bbox_inches='tight')
plt.close()
print("\nSaved: 10_Q3_pathway_enrichment.png")


# ============================================================
# Q4: ENVIRONMENTAL EFFECTS — PATHWAY ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("Q4: PATHWAY ENRICHMENT — ENVIRONMENTAL EFFECTS")
print("=" * 70)

# Seasonal change within each infection group at each stage
q4_env = {}
for grp in ['RB(+)', 'RB(-)']:
    for stage in stages:
        s21 = get_samples(meta, Infection_Group=grp, Sampling_Season=2021, Maturity_stage=stage)
        s22 = get_samples(meta, Infection_Group=grp, Sampling_Season=2022, Maturity_stage=stage)
        if len(s21) < 3 or len(s22) < 3:
            continue
        de = run_de(s22, s21)
        up, down = get_up_down(de)
        grp_label = 'pos' if grp == 'RB(+)' else 'neg'
        print(f"\n--- {grp} {stage}: 2022 vs 2021 (up={len(up)}, down={len(down)}) ---")
        for direction, genes in [('up', up), ('down', down)]:
            if len(genes) >= 10:
                enr = run_full_enrichment(genes)
                print(f"{direction.upper()}:\n{format_top_terms(enr, n=5)}")
                q4_env[(grp_label, stage, direction)] = enr

# Compare env pathways in infected vs uninfected
print("\n\nEnvironmental pathway comparison (RB+ vs RB- seasonal changes):")
for stage in stages:
    print(f"\n  {stage}:")
    for direction in ['up', 'down']:
        pos_enr = q4_env.get(('pos', stage, direction), {})
        neg_enr = q4_env.get(('neg', stage, direction), {})

        pos_terms = set()
        neg_terms = set()
        for cat in ['GO_BP', 'KEGG']:
            df_pos = pos_enr.get(cat, pd.DataFrame())
            df_neg = neg_enr.get(cat, pd.DataFrame())
            if len(df_pos) > 0:
                pos_terms.update(df_pos[df_pos['padj'] < 0.05]['term_name'].tolist())
            if len(df_neg) > 0:
                neg_terms.update(df_neg[df_neg['padj'] < 0.05]['term_name'].tolist())

        shared = pos_terms & neg_terms
        pos_only = pos_terms - neg_terms
        neg_only = neg_terms - pos_terms
        print(f"    {direction}: Shared={len(shared)}, RB(+)-only={len(pos_only)}, RB(-)-only={len(neg_only)}")

# Plot Q4
fig = plt.figure(figsize=(22, 16))
gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.5, wspace=0.6)

# Environmental response pathways per stage (RB+)
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[0, idx])
    enr = q4_env.get(('pos', stage, 'up'), {})
    plot_enrichment_bars(enr, ax, f"RB(+) {stage}\nUP in 2022", n=8)

# Environmental response pathways per stage (RB-)
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[1, idx])
    enr = q4_env.get(('neg', stage, 'up'), {})
    plot_enrichment_bars(enr, ax, f"RB(-) {stage}\nUP in 2022", n=8)

fig.suptitle("Q4: Environmental Effects — Pathway Enrichment (2022 vs 2021)", fontsize=14, y=1.01)
plt.savefig('11_Q4_pathway_enrichment.png', bbox_inches='tight')
plt.close()
print("\nSaved: 11_Q4_pathway_enrichment.png")


# ============================================================
# Q5: UNINFECTED ENVIRONMENT — PATHWAY ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("Q5: PATHWAY ENRICHMENT — UNINFECTED GRAPEVINE ENVIRONMENT")
print("=" * 70)

# Overall seasonal DEGs in RB(-)
neg_21 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2021)
neg_22 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2022)
de_neg = run_de(neg_22, neg_21)
up_neg, down_neg = get_up_down(de_neg)

print(f"\nOverall RB(-) 2022 vs 2021: up={len(up_neg)}, down={len(down_neg)}")
for direction, genes in [('up', up_neg), ('down', down_neg)]:
    if len(genes) >= 5:
        enr = run_full_enrichment(genes)
        print(f"\n{direction.upper()}:\n{format_top_terms(enr)}")

# Per-stage seasonal changes in RB(-)
q5_stage = {}
for stage in stages:
    s21 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2021, Maturity_stage=stage)
    s22 = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=2022, Maturity_stage=stage)
    if len(s21) < 3 or len(s22) < 3:
        continue
    de = run_de(s22, s21)
    up, down = get_up_down(de)
    print(f"\n--- RB(-) {stage}: 2022 vs 2021 (up={len(up)}, down={len(down)}) ---")
    for direction, genes in [('up', up), ('down', down)]:
        if len(genes) >= 10:
            enr = run_full_enrichment(genes)
            print(f"{direction.upper()}:\n{format_top_terms(enr, n=5)}")
            q5_stage[(stage, direction)] = enr

# Developmental program in RB(-)
print("\n\nDevelopmental pathway enrichment in healthy grapevines:")
q5_dev = {}
for i in range(len(stages) - 1):
    s1, s2 = stages[i], stages[i+1]
    samp1 = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=s1)
    samp2 = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=s2)
    de = run_de(samp2, samp1)
    up, down = get_up_down(de)
    print(f"\n--- {s1} -> {s2} (up={len(up)}, down={len(down)}) ---")
    for direction, genes in [('up', up), ('down', down)]:
        if len(genes) >= 10:
            enr = run_full_enrichment(genes)
            print(f"{direction.upper()}:\n{format_top_terms(enr, n=5)}")
            q5_dev[(f"{s1}->{s2}", direction)] = enr

# Plot Q5
fig = plt.figure(figsize=(22, 20))
gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.5, wspace=0.6)

# Row 1: Seasonal changes per stage (RB-)
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[0, idx])
    enr = q5_stage.get((stage, 'up'), {})
    plot_enrichment_bars(enr, ax, f"RB(-) {stage}\nUP in 2022", n=8)

# Row 2: Seasonal down per stage
for idx, stage in enumerate(stages):
    ax = fig.add_subplot(gs[1, idx])
    enr = q5_stage.get((stage, 'down'), {})
    plot_enrichment_bars(enr, ax, f"RB(-) {stage}\nDOWN in 2022", n=8)

# Row 3: Developmental transitions
transitions = ['Pre_veraison->Veraison', 'Veraison->Post_veraison', 'Post_veraison->Harvest']
for idx, trans in enumerate(transitions):
    ax = fig.add_subplot(gs[2, idx])
    enr = q5_dev.get((trans, 'up'), {})
    plot_enrichment_bars(enr, ax, f"RB(-) {trans}\nUP", n=8)

ax = fig.add_subplot(gs[2, 3])
# Summary: developmental transition heatmap
dev_terms = {}
for trans in transitions:
    for direction in ['up', 'down']:
        enr = q5_dev.get((trans, direction), {})
        for cat in ['GO_BP', 'KEGG']:
            df = enr.get(cat, pd.DataFrame())
            if len(df) > 0:
                for _, row in df[df['padj'] < 0.05].head(5).iterrows():
                    key = row['term_name']
                    if key not in dev_terms:
                        dev_terms[key] = {}
                    sign = 1 if direction == 'up' else -1
                    dev_terms[key][trans] = sign * (-np.log10(max(row['padj'], 1e-30)))

if dev_terms:
    hm_dev = pd.DataFrame(dev_terms).T.fillna(0)
    hm_dev = hm_dev.reindex(columns=[t for t in transitions if t in hm_dev.columns])
    hm_var = hm_dev.var(axis=1).nlargest(15).index
    hm_dev = hm_dev.loc[hm_var]
    sns.heatmap(hm_dev, ax=ax, cmap='RdBu_r', center=0,
                yticklabels=[textwrap.fill(str(t)[:40], 30) for t in hm_dev.index],
                cbar_kws={'label': 'Signed\n-log10(padj)'}, linewidths=0.5)
    ax.set_title("Developmental pathways\n(healthy grapevines)")
    ax.tick_params(axis='y', labelsize=6)

fig.suptitle("Q5: Environmental Changes in Uninfected Grapevines — Pathway Enrichment", fontsize=14, y=1.01)
plt.savefig('12_Q5_pathway_enrichment.png', bbox_inches='tight')
plt.close()
print("\nSaved: 12_Q5_pathway_enrichment.png")


# ============================================================
# Q6: TRANSCRIPTOME TRAJECTORY — PATHWAY ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("Q6: PATHWAY ENRICHMENT — TRANSCRIPTOME TRAJECTORY")
print("=" * 70)

infected_meta = meta[meta['Infection_Group'] == 'RB(+)']
age_groups = sorted(infected_meta['Actual_Infection_Age'].unique())
all_neg = get_samples(meta, Infection_Group='RB(-)')

q6_enrichments = {}
for age in age_groups:
    infected = list(infected_meta[infected_meta['Actual_Infection_Age'] == age].index)
    de = run_de(infected, all_neg)
    up, down = get_up_down(de)
    print(f"\n--- Age {age}yr vs controls (up={len(up)}, down={len(down)}) ---")
    for direction, genes in [('up', up), ('down', down)]:
        if len(genes) >= 5:
            enr = run_full_enrichment(genes)
            print(f"{direction.upper()}:\n{format_top_terms(enr, n=5)}")
            q6_enrichments[(age, direction)] = enr
        else:
            q6_enrichments[(age, direction)] = {}

# Trajectory heatmap
print("\n\nPathway trajectory over infection course:")
all_terms_q6 = {}
for age in age_groups:
    for direction in ['up', 'down']:
        enr = q6_enrichments.get((age, direction), {})
        for cat in ['GO_BP', 'GO_MF', 'KEGG']:
            df = enr.get(cat, pd.DataFrame())
            if len(df) > 0:
                for _, row in df[df['padj'] < 0.1].head(8).iterrows():
                    key = f"[{cat}] {row['term_name']}"
                    if key not in all_terms_q6:
                        all_terms_q6[key] = {}
                    sign = 1 if direction == 'up' else -1
                    all_terms_q6[key][f"{age}yr"] = sign * (-np.log10(max(row['padj'], 1e-30)))

# Consecutive age transitions
print("\n\nConsecutive transition pathways:")
q6_trans = {}
for i in range(len(age_groups) - 1):
    a1, a2 = age_groups[i], age_groups[i+1]
    s1 = list(infected_meta[infected_meta['Actual_Infection_Age'] == a1].index)
    s2 = list(infected_meta[infected_meta['Actual_Infection_Age'] == a2].index)
    de = run_de(s2, s1)
    up, down = get_up_down(de)
    if len(up) + len(down) >= 5:
        print(f"\n--- {a1}yr -> {a2}yr (up={len(up)}, down={len(down)}) ---")
        for direction, genes in [('up', up), ('down', down)]:
            if len(genes) >= 5:
                enr = run_full_enrichment(genes)
                print(f"{direction.upper()}:\n{format_top_terms(enr, n=5)}")
                q6_trans[(f"{a1}->{a2}", direction)] = enr

# Plot Q6
fig = plt.figure(figsize=(22, 18))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.5, wspace=0.55)

# Row 1: Enrichment for key ages (up)
for idx, age in enumerate([0, 2, 5]):
    if age in age_groups:
        ax = fig.add_subplot(gs[0, idx])
        enr = q6_enrichments.get((age, 'up'), {})
        plot_enrichment_bars(enr, ax, f"Age {age}yr UP (vs controls)", n=8)

# Row 2: Down at key ages
for idx, age in enumerate([0, 2, 5]):
    if age in age_groups:
        ax = fig.add_subplot(gs[1, idx])
        enr = q6_enrichments.get((age, 'down'), {})
        plot_enrichment_bars(enr, ax, f"Age {age}yr DOWN (vs controls)", n=8)

# Row 3: Pathway trajectory heatmap
if all_terms_q6:
    hm_q6 = pd.DataFrame(all_terms_q6).T.fillna(0)
    col_order = [f"{a}yr" for a in age_groups]
    hm_q6 = hm_q6.reindex(columns=[c for c in col_order if c in hm_q6.columns])
    hm_var = hm_q6.var(axis=1).nlargest(25).index
    hm_q6 = hm_q6.loc[hm_var]

    ax = fig.add_subplot(gs[2, :])
    sns.heatmap(hm_q6, ax=ax, cmap='RdBu_r', center=0,
                yticklabels=[textwrap.fill(str(t)[:55], 45) for t in hm_q6.index],
                xticklabels=True,
                cbar_kws={'label': 'Signed -log10(padj)\n(+Up, -Down vs controls)'},
                linewidths=0.5)
    ax.set_title("Pathway enrichment trajectory over infection course")
    ax.tick_params(axis='y', labelsize=7)
    ax.set_xlabel("Infection Age")

fig.suptitle("Q6: Host Transcriptome Trajectory — Pathway Enrichment", fontsize=14, y=1.01)
plt.savefig('13_Q6_pathway_enrichment.png', bbox_inches='tight')
plt.close()
print("\nSaved: 13_Q6_pathway_enrichment.png")


# ============================================================
# COMPREHENSIVE SUMMARY FIGURE
# ============================================================
print("\n" + "=" * 70)
print("GENERATING COMPREHENSIVE SUMMARY")
print("=" * 70)

# Build a master pathway landscape: stage x infection
# Focus on the most informative comparison: RB+ vs RB- across stages
summary_terms = {}
for stage in stages:
    for direction in ['up', 'down']:
        enr = q2_enrichments.get((stage, direction), {})
        for cat in ['GO_BP', 'GO_MF', 'KEGG']:
            df = enr.get(cat, pd.DataFrame())
            if len(df) > 0:
                for _, row in df[df['padj'] < 0.05].head(6).iterrows():
                    key = f"[{row['category']}] {row['term_name']}"
                    if key not in summary_terms:
                        summary_terms[key] = {}
                    sign = 1 if direction == 'up' else -1
                    summary_terms[key][stage] = sign * (-np.log10(max(row['padj'], 1e-30)))

fig = plt.figure(figsize=(18, 16))

if summary_terms:
    hm_summary = pd.DataFrame(summary_terms).T.fillna(0)
    hm_summary = hm_summary.reindex(columns=[s for s in stages if s in hm_summary.columns])
    hm_var = hm_summary.var(axis=1).nlargest(35).index
    hm_summary = hm_summary.loc[hm_var]

    ax = plt.gca()
    sns.heatmap(hm_summary, cmap='RdBu_r', center=0,
                yticklabels=[textwrap.fill(str(t)[:60], 50) for t in hm_summary.index],
                xticklabels=True,
                cbar_kws={'label': 'Signed -log10(padj): +Red=UP in RB(+), -Blue=DOWN in RB(+)'},
                linewidths=0.5)
    ax.set_title("Master Pathway Landscape: GRBV Response Across Maturity Stages\n"
                 "(GO Biological Process, Molecular Function, KEGG)", fontsize=13)
    ax.tick_params(axis='y', labelsize=7)

plt.tight_layout()
plt.savefig('14_master_pathway_landscape.png', bbox_inches='tight')
plt.close()
print("Saved: 14_master_pathway_landscape.png")


# ============================================================
# FINAL SUMMARY PRINTOUT
# ============================================================
print("\n" + "=" * 70)
print("COMPLETE PATHWAY-LEVEL SUMMARY")
print("=" * 70)

print("""
=========================================================================
Q1: PATHWAY CHANGES WITH YEARS OF INFECTION
=========================================================================
""")
for ysi in ysi_values:
    print(f"\n  YSI = {ysi}:")
    for direction in ['up', 'down']:
        enr = q1_enrichments.get((ysi, direction), {})
        combined = get_combined_top(enr, n_per_cat=3)
        if len(combined) > 0:
            sig = combined[combined['padj'] < 0.05]
            if len(sig) > 0:
                terms = '; '.join(sig['term_name'].head(5).tolist())
                print(f"    {direction.upper()}: {terms}")
            else:
                print(f"    {direction.upper()}: No significant pathways")
        else:
            print(f"    {direction.upper()}: No enrichment")

print("""
=========================================================================
Q2: MATURITY STAGE PATHWAY ENRICHMENT
=========================================================================
""")
for stage in stages:
    print(f"\n  {stage}:")
    for direction in ['up', 'down']:
        enr = q2_enrichments.get((stage, direction), {})
        combined = get_combined_top(enr, n_per_cat=3)
        if len(combined) > 0:
            sig = combined[combined['padj'] < 0.05]
            if len(sig) > 0:
                terms = '; '.join(sig['term_name'].head(5).tolist())
                print(f"    {direction.upper()}: {terms}")
            else:
                print(f"    {direction.upper()}: No significant pathways")
        else:
            print(f"    {direction.upper()}: No enrichment")

print("""
=========================================================================
Q3: SEASONAL PATHWAY EFFECTS
=========================================================================
""")
for season in [2021, 2022]:
    print(f"\n  Season {season} infection response:")
    for direction in ['up', 'down']:
        enr = q3_enrichments.get((season, 'infection', direction), {})
        combined = get_combined_top(enr, n_per_cat=3)
        if len(combined) > 0:
            sig = combined[combined['padj'] < 0.05]
            if len(sig) > 0:
                terms = '; '.join(sig['term_name'].head(5).tolist())
                print(f"    {direction.upper()}: {terms}")
            else:
                print(f"    {direction.upper()}: No significant pathways")
        else:
            print(f"    {direction.upper()}: No enrichment")

print("""
=========================================================================
Q4: ENVIRONMENTAL PATHWAY EFFECTS
=========================================================================
""")
for stage in stages:
    print(f"\n  {stage}:")
    for grp_label, grp_name in [('pos', 'RB(+)'), ('neg', 'RB(-)')]:
        for direction in ['up', 'down']:
            enr = q4_env.get((grp_label, stage, direction), {})
            combined = get_combined_top(enr, n_per_cat=2)
            if len(combined) > 0:
                sig = combined[combined['padj'] < 0.05]
                if len(sig) > 0:
                    terms = '; '.join(sig['term_name'].head(3).tolist())
                    print(f"    {grp_name} {direction.upper()}: {terms}")

print("""
=========================================================================
Q5: UNINFECTED ENVIRONMENT PATHWAYS
=========================================================================
""")
for stage in stages:
    print(f"\n  {stage} (2022 vs 2021 in RB-):")
    for direction in ['up', 'down']:
        enr = q5_stage.get((stage, direction), {})
        combined = get_combined_top(enr, n_per_cat=3)
        if len(combined) > 0:
            sig = combined[combined['padj'] < 0.05]
            if len(sig) > 0:
                terms = '; '.join(sig['term_name'].head(5).tolist())
                print(f"    {direction.upper()}: {terms}")

print("""
=========================================================================
Q6: TRANSCRIPTOME TRAJECTORY PATHWAYS
=========================================================================
""")
for age in age_groups:
    print(f"\n  Age {age}yr vs controls:")
    for direction in ['up', 'down']:
        enr = q6_enrichments.get((age, direction), {})
        combined = get_combined_top(enr, n_per_cat=3)
        if len(combined) > 0:
            sig = combined[combined['padj'] < 0.05]
            if len(sig) > 0:
                terms = '; '.join(sig['term_name'].head(5).tolist())
                print(f"    {direction.upper()}: {terms}")
            else:
                print(f"    {direction.upper()}: No significant pathways")
        else:
            print(f"    {direction.upper()}: No enrichment")

print("\n" + "=" * 70)
print("ALL PATHWAY PLOTS GENERATED:")
print("  08_Q1_pathway_enrichment.png")
print("  09_Q2_pathway_enrichment.png")
print("  10_Q3_pathway_enrichment.png")
print("  11_Q4_pathway_enrichment.png")
print("  12_Q5_pathway_enrichment.png")
print("  13_Q6_pathway_enrichment.png")
print("  14_master_pathway_landscape.png")
print("=" * 70)
