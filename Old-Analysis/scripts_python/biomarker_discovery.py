#!/usr/bin/env python3
"""
GRBV Biomarker Discovery Analysis.
Identifies robust biomarker genes that are:
  1. Consistently DE across both seasons (2021, 2022)
  2. Consistently DE across all maturity stages
  3. Consistently DE across infection age cohorts
  4. High classification accuracy (AUC-ROC)
  5. Functionally annotated (GO/KEGG)

Criteria for a "robust biomarker":
  - Significant (padj < 0.05, |log2FC| >= 1) in MULTIPLE conditions
  - Consistent direction (always up or always down in infected)
  - High effect size (large fold change, low p-value)
  - High discriminatory power (AUC > 0.8)
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
from scipy.stats import fisher_exact, rankdata, mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
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
print("BIOMARKER DISCOVERY: LOADING DATA")
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

meta['Actual_Infection_Age'] = meta.apply(
    lambda row: int(row['Sampling_Season']) - int(row['Years_Since_Infection'])
    if row['Years_Since_Infection'] > 0 else 0, axis=1
)

print(f"Genes: {counts_filtered.shape[0]}, Samples: {counts_filtered.shape[1]}")

# Load annotations
db_path = 'org.Vvinifera.eg.db/inst/extdata/org.Vvinifera.eg.sqlite'
conn = sqlite3.connect(db_path)
gene_info = pd.read_sql("SELECT SYMBOL, GENENAME FROM gene_info", conn)
gene_name_map = dict(zip(gene_info['SYMBOL'], gene_info['GENENAME']))

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
except Exception:
    kegg_names = {}

from goatools.obo_parser import GODag
godag = GODag('go-basic.obo')

# Build gene -> GO/KEGG mappings
gene2go = defaultdict(set)
gene2kegg_set = defaultdict(set)
for _, row in go_df.iterrows():
    gene2go[row['SYMBOL']].add(f"{row['GO']} ({row['ONTOLOGY']})")
for _, row in kegg_df.iterrows():
    gene2kegg_set[row['SYMBOL']].add(row['KEGG'])

# Enrichment setup
go2genes_bp, go2genes_mf, go2genes_cc, kegg2genes = defaultdict(set), defaultdict(set), defaultdict(set), defaultdict(set)
for _, row in go_df.iterrows():
    gene, term, ont = row['SYMBOL'], row['GO'], row['ONTOLOGY']
    if ont == 'BP': go2genes_bp[term].add(gene)
    elif ont == 'MF': go2genes_mf[term].add(gene)
    elif ont == 'CC': go2genes_cc[term].add(gene)
for _, row in kegg_df.iterrows():
    kegg2genes[row['KEGG']].add(row['SYMBOL'])

all_genes = set(counts_filtered.index)
bg_go_bp = all_genes & set(go_df[go_df['ONTOLOGY'] == 'BP']['SYMBOL'])
bg_go_mf = all_genes & set(go_df[go_df['ONTOLOGY'] == 'MF']['SYMBOL'])
bg_kegg = all_genes & set(kegg_df['SYMBOL'])

conn.close()

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
seasons = [2021, 2022]


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_samples(meta_df, **kwargs):
    mask = pd.Series(True, index=meta_df.index)
    for col, val in kwargs.items():
        mask &= meta_df[col].isin(val) if isinstance(val, list) else meta_df[col] == val
    return list(meta_df[mask].index)


def run_de(group1_samples, group2_samples):
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


def fisher_enrichment(query_genes, term2genes, background_genes, min_genes=3):
    query = set(query_genes) & background_genes
    n_query, n_bg = len(query), len(background_genes)
    results = []
    for term, term_genes in term2genes.items():
        term_in_bg = term_genes & background_genes
        if len(term_in_bg) < min_genes: continue
        a = len(query & term_in_bg)
        if a < 2: continue
        b, c = len(query - term_in_bg), len(term_in_bg - query)
        d = n_bg - len(query | term_in_bg)
        _, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        results.append({
            'term': term, 'overlap': a, 'term_size': len(term_in_bg),
            'fold_enrichment': (a/n_query)/(len(term_in_bg)/n_bg) if n_query > 0 and len(term_in_bg) > 0 else 0,
            'pvalue': pval
        })
    if not results: return pd.DataFrame()
    df = pd.DataFrame(results)
    df['padj'] = np.minimum(df['pvalue'] * len(df) / rankdata(df['pvalue']), 1.0)
    return df.sort_values('pvalue')


# ============================================================
# 2. MULTI-CONDITION DE ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: RUNNING DE ACROSS ALL CONDITIONS")
print("=" * 70)

# Collect DE results across all conditions
condition_de = {}

# A. Overall: RB(+) vs RB(-)
all_pos = get_samples(meta, Infection_Group='RB(+)')
all_neg = get_samples(meta, Infection_Group='RB(-)')
de_overall = run_de(all_pos, all_neg)
condition_de['Overall'] = de_overall
print(f"Overall: {len(get_degs(de_overall))} DEGs")

# B. Per season
for season in seasons:
    pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season)
    neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season)
    de = run_de(pos, neg)
    condition_de[f'Season_{season}'] = de
    print(f"Season {season}: {len(get_degs(de))} DEGs")

# C. Per maturity stage
for stage in stages:
    pos = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=stage)
    neg = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=stage)
    de = run_de(pos, neg)
    condition_de[f'Stage_{stage}'] = de
    print(f"Stage {stage}: {len(get_degs(de))} DEGs")

# D. Per season × stage
for season in seasons:
    for stage in stages:
        pos = get_samples(meta, Infection_Group='RB(+)', Sampling_Season=season, Maturity_stage=stage)
        neg = get_samples(meta, Infection_Group='RB(-)', Sampling_Season=season, Maturity_stage=stage)
        if len(pos) < 3 or len(neg) < 3:
            continue
        de = run_de(pos, neg)
        condition_de[f'{season}_{stage}'] = de
        print(f"{season}-{stage}: {len(get_degs(de))} DEGs")

# E. Per infection age
infected_meta = meta[meta['Infection_Group'] == 'RB(+)']
age_groups = sorted(infected_meta['Actual_Infection_Age'].unique())
for age in age_groups:
    infected = list(infected_meta[infected_meta['Actual_Infection_Age'] == age].index)
    de = run_de(infected, all_neg)
    condition_de[f'Age_{age}yr'] = de
    print(f"Age {age}yr: {len(get_degs(de))} DEGs")


# ============================================================
# 3. BIOMARKER SCORING: CONSISTENCY + EFFECT SIZE + DIRECTION
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: SCORING GENES FOR BIOMARKER POTENTIAL")
print("=" * 70)

# For each gene, compute biomarker scores
all_gene_ids = log_cpm.index.tolist()
biomarker_scores = []

for gene in all_gene_ids:
    score = {
        'gene': gene,
        'gene_name': gene_name_map.get(gene, ''),
        'n_conditions_sig': 0,
        'n_conditions_tested': 0,
        'directions': [],
        'log2FCs': [],
        'padjs': [],
        'mean_abs_fc': 0,
        'min_padj': 1.0,
        'direction_consistent': False,

        # Season consistency
        'sig_2021': False, 'sig_2022': False,
        'fc_2021': 0, 'fc_2022': 0,

        # Stage consistency
        'n_stages_sig': 0,
        'stage_fcs': {},

        # Season×Stage consistency
        'n_season_stage_sig': 0,

        # Age consistency
        'n_ages_sig': 0,
    }

    for cond_name, de_df in condition_de.items():
        gene_row = de_df[de_df['gene'] == gene]
        if len(gene_row) == 0:
            continue

        score['n_conditions_tested'] += 1
        fc = gene_row.iloc[0]['log2FC']
        padj = gene_row.iloc[0]['padj']

        if padj < 0.05 and abs(fc) >= 1.0:
            score['n_conditions_sig'] += 1
            score['directions'].append(1 if fc > 0 else -1)
            score['log2FCs'].append(fc)
            score['padjs'].append(padj)

            if cond_name == 'Season_2021':
                score['sig_2021'] = True; score['fc_2021'] = fc
            elif cond_name == 'Season_2022':
                score['sig_2022'] = True; score['fc_2022'] = fc

            if cond_name.startswith('Stage_'):
                score['n_stages_sig'] += 1
                stage = cond_name.replace('Stage_', '')
                score['stage_fcs'][stage] = fc

            if cond_name.startswith('202') and '_' in cond_name:
                parts = cond_name.split('_', 1)
                if parts[0] in ['2021', '2022']:
                    score['n_season_stage_sig'] += 1

            if cond_name.startswith('Age_'):
                score['n_ages_sig'] += 1
        else:
            # Track FC even if not significant (for overall direction)
            if cond_name.startswith('Season_'):
                if '2021' in cond_name: score['fc_2021'] = fc
                elif '2022' in cond_name: score['fc_2022'] = fc

    if score['directions']:
        score['direction_consistent'] = len(set(score['directions'])) == 1
        score['mean_abs_fc'] = np.mean(np.abs(score['log2FCs']))
        score['min_padj'] = min(score['padjs'])
        score['primary_direction'] = 'UP' if np.mean(score['directions']) > 0 else 'DOWN'
    else:
        score['primary_direction'] = 'NS'

    biomarker_scores.append(score)

bm_df = pd.DataFrame(biomarker_scores)

# Composite biomarker score:
# Weighted combination of consistency, effect size, and significance
bm_df['consistency_score'] = bm_df['n_conditions_sig'] / bm_df['n_conditions_tested'].clip(lower=1)
bm_df['effect_score'] = bm_df['mean_abs_fc'].clip(upper=5) / 5  # Normalize to 0-1
bm_df['significance_score'] = (-np.log10(bm_df['min_padj'].clip(lower=1e-50))) / 50  # Normalize

bm_df['season_robust'] = bm_df['sig_2021'] & bm_df['sig_2022'] & bm_df['direction_consistent']

# Overall biomarker score
bm_df['biomarker_score'] = (
    0.30 * bm_df['consistency_score'] +
    0.25 * bm_df['effect_score'] +
    0.20 * bm_df['significance_score'] +
    0.15 * (bm_df['direction_consistent'].astype(float)) +
    0.10 * (bm_df['season_robust'].astype(float))
)

bm_df = bm_df.sort_values('biomarker_score', ascending=False)

print(f"\nTotal genes scored: {len(bm_df)}")
print(f"Genes significant in >=1 condition: {(bm_df['n_conditions_sig'] >= 1).sum()}")
print(f"Genes significant in >=5 conditions: {(bm_df['n_conditions_sig'] >= 5).sum()}")
print(f"Genes significant in >=10 conditions: {(bm_df['n_conditions_sig'] >= 10).sum()}")
print(f"Season-robust (sig in both 2021 & 2022, same direction): {bm_df['season_robust'].sum()}")
print(f"Stage-consistent (sig in >=3 stages): {(bm_df['n_stages_sig'] >= 3).sum()}")


# ============================================================
# 4. AUC-ROC ANALYSIS FOR TOP CANDIDATES
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: AUC-ROC CLASSIFICATION ANALYSIS")
print("=" * 70)

# Calculate individual gene AUC for RB(+) vs RB(-)
labels = (meta.loc[log_cpm.columns, 'Infection_Group'] == 'RB(+)').astype(int).values

# Top 500 biomarker candidates
top_candidates = bm_df[bm_df['n_conditions_sig'] >= 1].head(500)['gene'].tolist()

auc_results = []
for gene in top_candidates:
    expr = log_cpm.loc[gene].values
    try:
        auc = roc_auc_score(labels, expr)
        # Flip AUC if < 0.5 (gene is downregulated)
        auc_directed = max(auc, 1 - auc)
        auc_results.append({'gene': gene, 'AUC': auc_directed, 'raw_AUC': auc})
    except Exception:
        pass

auc_df = pd.DataFrame(auc_results).sort_values('AUC', ascending=False)
print(f"\nAUC computed for {len(auc_df)} genes")
print(f"Genes with AUC > 0.90: {(auc_df['AUC'] > 0.90).sum()}")
print(f"Genes with AUC > 0.85: {(auc_df['AUC'] > 0.85).sum()}")
print(f"Genes with AUC > 0.80: {(auc_df['AUC'] > 0.80).sum()}")
print(f"Genes with AUC > 0.75: {(auc_df['AUC'] > 0.75).sum()}")

# Merge AUC into biomarker df
bm_df = bm_df.merge(auc_df[['gene', 'AUC']], on='gene', how='left')
bm_df['AUC'] = bm_df['AUC'].fillna(0.5)

# Update biomarker score with AUC
bm_df['biomarker_score_final'] = bm_df['biomarker_score'] * 0.7 + (bm_df['AUC'] - 0.5) * 2 * 0.3
bm_df = bm_df.sort_values('biomarker_score_final', ascending=False)


# ============================================================
# 5. SEASON-SPECIFIC AUC (robustness check)
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: SEASON-SPECIFIC AUC (Robustness check)")
print("=" * 70)

# Test if biomarkers work in each season independently
top50 = bm_df.head(50)['gene'].tolist()

season_auc = {}
for season in seasons:
    s_mask = meta['Sampling_Season'] == season
    s_samples = meta[s_mask].index.tolist()
    s_labels = (meta.loc[s_samples, 'Infection_Group'] == 'RB(+)').astype(int).values

    for gene in top50:
        expr = log_cpm.loc[gene, s_samples].values
        try:
            auc = roc_auc_score(s_labels, expr)
            auc = max(auc, 1 - auc)
            season_auc[(gene, season)] = auc
        except Exception:
            season_auc[(gene, season)] = 0.5

# Stage-specific AUC
stage_auc = {}
for stage in stages:
    st_mask = meta['Maturity_stage'] == stage
    st_samples = meta[st_mask].index.tolist()
    st_labels = (meta.loc[st_samples, 'Infection_Group'] == 'RB(+)').astype(int).values

    for gene in top50:
        expr = log_cpm.loc[gene, st_samples].values
        try:
            auc = roc_auc_score(st_labels, expr)
            auc = max(auc, 1 - auc)
            stage_auc[(gene, stage)] = auc
        except Exception:
            stage_auc[(gene, stage)] = 0.5


# ============================================================
# 6. MULTI-GENE PANEL EVALUATION
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: MULTI-GENE PANEL CLASSIFIER")
print("=" * 70)

# Build panels of 5, 10, 20, 50 genes and test with cross-validation
panel_sizes = [3, 5, 10, 20, 50]
panel_results = {}

for panel_size in panel_sizes:
    panel_genes = bm_df.head(panel_size)['gene'].tolist()
    X = log_cpm.loc[panel_genes].T.values
    y = labels

    # 5-fold stratified CV with logistic regression
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs_lr = []
    aucs_rf = []

    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Logistic Regression
        lr = LogisticRegression(max_iter=1000, random_state=42)
        lr.fit(X_train, y_train)
        prob_lr = lr.predict_proba(X_test)[:, 1]
        aucs_lr.append(roc_auc_score(y_test, prob_lr))

        # Random Forest
        rf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
        rf.fit(X_train, y_train)
        prob_rf = rf.predict_proba(X_test)[:, 1]
        aucs_rf.append(roc_auc_score(y_test, prob_rf))

    panel_results[panel_size] = {
        'lr_mean': np.mean(aucs_lr), 'lr_std': np.std(aucs_lr),
        'rf_mean': np.mean(aucs_rf), 'rf_std': np.std(aucs_rf),
        'genes': panel_genes
    }
    print(f"  Panel {panel_size} genes: LR AUC={np.mean(aucs_lr):.4f}±{np.std(aucs_lr):.4f}, "
          f"RF AUC={np.mean(aucs_rf):.4f}±{np.std(aucs_rf):.4f}")

# Season-specific panel validation
print("\nSeason-specific panel validation (top 10 genes):")
panel_10 = bm_df.head(10)['gene'].tolist()
for season in seasons:
    s_mask = meta['Sampling_Season'] == season
    s_samples = meta[s_mask].index.tolist()
    X_s = log_cpm.loc[panel_10, s_samples].T.values
    y_s = (meta.loc[s_samples, 'Infection_Group'] == 'RB(+)').astype(int).values

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []
    for train_idx, test_idx in cv.split(X_s, y_s):
        lr = LogisticRegression(max_iter=1000, random_state=42)
        lr.fit(X_s[train_idx], y_s[train_idx])
        prob = lr.predict_proba(X_s[test_idx])[:, 1]
        aucs.append(roc_auc_score(y_s[test_idx], prob))
    print(f"  Season {season}: LR AUC = {np.mean(aucs):.4f} ± {np.std(aucs):.4f}")

# Stage-specific panel validation
print("\nStage-specific panel validation (top 10 genes):")
for stage in stages:
    st_mask = meta['Maturity_stage'] == stage
    st_samples = meta[st_mask].index.tolist()
    X_st = log_cpm.loc[panel_10, st_samples].T.values
    y_st = (meta.loc[st_samples, 'Infection_Group'] == 'RB(+)').astype(int).values

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []
    for train_idx, test_idx in cv.split(X_st, y_st):
        lr = LogisticRegression(max_iter=1000, random_state=42)
        lr.fit(X_st[train_idx], y_st[train_idx])
        prob = lr.predict_proba(X_st[test_idx])[:, 1]
        aucs.append(roc_auc_score(y_st[test_idx], prob))
    print(f"  {stage}: LR AUC = {np.mean(aucs):.4f} ± {np.std(aucs):.4f}")


# ============================================================
# 7. CATEGORIZE BIOMARKERS
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: BIOMARKER CATEGORIES")
print("=" * 70)

# Category 1: Universal biomarkers (sig in both seasons AND >=3 stages)
universal = bm_df[(bm_df['season_robust']) & (bm_df['n_stages_sig'] >= 3) & (bm_df['AUC'] >= 0.7)]
print(f"\n1. UNIVERSAL BIOMARKERS (season-robust + ≥3 stages + AUC≥0.7): {len(universal)}")

# Category 2: Harvest-specific biomarkers (strongest at Harvest)
harvest_sig = set(get_degs(condition_de.get('Stage_Harvest', pd.DataFrame()))['gene'])
non_early = set(bm_df[(bm_df['n_stages_sig'] <= 1)]['gene'])
harvest_specific = bm_df[(bm_df['gene'].isin(harvest_sig)) & (bm_df['gene'].isin(non_early)) & (bm_df['AUC'] >= 0.7)]
print(f"2. HARVEST-SPECIFIC BIOMARKERS: {len(harvest_specific)}")

# Category 3: Early-detection biomarkers (sig at Pre_veraison or Veraison)
early_stages_genes = set()
for s in ['Pre_veraison', 'Veraison']:
    de = condition_de.get(f'Stage_{s}', pd.DataFrame())
    if len(de) > 0:
        early_stages_genes.update(get_degs(de)['gene'].tolist())
early_biomarkers = bm_df[(bm_df['gene'].isin(early_stages_genes)) & (bm_df['AUC'] >= 0.65)]
print(f"3. EARLY-DETECTION BIOMARKERS (Pre_veraison/Veraison, AUC≥0.65): {len(early_biomarkers)}")

# Category 4: Infection-age biomarkers (sig in >=3 age groups, correlated with age)
age_biomarkers = bm_df[(bm_df['n_ages_sig'] >= 3) & (bm_df['direction_consistent']) & (bm_df['AUC'] >= 0.7)]
print(f"4. INFECTION-AGE BIOMARKERS (≥3 age groups + consistent + AUC≥0.7): {len(age_biomarkers)}")


# ============================================================
# 8. FUNCTIONAL ANNOTATION OF TOP BIOMARKERS
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: FUNCTIONAL ANNOTATION OF TOP BIOMARKERS")
print("=" * 70)

# Top 30 overall biomarkers
top30 = bm_df.head(30)

print(f"\n{'='*120}")
print(f"{'Rank':<5} {'Gene':<18} {'Gene Name':<45} {'Score':<7} {'AUC':<7} {'Dir':<5} "
      f"{'#Sig':<5} {'Stages':<7} {'Season':<8} {'meanFC':<8}")
print(f"{'='*120}")

for rank, (_, row) in enumerate(top30.iterrows(), 1):
    name = str(row['gene_name'])[:42] if row['gene_name'] else ''
    season_str = '2021+2022' if row['season_robust'] else ('2021' if row['sig_2021'] else ('2022' if row['sig_2022'] else '-'))
    print(f"{rank:<5} {row['gene']:<18} {name:<45} {row['biomarker_score_final']:.4f} "
          f"{row['AUC']:.3f}  {row['primary_direction']:<4} "
          f"{row['n_conditions_sig']:<5} {row['n_stages_sig']:<7} {season_str:<8} "
          f"{row['mean_abs_fc']:.2f}")

# Enrichment of top biomarkers
print("\n\nPathway enrichment of top 100 biomarkers:")
top100_genes = bm_df.head(100)['gene'].tolist()
top100_up = bm_df[(bm_df.index.isin(bm_df.head(100).index)) & (bm_df['primary_direction'] == 'UP')]['gene'].tolist()
top100_down = bm_df[(bm_df.index.isin(bm_df.head(100).index)) & (bm_df['primary_direction'] == 'DOWN')]['gene'].tolist()

for direction, genes, label in [('UP', top100_up, 'Up-biomarkers'), ('DOWN', top100_down, 'Down-biomarkers')]:
    print(f"\n  {label} ({len(genes)} genes):")
    for name, t2g, bg in [('GO-BP', go2genes_bp, bg_go_bp), ('GO-MF', go2genes_mf, bg_go_mf), ('KEGG', kegg2genes, bg_kegg)]:
        df = fisher_enrichment(genes, t2g, bg)
        if len(df) > 0:
            if name.startswith('GO'):
                df['term_name'] = df['term'].apply(lambda x: godag[x].name if x in godag else x)
            else:
                df['term_name'] = df['term'].apply(lambda x: kegg_names.get(x, x))
            sig = df[df['padj'] < 0.05].head(5)
            if len(sig) > 0:
                for _, r in sig.iterrows():
                    print(f"    {name}: {r['term_name'][:55]:55s} FE={r['fold_enrichment']:.1f} padj={r['padj']:.2e} [{r['overlap']}/{r['term_size']}]")


# ============================================================
# PLOT 1: BIOMARKER OVERVIEW
# ============================================================
print("\n" + "=" * 70)
print("GENERATING PLOTS")
print("=" * 70)

fig = plt.figure(figsize=(24, 22))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.5, wspace=0.45)

# 1. Biomarker score distribution
ax = fig.add_subplot(gs[0, 0])
scores = bm_df[bm_df['n_conditions_sig'] >= 1]['biomarker_score_final']
ax.hist(scores, bins=50, color='steelblue', alpha=0.8, edgecolor='black', linewidth=0.3)
ax.axvline(scores.quantile(0.95), color='red', linestyle='--', linewidth=1, label=f'Top 5% ({scores.quantile(0.95):.3f})')
ax.axvline(scores.quantile(0.99), color='darkred', linestyle='--', linewidth=1, label=f'Top 1% ({scores.quantile(0.99):.3f})')
ax.set_xlabel("Biomarker Score")
ax.set_ylabel("Number of Genes")
ax.set_title("Biomarker Score Distribution\n(genes with ≥1 DE condition)")
ax.legend(fontsize=12)

# 2. Consistency vs AUC
ax = fig.add_subplot(gs[0, 1])
plot_df = bm_df[bm_df['n_conditions_sig'] >= 1].copy()
sc = ax.scatter(plot_df['n_conditions_sig'], plot_df['AUC'],
                c=plot_df['mean_abs_fc'], cmap='YlOrRd', s=15, alpha=0.6, edgecolors='none')
plt.colorbar(sc, ax=ax, label='Mean |log2FC|', shrink=0.8)
ax.set_xlabel("Number of Conditions Significant")
ax.set_ylabel("AUC-ROC")
ax.axhline(0.8, color='grey', linestyle='--', linewidth=0.5)
ax.set_title("Consistency vs Discriminatory Power")

# Annotate top 5
for _, row in bm_df.head(5).iterrows():
    if row['n_conditions_sig'] >= 1:
        ax.annotate(row['gene'], (row['n_conditions_sig'], row['AUC']),
                    fontsize=9, ha='left', va='bottom')

# 3. Multi-gene panel AUC
ax = fig.add_subplot(gs[0, 2])
sizes = list(panel_results.keys())
lr_means = [panel_results[s]['lr_mean'] for s in sizes]
lr_stds = [panel_results[s]['lr_std'] for s in sizes]
rf_means = [panel_results[s]['rf_mean'] for s in sizes]
rf_stds = [panel_results[s]['rf_std'] for s in sizes]

ax.errorbar(sizes, lr_means, yerr=lr_stds, marker='o', linewidth=2, capsize=5, label='Logistic Regression', color='#3498db')
ax.errorbar(sizes, rf_means, yerr=rf_stds, marker='s', linewidth=2, capsize=5, label='Random Forest', color='#e74c3c')
ax.set_xlabel("Panel Size (number of genes)")
ax.set_ylabel("5-fold CV AUC-ROC")
ax.set_title("Multi-Gene Panel Classifier\n(5-fold CV)")
ax.legend(fontsize=12)
ax.set_ylim(0.5, 1.05)
ax.axhline(0.9, color='grey', linestyle='--', linewidth=0.5)
ax.set_xticks(sizes)

# 4. ROC curves for top individual genes
ax = fig.add_subplot(gs[1, 0])
top5_genes = bm_df.head(5)['gene'].tolist()
for gene in top5_genes:
    expr = log_cpm.loc[gene].values
    auc_val = roc_auc_score(labels, expr)
    if auc_val < 0.5:
        expr = -expr
        auc_val = 1 - auc_val
    fpr, tpr, _ = roc_curve(labels, expr)
    gene_nm = gene_name_map.get(gene, '')[:25]
    ax.plot(fpr, tpr, linewidth=1.5, label=f'{gene} ({gene_nm})\nAUC={auc_val:.3f}')

ax.plot([0, 1], [0, 1], 'k--', linewidth=0.5)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("ROC Curves: Top 5 Biomarker Genes")
ax.legend(fontsize=10, loc='lower right')

# 5. Season-specific AUC heatmap for top 30
ax = fig.add_subplot(gs[1, 1])
top30_genes = bm_df.head(30)['gene'].tolist()
auc_matrix = pd.DataFrame(index=top30_genes, columns=['Overall', '2021', '2022'] + stages)

for gene in top30_genes:
    auc_matrix.loc[gene, 'Overall'] = bm_df[bm_df['gene'] == gene]['AUC'].values[0]
    for season in seasons:
        auc_matrix.loc[gene, str(season)] = season_auc.get((gene, season), 0.5)
    for stage in stages:
        auc_matrix.loc[gene, stage] = stage_auc.get((gene, stage), 0.5)

auc_matrix = auc_matrix.astype(float)
gene_labels = [f"{g} ({gene_name_map.get(g, '')[:20]})" for g in top30_genes]

sns.heatmap(auc_matrix, ax=ax, cmap='RdYlGn', vmin=0.5, vmax=1.0,
            yticklabels=gene_labels, xticklabels=True,
            cbar_kws={'label': 'AUC-ROC', 'shrink': 0.8},
            linewidths=0.5, annot=True, fmt='.2f', annot_kws={'fontsize': 5})
ax.tick_params(axis='y', labelsize=5.5)
ax.tick_params(axis='x', labelsize=7, rotation=45)
ax.set_title("AUC Across Conditions\n(Top 30 Biomarkers)")

# 6. Expression heatmap for top 20 biomarkers
ax = fig.add_subplot(gs[1, 2])
top20_genes = bm_df.head(20)['gene'].tolist()

# Order samples by infection status, then stage, then season
sample_order = meta.sort_values(['Infection_Group', 'Maturity_stage', 'Sampling_Season']).index
expr_matrix = log_cpm.loc[top20_genes, sample_order]

# Z-score normalize per gene
expr_z = (expr_matrix.T - expr_matrix.mean(axis=1)).T
expr_z = (expr_z.T / expr_matrix.std(axis=1).clip(lower=0.01)).T

# Truncate for visualization
expr_z = expr_z.clip(-3, 3)

sns.heatmap(expr_z, ax=ax, cmap='RdBu_r', center=0,
            yticklabels=[f"{g}" for g in top20_genes],
            xticklabels=False,
            cbar_kws={'label': 'Z-score', 'shrink': 0.6})
ax.set_title("Expression Heatmap: Top 20 Biomarkers\n(Z-scored, samples ordered by group)")
ax.tick_params(axis='y', labelsize=6)

# Add color bar for sample grouping
infection_colors = ['red' if meta.loc[s, 'Infection_Group'] == 'RB(+)' else 'blue' for s in sample_order]
for i, c in enumerate(infection_colors):
    ax.axvline(i, color=c, linewidth=0.1, alpha=0.3)

# 7. Volcano plot highlighting biomarkers
ax = fig.add_subplot(gs[2, 0])
de_overall_plot = de_overall.copy()
de_overall_plot['neg_log10p'] = -np.log10(de_overall_plot['padj'].clip(lower=1e-300))

# Color: grey=NS, light red/blue=DE, bright=biomarker
top100_set = set(bm_df.head(100)['gene'].tolist())
top20_set = set(top20_genes)
colors_vol = []
sizes_vol = []
for _, row in de_overall_plot.iterrows():
    if row['gene'] in top20_set:
        colors_vol.append('gold')
        sizes_vol.append(30)
    elif row['gene'] in top100_set:
        colors_vol.append('#e74c3c' if row['log2FC'] > 0 else '#3498db')
        sizes_vol.append(12)
    elif row['padj'] < 0.05 and abs(row['log2FC']) >= 1:
        colors_vol.append('#f5b7b1' if row['log2FC'] > 0 else '#aed6f1')
        sizes_vol.append(4)
    else:
        colors_vol.append('#d5d8dc')
        sizes_vol.append(2)

ax.scatter(de_overall_plot['log2FC'], de_overall_plot['neg_log10p'],
           c=colors_vol, s=sizes_vol, alpha=0.6, edgecolors='none')
ax.axvline(1, color='red', linestyle='--', linewidth=0.5)
ax.axvline(-1, color='blue', linestyle='--', linewidth=0.5)
ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)

# Label top 8
for _, row in bm_df.head(8).iterrows():
    gene_de = de_overall_plot[de_overall_plot['gene'] == row['gene']]
    if len(gene_de) > 0:
        gr = gene_de.iloc[0]
        ax.annotate(row['gene'], (gr['log2FC'], gr['neg_log10p']),
                    fontsize=8, fontweight='bold', ha='center', va='bottom',
                    xytext=(0, 5), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.7))

ax.set_xlabel("log2 Fold Change (RB+ vs RB-)")
ax.set_ylabel("-log10(padj)")
ax.set_title("Volcano: Top Biomarkers Highlighted\n(Gold=Top 20, Red/Blue=Top 100)")

# 8. Biomarker category summary
ax = fig.add_subplot(gs[2, 1])
categories = ['Universal\n(season+stage)', 'Early\ndetection', 'Harvest\nspecific', 'Infection\nage']
cat_counts = [len(universal), len(early_biomarkers), len(harvest_specific), len(age_biomarkers)]
cat_colors = ['#2ecc71', '#3498db', '#9b59b6', '#e67e22']
bars = ax.bar(categories, cat_counts, color=cat_colors, alpha=0.85, edgecolor='black', linewidth=0.5)
ax.set_ylabel("Number of Biomarker Genes")
ax.set_title("Biomarker Categories")
for b, v in zip(bars, cat_counts):
    ax.text(b.get_x() + b.get_width()/2, v + max(cat_counts)*0.02, str(v), ha='center', fontsize=15)

# 9. Expression boxplots for top 6 biomarkers
ax = fig.add_subplot(gs[2, 2])
top6 = bm_df.head(6)['gene'].tolist()

box_data = []
for gene in top6:
    for _, row_meta in meta.iterrows():
        sample = row_meta.name
        box_data.append({
            'Gene': gene,
            'Expression': log_cpm.loc[gene, sample],
            'Group': row_meta['Infection_Group']
        })

box_df = pd.DataFrame(box_data)
sns.boxplot(data=box_df, x='Gene', y='Expression', hue='Group',
            palette={'RB(+)': '#e74c3c', 'RB(-)': '#3498db'},
            ax=ax, fliersize=1, linewidth=0.8)
ax.set_title("Expression of Top 6 Biomarkers\n(RB+ vs RB-)")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
ax.legend(fontsize=12)

fig.suptitle("GRBV Biomarker Discovery: Overview", fontsize=22, y=1.02)
plt.savefig('20_biomarker_overview.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 20_biomarker_overview.png")


# ============================================================
# PLOT 2: DETAILED BIOMARKER PROFILES
# ============================================================

fig = plt.figure(figsize=(24, 24))
gs = gridspec.GridSpec(4, 4, figure=fig, hspace=0.6, wspace=0.5)

# Individual gene profiles for top 16 biomarkers
top16 = bm_df.head(16)

for idx, (_, bm_row) in enumerate(top16.iterrows()):
    gene = bm_row['gene']
    row_idx, col_idx = idx // 4, idx % 4
    ax = fig.add_subplot(gs[row_idx, col_idx])

    # Box plot: season × stage
    plot_data = []
    for season in seasons:
        for stage in stages:
            for grp in ['RB(+)', 'RB(-)']:
                samples = get_samples(meta, Infection_Group=grp, Sampling_Season=season, Maturity_stage=stage)
                for s in samples:
                    plot_data.append({
                        'Season': str(season),
                        'Stage': stage,
                        'Group': grp,
                        'Expression': log_cpm.loc[gene, s],
                        'x_label': f"{season}\n{stage[:4]}"
                    })

    pdf = pd.DataFrame(plot_data)

    # Plot means with error bars for each season×stage×group
    x_pos = 0
    xticks = []
    xticklabels = []
    for season in seasons:
        for s_idx, stage in enumerate(stages):
            for grp, color, offset in [('RB(-)', '#3498db', -0.15), ('RB(+)', '#e74c3c', 0.15)]:
                sub = pdf[(pdf['Season'] == str(season)) & (pdf['Stage'] == stage) & (pdf['Group'] == grp)]
                if len(sub) > 0:
                    mean_v = sub['Expression'].mean()
                    sem_v = sub['Expression'].sem()
                    ax.bar(x_pos + offset, mean_v, width=0.25, color=color, alpha=0.7)
                    ax.errorbar(x_pos + offset, mean_v, yerr=sem_v, color='black', capsize=2, linewidth=0.5)

            xticks.append(x_pos)
            xticklabels.append(f"{stage[:4]}")
            x_pos += 1
        x_pos += 0.5  # gap between seasons

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=8, rotation=45)
    ax.axvline(3.75, color='grey', linestyle=':', linewidth=0.5)  # season separator

    gene_nm = gene_name_map.get(gene, '')
    gene_nm_short = gene_nm[:30] if gene_nm else ''
    ax.set_title(f"{gene}\n{gene_nm_short}\nAUC={bm_row['AUC']:.3f} | {bm_row['primary_direction']} | "
                 f"Sig in {bm_row['n_conditions_sig']}/{bm_row['n_conditions_tested']} cond.",
                 fontsize=10)
    ax.set_ylabel("log2(CPM+1)", fontsize=10)

    # Add season labels
    ax.text(1.5, ax.get_ylim()[1]*0.95, "2021", ha='center', fontsize=12, fontweight='bold', color='#2c3e50')
    ax.text(6, ax.get_ylim()[1]*0.95, "2022", ha='center', fontsize=12, fontweight='bold', color='#2c3e50')

fig.suptitle("Top 16 GRBV Biomarkers: Expression Across Season × Stage × Infection\n"
             "(Red = RB+, Blue = RB-; bars = mean ± SEM)",
             fontsize=21, y=1.02)
plt.savefig('21_biomarker_profiles.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 21_biomarker_profiles.png")


# ============================================================
# PLOT 3: BIOMARKER ROBUSTNESS & VALIDATION
# ============================================================

fig = plt.figure(figsize=(22, 18))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.5, wspace=0.45)

# 1. FC consistency: 2021 vs 2022 for top biomarkers
ax = fig.add_subplot(gs[0, 0])
top50_genes = bm_df.head(50)['gene'].tolist()

de_21 = condition_de['Season_2021'].set_index('gene')
de_22 = condition_de['Season_2022'].set_index('gene')

fc21 = [de_21.loc[g, 'log2FC'] if g in de_21.index else 0 for g in top50_genes]
fc22 = [de_22.loc[g, 'log2FC'] if g in de_22.index else 0 for g in top50_genes]

ax.scatter(fc21, fc22, c='#e74c3c', s=40, alpha=0.7, edgecolors='black', linewidth=0.3)
r, p = stats.pearsonr(fc21, fc22)
ax.plot([-6, 6], [-6, 6], 'k--', linewidth=0.5)
ax.set_xlabel("log2FC in 2021")
ax.set_ylabel("log2FC in 2022")
ax.set_title(f"FC Consistency of Top 50 Biomarkers\n(r={r:.3f}, p={p:.2e})")

# Label extreme points
for i, gene in enumerate(top50_genes[:8]):
    ax.annotate(gene, (fc21[i], fc22[i]), fontsize=8, ha='left')

# 2. AUC stability: box plot of AUC across conditions
ax = fig.add_subplot(gs[0, 1])
auc_stability_data = []
for gene in top30_genes:
    aucs_gene = []
    aucs_gene.append(bm_df[bm_df['gene'] == gene]['AUC'].values[0])
    for season in seasons:
        aucs_gene.append(season_auc.get((gene, season), 0.5))
    for stage in stages:
        aucs_gene.append(stage_auc.get((gene, stage), 0.5))
    auc_stability_data.append({
        'gene': gene,
        'mean_auc': np.mean(aucs_gene),
        'min_auc': np.min(aucs_gene),
        'std_auc': np.std(aucs_gene)
    })

stab_df = pd.DataFrame(auc_stability_data).sort_values('min_auc', ascending=False)
x = np.arange(len(stab_df))
ax.barh(x, stab_df['mean_auc'].values, xerr=stab_df['std_auc'].values,
        color='#2ecc71', alpha=0.8, capsize=2)
ax.scatter(stab_df['min_auc'].values, x, color='red', s=20, zorder=5, label='Min AUC')
ax.set_yticks(x)
gene_labels_stab = [f"{g} ({gene_name_map.get(g, '')[:18]})" for g in stab_df['gene']]
ax.set_yticklabels(gene_labels_stab, fontsize=8)
ax.set_xlabel("AUC-ROC")
ax.set_title("AUC Stability: Top 30 Biomarkers\n(bar=mean±std, red=minimum)")
ax.axvline(0.7, color='grey', linestyle='--', linewidth=0.5)
ax.legend(fontsize=10)
ax.invert_yaxis()

# 3. Condition-specific signature: where each biomarker is significant
ax = fig.add_subplot(gs[0, 2])
sig_matrix_data = {}
conditions_ordered = ['Overall', 'Season_2021', 'Season_2022'] + \
                      [f'Stage_{s}' for s in stages] + \
                      [f'{s}_{st}' for s in seasons for st in stages]

for gene in top30_genes:
    sig_matrix_data[gene] = {}
    for cond in conditions_ordered:
        de = condition_de.get(cond, pd.DataFrame())
        if len(de) > 0:
            gene_row = de[de['gene'] == gene]
            if len(gene_row) > 0:
                padj = gene_row.iloc[0]['padj']
                fc = gene_row.iloc[0]['log2FC']
                if padj < 0.05 and abs(fc) >= 1.0:
                    sig_matrix_data[gene][cond] = np.sign(fc)
                else:
                    sig_matrix_data[gene][cond] = 0
            else:
                sig_matrix_data[gene][cond] = 0
        else:
            sig_matrix_data[gene][cond] = np.nan

sig_matrix = pd.DataFrame(sig_matrix_data).T
sig_matrix = sig_matrix.reindex(columns=[c for c in conditions_ordered if c in sig_matrix.columns])

# Shorten column names for display
col_map = {}
for c in sig_matrix.columns:
    if c == 'Overall': col_map[c] = 'All'
    elif c.startswith('Season_'): col_map[c] = c.replace('Season_', 'S')
    elif c.startswith('Stage_'): col_map[c] = c.replace('Stage_', '').replace('_veraison', 'V').replace('Pre_', 'PreV').replace('Post_', 'PostV')
    elif c.startswith('202'):
        parts = c.split('_', 1)
        col_map[c] = f"{parts[0][2:]}-{parts[1][:4]}"
    else:
        col_map[c] = c[:8]

sig_matrix_display = sig_matrix.rename(columns=col_map)

sns.heatmap(sig_matrix_display, ax=ax, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
            yticklabels=[f"{g}" for g in sig_matrix.index],
            xticklabels=True,
            cbar_kws={'label': 'Direction (+1=Up, -1=Down)', 'shrink': 0.6},
            linewidths=0.5)
ax.tick_params(axis='y', labelsize=5.5)
ax.tick_params(axis='x', labelsize=5, rotation=90)
ax.set_title("Significance Matrix\n(Top 30 Biomarkers × Conditions)")

# 4. Panel performance across seasons and stages
ax = fig.add_subplot(gs[1, 0])
panel_10_genes = bm_df.head(10)['gene'].tolist()

validation_results = {}
for label, mask_fn in [('All', lambda m: pd.Series(True, index=m.index)),
                        ('2021', lambda m: m['Sampling_Season'] == 2021),
                        ('2022', lambda m: m['Sampling_Season'] == 2022),
                        ('Pre_ver', lambda m: m['Maturity_stage'] == 'Pre_veraison'),
                        ('Veraison', lambda m: m['Maturity_stage'] == 'Veraison'),
                        ('Post_ver', lambda m: m['Maturity_stage'] == 'Post_veraison'),
                        ('Harvest', lambda m: m['Maturity_stage'] == 'Harvest')]:
    mask = mask_fn(meta)
    s_samples = meta[mask].index.tolist()
    X_v = log_cpm.loc[panel_10_genes, s_samples].T.values
    y_v = (meta.loc[s_samples, 'Infection_Group'] == 'RB(+)').astype(int).values

    if len(np.unique(y_v)) < 2:
        continue

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []
    for train_idx, test_idx in cv.split(X_v, y_v):
        lr = LogisticRegression(max_iter=1000, random_state=42)
        lr.fit(X_v[train_idx], y_v[train_idx])
        prob = lr.predict_proba(X_v[test_idx])[:, 1]
        aucs.append(roc_auc_score(y_v[test_idx], prob))
    validation_results[label] = {'mean': np.mean(aucs), 'std': np.std(aucs)}

v_labels = list(validation_results.keys())
v_means = [validation_results[l]['mean'] for l in v_labels]
v_stds = [validation_results[l]['std'] for l in v_labels]

bars = ax.bar(range(len(v_labels)), v_means, yerr=v_stds, capsize=3,
              color=['#2ecc71', '#3498db', '#e67e22', '#2ca02c', '#ff7f0e', '#d62728', '#9467bd'],
              alpha=0.85, edgecolor='black', linewidth=0.5)
ax.set_xticks(range(len(v_labels)))
ax.set_xticklabels(v_labels, rotation=45, ha='right', fontsize=12)
ax.set_ylabel("5-fold CV AUC")
ax.set_title("10-Gene Panel Performance\nAcross Conditions")
ax.set_ylim(0.5, 1.05)
ax.axhline(0.9, color='grey', linestyle='--', linewidth=0.5)

for b, m in zip(bars, v_means):
    ax.text(b.get_x() + b.get_width()/2, m + 0.02, f'{m:.3f}', ha='center', fontsize=10)

# 5. Feature importance from Random Forest
ax = fig.add_subplot(gs[1, 1])
panel_20_genes = bm_df.head(20)['gene'].tolist()
X_rf = log_cpm.loc[panel_20_genes].T.values
y_rf = labels

rf = RandomForestClassifier(n_estimators=500, random_state=42, max_depth=8)
rf.fit(X_rf, y_rf)
importances = rf.feature_importances_

imp_df = pd.DataFrame({'gene': panel_20_genes, 'importance': importances})
imp_df = imp_df.sort_values('importance', ascending=True)

ax.barh(range(len(imp_df)), imp_df['importance'].values, color='#9b59b6', alpha=0.85)
ax.set_yticks(range(len(imp_df)))
gene_labels_imp = [f"{g} ({gene_name_map.get(g, '')[:20]})" for g in imp_df['gene']]
ax.set_yticklabels(gene_labels_imp, fontsize=9)
ax.set_xlabel("Feature Importance (Random Forest)")
ax.set_title("RF Feature Importance\n(Top 20 Biomarkers)")

# 6. ROC curves for panel classifier
ax = fig.add_subplot(gs[1, 2])
for panel_size in [3, 5, 10, 20]:
    genes_p = bm_df.head(panel_size)['gene'].tolist()
    X_p = log_cpm.loc[genes_p].T.values

    lr = LogisticRegression(max_iter=1000, random_state=42)
    lr.fit(X_p, labels)
    prob = lr.predict_proba(X_p)[:, 1]
    fpr, tpr, _ = roc_curve(labels, prob)
    auc_val = roc_auc_score(labels, prob)
    ax.plot(fpr, tpr, linewidth=1.5, label=f'{panel_size}-gene panel (AUC={auc_val:.3f})')

ax.plot([0, 1], [0, 1], 'k--', linewidth=0.5)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("Panel Classifier ROC Curves\n(Logistic Regression)")
ax.legend(fontsize=12)

# 7-8-9: Category-specific analysis
categories_info = [
    ('Universal Biomarkers', universal.head(15), '#2ecc71'),
    ('Early-Detection', early_biomarkers.head(15), '#3498db'),
    ('Infection-Age', age_biomarkers.head(15), '#e67e22')
]

for c_idx, (cat_name, cat_df, cat_color) in enumerate(categories_info):
    ax = fig.add_subplot(gs[2, c_idx])

    if len(cat_df) == 0:
        ax.text(0.5, 0.5, f'No {cat_name}\nbiomarkers found', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(cat_name)
        continue

    cat_genes = cat_df['gene'].tolist()[:12]
    y = np.arange(len(cat_genes))

    aucs_cat = [cat_df[cat_df['gene'] == g]['AUC'].values[0] for g in cat_genes]
    fcs_cat = [cat_df[cat_df['gene'] == g]['mean_abs_fc'].values[0] for g in cat_genes]

    ax.barh(y, aucs_cat, color=cat_color, alpha=0.85)
    ax.set_yticks(y)
    gene_labels_cat = [f"{g} ({gene_name_map.get(g, '')[:22]})" for g in cat_genes]
    ax.set_yticklabels(gene_labels_cat, fontsize=8)
    ax.set_xlabel("AUC-ROC")
    ax.set_title(f"{cat_name}\n(top {len(cat_genes)} genes)")
    ax.axvline(0.7, color='grey', linestyle='--', linewidth=0.5)
    ax.invert_yaxis()

    # Add FC annotation
    for i, (a, f) in enumerate(zip(aucs_cat, fcs_cat)):
        ax.text(a + 0.01, i, f'FC={f:.1f}', va='center', fontsize=8)

fig.suptitle("GRBV Biomarker Robustness & Validation", fontsize=22, y=1.02)
plt.savefig('22_biomarker_robustness.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 22_biomarker_robustness.png")


# ============================================================
# FINAL OUTPUT: BIOMARKER TABLE
# ============================================================
print("\n" + "=" * 70)
print("FINAL BIOMARKER TABLE (Top 50)")
print("=" * 70)

# Build comprehensive output table
output_cols = ['gene', 'gene_name', 'primary_direction', 'biomarker_score_final', 'AUC',
               'n_conditions_sig', 'n_stages_sig', 'n_ages_sig', 'season_robust',
               'mean_abs_fc', 'min_padj']
output_df = bm_df[output_cols].head(50).copy()
output_df.columns = ['Gene_ID', 'Gene_Name', 'Direction', 'Biomarker_Score', 'AUC_ROC',
                      'N_Conditions_Sig', 'N_Stages_Sig', 'N_Ages_Sig', 'Season_Robust',
                      'Mean_|log2FC|', 'Min_padj']

# Add GO annotation summary
go_annotations = []
kegg_annotations = []
for gene in output_df['Gene_ID']:
    go_terms = gene2go.get(gene, set())
    bp_terms = [t for t in go_terms if '(BP)' in t]
    go_str = '; '.join(sorted(bp_terms)[:3]) if bp_terms else ''
    go_annotations.append(go_str)

    kegg_terms = gene2kegg_set.get(gene, set())
    kegg_str = '; '.join([kegg_names.get(k, k) for k in sorted(kegg_terms)][:2]) if kegg_terms else ''
    kegg_annotations.append(kegg_str)

output_df['GO_BP'] = go_annotations
output_df['KEGG'] = kegg_annotations

# Save to file
output_df.to_csv('GRBV_biomarkers_top50.csv', index=False)
print(f"\nSaved: GRBV_biomarkers_top50.csv")

# Print
print(f"\n{'='*140}")
print(f"{'#':<4} {'Gene':<18} {'Name':<35} {'Dir':<5} {'Score':<7} {'AUC':<6} "
      f"{'Cond':<5} {'Stg':<4} {'Age':<4} {'Ssn':<5} {'FC':<6} {'padj':<10}")
print(f"{'='*140}")

for i, (_, row) in enumerate(output_df.iterrows(), 1):
    name = str(row['Gene_Name'])[:32] if pd.notna(row['Gene_Name']) else ''
    ssn = 'Both' if row['Season_Robust'] else '-'
    print(f"{i:<4} {row['Gene_ID']:<18} {name:<35} {row['Direction']:<5} "
          f"{row['Biomarker_Score']:.4f} {row['AUC_ROC']:.3f} "
          f"{row['N_Conditions_Sig']:<5} {row['N_Stages_Sig']:<4} {row['N_Ages_Sig']:<4} "
          f"{ssn:<5} {row['Mean_|log2FC|']:.2f}  {row['Min_padj']:.1e}")

# Summary statistics
print(f"\n\n{'='*70}")
print("BIOMARKER DISCOVERY SUMMARY")
print(f"{'='*70}")

top50_summary = output_df.head(50)
print(f"\nTop 50 Biomarkers:")
print(f"  Upregulated: {(top50_summary['Direction'] == 'UP').sum()}")
print(f"  Downregulated: {(top50_summary['Direction'] == 'DOWN').sum()}")
print(f"  Season-robust: {top50_summary['Season_Robust'].sum()}")
print(f"  Mean AUC: {top50_summary['AUC_ROC'].mean():.3f}")
print(f"  Median conditions significant: {top50_summary['N_Conditions_Sig'].median():.0f}")

print(f"\nMulti-gene panel performance (5-fold CV):")
for ps in panel_sizes:
    pr = panel_results[ps]
    print(f"  {ps:2d}-gene panel: LR AUC={pr['lr_mean']:.4f}±{pr['lr_std']:.4f}  |  "
          f"RF AUC={pr['rf_mean']:.4f}±{pr['rf_std']:.4f}")

print(f"\nBiomarker Categories:")
print(f"  Universal (season+stage robust, AUC≥0.7): {len(universal)}")
print(f"  Early-detection (Pre_ver/Veraison, AUC≥0.65): {len(early_biomarkers)}")
print(f"  Harvest-specific (AUC≥0.7): {len(harvest_specific)}")
print(f"  Infection-age responsive (≥3 ages, AUC≥0.7): {len(age_biomarkers)}")

print(f"\n{'='*70}")
print("GENERATED FILES:")
print("  20_biomarker_overview.png        — Score distribution, AUC, panels, volcano, categories")
print("  21_biomarker_profiles.png        — Top 16 gene expression profiles across Season×Stage")
print("  22_biomarker_robustness.png      — FC consistency, AUC stability, significance matrix")
print("  GRBV_biomarkers_top50.csv        — Table of top 50 biomarkers with annotations")
print(f"{'='*70}")
