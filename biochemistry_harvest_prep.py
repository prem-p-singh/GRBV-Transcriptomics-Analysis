"""
Biochemistry x transcriptome prep (Harvest).

Reads Biochemical_Data.xlsx, tidies it, merges with normalized expression for
SPS1, MBF1c, LOC100266309, plus an sHSP module score and viral titer, then
computes three stats tables:

  (1) phenotypic validation correlations (global / partial-season / per-season)
  (2) cohort trajectory linear regression (biochem ~ YSI_years, per season)
  (3) Season x Infection two-way ANOVA + season-specific GRBV contrasts

No figures are generated here. See figures/scripts_R/supp_biochem.R for the
S61-S63 publication figures.

Outputs (all in repo root):
  biochemistry_harvest.csv
  biochemistry_transcriptome_harvest_merged.csv
  biochemistry_validation_correlations.csv
  biochemistry_cohort_trajectory.csv
  biochemistry_season_infection_anova.csv
"""

import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

ROOT = Path('/Users/prempratapsingh/Data_Analysis')

TARGET_GENES = {
    'SPS1': 'LOC100232974',
    'MBF1c': 'LOC100249249',
    'LOC100266309': 'LOC100266309',
}
BIOCHEM_VARS = ['Brix', 'pH', 'TA', 'Malic', 'NOPA', 'NH3', 'YAN']


def load_biochem() -> pd.DataFrame:
    frames = []
    for year in (2021, 2022):
        df = pd.read_excel(ROOT / 'Biochemical_Data.xlsx', sheet_name=str(year))
        df.columns = [str(c).strip() for c in df.columns]
        df = df.rename(columns={
            'Sample harvest': 'Sample',
            'NOPA (mg/L)': 'NOPA', 'NH3 (mg/L)': 'NH3', 'YAN (mg/L)': 'YAN',
        })
        df = df[df['Sample'].astype(str).str.match(r'^2[12]ME03PASO')]
        df['Sampling_Season'] = year
        frames.append(df)
    b = pd.concat(frames, ignore_index=True)

    pat = re.compile(r'^2[12]ME03PASO_RB\(([-+])\)(\d{2})?_(\d)$')
    info = b['Sample'].apply(lambda s: pat.match(s).groups())
    b['Infection_Group'] = info.apply(lambda t: 'RB(+)' if t[0] == '+' else 'RB(-)')
    b['Cohort'] = info.apply(lambda t: t[1])
    b['Biorep'] = info.apply(lambda t: int(t[2]))
    b['Infection_Year'] = pd.to_numeric('20' + b['Cohort'].fillna(''), errors='coerce')
    b['YSI_years'] = np.where(b['Infection_Group'] == 'RB(+)',
                              b['Sampling_Season'] - b['Infection_Year'],
                              np.nan)
    for c in BIOCHEM_VARS + ['Tare (Kg)', 'Total Weight (Kg)', 'Volume (Liter)', 'Conversion']:
        if c in b.columns:
            b[c] = pd.to_numeric(b[c], errors='coerce')
    keep = ['Sample', 'Sampling_Season', 'Infection_Group', 'Cohort',
            'Infection_Year', 'YSI_years', 'Biorep',
            'Total Weight (Kg)', 'Volume (Liter)'] + BIOCHEM_VARS
    return b[keep].rename(columns={'Total Weight (Kg)': 'Cluster_Weight_kg',
                                   'Volume (Liter)': 'Juice_Volume_L'})


def harvest_expression(meta_h: pd.DataFrame) -> pd.DataFrame:
    counts = pd.read_csv(ROOT / 'counts_matrix.txt', sep='\t', index_col=0)
    h_samples = [s for s in meta_h['Sample'] if s in counts.columns]
    sub = counts[h_samples]
    cpm = sub.divide(sub.sum(axis=0), axis=1) * 1e6
    log2cpm = np.log2(cpm + 1)

    expr = {}
    for label, gid in TARGET_GENES.items():
        if gid in log2cpm.index:
            expr[label] = log2cpm.loc[gid]

    hs = pd.read_csv(ROOT / 'heat_stress_DEGs_all_conditions.csv')
    shsp_ids = [g for g in hs[hs['Category'] == 'sHSP (small heat shock proteins)']['Gene'].unique()
                if g in log2cpm.index]
    if shsp_ids:
        m = log2cpm.loc[shsp_ids]
        z = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0)
        expr['sHSP_score'] = z.mean(axis=0)
    e = pd.DataFrame(expr)
    e.index.name = 'Sample'
    return e.reset_index()


def titer_harvest() -> pd.DataFrame:
    t = pd.read_csv(ROOT / 'Virus titer.txt', sep='\t')
    t.columns = [c.strip() for c in t.columns]
    harv_col = [c for c in t.columns if 'Harvest' in c][0]
    out = t[['Season', 'Year', 'Biorep', harv_col]].rename(columns={harv_col: 'Titer_Harvest'})
    out['Titer_log10'] = np.log10(out['Titer_Harvest'].clip(lower=1))
    return out


def merge_biochem_transcriptome(biochem: pd.DataFrame) -> pd.DataFrame:
    meta = pd.read_csv(ROOT / 'sample_metadata_v2.txt', sep='\t')
    meta_h = meta[meta['Maturity_stage'] == 'Harvest'].copy()
    meta_h['Cohort'] = meta_h['Years_Since_Infection'].astype(str).str[-2:]
    meta_h.loc[meta_h['Infection_Group'] == 'RB(-)', 'Cohort'] = None

    expr = harvest_expression(meta_h)
    meta_expr = meta_h.merge(expr, on='Sample', how='left')

    titer = titer_harvest()
    titer = titer.rename(columns={'Season': 'Sampling_Season', 'Year': 'Infection_Year'})
    meta_expr['Infection_Year'] = pd.to_numeric(meta_expr['Years_Since_Infection'], errors='coerce')
    meta_expr = meta_expr.merge(
        titer[['Sampling_Season', 'Infection_Year', 'Biorep', 'Titer_Harvest', 'Titer_log10']],
        left_on=['Sampling_Season', 'Infection_Year', 'biorep_id'],
        right_on=['Sampling_Season', 'Infection_Year', 'Biorep'],
        how='left').drop(columns=['Biorep'])

    merged = biochem.merge(
        meta_expr,
        left_on=['Sampling_Season', 'Infection_Group', 'Cohort', 'Biorep'],
        right_on=['Sampling_Season', 'Infection_Group', 'Cohort', 'biorep_id'],
        how='left', suffixes=('', '_tx'))
    return merged


def validation_correlations(m: pd.DataFrame) -> pd.DataFrame:
    predictors = [c for c in ['SPS1', 'MBF1c', 'LOC100266309', 'sHSP_score', 'Titer_log10']
                  if c in m.columns]
    rows = []
    for p in predictors:
        for v in BIOCHEM_VARS:
            sub = m[[p, v, 'Sampling_Season']].dropna()
            if len(sub) < 6:
                continue
            rho_g, p_g = stats.spearmanr(sub[p], sub[v])
            per_season = {}
            for season in (2021, 2022):
                ss = sub[sub['Sampling_Season'] == season]
                if len(ss) >= 4:
                    r, pv = stats.spearmanr(ss[p], ss[v])
                    per_season[season] = (r, pv, len(ss))
                else:
                    per_season[season] = (np.nan, np.nan, len(ss))
            xs = sub[p].values; ys = sub[v].values
            s = (sub['Sampling_Season'] == 2022).astype(float).values
            xr = xs - np.polyval(np.polyfit(s, xs, 1), s)
            yr = ys - np.polyval(np.polyfit(s, ys, 1), s)
            rho_p, p_p = stats.spearmanr(xr, yr)
            rows.append({'predictor': p, 'biochem': v, 'n_total': len(sub),
                         'rho_global': rho_g, 'p_global': p_g,
                         'rho_partial_season': rho_p, 'p_partial_season': p_p,
                         'rho_2021': per_season[2021][0], 'p_2021': per_season[2021][1], 'n_2021': per_season[2021][2],
                         'rho_2022': per_season[2022][0], 'p_2022': per_season[2022][1], 'n_2022': per_season[2022][2]})
    return pd.DataFrame(rows)


def cohort_trajectory(biochem: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for season, g_season in biochem.groupby('Sampling_Season'):
        rbp = g_season[g_season['Infection_Group'] == 'RB(+)']
        for v in BIOCHEM_VARS:
            sub = rbp[['YSI_years', v]].dropna()
            if len(sub) < 5:
                continue
            slope, intercept, r, p, se = stats.linregress(sub['YSI_years'], sub[v])
            rows.append({'season': season, 'biochem': v, 'n': len(sub),
                         'slope_per_year': slope, 'intercept': intercept,
                         'r': r, 'p': p, 'se': se})
    return pd.DataFrame(rows)


def season_infection_anova(biochem: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for v in BIOCHEM_VARS:
        sub = biochem[['Sampling_Season', 'Infection_Group', v]].dropna().copy()
        sub.columns = ['Season', 'Infection', 'y']
        sub['Season'] = sub['Season'].astype('category')
        sub['Infection'] = sub['Infection'].astype('category')
        model = ols('y ~ C(Season) * C(Infection)', data=sub).fit()
        aov = sm.stats.anova_lm(model, typ=2)

        def _r(term):
            if term not in aov.index:
                return dict(F=np.nan, p=np.nan, eta2=np.nan)
            F = aov.loc[term, 'F']; p = aov.loc[term, 'PR(>F)']
            ss = aov.loc[term, 'sum_sq']; sst = aov['sum_sq'].sum()
            return dict(F=F, p=p, eta2=ss / sst)

        s = _r('C(Season)'); i = _r('C(Infection)'); ix = _r('C(Season):C(Infection)')

        def grbv_effect(season):
            d = sub[sub['Season'] == season]
            a = d[d['Infection'] == 'RB(+)']['y']
            b = d[d['Infection'] == 'RB(-)']['y']
            if len(a) < 2 or len(b) < 2:
                return np.nan, np.nan
            t, p = stats.ttest_ind(a, b, equal_var=False)
            return float(a.mean() - b.mean()), float(p)

        d21, p21 = grbv_effect(2021); d22, p22 = grbv_effect(2022)
        rows.append({'biochem': v, 'n': len(sub),
                     'F_season': s['F'], 'p_season': s['p'], 'eta2_season': s['eta2'],
                     'F_infection': i['F'], 'p_infection': i['p'], 'eta2_infection': i['eta2'],
                     'F_interaction': ix['F'], 'p_interaction': ix['p'], 'eta2_interaction': ix['eta2'],
                     'GRBV_effect_2021': d21, 'p_GRBV_2021': p21,
                     'GRBV_effect_2022': d22, 'p_GRBV_2022': p22})
    return pd.DataFrame(rows)


def main():
    biochem = load_biochem()
    biochem.to_csv(ROOT / 'biochemistry_harvest.csv', index=False)
    print(f'[1] biochemistry_harvest.csv : {len(biochem)} rows')

    merged = merge_biochem_transcriptome(biochem)
    merged.to_csv(ROOT / 'biochemistry_transcriptome_harvest_merged.csv', index=False)
    n_with = merged.dropna(subset=['SPS1']).shape[0] if 'SPS1' in merged.columns else 0
    print(f'[2] biochemistry_transcriptome_harvest_merged.csv : {len(merged)} rows, {n_with} with expression')

    val = validation_correlations(merged)
    val.to_csv(ROOT / 'biochemistry_validation_correlations.csv', index=False)
    print(f'[3] biochemistry_validation_correlations.csv : {len(val)} pairs')

    traj = cohort_trajectory(biochem)
    traj.to_csv(ROOT / 'biochemistry_cohort_trajectory.csv', index=False)
    print(f'[4] biochemistry_cohort_trajectory.csv : {len(traj)} rows')

    aov = season_infection_anova(biochem)
    aov.to_csv(ROOT / 'biochemistry_season_infection_anova.csv', index=False)
    print(f'[5] biochemistry_season_infection_anova.csv : {len(aov)} rows')


if __name__ == '__main__':
    main()
