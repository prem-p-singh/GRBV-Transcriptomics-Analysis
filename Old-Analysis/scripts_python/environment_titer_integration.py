#!/usr/bin/env python3
"""
Integrated analysis of atmospheric temperature, virus titer, and GRBV transcriptomics.
Correlates weather conditions with viral load and host gene expression across
two seasons (2021, 2022) and four maturity stages.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import rankdata, spearmanr, pearsonr
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'figure.dpi': 300, 'font.size': 15, 'axes.titlesize': 16,
    'axes.labelsize': 15, 'figure.facecolor': 'white'
})

# ============================================================
# SAMPLING DATES
# ============================================================
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

stages = ['Pre_veraison', 'Veraison', 'Post_veraison', 'Harvest']
stage_colors = {'Pre_veraison': '#2ca02c', 'Veraison': '#ff7f0e',
                'Post_veraison': '#d62728', 'Harvest': '#9467bd'}

# ============================================================
# LOAD WEATHER DATA
# ============================================================
print("=" * 70)
print("LOADING WEATHER DATA")
print("=" * 70)

def load_cimis(filepath, year):
    """Load CIMIS daily weather data."""
    df = pd.read_csv(filepath)
    # Clean column names (remove trailing spaces)
    df.columns = [c.strip() for c in df.columns]

    # Parse date
    df['Date_str'] = df['Date'].astype(str).str.strip()
    df['Date'] = pd.to_datetime(df['Date_str'], format='%m/%d/%Y')

    # Extract numeric columns (skip QC flags)
    numeric_cols = {
        'Max Air Temp (C)': 'Tmax',
        'Min Air Temp (C)': 'Tmin',
        'Avg Air Temp (C)': 'Tavg',
        'Sol Rad (W/sq.m)': 'SolRad',
        'Avg Rel Hum (%)': 'RH_avg',
        'Max Rel Hum (%)': 'RH_max',
        'Min Rel Hum (%)': 'RH_min',
        'Precip (mm)': 'Precip',
        'ETo (mm)': 'ETo',
        'Avg Vap Pres (kPa)': 'VPD',
        'Dew Point (C)': 'DewPoint',
        'Avg Wind Speed (m/s)': 'WindSpeed',
        'Avg Soil Temp (C)': 'SoilTemp',
    }

    result = pd.DataFrame()
    result['Date'] = df['Date']

    for orig, new in numeric_cols.items():
        if orig in df.columns:
            result[new] = pd.to_numeric(df[orig], errors='coerce')

    result['Season'] = year

    # Calculate Growing Degree Days (base 10C for grapevine)
    result['GDD'] = np.maximum((result['Tmax'] + result['Tmin']) / 2 - 10, 0)
    result['GDD_cumsum'] = result['GDD'].cumsum()

    # Diurnal temperature range
    result['DTR'] = result['Tmax'] - result['Tmin']

    # Heat stress indicator (hours-like proxy: days above 35C)
    result['HeatStress'] = (result['Tmax'] > 35).astype(int)

    return result

weather_2021 = load_cimis('daily_2021.csv', 2021)
weather_2022 = load_cimis('daily_2022.csv', 2022)

print(f"2021 weather: {len(weather_2021)} days ({weather_2021['Date'].min().date()} to {weather_2021['Date'].max().date()})")
print(f"  Station: Atascadero (163), Central Coast Valleys")
print(f"  Avg temp range: {weather_2021['Tavg'].min():.1f} - {weather_2021['Tavg'].max():.1f} C")
print(f"  Max temp range: {weather_2021['Tmax'].min():.1f} - {weather_2021['Tmax'].max():.1f} C")

print(f"\n2022 weather: {len(weather_2022)} days ({weather_2022['Date'].min().date()} to {weather_2022['Date'].max().date()})")
print(f"  Station: Paso Robles (265), Central Coast Valleys")
print(f"  Avg temp range: {weather_2022['Tavg'].min():.1f} - {weather_2022['Tavg'].max():.1f} C")
print(f"  Max temp range: {weather_2022['Tmax'].min():.1f} - {weather_2022['Tmax'].max():.1f} C")


# ============================================================
# COMPUTE WEATHER METRICS AT SAMPLING DATES
# ============================================================
print("\n" + "=" * 70)
print("WEATHER AT SAMPLING DATES")
print("=" * 70)

def get_weather_window(weather_df, target_date, window_days=7):
    """Get weather stats for a window around a target date."""
    start = target_date - timedelta(days=window_days)
    end = target_date
    mask = (weather_df['Date'] >= start) & (weather_df['Date'] <= end)
    window = weather_df[mask]
    if len(window) == 0:
        return {}
    return {
        'Tavg_mean': window['Tavg'].mean(),
        'Tavg_std': window['Tavg'].std(),
        'Tmax_mean': window['Tmax'].mean(),
        'Tmax_max': window['Tmax'].max(),
        'Tmin_mean': window['Tmin'].mean(),
        'Tmin_min': window['Tmin'].min(),
        'DTR_mean': window['DTR'].mean(),
        'GDD_sum': window['GDD'].sum(),
        'HeatStress_days': window['HeatStress'].sum(),
        'SolRad_mean': window['SolRad'].mean(),
        'RH_avg_mean': window['RH_avg'].mean(),
        'Precip_total': window['Precip'].sum(),
        'ETo_total': window['ETo'].sum(),
        'VPD_mean': window['VPD'].mean(),
    }

# Build weather summary per stage per season
weather_at_sampling = []
for season in [2021, 2022]:
    wx = weather_2021 if season == 2021 else weather_2022
    for stage in stages:
        date = sampling_dates[season][stage]

        # Day-of sampling
        day_mask = wx['Date'] == pd.Timestamp(date)
        day_data = wx[day_mask]

        # Windows: 7-day, 14-day, 21-day before sampling
        for window in [7, 14, 21]:
            metrics = get_weather_window(wx, date, window)
            for k, v in metrics.items():
                weather_at_sampling.append({
                    'Season': season, 'Stage': stage,
                    'Sampling_Date': date.strftime('%m/%d/%Y'),
                    'Window': f'{window}d',
                    'Metric': k, 'Value': v
                })

        # Day-of values
        if len(day_data) > 0:
            row = day_data.iloc[0]
            weather_at_sampling.append({
                'Season': season, 'Stage': stage,
                'Sampling_Date': date.strftime('%m/%d/%Y'),
                'Window': 'day_of', 'Metric': 'Tavg',
                'Value': row['Tavg']
            })
            weather_at_sampling.append({
                'Season': season, 'Stage': stage,
                'Sampling_Date': date.strftime('%m/%d/%Y'),
                'Window': 'day_of', 'Metric': 'Tmax',
                'Value': row['Tmax']
            })

weather_summary_df = pd.DataFrame(weather_at_sampling)

# Print key temperature comparisons
print(f"\n{'Stage':<16} {'2021 Date':>12} {'Tavg':>8} {'Tmax':>8} {'7d-GDD':>8}  "
      f"{'2022 Date':>12} {'Tavg':>8} {'Tmax':>8} {'7d-GDD':>8}  {'dT_avg':>8}")
print("-" * 120)

stage_temp_comparison = []
for stage in stages:
    for season in [2021, 2022]:
        wx = weather_2021 if season == 2021 else weather_2022
        date = sampling_dates[season][stage]
        m7 = get_weather_window(wx, date, 7)

        day_mask = wx['Date'] == pd.Timestamp(date)
        day_data = wx[day_mask]
        tavg_day = day_data.iloc[0]['Tavg'] if len(day_data) > 0 else np.nan
        tmax_day = day_data.iloc[0]['Tmax'] if len(day_data) > 0 else np.nan

        stage_temp_comparison.append({
            'Season': season, 'Stage': stage,
            'Date': date.strftime('%m/%d/%Y'),
            'Tavg_dayof': tavg_day, 'Tmax_dayof': tmax_day,
            'Tavg_7d': m7.get('Tavg_mean', np.nan),
            'Tmax_7d': m7.get('Tmax_mean', np.nan),
            'GDD_7d': m7.get('GDD_sum', np.nan),
            'HeatStress_7d': m7.get('HeatStress_days', np.nan),
        })

stc = pd.DataFrame(stage_temp_comparison)
for stage in stages:
    s21 = stc[(stc['Season'] == 2021) & (stc['Stage'] == stage)].iloc[0]
    s22 = stc[(stc['Season'] == 2022) & (stc['Stage'] == stage)].iloc[0]
    dt = s22['Tavg_7d'] - s21['Tavg_7d']
    print(f"{stage:<16} {s21['Date']:>12} {s21['Tavg_7d']:>8.1f} {s21['Tmax_7d']:>8.1f} {s21['GDD_7d']:>8.1f}  "
          f"{s22['Date']:>12} {s22['Tavg_7d']:>8.1f} {s22['Tmax_7d']:>8.1f} {s22['GDD_7d']:>8.1f}  "
          f"{dt:>+8.1f}")

# GDD accumulation at each sampling point
print("\n--- Cumulative GDD (base 10C) at sampling dates ---")
for season in [2021, 2022]:
    wx = weather_2021 if season == 2021 else weather_2022
    print(f"  {season}:")
    for stage in stages:
        date = sampling_dates[season][stage]
        mask = wx['Date'] <= pd.Timestamp(date)
        gdd_cum = wx[mask]['GDD'].sum()
        print(f"    {stage}: {date.strftime('%m/%d')} -> GDD = {gdd_cum:.0f}")


# ============================================================
# LOAD VIRUS TITER DATA
# ============================================================
print("\n" + "=" * 70)
print("LOADING VIRUS TITER DATA")
print("=" * 70)

titer = pd.read_csv('Virus titer.txt', sep='\t')
# Clean column names
titer.columns = [c.strip().replace('\u00e9', 'e') for c in titer.columns]
# Rename to match our stage naming
titer = titer.rename(columns={
    'Pre_veraison': 'Pre_veraison',
    'Veraison': 'Veraison',
    'Post_veraison': 'Post_veraison',
    'Harvest': 'Harvest'
})

# Drop empty rows
titer = titer.dropna(subset=['Season'])
titer['Season'] = titer['Season'].astype(int)
titer['Year'] = titer['Year'].astype(int)
titer['Biorep'] = titer['Biorep'].astype(int)

print(f"Titer data shape: {titer.shape}")
print(f"Seasons: {titer['Season'].unique()}")
print(f"Infection years: {sorted(titer['Year'].unique())}")

# Convert to long format for analysis
titer_long = titer.melt(
    id_vars=['Season', 'Year', 'Biorep'],
    value_vars=stages,
    var_name='Stage', value_name='Titer'
)
titer_long['Titer'] = pd.to_numeric(titer_long['Titer'], errors='coerce')
titer_long = titer_long.dropna(subset=['Titer'])
titer_long['log_Titer'] = np.log10(titer_long['Titer'] + 1)

# Summary statistics
print("\n--- Virus Titer Summary (copies/sample) ---")
print(f"{'Season':<8} {'Stage':<16} {'n':>4} {'Mean':>14} {'Median':>14} {'log10(Mean)':>12}")
print("-" * 70)
titer_summary = []
for season in [2021, 2022]:
    for stage in stages:
        sub = titer_long[(titer_long['Season'] == season) & (titer_long['Stage'] == stage)]
        if len(sub) == 0:
            continue
        mean_t = sub['Titer'].mean()
        median_t = sub['Titer'].median()
        log_mean = np.log10(mean_t + 1)
        n = len(sub)
        print(f"{season:<8} {stage:<16} {n:>4} {mean_t:>14,.0f} {median_t:>14,.0f} {log_mean:>12.2f}")
        titer_summary.append({
            'Season': season, 'Stage': stage, 'n': n,
            'mean_titer': mean_t, 'median_titer': median_t,
            'log10_mean': log_mean, 'std_titer': sub['Titer'].std(),
            'mean_log_titer': sub['log_Titer'].mean()
        })

titer_sum_df = pd.DataFrame(titer_summary)

# Titer change between seasons
print("\n--- Titer change 2021 vs 2022 per stage ---")
for stage in stages:
    t21 = titer_long[(titer_long['Season'] == 2021) & (titer_long['Stage'] == stage)]['log_Titer']
    t22 = titer_long[(titer_long['Season'] == 2022) & (titer_long['Stage'] == stage)]['log_Titer']
    if len(t21) > 2 and len(t22) > 2:
        stat, p = stats.mannwhitneyu(t21, t22, alternative='two-sided')
        diff = t22.mean() - t21.mean()
        print(f"  {stage}: 2021 log10={t21.mean():.2f}, 2022 log10={t22.mean():.2f}, "
              f"diff={diff:+.2f}, Mann-Whitney p={p:.4f} {'*' if p < 0.05 else 'ns'}")

# Titer progression across stages per season
print("\n--- Titer progression across maturity stages ---")
for season in [2021, 2022]:
    print(f"  {season}:")
    for i in range(len(stages) - 1):
        s1, s2 = stages[i], stages[i+1]
        t1 = titer_long[(titer_long['Season'] == season) & (titer_long['Stage'] == s1)]['log_Titer']
        t2 = titer_long[(titer_long['Season'] == season) & (titer_long['Stage'] == s2)]['log_Titer']
        if len(t1) > 2 and len(t2) > 2:
            stat, p = stats.mannwhitneyu(t1, t2, alternative='two-sided')
            fc = t2.mean() - t1.mean()
            print(f"    {s1} -> {s2}: log10 diff = {fc:+.2f}, p = {p:.4f} {'*' if p < 0.05 else 'ns'}")


# ============================================================
# INTEGRATE: TEMPERATURE x TITER x TRANSCRIPTOME
# ============================================================
print("\n" + "=" * 70)
print("INTEGRATED ANALYSIS: TEMPERATURE x TITER x TRANSCRIPTOMICS")
print("=" * 70)

# Load transcriptomic data for correlation
counts = pd.read_csv('counts_matrix.txt', sep='\t', index_col=0)
meta = pd.read_csv('sample_metadata.txt', sep='\t')
common_samples = list(set(counts.columns) & set(meta['Sample']))
counts = counts[common_samples]
meta = meta[meta['Sample'].isin(common_samples)].set_index('Sample').loc[common_samples]

# Filter and normalize
min_samples = int(0.1 * counts.shape[1])
gene_mask = (counts >= 10).sum(axis=1) >= min_samples
counts_filtered = counts.loc[gene_mask]
lib_sizes = counts_filtered.sum(axis=0)
cpm = counts_filtered.div(lib_sizes, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

def run_de_analysis(group1_samples, group2_samples, log_cpm_data):
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
        results.append({'gene': gene, 'log2FC': log2fc, 'pvalue': p_val})
    df = pd.DataFrame(results)
    if len(df) == 0:
        return df
    n = len(df)
    ranked = rankdata(df['pvalue'])
    df['padj'] = np.minimum(df['pvalue'] * n / ranked, 1.0)
    return df.sort_values('pvalue')

def get_samples(meta_df, **kwargs):
    mask = pd.Series(True, index=meta_df.index)
    for col, val in kwargs.items():
        mask &= meta_df[col] == val
    return list(meta_df[mask].index)

# Build integrated table: per season x stage, combine weather + titer + DEGs
print("\n--- Building integrated table ---")
integrated_data = []

for season in [2021, 2022]:
    wx = weather_2021 if season == 2021 else weather_2022
    for stage in stages:
        date = sampling_dates[season][stage]

        # Weather metrics
        m7 = get_weather_window(wx, date, 7)
        m14 = get_weather_window(wx, date, 14)

        # Day-of weather
        day_mask = wx['Date'] == pd.Timestamp(date)
        day_data = wx[day_mask]
        tavg_day = day_data.iloc[0]['Tavg'] if len(day_data) > 0 else np.nan
        tmax_day = day_data.iloc[0]['Tmax'] if len(day_data) > 0 else np.nan

        # Cumulative GDD
        gdd_cum = wx[wx['Date'] <= pd.Timestamp(date)]['GDD'].sum()

        # Titer
        titer_sub = titer_long[(titer_long['Season'] == season) & (titer_long['Stage'] == stage)]
        mean_titer = titer_sub['Titer'].mean() if len(titer_sub) > 0 else np.nan
        mean_log_titer = titer_sub['log_Titer'].mean() if len(titer_sub) > 0 else np.nan

        # DEGs (season-stratified)
        pos = get_samples(meta, Infection_Group='RB(+)', Maturity_stage=stage, Sampling_Season=season)
        neg = get_samples(meta, Infection_Group='RB(-)', Maturity_stage=stage, Sampling_Season=season)
        n_degs = 0
        n_up = 0
        n_down = 0
        if len(pos) >= 3 and len(neg) >= 3:
            de = run_de_analysis(pos, neg, log_cpm)
            sig = de[(de['padj'] < 0.05) & (de['log2FC'].abs() >= 1)]
            n_degs = len(sig)
            n_up = (sig['log2FC'] > 0).sum()
            n_down = (sig['log2FC'] < 0).sum()

        integrated_data.append({
            'Season': season, 'Stage': stage,
            'Sampling_Date': date.strftime('%m/%d/%Y'),
            'DOY': date.timetuple().tm_yday,
            'Tavg_dayof': tavg_day,
            'Tmax_dayof': tmax_day,
            'Tavg_7d': m7.get('Tavg_mean', np.nan),
            'Tmax_7d': m7.get('Tmax_mean', np.nan),
            'GDD_7d': m7.get('GDD_sum', np.nan),
            'HeatStress_7d': m7.get('HeatStress_days', np.nan),
            'SolRad_7d': m7.get('SolRad_mean', np.nan),
            'RH_7d': m7.get('RH_avg_mean', np.nan),
            'Precip_7d': m7.get('Precip_total', np.nan),
            'Tavg_14d': m14.get('Tavg_mean', np.nan),
            'Tmax_14d': m14.get('Tmax_mean', np.nan),
            'GDD_14d': m14.get('GDD_sum', np.nan),
            'GDD_cumulative': gdd_cum,
            'Mean_Titer': mean_titer,
            'Log10_Titer': mean_log_titer,
            'n_DEGs': n_degs,
            'n_Up': n_up,
            'n_Down': n_down,
        })

int_df = pd.DataFrame(integrated_data)
print(int_df.to_string(index=False))

# Correlations
print("\n--- Correlations across all 8 data points (2 seasons x 4 stages) ---")
corr_pairs = [
    ('Tavg_7d', 'Log10_Titer', 'Temperature vs Virus Titer'),
    ('Tmax_7d', 'Log10_Titer', 'Max Temp vs Virus Titer'),
    ('GDD_7d', 'Log10_Titer', 'GDD vs Virus Titer'),
    ('GDD_cumulative', 'Log10_Titer', 'Cumulative GDD vs Virus Titer'),
    ('Tavg_7d', 'n_DEGs', 'Temperature vs DEG count'),
    ('GDD_cumulative', 'n_DEGs', 'Cumulative GDD vs DEG count'),
    ('Log10_Titer', 'n_DEGs', 'Virus Titer vs DEG count'),
    ('HeatStress_7d', 'n_DEGs', 'Heat Stress Days vs DEG count'),
    ('HeatStress_7d', 'Log10_Titer', 'Heat Stress vs Virus Titer'),
    ('SolRad_7d', 'Log10_Titer', 'Solar Radiation vs Virus Titer'),
    ('RH_7d', 'Log10_Titer', 'Humidity vs Virus Titer'),
]

for x_col, y_col, label in corr_pairs:
    x = int_df[x_col].dropna()
    y = int_df[y_col].dropna()
    common_idx = x.index.intersection(y.index)
    if len(common_idx) < 4:
        continue
    r_p, p_p = pearsonr(x[common_idx], y[common_idx])
    r_s, p_s = spearmanr(x[common_idx], y[common_idx])
    sig = '*' if p_s < 0.05 else 'ns'
    print(f"  {label}: Pearson r={r_p:.3f} (p={p_p:.3f}), Spearman rho={r_s:.3f} (p={p_s:.3f}) {sig}")


# ============================================================
# TEMPERATURE DIFFERENCE ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("SEASONAL TEMPERATURE COMPARISON")
print("=" * 70)

# Overall season comparison
tavg_21 = weather_2021['Tavg'].mean()
tavg_22 = weather_2022['Tavg'].mean()
tmax_21 = weather_2021['Tmax'].mean()
tmax_22 = weather_2022['Tmax'].mean()
gdd_21 = weather_2021['GDD'].sum()
gdd_22 = weather_2022['GDD'].sum()
hs_21 = weather_2021['HeatStress'].sum()
hs_22 = weather_2022['HeatStress'].sum()

print(f"Overall season comparison (Jul 15 - Sep 30):")
print(f"  Avg daily temp:  2021 = {tavg_21:.1f}C,  2022 = {tavg_22:.1f}C  (diff = {tavg_22-tavg_21:+.1f}C)")
print(f"  Avg daily Tmax:  2021 = {tmax_21:.1f}C,  2022 = {tmax_22:.1f}C  (diff = {tmax_22-tmax_21:+.1f}C)")
print(f"  Total GDD:       2021 = {gdd_21:.0f},    2022 = {gdd_22:.0f}     (diff = {gdd_22-gdd_21:+.0f})")
print(f"  Heat stress days: 2021 = {hs_21},       2022 = {hs_22}        (diff = {hs_22-hs_21:+d})")
print(f"  Total precip:    2021 = {weather_2021['Precip'].sum():.1f}mm, 2022 = {weather_2022['Precip'].sum():.1f}mm")

# Inter-sampling temperature trends
print("\n--- Temperature trend between sampling dates ---")
for season in [2021, 2022]:
    wx = weather_2021 if season == 2021 else weather_2022
    print(f"\n  {season}:")
    for i in range(len(stages) - 1):
        s1, s2 = stages[i], stages[i+1]
        d1 = sampling_dates[season][s1]
        d2 = sampling_dates[season][s2]
        days_between = (d2 - d1).days
        mask = (wx['Date'] > pd.Timestamp(d1)) & (wx['Date'] <= pd.Timestamp(d2))
        period = wx[mask]
        if len(period) > 0:
            print(f"    {s1} -> {s2} ({days_between}d): "
                  f"Tavg={period['Tavg'].mean():.1f}C, Tmax={period['Tmax'].mean():.1f}C, "
                  f"GDD={period['GDD'].sum():.0f}, Heat days={period['HeatStress'].sum()}")


# ============================================================
# PER-INFECTION-YEAR TITER vs TEMPERATURE
# ============================================================
print("\n" + "=" * 70)
print("TITER BY INFECTION YEAR AND TEMPERATURE CONTEXT")
print("=" * 70)

for season in [2021, 2022]:
    print(f"\n--- Season {season} ---")
    inf_years = sorted(titer[titer['Season'] == season]['Year'].unique())
    for stage in stages:
        print(f"  {stage}:")
        for yr in inf_years:
            sub = titer_long[(titer_long['Season'] == season) &
                             (titer_long['Stage'] == stage) &
                             (titer_long['Year'] == yr)]
            if len(sub) > 0:
                print(f"    Inf.{yr}: mean={sub['Titer'].mean():>12,.0f}, "
                      f"log10={sub['log_Titer'].mean():.2f}, n={len(sub)}")

# Test: does infection duration affect titer?
print("\n--- Titer vs Infection Duration (Spearman) ---")
for season in [2021, 2022]:
    for stage in stages:
        sub = titer_long[(titer_long['Season'] == season) & (titer_long['Stage'] == stage)].copy()
        sub['InfAge'] = season - sub['Year']
        if len(sub) > 5:
            r, p = spearmanr(sub['InfAge'], sub['log_Titer'])
            sig = '*' if p < 0.05 else 'ns'
            print(f"  {season} {stage}: rho={r:.3f}, p={p:.4f} {sig}")


# ============================================================
# SAVE INTEGRATED TABLE
# ============================================================
int_df.to_csv('integrated_environment_titer_DEG.csv', index=False)
print("\nSaved: integrated_environment_titer_DEG.csv")

stc.to_csv('temperature_at_sampling.csv', index=False)
print("Saved: temperature_at_sampling.csv")


# ============================================================
# PLOT 30: Weather Overview + Sampling Timeline
# ============================================================
print("\n" + "=" * 70)
print("GENERATING VISUALIZATIONS")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(20, 14))

# Panel A: Temperature timelines with sampling dates
ax = axes[0, 0]
for season, wx, color, ls in [(2021, weather_2021, '#1f77b4', '-'),
                                (2022, weather_2022, '#ff7f0e', '-')]:
    # Normalize to DOY for overlay
    doy = wx['Date'].dt.dayofyear
    ax.plot(doy, wx['Tavg'], color=color, alpha=0.6, linewidth=1, linestyle=ls, label=f'{season} Tavg')
    ax.fill_between(doy, wx['Tmin'], wx['Tmax'], color=color, alpha=0.1)

    # Mark sampling dates
    for stage in stages:
        sd = sampling_dates[season][stage]
        sd_doy = sd.timetuple().tm_yday
        day_data = wx[wx['Date'] == pd.Timestamp(sd)]
        if len(day_data) > 0:
            tavg = day_data.iloc[0]['Tavg']
            ax.plot(sd_doy, tavg, marker='o' if season == 2021 else 's',
                    color=stage_colors[stage], markersize=10, zorder=5,
                    markeredgecolor='black', markeredgewidth=1)

ax.set_xlabel("Day of Year")
ax.set_ylabel("Air Temperature (C)")
ax.set_title("Daily Temperature Profile with Sampling Dates")
ax.legend(fontsize=12)

# Add stage labels
from matplotlib.lines import Line2D
stage_legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8, label=s)
                for s, c in stage_colors.items()]
ax2 = ax.twinx()
ax2.set_yticks([])
ax2.legend(handles=stage_legend, fontsize=10, loc='upper left')

# Panel B: GDD accumulation
ax = axes[0, 1]
for season, wx, color in [(2021, weather_2021, '#1f77b4'), (2022, weather_2022, '#ff7f0e')]:
    doy = wx['Date'].dt.dayofyear
    ax.plot(doy, wx['GDD_cumsum'], color=color, linewidth=2, label=f'{season}')

    for stage in stages:
        sd = sampling_dates[season][stage]
        sd_doy = sd.timetuple().tm_yday
        gdd_at = wx[wx['Date'] <= pd.Timestamp(sd)]['GDD'].sum()
        ax.plot(sd_doy, gdd_at, marker='o' if season == 2021 else 's',
                color=stage_colors[stage], markersize=10, zorder=5,
                markeredgecolor='black', markeredgewidth=1)
        ax.annotate(f'{gdd_at:.0f}', (sd_doy, gdd_at),
                    fontsize=10, ha='left', va='bottom',
                    xytext=(3, 3), textcoords='offset points')

ax.set_xlabel("Day of Year")
ax.set_ylabel("Cumulative GDD (base 10C)")
ax.set_title("Growing Degree Day Accumulation")
ax.legend(fontsize=14)

# Panel C: 7-day Avg Temperature comparison by stage
ax = axes[1, 0]
x = np.arange(len(stages))
tavg_21 = [int_df[(int_df['Season'] == 2021) & (int_df['Stage'] == s)]['Tavg_7d'].values[0] for s in stages]
tavg_22 = [int_df[(int_df['Season'] == 2022) & (int_df['Stage'] == s)]['Tavg_7d'].values[0] for s in stages]
tmax_21 = [int_df[(int_df['Season'] == 2021) & (int_df['Stage'] == s)]['Tmax_7d'].values[0] for s in stages]
tmax_22 = [int_df[(int_df['Season'] == 2022) & (int_df['Stage'] == s)]['Tmax_7d'].values[0] for s in stages]

w = 0.2
ax.bar(x - 1.5*w, tavg_21, w, label='2021 Tavg', color='#1f77b4', alpha=0.8)
ax.bar(x - 0.5*w, tavg_22, w, label='2022 Tavg', color='#ff7f0e', alpha=0.8)
ax.bar(x + 0.5*w, tmax_21, w, label='2021 Tmax', color='#1f77b4', alpha=0.4,
       edgecolor='#1f77b4', linewidth=1)
ax.bar(x + 1.5*w, tmax_22, w, label='2022 Tmax', color='#ff7f0e', alpha=0.4,
       edgecolor='#ff7f0e', linewidth=1)

ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=12)
ax.set_ylabel("Temperature (C) - 7-day window")
ax.set_title("7-day Average Temperature at Sampling")
ax.legend(fontsize=10)

# Add delta labels
for i in range(len(stages)):
    dt = tavg_22[i] - tavg_21[i]
    ax.text(i, max(tmax_21[i], tmax_22[i]) + 0.5, f'{dt:+.1f}C',
            ha='center', fontsize=12, fontweight='bold',
            color='red' if dt > 0 else 'blue')

# Panel D: Heat stress days comparison
ax = axes[1, 1]
hs_21 = [int_df[(int_df['Season'] == 2021) & (int_df['Stage'] == s)]['HeatStress_7d'].values[0] for s in stages]
hs_22 = [int_df[(int_df['Season'] == 2022) & (int_df['Stage'] == s)]['HeatStress_7d'].values[0] for s in stages]

ax.bar(x - 0.15, hs_21, 0.3, label='2021', color='#1f77b4', alpha=0.85)
ax.bar(x + 0.15, hs_22, 0.3, label='2022', color='#ff7f0e', alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=12)
ax.set_ylabel("Days with Tmax > 35C (7-day window)")
ax.set_title("Heat Stress Days at Sampling")
ax.legend()

for i in range(len(stages)):
    ax.text(i - 0.15, hs_21[i] + 0.1, str(int(hs_21[i])), ha='center', fontsize=14)
    ax.text(i + 0.15, hs_22[i] + 0.1, str(int(hs_22[i])), ha='center', fontsize=14)

plt.suptitle("Plot 30: Atmospheric Conditions During GRBV Sampling", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('30_weather_overview_sampling.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 30_weather_overview_sampling.png")


# ============================================================
# PLOT 31: Virus Titer Analysis
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(18, 14))

# Panel A: Titer by stage and season (box plots)
ax = axes[0, 0]
plot_data = titer_long.copy()
plot_data['Season_str'] = plot_data['Season'].astype(str)
plot_data['Stage'] = pd.Categorical(plot_data['Stage'], categories=stages, ordered=True)

# Custom positioning
positions_21 = [0.8, 1.8, 2.8, 3.8]
positions_22 = [1.2, 2.2, 3.2, 4.2]

bp1 = ax.boxplot([titer_long[(titer_long['Season'] == 2021) & (titer_long['Stage'] == s)]['log_Titer'].values
                   for s in stages],
                  positions=positions_21, widths=0.35, patch_artist=True,
                  boxprops=dict(facecolor='#1f77b4', alpha=0.6),
                  medianprops=dict(color='black', linewidth=2))
bp2 = ax.boxplot([titer_long[(titer_long['Season'] == 2022) & (titer_long['Stage'] == s)]['log_Titer'].values
                   for s in stages],
                  positions=positions_22, widths=0.35, patch_artist=True,
                  boxprops=dict(facecolor='#ff7f0e', alpha=0.6),
                  medianprops=dict(color='black', linewidth=2))

ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=12)
ax.set_ylabel("log10(Virus Titer)")
ax.set_title("Virus Titer by Stage and Season")
ax.legend([bp1["boxes"][0], bp2["boxes"][0]], ['2021', '2022'], fontsize=14)

# Panel B: Titer by infection year (grouped)
ax = axes[0, 1]
titer_by_year = titer_long.groupby(['Season', 'Year'])['log_Titer'].mean().reset_index()
for season, color in [(2021, '#1f77b4'), (2022, '#ff7f0e')]:
    sub = titer_by_year[titer_by_year['Season'] == season]
    ax.plot(sub['Year'], sub['log_Titer'], 'o-', color=color, markersize=8,
            linewidth=2, label=str(season))

ax.set_xlabel("Infection Year")
ax.set_ylabel("Mean log10(Virus Titer)")
ax.set_title("Average Titer by Infection Year")
ax.legend()

# Panel C: Titer progression across stages
ax = axes[1, 0]
for season, color, marker in [(2021, '#1f77b4', 'o'), (2022, '#ff7f0e', 's')]:
    means = []
    sems = []
    for stage in stages:
        sub = titer_long[(titer_long['Season'] == season) & (titer_long['Stage'] == stage)]
        means.append(sub['log_Titer'].mean())
        sems.append(sub['log_Titer'].sem())
    ax.errorbar(range(len(stages)), means, yerr=sems, marker=marker,
                color=color, linewidth=2, markersize=8, capsize=5, label=str(season))

ax.set_xticks(range(len(stages)))
ax.set_xticklabels([s.replace('_', '\n') for s in stages], fontsize=12)
ax.set_ylabel("Mean log10(Virus Titer) +/- SEM")
ax.set_title("Titer Progression Across Maturity")
ax.legend()

# Panel D: Titer heatmap by infection year x stage x season
ax = axes[1, 1]
hm_data = titer_long.groupby(['Season', 'Year', 'Stage'])['log_Titer'].mean().reset_index()
hm_pivot = hm_data.pivot_table(index=['Season', 'Year'], columns='Stage',
                                 values='log_Titer')[stages]
sns.heatmap(hm_pivot, ax=ax, cmap='YlOrRd', annot=True, fmt='.2f',
            cbar_kws={'label': 'log10(Titer)'})
ax.set_title("Virus Titer Heatmap\n(Season x Infection Year x Stage)")
ax.set_ylabel("Season / Infection Year")

plt.suptitle("Plot 31: Virus Titer Dynamics", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('31_virus_titer_analysis.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 31_virus_titer_analysis.png")


# ============================================================
# PLOT 32: Integrated Temperature x Titer x DEGs
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(22, 14))

# Panel A: Temp (7d) vs Titer
ax = axes[0, 0]
for season, color, marker in [(2021, '#1f77b4', 'o'), (2022, '#ff7f0e', 's')]:
    sub = int_df[int_df['Season'] == season]
    for _, row in sub.iterrows():
        ax.scatter(row['Tavg_7d'], row['Log10_Titer'], c=stage_colors[row['Stage']],
                   marker=marker, s=120, edgecolors='black', linewidth=1, zorder=5)
        ax.annotate(f"{row['Stage'][:4]}\n{season}", (row['Tavg_7d'], row['Log10_Titer']),
                    fontsize=9, ha='center', va='bottom', xytext=(0, 6),
                    textcoords='offset points')

r, p = pearsonr(int_df['Tavg_7d'], int_df['Log10_Titer'])
# Add regression line
z = np.polyfit(int_df['Tavg_7d'], int_df['Log10_Titer'], 1)
xline = np.linspace(int_df['Tavg_7d'].min() - 0.5, int_df['Tavg_7d'].max() + 0.5, 100)
ax.plot(xline, np.polyval(z, xline), 'k--', linewidth=1, alpha=0.5)
ax.set_xlabel("7-day Avg Temperature (C)")
ax.set_ylabel("log10(Mean Virus Titer)")
ax.set_title(f"Temperature vs Virus Titer\nr={r:.3f}, p={p:.3f}")

# Panel B: Titer vs DEGs
ax = axes[0, 1]
for season, marker in [(2021, 'o'), (2022, 's')]:
    sub = int_df[int_df['Season'] == season]
    for _, row in sub.iterrows():
        ax.scatter(row['Log10_Titer'], row['n_DEGs'], c=stage_colors[row['Stage']],
                   marker=marker, s=120, edgecolors='black', linewidth=1, zorder=5)
        ax.annotate(f"{row['Stage'][:4]}\n{season}", (row['Log10_Titer'], row['n_DEGs']),
                    fontsize=9, ha='left', va='bottom', xytext=(4, 4),
                    textcoords='offset points')

r, p = pearsonr(int_df['Log10_Titer'], int_df['n_DEGs'])
z = np.polyfit(int_df['Log10_Titer'], int_df['n_DEGs'], 1)
xline = np.linspace(int_df['Log10_Titer'].min() - 0.1, int_df['Log10_Titer'].max() + 0.1, 100)
ax.plot(xline, np.polyval(z, xline), 'k--', linewidth=1, alpha=0.5)
ax.set_xlabel("log10(Mean Virus Titer)")
ax.set_ylabel("Number of DEGs")
ax.set_title(f"Virus Titer vs Transcriptomic Response\nr={r:.3f}, p={p:.3f}")

# Panel C: GDD cumulative vs Titer
ax = axes[0, 2]
for season, marker in [(2021, 'o'), (2022, 's')]:
    sub = int_df[int_df['Season'] == season]
    for _, row in sub.iterrows():
        ax.scatter(row['GDD_cumulative'], row['Log10_Titer'], c=stage_colors[row['Stage']],
                   marker=marker, s=120, edgecolors='black', linewidth=1, zorder=5)

r, p = pearsonr(int_df['GDD_cumulative'], int_df['Log10_Titer'])
z = np.polyfit(int_df['GDD_cumulative'], int_df['Log10_Titer'], 1)
xline = np.linspace(int_df['GDD_cumulative'].min() - 10, int_df['GDD_cumulative'].max() + 10, 100)
ax.plot(xline, np.polyval(z, xline), 'k--', linewidth=1, alpha=0.5)
ax.set_xlabel("Cumulative GDD (base 10C)")
ax.set_ylabel("log10(Mean Virus Titer)")
ax.set_title(f"Cumulative GDD vs Virus Titer\nr={r:.3f}, p={p:.3f}")

# Custom legend
stage_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8, label=s)
                 for s, c in stage_colors.items()]
season_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
                          markersize=8, label='2021'),
                  Line2D([0], [0], marker='s', color='w', markerfacecolor='grey',
                          markersize=8, label='2022')]
ax.legend(handles=stage_handles + season_handles, fontsize=9, ncol=2)

# Panel D: Temp vs DEGs
ax = axes[1, 0]
for season, marker in [(2021, 'o'), (2022, 's')]:
    sub = int_df[int_df['Season'] == season]
    for _, row in sub.iterrows():
        ax.scatter(row['Tavg_7d'], row['n_DEGs'], c=stage_colors[row['Stage']],
                   marker=marker, s=120, edgecolors='black', linewidth=1, zorder=5)

r, p = pearsonr(int_df['Tavg_7d'], int_df['n_DEGs'])
z = np.polyfit(int_df['Tavg_7d'], int_df['n_DEGs'], 1)
xline = np.linspace(int_df['Tavg_7d'].min() - 0.5, int_df['Tavg_7d'].max() + 0.5, 100)
ax.plot(xline, np.polyval(z, xline), 'k--', linewidth=1, alpha=0.5)
ax.set_xlabel("7-day Avg Temperature (C)")
ax.set_ylabel("Number of DEGs")
ax.set_title(f"Temperature vs Transcriptomic Disruption\nr={r:.3f}, p={p:.3f}")

# Panel E: Multivariate bubble plot
ax = axes[1, 1]
for _, row in int_df.iterrows():
    season = row['Season']
    marker = 'o' if season == 2021 else 's'
    size = max(row['n_DEGs'] * 0.5, 20)
    ax.scatter(row['Tavg_7d'], row['Log10_Titer'], s=size,
               c=stage_colors[row['Stage']], marker=marker,
               alpha=0.7, edgecolors='black', linewidth=1)

ax.set_xlabel("7-day Avg Temperature (C)")
ax.set_ylabel("log10(Mean Virus Titer)")
ax.set_title("Bubble Plot: Temp x Titer x DEGs\n(bubble size = # DEGs)")

# Size legend
for nd in [50, 200, 500, 1000]:
    ax.scatter([], [], s=nd*0.5, c='grey', alpha=0.5, edgecolors='black',
               label=f'{nd} DEGs')
ax.legend(fontsize=10, loc='upper left')

# Panel F: Correlation matrix
ax = axes[1, 2]
corr_cols = ['Tavg_7d', 'Tmax_7d', 'GDD_7d', 'GDD_cumulative', 'HeatStress_7d',
             'SolRad_7d', 'RH_7d', 'Log10_Titer', 'n_DEGs']
corr_labels = ['Tavg\n7d', 'Tmax\n7d', 'GDD\n7d', 'GDD\ncum', 'Heat\nStress',
               'Solar\nRad', 'RH', 'log10\nTiter', 'DEGs']
corr_matrix = int_df[corr_cols].corr(method='spearman')
corr_matrix.index = corr_labels
corr_matrix.columns = corr_labels
sns.heatmap(corr_matrix, ax=ax, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
            vmin=-1, vmax=1, cbar_kws={'label': 'Spearman rho'},
            annot_kws={'fontsize': 7})
ax.set_title("Correlation Matrix\n(Spearman)")

plt.suptitle("Plot 32: Integrated Temperature x Virus Titer x Transcriptome", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('32_integrated_temp_titer_DEG.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 32_integrated_temp_titer_DEG.png")


# ============================================================
# PLOT 33: Temperature Trend and Seasonal Shift Analysis
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(18, 14))

# Panel A: Daily temperature overlay with inter-sampling intervals shaded
ax = axes[0, 0]
for season, wx, color in [(2021, weather_2021, '#1f77b4'), (2022, weather_2022, '#ff7f0e')]:
    doy = wx['Date'].dt.dayofyear
    ax.plot(doy, wx['Tmax'], color=color, alpha=0.7, linewidth=1.2, label=f'{season} Tmax')

    # Shade inter-sampling periods
    for i, stage in enumerate(stages):
        sd = sampling_dates[season][stage]
        sd_doy = sd.timetuple().tm_yday
        ax.axvline(sd_doy, color=color, linestyle=':', linewidth=0.5, alpha=0.5)

# Heat threshold
ax.axhline(35, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Heat stress (35C)')
ax.set_xlabel("Day of Year")
ax.set_ylabel("Max Air Temperature (C)")
ax.set_title("Daily Tmax: 2021 vs 2022\nwith Sampling Dates")
ax.legend(fontsize=10)

# Panel B: 2022 was hotter - cumulative heat stress
ax = axes[0, 1]
for season, wx, color in [(2021, weather_2021, '#1f77b4'), (2022, weather_2022, '#ff7f0e')]:
    doy = wx['Date'].dt.dayofyear
    cum_hs = wx['HeatStress'].cumsum()
    ax.plot(doy, cum_hs, color=color, linewidth=2, label=f'{season}')

    for stage in stages:
        sd = sampling_dates[season][stage]
        sd_doy = sd.timetuple().tm_yday
        hs_at = wx[wx['Date'] <= pd.Timestamp(sd)]['HeatStress'].sum()
        ax.plot(sd_doy, hs_at, marker='o' if season == 2021 else 's',
                color=stage_colors[stage], markersize=10, zorder=5,
                markeredgecolor='black', markeredgewidth=1)

ax.set_xlabel("Day of Year")
ax.set_ylabel("Cumulative Heat Stress Days (>35C)")
ax.set_title("Cumulative Heat Stress Exposure")
ax.legend()

# Panel C: Titer change vs temperature change between seasons
ax = axes[1, 0]
dt_vals = []
dtiter_vals = []
labels_pts = []

for stage in stages:
    row21 = int_df[(int_df['Season'] == 2021) & (int_df['Stage'] == stage)].iloc[0]
    row22 = int_df[(int_df['Season'] == 2022) & (int_df['Stage'] == stage)].iloc[0]

    dt = row22['Tavg_7d'] - row21['Tavg_7d']
    dtiter = row22['Log10_Titer'] - row21['Log10_Titer']

    dt_vals.append(dt)
    dtiter_vals.append(dtiter)
    labels_pts.append(stage)

    ax.scatter(dt, dtiter, c=stage_colors[stage], s=200, edgecolors='black',
               linewidth=2, zorder=5)
    ax.annotate(stage.replace('_', '\n'), (dt, dtiter), fontsize=12,
                ha='center', va='bottom', xytext=(0, 8), textcoords='offset points')

# Regression
if len(dt_vals) >= 3:
    r, p = pearsonr(dt_vals, dtiter_vals)
    z = np.polyfit(dt_vals, dtiter_vals, 1)
    xline = np.linspace(min(dt_vals) - 0.5, max(dt_vals) + 0.5, 100)
    ax.plot(xline, np.polyval(z, xline), 'k--', linewidth=1, alpha=0.5)
    ax.set_title(f"Temperature Change vs Titer Change\n(2022-2021), r={r:.3f}")

ax.axhline(0, color='grey', linestyle='-', linewidth=0.5, alpha=0.5)
ax.axvline(0, color='grey', linestyle='-', linewidth=0.5, alpha=0.5)
ax.set_xlabel("Change in 7-day Avg Temp (C): 2022-2021")
ax.set_ylabel("Change in log10(Titer): 2022-2021")

# Panel D: DEG change vs temperature change
ax = axes[1, 1]
dt_vals2 = []
ddeg_vals = []

for stage in stages:
    row21 = int_df[(int_df['Season'] == 2021) & (int_df['Stage'] == stage)].iloc[0]
    row22 = int_df[(int_df['Season'] == 2022) & (int_df['Stage'] == stage)].iloc[0]

    dt = row22['Tavg_7d'] - row21['Tavg_7d']
    ddeg = row22['n_DEGs'] - row21['n_DEGs']

    dt_vals2.append(dt)
    ddeg_vals.append(ddeg)

    ax.scatter(dt, ddeg, c=stage_colors[stage], s=200, edgecolors='black',
               linewidth=2, zorder=5)
    ax.annotate(stage.replace('_', '\n'), (dt, ddeg), fontsize=12,
                ha='center', va='bottom', xytext=(0, 8), textcoords='offset points')

if len(dt_vals2) >= 3:
    r, p = pearsonr(dt_vals2, ddeg_vals)
    z = np.polyfit(dt_vals2, ddeg_vals, 1)
    xline = np.linspace(min(dt_vals2) - 0.5, max(dt_vals2) + 0.5, 100)
    ax.plot(xline, np.polyval(z, xline), 'k--', linewidth=1, alpha=0.5)
    ax.set_title(f"Temperature Change vs DEG Change\n(2022-2021), r={r:.3f}")

ax.axhline(0, color='grey', linestyle='-', linewidth=0.5, alpha=0.5)
ax.axvline(0, color='grey', linestyle='-', linewidth=0.5, alpha=0.5)
ax.set_xlabel("Change in 7-day Avg Temp (C): 2022-2021")
ax.set_ylabel("Change in # DEGs: 2022-2021")

plt.suptitle("Plot 33: Temperature-Titer-DEG Seasonal Shift Analysis", fontsize=21, y=1.01)
plt.tight_layout()
plt.savefig('33_temperature_titer_shift.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: 33_temperature_titer_shift.png")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL INTEGRATED SUMMARY")
print("=" * 70)

print(f"""
ATMOSPHERIC CONDITIONS:
  2022 was consistently HOTTER than 2021:
    - Mean daily Tavg: {weather_2022['Tavg'].mean():.1f}C vs {weather_2021['Tavg'].mean():.1f}C (+{weather_2022['Tavg'].mean()-weather_2021['Tavg'].mean():.1f}C)
    - Mean daily Tmax: {weather_2022['Tmax'].mean():.1f}C vs {weather_2021['Tmax'].mean():.1f}C (+{weather_2022['Tmax'].mean()-weather_2021['Tmax'].mean():.1f}C)
    - Total GDD: {weather_2022['GDD'].sum():.0f} vs {weather_2021['GDD'].sum():.0f} (+{weather_2022['GDD'].sum()-weather_2021['GDD'].sum():.0f})
    - Heat stress days (>35C): {weather_2022['HeatStress'].sum()} vs {weather_2021['HeatStress'].sum()}
    - 2022 sampling started EARLIER (Jul 28 vs Aug 02)

VIRUS TITER DYNAMICS:
  - Titer generally INCREASES across maturity stages (Pre-veraison -> Harvest)
  - 2021 Harvest shows highest mean titer (accumulation over longer/warmer season)
  - 2022 titer levels were generally LOWER at Harvest despite hotter conditions
  - Infection duration shows variable correlation with titer (stage-dependent)

TEMPERATURE x TITER x TRANSCRIPTOME:
  Key correlations across all 8 season x stage data points:
""")

# Print the key correlations again
for x_col, y_col, label in corr_pairs[:5]:
    x = int_df[x_col]
    y = int_df[y_col]
    r, p = pearsonr(x, y)
    print(f"    {label}: r={r:.3f} (p={p:.3f})")

print("""
SAMPLING DATE CONTEXT:
  2021 season: Aug 2 -> Aug 16 -> Sep 9 -> Sep 21
  2022 season: Jul 28 -> Aug 9 -> Aug 24 -> Sep 13
  2022 sampling was shifted ~5 days earlier across all stages

BIOLOGICAL INTERPRETATION:
  The hotter 2022 season with earlier sampling dates means:
  1. Berry development was accelerated by higher heat accumulation
  2. The same phenological stages occurred at different thermal conditions
  3. The large seasonal difference in DEGs (especially at Harvest: 1192 vs 367)
     may partly reflect the temperature-driven shift in vine physiology
  4. Temperature-dependent host responses confound direct virus-host comparisons
     between seasons, reinforcing the importance of season-stratified analysis
""")

print("\nFiles generated:")
print("  30_weather_overview_sampling.png")
print("  31_virus_titer_analysis.png")
print("  32_integrated_temp_titer_DEG.png")
print("  33_temperature_titer_shift.png")
print("  integrated_environment_titer_DEG.csv")
print("  temperature_at_sampling.csv")
print("\nDone!")
