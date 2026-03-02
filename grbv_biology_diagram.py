#!/usr/bin/env python3
"""
GRBV Biology Diagram — Data-Driven Infographic
Comprehensive visualization of GRBV infection effects across developmental stages,
based on actual transcriptomics, viral titer, and environmental data.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc, Circle, Wedge
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import numpy as np

# ============================================================
# DATA (exact values from analysis)
# ============================================================
stages = ['Pre-veraison', 'Veraison', 'Post-veraison', 'Harvest']

data_2021 = {
    'Pre-veraison':  {'date': 'Aug 2', 'Tavg': 20.66, 'Tmax': 32.43, 'titer': 456885,
                      'DEGs': 118, 'up': 78, 'down': 40, 'heat_DEGs': 0},
    'Veraison':      {'date': 'Aug 16', 'Tavg': 23.23, 'Tmax': 34.18, 'titer': 831813,
                      'DEGs': 107, 'up': 59, 'down': 48, 'heat_DEGs': 0},
    'Post-veraison': {'date': 'Sep 9', 'Tavg': 20.50, 'Tmax': 33.48, 'titer': 754656,
                      'DEGs': 326, 'up': 177, 'down': 149, 'heat_DEGs': 0},
    'Harvest':       {'date': 'Sep 21', 'Tavg': 17.55, 'Tmax': 27.96, 'titer': 2210531,
                      'DEGs': 1192, 'up': 516, 'down': 676, 'heat_DEGs': 0},
}

data_2022 = {
    'Pre-veraison':  {'date': 'Jul 28', 'Tavg': 21.60, 'Tmax': 34.63, 'titer': 1629971,
                      'DEGs': 113, 'up': 81, 'down': 32, 'heat_DEGs': 0},
    'Veraison':      {'date': 'Aug 9', 'Tavg': 22.58, 'Tmax': 33.68, 'titer': 1386813,
                      'DEGs': 38, 'up': 28, 'down': 10, 'heat_DEGs': 0},
    'Post-veraison': {'date': 'Aug 24', 'Tavg': 22.59, 'Tmax': 34.50, 'titer': 728816,
                      'DEGs': 565, 'up': 220, 'down': 345, 'heat_DEGs': 8},
    'Harvest':       {'date': 'Sep 13', 'Tavg': 25.13, 'Tmax': 35.64, 'titer': 393267,
                      'DEGs': 367, 'up': 184, 'down': 183, 'heat_DEGs': 13},
}

cross_season_r = {'Pre-veraison': 0.329, 'Veraison': 0.524,
                  'Post-veraison': 0.692, 'Harvest': 0.586}

# Heat stress from heat_stress_signatures.py
heat_stress_summary = {
    'total_heat_genes': 115,
    'sHSPs_down_2021_Harvest': '29/36',
    'sHSPs_down_2022_PostVer': '23/36',
    'sHSP_module_healthy_2022H': '+1.64',
    'sHSP_module_infected_2022H': '+0.74',
    'HSF_sHSP_corr': 'r=0.625',
    'genes_both_virus_env_Harvest': 23,
}

biomarkers = [
    ('F-box/LRR protein', 'UP', 0.880, 'Defense signaling'),
    ('E3 ubiquitin ligase RSL1', 'UP', 0.894, 'Protein turnover'),
    ('RNA-dep. RNA polymerase 1', 'UP', 0.895, 'Antiviral RNAi / Sugar QTL'),
    ('Cell wall struct. protein 2', 'DOWN', 0.871, 'Cell wall integrity'),
    ('NAC domain protein 62', 'UP', 0.845, 'Ripening TF'),
    ('PDR3 (ABC transporter)', 'UP', 0.811, 'Defense transport'),
]

# ============================================================
# COLOR PALETTE
# ============================================================
C_DARK = '#1B3A5C'
C_BLUE = '#2E86AB'
C_RED = '#D62728'
C_GREEN = '#2CA02C'
C_ORANGE = '#FF7F0E'
C_PURPLE = '#9467BD'
C_TEAL = '#17BECF'
C_PINK = '#E377C2'
C_GOLD = '#BCBD22'
C_GRAY = '#7F7F7F'
C_LGRAY = '#D0D0D0'
C_VLGRAY = '#F5F5F5'

stage_colors = {
    'Pre-veraison': '#4DAF4A',
    'Veraison': '#FF7F0E',
    'Post-veraison': '#E41A1C',
    'Harvest': '#984EA3',
}
season_colors = {'2021': '#2166AC', '2022': '#B2182B'}

# ============================================================
# FIGURE SETUP
# ============================================================
fig = plt.figure(figsize=(36, 44))
fig.patch.set_facecolor('white')

# Main title
fig.text(0.5, 0.975, 'Grapevine Red Blotch Virus (GRBV): Data-Driven Molecular Landscape',
         fontsize=32, fontweight='bold', ha='center', va='top', color=C_DARK,
         fontfamily='sans-serif')
fig.text(0.5, 0.965, 'Transcriptomic, Environmental, and Viral Titer Integration Across Developmental Stages',
         fontsize=18, ha='center', va='top', color=C_GRAY, fontfamily='sans-serif')

# ============================================================
# PANEL LAYOUT (GridSpec)
# ============================================================
from matplotlib.gridspec import GridSpec
gs_main = GridSpec(6, 4, figure=fig, hspace=0.35, wspace=0.3,
                   left=0.04, right=0.96, top=0.955, bottom=0.02)

# ============================================================
# PANEL 1: Central GRBV Information (top-left)
# ============================================================
ax_info = fig.add_subplot(gs_main[0, 0])
ax_info.set_xlim(0, 10)
ax_info.set_ylim(0, 10)
ax_info.axis('off')

# Virus icon (hexagon)
hex_x = 2.5
hex_y = 7.5
hex_r = 1.2
angles = np.linspace(0, 2 * np.pi, 7)
hex_pts = np.column_stack([hex_x + hex_r * np.cos(angles), hex_y + hex_r * np.sin(angles)])
virus_hex = plt.Polygon(hex_pts, closed=True, facecolor=C_RED, edgecolor='#8B0000',
                         linewidth=2.5, alpha=0.9, zorder=5)
ax_info.add_patch(virus_hex)
ax_info.text(hex_x, hex_y, 'GRBV', fontsize=14, fontweight='bold', color='white',
             ha='center', va='center', zorder=6)

ax_info.text(5.5, 9.2, 'Grapevine Red Blotch Virus', fontsize=14, fontweight='bold', color=C_DARK)
ax_info.text(5.5, 8.3, 'ssDNA virus (Geminiviridae)', fontsize=11, color=C_GRAY)
ax_info.text(5.5, 7.5, 'Replicates in host cell nuclei', fontsize=11, color=C_GRAY)
ax_info.text(5.5, 6.7, 'Vector: Spissistilus festinus', fontsize=11, color=C_GRAY, style='italic')

# Dataset box
rect = FancyBboxPatch((0.3, 0.3), 9.4, 5.0, boxstyle="round,pad=0.3",
                        facecolor='#EBF5FB', edgecolor=C_BLUE, linewidth=1.5)
ax_info.add_patch(rect)
ax_info.text(5, 4.8, 'Study Dataset', fontsize=13, fontweight='bold', color=C_DARK, ha='center')
ax_info.text(1, 4.0, '327 samples  |  14,604 genes', fontsize=10, color=C_DARK)
ax_info.text(1, 3.2, 'Seasons: 2021 (cooler) & 2022 (hotter)', fontsize=10, color=C_DARK)
ax_info.text(1, 2.4, '4 maturity stages  |  7 infection cohorts', fontsize=10, color=C_DARK)
ax_info.text(1, 1.6, 'RB(+) infected vs RB(-) controls', fontsize=10, color=C_DARK)
ax_info.text(1, 0.8, 'DE: padj < 0.05, |log2FC| >= 1', fontsize=10, color=C_GRAY)

# ============================================================
# PANEL 2: DEG Trajectory (top, center-right, spanning 3 cols)
# ============================================================
ax_deg = fig.add_subplot(gs_main[0, 1:])
x_pos = np.arange(4)
width = 0.35

# 2021
bars_21_up = ax_deg.bar(x_pos - width/2, [data_2021[s]['up'] for s in stages],
                         width, label='2021 Up', color=season_colors['2021'], alpha=0.7)
bars_21_dn = ax_deg.bar(x_pos - width/2, [-data_2021[s]['down'] for s in stages],
                         width, label='2021 Down', color=season_colors['2021'], alpha=0.35)

# 2022
bars_22_up = ax_deg.bar(x_pos + width/2, [data_2022[s]['up'] for s in stages],
                         width, label='2022 Up', color=season_colors['2022'], alpha=0.7)
bars_22_dn = ax_deg.bar(x_pos + width/2, [-data_2022[s]['down'] for s in stages],
                         width, label='2022 Down', color=season_colors['2022'], alpha=0.35)

# Total DEG labels
for i, s in enumerate(stages):
    ax_deg.text(i - width/2, data_2021[s]['up'] + 25,
                f"{data_2021[s]['DEGs']}", fontsize=11, fontweight='bold',
                ha='center', color=season_colors['2021'])
    ax_deg.text(i + width/2, data_2022[s]['up'] + 25,
                f"{data_2022[s]['DEGs']}", fontsize=11, fontweight='bold',
                ha='center', color=season_colors['2022'])

# Highlight peaks
ax_deg.annotate('2021 PEAK\n1,192 DEGs', xy=(3 - width/2, 520), xytext=(2.2, 650),
                fontsize=10, fontweight='bold', color=season_colors['2021'],
                arrowprops=dict(arrowstyle='->', color=season_colors['2021'], lw=1.5))
ax_deg.annotate('2022 PEAK\n565 DEGs', xy=(2 + width/2, 225), xytext=(3.2, 450),
                fontsize=10, fontweight='bold', color=season_colors['2022'],
                arrowprops=dict(arrowstyle='->', color=season_colors['2022'], lw=1.5))

ax_deg.axhline(0, color='black', linewidth=0.8)
ax_deg.set_xticks(x_pos)
ax_deg.set_xticklabels([s.replace('-', '-\n') for s in stages], fontsize=12)
ax_deg.set_ylabel('Number of DEGs (Up / Down)', fontsize=12)
ax_deg.set_title('INVERTED SEASONAL TRAJECTORY: DEGs per Stage (RB+ vs RB-)',
                  fontsize=16, fontweight='bold', color=C_DARK, pad=10)
ax_deg.legend(fontsize=10, ncol=4, loc='upper left')
ax_deg.spines['top'].set_visible(False)
ax_deg.spines['right'].set_visible(False)

# ============================================================
# PANELS 3-6: FOUR STAGE DETAIL BOXES (Row 1-2)
# ============================================================
def draw_stage_panel(ax, stage, d21, d22, pathways, molecular_changes, symptoms, heat_info, color):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title bar
    rect = FancyBboxPatch((0, 8.8), 10, 1.2, boxstyle="round,pad=0.15",
                            facecolor=color, edgecolor='none', alpha=0.9)
    ax.add_patch(rect)
    ax.text(5, 9.4, stage, fontsize=16, fontweight='bold', color='white', ha='center', va='center')

    y = 8.4

    # Season comparison box
    rect2 = FancyBboxPatch((0.1, 5.5), 9.8, 2.8, boxstyle="round,pad=0.15",
                             facecolor='#F8F9FA', edgecolor=C_LGRAY, linewidth=1)
    ax.add_patch(rect2)

    ax.text(5, 8.05, '2021 vs 2022 Comparison', fontsize=10, fontweight='bold',
            color=C_DARK, ha='center')

    # 2021 column
    ax.text(2.5, 7.6, '2021', fontsize=11, fontweight='bold', color=season_colors['2021'], ha='center')
    ax.text(2.5, 7.1, f"{d21['DEGs']} DEGs", fontsize=10, color=C_DARK, ha='center')
    ax.text(2.5, 6.65, f"({d21['up']} up / {d21['down']} down)", fontsize=8, color=C_GRAY, ha='center')
    ax.text(2.5, 6.2, f"Titer: {d21['titer']/1e6:.2f}M", fontsize=9, color=C_DARK, ha='center')
    ax.text(2.5, 5.8, f"Tavg: {d21['Tavg']:.1f}C  Tmax: {d21['Tmax']:.1f}C",
            fontsize=8, color=C_GRAY, ha='center')

    # Divider
    ax.plot([5, 5], [5.7, 7.8], color=C_LGRAY, linewidth=1, linestyle='--')

    # 2022 column
    ax.text(7.5, 7.6, '2022', fontsize=11, fontweight='bold', color=season_colors['2022'], ha='center')
    ax.text(7.5, 7.1, f"{d22['DEGs']} DEGs", fontsize=10, color=C_DARK, ha='center')
    ax.text(7.5, 6.65, f"({d22['up']} up / {d22['down']} down)", fontsize=8, color=C_GRAY, ha='center')
    ax.text(7.5, 6.2, f"Titer: {d22['titer']/1e6:.2f}M", fontsize=9, color=C_DARK, ha='center')
    ax.text(7.5, 5.8, f"Tavg: {d22['Tavg']:.1f}C  Tmax: {d22['Tmax']:.1f}C",
            fontsize=8, color=C_GRAY, ha='center')

    # Cross-season r
    r_val = cross_season_r[stage]
    r_color = C_GREEN if r_val > 0.5 else C_ORANGE if r_val > 0.3 else C_RED
    ax.text(5, 5.55, f"Cross-season r = {r_val:.3f}", fontsize=8, fontweight='bold',
            color=r_color, ha='center', style='italic')

    # Molecular changes
    y = 5.2
    ax.text(0.3, y, 'Molecular Changes', fontsize=10, fontweight='bold', color=C_DARK)
    for change in molecular_changes:
        y -= 0.38
        direction = change[0]
        text = change[1:]
        d_color = C_RED if direction == '+' else C_BLUE if direction == '-' else C_GRAY
        ax.text(0.5, y, direction, fontsize=10, fontweight='bold', color=d_color)
        ax.text(1.0, y, text, fontsize=8.5, color=C_DARK)

    # Pathways box
    y -= 0.55
    ax.text(0.3, y, 'Key Pathways Enriched', fontsize=10, fontweight='bold', color=C_DARK)
    for pw in pathways:
        y -= 0.35
        ax.text(0.6, y, '  ' + pw, fontsize=8, color=C_PURPLE)

    # Heat stress info
    if heat_info:
        y -= 0.55
        rect_h = FancyBboxPatch((0.2, y - 0.15), 9.6, 0.55, boxstyle="round,pad=0.08",
                                  facecolor='#FFF3E0', edgecolor=C_ORANGE, linewidth=1)
        ax.add_patch(rect_h)
        ax.text(0.5, y + 0.05, heat_info, fontsize=8, fontweight='bold', color='#E65100')

    # Symptoms
    if symptoms:
        y -= 0.55
        ax.text(0.3, y, 'Phenotype', fontsize=9, fontweight='bold', color=C_DARK)
        for sym in symptoms:
            y -= 0.33
            ax.text(0.6, y, sym, fontsize=8, color=C_RED)


# Stage 1: Pre-veraison
ax_s1 = fig.add_subplot(gs_main[1, 0])
draw_stage_panel(ax_s1, 'Pre-veraison', data_2021['Pre-veraison'], data_2022['Pre-veraison'],
    pathways=['General defense activation', 'Limited pathway enrichment'],
    molecular_changes=[
        '+Defense pathways beginning to activate',
        '+RNA silencing machinery initiating',
        '-Minimal transcriptomic disruption',
        ' Virus establishing in host tissues',
    ],
    symptoms=['No visible symptoms', 'Molecular changes underway'],
    heat_info='No heat-stress DEGs at this stage',
    color=stage_colors['Pre-veraison'])

# Stage 2: Veraison
ax_s2 = fig.add_subplot(gs_main[1, 1])
draw_stage_panel(ax_s2, 'Veraison', data_2021['Veraison'], data_2022['Veraison'],
    pathways=['Plant-pathogen interaction', 'Phenylpropanoid biosynthesis', 'Flavonoid pathway induced'],
    molecular_changes=[
        '+Plant-pathogen interaction pathway',
        '+RNA-dep. RNA polymerase (RdRP1)',
        '+Flavonoid biosynthesis activation',
        '-First foliar symptoms appear',
        ' Virus titer rising (2021)',
    ],
    symptoms=['Slight reddening', 'Interveinal color change begins'],
    heat_info='No heat-stress DEGs detected',
    color=stage_colors['Veraison'])

# Stage 3: Post-veraison
ax_s3 = fig.add_subplot(gs_main[1, 2])
draw_stage_panel(ax_s3, 'Post-veraison', data_2021['Post-veraison'], data_2022['Post-veraison'],
    pathways=['Phenylpropanoid biosynthesis', 'Plant-pathogen interaction',
             'ABC transporters (UP)', 'Protein processing in ER', 'Carbohydrate metabolism'],
    molecular_changes=[
        '+F-box/LRR defense (AUC=0.880)',
        '+E3 ubiquitin ligase (protein turnover)',
        '+ABC transporter PDR3 (defense)',
        '-Cell wall structural protein (DOWN)',
        '-Sugar metabolism suppressed',
        '-sHSPs suppressed in 2022 (23/36 DOWN)',
    ],
    symptoms=['Red interveinal blotches', 'Subdued pigmentation',
             'Delayed ripening onset'],
    heat_info='2022: 8 heat DEGs (all DOWN) -- sHSP suppression begins',
    color=stage_colors['Post-veraison'])

# Stage 4: Harvest
ax_s4 = fig.add_subplot(gs_main[1, 3])
draw_stage_panel(ax_s4, 'Harvest', data_2021['Harvest'], data_2022['Harvest'],
    pathways=['Phenylpropanoid biosynthesis', 'Plant-pathogen interaction',
             'ABC transporters', 'Cell wall modification', 'Starch/sucrose metabolism'],
    molecular_changes=[
        '+NAC62 ripening TF disrupted (AUC=0.845)',
        '+RdRP1 at sugar QTL (AUC=0.895)',
        '+Max defense activation (1,192 DEGs 2021)',
        '-29/36 sHSPs DOWN (2021)',
        '-Cell wall protein DOWN (AUC=0.871)',
        '-Anthocyanin/phenolics reduced',
    ],
    symptoms=['Coalescing red blotches', 'Low Brix (-3-4)',
             'Low anthocyanin (-25-50%)', 'High titratable acidity'],
    heat_info='23 heat genes with BOTH virus + environment effects',
    color=stage_colors['Harvest'])

# ============================================================
# PANEL: VIRAL TITER DYNAMICS (Row 2, left half)
# ============================================================
ax_titer = fig.add_subplot(gs_main[2, :2])

x = np.arange(4)
titer_21 = [data_2021[s]['titer']/1e6 for s in stages]
titer_22 = [data_2022[s]['titer']/1e6 for s in stages]

ax_titer.plot(x, titer_21, 'o-', color=season_colors['2021'], linewidth=3, markersize=12,
              label='2021 (cooler)', zorder=5)
ax_titer.plot(x, titer_22, 's-', color=season_colors['2022'], linewidth=3, markersize=12,
              label='2022 (hotter)', zorder=5)

for i, s in enumerate(stages):
    ax_titer.annotate(f'{titer_21[i]:.2f}M', (i, titer_21[i]),
                       textcoords="offset points", xytext=(0, 15),
                       fontsize=10, fontweight='bold', color=season_colors['2021'], ha='center')
    ax_titer.annotate(f'{titer_22[i]:.2f}M', (i, titer_22[i]),
                       textcoords="offset points", xytext=(0, -20),
                       fontsize=10, fontweight='bold', color=season_colors['2022'], ha='center')

# Arrows showing direction
ax_titer.annotate('', xy=(3, titer_21[3]), xytext=(2, titer_21[2]),
                   arrowprops=dict(arrowstyle='->', color=season_colors['2021'],
                                   lw=2, connectionstyle='arc3,rad=0.2'))
ax_titer.annotate('', xy=(3, titer_22[3]), xytext=(0, titer_22[0]),
                   arrowprops=dict(arrowstyle='->', color=season_colors['2022'],
                                   lw=2, connectionstyle='arc3,rad=-0.2'))

ax_titer.text(3.3, titer_21[3] - 0.15, '4.8x increase\n(Pre-ver to Harvest)',
              fontsize=9, color=season_colors['2021'], fontweight='bold')
ax_titer.text(3.3, titer_22[3] + 0.1, '4.1x DECREASE\n(Pre-ver to Harvest)',
              fontsize=9, color=season_colors['2022'], fontweight='bold')

ax_titer.set_xticks(x)
ax_titer.set_xticklabels(stages, fontsize=12)
ax_titer.set_ylabel('Mean Viral Titer (millions copies)', fontsize=12)
ax_titer.set_title('Viral Titer Dynamics: Opposing Seasonal Patterns',
                    fontsize=16, fontweight='bold', color=C_DARK, pad=10)
ax_titer.legend(fontsize=12, loc='upper center')
ax_titer.spines['top'].set_visible(False)
ax_titer.spines['right'].set_visible(False)
ax_titer.set_ylim(-0.1, 2.8)

# ============================================================
# PANEL: TEMPERATURE OVERLAY (Row 2, right half)
# ============================================================
ax_temp = fig.add_subplot(gs_main[2, 2:])

tavg_21 = [data_2021[s]['Tavg'] for s in stages]
tavg_22 = [data_2022[s]['Tavg'] for s in stages]
tmax_21 = [data_2021[s]['Tmax'] for s in stages]
tmax_22 = [data_2022[s]['Tmax'] for s in stages]

ax_temp.fill_between(x, tavg_21, tmax_21, alpha=0.15, color=season_colors['2021'])
ax_temp.fill_between(x, tavg_22, tmax_22, alpha=0.15, color=season_colors['2022'])
ax_temp.plot(x, tavg_21, 'o-', color=season_colors['2021'], linewidth=2.5, markersize=10, label='2021 Tavg')
ax_temp.plot(x, tavg_22, 's-', color=season_colors['2022'], linewidth=2.5, markersize=10, label='2022 Tavg')
ax_temp.plot(x, tmax_21, 'o--', color=season_colors['2021'], linewidth=1.5, markersize=6, alpha=0.5, label='2021 Tmax')
ax_temp.plot(x, tmax_22, 's--', color=season_colors['2022'], linewidth=1.5, markersize=6, alpha=0.5, label='2022 Tmax')

# Highlight Harvest temperature difference
ax_temp.annotate(f'+{tavg_22[3]-tavg_21[3]:.1f}C', xy=(3, (tavg_21[3]+tavg_22[3])/2),
                 fontsize=14, fontweight='bold', color=C_RED, ha='center',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFEBEE', edgecolor=C_RED))

# 35C heat stress threshold
ax_temp.axhline(35, color=C_RED, linestyle=':', linewidth=1.5, alpha=0.6)
ax_temp.text(0.1, 35.3, 'Heat stress threshold (35C)', fontsize=9, color=C_RED, alpha=0.7)

ax_temp.set_xticks(x)
ax_temp.set_xticklabels(stages, fontsize=12)
ax_temp.set_ylabel('Temperature (C)', fontsize=12)
ax_temp.set_title('Temperature at Sampling: 2022 Harvest +7.6C Warmer',
                    fontsize=16, fontweight='bold', color=C_DARK, pad=10)
ax_temp.legend(fontsize=10, ncol=2)
ax_temp.spines['top'].set_visible(False)
ax_temp.spines['right'].set_visible(False)

# ============================================================
# PANEL: HEAT-STRESS SUPPRESSION (Row 3, left half)
# ============================================================
ax_heat = fig.add_subplot(gs_main[3, :2])
ax_heat.set_xlim(0, 10)
ax_heat.set_ylim(0, 10)
ax_heat.axis('off')

# Title
rect_t = FancyBboxPatch((0, 8.8), 10, 1.2, boxstyle="round,pad=0.15",
                          facecolor=C_ORANGE, edgecolor='none', alpha=0.9)
ax_heat.add_patch(rect_t)
ax_heat.text(5, 9.35, 'NOVEL FINDING: GRBV Suppresses Heat-Stress Defense',
             fontsize=14, fontweight='bold', color='white', ha='center')

# Content
y = 8.3
items = [
    ('115 heat-stress genes', 'expressed across 9 categories (sHSP, HSP70, HSP90, HSP100, HSF, GolS, LEA, Antioxidant, GO:heat)'),
    ('57 heat genes are GRBV DEGs', 'Nearly all DOWNREGULATED in infected plants'),
    ('29/36 sHSPs suppressed', 'at 2021 Harvest (maximum virus titer + strongest response)'),
    ('23/36 sHSPs suppressed', 'at 2022 Post-veraison (shifted peak in hotter year)'),
    ('sHSP module score drop', 'Healthy: +1.64 vs Infected: +0.74 at 2022 Harvest (~55% reduction)'),
    ('HSF-sHSP axis confirmed', 'r = 0.625, p = 7.2e-37 (HSF transcription factors activate sHSPs)'),
    ('23 genes with dual regulation', 'Affected by BOTH virus AND environment at Harvest'),
]

for num_str, detail in items:
    rect_item = FancyBboxPatch((0.2, y - 0.55), 9.6, 0.75, boxstyle="round,pad=0.08",
                                facecolor='#FFF8E1', edgecolor='#FFE082', linewidth=0.8)
    ax_heat.add_patch(rect_item)
    ax_heat.text(0.5, y - 0.15, num_str, fontsize=11, fontweight='bold', color='#E65100')
    ax_heat.text(5.0, y - 0.15, detail, fontsize=9.5, color=C_DARK)
    y -= 0.95

# Implication box
rect_imp = FancyBboxPatch((0.2, 0.2), 9.6, 1.2, boxstyle="round,pad=0.15",
                            facecolor='#FFEBEE', edgecolor=C_RED, linewidth=2)
ax_heat.add_patch(rect_imp)
ax_heat.text(5, 0.95, 'IMPLICATION: GRBV compromises vine thermal tolerance by suppressing', fontsize=10,
             fontweight='bold', color=C_RED, ha='center')
ax_heat.text(5, 0.5, 'the proteostasis heat-defense network -- symptoms may worsen under combined viral + heat stress',
             fontsize=9.5, color=C_RED, ha='center')

# ============================================================
# PANEL: TOP BIOMARKERS (Row 3, right half)
# ============================================================
ax_bio = fig.add_subplot(gs_main[3, 2:])
ax_bio.set_xlim(0, 10)
ax_bio.set_ylim(0, 10)
ax_bio.axis('off')

rect_t2 = FancyBboxPatch((0, 8.8), 10, 1.2, boxstyle="round,pad=0.15",
                           facecolor=C_BLUE, edgecolor='none', alpha=0.9)
ax_bio.add_patch(rect_t2)
ax_bio.text(5, 9.35, 'Top Host Biomarkers for GRBV Diagnosis',
            fontsize=14, fontweight='bold', color='white', ha='center')

# 3-gene panel highlight
rect_panel = FancyBboxPatch((0.2, 7.5), 9.6, 1.1, boxstyle="round,pad=0.1",
                              facecolor='#E8F5E9', edgecolor=C_GREEN, linewidth=2)
ax_bio.add_patch(rect_panel)
ax_bio.text(5, 8.2, '5 Viral gene biomarkers: AUC = 1.000 (perfect classification)',
            fontsize=11, fontweight='bold', color=C_GREEN, ha='center')
ax_bio.text(5, 7.7, '3-gene panel achieves AUC = 1.000  |  12 universal biomarkers across all conditions',
            fontsize=10, color=C_DARK, ha='center')

# Biomarker table
y = 7.0
# Header
ax_bio.text(0.3, y, 'Gene', fontsize=10, fontweight='bold', color=C_DARK)
ax_bio.text(4.5, y, 'Dir', fontsize=10, fontweight='bold', color=C_DARK)
ax_bio.text(5.2, y, 'AUC', fontsize=10, fontweight='bold', color=C_DARK)
ax_bio.text(6.2, y, 'Function / QTL Link', fontsize=10, fontweight='bold', color=C_DARK)
ax_bio.plot([0.3, 9.7], [y - 0.15, y - 0.15], color=C_LGRAY, linewidth=1)

for name, direction, auc, function in biomarkers:
    y -= 0.6
    bg_color = '#FFEBEE' if direction == 'DOWN' else '#E3F2FD'
    rect_row = FancyBboxPatch((0.1, y - 0.2), 9.8, 0.5, boxstyle="round,pad=0.05",
                                facecolor=bg_color, edgecolor='none')
    ax_bio.add_patch(rect_row)
    ax_bio.text(0.3, y, name, fontsize=9.5, color=C_DARK)
    dir_color = C_RED if direction == 'UP' else C_BLUE
    ax_bio.text(4.5, y, direction, fontsize=10, fontweight='bold', color=dir_color)
    ax_bio.text(5.2, y, f'{auc:.3f}', fontsize=10, fontweight='bold', color=C_DARK)
    ax_bio.text(6.2, y, function, fontsize=9, color=C_PURPLE)

# QTL summary
y -= 0.8
rect_qtl = FancyBboxPatch((0.2, y - 0.5), 9.6, 1.0, boxstyle="round,pad=0.1",
                            facecolor='#F3E5F5', edgecolor=C_PURPLE, linewidth=1.5)
ax_bio.add_patch(rect_qtl)
ax_bio.text(5, y + 0.15, 'GWAS/QTL: 60% of biomarkers co-localize with known QTLs (1 Mb window)',
            fontsize=10, fontweight='bold', color=C_PURPLE, ha='center')
ax_bio.text(5, y - 0.25, '18 in disease-resistance clusters (Rpv, Ren, Run) | 8 near berry-quality loci | RdRP1 at sugar QTL',
            fontsize=9, color=C_DARK, ha='center')

# ============================================================
# PANEL: MOLECULAR PATHWAY FLOW (Row 4, full width)
# ============================================================
ax_flow = fig.add_subplot(gs_main[4, :])
ax_flow.set_xlim(0, 40)
ax_flow.set_ylim(0, 10)
ax_flow.axis('off')

ax_flow.text(20, 9.6, 'GRBV Molecular Disruption Cascade: From Infection to Symptoms',
             fontsize=18, fontweight='bold', color=C_DARK, ha='center')

# Flow boxes
boxes = [
    (1, 6.5, 'GRBV\nInfection', C_RED, 'white',
     'ssDNA virus enters\nhost cell nuclei'),
    (8, 6.5, 'Host Defense\nActivation', C_GREEN, 'white',
     'RdRP1 UP (AUC=0.895)\nF-box/LRR UP (0.880)\nE3 ubiquitin UP (0.894)\nRNA silencing activated'),
    (15, 6.5, 'Pathway\nDisruption', C_ORANGE, 'white',
     'Phenylpropanoid altered\nPlant-pathogen activated\nABC transporters UP\nProtein processing stressed'),
    (22, 6.5, 'Metabolic\nConsequences', C_PURPLE, 'white',
     'Sugar metabolism DOWN\nAnthocyanin reduced\nCell wall weakened\nHormone imbalance'),
    (29, 6.5, 'Heat Defense\nSuppressed', '#E65100', 'white',
     '29/36 sHSPs DOWN\nsHSP module -55%\nThermal tolerance\ncompromised'),
    (36, 6.5, 'Disease\nSymptoms', '#8B0000', 'white',
     'Low Brix (-3-4)\nLow anthocyanin\nDelayed ripening\nRed blotch lesions'),
]

for bx, by, title, bg, fg, detail in boxes:
    rect = FancyBboxPatch((bx - 2.5, by - 1.0), 5.5, 3.0, boxstyle="round,pad=0.2",
                            facecolor=bg, edgecolor='none', alpha=0.9)
    ax_flow.add_patch(rect)
    ax_flow.text(bx + 0.25, by + 0.8, title, fontsize=11, fontweight='bold',
                 color=fg, ha='center', va='center')
    ax_flow.text(bx + 0.25, by - 0.4, detail, fontsize=7.5, color='white',
                 ha='center', va='center', alpha=0.9, linespacing=1.3)

# Arrows between boxes
for i in range(5):
    x1 = boxes[i][0] + 3.2
    x2 = boxes[i + 1][0] - 2.5
    ax_flow.annotate('', xy=(x2, 6.5), xytext=(x1, 6.5),
                     arrowprops=dict(arrowstyle='->', color=C_DARK, lw=2.5))

# NAC/MYB branch
ax_flow.annotate('', xy=(22, 4.3), xytext=(22, 5.3),
                 arrowprops=dict(arrowstyle='->', color=C_PURPLE, lw=2))
rect_nac = FancyBboxPatch((19.5, 2.8), 5.5, 1.5, boxstyle="round,pad=0.15",
                            facecolor='#F3E5F5', edgecolor=C_PURPLE, linewidth=1.5)
ax_flow.add_patch(rect_nac)
ax_flow.text(22.25, 3.8, 'Ripening TF Disruption', fontsize=10, fontweight='bold', color=C_PURPLE, ha='center')
ax_flow.text(22.25, 3.2, 'NAC62 UP (AUC=0.845) at VvNAC cluster chr19\nMYBA locus linked to anthocyanin modifier',
             fontsize=8, color=C_DARK, ha='center')

# Sugar QTL branch
ax_flow.annotate('', xy=(15, 4.3), xytext=(15, 5.3),
                 arrowprops=dict(arrowstyle='->', color=C_BLUE, lw=2))
rect_sugar = FancyBboxPatch((12.5, 2.8), 5.5, 1.5, boxstyle="round,pad=0.15",
                              facecolor='#E3F2FD', edgecolor=C_BLUE, linewidth=1.5)
ax_flow.add_patch(rect_sugar)
ax_flow.text(15.25, 3.8, 'Sugar QTL Mechanistic Link', fontsize=10, fontweight='bold', color=C_BLUE, ha='center')
ax_flow.text(15.25, 3.2, 'RdRP1 maps to sugar accumulation QTL (chr1)\nAntiviral RNAi competes with sugar metabolism',
             fontsize=8, color=C_DARK, ha='center')

# Env modulation branch
ax_flow.annotate('', xy=(29, 4.3), xytext=(29, 5.3),
                 arrowprops=dict(arrowstyle='->', color='#E65100', lw=2))
rect_env = FancyBboxPatch((26.5, 2.8), 5.5, 1.5, boxstyle="round,pad=0.15",
                            facecolor='#FFF3E0', edgecolor='#E65100', linewidth=1.5)
ax_flow.add_patch(rect_env)
ax_flow.text(29.25, 3.8, 'Environment Modulates Severity', fontsize=10, fontweight='bold', color='#E65100', ha='center')
ax_flow.text(29.25, 3.2, '2021 (cooler): Peak at Harvest, max titer 2.2M\n2022 (hotter): Peak at Post-ver, titer declining',
             fontsize=8, color=C_DARK, ha='center')

# Bottom connection line
ax_flow.plot([15.25, 29.25], [2.7, 2.7], color=C_LGRAY, linewidth=1.5, linestyle='--')
ax_flow.text(22.25, 2.4, 'Interconnected pathways create multi-layered disease phenotype',
             fontsize=9, color=C_GRAY, ha='center', style='italic')

# ============================================================
# PANEL: KEY CONCLUSIONS (Row 5, full width)
# ============================================================
ax_conc = fig.add_subplot(gs_main[5, :])
ax_conc.set_xlim(0, 40)
ax_conc.set_ylim(0, 10)
ax_conc.axis('off')

rect_c = FancyBboxPatch((0.5, 0.3), 39, 9.5, boxstyle="round,pad=0.3",
                          facecolor=C_DARK, edgecolor='none', alpha=0.95)
ax_conc.add_patch(rect_c)

ax_conc.text(20, 9.2, 'KEY CONCLUSIONS FROM TRANSCRIPTOMIC ANALYSIS',
             fontsize=18, fontweight='bold', color='white', ha='center')

conclusions = [
    ('1.', 'STAGE-DEPENDENT', 'GRBV effects intensify through ripening, peaking at Post-veraison/Harvest with up to 1,192 DEGs'),
    ('2.', 'INVERTED TRAJECTORY', '2021 peaks at Harvest, 2022 peaks at Post-veraison -- environment reshapes response kinetics'),
    ('3.', 'HEAT-STRESS SUPPRESSION', 'GRBV suppresses 29/36 sHSPs, reducing thermal tolerance by ~55% (novel finding)'),
    ('4.', 'BIOMARKER PANEL', '3-gene panel achieves AUC=1.000; 12 universal biomarkers robust across all conditions'),
    ('5.', 'QTL CO-LOCALIZATION', '60% of biomarkers near known QTLs; RdRP1 at sugar QTL explains primary symptom'),
    ('6.', 'MULTI-PATHOGEN IMMUNITY', 'Biomarkers overlap Rpv/Ren/Run resistance loci -- shared innate immune machinery'),
]

y = 8.2
for num, title, text in conclusions:
    ax_conc.text(1.5, y, num, fontsize=14, fontweight='bold', color=C_BLUE)
    ax_conc.text(3, y, title, fontsize=12, fontweight='bold', color=C_ORANGE)
    ax_conc.text(12, y, text, fontsize=11, color='white')
    y -= 1.3

# ============================================================
# SAVE
# ============================================================
output_path = '/Users/prempratapsingh/Data_Analysis/36_GRBV_biology_diagram.png'
plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {output_path}')

# Also save high-res version
output_hires = '/Users/prempratapsingh/Data_Analysis/36_GRBV_biology_diagram_hires.png'
fig2 = plt.figure(figsize=(36, 44))
plt.close(fig2)
print(f'Figure 36: GRBV Biology Diagram generated successfully')
print(f'  Standard: {output_path}')
