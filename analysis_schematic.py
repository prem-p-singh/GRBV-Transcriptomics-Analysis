#!/usr/bin/env python3
"""
Creates a clean, publication-quality schematic diagram of the GRBV analysis pipeline.
Simplified layout with large readable text. No emojis.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(22, 32))
ax.set_xlim(0, 100)
ax.set_ylim(0, 120)
ax.axis('off')
fig.patch.set_facecolor('white')

# ============================================================
# COLORS
# ============================================================
c_input = '#dbeafe'
c_input_b = '#2563eb'
c_prep = '#fef3c7'
c_prep_b = '#d97706'
c_de = '#fce7f3'
c_de_b = '#db2777'
c_pw = '#d1fae5'
c_pw_b = '#059669'
c_seas = '#ede9fe'
c_seas_b = '#7c3aed'
c_bio = '#ffedd5'
c_bio_b = '#ea580c'
c_out = '#f0fdf4'
c_out_b = '#16a34a'
c_arrow = '#374151'


def box(x, y, w, h, title, lines, fc, ec, title_size=12, line_size=9.5):
    """Draw a rounded box with bold title + body lines."""
    patch = FancyBboxPatch((x - w/2, y - h/2), w, h,
                            boxstyle="round,pad=0.35", facecolor=fc, edgecolor=ec,
                            linewidth=2.2, zorder=2)
    ax.add_patch(patch)
    # Title
    n_lines = 1 + len(lines)
    spacing = min(h / (n_lines + 1.2), 1.4)
    top_y = y + h/2 - spacing * 0.9
    ax.text(x, top_y, title, ha='center', va='center',
            fontsize=title_size, fontweight='bold', color=ec, zorder=3)
    # Body lines
    for i, line in enumerate(lines):
        ax.text(x, top_y - spacing * (i + 1), line, ha='center', va='center',
                fontsize=line_size, color='#1f2937', zorder=3)


def arrow(x1, y1, x2, y2, color=c_arrow, lw=2.5):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='->', color=color,
                         linewidth=lw, mutation_scale=18, zorder=1)
    ax.add_patch(a)


def section_label(x, y, text, color, fs=15):
    ax.text(x, y, text, ha='center', va='center', fontsize=fs,
            fontweight='bold', color=color)


def section_bg(x, y, w, h, color, alpha=0.12):
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6",
                            facecolor=color, edgecolor='none', alpha=alpha, zorder=0)
    ax.add_patch(patch)


# ============================================================
# TITLE
# ============================================================
ax.text(50, 119, 'GRBV Transcriptomics Analysis Pipeline', ha='center',
        fontsize=22, fontweight='bold', color='#111827')
ax.text(50, 117.5, 'Grapevine Red Blotch Virus  |  Multi-factorial RNA-Seq Study', ha='center',
        fontsize=12, color='#6b7280', style='italic')

# ============================================================
# STEP 1: INPUT DATA
# ============================================================
section_bg(2, 109.5, 96, 7, c_input_b)
section_label(50, 115.5, 'STEP 1:  INPUT DATA', c_input_b)

box(18, 112, 24, 3.5, 'Raw Count Matrix',
    ['28,168 genes x 327 samples', '(RNA-Seq read counts)'],
    c_input, c_input_b)

box(50, 112, 24, 3.5, 'Sample Metadata',
    ['327 samples x 5 variables', 'Infection | Season | Stage | YSI'],
    c_input, c_input_b)

box(82, 112, 24, 3.5, 'OrgDb Annotations',
    ['org.Vvinifera.eg.sqlite', '59,449 GO + 13,801 KEGG terms'],
    c_input, c_input_b)

# ============================================================
# STEP 2: PREPROCESSING
# ============================================================
arrow(34, 110, 34, 107)
arrow(50, 110, 50, 107)

section_bg(8, 99, 60, 8.5, c_prep_b)
section_label(38, 106.5, 'STEP 2:  PREPROCESSING & NORMALIZATION', c_prep_b)

box(25, 103, 26, 3.2, 'Gene Filtering',
    ['Keep genes with >= 10 counts in >= 10% samples', '28,168  -->  14,604 genes retained'],
    c_prep, c_prep_b, line_size=9)

arrow(25, 101.2, 25, 100.2)

box(25, 99.8, 26, 1.2, '', [], c_prep, c_prep_b)
ax.text(25, 99.8, 'CPM normalization  +  log2(CPM + 1)', ha='center', va='center',
        fontsize=10, fontweight='bold', color=c_prep_b, zorder=3)

# Experimental design side box
box(78, 103, 30, 6.5, 'Experimental Design',
    ['RB(+): 264 infected | RB(-): 63 control',
     'Seasons: 2021, 2022',
     'Stages: Pre-veraison -> Veraison -> Post-veraison -> Harvest',
     'Infection cohorts: 2016, 2019, 2020, 2021, 2022'],
    '#f9fafb', '#9ca3af', title_size=11, line_size=9)

# ============================================================
# STEP 3: THREE PARALLEL ANALYSIS TRACKS
# ============================================================
arrow(18, 98.8, 18, 96.5)
arrow(38, 98.8, 50, 96.5)
arrow(50, 98.8, 82, 96.5)

# --- TRACK A: Differential Expression ---
section_bg(1, 72.5, 30, 24.5, c_de_b)
section_label(16, 96, 'STEP 3a:', c_de_b, fs=12)
section_label(16, 94.8, 'DIFFERENTIAL EXPRESSION', c_de_b, fs=11)

box(16, 91.5, 26, 3.2, "Welch's t-test",
    ['Per gene across 14,604 genes'],
    c_de, c_de_b)

arrow(16, 89.7, 16, 88.5)

box(16, 87, 26, 3, 'BH FDR Correction',
    ['padj < 0.05 significance threshold'],
    c_de, c_de_b)

arrow(16, 85.3, 16, 84)

box(16, 82.5, 26, 3, 'DEG Identification',
    ['|log2FC| >= 1  &  padj < 0.05'],
    c_de, c_de_b)

arrow(16, 80.8, 16, 79.5)

box(16, 77.5, 26, 4, '20 Comparisons',
    ['RB(+) vs RB(-) per stage',
     'Per season x stage (8 combos)',
     'Per infection age (6 cohorts)'],
    c_de, c_de_b, line_size=8.5)

arrow(16, 75.3, 16, 74)

box(16, 73, 26, 2, '', [], c_out, c_de_b)
ax.text(16, 73, 'PCA + Volcano + Heatmaps + Variance partition', ha='center',
        fontsize=8, fontweight='bold', color=c_de_b, zorder=3)


# --- TRACK B: Pathway Enrichment ---
section_bg(34, 72.5, 32, 24.5, c_pw_b)
section_label(50, 96, 'STEP 3b:', c_pw_b, fs=12)
section_label(50, 94.8, 'PATHWAY ENRICHMENT', c_pw_b, fs=11)

box(50, 91.5, 28, 3, 'Gene-to-Term Mapping',
    ['LOC IDs -> GO (BP/MF/CC) & KEGG (vvi)'],
    c_pw, c_pw_b)

arrow(50, 89.8, 50, 88.5)

box(50, 87, 28, 3, "Fisher's Exact Test",
    ['2x2 contingency per term'],
    c_pw, c_pw_b)

arrow(50, 85.3, 50, 84)

box(50, 82.5, 28, 3, 'BH FDR Correction',
    ['padj < 0.05 threshold'],
    c_pw, c_pw_b)

arrow(50, 80.8, 50, 79.5)

box(50, 77.5, 28, 4, 'Enrichment Per Condition',
    ['UP / DOWN DEGs separately',
     'Stage-wise, season-wise, age-wise',
     'Pathway landscape heatmaps'],
    c_pw, c_pw_b, line_size=8.5)

arrow(50, 75.3, 50, 74)

box(50, 73, 28, 2, '', [], c_out, c_pw_b)
ax.text(50, 73, 'Enrichment bars + Dotplots + Heatmaps', ha='center',
        fontsize=8, fontweight='bold', color=c_pw_b, zorder=3)


# --- TRACK C: Seasonal x Maturity ---
section_bg(69, 72.5, 30, 24.5, c_seas_b)
section_label(84, 96, 'STEP 3c:', c_seas_b, fs=12)
section_label(84, 94.8, 'SEASONAL x MATURITY', c_seas_b, fs=11)

box(84, 91.5, 26, 3, 'Season-Stratified DE',
    ['RB(+) vs RB(-) in 2021 & 2022'],
    c_seas, c_seas_b)

arrow(84, 89.8, 84, 88.5)

box(84, 87, 26, 3, 'Cross-Season Comparison',
    ['FC correlation per stage (Pearson r)'],
    c_seas, c_seas_b)

arrow(84, 85.3, 84, 84)

box(84, 82.5, 26, 3, 'DEG Overlap Analysis',
    ['Shared vs season-specific genes'],
    c_seas, c_seas_b)

arrow(84, 80.8, 84, 79.5)

box(84, 77.5, 26, 4, 'Stage Transitions',
    ['Change in infection effect',
     'Intensifying vs weakening genes',
     'Transition pathway enrichment'],
    c_seas, c_seas_b, line_size=8.5)

arrow(84, 75.3, 84, 74)

box(84, 73, 26, 2, '', [], c_out, c_seas_b)
ax.text(84, 73, 'Season x Stage DEG + pathway grids', ha='center',
        fontsize=8, fontweight='bold', color=c_seas_b, zorder=3)


# ============================================================
# STEP 4: BIOMARKER DISCOVERY
# ============================================================
# Convergence arrows
arrow(16, 71.8, 30, 67, color=c_bio_b, lw=3)
arrow(50, 71.8, 50, 67, color=c_bio_b, lw=3)
arrow(84, 71.8, 70, 67, color=c_bio_b, lw=3)

section_bg(4, 39, 92, 28.5, c_bio_b)
section_label(50, 66.5, 'STEP 4:  BIOMARKER DISCOVERY', c_bio_b)

# Scoring row
box(18, 62.5, 24, 4, 'Multi-Condition DE',
    ['Score each gene across',
     '20 DE comparisons'],
    c_bio, c_bio_b)

box(50, 62.5, 28, 4, 'Biomarker Score',
    ['Consistency (30%) + Effect size (25%)',
     'Significance (20%) + Direction (15%)',
     'Season robustness (10%)'],
    c_bio, c_bio_b, line_size=8.5)

box(82, 62.5, 22, 4, 'AUC-ROC',
    ['Individual gene',
     'discriminatory power'],
    c_bio, c_bio_b)

arrow(30, 62.5, 36, 62.5, lw=2, color=c_bio_b)
arrow(64, 62.5, 71, 62.5, lw=2, color=c_bio_b)
arrow(50, 60.3, 50, 58.8)

# Integration
box(50, 57, 70, 3.2, 'Integrated Score + Panel Classification',
    ['Final Score = 0.7 x Biomarker Score + 0.3 x AUC',
     '5-fold CV with Logistic Regression & Random Forest (3 to 50 gene panels)'],
    c_bio, c_bio_b, line_size=9)

arrow(50, 55.2, 50, 53.5)

# Categories
box(16, 50.5, 20, 5, 'UNIVERSAL',
    ['12 biomarkers', 'Season-robust', '>= 3 stages', 'AUC >= 0.7'],
    '#d1fae5', '#059669', title_size=11, line_size=9)

box(38, 50.5, 20, 5, 'EARLY DETECT',
    ['35 biomarkers', 'Sig at Pre-veraison', 'or Veraison', 'AUC >= 0.65'],
    '#dbeafe', '#2563eb', title_size=11, line_size=9)

box(62, 50.5, 20, 5, 'HARVEST',
    ['23 biomarkers', 'Strongest at', 'Harvest stage', 'AUC >= 0.7'],
    '#ede9fe', '#7c3aed', title_size=11, line_size=9)

box(84, 50.5, 20, 5, 'AGE MARKERS',
    ['33 biomarkers', '>= 3 age groups', 'Consistent dir.', 'AUC >= 0.7'],
    '#ffedd5', '#ea580c', title_size=11, line_size=9)

# Functional annotation
arrow(50, 47.5, 50, 45.5)

box(50, 43, 75, 4.5, 'Functional Annotation of Biomarkers',
    ['UP in RB(+): RNA silencing (RdRp1), Ubiquitin defense (F-box, E3 ligase), NAC TFs, ABC transporters',
     'DOWN in RB(+): ER chaperone suppression (HSP80), Protein processing in ER, Stress response (H2O2, salt)',
     'Core GRBV signature: virus hijacks host protein quality control and suppresses stress tolerance'],
    c_pw, c_pw_b, title_size=11, line_size=8.5)


# ============================================================
# STEP 5: KEY RESULTS
# ============================================================
arrow(50, 40.5, 50, 38)

section_bg(4, 14, 92, 24.5, '#1e293b')

section_label(50, 37.5, 'KEY RESULTS', '#f8fafc', fs=15)

# Three result columns
box(18, 33, 24, 6, 'Top Biomarkers',
    ['5 viral genes (N761_gp1-5)',
     '   AUC = 1.000 (perfect)',
     '7 host genes (LOC series)',
     '   AUC = 0.81 - 0.90',
     '   Season-robust, >= 3 stages'],
    '#ecfdf5', '#059669', title_size=11, line_size=8.5)

box(50, 33, 28, 6, 'Classifier Performance',
    ['3-gene panel: AUC = 1.000 (5-fold CV)',
     '10-gene panel: AUC = 1.000 per season',
     '10-gene panel: AUC = 1.000 per stage',
     'Top 50 biomarker FC r = 0.960',
     '   (2021 vs 2022 consistency)'],
    '#eff6ff', '#2563eb', title_size=11, line_size=8.5)

box(82, 33, 24, 6, 'GRBV Biology',
    ['Core viral signature:',
     '  - ER chaperone suppression',
     '  - Antiviral RNA silencing',
     '  - Cell wall remodeling',
     '  - Delayed ripening at Harvest'],
    '#fff7ed', '#ea580c', title_size=11, line_size=8.5)

# 6 questions box
box(50, 25, 80, 6, '6 Research Questions Addressed',
    ['Q1: Infection duration: DEGs accumulate with age; acute response at Year 0, stabilization after',
     'Q2: Maturity stages: Dominant variance driver (68.5%); Harvest has most DEGs',
     'Q3: Seasonal effects: FC correlation r = 0.77; consistent core + season-specific response',
     'Q4: Environmental effects: Season = 5.2% variance; modulates timing of GRBV disruption',
     'Q5: Uninfected environment: 15% overlap with GRBV DEGs; 60% concordant direction',
     'Q6: Transcriptome trajectory: 6 temporal clusters; PCA drift with infection age'],
    '#f1f5f9', '#475569', title_size=11, line_size=8.5)

# Outputs box
box(50, 17, 80, 4.5, 'Analysis Outputs',
    ['22 publication-ready plots + GRBV_biomarkers_top50.csv',
     'Plots 01-07: DE analysis | 08-14: Pathway enrichment | 15-19: Seasonal progression | 20-22: Biomarkers',
     'Scripts: analysis.py | pathway_enrichment.py | seasonal_maturity_progression.py | biomarker_discovery.py'],
    '#f8fafc', '#94a3b8', title_size=11, line_size=8.5)

# Method tags at bottom
methods_text = ('Methods: Python 3  |  Welch t-test + BH FDR  |  Fisher exact test  |  '
                'PCA (sklearn)  |  Logistic Regression + Random Forest  |  '
                'goatools (GO OBO 2024-06-17)  |  KEGG REST API  |  OrgDb SQLite')
ax.text(50, 14.8, methods_text, ha='center', va='center', fontsize=7.5,
        color='#9ca3af', style='italic')

plt.savefig('23_analysis_schematic.png', bbox_inches='tight', dpi=200, facecolor='white')
plt.close()
print("Saved: 23_analysis_schematic.png")
