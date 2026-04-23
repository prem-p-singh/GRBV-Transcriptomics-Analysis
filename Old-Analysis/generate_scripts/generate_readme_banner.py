#!/usr/bin/env python3
"""
Generate professional README banner and infographic images for the
GRBV Transcriptomics Analysis GitHub repository.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle
import matplotlib.patheffects as pe
import numpy as np

# ============================================================
# 1. HERO BANNER — Project Title Banner
# ============================================================

fig, ax = plt.subplots(figsize=(16, 5))
ax.set_xlim(0, 16)
ax.set_ylim(0, 5)
ax.axis('off')

# Gradient background using imshow
gradient = np.linspace(0, 1, 256).reshape(1, -1)
gradient = np.vstack([gradient] * 50)
ax.imshow(gradient, aspect='auto', cmap='RdYlGn_r',
          extent=[0, 16, 0, 5], alpha=0.15, zorder=0)

# Add subtle DNA helix pattern
theta = np.linspace(0, 8 * np.pi, 500)
x_dna = np.linspace(-1, 17, 500)
y_dna1 = 2.5 + 1.8 * np.sin(theta)
y_dna2 = 2.5 - 1.8 * np.sin(theta)
ax.plot(x_dna, y_dna1, color='#8B0000', alpha=0.06, linewidth=8, zorder=1)
ax.plot(x_dna, y_dna2, color='#006400', alpha=0.06, linewidth=8, zorder=1)

# Connecting rungs
for i in range(0, 500, 15):
    if 0 <= x_dna[i] <= 16:
        ax.plot([x_dna[i], x_dna[i]], [y_dna1[i], y_dna2[i]],
                color='#555555', alpha=0.04, linewidth=3, zorder=1)

# Grape cluster icon (circles)
grape_colors = ['#722F37', '#8B2252', '#6B2D5B', '#5C1A33', '#7B3F5E']
grape_positions = [
    (1.3, 3.2), (1.7, 3.2), (1.1, 2.8), (1.5, 2.8), (1.9, 2.8),
    (1.3, 2.4), (1.7, 2.4), (1.5, 2.0)
]
for idx, (gx, gy) in enumerate(grape_positions):
    circle = Circle((gx, gy), 0.18, facecolor=grape_colors[idx % len(grape_colors)],
                     edgecolor='#4A0E1E', linewidth=1, alpha=0.85, zorder=5)
    ax.add_patch(circle)

# Grape stem
ax.plot([1.5, 1.5, 1.8], [3.5, 3.8, 4.2], color='#2E8B57', linewidth=2.5, alpha=0.8, zorder=4)
# Leaf
leaf_x = [1.8, 2.3, 2.6, 2.3, 1.8]
leaf_y = [4.2, 4.5, 4.2, 3.9, 4.2]
ax.fill(leaf_x, leaf_y, color='#228B22', alpha=0.6, zorder=4)
ax.plot(leaf_x, leaf_y, color='#006400', linewidth=1.5, alpha=0.8, zorder=4)

# Virus particle icon (right side)
virus_cx, virus_cy = 14.5, 2.5
circle_virus = Circle((virus_cx, virus_cy), 0.7, facecolor='#B22222',
                       edgecolor='#8B0000', linewidth=2, alpha=0.3, zorder=5)
ax.add_patch(circle_virus)
# Spikes
for angle in np.linspace(0, 2 * np.pi, 12, endpoint=False):
    sx = virus_cx + 0.7 * np.cos(angle)
    sy = virus_cy + 0.7 * np.sin(angle)
    ex = virus_cx + 1.0 * np.cos(angle)
    ey = virus_cy + 1.0 * np.sin(angle)
    ax.plot([sx, ex], [sy, ey], color='#8B0000', linewidth=2, alpha=0.3, zorder=5)
    spike_circle = Circle((ex, ey), 0.08, facecolor='#B22222', alpha=0.3, zorder=5)
    ax.add_patch(spike_circle)

# Inner pentagon pattern for virus
for angle in np.linspace(0, 2 * np.pi, 5, endpoint=False):
    ix = virus_cx + 0.35 * np.cos(angle)
    iy = virus_cy + 0.35 * np.sin(angle)
    small_c = Circle((ix, iy), 0.1, facecolor='#FF6347', alpha=0.3, zorder=5)
    ax.add_patch(small_c)

# Main title
ax.text(8, 3.3, 'GRBV Transcriptomics Analysis', fontsize=32, fontweight='bold',
        ha='center', va='center', color='#1a1a2e',
        path_effects=[pe.withStroke(linewidth=4, foreground='white')], zorder=10)

# Subtitle
ax.text(8, 2.2, 'Grapevine Red Blotch Virus × Vitis vinifera cv. Merlot',
        fontsize=16, ha='center', va='center', color='#333333', style='italic',
        path_effects=[pe.withStroke(linewidth=3, foreground='white')], zorder=10)

# Stats bar at bottom
stats = ['327 Samples', '14,604 Genes', '55 Figures', '20 Scripts', '2 Seasons']
stat_colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0', '#F44336']
bar_y = 0.6
for i, (stat, sc) in enumerate(zip(stats, stat_colors)):
    x_pos = 2.0 + i * 3.0
    fancy = FancyBboxPatch((x_pos - 1.1, bar_y - 0.35), 2.2, 0.7,
                            boxstyle="round,pad=0.1", facecolor=sc, alpha=0.85,
                            edgecolor='white', linewidth=1.5, zorder=8)
    ax.add_patch(fancy)
    ax.text(x_pos, bar_y, stat, fontsize=11, fontweight='bold',
            ha='center', va='center', color='white', zorder=9)

# Decorative line
ax.axhline(y=1.2, xmin=0.05, xmax=0.95, color='#cccccc', linewidth=0.5, zorder=3)

plt.tight_layout(pad=0)
plt.savefig('assets/banner.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("Saved assets/banner.png")


# ============================================================
# 2. PIPELINE OVERVIEW — Visual Workflow
# ============================================================

fig, ax = plt.subplots(figsize=(16, 6))
ax.set_xlim(0, 16)
ax.set_ylim(0, 6)
ax.axis('off')

# Light background
ax.add_patch(FancyBboxPatch((0.1, 0.1), 15.8, 5.8,
             boxstyle="round,pad=0.2", facecolor='#fafafa',
             edgecolor='#e0e0e0', linewidth=1))

# Title
ax.text(8, 5.5, 'Analysis Pipeline Overview', fontsize=20, fontweight='bold',
        ha='center', va='center', color='#1a1a2e')

# Pipeline stages
stages = [
    ('RNA-Seq\nData', '#E3F2FD', '#1565C0', '327 samples\n28,168 genes'),
    ('Quality\nFiltering', '#E8F5E9', '#2E7D32', '≥10 counts\n≥10% samples'),
    ('Differential\nExpression', '#FFF3E0', '#E65100', "Welch's t-test\nBH FDR"),
    ('Pathway\nEnrichment', '#F3E5F5', '#6A1B9A', 'GO + KEGG\nGSEA'),
    ('Biomarker\nDiscovery', '#FFEBEE', '#C62828', '3-gene panel\nAUC=1.000'),
    ('Integration\n& Networks', '#E0F7FA', '#00695C', 'WGCNA\nANOVA'),
]

box_w = 2.0
box_h = 2.2
start_x = 0.5
y_center = 2.8

for i, (title, bg_color, border_color, detail) in enumerate(stages):
    x = start_x + i * 2.55

    # Box
    box = FancyBboxPatch((x, y_center - box_h/2), box_w, box_h,
                          boxstyle="round,pad=0.15", facecolor=bg_color,
                          edgecolor=border_color, linewidth=2.5)
    ax.add_patch(box)

    # Step number circle
    num_circle = Circle((x + box_w/2, y_center + box_h/2 - 0.2), 0.22,
                         facecolor=border_color, edgecolor='white', linewidth=2, zorder=6)
    ax.add_patch(num_circle)
    ax.text(x + box_w/2, y_center + box_h/2 - 0.2, str(i+1),
            fontsize=11, fontweight='bold', ha='center', va='center',
            color='white', zorder=7)

    # Title
    ax.text(x + box_w/2, y_center + 0.3, title, fontsize=11, fontweight='bold',
            ha='center', va='center', color=border_color)

    # Detail
    ax.text(x + box_w/2, y_center - 0.5, detail, fontsize=9,
            ha='center', va='center', color='#555555')

    # Arrow
    if i < len(stages) - 1:
        ax.annotate('', xy=(x + box_w + 0.55, y_center),
                    xytext=(x + box_w + 0.05, y_center),
                    arrowprops=dict(arrowstyle='->', color='#888888',
                                    lw=2.5, connectionstyle='arc3'))

# Output summary
output_items = [
    ('55 Figures', '#FF9800'),
    ('26 CSV Tables', '#4CAF50'),
    ('PDF Report', '#2196F3'),
    ('PPTX Deck', '#9C27B0'),
]

for j, (item, color) in enumerate(output_items):
    ox = 2.5 + j * 3.0
    fancy = FancyBboxPatch((ox - 0.9, 0.3), 1.8, 0.6,
                            boxstyle="round,pad=0.08", facecolor=color, alpha=0.85,
                            edgecolor='white', linewidth=1.5)
    ax.add_patch(fancy)
    ax.text(ox, 0.6, item, fontsize=10, fontweight='bold',
            ha='center', va='center', color='white')

ax.text(8, 1.15, 'Outputs', fontsize=12, fontweight='bold',
        ha='center', va='center', color='#666666')

plt.tight_layout(pad=0)
plt.savefig('assets/pipeline_overview.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("Saved assets/pipeline_overview.png")


# ============================================================
# 3. KEY FINDINGS — Visual Summary Card
# ============================================================

fig, ax = plt.subplots(figsize=(16, 7))
ax.set_xlim(0, 16)
ax.set_ylim(0, 7)
ax.axis('off')

# Background
ax.add_patch(FancyBboxPatch((0.1, 0.1), 15.8, 6.8,
             boxstyle="round,pad=0.2", facecolor='#fafafa',
             edgecolor='#e0e0e0', linewidth=1))

# Title
ax.text(8, 6.5, 'Key Findings at a Glance', fontsize=22, fontweight='bold',
        ha='center', va='center', color='#1a1a2e')

findings = [
    ('Inverted Seasonal\nTrajectory', '2021 peaks at Harvest\n2022 peaks at Post-vér',
     '#E53935', '🔄'),
    ('sHSP Suppression', '29/36 sHSPs down\n~55% module reduction',
     '#FF6F00', '🌡'),
    ('3-Gene Diagnostic\nPanel', 'AUC = 1.000\nacross 5-fold CV',
     '#2E7D32', '🎯'),
    ('908 Universal\nDEGs', 'Robust across pooled\n& stratified approaches',
     '#1565C0', '🧬'),
    ('270 GSEA\nEnrichments', '77% down-regulated\nRibosome top pathway',
     '#6A1B9A', '📊'),
    ('13 Co-expression\nModules', 'Stage dominates\nr=0.941 (turquoise)',
     '#00695C', '🕸'),
]

cols = 3
rows = 2
card_w = 4.5
card_h = 2.2
x_start = 0.6
y_start = 3.8

for idx, (title, detail, color, icon) in enumerate(findings):
    row = idx // cols
    col = idx % cols
    x = x_start + col * 5.1
    y = y_start - row * 2.7

    # Card
    card = FancyBboxPatch((x, y), card_w, card_h,
                           boxstyle="round,pad=0.12", facecolor='white',
                           edgecolor=color, linewidth=2.5)
    ax.add_patch(card)

    # Left color bar
    bar = FancyBboxPatch((x, y), 0.15, card_h,
                          boxstyle="round,pad=0", facecolor=color, alpha=0.8)
    ax.add_patch(bar)

    # Title
    ax.text(x + 0.5, y + card_h - 0.55, title, fontsize=12, fontweight='bold',
            ha='left', va='center', color=color)

    # Detail
    ax.text(x + 0.5, y + 0.55, detail, fontsize=10,
            ha='left', va='center', color='#555555')

plt.tight_layout(pad=0)
plt.savefig('assets/key_findings.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("Saved assets/key_findings.png")


print("\nAll README images generated successfully!")
print("  - assets/banner.png")
print("  - assets/pipeline_overview.png")
print("  - assets/key_findings.png")
