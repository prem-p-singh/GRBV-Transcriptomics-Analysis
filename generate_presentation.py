#!/usr/bin/env python3
"""
Generate a publication-quality PowerPoint presentation for GRBV Transcriptomics study.
Embeds all 36 figures with organized section slides, speaker notes, and key findings.
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
import os

# ============================================================
# CONFIGURATION
# ============================================================
WORK_DIR = '/Users/prempratapsingh/Data_Analysis'
OUTPUT = os.path.join(WORK_DIR, 'GRBV_Transcriptomics_Presentation.pptx')

# Color palette
DARK_BLUE = RGBColor(0x1B, 0x3A, 0x5C)
ACCENT_BLUE = RGBColor(0x2E, 0x86, 0xAB)
ACCENT_RED = RGBColor(0xD6, 0x27, 0x28)
ACCENT_GREEN = RGBColor(0x2C, 0xA0, 0x2C)
ACCENT_ORANGE = RGBColor(0xFF, 0x7F, 0x0E)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT_GRAY = RGBColor(0xF0, 0xF0, 0xF0)
MED_GRAY = RGBColor(0x66, 0x66, 0x66)
DARK_GRAY = RGBColor(0x33, 0x33, 0x33)

SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)

prs = Presentation()
prs.slide_width = SLIDE_W
prs.slide_height = SLIDE_H

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def add_background(slide, color=WHITE):
    """Set slide background color."""
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_shape_bg(slide, left, top, width, height, color, alpha=None):
    """Add a colored rectangle shape."""
    shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape


def add_text_box(slide, left, top, width, height, text, font_size=18,
                 color=DARK_GRAY, bold=False, alignment=PP_ALIGN.LEFT, font_name='Calibri'):
    """Add a text box with formatted text."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.color.rgb = color
    p.font.bold = bold
    p.font.name = font_name
    p.alignment = alignment
    return txBox


def add_bullet_slide_content(slide, left, top, width, height, bullets, font_size=16,
                              color=DARK_GRAY, spacing=Pt(6)):
    """Add bulleted text to a slide."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, bullet in enumerate(bullets):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.text = bullet
        p.font.size = Pt(font_size)
        p.font.color.rgb = color
        p.font.name = 'Calibri'
        p.space_after = spacing
        p.level = 0
    return txBox


def add_image_safe(slide, img_path, left, top, width=None, height=None):
    """Add image if it exists."""
    full = os.path.join(WORK_DIR, img_path) if not os.path.isabs(img_path) else img_path
    if os.path.exists(full):
        kwargs = {'left': left, 'top': top}
        if width:
            kwargs['width'] = width
        if height:
            kwargs['height'] = height
        slide.shapes.add_picture(full, **kwargs)
        return True
    return False


def make_title_slide(title, subtitle="", notes=""):
    """Create a styled title/section slide."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
    add_background(slide, WHITE)

    # Left accent bar
    add_shape_bg(slide, Inches(0), Inches(0), Inches(0.4), SLIDE_H, DARK_BLUE)

    # Title
    add_text_box(slide, Inches(0.8), Inches(2.0), Inches(11.5), Inches(2.0),
                 title, font_size=36, color=DARK_BLUE, bold=True)

    if subtitle:
        add_text_box(slide, Inches(0.8), Inches(4.0), Inches(11.5), Inches(1.5),
                     subtitle, font_size=20, color=MED_GRAY)

    # Bottom line
    add_shape_bg(slide, Inches(0.8), Inches(6.8), Inches(11.5), Inches(0.04), ACCENT_BLUE)

    if notes:
        slide.notes_slide.notes_text_frame.text = notes
    return slide


def make_section_slide(section_title, section_number=""):
    """Create a section divider slide."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    add_background(slide, DARK_BLUE)

    if section_number:
        add_text_box(slide, Inches(0.8), Inches(2.0), Inches(2), Inches(1.5),
                     section_number, font_size=60, color=ACCENT_BLUE, bold=True)

    add_text_box(slide, Inches(0.8), Inches(3.2), Inches(11.5), Inches(2.0),
                 section_title, font_size=36, color=WHITE, bold=True)

    add_shape_bg(slide, Inches(0.8), Inches(5.5), Inches(4.0), Inches(0.05), ACCENT_BLUE)
    return slide


def make_figure_slide(title, img_file, bullets=None, notes=""):
    """Create a slide with a figure and optional bullet points."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    add_background(slide, WHITE)

    # Top bar
    add_shape_bg(slide, Inches(0), Inches(0), SLIDE_W, Inches(0.9), DARK_BLUE)
    add_text_box(slide, Inches(0.5), Inches(0.12), Inches(12), Inches(0.7),
                 title, font_size=22, color=WHITE, bold=True)

    if bullets:
        # Image on left, bullets on right
        add_image_safe(slide, img_file, Inches(0.3), Inches(1.1), width=Inches(8.0))
        add_bullet_slide_content(slide, Inches(8.6), Inches(1.3), Inches(4.4), Inches(5.5),
                                  bullets, font_size=14, color=DARK_GRAY)
    else:
        # Image centered
        add_image_safe(slide, img_file, Inches(0.5), Inches(1.1), width=Inches(12.0))

    if notes:
        slide.notes_slide.notes_text_frame.text = notes
    return slide


def make_two_figure_slide(title, img1, img2, notes=""):
    """Create a slide with two figures side by side."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    add_background(slide, WHITE)

    add_shape_bg(slide, Inches(0), Inches(0), SLIDE_W, Inches(0.9), DARK_BLUE)
    add_text_box(slide, Inches(0.5), Inches(0.12), Inches(12), Inches(0.7),
                 title, font_size=22, color=WHITE, bold=True)

    add_image_safe(slide, img1, Inches(0.2), Inches(1.1), width=Inches(6.3))
    add_image_safe(slide, img2, Inches(6.8), Inches(1.1), width=Inches(6.3))

    if notes:
        slide.notes_slide.notes_text_frame.text = notes
    return slide


def make_key_numbers_slide(title, numbers_dict):
    """Create a slide with key statistics in colored boxes."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    add_background(slide, WHITE)

    add_shape_bg(slide, Inches(0), Inches(0), SLIDE_W, Inches(0.9), DARK_BLUE)
    add_text_box(slide, Inches(0.5), Inches(0.12), Inches(12), Inches(0.7),
                 title, font_size=22, color=WHITE, bold=True)

    colors = [ACCENT_BLUE, ACCENT_RED, ACCENT_GREEN, ACCENT_ORANGE,
              RGBColor(0x94, 0x67, 0xBD), RGBColor(0x8C, 0x56, 0x4B)]

    items = list(numbers_dict.items())
    n_cols = min(len(items), 3)
    n_rows = (len(items) + n_cols - 1) // n_cols

    box_w = Inches(3.8)
    box_h = Inches(1.8)
    start_x = Inches(0.5)
    start_y = Inches(1.3)
    gap_x = Inches(0.35)
    gap_y = Inches(0.3)

    for idx, (number, label) in enumerate(items):
        row = idx // n_cols
        col = idx % n_cols
        x = start_x + col * (box_w + gap_x)
        y = start_y + row * (box_h + gap_y)

        box = add_shape_bg(slide, x, y, box_w, box_h, colors[idx % len(colors)])

        # Number
        add_text_box(slide, x + Inches(0.2), y + Inches(0.2), box_w - Inches(0.4), Inches(0.9),
                     str(number), font_size=36, color=WHITE, bold=True, alignment=PP_ALIGN.CENTER)
        # Label
        add_text_box(slide, x + Inches(0.2), y + Inches(1.0), box_w - Inches(0.4), Inches(0.7),
                     label, font_size=14, color=WHITE, alignment=PP_ALIGN.CENTER)

    return slide


# ============================================================
# BUILD PRESENTATION
# ============================================================
print("Building GRBV Transcriptomics Presentation...")

# ---- SLIDE 1: Title ----
make_title_slide(
    "Comprehensive Transcriptomic Analysis of\nGrapevine Red Blotch Virus (GRBV) Infection",
    "Vitis vinifera cv. Merlot | RNA-Seq Across Infection Duration,\nMaturity Stages, Seasons, and Environmental Context",
    notes="Presenter: [Your Name]\nInstitution: [Your Institution]\nDate: [Date]"
)

# ---- SLIDE 2: Overview / Key Numbers ----
make_key_numbers_slide("Study at a Glance", {
    "327": "Samples (RB+ infected & RB- controls)",
    "14,604": "Genes analyzed after quality filtering",
    "6": "Core research questions addressed",
    "36": "Publication-quality figures generated",
    "115": "Heat-stress genes characterized",
    "12": "Universal biomarkers identified",
})

# ---- SLIDE 3: Objectives ----
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_background(slide, WHITE)
add_shape_bg(slide, Inches(0), Inches(0), SLIDE_W, Inches(0.9), DARK_BLUE)
add_text_box(slide, Inches(0.5), Inches(0.12), Inches(12), Inches(0.7),
             "Study Objectives", font_size=22, color=WHITE, bold=True)
add_bullet_slide_content(slide, Inches(0.8), Inches(1.3), Inches(11.5), Inches(5.5), [
    "Q1: How does infection duration (years since infection) affect the grapevine transcriptome?",
    "Q2: How does berry maturity stage modulate the GRBV response?",
    "Q3: Are GRBV-induced changes reproducible across seasons (2021 vs 2022)?",
    "Q4: Does the environment modulate GRBV transcriptomic effects?",
    "Q5: What are baseline environmental changes in uninfected vines?",
    "Q6: What is the transcriptomic trajectory as infection ages?",
    "PLUS: Biomarker discovery, GWAS/QTL cross-referencing, and heat-stress signature analysis",
], font_size=18)

# ---- SLIDE 4: Analysis Pipeline ----
make_figure_slide("Analysis Pipeline", "23_analysis_schematic.png",
    notes="10 Python scripts, 35 figures, 13 CSV outputs. Pipeline: raw counts → filtering → "
          "normalization → DE analysis → pathway enrichment → biomarker discovery → QTL cross-ref → heat-stress analysis.")

# ======== SECTION: GLOBAL OVERVIEW ========
make_section_slide("Global Transcriptome Overview", "01")

make_figure_slide("Global PCA: Variance Hierarchy", "01_global_pca_overview.png", [
    "Maturity stage is the dominant source of variance",
    "Season is the second major axis",
    "Infection status (RB+ vs RB-) shows subtle but consistent separation",
    "Infection duration has the smallest effect",
    "Hierarchy: Stage > Season > Infection > Duration",
], notes="PCA of 14,604 genes across 327 samples. log2(CPM+1) normalized.")

# ======== SECTION: CORE QUESTIONS ========
make_section_slide("Core Research Questions (Q1-Q6)", "02")

# Q1
make_figure_slide("Q1: Infection Duration Effects", "02_Q1_infection_years.png", [
    "Biphasic response pattern discovered",
    "Early responders: defense genes (NB-LRR, PR)",
    "Late responders: metabolism genes",
    "Year cohort analysis across 2016-2022",
    "Infection age correlates with transcriptomic shift",
])

make_figure_slide("Q1: Pathway Enrichment by Duration", "08_Q1_pathway_enrichment.png",
    notes="GO and KEGG enrichment for infection duration comparisons.")

# Q2
make_figure_slide("Q2: Maturity Stage Effects", "03_Q2_maturity_stages.png", [
    "Response intensity increases with ripening",
    "Post-veraison and Harvest show strongest effects",
    "Stage-specific DEG sets with limited overlap",
    "Developmental context shapes viral response",
])

make_figure_slide("Q2: Pathway Enrichment by Stage", "09_Q2_pathway_enrichment.png",
    notes="GO and KEGG enrichment for maturity stage comparisons.")

# Q3
make_figure_slide("Q3: Seasonal Reproducibility", "04_Q3_seasonal_effects.png", [
    "Cross-season fold-change correlation: r = 0.33-0.69",
    "Core response is conserved across years",
    "Season-specific genes reveal environmental interaction",
    "2021 cooler, 2022 hotter (+7.6C at Harvest)",
])

# Q4
make_figure_slide("Q4: Environmental Modulation", "05_Q4_environmental_effects.png", [
    "Environment MODULATES (not just adds to) GRBV effects",
    "Interaction terms significant for many pathways",
    "Hotter conditions alter both direction and magnitude",
    "Not a simple additive model — true interaction",
])

# Q5-Q6
make_two_figure_slide("Q5: Baseline Environment & Q6: Infection Trajectory",
    "06_Q5_uninfected_environment.png", "07_Q6_transcriptome_trajectory.png",
    notes="Q5: Environmental changes in RB- controls. Q6: Biphasic trajectory — acute defense to chronic metabolism.")

# ======== SECTION: SEASONAL MATURITY PROGRESSION ========
make_section_slide("KEY FINDING: Inverted Seasonal Trajectory", "03")

make_figure_slide("Seasonal x Maturity DEG Landscape", "15_seasonal_maturity_DEG_comparison.png", [
    "KEY FINDING: Inverted response kinetics",
    "2021 (cooler): Peak at Harvest (1,192 DEGs)",
    "2022 (hotter): Peak at Post-veraison (565 DEGs)",
    "Diagnostic timing must account for season",
    "Environment reshapes GRBV response trajectory",
], notes="This is arguably the most important finding — the seasonal inversion of response timing.")

make_figure_slide("Pathway Dynamics Across Season x Stage", "16_seasonal_maturity_pathway_enrichment.png",
    notes="Pathway enrichment grid across all 8 season x stage combinations.")

make_two_figure_slide("Seasonal Pathway Landscape & PCA Transitions",
    "17_seasonal_maturity_pathway_landscape.png", "18_seasonal_maturity_transitions_PCA.png")

make_figure_slide("Stage-by-Stage Pathway Detail", "19_seasonal_maturity_pathway_detail.png",
    notes="Detailed pathway progression through ripening for each season.")

# ======== SECTION: PATHWAY ENRICHMENT ========
make_section_slide("Systems-Level Pathway Analysis", "04")

make_figure_slide("Master Pathway Landscape", "14_master_pathway_landscape.png", [
    "Phenylpropanoid biosynthesis: consistently enriched",
    "Plant-pathogen interaction: defense activation",
    "ABC transporters: altered metabolite transport",
    "Starch/sucrose metabolism: sugar accumulation disrupted",
    "Photosynthesis affected at specific stages",
], notes="Master heatmap showing enrichment significance across all 20+ conditions.")

make_two_figure_slide("Q3 & Q4 Pathway Enrichment",
    "10_Q3_pathway_enrichment.png", "11_Q4_pathway_enrichment.png")

make_two_figure_slide("Q5 & Q6 Pathway Enrichment",
    "12_Q5_pathway_enrichment.png", "13_Q6_pathway_enrichment.png")

# ======== SECTION: ENVIRONMENT & TITER ========
make_section_slide("Environmental Context & Viral Titer Integration", "05")

make_figure_slide("Weather Conditions at Sampling", "30_weather_overview_sampling.png", [
    "2021: cooler growing season",
    "2022: significantly hotter (+7.6C at Harvest)",
    "Temperature tracked with daily weather data",
    "GDD accumulation differed between years",
])

make_figure_slide("Virus Titer Patterns", "31_virus_titer_analysis.png", [
    "Viral load varies by maturity stage",
    "Season affects titer dynamics",
    "Titer correlates with transcriptomic magnitude",
])

make_two_figure_slide("Integrated Temperature-Titer-DEG Analysis",
    "32_integrated_temp_titer_DEG.png", "33_temperature_titer_shift.png",
    notes="Correlations between temperature, viral titer, and DEG counts across conditions.")

# ======== SECTION: HEAT-STRESS ========
make_section_slide("Heat-Stress Signatures & Virus-Heat Disentanglement", "06")

make_key_numbers_slide("Heat-Stress Analysis: Key Numbers", {
    "115": "Heat-stress genes expressed (9 categories)",
    "57": "Heat genes differentially expressed due to GRBV",
    "29/36": "sHSPs downregulated at 2021 Harvest (RB+)",
    "23": "Heat genes with BOTH virus + environment effects",
    "+7.6C": "Temperature difference at Harvest (2022 vs 2021)",
    "r = 0.625": "HSF-sHSP module score correlation",
})

make_figure_slide("Heat-Stress Gene Expression Landscape", "34_heat_stress_gene_landscape.png", [
    "9 categories: sHSP, HSP70, HSP90, HSP100, HSF, GolS, LEA, antioxidant, GO:heat",
    "sHSPs show strongest stage- and season-dependent expression",
    "2022 Harvest: highest sHSP module score in controls (+1.64)",
    "GRBV infection REDUCES sHSP module by ~50%",
    "Virus suppresses host proteostasis network",
], notes="Plot 34: 4-panel overview of heat gene categories, DEG counts, correlations, module scores.")

make_figure_slide("Heat vs Virus Effect Disentanglement", "35_heat_vs_virus_disentanglement.png", [
    "At Harvest: 23 genes regulated by BOTH virus AND environment",
    "At Post-veraison: virus effect dominates (42 virus-only DEGs)",
    "FC scatter shows virus OPPOSES environmental heat induction",
    "HSF → sHSP regulatory axis confirmed (r=0.625, p=7.2e-37)",
    "GRBV may compromise vine heat tolerance",
], notes="Plot 35: Classification, FC scatter, gene profiles, HSF-sHSP module correlation.")

# ======== SECTION: BIOMARKERS ========
make_section_slide("Biomarker Discovery & Diagnostics", "07")

make_figure_slide("Biomarker Scoring & Classification", "20_biomarker_overview.png", [
    "Composite score: consistency + effect size + significance + direction + robustness",
    "12 universal biomarkers robust across ALL conditions",
    "3-gene panel achieves AUC = 1.000 (perfect classification)",
    "Top host biomarker: F-box/LRR protein (AUC = 0.880)",
    "35 early-detection + 23 harvest-specific biomarkers",
])

make_figure_slide("Biomarker Expression Profiles", "21_biomarker_profiles.png",
    notes="Expression profiles of top 16 biomarkers across all season x stage conditions.")

make_figure_slide("Biomarker Robustness Analysis", "22_biomarker_robustness.png", [
    "Fold-change consistency across seasons verified",
    "AUC stability across cross-validation folds",
    "Significance matrix confirms multi-condition robustness",
    "RT-qPCR validation candidates identified",
])

# ======== SECTION: GWAS/QTL ========
make_section_slide("GWAS/QTL Cross-Reference & Breeding Implications", "08")

make_figure_slide("Biomarker-QTL Chromosome Map", "24_biomarker_QTL_chromosome_map.png", [
    "60% of biomarkers co-localize with known QTLs (1 Mb window)",
    "40% fall directly inside QTL intervals",
    "18 biomarkers in disease-resistance clusters (Rpv, Ren, Run)",
    "8 biomarkers near berry-quality loci",
    "Sugar QTL link: RNA-dep RNA polymerase → sugar accumulation",
])

make_figure_slide("QTL Cross-Reference Summary", "25_biomarker_QTL_crossref_summary.png",
    notes="Overlap statistics, category breakdown, functional interpretation of co-localization.")

# ======== SECTION: ROBUSTNESS ========
make_section_slide("Analytical Robustness & Integration", "09")

make_two_figure_slide("Pooled vs Stratified Comparison & Overlap Heatmap",
    "26_pooled_vs_stratified_comparison.png", "27_overlap_connectivity_heatmap.png")

make_two_figure_slide("Cross-Question Integration & Summary Dashboard",
    "28_Q1Q2_vs_Q3Q4_crossreference.png", "29_comprehensive_summary_dashboard.png",
    notes="Cross-referencing findings across all research questions. 908 universal DEGs confirmed.")

# ======== SECTION: INTEGRATED BIOLOGY DIAGRAM ========
make_section_slide("Integrated Molecular Model", "10")

make_figure_slide("GRBV Molecular Disruption Cascade", "36_GRBV_biology_diagram.png", [
    "Data-driven infographic integrating all findings",
    "DEG counts per season x stage with viral titer overlay",
    "Phenylpropanoid, flavonoid & stilbenoid pathways suppressed",
    "sHSP heat-protection network reduced by ~55%",
    "3-gene diagnostic panel: AUC = 1.000",
    "60% biomarker-QTL co-localization with berry quality loci",
], notes="Figure 36: Comprehensive biology diagram synthesizing all analyses into a molecular disruption model.")

# ======== SECTION: CONCLUSIONS ========
make_section_slide("Conclusions & Future Directions", "11")

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_background(slide, WHITE)
add_shape_bg(slide, Inches(0), Inches(0), SLIDE_W, Inches(0.9), DARK_BLUE)
add_text_box(slide, Inches(0.5), Inches(0.12), Inches(12), Inches(0.7),
             "Key Conclusions", font_size=22, color=WHITE, bold=True)
add_bullet_slide_content(slide, Inches(0.5), Inches(1.2), Inches(6.0), Inches(5.8), [
    "1. GRBV effects are profoundly stage-dependent, peaking at Post-veraison and Harvest",
    "2. Seasonal response trajectory is INVERTED: 2021 peaks at Harvest, 2022 peaks at Post-veraison",
    "3. GRBV suppresses heat-stress genes (sHSPs), potentially compromising vine thermal tolerance",
    "4. 12 universal biomarkers + 3-gene panel (AUC=1.0) enable robust diagnostics",
    "5. 60% of biomarkers co-localize with known disease-resistance and berry-quality QTLs",
    "6. Sugar accumulation QTL link explains primary GRBV symptom (delayed ripening, reduced Brix)",
], font_size=16, spacing=Pt(10))

# Right side: future directions
add_text_box(slide, Inches(7.0), Inches(1.2), Inches(5.5), Inches(0.5),
             "Future Directions", font_size=20, color=DARK_BLUE, bold=True)
add_bullet_slide_content(slide, Inches(7.0), Inches(1.8), Inches(5.8), Inches(4.5), [
    "RT-qPCR validation of top biomarkers in independent vineyards",
    "Functional validation of sHSP suppression mechanism (e.g., VIGS knockdown)",
    "Expand to additional cultivars (Cabernet Sauvignon, Zinfandel, Pinot Noir)",
    "Integrate with metabolomics (anthocyanin, sugar profiling)",
    "Develop field-deployable molecular diagnostic assay",
    "Test biomarker-QTL candidates in breeding populations",
], font_size=14, spacing=Pt(6))

# ---- Final acknowledgments slide ----
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_background(slide, DARK_BLUE)
add_text_box(slide, Inches(0.8), Inches(2.5), Inches(11.5), Inches(2.0),
             "Thank You", font_size=48, color=WHITE, bold=True, alignment=PP_ALIGN.CENTER)
add_text_box(slide, Inches(0.8), Inches(4.5), Inches(11.5), Inches(1.0),
             "Questions & Discussion", font_size=24, color=ACCENT_BLUE, alignment=PP_ALIGN.CENTER)
add_shape_bg(slide, Inches(4.5), Inches(5.8), Inches(4.3), Inches(0.04), ACCENT_BLUE)
add_text_box(slide, Inches(0.8), Inches(6.0), Inches(11.5), Inches(1.0),
             "[Your Name]  |  [Institution]  |  [Email]",
             font_size=16, color=LIGHT_GRAY, alignment=PP_ALIGN.CENTER)

# ============================================================
# SAVE
# ============================================================
prs.save(OUTPUT)
print(f"\nPresentation saved: {OUTPUT}")
print(f"Total slides: {len(prs.slides)}")
print("\nSlide summary:")
print("  1. Title slide")
print("  2. Study at a glance (key numbers)")
print("  3. Study objectives")
print("  4. Analysis pipeline (Fig 23)")
print("  5. Section: Global Overview")
print("  6. Global PCA (Fig 1)")
print("  7. Section: Core Questions Q1-Q6")
print("  8-13. Q1-Q6 figures with key findings")
print("  14. Section: Inverted Seasonal Trajectory")
print("  15-18. Seasonal maturity progression (Figs 15-19)")
print("  19. Section: Pathway Analysis")
print("  20-22. Pathway enrichment (Figs 10-14)")
print("  23. Section: Environment & Titer")
print("  24-27. Weather, titer, integration (Figs 30-33)")
print("  28. Section: Heat-Stress")
print("  29-31. Heat-stress key numbers & figures (Figs 34-35)")
print("  32. Section: Biomarkers")
print("  33-35. Biomarker scoring, profiles, robustness (Figs 20-22)")
print("  36. Section: GWAS/QTL")
print("  37-38. Chromosome map & cross-ref (Figs 24-25)")
print("  39. Section: Robustness")
print("  40-41. Pooled vs stratified, summary dashboard (Figs 26-29)")
print("  42. Section: Integrated Molecular Model")
print("  43. GRBV biology diagram (Fig 36)")
print("  44. Section: Conclusions")
print("  45. Key conclusions & future directions")
print("  46. Thank you / Q&A")
