#!/usr/bin/env python
"""Assemble the full GRBV figure atlas: 10 Main composites (R/Arial,
DESeq2-current) + 55 source figures as Supplementary.

Each Supplementary figure is wrapped in a single-page PDF with the figure
PNG, its title, and a status badge indicating whether the underlying data
has been refreshed under the DESeq2 reanalysis or still reflects the
original Welch-t analysis.
"""
from __future__ import annotations

import io
from pathlib import Path

from pypdf import PdfWriter
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.units import inch
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas

HERE       = Path(__file__).resolve().parent
DATA_ROOT  = Path("/Users/prempratapsingh/Data_Analysis")
MAIN_DIR   = DATA_ROOT / "figures" / "main"
SUPP_DIR   = DATA_ROOT / "figures" / "supp"
SOURCE_DIR = DATA_ROOT
OUT        = DATA_ROOT / "figures" / "GRBV_Full_Figure_Atlas.pdf"

def resolve_supp(n: int, fname: str) -> Path:
    """Return the regenerated supp PNG if it exists in SUPP_DIR, else the
    original source PNG in DATA_ROOT. First try direct fname match, then
    fall back to leading two-digit-number prefix glob."""
    direct = SUPP_DIR / fname
    if direct.exists():
        return direct
    prefix = f"S{n:02d}_"
    for p in SUPP_DIR.glob(f"{prefix}*.png"):
        return p
    return SOURCE_DIR / fname

MAINS = [
    "Main01_overview", "Main02_de_landscape", "Main03_inverted_trajectory",
    "Main04_heat_stress", "Main05_mediation",  "Main06_biomarker",
    "Main07_wgcna_tf",   "Main08_chr_qtl",     "Main09_gsea",
    "Main10_biology_cascade",
]

# --- 55 source figures with title + data-currency status -------------
# DATA STATUS:
#   current  = data underlying figure unchanged by fix passes
#   DESeq2   = figure should be refreshed using de_redo/de_deseq2_extras DESeq2 results
#   reframed = figure is superseded by a Main composite (here for completeness)
SOURCE_FIGS = [
    ( 1, "01_global_pca_overview.png",                "Global PCA overview",                                "current",  "See Main 1 (panel B)"),
    ( 2, "02_Q1_infection_years.png",                 "Q1 infection-duration DEGs",                         "DESeq2",   "See Main 2; refit with de_deseq2_extras Q1a/Q1b"),
    ( 3, "03_Q2_maturity_stages.png",                 "Q2 maturity-stage DEGs",                             "DESeq2",   "See Main 2 + de_deseq2_voom Q2"),
    ( 4, "04_Q3_seasonal_effects.png",                "Q3 seasonal effects in RB+",                         "DESeq2",   "See Main 2 + de_deseq2_voom Q3"),
    ( 5, "05_Q4_environmental_effects.png",           "Q4 environmental modulation",                        "DESeq2",   "See Main 2 + de_deseq2_extras Q4d"),
    ( 6, "06_Q5_uninfected_environment.png",          "Q5 baseline seasonal change in RB-",                 "DESeq2",   "See de_deseq2_voom Q5"),
    ( 7, "07_Q6_transcriptome_trajectory.png",        "Q6 pooled transcriptome trajectory",                 "DESeq2",   "See Main 2 panel B"),
    ( 8, "08_Q1_pathway_enrichment.png",              "Q1 pathway enrichment",                              "DESeq2",   "Recompute on Q1a/Q1b DEG sets"),
    ( 9, "09_Q2_pathway_enrichment.png",              "Q2 pathway enrichment",                              "DESeq2",   "Recompute on Q2 DEG sets"),
    (10, "10_Q3_pathway_enrichment.png",              "Q3 pathway enrichment",                              "DESeq2",   "Recompute on Q3 DEG sets"),
    (11, "11_Q4_pathway_enrichment.png",              "Q4 pathway enrichment",                              "DESeq2",   "Recompute on Q4d DEG sets"),
    (12, "12_Q5_pathway_enrichment.png",              "Q5 pathway enrichment",                              "DESeq2",   "Recompute on Q5 DEG sets"),
    (13, "13_Q6_pathway_enrichment.png",              "Q6 pathway enrichment",                              "DESeq2",   "Recompute on Q6 DEG sets"),
    (14, "14_master_pathway_landscape.png",           "Master pathway landscape (all questions)",           "DESeq2",   "Aggregate of S8-S13"),
    (15, "15_seasonal_maturity_DEG_comparison.png",   "Inverted seasonal trajectory (key)",                 "DESeq2",   "See Main 3 panel A"),
    (16, "16_seasonal_maturity_pathway_enrichment.png","Stage \u00d7 season pathway enrichment grid",       "DESeq2",   "Recompute from de_deseq2_extras Q2x"),
    (17, "17_seasonal_maturity_pathway_landscape.png","Per-season pathway landscape heatmaps",              "DESeq2",   "Recompute"),
    (18, "18_seasonal_maturity_transitions_PCA.png",  "Stage transitions PCA",                              "current",  ""),
    (19, "19_seasonal_maturity_pathway_detail.png",   "Pathway detail per stage",                           "DESeq2",   "Recompute"),
    (20, "20_biomarker_overview.png",                 "Biomarker overview (scores, AUC, panel)",            "reframed", "See Main 6"),
    (21, "21_biomarker_profiles.png",                 "Top-16 biomarker expression profiles",               "current",  ""),
    (22, "22_biomarker_robustness.png",               "Biomarker FC consistency + AUC stability",           "reframed", "See Main 6 panel A/B"),
    (23, "23_analysis_schematic.png",                 "Pipeline analysis schematic",                        "current",  ""),
    (24, "24_biomarker_QTL_chromosome_map.png",       "Biomarker positions on chromosome map",              "current",  "See Main 8 panel B"),
    (25, "25_biomarker_QTL_crossref_summary.png",     "QTL overlap summary (no permutation null)",          "reframed", "See Main 8 panel C"),
    (26, "26_pooled_vs_stratified_comparison.png",    "Pooled vs stratified DE concordance",                "DESeq2",   "Recompute with de_deseq2_voom"),
    (27, "27_overlap_connectivity_heatmap.png",       "DEG overlap across questions",                       "DESeq2",   "Recompute with de_redo/ + de_extras"),
    (28, "28_Q1Q2_vs_Q3Q4_crossreference.png",        "Cross-question DEG cross-reference",                 "DESeq2",   "Recompute"),
    (29, "29_comprehensive_summary_dashboard.png",    "Integrated findings dashboard",                      "DESeq2",   "Recompute"),
    (30, "30_weather_overview_sampling.png",          "Weather overview at sampling",                       "current",  ""),
    (31, "31_virus_titer_analysis.png",               "Viral titer dynamics",                               "current",  "See Main 3 panel C"),
    (32, "32_integrated_temp_titer_DEG.png",          "Integrated temperature-titer-DEG",                   "DESeq2",   "Recompute with de_deseq2_extras"),
    (33, "33_temperature_titer_shift.png",            "Temperature-titer correlations",                     "current",  ""),
    (34, "34_heat_stress_gene_landscape.png",         "Heat-gene landscape across categories",              "current",  "See Main 4 panel A"),
    (35, "35_heat_vs_virus_disentanglement.png",      "Heat vs virus disentanglement (HSF-sHSP axis)",      "current",  "See Main 4 panel B"),
    (36, "36_GRBV_biology_diagram.png",               "Integrated biology infographic",                     "current",  "See Main 10 (cleaner cascade)"),
    (37, "37_gsea_enrichment_plots.png",              "GSEA running enrichment plots",                      "current",  ""),
    (38, "38_gsea_vs_fisher_concordance.png",         "GSEA vs Fisher ORA concordance",                     "current",  ""),
    (39, "39_gsea_nes_heatmap.png",                   "GSEA NES heatmap",                                   "current",  "See Main 9 panel A"),
    (40, "40_variance_decomposition.png",             "ANOVA \u03b7\u00b2 variance decomposition",          "current",  "See Main 1 panels C/D"),
    (41, "41_interaction_gene_sets.png",              "Interaction gene-set counts and overlap",            "current",  ""),
    (42, "42_interaction_profiles.png",               "Top Inf \u00d7 Sea interaction profiles",            "current",  ""),
    (43, "43_tf_family_heatmap.png",                  "TF family %DE and mean FC heatmap",                  "current",  "See Main 7 panel B"),
    (44, "44_tf_enrichment_balloon.png",              "TF family DEG enrichment balloon plot",              "current",  ""),
    (45, "45_tf_expression_boxplots.png",             "Top TF gene expression boxplots",                    "current",  ""),
    (46, "46_ml_method_comparison.png",               "Multi-method ML AUC comparison",                     "reframed", "See Main 6"),
    (47, "47_stability_selection.png",                "LASSO stability selection (100 subsamples)",         "reframed", "See Main 6 panel C"),
    (48, "48_consensus_panel.png",                    "Consensus biomarker panel + ROC",                    "reframed", "See Main 6"),
    (49, "49_mediation_diagrams.png",                 "Mediation path diagrams (top 6)",                    "current",  "See Main 5 panel A (MBF1c)"),
    (50, "50_direct_vs_indirect_effects.png",         "Direct vs indirect effect scatter",                  "current",  "See Main 5 panels B-D"),
    (51, "51_chromosome_deg_ideogram.png",            "Chromosome DEG-density ideogram",                    "current",  ""),
    (52, "52_deg_manhattan_plot.png",                 "DEG enrichment Manhattan",                           "current",  "See Main 8 panel A"),
    (53, "53_module_trait_correlations.png",          "WGCNA module-trait correlation heatmap",             "current",  "See Main 7 panel A"),
    (54, "54_module_eigengene_profiles.png",          "WGCNA module eigengene profiles",                    "current",  ""),
    (55, "55_hub_gene_heatmaps.png",                  "WGCNA hub gene expression heatmaps",                 "current",  "See Main 7 panel C"),
    (56, "S56_viral_expression_atlas.png",            "Viral transcriptome atlas (5 ORFs \u00d7 8 cells \u00d7 cohorts)", "regenerated", "Per-ORF log-CPM across stage/season/cohort"),
    (57, "S57_viral_titer_coherence.png",             "Virus-titer coherence check",                        "regenerated", "RNA-Seq CPM vs RT-qPCR titer, r = 0.07 overall"),
    (58, "S58_per_ORF_host_correlates.png",           "Per-ORF host-gene correlate landscape",              "regenerated", "Top 15 host Spearman correlates per viral gene"),
    (59, "S59_viral_host_scatter_grid.png",           "Viral-host co-expression grid",                      "regenerated", "5 viral ORFs \u00d7 5 host anchors (MBF1c, SPS1, etc.)"),
    (60, "S60_viral_biology_diagram.png",             "Viral-gene function atlas (biology explainer)",      "regenerated", "Genome schematic + function cards + virus-host impact map"),
    (61, "S61_biochem_validation.png",                "Harvest biochem x transcriptome validation",         "regenerated", "Per-season + partial rho for SPS1, MBF1c, LOC100266309, sHSP, titer vs Brix / Malic"),
    (62, "S62_biochem_cohort_trajectory.png",         "Harvest biochem vs years-since-infection",           "regenerated", "Brix / pH / TA / Malic vs YSI per season; 2021 Brix recovery +0.10 \u00b0Bx / yr (r=0.76, p=0.0045)"),
    (63, "S63_biochem_season_infection.png",          "Season x Infection ANOVA on juice chemistry",        "regenerated", "TA interaction \u03b7\u00b2=0.77; cool 2021 +0.62 g/L vs hot 2022 -0.79 g/L; Brix 2.6x stronger in 2022"),
    (64, "S64_tier1_overlap_summary.png",             "Tier 1 DEG overlap: size, intersections, enrichment","regenerated", "Signature sizes (log-x), GRBV intersections (35 unique / 13 Botrytis-shared / 7 multi), Fisher q=2.1e-5 vs Blanco-Ulate Zinfandel"),
    (65, "S65_tier1_fc_concordance.png",              "Tier 1 log2FC Spearman concordance heatmap",         "regenerated", "GRBV (Merlot) x Blanco-Ulate GRBV (Zinfandel) rho=0.47; Blanco-Ulate x Song GLRaV-3 rho=-0.21; other cells \u2248 0"),
    (66, "S66_disease_progression_cube.png",          "Disease progression cube (RB+)",                     "regenerated", "Cohort x season x stage cube for titer / sHSP / GRBV signature / Brix"),
    (67, "S67_rbminus_baseline_cube.png",             "RB(-) baseline cube",                                "regenerated", "Healthy comparator that makes the passive-failure reframe visible"),
    (68, "S68_enrichment_season_stage.png",           "Season x stage GSEA + season-flip catalog",          "regenerated", "Top variable pathways + 20-pathway sign-flip catalog + 8-family trajectory"),
    (69, "S69_enrichment_cohort.png",                 "Per-cohort hyper-enrichment",                        "regenerated", "65x sHSP at 2016 (chronic) vs 5x sHSP at 2022 (acute); 7-year temporal architecture"),
    (70, "S70_literature_curated_pathways.png",       "Literature-curated 90-pathway NES heatmap",          "regenerated", "Ribosome / heat-shock / sugar / phenylpropanoid / immunity per stage x season"),
    (71, "S71_transcriptional_response_strength.png", "PCA-based response strength per cell",               "regenerated", "Peak phase 2020 cohort reaches PC-distance 21.9 at 2021 Harvest"),
    (72, "S72_variance_decomposition_ternary.png",    "Per-gene variance ternary (Singh 2025 framework)",   "regenerated", "Dev / Viral / Env split; sHSPs lie on Dev-Env edge with viral share \u2248 baseline"),
    (73, "S73_viral_orf_extensions.png",              "Viral ORF extensions: phenotype + module",           "regenerated", "V2 adds Brix/Malic prediction beyond titer; per-ORF x 13-module heatmap"),
    (74, "S74_pathway_per_gene.png",                  "Phenylpropanoid + PR + cell-wall per-gene model",    "regenerated", "31% Phenylpropanoid / 20% PR / 26% Cell-wall infection-significant; 45-51% 3-way"),
    (75, "S75_carbon_translocation.png",              "Sugar-transport module (SWEET/SUT/STP/CWINV)",       "regenerated", "77 genes; module rho_stage=0.94 but no main season or cohort effect (null hypothesis test of Martinez-Luscher 2018)"),
    (76, "S76_vine_heterogeneity.png",                "Vine-level tolerant vs susceptible",                 "regenerated", "44 vines: 13 tolerant (24.2 Bx) vs 13 susceptible (23.6 Bx) within one Merlot vineyard"),
    (77, "S77_tf_activity.png",                       "Regulon-based TF activity DE",                       "regenerated", "236 TFs scored; LOC100262652 padj=1.4e-64; MYB30 down padj=8.9e-10"),
    (78, "S78_segmented_cohort.png",                  "Segmented regression on cohort age (null result)",   "regenerated", "best_k=0 (linear) for all 15 responses; '3 phases' framing is descriptive only"),
    (79, "S79_bayesian_posterior.png",                "Bayesian posteriors on infection effect",            "regenerated", "RdRp1 +1.78, F-box +2.01, glycine-rich -1.33; MBF1c -0.25 [-0.81, +0.33], 19% in ROPE"),
    (80, "S80A_phenylpropanoid_PR_heatmap.png",       "S80A — Phenylpropanoid + PR proteins (91 genes)",    "regenerated", "Hierarchical clustering (Ward.D2); facet_nested season > stage; cohort years on x"),
    (81, "S80B_cell_wall_heatmap.png",                "S80B — Cell wall (161 genes, 5 sub-categories)",     "regenerated", "Pectin / Hemicellulose / Cellulose / Lignin / Expansin; only |log2FC|>=1 genes labelled"),
    (82, "S80C_cell_wall_focused.png",                "S80C — Cell wall genes with consistent GRBV response","regenerated", "21 genes (13 strict + 8 broad at |FC|>=0.45 in >=6/8); axis_nested by sub-category; cell text"),
    (83, "S81_sugar_full_heatmap.png",                "Carbon translocation full FC heatmap",               "regenerated", "33 sugar-transport genes x (season x cohort x stage); functional subclass ordering"),
    (84, "S84_synteny_pinot_merlot.png",              "Pinot Noir <-> Merlot cl. 343 synteny",              "regenerated", "JCVI MCScan: 41,248 anchors / 125 blocks; both Merlot haplotypes preserve colinearity with Pinot"),
    (85, "S85_salmon_reference_concordance.png",      "Direct reference-bias test: salmon re-quantification","regenerated", "Same 327 reads vs Pinot + Merlot; sig-in-both rho=0.95 / 100% sign-agree; marginal Jaccard 0.18"),
    (86, "S86_hormone_landscape.png",                 "Hormone-pathway response (7 hormones, 519 genes)",   "regenerated", "Q6 pooled covariate-adj DESeq2; no class enrichment, GA monotonic-down; cross-paper validation vs BU 2017"),
    (87, "S87_phenylpropanoid_kegg.png",              "Phenylpropanoid biosynthesis (KEGG vvi00940)",       "regenerated", "Pathview 8-cell split-box (pure stage order: 2021 Pre/Vér/Post/Harv | 2022); paralog mean log2FC; node.sum=mean"),
    (88, "S88_starch_sucrose_kegg.png",               "Starch and sucrose metabolism (vvi00500)",           "regenerated", "Sugar/Brix story; includes SPS1; GSEA 5/8 cond NES +0.89"),
    (89, "S89_plant_hormone_kegg.png",                "Plant hormone signal transduction (vvi04075)",       "regenerated", "Narrative tie-in to S86 hormone landscape (ethylene cross-paper / GA-down)"),
    (90, "S90_ribosome_kegg.png",                     "Ribosome (vvi03010)",                                "regenerated", "Translational shutdown; GSEA 6/8 cond NES -1.33 (most robust DOWN)"),
    (91, "S91_amino_acid_biosynth_kegg.png",          "Biosynthesis of amino acids (vvi01230)",             "regenerated", "GSEA 5/8 cond NES -1.95"),
    (92, "S92_proteasome_kegg.png",                   "Proteasome (vvi03050)",                              "regenerated", "Strongest |NES| in pathway set; GSEA 4/8 cond NES -2.43"),
    (93, "S93_hsf_mbf1c_shsp_axis.png",               "HSF -> MBF1c -> sHSP axis (custom flow)",            "new",         "Three-column ggplot diagram; passive-failure model (Thread 10); per-node 8-cell Q2x log2FC mini-bars"),
    (94, "S94_temp_titer_mediator_cascade.png",       "Temperature -> titer -> mediator cascade",           "new",         "MBF1c / SPS1 / LOC100266309 mediator nodes; edges from mediation_results.csv (tv=2000)"),
    (95, "S95_viral_orf_module_rewiring.png",         "Viral ORFs -> WGCNA module rewiring",                "new",         "Bipartite: 6 ORFs x 13 modules; rule 13 Spearman; FDR<0.05 or |rho|>=0.3"),
    (96, "S96_rumbaugh2022_cross_cultivar.png",       "Cross-cultivar replication vs Rumbaugh 2022 CS",     "new",         "3 panels: DCL2 veraison-peak trajectory (their headline mechanism, padj 5.5e-11 in our 2021); per-network Spearman concordance (S1 Ca-binding, S2 DCL2-gain, S3 DCL2-loss; rootstocks separate); Q1b cohort gradient (peak match at cohort2020 = 1-2 y post-infection)"),
    (97, "S97_rumbaugh2023_within_cultivar.png",      "Within-cultivar replication vs Rumbaugh 2023 Merlot","new",         "4 panels: per-gene harvest scatter (rho = 0.79 across 37 mapped genes); per-category Spearman bars (cell wall rho = 0.69 with 100% sign agreement); per-(category, stage) heatmap; Q1b cohort gradient. Cell wall + phenylpropanoid biology is broadly cohort-age-portable; PR proteins dip at cohort2020 (rho = 0.10)"),
]

PAGE_W, PAGE_H = LETTER  # 8.5 x 11 inches
MARGIN = 0.5 * inch
STATUS_COLOR = {
    "current":     (0.15, 0.55, 0.30),  # green
    "DESeq2":      (0.85, 0.30, 0.20),  # red (pending)
    "regenerated": (0.10, 0.40, 0.70),  # blue (done)
    "reframed":    (0.50, 0.30, 0.65),  # purple
    "new":         (0.10, 0.40, 0.70),  # blue (newly created under DESeq2 pipeline)
}

def render_supp_page(num: int, png_path: Path, title: str,
                     status: str, note: str) -> bytes:
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=LETTER)
    c.setTitle(title)

    # Header
    c.setFont("Helvetica-Bold", 14)
    c.drawString(MARGIN, PAGE_H - MARGIN - 14,
                 f"Supplementary Figure S{num}.  {title}")

    # Status badge
    badge = {
        "current":     "Data current (kept as-is, restyle pending)",
        "DESeq2":      "DATA UPDATE NEEDED \u2014 refit with DESeq2 (de_redo/ + de_extras)",
        "regenerated": "Regenerated under DESeq2, R/ggplot2, Arial",
        "reframed":    "Reframed in Main composite (see note)",
        "new":         "New figure produced under DESeq2 pipeline",
    }[status]
    c.setFont("Helvetica-Bold", 9)
    r, g, b = STATUS_COLOR[status]
    c.setFillColorRGB(r, g, b)
    c.drawString(MARGIN, PAGE_H - MARGIN - 30, f"Status: {badge}")
    c.setFillColorRGB(0, 0, 0)

    if note:
        c.setFont("Helvetica-Oblique", 9)
        c.drawString(MARGIN, PAGE_H - MARGIN - 44, f"Note: {note}")

    # Figure image
    if png_path.exists():
        img = ImageReader(str(png_path))
        iw, ih = img.getSize()
        avail_w = PAGE_W - 2 * MARGIN
        avail_h = PAGE_H - 2 * MARGIN - 70
        scale = min(avail_w / iw, avail_h / ih)
        w = iw * scale; h = ih * scale
        x = (PAGE_W - w) / 2
        y = MARGIN + 10
        c.drawImage(img, x, y, width=w, height=h, preserveAspectRatio=True,
                    anchor='c')
    else:
        c.setFont("Helvetica", 11)
        c.setFillColorRGB(0.5, 0, 0)
        c.drawString(MARGIN, PAGE_H / 2,
                     f"[file not found: {png_path.name}]")

    c.setFont("Helvetica", 7)
    c.setFillColorRGB(0.4, 0.4, 0.4)
    c.drawString(MARGIN, MARGIN - 4,
                 f"GRBV transcriptomics atlas \u2014 page S{num}")
    c.showPage()
    c.save()
    return buf.getvalue()

# --------- assemble -------------------------------------------------
# Main-figure page renderer — builds a fresh PDF page from the PNG.
# (Main figures are no longer shipped as PDF; PNG + SVG only.)
def render_main_page(stem: str) -> bytes:
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=LETTER)
    c.setTitle(stem.replace("_", " "))
    png = MAIN_DIR / f"{stem}.png"
    if not png.exists():
        raise SystemExit(f"missing main composite PNG: {png}")
    img = ImageReader(str(png))
    iw, ih = img.getSize()
    avail_w = PAGE_W - 2 * MARGIN
    avail_h = PAGE_H - 2 * MARGIN - 20
    scale = min(avail_w / iw, avail_h / ih)
    w = iw * scale; h = ih * scale
    x = (PAGE_W - w) / 2
    y = (PAGE_H - h) / 2
    c.drawImage(img, x, y, width=w, height=h, preserveAspectRatio=True, anchor='c')
    c.setFont("Helvetica", 8)
    c.setFillGray(0.4)
    c.drawString(MARGIN, MARGIN, stem.replace("_", " "))
    c.showPage()
    c.save()
    return buf.getvalue()

writer = PdfWriter()

# 10 Main composites — rendered from PNG
for stem in MAINS:
    page_pdf = render_main_page(stem)
    writer.append(io.BytesIO(page_pdf))

# 55 Supplementary pages
for num, fname, title, status, note in SOURCE_FIGS:
    img_path = resolve_supp(num, fname)
    regen = img_path.parent.name == "supp"
    if regen and status == "DESeq2":
        status = "regenerated"  # upgrade badge
    page_pdf = render_supp_page(num, img_path, title, status, note)
    writer.append(io.BytesIO(page_pdf))

with open(OUT, "wb") as fh:
    writer.write(fh)
total = len(MAINS) + len(SOURCE_FIGS)
print(f"wrote {OUT} ({total} pages: {len(MAINS)} main + {len(SOURCE_FIGS)} supp)")
