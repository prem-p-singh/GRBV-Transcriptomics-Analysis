# GRBV Transcriptomics Analysis

![Python](https://img.shields.io/badge/Python-3.6%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen)

A comprehensive RNA-Seq transcriptomics pipeline for studying **Grapevine Red Blotch Virus (GRBV)** effects on *Vitis vinifera* (Merlot). The analysis covers 327 field samples across 2 growing seasons and 4 maturity stages, producing 55 publication-quality figures, 26+ data tables, and automated reports.

---

## Table of Contents

- [Study Overview](#study-overview)
- [Repository Structure](#repository-structure)
- [Input Data Requirements](#input-data-requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output Files](#output-files)
- [Key Findings](#key-findings)
- [Citation](#citation)
- [License](#license)

---

## Study Overview

| Parameter | Details |
|-----------|---------|
| **Species** | *Vitis vinifera* cv. Merlot |
| **Pathogen** | Grapevine Red Blotch Virus (GRBV) |
| **Samples** | 327 total (264 RB+ infected, 63 RB- control) |
| **Genes** | 28,168 raw &rarr; 14,604 after filtering |
| **Seasons** | 2021 (cooler) and 2022 (hotter, +7.6 &deg;C at Harvest) |
| **Maturity Stages** | Pre-veraison, Veraison, Post-veraison, Harvest |
| **Infection Cohorts** | 7 years (2016-2022) |
| **DE Thresholds** | adj. p-value < 0.05, \|log2FC\| &ge; 1 |

The pipeline addresses six core research questions:
1. **Q1** - How does infection duration affect transcriptome responses?
2. **Q2** - Are GRBV effects maturity-stage dependent?
3. **Q3** - Do seasonal conditions modulate GRBV responses?
4. **Q4** - Does environment amplify or dampen viral effects?
5. **Q5** - What seasonal variation exists in uninfected vines?
6. **Q6** - How does transcriptome trajectory shift with infection age?

---

## Repository Structure

### Analysis Scripts (run in order)

| # | Script | Description |
|---|--------|-------------|
| 1 | `analysis.py` | Core differential expression analysis for 6 research questions (Welch's t-test, BH FDR) |
| 2 | `pathway_enrichment.py` | GO (BP/MF/CC) and KEGG pathway enrichment for all 6 questions |
| 3 | `seasonal_maturity_progression.py` | Season-stratified DE across 4 maturity stages; reveals inverted response trajectory |
| 4 | `biomarker_discovery.py` | Multi-condition biomarker scoring; 3-gene panel with AUC = 1.0 |
| 5 | `gwas_qtl_crossref.py` | Cross-references biomarkers with 121 published GWAS/QTL loci |
| 6 | `stratified_vs_pooled_comparison.py` | Robustness validation: season-stratified vs pooled DE consistency |
| 7 | `environment_titer_integration.py` | Integrates weather data, viral titer, and DEGs across seasons/stages |
| 8 | `heat_stress_signatures.py` | Heat-stress gene analysis (115 genes, 9 categories); sHSP suppression by GRBV |

### Extended Analysis Scripts

| # | Script | Description |
|---|--------|-------------|
| 9 | `analysis_schematic.py` | Publication-quality pipeline flowchart diagram |
| 10 | `grbv_biology_diagram.py` | Data-driven biological infographic of molecular disruption cascade |
| 11 | `gsea_analysis.py` | Permutation-based GSEA (500 perm) across 8 season × stage conditions |
| 12 | `multifactor_interaction_model.py` | Infection × Season × Stage ANOVA with variance decomposition |
| 13 | `tf_family_analysis.py` | Transcription factor family profiling (WRKY, MYB, NAC, bHLH, ERF/AP2, HSF) |
| 14 | `advanced_ml_biomarkers.py` | Multi-method ML comparison (LR, RF, LASSO, SVM, GBM) + stability selection |
| 15 | `causal_mediation_analysis.py` | Temperature → Titer → Gene expression mediation analysis |
| 16 | `chromosome_deg_hotspots.py` | Genome-wide DEG density and chromosomal hotspot analysis |
| 17 | `coexpression_network.py` | WGCNA co-expression module analysis (Python wrapper for R outputs) |

### Report Generation Scripts

| # | Script | Description |
|---|--------|-------------|
| 18 | `generate_report_pdf.py` | Automated PDF report with embedded figures and discussion |
| 19 | `generate_presentation.py` | Conference PowerPoint deck with all figures |
| 20 | `generate_workflow_pdf.py` | Paper-writing workflow guide (6 phases, journal targets) |

---

## Input Data Requirements

The following files are required in the working directory but are **not included** in this repository due to size:

| File | Description | Size |
|------|-------------|------|
| `counts_matrix.txt` | Raw RNA-Seq count matrix (28,168 genes x 327 samples), tab-separated | ~23 MB |
| `sample_metadata.txt` | Sample metadata with columns: Sample, Infection_Group, Sampling_Season, Maturity_stage, Infection_YSI | ~12 KB |
| `org.Vvinifera.eg.db/` | SQLite annotation database (GO terms, KEGG pathways, gene symbols) | ~12 MB |
| `go-basic.obo` | Gene Ontology OBO file (required by goatools) | ~30 MB |
| `daily_2021.csv` | Daily weather data for 2021 growing season | - |
| `daily_2022.csv` | Daily weather data for 2022 growing season | - |
| `Virus titer.txt` | Viral titer measurements across seasons and stages | ~4 KB |
| `MEs.csv` | WGCNA module eigengenes (from R analysis) | - |
| `Module_genes.xlsx` | WGCNA module gene assignments (14 sheets) | - |

---

## Installation

### Prerequisites

- Python 3.6 or higher
- pip package manager

### Setup

```bash
# Clone the repository
git clone https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis.git
cd GRBV-Transcriptomics-Analysis

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

| Package | Purpose |
|---------|---------|
| numpy | Numerical computing |
| pandas | Data manipulation |
| scipy | Statistical tests (Welch's t-test, Fisher's exact, Spearman) |
| scikit-learn | PCA, ROC-AUC, cross-validation, Random Forest |
| matplotlib | Figure generation |
| seaborn | Statistical visualizations |
| goatools | Gene Ontology enrichment analysis |
| fpdf2 | PDF report generation |
| python-pptx | PowerPoint presentation generation |
| openpyxl | Reading Excel files (WGCNA module gene lists) |

---

## Usage

### Full Pipeline Execution

Run scripts sequentially in the numbered order. Each script reads its required inputs and produces outputs used by downstream scripts.

```bash
# Step 1: Core differential expression (produces plots 01-07, base DE results)
python analysis.py

# Step 2: Pathway enrichment (produces plots 08-14)
python pathway_enrichment.py

# Step 3: Seasonal x maturity analysis (produces plots 15-19)
python seasonal_maturity_progression.py

# Step 4: Biomarker discovery (produces plots 20-22, biomarker CSVs)
python biomarker_discovery.py

# Step 5: GWAS/QTL cross-reference (produces plots 24-25)
python gwas_qtl_crossref.py

# Step 6: Robustness validation (produces plots 26-29)
python stratified_vs_pooled_comparison.py

# Step 7: Environment-titer integration (produces plots 30-33)
python environment_titer_integration.py

# Step 8: Heat-stress analysis (produces plots 34-35)
python heat_stress_signatures.py

# Step 9: Pipeline schematic + biology diagram (produces plots 23, 36)
python analysis_schematic.py
python grbv_biology_diagram.py

# Step 10: GSEA analysis (produces plots 37-39)
python gsea_analysis.py

# Step 11: Multi-factor ANOVA (produces plots 40-42)
python multifactor_interaction_model.py

# Step 12: TF family analysis (produces plots 43-45)
python tf_family_analysis.py

# Step 13: Advanced ML biomarkers (produces plots 46-48)
python advanced_ml_biomarkers.py

# Step 14: Causal mediation (produces plots 49-50)
python causal_mediation_analysis.py

# Step 15: Chromosome DEG hotspots (produces plots 51-52)
python chromosome_deg_hotspots.py

# Step 16: Co-expression network (produces plots 53-55)
python coexpression_network.py

# Step 17: Generate outputs
python generate_report_pdf.py        # PDF report
python generate_presentation.py      # PowerPoint deck
python generate_workflow_pdf.py      # Paper-writing guide
```

### Individual Script Execution

Each analysis script can be run independently if its required input files are present. See the input data requirements above.

---

## Output Files

### Figures (55 PNGs)

| Plots | Source Script | Content |
|-------|--------------|---------|
| 01-07 | `analysis.py` | PCA overview, DE results for Q1-Q6 |
| 08-14 | `pathway_enrichment.py` | Pathway enrichment heatmaps for Q1-Q6, master landscape |
| 15-19 | `seasonal_maturity_progression.py` | Season x stage DEG grids, PCA transitions, pathway detail |
| 20-22 | `biomarker_discovery.py` | Biomarker scores, expression profiles, robustness |
| 23 | `analysis_schematic.py` | Pipeline flowchart |
| 24-25 | `gwas_qtl_crossref.py` | Chromosome ideogram, QTL cross-reference summary |
| 26-29 | `stratified_vs_pooled_comparison.py` | Pooled vs stratified, overlap heatmap, cross-question, dashboard |
| 30-33 | `environment_titer_integration.py` | Weather, titer dynamics, integrated multi-parameter plots |
| 34-35 | `heat_stress_signatures.py` | Heat-stress gene landscape, heat vs virus disentanglement |
| 36 | `grbv_biology_diagram.py` | Data-driven molecular disruption infographic |
| 37-39 | `gsea_analysis.py` | GSEA enrichment scores, Fisher concordance, NES heatmap |
| 40-42 | `multifactor_interaction_model.py` | Variance decomposition, interaction sets, expression profiles |
| 43-45 | `tf_family_analysis.py` | TF family heatmap, enrichment balloon, expression boxplots |
| 46-48 | `advanced_ml_biomarkers.py` | ML AUC comparison, stability selection, consensus panel |
| 49-50 | `causal_mediation_analysis.py` | Mediation diagrams, direct vs indirect effects |
| 51-52 | `chromosome_deg_hotspots.py` | Chromosome ideogram DEG density, Manhattan plot |
| 53-55 | `coexpression_network.py` | Module-trait correlations, eigengene profiles, hub gene heatmaps |

### Data Tables (26+ CSVs)

| File | Content |
|------|---------|
| `GRBV_biomarkers_top50.csv` | Top 50 biomarkers with multi-condition scores and AUC values |
| `GRBV_biomarker_QTL_crossref.csv` | Biomarker-QTL co-localization details |
| `biomarker_genomic_locations.csv` | Genomic coordinates for biomarker genes |
| `heat_stress_DEGs_all_conditions.csv` | Heat-stress gene DE results across all conditions |
| `heat_gene_temperature_correlations.csv` | Spearman correlations of heat genes with temperature |
| `heat_DEGs_recurrence.csv` | Recurrent heat DEGs across conditions |
| `universal_DEGs_across_approaches.csv` | 908 robust consensus DEGs |
| `integrated_environment_titer_DEG.csv` | Temperature-titer-DEG integration table |
| `temperature_at_sampling.csv` | Temperature data at sampling timepoints |
| `Q1_pooled_vs_stratified.csv` | Q1 pooled vs stratified DE comparison |
| `Q2_pooled_vs_stratified.csv` | Q2 pooled vs stratified DE comparison |

### Reports

| File | Content |
|------|---------|
| `GRBV_Transcriptomics_Report.pdf` | Comprehensive report with 25+ figures and discussion |
| `GRBV_Transcriptomics_Presentation.pptx` | 44-slide conference-ready presentation |
| `PAPER_WRITING_WORKFLOW.pdf` | 6-phase paper-writing guide with journal targets |

---

## Key Findings

- **Stage-dependent viral impact**: GRBV effects are strongest at Post-veraison and Harvest, with up to 1,873 DEGs in a single condition
- **Inverted seasonal trajectory**: 2021 DEG counts peak at Harvest; 2022 peaks at Post-veraison, driven by a +7.6 &deg;C temperature differential
- **Heat-stress suppression**: GRBV downregulates 29/36 small heat-shock proteins (sHSPs), reducing the sHSP module score by ~55%
- **Diagnostic biomarkers**: A 3-gene panel achieves AUC = 1.000 across 5-fold cross-validation; 12 universal biomarkers significant in all conditions
- **Genomic validation**: 60% of biomarker genes co-localize within 1 Mb of known GWAS/QTL loci for berry quality traits
- **Pathway disruption**: Consistent downregulation of phenylpropanoid, flavonoid, and stilbenoid biosynthesis across seasons and stages
- **13 co-expression modules**: WGCNA identifies stage-dominated networks; black module most infection-correlated (r=+0.126, p=0.022); 130 hub genes

---

## Citation

If you use this pipeline or its results, please cite:

> Singh, P.P. (2026). Comprehensive transcriptomic analysis of Grapevine Red Blotch Virus effects on *Vitis vinifera* cv. Merlot across seasons and maturity stages. *[Journal TBD]*.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
