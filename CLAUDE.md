# GRBV Transcriptomics Analysis — Project Context

## Study Overview
- **Species**: *Vitis vinifera* cv. Merlot
- **Pathogen**: Grapevine Red Blotch Virus (GRBV), family Geminiviridae
- **Samples**: 327 total (264 RB+ infected, 63 RB- control)
- **Genes**: 28,168 raw → 14,604 after filtering (≥10 counts in ≥10% samples)
- **Seasons**: 2021 (cooler) and 2022 (hotter, +7.6°C at Harvest)
- **Maturity Stages**: Pre-veraison, Veraison, Post-veraison, Harvest
- **Infection Cohorts**: 7 years (2016-2022)
- **DE Thresholds**: Welch's t-test, BH FDR, adj. p-value < 0.05, |log2FC| ≥ 1
- **GitHub**: https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis (private)

## Working Directory
`/Users/prempratapsingh/Data_Analysis/`

## 6 Core Research Questions
- Q1: Effect of infection duration (years since infection)
- Q2: Effect of berry maturity stage
- Q3: Seasonal reproducibility (2021 vs 2022)
- Q4: Environmental modulation of GRBV effects
- Q5: Baseline environmental changes in uninfected vines
- Q6: Infection age trajectory

## Key Findings
1. **Inverted seasonal trajectory**: 2021 DEGs peak at Harvest (1,192); 2022 peaks at Post-veraison (565)
2. **sHSP suppression**: GRBV downregulates 29/36 sHSPs, reducing module score by ~55% (healthy +1.64 → infected +0.74)
3. **HSF-sHSP axis**: r = 0.625, p = 7.2e-37
4. **3-gene diagnostic panel**: AUC = 1.000 across 5-fold CV
5. **12 universal biomarkers**: significant in all conditions
6. **60% QTL co-localization**: biomarkers near known GWAS/QTL loci (1 Mb window)
7. **Sugar QTL link**: Top host biomarker (RNA-dep RNA polymerase 1) maps to sugar accumulation QTL on chr1
8. **908 universal DEGs**: confirmed across pooled and stratified approaches
9. **Cross-season FC correlation**: r = 0.33 (Pre-veraison) to r = 0.69 (Post-veraison)
10. **115 heat-stress genes**: 9 categories, 57 GRBV DEGs, 23 dual-regulated (virus + environment)
11. **Stage dominates variance**: η²=0.426 (ANOVA), Stage > Season×Stage > Infection×Stage > Season > Infection
12. **773 Infection×Season interaction genes**: formally confirms inverted seasonal trajectory
13. **4,146 three-way interaction genes**: Infection × Season × Stage
14. **All 5 ML methods achieve AUC=1.000**: LR, RF, LASSO, SVM-RBF, Gradient Boosting (viral genes N761_gp1-3 dominate)
15. **267 TF genes across 6 families**: ERF/AP2 (79) most disrupted, HSF (13) most consistently affected
16. **90 temperature-responsive genes**: among top 500 variable genes in infected vines
17. **13 WGCNA co-expression modules**: 10,221 genes; Stage dominates module-trait correlations (turquoise r=0.941); black module most infection-correlated (r=+0.126, p=0.022)
18. **130 hub genes identified**: top 10 per module by kME (module eigengene correlation)
19. **GSEA: 270 significant gene set-condition pairs** (FDR<0.05): 375 gene sets × 8 conditions, 500 permutations; down-regulation dominates (77%); KEGG:Ribosome and GO:translation initiation factor activity significant in 6/8 conditions; cross-season NES correlation mirrors inverted trajectory (Pre-vér r=0.09, Post-vér r=0.63)

## Analysis Scripts (run in order)

| # | Script | Purpose | Outputs |
|---|--------|---------|---------|
| 1 | `analysis.py` | Core DE for 6 research questions | Figs 01-07 |
| 2 | `pathway_enrichment.py` | GO (BP/MF/CC) + KEGG enrichment | Figs 08-14 |
| 3 | `seasonal_maturity_progression.py` | Season-stratified DE across 4 stages | Figs 15-19 |
| 4 | `biomarker_discovery.py` | Multi-condition scoring, 3-gene panel | Figs 20-22 |
| 5 | `gwas_qtl_crossref.py` | Cross-ref with 121 GWAS/QTL loci | Figs 24-25 |
| 6 | `stratified_vs_pooled_comparison.py` | Robustness: pooled vs stratified | Figs 26-29 |
| 7 | `environment_titer_integration.py` | Weather + titer + DEG integration | Figs 30-33 |
| 8 | `heat_stress_signatures.py` | 115 heat genes, sHSP suppression | Figs 34-35 |
| 9 | `analysis_schematic.py` | Pipeline flowchart | Fig 23 |
| 10 | `grbv_biology_diagram.py` | Integrated biology infographic | Fig 36 |
| 11 | `gsea_analysis.py` | NumPy-vectorized GSEA (500 perm, HPC) | Figs 37-39 |
| 12 | `multifactor_interaction_model.py` | Infection × Season × Stage ANOVA | Figs 40-42 |
| 13 | `tf_family_analysis.py` | TF family profiling (6 families) | Figs 43-45 |
| 14 | `advanced_ml_biomarkers.py` | Multi-method ML biomarker comparison | Figs 46-48 |
| 15 | `causal_mediation_analysis.py` | Temp → Titer → Gene mediation | Figs 49-50 |
| 16 | `chromosome_deg_hotspots.py` | Genome-wide DEG density & hotspots | Figs 51-52 |
| 17 | `coexpression_network.py` | WGCNA co-expression module analysis | Figs 53-55 |
| 18 | `generate_report_pdf.py` | Comprehensive PDF report | GRBV_Transcriptomics_Report.pdf |
| 19 | `generate_presentation.py` | Conference PPTX deck | GRBV_Transcriptomics_Presentation.pptx |
| 20 | `generate_workflow_pdf.py` | Paper-writing guide PDF | PAPER_WRITING_WORKFLOW.pdf |

## Figures (55 PNGs)

| Fig | File | Content |
|-----|------|---------|
| 01 | `01_global_pca_overview.png` | PCA: Stage > Season > Infection > Duration |
| 02 | `02_Q1_infection_years.png` | Infection duration: biphasic response |
| 03 | `03_Q2_maturity_stages.png` | Stage-dependent DEGs, Post-ver/Harvest strongest |
| 04 | `04_Q3_seasonal_effects.png` | Cross-season FC correlation r=0.33-0.69 |
| 05 | `05_Q4_environmental_effects.png` | Environment modulates (not just adds to) GRBV |
| 06 | `06_Q5_uninfected_environment.png` | Baseline seasonal changes in RB- controls |
| 07 | `07_Q6_transcriptome_trajectory.png` | Biphasic: acute defense → chronic metabolism |
| 08 | `08_Q1_pathway_enrichment.png` | Pathway enrichment for infection duration |
| 09 | `09_Q2_pathway_enrichment.png` | Pathway enrichment for maturity stage |
| 10 | `10_Q3_pathway_enrichment.png` | Pathway enrichment for seasonal effects |
| 11 | `11_Q4_pathway_enrichment.png` | Pathway enrichment for environment |
| 12 | `12_Q5_pathway_enrichment.png` | Pathway enrichment for uninfected baseline |
| 13 | `13_Q6_pathway_enrichment.png` | Pathway enrichment for infection trajectory |
| 14 | `14_master_pathway_landscape.png` | Master heatmap: all conditions, phenylpropanoid top |
| 15 | `15_seasonal_maturity_DEG_comparison.png` | KEY: Inverted trajectory 2021 vs 2022 |
| 16 | `16_seasonal_maturity_pathway_enrichment.png` | Pathway grid across 8 season x stage combos |
| 17 | `17_seasonal_maturity_pathway_landscape.png` | Pathway landscape heatmaps by season |
| 18 | `18_seasonal_maturity_transitions_PCA.png` | PCA stage transitions |
| 19 | `19_seasonal_maturity_pathway_detail.png` | Detailed pathway progression per stage |
| 20 | `20_biomarker_overview.png` | Scores, AUC, 3-gene panel, volcano |
| 21 | `21_biomarker_profiles.png` | Top 16 biomarker expression profiles |
| 22 | `22_biomarker_robustness.png` | FC consistency, AUC stability |
| 23 | `23_analysis_schematic.png` | Pipeline flowchart |
| 24 | `24_biomarker_QTL_chromosome_map.png` | Chromosome ideogram with biomarker positions |
| 25 | `25_biomarker_QTL_crossref_summary.png` | QTL overlap statistics |
| 26 | `26_pooled_vs_stratified_comparison.png` | Pooled vs stratified DE comparison |
| 27 | `27_overlap_connectivity_heatmap.png` | DEG overlap across all questions |
| 28 | `28_Q1Q2_vs_Q3Q4_crossreference.png` | Cross-question integration |
| 29 | `29_comprehensive_summary_dashboard.png` | Integrated findings dashboard |
| 30 | `30_weather_overview_sampling.png` | Temperature/weather at sampling |
| 31 | `31_virus_titer_analysis.png` | Viral load dynamics |
| 32 | `32_integrated_temp_titer_DEG.png` | Multi-parameter integration |
| 33 | `33_temperature_titer_shift.png` | Temperature-titer correlations |
| 34 | `34_heat_stress_gene_landscape.png` | Heat gene categories, DEG counts, module scores |
| 35 | `35_heat_vs_virus_disentanglement.png` | Virus vs heat classification, HSF-sHSP axis |
| 36 | `36_GRBV_biology_diagram.png` | Integrated molecular disruption cascade |
| 37 | `37_gsea_enrichment_plots.png` | GSEA running ES across season × stage |
| 38 | `38_gsea_vs_fisher_concordance.png` | GSEA vs Fisher ORA concordance scatter |
| 39 | `39_gsea_nes_heatmap.png` | GSEA NES heatmap across conditions |
| 40 | `40_variance_decomposition.png` | ANOVA η² variance decomposition |
| 41 | `41_interaction_gene_sets.png` | Interaction term gene counts & overlap |
| 42 | `42_interaction_profiles.png` | Top Infection×Season interaction profiles |
| 43 | `43_tf_family_heatmap.png` | TF family %DE and mean FC across conditions |
| 44 | `44_tf_enrichment_balloon.png` | TF family DEG enrichment balloon plot |
| 45 | `45_tf_expression_boxplots.png` | Top TF gene expression boxplots |
| 46 | `46_ml_method_comparison.png` | Multi-method AUC comparison + ROC |
| 47 | `47_stability_selection.png` | LASSO stability selection (100 subsamples) |
| 48 | `48_consensus_panel.png` | Consensus biomarker panel + ROC |
| 49 | `49_mediation_diagrams.png` | Causal mediation path diagrams |
| 50 | `50_direct_vs_indirect_effects.png` | Direct vs titer-mediated temperature effects |
| 51 | `51_chromosome_deg_ideogram.png` | Chromosome ideogram with DEG density |
| 52 | `52_deg_manhattan_plot.png` | DEG enrichment Manhattan + per-chr counts |
| 53 | `53_module_trait_correlations.png` | WGCNA module-trait correlation heatmap |
| 54 | `54_module_eigengene_profiles.png` | Module eigengene profiles across conditions |
| 55 | `55_hub_gene_heatmaps.png` | Hub gene expression heatmaps per module |

## CSV Data Tables (26 files)

| File | Content |
|------|---------|
| `GRBV_biomarkers_top50.csv` | Top 50 biomarkers with composite scores, AUC, direction |
| `GRBV_biomarker_QTL_crossref.csv` | Biomarker-QTL co-localization (27/45 overlap) |
| `biomarker_genomic_locations.csv` | Chromosomal coordinates for biomarker genes |
| `universal_DEGs_across_approaches.csv` | 908 robust consensus DEGs |
| `Q1_pooled_vs_stratified.csv` | Q1 pooled vs stratified DE comparison |
| `Q2_pooled_vs_stratified.csv` | Q2 pooled vs stratified DE comparison |
| `integrated_environment_titer_DEG.csv` | Temperature-titer-DEG integration table |
| `temperature_at_sampling.csv` | Temperature at sampling timepoints |
| `heat_stress_DEGs_all_conditions.csv` | Heat gene DE results across all conditions |
| `heat_gene_temperature_correlations.csv` | Spearman correlations with temperature |
| `heat_DEGs_recurrence.csv` | Recurrent heat DEGs across conditions |
| `daily_2021.csv` | Daily weather data 2021 season |
| `daily_2022.csv` | Daily weather data 2022 season |
| `gsea_results_all_conditions.csv` | GSEA NES, FDR for all gene sets × conditions |
| `multifactor_anova_results.csv` | Per-gene ANOVA F-stats, p-values, η² (14,604 genes) |
| `interaction_gene_categories.csv` | Genes with significant interaction terms |
| `tf_family_DE_summary.csv` | TF family DE enrichment across conditions |
| `tf_family_gene_details.csv` | Individual TF gene DE results |
| `advanced_ml_comparison.csv` | 5-method ML AUC comparison |
| `stability_selection_results.csv` | LASSO stability selection frequencies |
| `consensus_biomarker_panel.csv` | Cross-method consensus biomarker genes |
| `mediation_results.csv` | Temperature → Titer → Gene mediation results |
| `chromosome_deg_density.csv` | 1 Mb window DEG density per chromosome |
| `deg_hotspot_qtl_overlap.csv` | DEG hotspot × QTL region overlaps |
| `coexpression_modules.csv` | WGCNA module summary (size, hub gene count, trait correlations) |
| `module_trait_correlations.csv` | Module-trait Pearson correlations and p-values |
| `hub_genes_per_module.csv` | Top 10 hub genes per module (by kME) |

## Input Data (not in repo — too large)

| File | Size | Content |
|------|------|---------|
| `counts_matrix.txt` | 23 MB | Raw RNA-Seq counts (28,168 genes x 327 samples) |
| `sample_metadata.txt` | 12 KB | Sample, Infection_Group, Sampling_Season, Maturity_stage, Infection_YSI |
| `org.Vvinifera.eg.db/` | 12 MB | SQLite annotation DB (59,449 GO, 13,801 KEGG) |
| `go-basic.obo` | 30 MB | Gene Ontology OBO file |
| `Virus titer.txt` | 4 KB | Viral titer measurements |
| `MEs.csv` | - | WGCNA module eigengenes (from R analysis) |
| `Module_genes.xlsx` | - | WGCNA module gene assignments (14 sheets) |

## Reports & Outputs

| File | Description |
|------|-------------|
| `GRBV_Transcriptomics_Report.pdf` | 55-page comprehensive report (8 sections, 55 figures) |
| `GRBV_Transcriptomics_Presentation.pptx` | 64-slide conference deck with all figures |
| `PAPER_WRITING_WORKFLOW.pdf` | 6-phase paper-writing guide |
| `GRBV_Transcriptomics_Report.md` | Markdown version of report |
| `PAPER_WRITING_WORKFLOW.md` | Markdown version of workflow |

## Top Host Biomarkers

| Rank | Gene | Name | Score | AUC | Direction |
|------|------|------|-------|-----|-----------|
| 1 | LOC100854143 | F-box/LRR-repeat protein | 0.741 | 0.880 | UP |
| 2 | LOC100256964 | E3 ubiquitin ligase RSL1 | 0.731 | 0.894 | UP |
| 3 | LOC104881857 | Protein NTM1-like 9 | 0.725 | 0.869 | UP |
| 4 | LOC100268011 | RNA-dep. RNA polymerase 1 | 0.716 | 0.895 | UP |
| 5 | LOC100257913 | Glycine-rich cell wall prot. | 0.705 | 0.871 | DOWN |

## Heat-Stress Gene Categories

| Category | Expressed | GRBV DEGs | Key Finding |
|----------|-----------|-----------|-------------|
| sHSP | 36 | 29 DOWN | Strongest suppression |
| HSP70 | 18 | 8 | Moderate effect |
| HSP90 | 7 | 3 | Selective |
| HSP100/ClpB | 5 | 2 | Minimal |
| HSF | 15 | 6 | Regulatory hub |
| GolS | 5 | 3 | Osmotic link |
| LEA | 8 | 4 | Desiccation |
| Antioxidant | 12 | 5 | ROS defense |
| GO:heat | 9 | 3 | General heat |

## Dependencies
numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib, seaborn, goatools, fpdf2, python-pptx, openpyxl
