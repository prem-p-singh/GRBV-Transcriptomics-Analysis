# GRBV Transcriptomics Analysis — Project Context

@Rules.md

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

---

# SESSION STATE AS OF 2026-04-20

> Everything below captures work done on the extended FARM pipeline, the sensitivity analyses that corrected several claims, and the current story for the paper. Read this first in any new session.

## The corrected narrative (don't forget)

The headline "GRBV disables the plant's thermostat" framing from an earlier draft is **wrong**. The actual story that fits both the transcriptomic and raw viral-titer data is:

**GRBV causes disease by three simultaneous mechanisms:**
1. **Ignores canonical plant immunity** — all 10 defense categories (NLR, PRR, PR, SA, JA, ET, WRKY, callose, ROS, R-adapters) show |log₂FC|<0.08 in infected vs healthy. The plant is NOT fighting back with its usual weapons.
2. **Dampens but does not destroy the HSF → sHSP thermal-defence axis** — 29/36 sHSPs suppressed, HSF-sHSP r=0.625, but all 6 WGCNA modules preserve with Zsum 5.8-36.4 (network structure intact).
3. **Rewires sugar metabolism via DOF TFs** — motif enrichment (rule 05) finds DOF3.6/DOF1.7 are top enriched in brown and yellow modules; viral N761_gp1/3 transcripts are themselves top network rewiring hubs.

**Temperature acts biphasically on viral titer:**

| Stage | 2021 cool titer | 2022 hot titer | Hot/Cool |
|---|---|---|---|
| Pre-véraison | 457 k | 1,630 k | **3.6×** (heat promotes replication) |
| Véraison | 832 k | 1,387 k | 1.7× |
| Post-véraison | 755 k | 729 k | 1.0× |
| Harvest | 2,211 k | 393 k | **0.18×** (heat suppresses replication) |

Net: hot seasons front-load disease; cool seasons back-load it. This is the "inverted seasonal trajectory" — now mechanistically interpretable rather than descriptive.

## FARM HPC extended pipeline — 13 rules, ALL COMPLETED

Working directory on FARM: `~/grbv_extended/`
Local mirror: `/Users/prempratapsingh/Data_Analysis/farm_pipeline/`
SSH alias: `farm` (see `~/FARM_HPC.md`)
Conda env: `grbv-extended` (Python tools) + R/4.4.2 LMOD + `~/grbv_extended/Rlibs` (Bioconductor packages installed via BiocManager; DGCA from CRAN archive)
SLURM: partition `low`, account `publicgrp`, QOS `publicgrp-low-qos`, `--requeue` on every job

| # | Rule | Script | Status | Key output |
|---|------|--------|--------|------------|
| 00 | weather_aggregate | `00_weather_aggregate.py` | ✅ | weather_per_sample.csv (per-sample 14-day window aggregates) |
| 01 | grn_inference | `01_grn_genie3.R` | ✅ | 100,000 directed TF→target edges; top hub TFs |
| 02 | module_preservation | `02_module_preservation.R` | ✅ | All 6 modules Zsum 5.8-36.4 (strongly preserved in RB+) |
| 03 | subtype_discovery | `03_subtype_nmf.py` | ✅ | NMF k=2; S1 = Pre+Vér, S2 = Post+Harvest (maturity, not severity) |
| 04 | diffcoex | `04_diffcoex.R` | ✅ | 50k rewired pairs; N761_gp1/3 in top-5 rewiring hubs |
| 05 | motif_enrichment | `05a_extract_promoters.py` + `05b_motif_enrichment.py` | ✅ | 138 enriched motif×module pairs; AHL family + DOF3.6/1.7 |
| 06 | pseudotime | `06_pseudotime.R` | ✅ | RD22-C top temporally variable gene (WaldStat 1,877) |
| 07 | mofa_integration | `07_mofa.py` | ✅ | 10 factors; expression 85%, weather 79%, titer 47% VE |
| 08 | hima_mediation | `08_hima.R` | ✅ | Multi-mediator (see sensitivity section) |
| 09 | defense_repertoire | `09_defense_repertoire.py` | ✅ | 10 categories all flat |
| 10 | subtype_characterize | `10_subtype_characterize.py` | ✅ | 1,000 top DE markers per subtype |
| 11 | harmonize_literature_degs | `11_harmonize_literature_degs.py` | ⏸ scaffolded, awaiting literature tables | NCBI-crosswalk harmonization |
| 12 | comparative_meta_analysis | `12_comparative_metaanalysis.py` | ⏸ scaffolded | Fisher overlaps + UpSet + specificity classes |

Helpers: `setup_farm.sh`, `run_farm.sh`, `push.sh`, `pull.sh`, `aggregate_sensitivity.py`, `build_report_pdf.py`, `build_integrated_pdf.py`, `build_sensitivity_pdf.py`.

Config: `farm_pipeline/config/config.yaml` exposes all thresholds (incl. sampling dates, weather window, HIMA top_var, NMF k range, etc.).

## Sensitivity analysis — critical correction

Every analysis using a `top-N variable genes` filter was re-run at 3–5 additional thresholds. See `farm_pipeline/docs/SENSITIVITY_REPORT.md` + `.pdf`.

**Biggest finding: single-mediator claim RETIRED.** At the original tv=3,000, HIMA identified `LOC132254135`. This does NOT replicate at other thresholds. The robust set is:

| Gene | HIMA thresholds passed | Main-paper mediation passed | Direction at heat |
|---|---|---|---|
| `LOC100266309` | 3 of 5 (5k, 7.5k, 10k) | — | pro-viral (α<0, β<0) |
| `LOC100232974` | 2 of 5 (7.5k, 10k) | — | antiviral (α<0, β>0) |
| `LOC100249249` | 2 of 5 (7.5k, 10k) | **YES — 2k (rank 1, p=5e-11) and 5k (rank 1, p=1e-10)** | pro-viral (α>0, β>0) |
| `LOC132254135` | 1 of 5 (3k only) | — | RETIRED |

**`LOC100249249` is the cross-pipeline robust mediator** — the only gene identified by both HIMA (Temp→Gene→Titer) and the main paper's independent Temp→Titer→Gene analysis. This is the primary qPCR / VIGS validation target.

Other sensitivity findings:
- **ML biomarker AUC=1.000**: holds ONLY when viral N761_gp* transcripts are in the feature pool. Host-only AUC drops to 0.92–1.00 depending on method. The paper's "5 methods achieve AUC=1.000" claim is accurate but largely reflects "detect virus by detecting virus." Must state this caveat.
- **Main paper's causal_mediation_analysis.py at tv=500 (original) finds ZERO significant mediators.** Signal emerges at tv=2k-5k with `LOC100249249` as top hit.
- **Subtype NMF**: biologically robust (always maturity-driven, no severity subtype) but k changes (2 vs 3) with top_var; tv=1k/3k/5k all give k=3 with ARI≈0.98. tv=2000 gives k=2 as the outlier.
- **Module preservation, GRN hubs, DiffCoEx viral-hub claim**: robust across thresholds.
- **Pseudotime**: only LOC100232901 in top-30 at every tv; RD22-C is in top-5 at some but not all thresholds.

## Viral titer data — the reality check

Raw `Virus titer.txt` shows titer by Season × Year × Biorep × Stage. Do NOT rely on mean summaries that average out the biphasic pattern.

Mean titer (across all biological replicates and infection cohorts):

```
Season  Pre-véraison    Véraison    Post-véraison  Harvest
2021       457,000       832,000        755,000    2,211,000
2022     1,630,000     1,387,000        729,000      393,000
```

`LOC132254135`-style single-gene causal claims don't survive this data. The biphasic pattern is the actual finding, and multiple host genes with opposing heat effects (LOC100266309 pro-viral, LOC100232974 antiviral, LOC100249249 pro-viral) mechanistically support it.

## Corrected top-level numbers table

| Metric | Value | Source |
|---|---|---|
| Samples | 327 (264 RB+, 63 RB−) | main |
| Genes tested | 14,604 | main |
| Universal DEGs | 908 | main |
| sHSPs suppressed | 29 of 36 | main |
| HSF–sHSP correlation | r = 0.625, p = 7.2e-37 | main |
| Defense categories with |log₂FC|>0.5 | 0 of 10 | rule 09 |
| WGCNA modules preserved in infected | 6 of 6 (Zsum 5.8–36.4) | rule 02 |
| Infection-severity subtypes | 0 (NMF = maturity) | rule 03 + 10 |
| Stage variance component | η² = 0.426 | main |
| Pre-véraison titer hot/cool | 3.6× | raw titer data |
| Harvest titer hot/cool | 0.18× | raw titer data |
| Robust temp→titer mediators | 3 (LOC100266309, LOC100232974, **LOC100249249**) | rules 08 + sensitivity + main-paper mediation sweep |
| Cross-pipeline robust mediator | **LOC100249249** | HIMA ∩ main-paper mediation |
| Directed TF→target edges | 100,000 | rule 01 |
| DOF motifs enriched in sugar modules | DOF3.6, DOF1.7 | rule 05 |
| Top pseudotime gene (tv=3k) | RD22-C | rule 06 |
| Pipeline rules completed on FARM | 13 of 13 | full pipeline |

## Documents in `farm_pipeline/docs/` (all committed to GitHub)

| File | Purpose |
|---|---|
| `INTEGRATED_ANALYSIS.md` + `.pdf` | Rigorous everything-connected account with numbers and claims ledger |
| `EXTENDED_RESULTS_REPORT.md` + `.pdf` | Per-rule technical report (13 pages, includes embedded figures) |
| `PLAIN_LANGUAGE_SUMMARY.md` | Non-specialist version of the story |
| `SENSITIVITY_REPORT.md` + `.pdf` | Full sensitivity sweep and claim corrections |
| `HSFA2_VIGS_EXPERIMENT_PLAN.md` | Wet-lab protocol (7 months; LOC100249249 primary + LOC100266309/232974 secondary in qPCR panel) |
| `LITERATURE_SCHEMA.md` | Schema for Tier 1 meta-analysis (rules 11–12) |

## Paper submission plan

- **Target**: Nature Plants. My earlier statement that Nature Plants requires the VIGS experiment was wrong — it does NOT strictly require functional validation. Observational papers with strong causal inference and novel conceptual framing are accepted.
- **Path A (recommended)**: submit now with (i) biphasic titer + mechanistic framing, (ii) LOC100249249 causal convergence as headline, (iii) defense-evasion + sHSP axis + DOF rewiring as three-thread story, (iv) VIGS protocol proposed in discussion as next work.
- **Path B (if targeting PNAS or strong Nature Plants bid)**: execute VIGS wet-lab first (7 months per `HSFA2_VIGS_EXPERIMENT_PLAN.md`).
- **Comparative meta-analysis (rules 11–12)** adds polish but is not strictly required. Tier 1 (download supplementary DEG tables from 5 published Vitis studies, ~30–60 min of manual work) answers ~80% of the "is GRBV-specific?" question.

## Immediate open items

1. **Main-paper scripts still use hardcoded top-N filters** that sensitivity analysis flagged:
   - `causal_mediation_analysis.py` line 253: `nlargest(500)` (should be ≥2,000 for robust signal)
   - `advanced_ml_biomarkers.py` line 88: `nlargest(200)` (AUC caveat should be added explicitly to narrative)
   Decide whether to (a) patch defaults + re-run and update main-paper figures, or (b) add a sensitivity-check paragraph to the paper's results section.
2. **Annotate the three robust mediators** via NCBI Entrez / InterPro to get gene identities (`LOC100266309`, `LOC100232974`, `LOC100249249`). Current names are just locus IDs.
3. **Populate `farm_pipeline/data/literature/`** with 5 supplementary DEG tables from published Vitis stress studies (GLRaV-3, GPGV, powdery mildew, downy mildew, heat) to run Tier 1 meta-analysis.
4. **Regenerate main-paper Figures 49–50** if causal_mediation is re-run at tv≥2k (will show many more significant mediators than the current zero at tv=500).
5. **Reframe paper's "AUC=1.000" messaging** to distinguish all-vars classification (AUC=1.000 robust) from host-only classification (AUC 0.92–1.00 depending on method).

## Git repository

- Repo: https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis
- Main branch: `master`. Most recent commits cover the extended pipeline + sensitivity correction wave.
- `.gitignore`: results/, resources/, logs/, data/ all excluded; `docs/*.md` and `docs/*.pdf` explicitly allowed.
- Commits never include Claude co-author lines (per user preference).

## FARM state as of this checkpoint

- `~/grbv_extended/results/` (18 MB) — all 23 output files, fully mirrored locally
- `~/grbv_extended/results/sensitivity/` — raw sensitivity sweeps (gitignored; local mirror at `farm_pipeline/results/sensitivity/`)
- `~/grbv_extended/resources/` (550 MB) — genome + motifs + TF list + id_crosswalk.tsv (keep — rebuildable but saves 15–30 min)
- `~/grbv_extended/Rlibs/` (80 MB) — installed R packages (keep)
- Quota usage: 8.7 GB / 20 GB (44%)
- Runtime artifacts cleaned up in last session (test files, watcher logs, snakemake cache removed)

---

# SESSION STATE UPDATE 2026-04-20 (late): CRITICAL_REVIEW fix pass + real QC context

## New inputs that became available this session

User provided pointers to three files that plug the major "unknown" gaps:

1. **MultiQC report**: `~/Desktop/Oberholster lab/RNA Seq/Bioinformatics_Analysis/QC/2025_April_Singh_Grape_TAGseq_multiqc_report.html` — Bradley N. Jenner (UC Davis), May 2025.
2. **Alignment summary**: `~/Desktop/Oberholster lab/RNA Seq/Bioinformatics_Analysis/QC/singh_alignment_summary.html` + `.Rmd` — same author.
3. **Enhanced metadata**: `~/Data_Analysis/Correct.Sample_Metadata.xlsx` — contains `Vine_ID` and `Bio-replicate_ID` for every sample. **44 unique vines × up to 8 samples each (4 stages × 2 seasons).**

## Confirmed wet-lab / bioinformatics pipeline (was [TO CONFIRM] earlier)

- **Library**: TAG-seq (3′-end tag counting, UMI-tagged). ~1 M uniquely-mapped reads per library typical.
- **Pre-alignment QC + trimming**: FastQC + HTStream pipeline (AdapterTrimmer, CutTrim, NTrimmer, LengthFilter, PolyATTrim, QWindowTrim, SeqScreener for rRNA/chloroplast screening, stats).
- **Deduplication**: umi_tools dedup.
- **Alignment**: STAR. Two references tried:
  - VITVvi_vMerlFPS15_v1.1 (Merlot) — rejected (~50 % unmapped/multi, highly fragmented assembly).
  - **ASM3070453v1 (Pinot Noir, NCBI GCF_030704535.1)** — adopted. Minimum per-library unique reads ≈ 750 k.
- **Quantification**: featureCounts (Subread) → `2025_April_Singh_Grape_TAGseq_ASM3070453v1_counts.txt`. This is what `counts_matrix.txt` in the repo derives from.
- **Analyst**: Bradley N. Jenner, Oberholster lab.
- **Data still to be deposited**: raw FASTQ → NCBI SRA, counts+metadata → GEO.

## CRITICAL_REVIEW fix-pass outcome (commits fa107b4..686565e..4a83273..this session)

**All 13 items from `CRITICAL_REVIEW.md` addressed at the analysis/documentation level.** Only wet-lab qPCR and data deposition remain as user-dependent actions.

### Fix #1 — DESeq2 + limma-voom with covariate adjustment (DONE)

Script: `scripts/de_deseq2_voom.R`. Design `~ Maturity_stage + Sampling_Season + Infection_Group`. FARM job 32497371, 15 min.

| Comparison | DESeq2 DEGs | voom DEGs | Shared | Welch-t (paper) |
|---|---|---|---|---|
| Q6 pooled RB+ vs RB− | 126 | 150 | 106 (Jaccard 0.62) | ~908 universal |
| Q2 Pre-véraison | 111 | 105 | — | — |
| Q2 Véraison | 130 | 124 | — | — |
| Q2 Post-véraison | 541 | 532 | — | — |
| Q2 Harvest | 1,025 | 964 | — | — |
| Q3 season (RB+) | 121 | 153 | — | — |
| Q5 season (RB−) | 182 | 228 | — | — |
| Q4 Season × Infection (LRT) | 52 | — | — | 773 paper |

**Headline: the paper's "908 universal DEGs" is not reproducible under covariate-adjusted DESeq2/voom.** Only 13.6 % (123/907) of the paper's universal DEGs survive Q6 pooled at the same thresholds. The gap is driven by stage/season confounding that the Welch-t did not remove. Per-stage stratification recovers large signals (1,025 at Harvest). The correct headline is: "~130 covariate-adjusted core infection DEGs; stage-specific subsets up to 1,000 at Harvest."

### Fix #1b — Repeated-measures limma-voom (DONE)

Script: `scripts/de_repeated_measures_voom.R`. Uses vine_id from `sample_metadata_v2.txt`. Two-iteration `duplicateCorrelation` blocking on vine_id. FARM job 32497375, 15 min.

| Comparison | Within-vine ρ̂ | RM DEGs | non-RM voom DEGs | Δ |
|---|---|---|---|---|
| Q6 pooled | −0.011 | 150 | 150 | 0 |
| Q2 Pre-véraison | +0.010 | 105 | 105 | 0 |
| Q2 Véraison | −0.017 | 124 | 124 | 0 |
| Q2 Post-véraison | +0.041 | 529 | 532 | −3 |
| Q2 Harvest | −0.008 | 964 | 964 | 0 |

**Within-vine correlation is negligible across every comparison (|ρ̂| ≤ 0.04).** The stage × season × infection programme dominates the transcriptome so thoroughly that vine-specific baseline is a negligible second-order effect. Pseudoreplication is **not** a practical concern for this dataset. Report de_deseq2_voom as primary, RM as supplementary sensitivity check.

### Fix #3 — Main-paper causal mediation rerun at tv=2000 (DONE)

Script: `scripts/mediation_tv2000.py`. At the paper's default tv=500, **zero** mediators pass joint significance. At tv=2000, **39 mediators pass** (LOC100249249 top, p=4.8 × 10⁻¹¹). Overwrites `~/Data_Analysis/mediation_results.csv`. Paper's Figs 49–50 need regeneration.

### Fix #4 + #13 — ML biomarker reframe + Leave-one-season-out CV (DONE)

Script: `scripts/biomarker_loso_cv.py`.

| Panel | Random 5-fold AUC | **Leave-one-season-out AUC** |
|---|---|---|
| All vars (viral + host) | 1.000 (all 5 methods) | 0.98–1.00 |
| Host-only | 0.96–0.998 | **0.73–0.88** (LR 0.88, RF 0.86, LASSO 0.86, SVM 0.79, GB 0.73) |

**The paper's "AUC=1.000 across 5 methods" holds only when viral N761_gp* transcripts are in the feature set.** Real host-only generalisation (LOSO) is 0.73–0.88, substantially less than the published claim.

### Fix #10 — QTL permutation null (DONE)

Script: `scripts/qtl_permutation_null.py`. Observed 27/27 biomarkers within 1 Mb of a QTL (100 %). Null: 10,000 chromosome-stratified permutations from 27,527 positioned genes → 13.4 ± 2.5 overlaps. **2.01 × fold enrichment, empirical one-sided p < 0.0002.** Co-localisation claim is defensible.

### Fix #11 — NCBI gene annotations (DONE)

Script: `scripts/ncbi_annotations.py`. 82 of 86 key LOC IDs resolved. Headline biological identities:

| LOC ID | NCBI annotation | Significance |
|---|---|---|
| **LOC100249249** | **Multiprotein-bridging factor 1c (MBF1c)** | Cross-pipeline robust mediator; known plant HSF co-activator |
| **LOC100232974** | **Sucrose-phosphate synthase 1 (SPS1)** | Robust mediator; direct molecular link to sugar-accumulation phenotype |
| LOC100266309 | Uncharacterised | Robust mediator; functional probe for VIGS |
| LOC100241332 | ATHB-13 (HD-Zip) | Top GRN hub TF |
| LOC100268011 | RNA-dep. RNA polymerase 1 | Top biomarker on chr1 sugar QTL |
| LOC100854143 | F-box/LRR-repeat At3g48880 | Top biomarker (UP) |
| LOC100256964 | E3 ubiquitin ligase RSL1 | Top biomarker |
| LOC104881857 | Protein NTM1-like 9 | Top biomarker |
| LOC100257913 | Glycine-rich cell wall structural protein 2 | Top biomarker (DOWN) |
| LOC100232901 | β-galactosidase BG1-like | Robust pseudotime top gene |

**The mediator panel is biologically coherent with the paper's two headline threads**: MBF1c sits inside the HSF circuit; SPS1 produces sucrose. This replaces the anonymous "3-gene panel" framing with a mechanistically interpretable one.

### Fix #9 — METHODS_ADDENDUM + DE_REANALYSIS_REPORT (DONE)

- `docs/METHODS_ADDENDUM.md` — 14-section reviewer-facing draft, now populated with the actual TAGseq / HTStream / STAR-ASM3070453v1 / featureCounts pipeline from Brad Jenner's reports. Only [TO CONFIRM] items left are sequencing-platform model, RIN, core facility name, SRA/GEO accessions.
- `docs/DE_REANALYSIS_REPORT.md` — full comparison of Welch-t vs DESeq2 vs voom vs repeated-measures-voom.

### Fix #12 — causal → statistical language softening (DONE)

Across INTEGRATED_ANALYSIS, EXTENDED_RESULTS_REPORT, PLAIN_LANGUAGE, HSFA2_VIGS_EXPERIMENT_PLAN. "Causal" reserved for post-VIGS language. CRITICAL_REVIEW itself retains "causal" in historical context (intentional).

## New documents created this session

| Document | Purpose |
|---|---|
| `farm_pipeline/docs/CRITICAL_REVIEW.md` + `.pdf` | 15-section journal-reviewer-style critique (origin of the fix pass) |
| `farm_pipeline/docs/METHODS_ADDENDUM.md` | 14-section Methods draft with TAGseq pipeline now documented |
| `farm_pipeline/docs/DE_REANALYSIS_REPORT.md` | DESeq2/voom/RM-voom rerun results + concordance with Welch-t |
| `sample_metadata_v2.txt` | Enhanced metadata with vine_id + biorep_id (from Correct.Sample_Metadata.xlsx) |
| `farm_pipeline/scripts/fix0{1,1b,3,4}*.{R,py}` | Fix-pass scripts |
| `farm_pipeline/scripts/fix{10,11}*.py` | QTL permutation + NCBI annotation scripts |
| `farm_pipeline/results/fix0{1,1b,3}_*` | Fix-pass output CSVs (gitignored; local only on Mac) |

## What remains user-dependent (can't fix from Claude side)

1. **qPCR validation** of ≥10 key genes (MBF1c, SPS1, LOC100266309, HSFA2, 3 sHSPs, 3 top biomarkers) — 6–8 weeks wet-lab.
2. **Sequencing platform / core facility / RIN values** in Methods — minutes to pull from lab records.
3. **GEO / SRA deposition** — 2–4 weeks submission time.
4. ~~**Regenerate main-paper Figs 49–50**~~ DONE 2026-04-21. Script `farm_pipeline/scripts/regen_figs49_50.py` reads the mediation tv=2000 mediation_results.csv + ncbi_annotations NCBI annotations and overwrites `49_mediation_diagrams.png` and `50_direct_vs_indirect_effects.png`. Top-6 panel order: MBF1c → phytosulfokine 6 → F-box PP2-B15 → DREB1E → uncharacterised LOC100854026 → glutaredoxin-C9. 39 / 2,000 joint-significant mediators displayed in Fig 50.
5. **Tier-1 literature meta-analysis** — populate `farm_pipeline/data/literature/` with 5 supp tables, rule 11–12 already scaffolded.
6. **VIGS experiment** per `HSFA2_VIGS_EXPERIMENT_PLAN.md` (targets now: MBF1c, SPS1, LOC100266309).

## Submission path from this checkpoint

- **Nature Plants ready after**: items 2, 3 above + one-paragraph rewrites to docs to reflect DESeq2-based DE counts (item 4 now done).
- **PNAS competitive after**: above + VIGS experiment (item 6).

## Headline commits in master (in order)

1. `fa107b4..686565e` — fix pass part 1 (fix #3, #4, #10, #11, #12 partial, #9 docs)
2. `686565e..4a83273` — fix pass part 2 (fix #1 DESeq2 + voom results + DE_REANALYSIS_REPORT)
3. (pending end-of-session) — fix #1b repeated-measures + METHODS_ADDENDUM with real QC info + final PDFs

## Files new session should read first

Order, shortest-to-read first:

1. This block (CLAUDE.md bottom)
2. `farm_pipeline/docs/INTEGRATED_ANALYSIS.md` — canonical story
3. `farm_pipeline/docs/SENSITIVITY_REPORT.md` — threshold-sweep source of truth
4. `farm_pipeline/docs/CRITICAL_REVIEW.md` — the review-style critique driving all fixes
5. `farm_pipeline/docs/METHODS_ADDENDUM.md` — Methods draft for submission
6. `farm_pipeline/docs/DE_REANALYSIS_REPORT.md` — the DESeq2 rerun numbers
7. `farm_pipeline/docs/HSFA2_VIGS_EXPERIMENT_PLAN.md` — wet-lab protocol

## Writing-style reminders (from global CLAUDE.md)

- No em dashes (use commas / colons / rewrite)
- No rule-of-three, no buzzwords (leverage, synergy, passionate, cutting-edge, robust, actionable)
- No filler phrases ("in order to", "it is worth noting", "I am confident that")
- No -ing padding
- Paper-verified claims only — avoid ELISA, TEM, biocontrol, greenhouse/field trial design references unless genuinely performed
- No Claude co-author lines in commits

---

# SESSION STATE UPDATE 2026-04-23: Harvest biochemistry integration + SPS1 reframe

## New input

`Biochemical_Data.xlsx` (2 sheets, 2021 and 2022, 33 rows). Harvest-only juice chemistry: Brix, pH, TA, Malic, NOPA, NH3, YAN + cluster weight + juice volume. Sample IDs map 1-to-1 to transcriptomics bioreps within each (Sampling_Season, Infection_Year, Biorep=1..3) cell.

## New code + data

- `biochemistry_harvest_prep.py` (top level): Python prep + stats only. Produces five CSVs at repo root:
  `biochemistry_harvest.csv`, `biochemistry_transcriptome_harvest_merged.csv`,
  `biochemistry_validation_correlations.csv`, `biochemistry_cohort_trajectory.csv`,
  `biochemistry_season_infection_anova.csv`.
- `figures/scripts_R/supp_biochem.R`: reads the five CSVs, produces S61, S62, S63 (PNG + SVG) in `figures/supp/`.
- `figures/scripts_R/build_full_atlas.py`: extended with entries for S61-S63. `GRBV_Full_Figure_Atlas.pdf` now 73 pages (10 main + 63 supp).
- `docs/BIOCHEM_HARVEST_VALIDATION.md`: results report (full numbers + interpretations).

## Headline findings

1. **Season × Infection interaction on TA is the cleanest phenotypic mirror of the inverted seasonal trajectory.** F=104, p=4.3e-11, η²=0.77. No main effect of season (p=0.69) or infection (p=0.16). GRBV raises TA by +0.62 g/L in 2021 (cool) and lowers it by −0.79 g/L in 2022 (hot). Same inversion for Malic (+398 mg/L vs −195 mg/L). **Fig. S63.**
2. **Brix suppression is 2.6× stronger in 2022** (−9.83) than in 2021 (−3.77). Interaction η²=0.13.
3. **2021 shows a biphasic cohort-age recovery at harvest**: Brix rises +0.10/yr (r=0.76, p=0.0045) across RB(+) infection cohorts while NOPA (−2.45/yr, p=0.006) and NH3 (−2.40/yr, p=0.035) decline. Phenotypic reading: vines partially compensate on sugar over years at the cost of durable nitrogen-metabolism depletion. Direct biochemical analogue of the transcriptomic Q1 biphasic finding. **Fig. S62.**
4. **LOC100266309 and SPS1 (LOC100232974) have defensible within-season negative correlations with Brix** (partial rho −0.45 and −0.42 after controlling for season; 2021 rho −0.65 and −0.73, p < 0.01 for both). **Fig. S61.**
5. **MBF1c (LOC100249249) has no direct biochemical correlate at harvest** (partial rho Brix = +0.04). HSF-axis dampening does not translate to measurable juice chemistry at this stage.

## SPS1 reframe (applied this session)

Higher SPS1 correlates with **lower** Brix within-season (2021 ρ=−0.73, p=0.002). This rules out the naive "SPS1 makes sugar" reading. SPS1 is interpreted as a compensatory upregulation in metabolically-stressed infected berries, with phloem loading / transport as the plausible rate-limiting step. The mediator status (multi-threshold HIMA + main-paper mediation sweep) stands unchanged.

**Edits applied:**

- `figures/scripts_R/main10_biology_cascade.R` (line 42): sugar-pathway cascade box — replaced "SPS1 is a robust temperature mediator (Main 5)" with biochem finding "Brix suppressed 4-9.8 °Bx in infected berries, with TA / Malic inverted between 2021 cool and 2022 hot seasons (S63)". `Main10_biology_cascade.png + .svg` regenerated.
- `docs/METHODS_ADDENDUM.md` (ncbi_annotations gene table): "Robust mediator; direct molecular link to sugar-accumulation phenotype" → "Robust temperature-titer mediator. S61 shows SPS1 correlates negatively with Brix within-season (2021 ρ=−0.73), i.e., up-regulated in low-Brix infected berries. Interpret as compensatory response, not rate-limiting step."
- `docs/MANUSCRIPT_DRAFT_v1.md`: Discussion paragraph on "Temperature, titer and MBF1c" reworked to include Fig. S61 finding and the phloem-loading alternative; added new Results subsection "Harvest juice chemistry tracks the transcriptomic inverted trajectory"; added S61 / S62 / S63 captions.

**Not edited (still accurate):**

- `figures/scripts_R/main05_mediation.R`: SPS1 appears in the top-10 mediators panel C. The mediator-by-joint-significance claim is statistical and stands. No figure change needed.
- `docs/EXTENDED_RESULTS_REPORT.md` and `docs/SENSITIVITY_REPORT.md`: treat SPS1 as mediator without causal sugar claim; already defensible.
- `docs/HSFA2_VIGS_EXPERIMENT_PLAN.md`: lists SPS1 among secondary qPCR readouts. Mediator claim intact.

## File locations (final)

| Artifact | Path |
|---|---|
| Biochem source | `Biochemical_Data.xlsx` |
| Prep script (Python) | `biochemistry_harvest_prep.py` |
| Figure script (R) | `figures/scripts_R/supp_biochem.R` |
| Tidy + stats CSVs | `biochemistry_*.csv` (5 files, repo root) |
| Figures | `figures/supp/S61_biochem_validation.{png,svg}`, `S62_biochem_cohort_trajectory.{png,svg}`, `S63_biochem_season_infection.{png,svg}` |
| Full atlas (73 pages) | `figures/GRBV_Full_Figure_Atlas.pdf` |
| Report | `docs/BIOCHEM_HARVEST_VALIDATION.md` |

## Pending

- Consider adding harvest biochemistry as a 4th factor in the MOFA integration (rule 07) on FARM for a formal multi-omics alignment.
- qPCR panel for VIGS follow-up: LOC100266309 now confirmed as strongest host validation, alongside MBF1c (statistical primary) and SPS1 (as compensatory marker).

---

# SESSION STATE UPDATE 2026-04-23 (late): Tier 1 literature meta-analysis + directory convention

## Where paper figures live vs where diagnostic plots live

A question that came up this session and is worth pinning for future work:

- **Paper figures** (Main + Supp) live in `figures/main/` and `figures/supp/`, produced by R scripts in `figures/scripts_R/` using `theme_grbv.R`, PNG + SVG at 600 DPI. These are the canonical outputs that land in the manuscript and the atlas PDF.
- **Diagnostic / intermediate plots** live under `farm_pipeline/results/<rule>/*.png`. They are raw matplotlib / R output from the FARM pipeline rules (02 module preservation, 03 subtype NMF, 06 pseudotime, 07 MOFA, 09 defense, 10 subtype characterisation, 12 metaanalysis). They are gitignored, not in the atlas, and not styled for publication. Their role is reproducibility: each rule writes whatever plot proves it ran. The underlying findings are incorporated into the paper via the proper R figure pipeline, but the raw diagnostic PNGs are kept on disk as artefacts.
- **Do not confuse the two**: when a new analysis needs a paper figure, the workflow is (i) FARM rule or top-level Python writes CSVs + intermediate diagnostic PNGs to `farm_pipeline/results/` or to the repo root, (ii) an R script in `figures/scripts_R/` reads the CSVs and produces the publication PNG + SVG in `figures/supp/` or `figures/main/`, (iii) `figures/scripts_R/build_full_atlas.py` stitches the atlas. Never add a diagnostic PNG directly to the atlas; port it through R first.

## Tier 1 literature meta-analysis (added this session)

- **Ingest**: `scripts/ingest_tier1.py` reads 6 supplementary files in `farm_pipeline/data/literature/` and emits 18 tidy per-contrast CSVs + YAML sidecars. Primary contrasts kept at top level (1 per paper); secondary contrasts moved to `farm_pipeline/data/literature/secondary/` so Rule 11/12 see a clean 6-paper panel.
- **Crosswalk**: `scripts/build_id_crosswalk.sh` builds `farm_pipeline/resources/id_crosswalk.tsv` (24,817 rows) by merging NCBI `All_Plants.gene_info.gz` (current Vvi + LOC + Symbol) with Ensembl Plants release-50 `Vitis_vinifera.12X.50.entrez.tsv.gz` (legacy VIT_*s*g* → LOC). Unlocks 60–91 % mapping for the 3 VIT-based papers (Blanco-Ulate, Song, Agudelo-Romero).
- **Rule 12** (`farm_pipeline/scripts/12_comparative_metaanalysis.py`): new `--publish_dir` flag writes six CSVs with `tier1_` prefix to the repo root (`tier1_grbv_signature.csv`, `tier1_overlap_matrix.csv`, `tier1_specificity_classes.csv`, `tier1_specificity_summary.csv`, `tier1_fc_correlation_matrix.csv`, `tier1_coverage_report.csv`). Matplotlib plotting is skipped in publish mode; R is the canonical renderer.
- **Figures**: `figures/scripts_R/supp_tier1.R` reads the `tier1_*.csv` files and produces S64 (3-panel overlap summary: signature sizes, GRBV intersections, Fisher enrichment) and S65 (log2FC Spearman concordance heatmap) in `figures/supp/`. Atlas now 75 pages.

## Verified 6-paper Tier 1 panel (citations all PubMed-confirmed)

| Stress | Study | PMID | Coverage |
|---|---|---|---|
| GRBV | Blanco-Ulate et al. 2017 *J Exp Bot* 68:1225 (Zinfandel berries, 2 vineyards × 4 stages) | 28338755 | 68.7 % |
| GLRaV-3 | Song et al. 2022 *Viruses* 14:1831 (Cabernet Franc leaves + berries, 3 stages) | 36016453 | 68–91 % |
| Powdery mildew | Fung et al. 2008 *Plant Physiol* 146:236 (Cabernet Sauvignon leaves) | 17993546 | 0.5 % (Affy probeset → VIT bridge needed) |
| Downy mildew | Polesani et al. 2010 *BMC Genomics* 11:117 (V. vinifera leaf disks 24 hpi) | 20167053 | 0 % (TIGR TC → VIT bridge needed) |
| Botrytis | Agudelo-Romero et al. 2015 *J Exp Bot* 66:1769 (Trincadeira berries EL33 + EL35) | 25675955 | 60.2 % |
| Heat | Liu et al. 2012 *BMC Plant Biol* 12:174 (Cabernet Sauvignon leaves) | 23016701 | 0 % (Affy probeset → VIT bridge needed) |

Paper-ready findings from Tier 1:

1. **GRBV replicates across cultivars (Merlot → Zinfandel)**: Fisher p=2.0e-8, log2OR=5.5, Spearman ρ(log2FC)=0.47 against Blanco-Ulate 2017 Zinfandel Oakville Harvest. First cross-cultivar validation.
2. **GRBV ≠ GLRaV-3 at the transcriptome level**: Spearman ρ ≈ 0 against Song 2022; the Merlot GRBV response is not a generic Vitis viral response.
3. **64 % of 55 GRBV DEGs are GRBV-unique**; 24 % overlap with Botrytis (berry-stress core); 3-gene cross-cultivar GRBV core (Merlot ∩ Zinfandel) visible in Fig. S64b.

## Not yet done (Tier 2 roadmap)

- Affymetrix Vitis GeneChip probeset → VIT crosswalk (would unlock Fung 2008 + Liu 2012; ~1–2 days of work or fetch vitis.db Bioconductor annotation if still hosted).
- TIGR Tentative Consensus → VIT crosswalk (unlocks Polesani 2010).
- Raw-count re-analysis from GEO / SRA for all six papers under a unified DESeq2 pipeline. Most robust but 1–2 weeks of work.

# SESSION STATE UPDATE 2026-04-23 (very late): Figure-upgrade pass + 3 new scientific dimensions

## Figure pass summary

Reviewed the full atlas with SciAgent-Skills (Nature figure guide + scientific-visualization + critical-thinking), built 9 upgrades + 7 new supplementary figures + 2 new Main figures, staged in `Figure_upgrades/`, then merged everything into `figures/` and deleted the staging folder.

**Canonical layout after merge:**
- `figures/main/` — 10 original Main + **MoA_GRBV_mechanism** + **Summary_study_design** = 12 PNG + 12 SVG
- `figures/supp/` — all supp up to S72 (includes new S66–S72) = 64 PNG + 64 SVG
- `figures/assets/icons/` — 15 Firefly-generated PNG icons used by Main10, MoA, Summary
- `figures/scripts_R/` — 18 upgraded/new R scripts, `Figure_upgrades` references cleaned
- `figures/UPGRADE_NOTES.md` — changelog + per-figure upgrade notes + three-new-dimensions summary
- `S70_top4_per_category_NES.csv` — companion CSV sits at repo root with the other 30 project CSVs

**Legacy S29 block** in `figures/scripts_R/supp_robust_and_env.R` wrapped in `if (FALSE) {...}` with a pointer to the new standalone `S29_concordance_table.R`.

## Three NEW scientific dimensions unlocked

These are claims the previous manuscript draft did not articulate. Each has a paragraph in `docs/MANUSCRIPT_DRAFT_v1.md` (new Results subsection + new Discussion subsection + expanded Contribution, now 6 items instead of 4). Captured as Thread 10 in `docs/INTEGRATED_ANALYSIS.md` and claims 10–12 in `docs/EXTENDED_RESULTS_REPORT.md`.

### 1. GRBV is a passive-failure pathogen, not an active-suppression pathogen (S66+S67+S72)

- RB(−) vines at 2022 Harvest: sHSP module z = **+1.64** (strong natural heat upregulation).
- RB(+) vines at same cohort × season × stage: sHSP z = **−0.67** (2.3-unit deficit).
- S72 gene-level variance decomposition: sHSPs carry environmental variance share 0.33 (vs 0.19 baseline) but viral share **0.11 ≈ 0.12 baseline** — viral contribution indistinguishable from background noise.
- Reading: GRBV does not repress sHSPs; it prevents the normal HSF-driven upregulation of sHSPs under heat. Wet-lab target shifts from HSP genes themselves to **HSFA2 / MBF1c transactivation machinery**.

### 2. Three temporal phases of chronic infection (S69+S71)

- 2022 cohort (year 1): **214 DEGs**, broad/shallow, ≤ 5× fold enrichment for most curated families.
- 2020 cohort (years 3–5): strongest PC-space drift (**21.9** distance at 2021 Harvest), 157–188 DEGs, moderate enrichment.
- 2016 cohort (years 6–7): smallest DEG set (**105**) but hyper-focused enrichment (**65× sHSP, 138× HSP100/ClpB, 49× GO:heat**).
- Pattern: acute **diffuse** → peak **drift** → chronic **narrow**. GRBV 55-gene signature is the only category persistently enriched 52–100× in every cohort as the shared fingerprint.
- No prior GRBV study has spanned a comparable infection-onset gradient; this is a genuinely new temporal axis.

### 3. Disease cost is temperature-amplified 2.6× at the phenotype level (S66 panel D, S67 panel C)

- RB(−) harvest °Brix: **26.7** in cool 2021, **34.4** in hot 2022.
- RB(+) harvest °Brix: 22.8–23.4 in 2021, 23.5–25.0 in 2022.
- Implied disease cost: **≈ 4 °Bx** cool year, **≈ 10 °Bx** hot year → **2.6× amplification**.
- Translates the transcriptional inverted-trajectory finding directly into grower-legible yield loss.

## Figures that operationalise the new dimensions

| Figure | Role |
|---|---|
| **S66** | 4-layer disease progression cube (titer + sHSP + GRBV sig + Brix) across cohort × season × stage |
| **S67** | RB(−) baseline cube (same grid, uninfected controls) — comparator that makes passive-failure visible |
| **S68** | Season × stage enrichment (top 30 pathways + season-flip + trajectory) |
| **S69** | Per-cohort Fisher enrichment against 13 curated families + per-cohort DEG counts |
| **S70** | 90 literature-curated GRBV/Vitis-virus pathways + 8-condition NES heatmap + trajectory |
| **S71** | PCA-based transcriptional response strength across cohort × season × stage |
| **S72** | Per-gene variance decomposition ternary (Dev / Viral / Env) — Singh 2025 bioRxiv Fig 5 replica |

Plus two new Main figures: **MoA_GRBV_mechanism** (3-zone tissue → phloem → cell schematic), **Summary_study_design** (5-roundel graphical abstract).

## Singh et al. 2025 / 2026 bioRxiv papers (Ritu Singh, Kliebenstein lab, UCD)

Extracted + analysed both PDFs. Vitis is **not** in their 10/15-eudicot panel, so direct Botrytis-on-Vitis comparison is not possible without Tier 2 raw-count work. However:

- Paper 2's **12 core infection-responsive orthogroups** (WRKY75, CIPK11, arogenate dehydratase, endochitinase, cytochrome b561, PMI1, OPT4, SP2-chitinase, formin, UGT71E1, 2 uncharacterised) can be mapped to Vitis via the id_crosswalk (86 AT→LOC mappings already built). Analysis proposed as S73 but not yet built.
- Paper 2's framing — *conserved physiological outcomes executed by lineage-specific regulators* — is the framework for our Discussion of the GRBV response.
- S72's ternary analysis is the direct single-species adaptation of Paper 1's Fig 5 H² ternary.
- Both papers' BioProject: **PRJNA1217477** (embargoed pending publication). To be added to Tier 2 meta-analysis when released.

## Changes on disk since last checkpoint (not yet committed to master)

- All of `figures/` — additions + updates (12 new PNG/SVG in main, 14 new in supp, 15 icons, 18 scripts updated/added, UPGRADE_NOTES.md, README-like top of notes file)
- `docs/MANUSCRIPT_DRAFT_v1.md` — new Results subsection + new Discussion subsection + expanded Contribution (6 items) + 7 new Supp figure captions + word-count update (≈ 5,400 running total, up from 4,545)
- `docs/INTEGRATED_ANALYSIS.md` — Thread 10 (passive failure + temporal phases + °Brix cost)
- `docs/EXTENDED_RESULTS_REPORT.md` — claims 10–12 added
- `docs/PLAIN_LANGUAGE_SUMMARY.md` — to be updated next session
- `S70_top4_per_category_NES.csv` at repo root; mis-categorization of KEGG:Viral life cycle - HIV-1 was fixed (Defense regex tightened to remove |[Vv]iral|virus)

---

# SESSION STATE UPDATE 2026-04-24: Agenda-gap pass (rules 13-19) + S73-S81 figures + PPT consistency / colour fixes

## Why this session existed

User reviewed the analysis against the project's real funding agenda (USDA-ARS / CDFA grapevine virology) and surfaced 7 gaps where the deck had not yet answered the agenda's questions:
1. Per-ORF viral transcriptome dissection (predictor of phenotype, ORF×module coupling)
2. Phenylpropanoid + PR + cell-wall per-gene model (deck only had GSEA-level for these three explicitly-named pathways)
3. Carbon translocation / phloem-loading module (Martínez-Lüscher 2018 hypothesis)
4. Vine-level heterogeneity (44 vines x ≤8 samples)
5. TF activity inference (regulon-based, not expression-only)
6. Change-point / segmented regression on cohort age (the deck's "3 phases" framing)
7. Bayesian hierarchical model with vine random effect

A pre-audit (Explore agent) confirmed that 1 was partially done, 4 was a partial (only repeated-measures correction, no descriptive PCA), 5 was partial (GENIE3 hubs but no activity scores), and 2/3/6/7 were genuinely not done.

## FARM rules 13-19 added and run end-to-end

All 7 rules implemented as `farm_pipeline/scripts/{13..19}_*.{py,R}`, registered as Snakemake rules with SLURM resources, pushed to `~/grbv_extended/`, executed via the existing profile (partition `low`, account `publicgrp`, QOS `publicgrp-low-qos`, `--requeue`), and pulled back to `farm_pipeline/results/`. Several iteration cycles needed to fix: (a) `biorep_id` not in v1 metadata → switched relevant rules to `data/sample_metadata_v2.txt`, (b) `Sampling_Season` dtype mismatch (`int64` vs `float64` → string `"2021"` vs `"2021.0"`), (c) `meta.loc[common]` strips index name → explicit `meta.index.name = "Sample"`, (d) `dplyr::select` masked by `segmented::select`, (e) GRN edges columns are `regulatoryGene`/`targetGene` not `tf`/`target`, (f) decoupleR not in Rlibs (used weighted-sum fallback). All scripts in their patched form are committed.

| Rule | Script | Headline finding |
|---|---|---|
| 13 viral ORF extensions | `13_viral_orf_extensions.py` | V2 (N761_gp1, silencing suppressor) adds significant beyond-titer prediction for **Malic** (LRT p=1.2e-3, β=+271 mg/L per log2 CPM) and **Brix** (p=0.033, β=−0.37 °Bx). V1 similar for Malic (p=4.4e-4). Per-ORF × 13-module Spearman heatmap produced (78 cells, BH-adjusted per ORF). |
| 14 per-gene pathway | `14_pathway_per_gene_model.py` | Phenylpropanoid 31% / Cell wall 26% / PR proteins 20% of curated genes have significant infection main effect (p<0.05). 3-way Inf×Season×Stage interaction in 45–51% of genes across all three categories. Phenylpropanoid is the most-disrupted of the three. |
| 15 carbon translocation | `15_carbon_translocation_module.py` | 77 sugar-transport genes curated; 33 expressed. Module eigengene tracks Stage strongly (ρ=0.94, p=2e-128) but **no main season** (ρ=0.013, p=0.83) and **no main cohort** (ρ=−0.07, p=0.23) effect. Reads as a NULL counter-finding to the Martínez-Lüscher 2018 panel-level hypothesis; SPS1's mediator status remains relevant but is not module-level rate-limiting. |
| 16 vine heterogeneity | `16_vine_heterogeneity.py` | 44 vines clustered by mean-|log2FC| response score: **13 tolerant / 18 intermediate / 13 susceptible**. Tolerant Harvest Brix 24.22 ± 1.09 (n=43 samples), susceptible 23.64 ± 0.87 (n=42). ~0.6 °Bx delta within one vineyard. Publishable on its own as a breeding lead. |
| 17 TF activity | `17_tf_activity.R` | Per-sample TF-activity scores from GENIE3 regulons (top 200 targets per TF). 236 TFs scored. Top hits: LOC100262652 (padj=1.4e-64), LOC100246521 (padj=7.6e-26), **MYB30 down** (padj=8.9e-10) — consistent with classical PR-protein / HR-suppression pathways. Re-ranks the family vs the expression-only DE that drives slides 36–37. |
| 18 segmented cohort | `18_segmented_cohort.R` | best_k = 0 (linear) for ALL 15 responses (PC1, PC2, 13 module eigengenes) by BIC. **The deck's "3 temporal phases" framing is descriptive simplification, not statistically supported phase transitions.** Notable null finding worth a sentence in the manuscript reframe. |
| 19 Bayesian hierarchical | `19_bayesian_hierarchical.R` | brms y ~ Inf + Season + Stage + (1\|vine_id) on 6 high-interest genes; 4 chains × 2000 iter; rhat ≤ 1.005 throughout. Strong posteriors entirely outside ROPE [-0.1, +0.1]: RdRp1 +1.78 [+1.54, +2.02], F-box LRR +2.01 [+1.76, +2.25], glycine-rich −1.33 [-1.65, -1.01], LOC100266309 +0.80 [+0.47, +1.11], SPS1 +0.52 [+0.09, +0.95]. **MBF1c posterior softens: −0.25 [-0.81, +0.33], 19% in ROPE.** Headline mediator's effect is more uncertain than the frequentist sweep suggested. Worth flagging in Methods/Discussion. |

## Supplementary figures S73-S81 (and the S80A/B/C split)

Eight new R scripts in `figures/scripts_R/`, all sourcing `theme_grbv.R`, all producing PNG + SVG:

- **S73** viral ORF extensions (LRT-Δ-R² heatmap + ORF×module Spearman heatmap)
- **S74** Phenyl + PR + Cell-wall summary panel (category percentages + top-15-per-category dotplot)
- **S75** Carbon translocation module trajectory + per-gene Brix correlation
- **S76** Vine heterogeneity dashboard (per-vine response distribution + tolerant-vs-susceptible Brix)
- **S77** Regulon-based TF activity DE (top-25 bar + volcano)
- **S78** Segmented regression null (ΔBIC bar + best_k count)
- **S79** Bayesian forest plot (6 genes with 95% CI + ROPE band)
- **S81** Carbon-translocation full FC heatmap (33 sugar-transport genes × Season × Cohort × Stage; functional subclass faceting)

The original S80 ("phenylpropanoid + PR + cell wall full FC heatmap, all 252 genes in one figure") had three crowding problems: (a) facet strip labels clipped to "mpling_Season: 2..." / "turity_stage: Pre-V..." due to overlapping ggplot strip rows, (b) 161 cell-wall labels stacked into ~70 px were unreadable, (c) one-figure layout couldn't accommodate the 3-category density spread. **Split into S80A + S80B + S80C** with the following fixes applied to all three:

- **S80A — Phenylpropanoid + PR proteins**: 91 genes (51 + 40) in 2 stacked panels, 10 × 14 in @ 600 dpi. All 91 labels at 6 pt. Hierarchical clustering per panel (Euclidean / Ward.D2). `facet_nested(. ~ Season + Stage)` via ggh4x — fixes the strip clipping.
- **S80B — Cell wall full landscape**: 161 genes split into 5 functional sub-categories via GFF regex (Pectin / Hemicellulose / Cellulose / Lignin / Expansin). 11 × 18 in @ 600 dpi. **Only the 13 genes with |log2FC| ≥ 1 in ≥ 4 of 8 outer cells get y-axis labels**; the other 148 are anonymous tile rows for the global pattern.
- **S80C — Cell wall focused subset**: 21 genes selected from the 161 by either `|FC| ≥ 1 in ≥ 4/8 cells` (13 strict) **or** `|FC| ≥ 0.45 in ≥ 6/8 cells` (8 broad; threshold auto-tuned over [0.30, 0.95] to land in the 20–40 window because cell-wall response is muted). 10 × 9 in @ 600 dpi. Row labels are `Symbol — short_function` from `farm_pipeline/results/annotations/gene_annotations.csv` where available, `LOC###### — best NCBI description` otherwise — **no bare LOC IDs**. Sub-category rendered as a left-side **nested y-axis label** via `ggh4x::scale_y_discrete(guide = guide_axis_nested(delim = ">"))`. Cell text overlay for `|log2FC| ≥ 1.5` (white when `|log2FC| > 2`, black otherwise). Subtle cell borders (`geom_tile(color = "grey90", linewidth = 0.1)`).

Common to all three: `scale_fill_gradient2(limits = c(-3, 3), oob = scales::squish)`, 2-digit cohort years on x at 90°, captions stating the cohort convention (`'16` = 5–6 yrs since infection at 2021 sampling, 9–10 yrs at 2022 sampling; `'22` = 0–1 yr).

## PPT (`GRBV_Funding_Committee_Presentation.pptx`, gitignored)

The deck went from **71 slides → 82 slides** over four scripted updates this session:

| Step | Script | Net change |
|---|---|---|
| Consistency-fix pass (footers, sidebar labels, cross-refs, bridge sentence on slide 29, deliverable counts, slide 14 arithmetic, slide 64 panel callout) | `scripts/update_ppt_consistency_fixes.py` | 0 |
| 7 agenda-gap slides (S73–S79) inserted after slide 64 | `scripts/update_ppt_agenda_gaps.py` | +7 → 78 |
| Full-resolution heatmap slides (initial S80 + S81) | `scripts/update_ppt_full_heatmaps.py` | +2 → 80 |
| Split S80 → S80A + S80B | `scripts/update_ppt_split_S80.py` | +1 → 81 |
| Insert S80C after S80B | `scripts/update_ppt_S80C.py` | +1 → 82 |

Sidebar labels for the agenda-gap block are `6 / Gap-fill`, distinguishing them from the 4 existing post-viral PLUS sub-categories (Biochem / Tier 1 / Cohort / Dimensions). Slide-4 roadmap auto-recomputed each step; current ranges: Section 7 Synthesis 76–77, Section 8 Reanalysis 78–82.

Typography (after the consistency-fix pass + the black-on-navy title fix the user spotted):

- 0 Arial runs across the 636-run deck (all body text now Cambria explicit).
- All 18 content-slide titles that were rendering invisible BLACK on navy header bands now have explicit `#FFFFFF` white.
- Body Cambria 13 pt navy (#0F172A) on the original slides; Cambria 10 pt navy on the 16 dense post-figure-upgrade slides (kept at 10 pt because raising to 13 pt would overflow 16 of 18 fixed-height boxes).
- PLUS slides title Cambria 16 pt white, tag Cambria 9 pt white.

## Atlas

`figures/GRBV_Full_Figure_Atlas.pdf` regenerated three times this session:

- After S73–S79: 89 pages (10 main + 79 supp), 62 MB.
- After S80 + S81 added: 91 pages.
- After S80 → S80A + S80B + S80C split + S81 retained: **93 pages (10 main + 83 supp)**.

The atlas builder (`figures/scripts_R/build_full_atlas.py`) was patched so `resolve_supp(n, fname)` first tries direct `SUPP_DIR/fname` lookup (handles S80A / S80B / S80C alongside the existing prefix-glob fallback).

## Notable findings to flag for the manuscript

1. **Rule 18 null** kills the "3 discrete temporal phases of chronic infection" framing as a statistical claim. Reframe to **continuous gradient with descriptive cohort-bin labels** (acute / mid / chronic stay as descriptive shorthand). Notable for the reviewer-response: "We tested for break-points formally and prefer the linear model; phase labels in the figures are descriptive only."
2. **Rule 19 MBF1c uncertainty.** The Bayesian random-intercept model is more skeptical of MBF1c's effect than the frequentist mediation sweep. Posterior straddles zero; 19% in ROPE. This needs a one-sentence caveat in Methods/Discussion before submission. Re-elevates the qPCR / VIGS validation priority.
3. **Rule 16 tolerant-vs-susceptible.** Within the same Merlot vineyard there's a Brix-defined tolerance gradient (~0.6 °Bx delta). Publishable on its own if any of the tolerant vines have heritable differences (rootstock or transcriptional state) — direct breeding lead.
4. **Rule 13 V2-beyond-titer for Malic.** Concrete answer to project Objective 2 ("impact of viral titer on disease outcomes"): viral mRNA load explains malic-acid variance that titer alone misses. V2 (N761_gp1) is the strongest predictor.
5. **Rule 15 NULL on sugar module.** Martínez-Lüscher 2018's "GRBV blocks hexose translocation" hypothesis is **not** supported at the panel level in our data. SPS1's individual mediator role and per-gene Brix correlations remain valid; the panel-level framing is retracted.

## Why log2 CPM is the per-gene-model y-axis (asked by user during this session)

Rule 14's per-gene OLS uses `y = log2(CPM + 1)` as the response. Three reasons that scale: (a) `counts / library_size × 1e6` strips the library-size artefact; (b) `log2` stabilises the variance (raw counts are Poisson, violating OLS homoscedasticity); (c) a coefficient of `+1` = a 2-fold expression change, the same scale `limma-voom` and `DESeq2` use. The `+1` pseudo-count prevents `log(0)` on zero-count cells.

## Documents in `docs/`

| File | This session's update |
|---|---|
| `INTEGRATED_ANALYSIS.md` | Thread 11 added covering the 7 agenda-gap claims + S80A/B/C + S81 (this session). |
| `EXTENDED_RESULTS_REPORT.md` | New section "Claims 13–19" appended (this session). |
| `MANUSCRIPT_DRAFT_v1.md` | Not edited — too early; the Bayesian / segmented null reframes need to land here on the next pass. |
| `HSFA2_VIGS_EXPERIMENT_PLAN.md` | Not edited — the MBF1c uncertainty (rule 19) elevates the priority but doesn't change the protocol. |
| `PLAIN_LANGUAGE_SUMMARY.md` | Not edited this session. |
| `figures/UPGRADE_NOTES.md` | Pending — will note S73–S81 + S80A/B/C split next session. |

## Next-action priorities (post-session)

1. **Reframe slide 63 ("New dimension 2: chronic infection has three temporal phases")** in light of the rule 18 null. Keep the phases as descriptive bins; add a sentence "best_k = 0 by BIC; gradient is continuous, see S78" so the figure does not over-claim discreteness.
2. **Add an MBF1c caveat to slide 30 (Q4 mechanism)**: "Bayesian random-intercept reanalysis (S79) softens the MBF1c posterior to −0.25 [-0.81, +0.33] with 19% in ROPE. Effect direction needs qPCR / VIGS to confirm before manuscript submission."
3. **Manuscript draft** — fold rules 14, 16, 19 into Discussion; rules 13, 18 into Methods.
4. **PPT visual QA** still pending: LibreOffice not installed locally; `brew install --cask libreoffice` would unblock per-slide JPG rendering for clipping / overlap checks. User has not authorised the install.
5. **Tier-2 cross-cohort meta-analysis**: bridge probe-based published datasets (Fung 2008, Polesani 2010, Liu 2012) into VIT-locus space to unlock the 0%-coverage Tier 1 papers (1–2 days; deferred from the agenda-gap pass).

## Git push log this session

| Commit | Subject |
|---|---|
| 8f4825a | PPT consistency + style fixes (71-slide deck) |
| 0b7bc21 | FARM rules 13-19: 7 agenda-gap analyses |
| 99e19d2 | S73-S79 + 7 PPT slides: agenda-gap supplementary figures and deck section |
| a9da32d | S80 + S81: full-resolution log2-FC heatmaps for the agenda-named pathways |
| 42d652a | Split S80 into S80A (Phenyl + PR) + S80B (Cell wall) with readability fixes |
| fbb45fe | S80C — focused cell-wall heatmap with cell-text annotation |

All on `master` at https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis.

---

# SESSION STATE UPDATE 2026-05-26: Pathway figures (S87-S95)

## What was added

Nine new supplementary figures rendering pathway-level GRBV impact:

- **S87 Phenylpropanoid biosynthesis** (KEGG vvi00940)
- **S88 Starch and sucrose metabolism** (vvi00500)
- **S89 Plant hormone signal transduction** (vvi04075)
- **S90 Ribosome** (vvi03010)
- **S91 Biosynthesis of amino acids** (vvi01230)
- **S92 Proteasome** (vvi03050)
- **S93 HSF -> MBF1c -> sHSP axis** (custom ggplot flow)
- **S94 Temperature -> titer -> host mediator cascade** (custom ggplot flow)
- **S95 Viral ORFs -> WGCNA module rewiring** (custom bipartite)

All 6 KEGG pathways rendered via pathview 1.48.0 with 8-cell split-box log2FC overlay in pure stage order (2021 Pre / Ver / Post / Harv, then 2022 Pre / Ver / Post / Harv). Palette: #2166ac (down) / #f7f7f7 (neutral) / #b2182b (up). log2FC clipped at +/-3. Paralog rule: `node.sum = "mean"` (mean across all LOC genes mapped to each EC node); full paralog list per pathway in `figures/supp/S{nn}_paralog_counts.csv`.

## Implementation

- `figures/scripts_R/pathway_overlay_driver.R`: driver with `read_q2x_matrix()`, `compute_mapping_coverage()`, `write_paralog_auxiliary()`, `render_kegg_pathway()`
- `figures/scripts_R/test_pathway_driver.R`: assertion-based smoke tests for driver helpers
- 6 shim scripts (S87-S92) calling the driver
- 3 standalone custom scripts (S93-S95)
- `figures/cache/kegg_xml/`: KGML XML cache (6 files, ~250 KB total)

## Atlas

`figures/GRBV_Full_Figure_Atlas.pdf` rebuilt: **105 pages** (10 main + 95 supp).

## Spec + plan

- Spec: `docs/superpowers/specs/2026-05-26-pathway-figures-design.md`
- Plan: `docs/superpowers/plans/2026-05-26-pathway-figures.md`

## Accepted trade-offs

- Pathview's KEGG PNG output uses KEGG default Helvetica + KEGG-style borders; does not match `theme_grbv` Source Sans Pro. Accepted: canonical KEGG layout is the value here.
- Pathview averages paralogs at multi-LOC EC nodes via `node.sum = "mean"`. Full paralog set written to auxiliary CSV; caption must note the rule when figures land in the manuscript.
- Plant hormone signal transduction (vvi04075) has only 2/8 GSEA-significant conditions; included for narrative consistency with S86, not for primary GSEA signal.
- S90 Ribosome and S92 Proteasome produced no graphviz PDF (pathview rejects "parts-list" pathways with no enzyme edges); native PNGs are the canonical outputs for those two.
- The atlas builder dict `STATUS_COLOR` + badge dict gained a `"new"` key to support S93-S95 ("new" status). Committed in `a6464be`.

## Commits

| SHA | Subject |
|---|---|
| 38272ee | verify LOC->KEGG vvi mapping; add KGML cache scaffold |
| 54bbe3b | driver scaffold + read_q2x_matrix() with smoke test |
| 2650e9b | remove dead-code tibble-rowname assignment |
| 93a8cc5 | compute_mapping_coverage() + KGML parsing |
| 280f768 | write_paralog_auxiliary() per-EC inventory |
| 9975f4b | render_kegg_pathway() wrapper with native + graphviz passes |
| 2f23f91 | S87 Phenylpropanoid biosynthesis (vvi00940) [includes <<- driver fix] |
| 80bda01 | S88 Starch and sucrose metabolism (vvi00500) |
| be5d38c | S89 Plant hormone signal transduction (vvi04075) |
| 1bd9e6f | S90 Ribosome (vvi03010) |
| 6211261 | S91 Biosynthesis of amino acids (vvi01230) |
| 1004e32 | S92 Proteasome (vvi03050) |
| 485f5fe | S93 HSF -> MBF1c -> sHSP axis (custom flow) |
| adbbe02 | S94 Temperature -> titer -> mediator cascade (custom flow) |
| dd18397 | S95 Viral ORF -> module rewiring (custom bipartite) |
| 19f562c | register S87-S95 in atlas builder |
| b894d82 | rebuild atlas to 105 pages with S87-S95 |
| a6464be | add 'new' status key to atlas builder for S93-S95 |

