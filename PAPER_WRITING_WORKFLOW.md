# Paper-Writing Workflow: GRBV Transcriptomics Study

## Quick Reference

| Item | Status |
|------|--------|
| **Dataset** | 327 samples, 14,604 genes, 2 seasons, 4 maturity stages |
| **Figures ready** | 35 publication-quality PNGs |
| **CSV tables** | 13 supplementary data files |
| **Draft report** | GRBV_Transcriptomics_Report.md (500 lines) |
| **PPT** | GRBV_Transcriptomics_Presentation.pptx (auto-generated) |

---

## Phase 1: Organize Your Story Arc (Day 1)

### 1.1 Define the Central Narrative

Your data supports this story:

> **GRBV hijacks host transcription in a development- and environment-dependent manner, suppressing heat-stress defenses, inverting seasonal response kinetics, and targeting genes that co-localize with known disease-resistance and berry-quality QTLs — with implications for diagnostics and breeding.**

### 1.2 Recommended Paper Structure

| Section | Key Content | Figures | Tables |
|---------|-------------|---------|--------|
| **Abstract** | Multi-dimensional RNA-Seq; 6 research questions; inverted seasonal trajectory; biomarkers; QTL co-localization; heat-stress suppression | — | — |
| **Introduction** | GRBV biology, RBD symptoms, knowledge gaps, study objectives | Fig 23 (schematic) | — |
| **Methods** | RNA-Seq pipeline, DE analysis, pathway enrichment, biomarker scoring, QTL cross-ref, heat-stress analysis | — | Table 1 (dataset summary) |
| **Results 1**: Global Overview | PCA, variance hierarchy (stage > season > infection) | Fig 1 | — |
| **Results 2**: Infection Duration (Q1) | Biphasic pattern, early vs late responders | Fig 2, Fig 8 | Supp Table S1 |
| **Results 3**: Maturity Stage (Q2) | Stage-dependent response, Post-veraison/Harvest peaks | Fig 3, Fig 9 | Supp Table S2 |
| **Results 4**: Seasonal & Environmental (Q3-Q5) | Inverted trajectory (KEY FINDING), environment modulation | Fig 4-6, Fig 10-12, Fig 15-19 | Supp Table S3-S5 |
| **Results 5**: Infection Trajectory (Q6) | Biphasic: acute defense → chronic metabolism | Fig 7, Fig 13 | Supp Table S6 |
| **Results 6**: Pathway Landscape | Phenylpropanoid, plant-pathogen interaction, ABC transporters | Fig 14 | Supp Table S7 |
| **Results 7**: Heat-Stress Signatures | 115 heat genes, sHSP suppression by GRBV, virus-heat disentanglement | Fig 34, Fig 35 | heat_stress_DEGs_all_conditions.csv |
| **Results 8**: Environment-Titer Integration | Temperature-titer-DEG correlations | Fig 30-33 | integrated_environment_titer_DEG.csv |
| **Results 9**: Biomarker Discovery | Top 50 biomarkers, 3-gene panel (AUC=1.0), classification | Fig 20-22 | GRBV_biomarkers_top50.csv |
| **Results 10**: GWAS/QTL Cross-Reference | 60% co-localization, disease resistance clusters, sugar QTL | Fig 24-25 | GRBV_biomarker_QTL_crossref.csv |
| **Discussion** | Developmental modulation, seasonal inversion mechanism, heat-stress suppression, breeding implications, diagnostic timing | Fig 26-29 (summary dashboards) | — |
| **Conclusions** | 6 key conclusions | — | — |

### 1.3 Suggested Journal Targets

| Journal | Impact | Fit | Notes |
|---------|--------|-----|-------|
| *Plant Cell* | ~17 | Excellent | Multi-dimensional systems biology, mechanistic insight |
| *New Phytologist* | ~10 | Excellent | Plant-pathogen interaction, environmental modulation |
| *Plant Physiology* | ~8 | Very Good | Comprehensive transcriptomics, heat-stress biology |
| *Molecular Plant-Microbe Interactions* | ~5 | Very Good | Virus-host interaction focus |
| *Frontiers in Plant Science* | ~5 | Good | Broad plant biology, open access |
| *Phytopathology* | ~3 | Good | Plant disease focus |

---

## Phase 2: Figure Selection & Refinement (Day 2-3)

### 2.1 Main Text Figures (recommended 6-8 for most journals)

| Priority | Figure | File | Why Include |
|----------|--------|------|-------------|
| **1** | Analysis Pipeline | `23_analysis_schematic.png` | Orients the reader |
| **2** | Global PCA | `01_global_pca_overview.png` | Shows variance hierarchy |
| **3** | Seasonal Maturity DEG Comparison | `15_seasonal_maturity_DEG_comparison.png` | KEY FINDING: inverted trajectory |
| **4** | Master Pathway Landscape | `14_master_pathway_landscape.png` | Systems-level view |
| **5** | Heat vs Virus Disentanglement | `35_heat_vs_virus_disentanglement.png` | Novel biology |
| **6** | Biomarker Overview | `20_biomarker_overview.png` | Translational impact |
| **7** | QTL Chromosome Map | `24_biomarker_QTL_chromosome_map.png` | Breeding relevance |
| **8** | Summary Dashboard | `29_comprehensive_summary_dashboard.png` | Integrative conclusion |

### 2.2 Supplementary Figures (remaining 27)

All other figures (02-07, 08-13, 16-19, 21-22, 25-28, 30-34) go as supplementary.

### 2.3 Figure Refinement Checklist

For each main-text figure:
- [ ] Resolution >= 300 DPI (current: 150 DPI — re-render at 300)
- [ ] Font size >= 8pt in all panels
- [ ] Consistent color scheme across figures
- [ ] Panel labels (A, B, C, D) clearly visible
- [ ] Axis labels fully spelled out (no abbreviations)
- [ ] Legend positioned without overlapping data
- [ ] Color-blind friendly palette (check with Coblis simulator)

**To re-render at 300 DPI**: In each .py script, change `'figure.dpi': 150` to `'figure.dpi': 300`.

---

## Phase 3: Write the Manuscript (Day 3-10)

### 3.1 Writing Order (most efficient)

1. **Methods** — straightforward, sets the foundation (use GRBV_Transcriptomics_Report.md Section 2)
2. **Results** — present findings in logical order (use Report Section 3 as draft)
3. **Discussion** — interpret results, compare to literature
4. **Introduction** — frame the problem knowing what you found
5. **Abstract** — summarize after everything else is written
6. **Title** — finalize last

### 3.2 Key Numbers to Cite

```
Dataset:         327 samples, 14,604 genes (from 28,168 raw)
Infection groups: RB(+) infected vs RB(-) uninfected controls
Seasons:         2021 (cooler) and 2022 (hotter, +7.6C at Harvest)
Stages:          Pre-veraison, Veraison, Post-veraison, Harvest
Cohorts:         7 infection-year cohorts (2016-2022)

DE thresholds:   padj < 0.05, |log2FC| >= 1
Enrichment:      Fisher's exact, BH FDR correction

Key DEG counts:
  - 2021 Harvest: 1,192 DEGs (peak in 2021)
  - 2022 Post-veraison: 565 DEGs (peak in 2022)
  - Cross-season correlation: Pearson r = 0.33-0.69

Heat-stress:
  - 115 heat-stress genes expressed (9 categories)
  - 57 heat genes are GRBV DEGs
  - 29/36 sHSPs downregulated at 2021 Harvest in RB(+)
  - 23 heat genes with BOTH virus and environment effects at Harvest
  - HSF-sHSP module correlation: r=0.625, p=7.2e-37

Biomarkers:
  - 12 universal biomarkers (robust across all conditions)
  - 3-gene panel AUC = 1.000
  - 60% of biomarkers near known QTLs (within 1 Mb)
  - 18 biomarkers in disease-resistance clusters
  - Top host biomarker: F-box/LRR protein (AUC = 0.880)

Environmental:
  - 2022 Harvest was +7.6C warmer than 2021
  - sHSP module score: +1.64 (2022 Harvest RB-) vs +0.46 (2021 Harvest RB-)
  - GRBV reduces sHSP induction by ~50% (RB+ vs RB- at Harvest)
```

### 3.3 Manuscript Draft Sources

| Section | Source File | Lines |
|---------|-----------|-------|
| Introduction | GRBV_Transcriptomics_Report.md | 29-41 |
| Methods | GRBV_Transcriptomics_Report.md | 44-94 |
| Results (all Q1-Q6) | GRBV_Transcriptomics_Report.md | 98-394 |
| Discussion | GRBV_Transcriptomics_Report.md | 396-425 |
| Conclusions | GRBV_Transcriptomics_Report.md | 435-448 |
| Heat stress results | heat_stress_signatures.py output | (run script) |
| Env-titer results | environment_titer_integration.py output | (run script) |

### 3.4 Literature to Cite

Search for these key references:
- GRBV original characterization (Al Rwahnih et al., 2013)
- GRBV genome and classification (Varsani et al., 2017)
- GRBV transcriptomics (Blanco-Ulate et al., 2017; Rumbaugh et al., 2021)
- GRBV impacts on berry metabolism (Girardello et al., 2020)
- Grapevine heat-stress HSPs (Liu et al., 2012; Wang et al., 2014)
- Grapevine QTL/GWAS resistance loci (Zyprian et al., 2016; Possamai et al., 2020)
- Biomarker RT-qPCR diagnostics in plants (various)

---

## Phase 4: Supplementary Materials (Day 8-12)

### 4.1 Supplementary Tables

| Table | Source CSV | Description |
|-------|-----------|-------------|
| Table S1 | GRBV_biomarkers_top50.csv | Complete biomarker list with scores |
| Table S2 | GRBV_biomarker_QTL_crossref.csv | QTL cross-reference details |
| Table S3 | universal_DEGs_across_approaches.csv | 908 robust DEGs |
| Table S4 | heat_stress_DEGs_all_conditions.csv | Heat-gene DE results |
| Table S5 | heat_gene_temperature_correlations.csv | Temperature correlations |
| Table S6 | heat_DEGs_recurrence.csv | Recurrent heat DEGs |
| Table S7 | Q1_pooled_vs_stratified.csv | Analytical robustness |
| Table S8 | Q2_pooled_vs_stratified.csv | Analytical robustness |
| Table S9 | integrated_environment_titer_DEG.csv | Env-titer-DEG integration |

### 4.2 Data Availability Statement

```
Raw RNA-Seq data have been deposited in [NCBI GEO/SRA under accession GSEXXXXX].
All analysis scripts, processed count matrices, and supplementary data files are
available at [GitHub repository URL]. The Vitis vinifera functional annotation
database (org.Vvinifera.eg.db) was obtained from Bioconductor.
```

---

## Phase 5: Review & Polish (Day 12-15)

### 5.1 Internal Review Checklist

- [ ] All figure references in text match actual figure numbers
- [ ] All p-values and fold-changes verified against script outputs
- [ ] Methods section matches what scripts actually do (check thresholds)
- [ ] Supplementary table column headers explained
- [ ] No orphan references (every citation in text is in bibliography)
- [ ] Author contributions drafted
- [ ] Acknowledgments and funding sources listed
- [ ] Conflict of interest statement prepared

### 5.2 Common Reviewer Concerns to Address Proactively

| Concern | Pre-emptive Response |
|---------|---------------------|
| "Why Welch's t-test instead of DESeq2/edgeR?" | Discuss in Methods; note large sample size (n=327) makes parametric tests robust; show supplementary comparison if possible |
| "Multiple testing burden across 6 questions?" | Each question uses independent BH correction; universal DEGs (Table S3) validated across approaches |
| "Batch effects between seasons?" | Addressed by Q5 (baseline environmental analysis in RB- controls); PCA shows no obvious batch clustering |
| "Small sample size for some sub-groups?" | Minimum n=3 per group; acknowledge as limitation; stratified analysis confirms pooled results |
| "Biomarker validation?" | Discuss RT-qPCR validation as future work; cross-validation already performed (5-fold stratified) |

---

## Phase 6: Generate Presentation (Automated)

Run: `python3 generate_presentation.py`

This produces `GRBV_Transcriptomics_Presentation.pptx` with all 35 figures organized into a conference-ready slide deck.

---

## File Inventory

### Analysis Scripts (10)
1. `analysis.py` — Core DE analysis (Q1-Q6), PCA, Figs 1-7
2. `pathway_enrichment.py` — GO/KEGG enrichment, Figs 8-14
3. `seasonal_maturity_progression.py` — Season x Stage analysis, Figs 15-19
4. `biomarker_discovery.py` — Biomarker scoring & classification, Figs 20-22
5. `gwas_qtl_crossref.py` — QTL cross-reference, Figs 24-25
6. `stratified_vs_pooled_comparison.py` — Robustness check, Figs 26-29
7. `environment_titer_integration.py` — Weather + titer, Figs 30-33
8. `heat_stress_signatures.py` — Heat-stress signatures, Figs 34-35
9. `analysis_schematic.py` — Pipeline diagram, Fig 23
10. `generate_report_pdf.py` — PDF report generator

### Key Output Files
- `GRBV_Transcriptomics_Report.md` — Full draft report
- `GRBV_Transcriptomics_Presentation.pptx` — Slide deck (auto-generated)
- 35 PNG figures (01-35)
- 13 CSV supplementary data files
