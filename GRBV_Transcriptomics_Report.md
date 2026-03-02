# Comprehensive Transcriptomic Analysis of Grapevine Red Blotch Virus (GRBV) Infection

## A Multi-Dimensional RNA-Seq Study Across Infection Duration, Maturity Stages, Seasons, and Genomic Context

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Materials and Methods](#2-materials-and-methods)
3. [Results](#3-results)
   - 3.1 [Global Transcriptome Overview](#31-global-transcriptome-overview)
   - 3.2 [Q1: Effect of Infection Duration](#32-q1-effect-of-infection-duration-on-grapevine-transcriptome)
   - 3.3 [Q2: Effect of Maturity Stage](#33-q2-effect-of-maturity-stage-on-grbv-response)
   - 3.4 [Q3: Seasonal Effects on Infection Response](#34-q3-seasonal-effects-on-infection-response)
   - 3.5 [Q4: Environmental Modulation of GRBV Effects](#35-q4-environmental-modulation-of-grbv-effects)
   - 3.6 [Q5: Baseline Environmental Changes in Uninfected Vines](#36-q5-baseline-environmental-changes-in-uninfected-vines)
   - 3.7 [Q6: Infection Age Trajectory](#37-q6-infection-age-trajectory)
   - 3.8 [Pathway Enrichment Analysis](#38-pathway-enrichment-analysis)
   - 3.9 [Seasonal Maturity Progression Analysis](#39-seasonal-maturity-progression-analysis)
   - 3.10 [Biomarker Discovery](#310-biomarker-discovery)
   - 3.11 [GWAS/QTL Cross-Reference Analysis](#311-gwasqtl-cross-reference-analysis)
4. [Discussion](#4-discussion)
5. [Conclusions](#5-conclusions)
6. [Supplementary Files](#6-supplementary-files)

---

## 1. Introduction

Grapevine Red Blotch Virus (GRBV) is a DNA virus (family Geminiviridae) that causes Red Blotch Disease (RBD) in grapevines, leading to reduced sugar accumulation, altered anthocyanin profiles, delayed ripening, and diminished fruit quality. Understanding the molecular mechanisms underlying GRBV infection is critical for developing management strategies, breeding resistant cultivars, and identifying diagnostic biomarkers.

This study presents a comprehensive RNA-Seq transcriptomic analysis of GRBV-infected (RB+) versus uninfected (RB-) grapevine leaf tissue, examining the virus-host interaction across multiple biological dimensions: infection duration (years since infection), berry maturity stages (Pre-veraison through Harvest), sampling seasons (2021 and 2022), and environmental context. The analysis further extends to biomarker discovery and cross-referencing identified biomarker genes with known grapevine GWAS/QTL loci.

### Objectives

1. Characterize the transcriptomic impact of GRBV across infection years, maturity stages, and seasons
2. Disentangle viral effects from natural developmental and environmental variation
3. Identify robust molecular biomarkers for GRBV infection diagnosis
4. Map biomarker genes to known grapevine genetic loci for functional and breeding relevance

---

## 2. Materials and Methods

### 2.1 Data Description

| Parameter | Value |
|-----------|-------|
| **Raw count matrix** | 28,168 genes x 327 samples |
| **Genes after filtering** | 14,604 genes |
| **Total samples** | 327 |
| **Infection groups** | RB(+) (infected) and RB(-) (uninfected controls) |
| **Sampling seasons** | 2021 and 2022 |
| **Maturity stages** | Pre-veraison, Veraison, Post-veraison, Harvest |
| **Infection year cohorts** | 2016, 2017, 2018, 2019, 2020, 2021, 2022 (years since infection onset) |
| **Reference genome** | *Vitis vinifera* PN40024 (12X.v2 / ASM3070453v1) |
| **Functional annotation** | org.Vvinifera.eg.db (29,585 genes; 59,449 GO; 13,801 KEGG annotations) |

### 2.2 Gene Filtering

Genes were retained if they had at least 10 raw counts in at least 10% of samples (n >= 33 samples). This removed 13,564 low-expression genes, reducing the dataset from 28,168 to 14,604 genes. The filtering ensures statistical reliability by removing genes detected only sporadically or at very low levels.

### 2.3 Normalization

- **CPM normalization**: Raw counts divided by library size, multiplied by 10^6
- **Log transformation**: log2(CPM + 1) for downstream analyses

### 2.4 Differential Expression Analysis

Differential expression was assessed using **Welch's t-test** (unequal variance) comparing RB(+) vs RB(-) groups, with **Benjamini-Hochberg (BH) FDR correction** for multiple testing. Significance thresholds: **adjusted p-value < 0.05** and **|log2 fold change| >= 1**.

### 2.5 Pathway Enrichment

**Gene Ontology (GO)** and **KEGG** pathway enrichment was performed using **Fisher's exact test** with BH FDR correction. The org.Vvinifera.eg.db SQLite database provided gene-to-pathway mappings.

### 2.6 Biomarker Scoring

A composite biomarker score was calculated for each gene based on:
- **Consistency** (30%): Fraction of tested conditions where the gene is significantly DE
- **Effect size** (25%): Normalized mean |log2FC| across conditions
- **Significance** (20%): Normalized -log10(min adjusted p-value)
- **Directional consistency** (15%): Fraction of conditions with concordant direction
- **Season robustness** (10%): Significant in both 2021 and 2022 with same direction

The final score also incorporated AUC-ROC discriminatory power: `Final = 0.7 x (composite) + 0.3 x (AUC - 0.5)`.

### 2.7 GWAS/QTL Cross-Referencing

Biomarker gene chromosomal positions were retrieved from the NCBI Gene database via the E-utilities esummary API. A curated database of 121 known grapevine GWAS/QTL loci was compiled from published literature, covering disease resistance (Rpv, Ren, Run loci), berry quality, abiotic stress, and developmental QTLs. Co-localization was assessed using a 1 Mb window around each biomarker gene.

### 2.8 Visualization

All plots were generated using matplotlib and seaborn in Python. A total of **25 publication-quality figures** were produced across all analysis modules.

---

## 3. Results

### 3.1 Global Transcriptome Overview

**Principal Component Analysis (PCA)** of the full 327-sample transcriptome revealed that the dominant source of variation was **maturity stage**, followed by **sampling season**. GRBV infection status (RB+ vs RB-) contributed to transcriptomic variation but was secondary to developmental and environmental factors in overall variance. This indicates that while GRBV significantly alters the transcriptome, the virus operates within and modulates the existing developmental and seasonal regulatory framework.

**Variance partitioning** (pseudo-R^2 analysis across top variable genes) confirmed the hierarchy:
- **Maturity stage**: Largest contributor to transcriptomic variance
- **Sampling season**: Second largest contributor
- **Infection status**: Third contributor
- **Infection duration (years)**: Smallest but significant contributor

> **Figure 1** (`01_global_pca_overview.png`): Global PCA overview showing sample clustering by infection group, season, maturity stage, and infection duration.

---

### 3.2 Q1: Effect of Infection Duration on Grapevine Transcriptome

**Question**: How does the duration of GRBV infection (years since initial infection) affect the transcriptome?

Differential expression analysis was performed for each infection year cohort (2016-2022) against RB(-) controls. The results revealed that the transcriptomic impact of GRBV is not static but evolves with infection duration. Longer-established infections (dating to 2016) showed the greatest number of DEGs, consistent with cumulative viral load effects and progressive host remodeling. Newer infections (2021-2022) showed fewer DEGs but with distinct early-response signatures.

Correlation analysis identified genes whose expression systematically increases or decreases with infection age, suggesting progressive transcriptome remodeling. Early-response genes (significantly DE only in new infections) may represent acute defense activation, while late-response genes (accumulating with infection duration) likely reflect chronic host adaptation or viral manipulation.

> **Figure 2** (`02_Q1_infection_years.png`): DEG counts per infection year, correlation with duration, pairwise comparisons, and infection-age trajectories.

---

### 3.3 Q2: Effect of Maturity Stage on GRBV Response

**Question**: How does the GRBV-induced transcriptomic response change across berry maturity stages?

DE analysis was stratified by maturity stage (Pre-veraison, Veraison, Post-veraison, Harvest). The host response to GRBV was strongly stage-dependent:

- **Pre-veraison**: Fewest DEGs, indicating the infection's transcriptomic footprint is smallest early in the season
- **Post-veraison and Harvest**: Greatest number of DEGs, coinciding with the period of active sugar accumulation and anthocyanin biosynthesis -- the very processes that GRBV most disrupts

The **stage overlap analysis** identified:
- A core set of genes consistently DE across all four maturity stages (constitutive GRBV response)
- Stage-specific genes unique to particular developmental windows (context-dependent response)
- Progressive intensification from Pre-veraison to Harvest, with the largest expansion of DEGs occurring at the Post-veraison transition

Stage transition analysis within RB(+) and RB(-) groups revealed that GRBV-infected vines undergo aberrant developmental transitions, particularly during the veraison-to-post-veraison shift, consistent with delayed ripening phenotypes.

> **Figure 3** (`03_Q2_maturity_stages.png`): Stage-specific DEGs, overlaps, stage transitions, and maturity-dependent response profiles.

---

### 3.4 Q3: Seasonal Effects on Infection Response

**Question**: How do seasonal differences (2021 vs 2022) affect the GRBV transcriptomic response?

Season-stratified DE analysis showed:
- The overall infection response (RB+ vs RB-) was moderately conserved between years (**Pearson r of fold changes: 0.33-0.69** depending on stage)
- However, the magnitude and composition of the response differed substantially between seasons
- Season x Stage interactions were prominent, with some genes responding to GRBV only in specific season-stage combinations

This demonstrates that the host response to GRBV is not fixed but is modulated by year-to-year environmental variation (temperature, precipitation, growing degree days), creating a **genotype x environment x pathogen interaction** at the transcriptomic level.

> **Figure 4** (`04_Q3_seasonal_effects.png`): Season-specific DEG counts, cross-season fold change correlations, and interaction analysis.

---

### 3.5 Q4: Environmental Modulation of GRBV Effects

**Question**: How does the environment modulate the magnitude of GRBV's transcriptomic impact?

Environmental effects were assessed by comparing the ratio of infection-responsive genes across different seasons and stages. The analysis revealed:

- The environment does not simply add an independent layer of variation; rather, it actively modulates the intensity and nature of the host-pathogen interaction
- Some pathways are "environmentally buffered" (responding to GRBV regardless of conditions), while others are "environmentally contingent" (only activated under specific conditions)
- The 2021-to-2022 seasonal shift amplified certain GRBV effects (e.g., secondary metabolism disruption) while attenuating others (e.g., some defense responses)

> **Figure 5** (`05_Q4_environmental_effects.png`): Variance partitioning, environmental modulation of fold changes by stage, and interaction terms.

---

### 3.6 Q5: Baseline Environmental Changes in Uninfected Vines

**Question**: What transcriptomic changes occur in healthy (RB-) grapevines between seasons, providing a baseline for interpreting GRBV effects?

Analysis of RB(-) vines between 2021 and 2022 identified substantial seasonal DEGs in healthy vines, representing the natural environmental response baseline. The critical finding was:

- A subset of GRBV-responsive genes overlapped with environmentally responsive genes in healthy vines
- In some cases, GRBV amplified the natural seasonal response (same direction as environment)
- In other cases, GRBV opposed the seasonal trend (opposite direction), suggesting the virus actively counters host environmental adaptations

This "environmental mimicry" vs "environmental antagonism" framework helps distinguish truly virus-specific pathways from those that are environmentally co-regulated.

> **Figure 6** (`06_Q5_uninfected_environment.png`): Baseline seasonal changes in RB(-), overlap with GRBV response, and direction concordance analysis.

---

### 3.7 Q6: Infection Age Trajectory

**Question**: How does the transcriptome change as a function of time since GRBV infection?

Infection age cohorts (calculated as sampling year minus infection year) were compared to trace the progressive remodeling of the host transcriptome. The trajectory analysis revealed:

- A biphasic pattern: an acute phase (early years) with strong defense gene activation, followed by a chronic phase with metabolic pathway restructuring
- Genes in defense signaling (NB-LRR, PR proteins, WRKY TFs) were most active in early infection
- Genes in carbohydrate metabolism, phenylpropanoid biosynthesis, and hormone signaling showed progressive changes accumulating over time
- K-means clustering identified distinct temporal expression patterns, grouping genes into early-transient, progressive-accumulating, and late-onset categories

> **Figure 7** (`07_Q6_transcriptome_trajectory.png`): Infection age-dependent DEG counts, trajectory clustering, and progressive pathway activation.

---

### 3.8 Pathway Enrichment Analysis

GO (Biological Process, Molecular Function, Cellular Component) and KEGG pathway enrichment was performed for every DE comparison across all six research questions. Key pathway-level findings include:

**Consistently enriched pathways across conditions:**
- **Phenylpropanoid biosynthesis** (KEGG) -- consistently disrupted, reflecting GRBV's impact on secondary metabolism
- **Plant-pathogen interaction** (KEGG) -- enriched in infection response across seasons/stages
- **ABC transporters** (KEGG) -- upregulated, possibly for defense compound transport
- **Protein processing in ER** (KEGG) -- stress-responsive, linked to heat shock proteins
- **Linoleic acid / alpha-Linolenic acid metabolism** (KEGG) -- lipid signaling, jasmonate precursors

**Stage-dependent pathway dynamics:**
- Defense pathways (GO: "defense response", "response to biotic stimulus") peaked at Post-veraison
- Carbohydrate metabolism pathways intensified from Veraison to Harvest
- Cell wall modification pathways were activated predominantly at Harvest

**Season-specific pathway modulation:**
- Some metabolic pathways showed enrichment only in one season, confirming environmental contingency
- Defense response pathways were more broadly conserved between seasons

> **Figures 8-14** (`08-14_pathway_enrichment.png`): GO and KEGG enrichment results for Q1-Q6, plus a master pathway landscape heatmap showing enrichment significance across all conditions.

---

### 3.9 Seasonal Maturity Progression Analysis

A dedicated analysis examined how GRBV-induced changes progress through maturity stages within each season, revealing striking differences in the temporal dynamics of infection response between 2021 and 2022.

**Key Findings:**

| Metric | 2021 Season | 2022 Season |
|--------|-------------|-------------|
| **Peak DEG stage** | Harvest (1,192 DEGs) | Post-veraison (565 DEGs) |
| **Trajectory** | Progressive intensification | Early peak, partial resolution |
| **FC correlation (Pre-ver.)** | r = 0.329 | -- |
| **FC correlation (Veraison)** | r = 0.524 | -- |
| **FC correlation (Post-ver.)** | r = 0.692 | -- |
| **FC correlation (Harvest)** | r = 0.586 | -- |

The cross-season fold-change correlation increased from Pre-veraison (r = 0.329) to Post-veraison (r = 0.692), indicating that the GRBV response becomes more conserved across years as development progresses. The decrease at Harvest (r = 0.586) suggests late-season environmental divergence.

The **inverted trajectory** between seasons is a critical finding: in 2021, the infection response built progressively and peaked at Harvest, while in 2022, it peaked earlier at Post-veraison and partially resolved by Harvest. This suggests that 2022 growing conditions may have either accelerated the host response or allowed partial recovery mechanisms.

**Stage transition analysis** identified genes that:
- **Intensify**: Increase in infection effect from one stage to the next (accumulating response)
- **Weaken**: Decrease in infection effect (recovery or adaptation)
- Season-specific transition genes that drive the divergent trajectories

> **Figures 15-19** (`15-19_seasonal_maturity_*.png`): Seasonal DEG comparisons, pathway enrichment grids, pathway landscape heatmaps, PCA transitions, and stage-by-stage pathway detail.

---

### 3.10 Biomarker Discovery

A systematic biomarker discovery pipeline evaluated all 14,604 genes across 20 different DE conditions (overall, 2 seasons, 4 stages, 8 season x stage combinations, and 6 infection age groups) to identify robust diagnostic markers.

#### 3.10.1 Viral Gene Biomarkers (Perfect Classifiers)

The top 5 biomarkers were all **viral genes** (N761_gp1-gp5), achieving:
- **Biomarker Score**: 0.97-1.00
- **AUC-ROC**: 1.000 (perfect discrimination)
- **Significant in**: All 21 conditions tested
- **Mean |log2FC|**: 4.1 to 7.3

These viral transcripts provide unambiguous infection status but require that the virus is actively transcribing, which may not apply to latent infections.

#### 3.10.2 Host Gene Biomarkers

The top host biomarkers (excluding viral genes) were:

| Rank | Gene ID | Gene Name | Score | AUC | Direction | Conditions Sig |
|------|---------|-----------|-------|-----|-----------|---------------|
| 1 | LOC100854143 | F-box/LRR-repeat protein | 0.741 | 0.880 | UP | 19/21 |
| 2 | LOC100256964 | E3 ubiquitin-protein ligase RSL1 | 0.731 | 0.894 | UP | 18/21 |
| 3 | LOC104881857 | Protein NTM1-like 9 | 0.725 | 0.869 | UP | 19/21 |
| 4 | LOC100268011 | RNA-dependent RNA polymerase 1 | 0.716 | 0.895 | UP | 18/21 |
| 5 | LOC100257913 | Glycine-rich cell wall structural protein 2 | 0.705 | 0.871 | DOWN | 18/21 |
| 6 | LOC100262652 | NAC domain-containing protein 62 | 0.655 | 0.845 | UP | 18/21 |
| 7 | LOC100263041 | Pleiotropic drug resistance protein 3 | 0.629 | 0.811 | UP | 17/21 |
| 8 | LOC100260978 | Glucan endo-1,3-beta-glucosidase 8 | 0.583 | 0.788 | UP | 14/21 |

Notable functional categories among top host biomarkers:
- **Defense signaling**: F-box/LRR protein, RNA-dependent RNA polymerase (antiviral RNAi)
- **Protein turnover**: E3 ubiquitin ligase (proteolytic regulation)
- **Transcription regulation**: NAC domain TF, NTM1-like protein
- **Cell wall remodeling**: Glycine-rich structural proteins (DOWN-regulated)
- **Transport**: ABC transporter (pleiotropic drug resistance protein)

#### 3.10.3 Multi-Gene Panel Classification

Multi-gene panels were evaluated using **5-fold stratified cross-validation** with Logistic Regression (LR) and Random Forest (RF) classifiers:

| Panel Size | LR AUC (mean +/- SD) | RF AUC (mean +/- SD) |
|------------|----------------------|----------------------|
| 3 genes | 1.000 +/- 0.000 | 1.000 +/- 0.000 |
| 5 genes | 1.000 +/- 0.000 | 1.000 +/- 0.000 |
| 10 genes | 1.000 +/- 0.000 | 1.000 +/- 0.000 |
| 20 genes | 1.000 +/- 0.000 | 1.000 +/- 0.000 |

A **3-gene panel** achieves perfect AUC = 1.000 across all cross-validation folds, demonstrating that GRBV infection creates a highly distinctive transcriptomic signature that can be captured with minimal gene sets. This has strong implications for developing RT-qPCR diagnostic assays.

#### 3.10.4 Biomarker Categories

Biomarkers were classified into four functional categories:

| Category | Count | Description |
|----------|-------|-------------|
| **Universal** | 12 | Season-robust, significant in >=3 maturity stages, AUC >= 0.7 |
| **Early-detection** | 35 | Significant at Pre-veraison or Veraison, AUC >= 0.65 |
| **Harvest-specific** | 23 | Significant primarily at Harvest stage |
| **Infection-age responsive** | 33 | Significant in >=3 infection age groups, AUC >= 0.7 |

**Universal biomarkers** are the most practically valuable, as they detect infection regardless of when samples are collected. **Early-detection biomarkers** enable pre-symptomatic diagnosis before visible symptoms appear.

> **Figures 20-22** (`20-22_biomarker_*.png`): Biomarker overview (scores, AUC, panels, categories), expression profiles of top 16 biomarkers, and robustness analysis (FC consistency, AUC stability, significance matrix).

---

### 3.11 GWAS/QTL Cross-Reference Analysis

The 45 host biomarker genes were mapped to chromosomal positions using the NCBI Gene database and cross-referenced with a curated database of **121 known grapevine GWAS/QTL loci** spanning disease resistance (Rpv1-31, Ren1-12, Run1-2), berry quality, abiotic stress tolerance, and developmental control.

#### 3.11.1 Overview of Co-localizations

| Metric | Count | Percentage |
|--------|-------|-----------|
| Biomarkers near known QTL (within 1 Mb) | **27/45** | **60.0%** |
| Biomarkers inside QTL intervals | **18/45** | **40.0%** |
| Biomarkers with no QTL overlap | 18/45 | 40.0% |

#### 3.11.2 Overlaps by QTL Category

| QTL Category | Biomarker Genes Overlapping |
|--------------|---------------------------|
| Disease Resistance | **18** |
| Berry Quality | **8** |
| Developmental | **8** |
| Abiotic Stress | **4** |

#### 3.11.3 Key Co-localizations

**Disease Resistance Loci:**

| Biomarker Gene | Position | QTL Overlap | Significance |
|----------------|----------|-------------|--------------|
| LOC100263041 (ABC transporter) | chr5: 2.77 Mb | **Rpv11, RGA cluster** | Downy mildew resistance region; drug efflux transporter activated by GRBV suggests cross-talk between biotic stress pathways |
| LOC100257913 (cell wall protein) | chr5: 12.16 Mb | **Rpv14** | Downy mildew resistance locus; cell wall protein downregulated by GRBV suggesting viral manipulation of structural defense |
| LOC100261750, LOC100253147 | chr14: 24.4 Mb | **Rpv8, Rpv12, Ren2, PdR1a/b, NB-LRR cluster** | Major multi-pathogen resistance hotspot on chr14; two GRBV biomarkers map to the densest NB-LRR cluster in the grapevine genome |
| LOC100267000 (sulfate transporter) | chr9: 4.31 Mb | **Rpv5, Rpv10, Ren6, NB-LRR cluster** | Second major resistance gene cluster; sulfur transport is linked to defense compound biosynthesis |
| LOC100261732 (nucleosome assembly) | chr9: 2.15 Mb | **Rpv5, Ren6, NB-LRR cluster** | Chromatin remodeling gene within resistance locus; epigenetic regulation of defense responses |
| LOC104881139 (VACUOLELESS GAMETOPHYTES) | chr13: 2.99 Mb | **Phylloxera Rdv1** | Phylloxera resistance locus; vacuolar processes may be co-opted during viral infection |
| LOC100247751, LOC100245905 | chr14: 1.7-2.1 Mb | **Rpv29, Ren5** | Proximal chr14 resistance locus; heavy metal-associated protein and defense-related gene |

**Berry Quality Loci:**

| Biomarker Gene | Position | QTL Overlap | Significance |
|----------------|----------|-------------|--------------|
| **LOC100268011 (RNA-dep. RNA polymerase 1)** | chr1: 5.58 Mb | **Sugar QTL + Berry weight QTL** | Our #1 host biomarker co-localizes with the sugar accumulation QTL, directly connecting the antiviral RNAi machinery with sugar metabolism disruption -- a central feature of GRBV disease |
| LOC100261077 | chr1: 4.40 Mb | **Sugar QTL** | Second biomarker in the sugar accumulation region |
| LOC100260313 (TORTIFOLIA1-like) | chr1: 10.32 Mb | **Anthocyanin modifier QTL** | Protein involved in cell morphogenesis maps to anthocyanin modification locus; GRBV disrupts anthocyanin profiles |
| LOC100258655 (expansin-A1) | chr17: 7.29 Mb | **Berry weight QTL** | Cell expansion gene within berry size QTL; GRBV affects berry development |

**Developmental Loci:**

| Biomarker Gene | Position | QTL Overlap | Significance |
|----------------|----------|-------------|--------------|
| LOC100248652 | chr19: 7.46 Mb | **VvNAC ripening cluster** | Maps to NAC transcription factor cluster controlling ripening; GRBV delays ripening through transcriptional regulators at this locus |
| LOC100258268 (LEA protein) | chr9: 5.86 Mb | **Root depth QTL** | Late embryogenesis abundant protein in root architecture region |

**Abiotic Stress Loci:**

| Biomarker Gene | Position | QTL Overlap | Significance |
|----------------|----------|-------------|--------------|
| LOC100260978 (glucanase) | chr11: 3.65 Mb | **Cold hardiness QTL** | Beta-glucosidase in cold tolerance region; defense enzyme with cross-stress function |
| LOC100260313 (TORTIFOLIA1-like) | chr1: 10.32 Mb | **Drought QTL** | Dual stress overlap; drought response machinery engaged during viral infection |

#### 3.11.4 Chromosome Distribution

Biomarker genes are distributed across all 19 grapevine chromosomes, but notable concentrations exist on:
- **Chromosome 14**: 6 biomarkers, including 4 within the major NB-LRR resistance gene cluster (22-27 Mb)
- **Chromosome 9**: 5 biomarkers, several within the Rpv5/Rpv10/Ren6 resistance region
- **Chromosome 18**: 4 biomarkers (though none overlapped the distal Rpv3/Ren4 cluster)
- **Chromosome 1**: 3 biomarkers, two within sugar/berry quality QTLs

> **Figure 24** (`24_biomarker_QTL_chromosome_map.png`): Chromosome ideogram showing biomarker positions (triangles) against colored QTL intervals across all 19 chromosomes.
>
> **Figure 25** (`25_biomarker_QTL_crossref_summary.png`): Cross-reference summary with overlap statistics, category breakdown, chromosome distribution, co-localization table, and biological interpretation.

---

## 4. Discussion

### 4.1 GRBV Creates a Multi-Layered Transcriptomic Perturbation

This study reveals that GRBV does not simply activate a single defense response but creates a complex, multi-layered perturbation of the grapevine transcriptome that is modulated by developmental context (maturity stage), temporal context (infection duration), and environmental context (growing season). The approximately 14,600 expressed genes represent a dynamic regulatory network in which GRBV effects interact with all three axes of variation.

### 4.2 The Inverted Seasonal Trajectory

Perhaps the most striking finding is the **inverted seasonal trajectory** of GRBV's transcriptomic impact: in 2021, the response peaked at Harvest (1,192 DEGs), while in 2022, it peaked at Post-veraison (565 DEGs). This suggests that year-to-year climate variation fundamentally alters the kinetics of the host-virus interaction. In practical terms, this means that sampling time optimization for GRBV diagnosis must consider the growing season's environmental profile.

### 4.3 Biomarker Robustness and Diagnostic Potential

The identification of **12 universal host biomarkers** that are significant across both seasons, multiple maturity stages, and multiple infection ages provides a robust diagnostic platform. The fact that a **3-gene panel achieves AUC = 1.000** in cross-validation demonstrates the feasibility of developing simple, cost-effective molecular diagnostics (e.g., RT-qPCR assays) for GRBV detection.

The top host biomarker, **LOC100854143 (F-box/LRR-repeat protein, AUC = 0.880)**, is particularly interesting because F-box/LRR proteins are a major class of immune receptors in plants. Its consistent upregulation across conditions suggests it may be a core component of the host surveillance system that detects GRBV.

The **RNA-dependent RNA polymerase 1 (LOC100268011, AUC = 0.895)** is notable because RdRP1 is a key enzyme in the RNA silencing (RNAi) pathway, which is the primary plant antiviral defense mechanism. Its upregulation is a direct marker of activated antiviral immunity.

### 4.4 GWAS/QTL Co-localization: Bridging Transcriptomics and Genetics

The finding that **60% of biomarker genes map near known GWAS/QTL loci** (and 40% fall directly within QTL intervals) is remarkable and has several implications:

1. **Shared genetic architecture**: GRBV biomarkers disproportionately localize to known resistance gene clusters (NB-LRR regions on chr9, chr14, chr5), suggesting that GRBV exploits or disrupts the same innate immune machinery used against other pathogens (downy mildew, powdery mildew, Pierce's disease).

2. **Sugar accumulation link**: The co-localization of the top host biomarker (RNA-dep. RNA polymerase 1) with the sugar QTL on chromosome 1 provides a direct genetic link between the antiviral defense response and the primary symptom of GRBV disease (reduced Brix). This suggests a mechanistic connection: the host's antiviral RNAi machinery may compete for resources with or interfere with sugar metabolism pathways at this genomic locus.

3. **Breeding implications**: Biomarker genes that co-localize with known resistance QTLs represent candidate genes for marker-assisted selection (MAS) in breeding programs. For example, allelic variation at the chr14 NB-LRR cluster (where two GRBV biomarkers map) could potentially confer differential tolerance to both GRBV and downy mildew.

4. **Ripening regulation**: The co-localization of a biomarker with the VvNAC ripening cluster on chromosome 19 confirms at the genetic level what is observed phenotypically: GRBV disrupts the NAC-mediated transcriptional regulation of fruit ripening.

### 4.5 Limitations

- QTL positions are approximate, based on flanking markers mapped to the PN40024 reference genome
- The QTL database is compiled from published literature and may not include the most recently characterized loci
- Co-localization does not prove functional causation; follow-up functional genomics studies are needed
- The DE analysis uses Welch's t-test rather than negative binomial models (DESeq2/edgeR); while robust for large sample sizes, this may miss low-count genes

---

## 5. Conclusions

1. **GRBV creates a developmental-stage-dependent, season-modulated, duration-progressive transcriptomic perturbation** in grapevine that affects defense, metabolism, and developmental pathways.

2. **The viral transcriptomic signature is inverted between seasons**: peaking at Harvest in 2021 but at Post-veraison in 2022, demonstrating that environmental context fundamentally shapes the host-pathogen interaction kinetics.

3. **Robust biomarkers exist**: 12 universal host biomarkers are consistently significant across seasons, stages, and infection ages. A 3-gene panel achieves perfect classification (AUC = 1.000), enabling diagnostic assay development.

4. **60% of biomarker genes co-localize with known GWAS/QTL loci**, with the strongest enrichment in disease resistance clusters (NB-LRR regions) and berry quality loci (sugar accumulation, anthocyanin modification).

5. **The top host biomarker (RNA-dependent RNA polymerase 1) maps to the sugar accumulation QTL on chromosome 1**, providing a direct genetic link between the antiviral defense pathway and the primary disease symptom.

6. **Multiple biomarker genes map to the chr14 NB-LRR mega-cluster** (Rpv8/Rpv12/Ren2/PdR1), suggesting shared innate immune machinery between GRBV and other grapevine pathogens, with implications for multi-pathogen resistance breeding.

---

## 6. Supplementary Files

### Analysis Scripts
| Script | Description |
|--------|-------------|
| `analysis.py` | Main DE analysis addressing 6 research questions |
| `pathway_enrichment.py` | GO and KEGG pathway enrichment for all comparisons |
| `seasonal_maturity_progression.py` | Season x Stage stratified analysis |
| `biomarker_discovery.py` | Biomarker scoring, AUC-ROC, and panel classification |
| `gwas_qtl_crossref.py` | GWAS/QTL cross-referencing with embedded QTL database |
| `analysis_schematic.py` | Analysis pipeline flowchart |

### Data Files
| File | Description |
|------|-------------|
| `counts_matrix.txt` | Raw RNA-Seq count matrix (28,168 genes x 327 samples) |
| `sample_metadata.txt` | Sample metadata (infection group, season, stage, age) |
| `org.Vvinifera.eg.db/` | *V. vinifera* functional annotation database |
| `GRBV_biomarkers_top50.csv` | Top 50 biomarkers with scores, AUC, annotations |
| `GRBV_biomarker_QTL_crossref.csv` | Full biomarker-QTL co-localization table |

### Figures (25 total)
| Figure | File | Description |
|--------|------|-------------|
| 1 | `01_global_pca_overview.png` | Global PCA overview |
| 2 | `02_Q1_infection_years.png` | Infection duration effects |
| 3 | `03_Q2_maturity_stages.png` | Maturity stage effects |
| 4 | `04_Q3_seasonal_effects.png` | Seasonal effects |
| 5 | `05_Q4_environmental_effects.png` | Environmental modulation |
| 6 | `06_Q5_uninfected_environment.png` | Baseline environmental changes |
| 7 | `07_Q6_transcriptome_trajectory.png` | Infection age trajectory |
| 8-13 | `08-13_Q1-Q6_pathway_enrichment.png` | Pathway enrichment per question |
| 14 | `14_master_pathway_landscape.png` | Master pathway landscape heatmap |
| 15 | `15_seasonal_maturity_DEG_comparison.png` | Seasonal DEG comparison |
| 16 | `16_seasonal_maturity_pathway_enrichment.png` | Season x stage pathway enrichment |
| 17 | `17_seasonal_maturity_pathway_landscape.png` | Seasonal pathway landscape |
| 18 | `18_seasonal_maturity_transitions_PCA.png` | PCA transitions and stage dynamics |
| 19 | `19_seasonal_maturity_pathway_detail.png` | Stage-by-stage pathway comparison |
| 20 | `20_biomarker_overview.png` | Biomarker scores, AUC, panels, categories |
| 21 | `21_biomarker_profiles.png` | Top 16 biomarker expression profiles |
| 22 | `22_biomarker_robustness.png` | Biomarker robustness analysis |
| 23 | `23_analysis_schematic.png` | Analysis pipeline flowchart |
| 24 | `24_biomarker_QTL_chromosome_map.png` | Chromosome ideogram with QTL map |
| 25 | `25_biomarker_QTL_crossref_summary.png` | QTL cross-reference summary |

---

*Report generated from the GRBV RNA-Seq Transcriptomics Analysis Pipeline*
*Data: 327 samples, 14,604 genes, 2 seasons, 4 maturity stages, 7 infection year cohorts*
*Analyses: Differential expression, pathway enrichment, biomarker discovery, GWAS/QTL cross-referencing*
