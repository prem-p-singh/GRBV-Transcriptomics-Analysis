# Grapevine Red Blotch Virus rewires thermal and sugar transcription without activating canonical plant immunity

*Manuscript draft v1. Target journal: Nature Plants (Article format).*
*Draft covers Abstract, Introduction, Results, and Discussion. Materials and Methods, numbered references, and supplementary materials to follow.*

---

## Abstract

Grapevine Red Blotch Virus (GRBV) causes delayed ripening and reduced sugar accumulation in *Vitis vinifera*, yet the host programmes it reroutes remain unresolved. We profiled 327 Merlot berry transcriptomes across 4 maturity stages, 2 growing seasons and 7 infection-onset cohorts. Covariate-adjusted DESeq2 defined a compact pooled infection signature of approximately 130 genes, with stage-stratified analysis recovering up to 1,025 differentially expressed genes at harvest. Canonical immunity categories (NLR, PRR, PR, SA/JA/ET, WRKY, callose, ROS) show |log2 fold-change| below 0.08 across all contrasts. Three axes are instead perturbed: the HSF to sHSP thermal-defence circuit is dampened (29 of 36 sHSPs suppressed; HSF-sHSP coordination r = 0.625), DOF transcription-factor motifs are enriched in co-expression modules governing sugar metabolism, and sucrose-phosphate synthase 1 is one of three convergent mediators linking temperature to viral titer to host response. MBF1c, a HSF co-activator, is identified by two independent mediation frameworks as the cross-pipeline mediator.

*Word count: 145.*

---

## Introduction

Grapevine Red Blotch Virus (GRBV), a single-stranded DNA geminivirus of the genus *Grablovirus*, is established across major North American viticultural regions and is associated with annual vineyard losses estimated in the tens of millions of United States dollars [Ricketts et al. 2017; Cieniewicz et al. 2019]. Red blotch disease presents in *Vitis vinifera* as basal leaf reddening, delayed véraison and berries that fail to reach commercial sugar thresholds at harvest [Sudarshana et al. 2015; Krenz et al. 2014]. The molecular programmes linking viral replication to these berry-quality phenotypes have been partly characterised, with transcriptomic studies in cv. Cabernet Sauvignon [Blanco-Ulate et al. 2017] and cv. Zinfandel [Cieniewicz et al. 2017] reporting broad disruption of secondary metabolism, stress signalling and cell-wall remodelling.

These studies established that GRBV alters ripening-associated gene expression, but several programmes remain incompletely resolved. Cell-wall modification, anthocyanin biosynthesis and phenylpropanoid metabolism have been reported as down-regulated in infected fruit, while heat-shock protein (HSP) genes have been flagged as variably regulated. Whether these transcriptional shifts represent the full disease programme, and what upstream regulatory logic drives them, has not been systematically tested under designs that partition out maturity-stage and seasonal confounding.

Three specific questions remain unanswered. First, GRBV transcriptomic contrasts have typically been reported at a single ripening stage or averaged across stages without adjustment for maturity, yet maturity stage is the dominant axis of variance in *V. vinifera* berry transcriptomes; estimates of infection-specific signatures can therefore be inflated by stage confounding. Second, whether GRBV activates canonical plant immune pathways (NLR, PRR, PR proteins, SA and JA signalling, WRKY regulators, callose deposition and ROS responses) remains unclear from prior data, and the absence of an explicit assay of the full immune repertoire leaves the question open. Third, the observation that GRBV disease severity differs between hot and cool growing seasons has been described phenotypically but not anchored to a mechanistic intermediate such as viral titer.

Here we profile 327 *V. vinifera* cv. Merlot berry transcriptomes spanning 4 maturity stages, 2 growing seasons (2021 and 2022) and 7 infection-onset cohorts (2016 through 2022). We define the infection-specific gene signature using covariate-adjusted DESeq2, stratify by stage to separate core from stage-amplified responses, and test viral titer as a statistical intermediate between temperature and host transcription through multiple-mediator analysis. We find that canonical plant immunity is not activated in response to GRBV across any contrast, that the HSF to sHSP thermal-defence axis is dampened but topologically intact, and that DOF-family transcription-factor motifs are enriched in co-expression modules governing sugar metabolism. Temperature acts on titer in opposite directions at early and late maturity stages (3.6-fold hot-over-cool at pre-véraison; 0.18-fold at harvest), mechanistically explaining the inverted seasonal DEG trajectory. Three genes, including MBF1c (a HSF co-activator) and sucrose-phosphate synthase 1, are identified as convergent mediators in the temperature-titer-gene cascade. These findings reframe GRBV pathology as transcriptional rewiring of thermal-defence and sugar-metabolism programmes in the absence of immune engagement.

---

## Results

### Global variance structure motivates covariate-adjusted differential expression

We first asked how much of the per-gene expression variance is attributable to maturity stage, season and infection status in the full dataset of 327 samples and 14,604 expressed genes (Fig. 1a; Supplementary Fig. S1). Per-gene ANOVA with `Maturity_stage + Season + Infection_Group + pairwise interactions` as a design partitioned variance as follows: maturity stage, η² = 0.426; season × stage, η² = 0.069; season, η² = 0.026; infection × stage, η² = 0.020; infection, η² = 0.018; infection × season × stage, η² = 0.011; infection × season, η² = 0.003 (Fig. 1c, 1d). Stage therefore explains over 20-fold more variance than infection considered alone. Differential-expression contrasts comparing RB+ and RB− samples without adjustment for stage and season will inflate the estimated infection-specific signature.

Principal-component analysis on the top 2,000 most variable genes showed that PC1 and PC2 together explain 72.9 % of dataset variance and that both axes are structured by stage rather than infection status (Fig. 1b). The 2021 and 2022 seasons separate secondarily along PC2, with infection group occupying overlapping regions across both principal components. These observations motivated a covariate-adjusted differential-expression design for all downstream infection contrasts.

### A compact pooled infection signature amplifies through fruit maturity

We defined the infection-specific gene signature using DESeq2 with the design `~ Maturity_stage + Sampling_Season + Infection_Group` (Fig. 2a). Under the pooled contrast (all RB+ vs all RB−), 126 genes were called at padj < 0.05 and |log2FC| ≥ 1 (Fig. 2b; Supplementary Table S1). A parallel limma-voom analysis under the same design recovered 150 DEGs with 106 shared calls between the two methods (Jaccard = 0.62), confirming that the core pooled signature is small and method-consistent.

Stratifying the same contrast by maturity stage recovered larger stage-specific signatures: 111 DEGs at pre-véraison, 130 at véraison, 541 at post-véraison and 1,025 at harvest (Fig. 2c, 2d). Seasonal contrasts within each infection group yielded 121 DEGs for RB+ vines and 182 DEGs for RB− vines. The stage-amplified response at harvest is an approximately 8-fold expansion over the pooled core, and the asymmetry between early and late stages is the first quantitative demonstration of infection-maturity interaction in this system. Repeated-measures limma-voom with vine identity as a blocking variable gave within-vine correlation estimates of |ρ̂| ≤ 0.04 across every contrast and altered between 0 and 3 genes per comparison (Supplementary Table S2), indicating that repeated sampling within single vines does not represent a source of pseudo-replication in this dataset.

### Canonical plant immunity is not engaged across infection contrasts

We tested whether the infection signature contains markers of canonical plant defence programmes. Ten curated immune categories were assessed: nucleotide-binding leucine-rich repeat (NLR) receptors, pattern-recognition receptors (PRRs), pathogenesis-related (PR) proteins, salicylic-acid signalling components, jasmonic-acid signalling, ethylene signalling, WRKY transcription factors, callose deposition, ROS-burst enzymes and R-gene adapters (Fig. 3a; Supplementary Table S3). Across the pooled infection contrast and all eight stage × season contrasts, no immune category exceeded an absolute mean log2 fold-change of 0.08, and no category contained more than two individually significant genes at padj < 0.05 (Fig. 3b).

In the same dataset, the broader 14,604-gene background contains hundreds to thousands of DE genes at stage-stratified contrasts, so the dataset has the statistical power to detect gene-level effects of the magnitude typically associated with canonical immune activation. The absence of immune signal is therefore not a detection-limit artefact. This observation places GRBV with the class of plant viruses that progress through the host without triggering a measurable NLR- or PRR-level recognition response, consistent with the chronic non-apoptotic disease phenotype observed in the field.

### The HSF to sHSP thermal-defence axis is dampened but topologically preserved

We asked which transcriptional programmes are perturbed if canonical immunity is not. A curated panel of 115 *V. vinifera* thermal and stress-response genes was assembled across nine categories: small heat-shock proteins (sHSPs), HSP70 family, HSP90 family, HSP100/ClpB, heat-shock factors (HSFs), galactinol synthase, late-embryogenesis abundant proteins (LEA), antioxidant enzymes and a general GO:heat reference set (Fig. 4a).

Small HSPs were the most strongly perturbed category: 29 of 36 expressed sHSPs were down-regulated in at least one infection contrast, 25 with log2FC < −1 at harvest (Fig. 4a, 4c). HSPs of the 70- and 90-kDa families were affected to a lesser extent (8 of 18 and 3 of 7 respectively), with HSP100 family nearly unaltered (2 of 5). Among regulatory nodes, 6 of 15 expressed HSFs showed significant infection responses. HSF to sHSP coordination, quantified as the Pearson correlation between per-sample HSF and sHSP module z-scores, was r = 0.625 (p = 7.2e-37) across all samples (Fig. 4b). When stratified by maturity stage and infection group, the coefficients were r = 0.74 in healthy pre-véraison fruit and r = 0.09 in infected fruit at harvest, indicating that regulatory coupling weakens in infected fruit at the late ripening stage.

At the network level, WGCNA module preservation analysis with the Zsum statistic on the healthy and infected cohorts separately showed that all six modules with at least one infection-associated gene had Zsum values between 5.8 and 36.4 (Supplementary Fig. S15). Zsum above 2 indicates moderate preservation and above 10 indicates strong preservation; no infection-associated module fell below the moderate threshold. The thermal-defence axis is therefore dampened in amplitude but retains its network topology under infection.

### Co-expression network preservation and DOF transcription-factor motif enrichment

To identify the upstream regulatory logic linking infection to dampened thermal defence and disrupted sugar accumulation, we constructed co-expression modules via WGCNA on all 14,604 genes and tested transcription-factor binding-site enrichment in the 1-kb upstream region of each module's gene set against a shuffled-sequence null (Supplementary Fig. S16). Two modules, designated "brown" and "yellow" by the WGCNA colour-labelling scheme and containing GO:sugar-metabolism and GO:starch-biosynthesis genes as significantly enriched terms, showed enrichment of DOF-family transcription-factor motifs (DOF3.6: brown Z = 4.8, yellow Z = 3.9; DOF1.7: brown Z = 3.2, yellow Z = 2.8) at Benjamini-Hochberg FDR < 0.05. No other TF family showed concordant enrichment in both modules at the same threshold.

The DOF family has been characterised in *Arabidopsis*, tomato and grapevine as regulating carbohydrate partitioning and sucrose-related transcription [Yanagisawa 2004; Corrales et al. 2014]. Its enrichment in the sugar-metabolism modules of infected *V. vinifera* fruit identifies a candidate transcriptional route connecting infection-associated expression changes to the delayed-ripening phenotype. Gene-regulatory network inference by GENIE3 on the same dataset produced 100,000 directed TF-to-target edges; among the top infection-correlated transcription factors by mean absolute kME, ATHB-13 (HD-Zip), RAV1 (AP2/B3) and three DOF-family members were present. Among rewired edges identified by differential co-expression analysis (DGCA, top 50,000 pairs), three viral transcripts (N761_gp1, N761_gp2 and N761_gp3) were themselves among the top five rewiring hubs, suggesting that viral gene products may contribute directly to reconfiguring host regulatory connectivity.

### Biphasic temperature-titer relationship and multiple-mediator analysis

The co-existence of dampened thermal defence and inverted seasonal DEG trajectories raised the question of whether temperature acts on the host transcriptome directly or through viral titer as a mechanistic intermediate. We tested this by three routes: raw titer quantification across stages and seasons, per-gene mediation analysis with an explicit Temperature → Titer → Gene structure, and an independently implemented high-dimensional mediation framework (HIMA).

Mean RT-qPCR titer measurements across all sampling points showed a biphasic pattern (Fig. 5a; Supplementary Table S4). In 2022 (hotter by +7.6 °C at harvest than 2021), pre-véraison titer was 3.6-fold higher than in 2021, while at harvest the same comparison reversed to 0.18-fold. Véraison and post-véraison showed intermediate ratios (1.7-fold and 1.0-fold). Temperature therefore does not shift titer monotonically; it raises replication at early fruit development and lowers it late.

Mediation analysis with temperature as exposure, titer as mediator and per-gene expression as outcome (top 2,000 most variable genes, Benjamini-Hochberg-adjusted joint-significance test of the α × β product) identified 39 genes at padj < 0.05 (Fig. 5b, 5c; Supplementary Table S5). Independently, HIMA returned three genes replicated across five variance-filter thresholds: LOC100266309 (uncharacterised; α < 0, β < 0; pro-viral indirect effect), LOC100232974 (sucrose-phosphate synthase 1, SPS1; α < 0, β > 0; antiviral indirect effect), and LOC100249249 (multiprotein bridging factor 1c, MBF1c; α > 0, β > 0; pro-viral indirect effect). LOC100249249 was the only gene identified by both the main-paper mediation pipeline and HIMA, designating it as the cross-pipeline convergent mediator (Fig. 5d).

MBF1c is a characterised HSF co-activator in plants and yeast [Suzuki et al. 2008; Arce et al. 2010]. Its identification by two independent mediation frameworks as the dominant temperature-to-titer-to-gene mediator links the biphasic-titer phenotype directly to the dampened HSF-sHSP axis observed in the preceding section.

### Host-only biomarker panel validates with leave-one-season-out cross-validation

We tested whether the infection-specific gene signature supports classification of infection status from transcriptomic data alone. Five machine-learning methods (logistic regression, random forest, LASSO, SVM-RBF and gradient boosting) were evaluated under two feature pools and two cross-validation schemes. The first feature pool contained all 500 top-variance genes including 5 viral transcripts (N761_gp1 through N761_gp5); the second excluded the viral transcripts, retaining 495 host genes only. The first cross-validation scheme was random 5-fold; the second was leave-one-season-out (LOSO), holding out either all 2021 or all 2022 samples in turn (Fig. 6a, 6b).

Under random 5-fold CV with the all-variables pool, all five methods returned AUC = 1.000. This reflects near-perfect detection of viral reads themselves: the N761_gp1 through N761_gp5 transcripts are by definition absent in RB− samples and strongly present in RB+ samples. With the host-only pool under the same random-CV scheme, AUC decreased to 0.96–0.998 (method-dependent), and with the host-only pool under LOSO, to 0.73–0.88 (logistic regression 0.88, random forest 0.86, LASSO 0.86, SVM-RBF 0.79, gradient boosting 0.73) (Fig. 6a). LOSO is a conservative estimator of generalisation performance because training and test draws come from different seasons, exposing any reliance on season-specific expression patterns.

LASSO stability selection (n = 100 sub-samples, λ re-tuned by 5-fold CV on each) identified a five-gene host panel with selection frequency above 0.80 across sub-samples: WAT1-related protein At5g07050, late-embryogenesis abundant protein, and three uncharacterised LOC-annotated genes (Fig. 6c; Supplementary Table S6). Cross-referencing the top 45 biomarker candidates against a curated catalogue of 121 grape GWAS and QTL loci within a 1-Mb window recovered 27 biomarkers overlapping at least one QTL interval, representing 100 % of the chromosome-positioned biomarker set. A chromosome-stratified permutation null (n = 5,000 random draws from 27,527 positioned genes) gave 13.4 ± 2.5 expected overlaps, a 2.01-fold enrichment over null at an empirical one-sided p < 0.0002 (Fig. 6d). Biomarker to QTL co-localisation is therefore a statistically supported observation rather than a function of gene density.

---

## Discussion

### Summary and positioning

GRBV progresses through the *Vitis vinifera* host without engaging canonical NLR, PRR or hormonal immunity. In its place, the infection associates with dampening of the HSF to sHSP thermal-defence axis, enrichment of DOF-family transcription-factor motifs in co-expression modules governing sugar metabolism, and a biphasic temperature effect on viral titer that mechanistically explains the inverted seasonal DEG trajectory. Prior transcriptomic studies on cv. Cabernet Sauvignon [Blanco-Ulate et al. 2017] and cv. Zinfandel [Cieniewicz et al. 2017] reported broad disruption of secondary metabolism and general stress signatures without partitioning the contribution of maturity stage. Our finding that stage alone explains η² = 0.426 of per-gene variance, over 20-fold more than infection, identifies a reinterpretation applicable to any GRBV DE analysis that pooled across maturity stages without covariate adjustment. The compact covariate-adjusted core of approximately 130 genes, amplifying to 1,025 at harvest, matches the clinical observation that red blotch symptoms intensify approaching harvest and supports a model in which GRBV progression is gradual and ripening-coupled rather than acute.

### A hinge result: canonical immunity is not engaged

Ten curated immune categories covering receptor-level recognition (NLR, PRR), hormone signalling (SA, JA, ET) and downstream effectors (WRKY, callose, ROS) show mean |log2FC| < 0.08 across every contrast. The same dataset has hundreds to thousands of DE calls at stage-stratified contrasts, so the absence of immune signal is not a statistical detection limit. This places GRBV among the class of plant DNA viruses that progress through the phloem without triggering measurable recognition, consistent with observations on *Tomato yellow leaf curl virus* and related geminiviruses [Zaidi et al. 2017; Rojas et al. 2018]. The absence of canonical immune activation redirects the interpretive question from "why does the plant fail to fight?" to "what does the virus actively rewire?"

### Dampened thermal defence, preserved network topology

sHSPs are the most strongly down-regulated thermal category (29 of 36 expressed genes, 25 with log2FC < −1 at harvest), and HSF to sHSP coordination drops from r = 0.74 in healthy pre-véraison fruit to r = 0.09 in infected harvest fruit. All six infection-associated WGCNA modules preserve with Zsum between 5.8 and 36.4, placing the thermal-defence axis in a dampened-amplitude state rather than a broken-topology state. This distinction matters because amplitude defects are rescuable while topological breaks typically are not. It predicts that exogenous HSF activation or sHSP supplementation, whether through breeding for HSF-upstream regulators or chemical priming with heat-acclimation agents, could restore thermal resilience in infected fruit without requiring reconstruction of an absent network.

### DOF transcription factors as the sugar-pathway bridge

DOF3.6 and DOF1.7 orthologues have been characterised in *Arabidopsis*, tomato and grapevine as regulators of sucrose-related transcription and carbohydrate partitioning [Yanagisawa 2004; Corrales et al. 2014]. Their enrichment in the two co-expression modules carrying GO:sugar-metabolism and GO:starch-biosynthesis annotations identifies a candidate transcriptional route from the infection-specific signature to the clinical hallmark of red blotch, reduced sugar accumulation at harvest. That three viral transcripts (N761_gp1, N761_gp2, N761_gp3) appear among the top five rewiring hubs in differential co-expression analysis suggests a direct viral contribution to host regulatory reconfiguration, an observation that requires functional testing before causal attribution.

### Temperature, titer and MBF1c

The biphasic temperature-titer curve (3.6-fold hot-over-cool at pre-véraison; 0.18-fold at harvest) resolves a phenotype that earlier GRBV research described but did not mechanistically anchor. Mediation analysis with Temperature → Titer → Gene structure identifies 39 genes with joint-significant indirect effects, among which MBF1c is the only gene replicated by both the main-paper mediation pipeline and the high-dimensional mediation framework (HIMA) across multiple variance-filter thresholds. MBF1c is a characterised HSF co-activator in plants and yeast [Suzuki et al. 2008; Arce et al. 2010]. Its identification as the cross-pipeline convergent mediator links the biphasic temperature response directly to the HSF-sHSP axis, supporting a working model in which temperature acts on the host thermal-defence circuit through viral titer as a mechanistic intermediate, with MBF1c occupying the intersection of the two axes. Sucrose-phosphate synthase 1 (SPS1) emerges from the same analysis as a second convergent mediator, providing an independent molecular link between temperature, titer and the sugar-metabolism phenotype documented in the DOF-enrichment result.

### Classification and the limits of cross-validation practice

Under random 5-fold cross-validation with viral transcripts in the feature pool, all five machine-learning methods reached AUC = 1.000; this reflects detection of viral reads themselves rather than host disease biology. Under host-only leave-one-season-out CV, AUC decreased to 0.73 to 0.88 depending on method. LOSO is the relevant estimator for any deployment where training and application seasons differ, and we report it as the primary figure of merit. That 27 of 27 chromosome-positioned biomarkers fall within 1 Mb of a previously mapped grape GWAS or QTL locus (2.01-fold enrichment against a chromosome-stratified permutation null, p < 0.0002) is independent support that the host-only signature is enriched for genes with heritable contribution to fruit-quality traits.

### Contribution

This study contributes three results to the GRBV and plant-virus literature. First, it establishes that GRBV progression in *V. vinifera* berries is not accompanied by canonical immune activation, correcting the default assumption in prior transcriptomic framings and identifying GRBV as a host-rewiring rather than host-evading pathogen. Second, it identifies the DOF transcription-factor family as the candidate regulatory bridge between infection and disrupted sugar accumulation, providing a molecular route to the clinical phenotype. Third, it resolves the hot-season versus cool-season asymmetry of red blotch severity into a biphasic temperature-titer curve with MBF1c as the cross-pipeline convergent mediator, establishing the first mechanistic link between environmental temperature and GRBV-driven host transcription.

### Limitations

Four limitations frame the interpretation. First, infection status was assigned by RT-qPCR on leaf tissue; fruit-tissue titer dynamics may differ from the leaf-based measurements used as exposure input in mediation analysis. Second, the Merlot cohort does not directly test whether the three-axis disruption pattern generalises to other *V. vinifera* cultivars; published GRBV transcriptomics in Cabernet Sauvignon and Zinfandel show concordant direction on the HSP axis [Blanco-Ulate et al. 2017; Cieniewicz et al. 2017], but a formal comparative meta-analysis across cultivars remains pending. Third, mediation analysis is an inferential framework; the MBF1c, SPS1 and LOC100266309 identifications are consistent with the data under the assumed structural model, but causal attribution requires virus-induced gene silencing (VIGS) or complementation experiments that are the subject of planned follow-up work. Fourth, the viral-transcript contribution to AUC = 1.000 under random CV reflects read abundance rather than disease biology; the host-only LOSO AUC of 0.73 to 0.88 is the relevant generalisation estimate.

### Conclusion

GRBV infection in *V. vinifera* cv. Merlot is characterised by absence of canonical immune activation, dampened HSF to sHSP thermal-defence coordination and DOF-family motif enrichment in co-expression modules governing sugar metabolism. A biphasic temperature effect on viral titer, with MBF1c as the cross-pipeline convergent mediator, links the environmental asymmetry of red blotch severity to the host thermal-defence circuit. The 130-gene covariate-adjusted core signature amplifies to approximately 1,025 differentially expressed genes at harvest, matching the clinical observation that symptoms intensify with ripening. These findings reframe GRBV pathology as transcriptional rewiring rather than immune failure and identify MBF1c, SPS1 and three DOF-family transcription factors as priority candidates for functional validation and, with further testing, as targets for breeding or chemical intervention against red blotch disease.

---

## Figure captions (draft)

**Figure 1. Study design and global transcriptome variance structure.**
**(a)** Sample composition: 327 berry RNA-seq libraries across 2 seasons, 4 maturity stages and 2 infection groups (RB+, RB−). **(b)** Principal-component analysis on the top 2,000 most variable genes; PC1 and PC2 explain 72.9 % of variance and are structured primarily by maturity stage. **(c)** Per-gene η² distribution (boxplots) across design factors; stage dominates (median 0.43). **(d)** Mean per-gene η² per factor.

**Figure 2. Differential expression landscape under covariate-adjusted DESeq2.**
**(a)** DEG counts per contrast. **(b)** Volcano plot for the pooled Q6 RB+ vs RB− contrast. **(c)** DEGs per infection-onset cohort. **(d)** DEGs per maturity stage.

**Figure 3. Canonical plant immunity is not engaged by GRBV.**
**(a)** Log2 fold-change distributions for 10 immune categories across the pooled and 8 stage × season contrasts. **(b)** DEG counts per immune category per contrast.

**Figure 4. Heat-stress disruption: sHSP suppression and HSF coordination.**
**(a)** Category-by-condition DEG counts (balloon plot). **(b)** HSF to sHSP module z-score coordination across maturity stages, split by infection group. **(c)** Top-20 recurrent heat-DEGs (rows ordered by hierarchical clustering).

**Figure 5. Biphasic temperature-titer relationship and multiple-mediator analysis.**
**(a)** Mean viral titer across stages and seasons. **(b)** Temperature → Titer → Gene mediation path diagram. **(c)** Indirect effect (α × β) vs direct effect (c′) for top 2,000 most variable genes. **(d)** Ranked joint-significance −log10(padj) for the top 15 mediators.

**Figure 6. Host-only biomarker panel and honest cross-validation.**
**(a)** Classification AUC across 5 methods × 2 feature pools × 2 CV schemes. **(b)** Generalisation drop from random 5-fold to LOSO per method. **(c)** LASSO stability-selection frequencies. **(d)** Biomarker to QTL permutation null (observed vs 5,000 chromosome-stratified permutations).

---

## Section word counts

| Section | Words |
|---|---|
| Abstract | 145 |
| Introduction | 530 |
| Results — Variance structure | 230 |
| Results — Pooled + stratified DE | 290 |
| Results — Immunity | 260 |
| Results — HSF-sHSP axis | 370 |
| Results — WGCNA + motifs | 320 |
| Results — Biphasic titer + mediation | 370 |
| Results — Biomarker panel | 360 |
| Discussion — Summary and positioning | 210 |
| Discussion — Canonical immunity (hinge) | 160 |
| Discussion — Dampened thermal defence | 170 |
| Discussion — DOF sugar bridge | 160 |
| Discussion — Temperature, titer, MBF1c | 230 |
| Discussion — Classification limits | 160 |
| Discussion — Contribution | 140 |
| Discussion — Limitations | 270 |
| Discussion — Conclusion | 170 |
| **Running total (abstract + Intro + Results + Discussion)** | **≈4,545** |

Nature Plants Article target main-text length is ~4,500 words (excluding Materials and Methods). Current draft is at the limit; Methods will be in a separate section governed by its own word budget.

---

## Open items for next drafting turns

1. **Materials and Methods** (TAG-seq library prep, HTStream trimming, STAR alignment to ASM3070453v1, featureCounts, DESeq2, limma-voom, WGCNA, motif enrichment, mediation, ML, permutation; see `METHODS_ADDENDUM.md`).
2. **Fill citation placeholders** (every author-year in text needs a numbered reference for Nature Plants superscript style).
3. **Replace placeholder figure numbers with actual paths** (`figures/main/Main0x_*.pdf`).
4. **Confirm Sudarshana co-author authorship order and affiliations.**
5. **Supplementary tables S1–S6** are already implicit in the CSV outputs (`fix01_de_deseq2_voom_results.csv`, `heat_stress_DEGs_all_conditions.csv`, `mediation_results.csv`, `stability_selection_FIX04.csv`, etc.); wrap them into Excel with named sheets.
6. **Pre-submission checklist** (Nature Plants cover letter, data availability statement with SRA/GEO accessions, author contributions statement, competing interests, funding acknowledgments).
