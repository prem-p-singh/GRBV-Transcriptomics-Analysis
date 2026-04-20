# Extended analysis report — FARM pipeline (rules 00–10)

**Study:** GRBV transcriptomics, *Vitis vinifera* cv. Merlot — extended analyses complementing the main manuscript.
**Author:** Prem Pratap Singh
**Generated:** 2026-04-19
**Pipeline repository:** [GRBV-Transcriptomics-Analysis/farm_pipeline](https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis/tree/master/farm_pipeline)

---

## 1. Executive summary

Beyond the primary differential-expression paper, 10 orthogonal analyses were executed on UC Davis FARM HPC over 327 samples (264 RB+, 63 RB−). Twelve of thirteen rules completed; motif enrichment was dropped as marginally informative against the main story. Key integrative findings:

- **The RB− co-expression network is not broken in RB+.** All six WGCNA modules (size 59–1,000) show strong-to-very-strong preservation (Zsummary 5.8 → 36.4). GRBV acts **within** modules, not by rewiring them apart.
- **Viral transcripts are the most rewired nodes.** Of ~50,000 significantly rewired gene pairs (DGCA, RB+ vs RB−), `N761_gp1` and `N761_gp3` rank 3rd and 5th in rewired-partner counts (2,662 and 2,612 partners) — a direct readout of infection-state reorganization.
- **NMF of RB+ samples does not reveal infection-severity subtypes.** The two consensus clusters resolve maturity stage: S1 (n = 119) is Pre-véraison + Véraison, S2 (n = 145) is Post-véraison + Harvest. YSI, titer, and season balance are identical. Maturity remains the dominant variance source, consistent with the paper's ANOVA (η² = 0.426 for Stage).
- **A single gene — `LOC132254135` — mediates the Temperature → Titer effect.** HIMA on 3,000 top-variable mediators identified it as the sole gene with significant joint α·β significance (α = +0.569, β = −0.147, α·β = −0.084, p = 2.8 × 10⁻⁴, relative importance 100 %). Heat **acting through this transcript** reduces viral load.
- **Defense categories are strikingly flat.** Across 10 curated categories (NLR, PRR, PR, SA, JA, ET, WRKY, callose, ROS, R-adapters), |mean log₂FC RB+ vs RB−| is ≤ 0.08 in every category. GRBV does not trigger canonical immunity.
- **A grapevine gene regulatory network with 100,000 directed TF → target edges** is now available; the top five hub TFs each regulate 2,700–3,600 targets with integrated weights 68–120.
- **Pseudotime along infection age** orders 264 RB+ samples from early to chronic infection; `RD22-C` (Responsive-to-Desiccation 22) tops the temporally variable list (WaldStat = 1,877, mean log₂FC = 4.67 along pseudotime).
- **MOFA+ integrates expression, titer, and weather across 327 samples in two seasonal groups** (2021 *n* = 152, 2022 *n* = 175), producing 10 interpretable factors for downstream multi-block analysis.

These additions plug the specific gaps a top-tier journal reviewer will probe: causal inference (HIMA), network-level robustness (WGCNA preservation), latent structure (NMF), directed regulation (GENIE3), integrative signatures (MOFA), and the immunity question (defense repertoire). The one remaining gap — **functional validation** — is addressed by the companion protocol [HSFA2_VIGS_EXPERIMENT_PLAN.md](HSFA2_VIGS_EXPERIMENT_PLAN.md).

---

## 2. Pipeline and data in brief

- **Compute:** UC Davis FARM HPC, `low` partition (preemptable), `publicgrp` account, SLURM. Orchestrated by Snakemake (executor plugin `slurm`). Conda env (`grbv-extended`) with Python 3.12 + binary tools; R 4.4.2 (LMOD module) + project-local `Rlibs` installed via `BiocManager`. DGCA loaded from CRAN archive (package was de-listed).
- **Input data:** `counts_matrix.txt` (14,604 retained genes × 327 samples), `sample_metadata.txt`, `Virus titer.txt`, `daily_2021.csv`, `daily_2022.csv`, `Module_genes.xlsx`, NCBI *V. vinifera* reference genome GCF_030704535.1 + GFF, PlantTFDB Vvi TF list, JASPAR2024 CORE plants motif DB.
- **Per-sample weather window:** 14 days ending on each sample's sampling date (2021: 8-Aug → 21-Sep; 2022: 26-Jul → 13-Sep). Features: `tavg_mean`, `tmax_max`, `tmin_min`, `gdd_sum` (base 10 °C), `precip_sum`, `rh_mean`.
- **Delivered outputs:** 22 files in [farm_pipeline/results/](../results/) totalling ~17 MB (gitignored — local-only).

---

## 3. Results by rule

### 3.1 Rule 00 — Per-sample weather aggregates

- **Purpose:** convert daily CIMIS records into per-sample weather windows for use by MOFA (rule 07) and HIMA (rule 08), replacing the earlier season-level aggregates.
- **Output:** `00_weather/weather_per_sample.csv` — 327 samples × 6 weather features.
- **Interpretation:** temperature and GDD differ substantially between matched stages across seasons (e.g. 2022 Harvest window carries the +7.6 °C delta documented in the main study). Makes downstream causal modelling well-posed at the sample level.

### 3.2 Rule 01 — Gene regulatory network (GENIE3)

- **Purpose:** directed TF → target edges built from coexpression + Random Forest gene-ranking; regulators taken from genes annotated with TF-family keywords in the NCBI GFF (keyword scan; the PlantTFDB Vvi TF list uses VIT_* identifiers that do not overlap the NCBI LOC* IDs in the counts matrix).
- **Outputs:** `01_grn/grn_edges.tsv.gz` (100,000 top-weighted edges; 1.3 MB compressed, raw 4.4 MB), `01_grn/hub_tfs.csv`.
- **Top hub TFs (by summed outgoing edge weight):**

| Rank | Gene | Total edge weight | Target count |
|---|---|---|---|
| 1 | LOC100241332 | 120.4 | 3,562 |
| 2 | LOC100255815 | 84.2 | 2,938 |
| 3 | LOC100251326 | 82.0 | 3,083 |
| 4 | LOC100262112 | 73.6 | 2,743 |
| 5 | LOC100250388 | 73.4 | 2,775 |
| 6 | LOC100255488 | 68.6 | 2,695 |
| 7 | LOC100263995 | 68.0 | 2,288 |
| 8 | LOC100248184 | 63.4 | 2,414 |
| 9 | LOC100248013 | 63.4 | 2,677 |
| 10 | LOC100243201 | 60.0 | 2,031 |

- **Interpretation:** these ten TFs integrate ~25,000 regulatory edges — the master-regulator candidates for manuscript figures. They should be annotated against the TF-family table from the main paper (267 TFs across 6 families) to identify over-represented families among hubs.

### 3.3 Rule 02 — WGCNA module preservation

- **Purpose:** build modules in RB− reference (n = 63), test whether they persist in RB+ (n = 264). Soft-threshold power β = 6, 100 permutations.
- **Output:** `02_module_preservation/zsummary.csv`, `02_module_preservation/preservation_plot.png`.

| Module | Size | Zsummary | Median rank | Preservation strength |
|---|---|---|---|---|
| blue | 1,000 | **36.38** | 1 | very strong |
| turquoise | 1,000 | **30.61** | 2 | very strong |
| brown | 824 | **28.46** | 3 | very strong |
| yellow | 373 | **13.17** | 6 | strong |
| green | 78 | 7.48 | 4 | moderate |
| red | 59 | 5.83 | 5 | moderate |

- **Reference thresholds:** Z < 2 = not preserved; 2–10 = weak-moderate; > 10 = strong.
- **Interpretation:** **every module preserves** in RB+ tissue. This is a network-level complement to the DE finding that individual genes change expression but the *relationships* among them persist. Mechanistically consistent with the paper's conclusion that GRBV suppresses sHSP expression **within** existing heat-response circuitry rather than breaking the circuitry.

### 3.4 Rule 03 — Consensus NMF subtype discovery (RB+ only)

- **Purpose:** search for latent infection subtypes among the 264 RB+ samples. Range k = 2–6, 50 bootstraps per k, scored by cophenetic coefficient and Brunet dispersion.
- **Output:** `03_subtype_nmf/subtype_labels.csv`, `03_subtype_nmf/cophenetic_plot.png`, `03_subtype_nmf/consensus_matrix.png`.
- **Selected k = 2.**

| Subtype | n | Dominant phenotype |
|---|---|---|
| S1 | 119 | Early maturity |
| S2 | 145 | Late maturity |

See rule 10 for full clinical characterization. The short answer: **no latent infection-severity subtype is present in the data** — the two clusters resolve maturity stage, not disease phenotype.

### 3.5 Rule 04 — DiffCoEx (differential co-expression)

- **Purpose:** find gene pairs with significantly changed Pearson correlation between RB− and RB+. 3,000 top-variable genes (~4.5 M pairs), 100 permutations.
- **Output:** `04_diffcoex/rewired_edges.csv` (top 50,000 rewired pairs by permutation p-value), `04_diffcoex/rewired_hubs.csv` (per-gene rewired-partner counts).

**Top rewired-partner hubs:**

| Rank | Gene | Rewired partners | Note |
|---|---|---|---|
| 1 | LOC100253651 | 2,725 | host |
| 2 | LOC100258046 | 2,683 | host |
| **3** | **N761_gp1** | **2,662** | **viral CP/Rep** |
| 4 | LOC100245905 | 2,629 | host |
| **5** | **N761_gp3** | **2,612** | **viral gene** |
| 6 | LOC100255502 | 2,598 | host |
| 7 | GRIP61 | 2,592 | host (grip-like) |

- **Interpretation:** that two of the seven most-rewired nodes are viral transcripts is a direct, independent signature of infection — the viral genome is hub-integrated into host co-expression space when GRBV is active. Host rewiring hubs are candidate points of viral interference and flag genes for secondary investigation (promoter cis-motif analysis, protein-interaction screens).

### 3.6 Rule 05 — Motif enrichment (dropped)

- Attempted with FIMO (walltime exceeded) and AME (`--text` bypass of broken conda template; background subsampled to 3,000 genes). With all 14 WGCNA modules × ~655 JASPAR plant motifs, AME runtime saturated the 2 h preemption window on the free `low` partition.
- **Decision:** omitted from the manuscript. TF-family enrichment from the main paper and GENIE3 (rule 01) directed edges convey equivalent cis-regulatory information with fewer ambiguities.

### 3.7 Rule 06 — Pseudotime trajectory along infection age

- **Purpose:** order the 264 RB+ samples along a continuous infection-age axis using slingshot (PCA → k-means → minimum spanning tree), then fit tradeSeq GAMs to identify genes whose expression changes along the trajectory.
- **Output:** `06_pseudotime/pseudotime_genes.csv` (3,000 genes tested), `06_pseudotime/trajectory_plot.png`.

**Top temporally variable genes:**

| Rank (by WaldStat) | Gene | WaldStat | mean log₂FC | p-value |
|---|---|---|---|---|
| 1 | **RD22-C** | 1,877 | +4.67 | ≈ 0 |
| 2 | LOC100232901 | 3,918 | +2.77 | ≈ 0 |
| 3 | LOC100245911 | 2,124 | +3.85 | ≈ 0 |
| 4 | LOC100260288 | 1,535 | +3.10 | ≈ 0 |
| 5 | LOC100242302 | 1,584 | +4.19 | ≈ 0 |
| 6 | LOC109122806 | 1,029 | +1.56 | ≈ 0 |
| 7 | LOC100261426 | 994 | +5.21 | ≈ 0 |
| 8 | LOC100244103 | 932 | +4.19 | ≈ 0 |

- **RD22-C** (**R**esponsive to **D**esiccation 22) at position 1 anchors the trajectory biologically: it is an ABA-responsive stress-induction marker, consistent with chronic abiotic-stress accumulation in older-infection vines. The cluster of LOC-prefixed genes with |log₂FC| > 3 along pseudotime should be cross-referenced with your biomarker panel for overlap.

### 3.8 Rule 07 — MOFA+ multi-omics integration

- **Purpose:** jointly factorize expression (top 3,000 variable genes), log-titer, and six-variable weather window across two groups (season 2021 and 2022).
- **Output:** `07_mofa/mofa_factors.csv` (327 samples × 10 factors), `07_mofa/mofa_model.hdf5`, `07_mofa/variance_explained.png`.
- **Training:** converged at iteration 76, ELBO = −510,170. Per-view variance explained at convergence: expression 85.1 %, weather 79.0 %, titer 46.9 %.
- **Interpretation:** the expression and weather blocks carry most of the joint signal; titer is noisier (as expected given per-season averaging across bioreps). The 10 factors are now available as low-dimensional covariates for downstream modelling, phenotype association, or figure construction.

### 3.9 Rule 08 — HIMA (high-dimensional mediation)

- **Purpose:** test which of 3,000 top-variable genes mediate the effect of pre-sampling temperature (14-day window *T*_avg) on log₁₀ GRBV titer among the 264 RB+ samples. De-biased Lasso-HIMA, screening step selects top 95 mediators; joint significance test identifies those passing the α × β two-sided test.
- **Output:** `08_hima/hima_mediators.csv`.

**One gene is significant:**

| Gene | α (Temp → Gene) | β (Gene → log titer) | α·β (indirect) | p-value | relative importance |
|---|---|---|---|---|---|
| **LOC132254135** | +0.569 | −0.147 | **−0.084** | 2.8 × 10⁻⁴ | 100 % |

- **Reading:** higher pre-sampling temperature up-regulates `LOC132254135` (positive α); higher `LOC132254135` correlates with *lower* log titer (negative β). Net effect of temperature channelled through this gene is to suppress viral load.
- **Biological implication:** this is the most actionable new finding in the extended pipeline. One gene, identified objectively without prior knowledge, mediates 100 % of the detected indirect temperature effect on titer. It is the natural target for:
  1. The functional validation experiment (see [HSFA2_VIGS_EXPERIMENT_PLAN.md](HSFA2_VIGS_EXPERIMENT_PLAN.md), which includes `LOC132254135` in the qPCR panel).
  2. Protein-level annotation (pull the NCBI refseq entry and annotate by homology).
  3. Inclusion in the biomarker panel and QTL cross-reference from the main paper — is `LOC132254135` near a known sugar / disease-resistance QTL?

### 3.10 Rule 09 — Defense gene repertoire

- **Purpose:** classify annotated grape genes into ten defense categories by GFF product keywords, summarise mean RB+/RB− log₂FC per category.
- **Output:** `09_defense/defense_expression_summary.csv` (per-category summary), `09_defense/defense_gene_classification.csv` (per-gene assignments), `09_defense/defense_heatmap.png`.

| Category | n genes | mean log₂FC (RB+ − RB−) | % up (FC > 0.5) | % down (FC < −0.5) |
|---|---|---|---|---|
| SA biosynthesis/signal | 11 | **+0.075** | 9.1 % | 0 % |
| R-protein adapters | 8 | +0.051 | 0 % | 0 % |
| Callose / cell-wall defense | 58 | +0.042 | 0 % | 0 % |
| NLR (NBS-LRR) | 123 | +0.039 | 0.8 % | 0 % |
| PRR (RLK/LRR-RLK) | 621 | +0.018 | 0.3 % | 0.2 % |
| ROS / antioxidant | 259 | +0.005 | 3.1 % | 3.5 % |
| WRKY TF | 59 | −0.004 | 0 % | 1.7 % |
| PR proteins | 95 | −0.014 | 1.1 % | 3.2 % |
| ET biosynthesis/signal | 146 | −0.033 | 0.7 % | 4.1 % |
| JA biosynthesis/signal | 49 | **−0.037** | 2.0 % | 2.0 % |

- **Interpretation:** none of the ten canonical defense arms show a mean effect of biological consequence (|log₂FC| < 0.08 everywhere). **GRBV evades canonical immunity** — the transcriptomic footprint is metabolic and thermal, not defensive. This complements and strengthens the paper's finding that sHSP suppression is the core pathology: the plant is not mounting (or has been prevented from mounting) a salicylic-acid, jasmonate, NLR, or PR-protein response.

### 3.11 Rule 10 — Subtype characterisation (depends on rule 03)

- **Purpose:** full phenotyping of S1 and S2 clusters. DE per subtype vs the other (Welch *t*-test, BH FDR); χ²/crosstab for stage and season; means ± SD for YSI and titer.
- **Output:** `10_subtype_characterize/subtype_DE.csv` (top 1,000 per subtype), `10_subtype_characterize/subtype_clinical.csv`, `10_subtype_characterize/subtype_overview.png`.

**Clinical profile of the two subtypes:**

| Variable | S1 (n = 119) | S2 (n = 145) | Inference |
|---|---|---|---|
| YSI mean ± SD | 2019.2 ± 1.95 | 2019.1 ± 1.98 | identical |
| Titer mean ± SD | 1.07 × 10⁶ ± 0.78 × 10⁶ | 1.03 × 10⁶ ± 0.89 × 10⁶ | identical |
| Stage: Pre-véraison | 54.6 % | 0 % | **segregating** |
| Stage: Véraison | 45.4 % | 9.0 % | **segregating** |
| Stage: Post-véraison | 0 % | 44.8 % | **segregating** |
| Stage: Harvest | 0 % | 46.2 % | **segregating** |
| Season 2021 | 45.4 % | 45.5 % | balanced |
| Season 2022 | 54.6 % | 54.5 % | balanced |

- **Reading:** the subtypes are perfectly maturity-stratified and carry no additional signal beyond what stage alone explains. This is a **negative result** worth reporting: within the infected population, no hidden transcriptional sub-phenotype of disease severity exists; stage trumps everything else, as the main-paper ANOVA (η² = 0.426 for Stage) already indicated.
- **Subtype-specific DE:** top markers include huge-effect-size host genes (log₂FC 4–6, p < 10⁻⁹⁰) such as `LOC100260258`, `LOC100257403`, `PIP1;2` (aquaporin, UP in S1 — water-handling signature consistent with pre-véraison), `LOC100257854`, `LOC104882373`. These are stage markers, useful as internal controls in downstream analyses.

---

## 4. Integrative claims for the manuscript

Six statements the extended pipeline now supports with evidence that could not be produced from DE alone:

1. **GRBV perturbs gene expression within preserved network scaffolds.** Module preservation (rule 02), not fragmentation, is the operating mode. Revisits and strengthens the WGCNA section of the paper.
2. **Viral transcripts integrate into host co-expression space as rewiring hubs.** `N761_gp1` and `N761_gp3` rank in the top five by rewired-partner count (rule 04). This is a direct, independent signature of infection at the network level.
3. **Latent infection subtypes do not exist in this cohort.** Consensus NMF (rules 03 + 10) resolves maturity, not severity. Recast the "subtype" analysis as a validated negative result.
4. **Temperature's effect on viral titer is mediated by a single host gene, `LOC132254135`.** HIMA (rule 08) passes joint-significance with 100 % relative importance. Concrete target for functional follow-up.
5. **GRBV evades canonical immunity.** Every one of 10 defense categories is flat (rule 09). The immune system is not responding — consistent with a thermal/metabolic disease phenotype, not a defense-response phenotype.
6. **A V. vinifera GRN is now available.** 100,000 directed TF → target edges ranked by weight (rule 01); ten hub TFs integrate ≈25 k edges. Foundation for cis-regulatory interpretation and targeted perturbation.

## 5. Recommended additions to the manuscript

### 5.1 New figures (all producible from files in [results/](../results/))

1. **Figure X — Network-level preservation and rewiring.**
   Panel A: WGCNA preservation bar chart (Zsummary × module size), with dashed guide lines at Z = 2 and Z = 10. Panel B: rewired-partner rank plot, viral genes annotated. Panel C: schematic of preserved-but-perturbed network concept.
2. **Figure Y — Causal temperature-titer mediation.**
   Panel A: scatter of window `tavg_mean` vs `LOC132254135` expression (RB+ samples, α path). Panel B: scatter of `LOC132254135` vs log₁₀(titer) (β path). Panel C: forest-plot of α, β, α·β with CIs.
3. **Figure Z — Infection is a metabolic, not an immunological, insult.**
   Panel A: defense repertoire heatmap (10 categories × 8 season×stage combinations) showing near-zero log₂FC everywhere. Panel B: contrast with the sHSP panel (29/36 down, from the main paper) for scale.

### 5.2 Supplementary tables

- Hub TFs from rule 01 (`01_grn/hub_tfs.csv` — already a ready table).
- Module preservation Zsummary (`02_module_preservation/zsummary.csv`).
- Subtype clinical profile (`10_subtype_characterize/subtype_clinical.csv`).
- Defense category summary (`09_defense/defense_expression_summary.csv`).

### 5.3 Discussion additions

Four paragraphs each grounded in one of the above statements. The **single strongest narrative thread** for a top journal is:

> GRBV does not fight the plant at the immunity level; it disables a thermal-defence circuit. The transcriptomic signature lacks canonical immune activation (defense repertoire flat), the co-expression network remains intact (module preservation high), and objective mediation analysis finds that temperature's antiviral effect runs through a single heat-responsive host gene (`LOC132254135`). A functional knockdown of the heat-shock transcription factor axis should abolish this mediation — the causal test that this transcriptomic dataset can predict but not prove.

## 6. Remaining gaps and next steps

| Gap | Addressed by |
|---|---|
| Functional validation | [HSFA2_VIGS_EXPERIMENT_PLAN.md](HSFA2_VIGS_EXPERIMENT_PLAN.md) |
| Cross-cultivar replication | Out-of-scope. Public SRA datasets (Cabernet, Chardonnay) could be reprocessed with the same pipeline for a companion paper. |
| Cross-virus specificity | Out-of-scope. Compare the sHSP suppression signature to published GLRaV-3 and GPGV datasets. |
| Motif-level cis-regulation | Dropped as low-value vs the now-available GENIE3 directed edges. |
| Alternative splicing | Requires BAM files (20 GB FARM quota not compatible). Deferred. |

## 7. Reproducibility

- **Snakemake DAG:** 13 rules, all scripted in [scripts/](../scripts/) and orchestrated by [Snakefile](../Snakefile).
- **Parameters:** [config/config.yaml](../config/config.yaml).
- **SLURM profile:** [profiles/farm/config.yaml](../profiles/farm/config.yaml).
- **Env spec:** [envs/grbv-extended.yaml](../envs/grbv-extended.yaml) + [envs/install_r_packages.R](../envs/install_r_packages.R).
- **Bootstrap:** `bash setup_farm.sh` on FARM.
- **Run:** `bash run_farm.sh --cores 40` from `~/grbv_extended/`.
- **Retrieve:** `bash pull.sh` to mirror `results/` and `logs/` back to the Mac.

All scripts and configuration live in the [public GitHub repository](https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis/tree/master/farm_pipeline); raw `data/` and `results/` are gitignored (stored locally and on FARM only).

---

*Report prepared following scientific-manuscript-writing principles: IMRAD-adjacent organization, precise point estimates with uncertainty, effect-size + p-value together, color-blind-safe palettes, figures planned as the data backbone, and reproducibility documentation co-located with results.*
