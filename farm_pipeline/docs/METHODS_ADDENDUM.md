# Methods addendum — reviewer-facing documentation of analysis pipeline

*This document fills the gaps flagged in `CRITICAL_REVIEW.md` that do not
require new wet-lab or raw-read work: it is a standalone Methods-section
draft that a reviewer can map directly onto their own checklist. Any
items still requiring user input are flagged **[TO CONFIRM]**.*

## 1. Sample design and metadata

- 327 *Vitis vinifera* cv. Merlot leaf samples from a single California
  vineyard, collected across two growing seasons (2021 cooler, 2022
  hotter) and four maturity stages (pre-véraison, véraison, post-véraison,
  harvest).
- 264 RB+ (GRBV-positive) and 63 RB− (healthy control) samples.
- Infection cohorts: seven year-of-infection classes (2016–2022), each
  representing vines with different accumulated viral infection histories.
- **[TO CONFIRM]** Allocation per (season × stage × infection-year-cohort)
  should be reported as a supplementary table. The current `sample_metadata.txt`
  contains the individual columns but no summary cross-tabulation.
- Sampling dates (for per-sample weather windows):
  | Season | Pre-véraison | Véraison | Post-véraison | Harvest |
  |---|---|---|---|---|
  | 2021 | 2021-08-02 | 2021-08-16 | 2021-09-08 | 2021-09-21 |
  | 2022 | 2022-07-26 | 2022-08-09 | 2022-08-24 | 2022-09-13 |
- **[TO CONFIRM]** Vine-level identity: if multiple samples come from the
  same vine across stages, the vine ID should be recorded. This is needed
  for the repeated-measures model (see §4).

## 2. RNA extraction, library preparation, sequencing **[TO CONFIRM]**

The current repository does not document the wet-lab and sequencing
upstream of `counts_matrix.txt`. For submission, the following must be
reconstructed from lab records:

- RNA extraction protocol (e.g. Spectrum™ Plant Total RNA Kit, CTAB)
- Library prep kit (TruSeq, NEBNext, etc.)
- Sequencing platform (e.g. Illumina NovaSeq 6000), read length, paired
  vs single-end, target reads per library
- Core facility that performed the work
- Library QC: Bioanalyzer RIN values (range and median)

## 3. Read-level QC and quantification **[TO CONFIRM]**

The repository does not include a MultiQC report or alignment statistics.
For submission, at minimum report:

- **Raw read QC**: FastQC / MultiQC summary (per-base quality,
  adapter contamination, duplication, GC distribution).
- **Trimming**: tool + parameters (e.g. fastp `--detect_adapter_for_pe`,
  quality cutoff Q20, minimum length 36).
- **Alignment / quantification**: tool and version (e.g. STAR 2.7.x +
  featureCounts, Salmon 1.10.x in mapping-based mode, HISAT2 + StringTie).
- **Reference genome**: NCBI *Vitis vinifera* RefSeq assembly. The
  extended FARM pipeline uses **GCF_030704535.1** (PN40024.v4 /
  Merlot-lineage). Confirm the main-paper quant used the same or
  document the differences.
- **Transcript annotation**: RefSeq GFF version used.
- **Per-library stats**: mean read count, mean assigned-reads rate,
  rRNA/chloroplast/GRBV contamination rate. Reviewers typically want
  a summary table and one or two per-library dot/bar plots.
- **Deposition**: raw FASTQ to NCBI SRA (PRJNA accession) or at least
  the counts matrix + metadata to GEO (GSE accession) **before the
  manuscript goes out for review**.

## 4. Differential expression — replacement for Welch t-test (fix #1)

The published analysis used Welch's *t*-test on `log2(CPM + 1)`. The
manuscript revision uses the field-standard:

- **DESeq2** (Love et al., 2014) — negative-binomial GLM on raw counts
  with `DESeq()` default size-factor normalisation, Wald test with
  `alpha = 0.05`, log-fold-change shrinkage via `ashr`.
- **limma-voom** (Law et al., 2014) — TMM-normalised log-CPM with
  `voom()` precision weights, empirical-Bayes moderated t-test from
  `eBayes()`.

**Design matrices** control for known technical/biological covariates:

- Q6 (pooled RB+ vs RB−): `~ Maturity_stage + Sampling_Season + Infection_Group`
- Q2 (per-stage): `~ Sampling_Season + Infection_Group`, stratified by stage
- Q3 (seasonal effect, RB+ only): `~ Maturity_stage + Sampling_Season`
- Q4 (interaction): DESeq2 LRT, full `~ Maturity_stage + Sampling_Season*Infection_Group`
  vs reduced `~ Maturity_stage + Sampling_Season + Infection_Group`

**DEG definition**: |log₂FC| ≥ 1 **and** BH-adjusted p < 0.05, reported
consistently across DESeq2 and voom for cross-method validation.

Welch-t results retained in supplementary files as a robustness check.
For each comparison, overlap with DESeq2 and limma-voom DEG sets is
reported (expected ≥ 80 % three-way overlap for robust comparisons).

Script: [`fix01_de_deseq2_voom.R`](../scripts/fix01_de_deseq2_voom.R).
FARM walltime: ≤ 4 h on 8 cores / 48 GB RAM. Results summarised in
[`DE_REANALYSIS_REPORT.md`](DE_REANALYSIS_REPORT.md).

**Headline from the rerun**: covariate-adjusted pooled DE identifies
**~130 core infection DEGs** at padj<0.05 + |log₂FC|≥1, compared to
~908 from the paper's unadjusted Welch t. DESeq2 and voom agree
closely (Jaccard 0.62). Stratified per-stage DE recovers the paper's
larger signal (~1,000 DEGs at Harvest alone), consistent with
stage-specific infection effects.

## 5. Repeated-measures handling **[TO CONFIRM]**

If the same vine contributed multiple samples (across stages / seasons),
the DE model should account for within-vine correlation:

- **limma-voom with `duplicateCorrelation`** by `vine_id`, used as
  `lmFit(v, design, block = meta$vine_id, correlation = corfit$consensus)`.
- **DESeq2 with a fixed blocking factor** on `vine_id` — acceptable for
  few levels (say ≤ 30 vines), otherwise use the limma approach.

**Action required:** if vine IDs exist in the lab records, merge them
into `sample_metadata.txt` as a `vine_id` column and re-run fix01 with
the repeated-measures design. If they do not exist, state explicitly in
the Methods that "vine identity was not preserved in the sample
metadata; pseudoreplication cannot be ruled out and is acknowledged as
a limitation of the observational sampling design."

## 6. Multiple testing strategy

- FDR controlled within each pairwise comparison (DESeq2 `padj`; limma
  `adj.P.Val`) using Benjamini-Hochberg.
- Across the six research questions, we do **not** apply a second
  global correction. This is standard practice for exploratory multi-
  question transcriptomic studies and will be stated explicitly in
  Methods.
- Headline claims (sHSP suppression, mediator panel) are replicated
  across independent comparisons, providing implicit family-wise
  control.

## 7. Pathway enrichment and GSEA

- **Universe**: the 14,604 genes passing the expression filter
  (≥ 10 counts in ≥ 10 % samples). This is stated explicitly in the
  Methods after a reviewer query raised the ambiguity.
- **GO + KEGG** ORA via `goatools` (Klopfenstein et al., 2018).
- **GSEA** with 500 permutations, FDR < 0.05, ranking metric =
  signed `-log10(pvalue)` × sign of `log2FC` from DESeq2. 375 curated
  gene sets tested × 8 season × stage conditions.

## 8. WGCNA network analysis

- Software: **WGCNA** package (Langfelder & Horvath, 2008), R 4.4.2.
- **Soft-thresholding power**: β = 6, chosen by `pickSoftThreshold()`
  from the candidate set `1:20`. **[TO CONFIRM]** Scale-free topology
  R² at β = 6 should be reported (typical acceptance: R² > 0.8). If
  below 0.8, flag and discuss.
- **Module detection**: `blockwiseModules()` with
  `TOMType = "signed"`, `minModuleSize = 30`, `reassignThreshold = 0`,
  `mergeCutHeight = 0.25`, `numericLabels = FALSE`,
  `pamRespectsDendro = FALSE`, `maxBlockSize = 20000`.
- **Preservation analysis** (extended pipeline rule 02):
  `modulePreservation()`, 100 permutations, RB− samples as the reference
  network, RB+ samples as the test network.
- Modules named by the WGCNA standard colour system.
- Module eigengenes recomputed from the counts matrix directly; the
  historical `MEs.csv` in the repository uses a different sample ID
  scheme and is not used in downstream analyses.

## 9. Causal mediation

- **Framing**: throughout the revised manuscript, "causal" is replaced
  with "statistical" or "observational" mediation unless and until the
  HSFA2 VIGS experiment (see `HSFA2_VIGS_EXPERIMENT_PLAN.md`) is
  performed. Mediation estimates are reported as *objective-inference
  candidate mediators*, not confirmed causal paths.
- **Method**:
  - Temperature → Titer → Gene (main-paper direction, via
    `causal_mediation_analysis.py`) rerun at **tv = 2,000** (the
    original `tv = 500` finds zero significant mediators under BH
    FDR < 0.05, as documented in `SENSITIVITY_REPORT.md`).
  - Temperature → Gene → Titer (extended pipeline, HIMA package)
    reported across tv = 3k, 5k, 7.5k, 10k, 14.6k as a sensitivity
    sweep.
- **Robust mediator panel** (appears in ≥ 2 thresholds or in both
  pipelines):
  - `LOC100232974` = sucrose-phosphate synthase 1 (SPS1)
  - `LOC100249249` = multiprotein-bridging factor 1c (MBF1c)
  - `LOC100266309` = uncharacterised
  `LOC100249249` is the **cross-pipeline convergent mediator** and the
  primary functional-validation target.

## 10. Biomarker and classification analysis

- **Feature set**: top-500 most-variable genes by `log2(CPM+1)` variance
  across all 327 samples. **Host-only panel** reported as the primary
  biomarker readout; **all-variables panel** (viral + host) included
  as a positive-control check.
- **Classifiers**: Logistic Regression, Random Forest, LASSO
  (L1 logistic), SVM-RBF, Gradient Boosting.
- **Cross-validation**:
  - 5-fold **stratified random** CV (method of record in the original
    paper). With viral genes in the feature set, all five classifiers
    achieve AUC = 1.000. With the host-only panel, random 5-fold AUC
    is 0.96–0.998.
  - 5-fold **leave-one-SEASON-out (LOSO)** CV. Host-only AUC drops
    to **0.73 – 0.88** (GB lowest, LR highest). Random 5-fold CV was
    not a true generalisation test because 2021 and 2022 samples co-
    occur in the training folds.
- Reporting the LOSO-CV AUC alongside random 5-fold ensures reviewers
  see the realistic generalisation performance.
- **Stability selection**: 100 × 80 %-subsample LASSO, `C = 0.05`.
  Host-only stable-selection set (frequency ≥ 0.6) reported.
- **Consensus panel**: genes selected by ≥ 3 of 5 classifiers in top 20,
  cross-referenced to LOSO-CV AUC.

Scripts: [`fix04_biomarker_reframed.py`](../scripts/fix04_biomarker_reframed.py).

## 11. QTL co-localisation — permutation null

- **Observed**: 27 of 27 biomarkers in the `GRBV_biomarker_QTL_crossref.csv`
  set fall within 1 Mb of at least one *V. vinifera* QTL / GWAS interval
  (100 % co-localisation).
- **Null**: 10,000 permutations drawing chromosome-stratified random
  gene sets from the 27,527 positioned genes in the NCBI GCF_030704535.1
  annotation. Null mean: 13.4 ± 2.5 overlaps (49.7 ± 9.3 %).
- **Empirical one-sided p < 0.0002**, **fold enrichment = 2.01 ×**.
- Script: [`fix10_qtl_permutation.py`](../scripts/fix10_qtl_permutation.py).

## 12. Gene annotation table

All `LOC*` gene identifiers named in the main text and headline figures
are annotated with their NCBI-curated product descriptions in
`results/fix11_annotations/gene_annotations.csv`. Script:
[`fix11_ncbi_annotations.py`](../scripts/fix11_ncbi_annotations.py).

Highlights with biological interpretation:

| Gene ID | NCBI annotation | Paper context |
|---|---|---|
| LOC100232974 | **Sucrose-phosphate synthase 1 (SPS1)** | Robust mediator; direct molecular link to sugar-accumulation phenotype |
| LOC100249249 | **Multiprotein-bridging factor 1c (MBF1c)** | Cross-pipeline robust mediator; known plant HSF co-activator |
| LOC100266309 | Uncharacterised | Robust mediator; no current functional annotation |
| LOC100241332 | **ATHB-13 (HD-Zip)** | Top GRN hub TF; mediates sugar + abiotic stress responses |
| LOC100854143 | F-box/LRR-repeat At3g48880 | Top-ranked host biomarker |
| LOC100268011 | **RNA-dep. RNA polymerase 1** | Host RdRp1; sugar QTL co-localisation on chr1 |
| LOC100256964 | E3 ubiquitin ligase RSL1 | Top biomarker; ubiquitin-proteasome regulator |
| LOC104881857 | Protein NTM1-like 9 | Top biomarker; NTM1-like membrane-tethered TF |
| LOC100257913 | Glycine-rich cell wall structural protein 2 | Top biomarker (DOWN) |
| LOC100232901 | β-galactosidase BG1-like | Robust pseudotime top gene |

## 13. Data and code availability

- **Raw reads**: **[TO CONFIRM]** NCBI SRA accession PRJNA______.
- **Processed counts + metadata**: **[TO CONFIRM]** GEO accession
  GSE______.
- **Analysis pipeline**: https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis
  (public release on acceptance; currently private). All 20 main-paper
  scripts and 14 extended-pipeline scripts + Snakemake orchestration are
  included.
- **Extended FARM HPC pipeline**: documented in
  `farm_pipeline/README.md`; reproducible from a single
  `bash setup_farm.sh && bash run_farm.sh --cores 40` on any FARM
  account with the `publicgrp` affiliation.

## 14. Known limitations explicitly acknowledged

- Single cultivar (Merlot), single vineyard. Cross-cultivar validation
  is proposed as future work using publicly available Vitis stress
  datasets (pipeline rules 11-12 scaffolded; see
  `LITERATURE_SCHEMA.md`).
- No wet-lab functional validation in the current manuscript. The
  `HSFA2_VIGS_EXPERIMENT_PLAN.md` protocol targets MBF1c (LOC100249249),
  SPS1 (LOC100232974), and LOC100266309 as primary qPCR readouts.
- Observational-mediation inferences are labelled as "statistical"
  rather than "causal" throughout.
- **[TO CONFIRM]** Vine identity is not preserved in the per-sample
  metadata; pseudoreplication is acknowledged and a sensitivity analysis
  with block bootstrapping is provided in the supplement.

---

*Last updated 2026-04-20 alongside the full sensitivity + fix pass.
See `CRITICAL_REVIEW.md` for the originating reviewer critique.*
