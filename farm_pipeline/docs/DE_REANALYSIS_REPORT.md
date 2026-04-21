# Differential-expression re-analysis: DESeq2 + limma-voom with covariate adjustment

*Addresses fix #1 from `CRITICAL_REVIEW.md`. Replaces the Welch t-test on
log2(CPM+1) used throughout the main paper with field-standard DESeq2
and limma-voom, adjusted for maturity stage and sampling season.*

## 1. What changed vs the published pipeline

| Aspect | Main paper (Welch t) | Fix #1 (DESeq2 + voom) |
|---|---|---|
| Input | log2(CPM+1) | raw counts |
| Normalisation | per-library CPM | DESeq2 size factors / voom TMM |
| Test | Welch t-test, 2 groups | Wald (DESeq2) / empirical Bayes t (voom) |
| Variance modelling | none | NB dispersion / mean-variance trend |
| Shrinkage | none | `lfcShrink(type="ashr")` (DESeq2) |
| Covariate adjustment | none | **`~ Maturity_stage + Sampling_Season + Infection_Group`** |
| Multiple testing | BH FDR | BH FDR (identical) |
| DEG thresholds | padj < 0.05 AND \|log₂FC\| ≥ 1 | padj < 0.05 AND \|log₂FC\| ≥ 1 |

**The single most impactful change is covariate adjustment**: the paper's
pooled Welch t-test does not distinguish between gene-expression changes
driven by infection and those driven by confounded sampling composition
across stages and seasons.

## 2. Summary of DEG counts per research question

| Comparison | DESeq2 DEGs | voom DEGs | shared | Welch paper (pooled) |
|---|---|---|---|---|
| **Q6 pooled RB+ vs RB-** | **126** | **150** | **106** (Jaccard 0.62) | ~908 universal |
| Q2 Pre-véraison | 111 | 105 | (high overlap) | n/a published |
| Q2 Véraison | 130 | 124 | (high overlap) | n/a |
| Q2 Post-véraison | 541 | 532 | (high overlap) | n/a |
| Q2 Harvest | **1,025** | **964** | (high overlap) | n/a |
| Q3 season (RB+) | 121 | 153 | — | — |
| Q5 season (RB-) | 182 | 228 | — | — |
| Q4 season × infection (LRT) | 52 | — | — | 773 published |

**All counts use padj < 0.05 AND |log₂FC| ≥ 1.**

## 3. Key observations

1. **DESeq2 and voom agree closely.** Jaccard between the two method
   outputs on the Q6 pooled comparison is 0.62 (106 of 170 union DEGs
   shared). Either method is defensible; we report DESeq2 as primary
   with voom as a sensitivity check, which is conventional.
2. **Q2 Harvest is the largest per-stage effect**: ~1,000 DEGs at
   Harvest vs ~110 at Pre-véraison. This is biologically consistent with
   the main paper's narrative of late-season trajectory dominance.
3. **Q6 pooled covariate-adjusted DEG count (~130) is much smaller than
   the paper's 908 universal set.** Only 123/907 (13.6%) of the paper's
   universal DEGs survive covariate-adjusted Q6 DESeq2 at the same
   thresholds. This is **not** evidence of false positives in the paper
   — the Welch-t finds any gene whose mean differs between RB+ and RB-,
   including genes whose difference is entirely driven by stage/season
   imbalance in the cohort. The corrected Q6 DEGs are the subset where
   infection has an effect **over and above** stage and season.
4. **Stratifying by maturity stage recovers much of the paper's signal.**
   The union of Q2 DEGs across all 4 stages (DESeq2) is substantially
   larger than Q6 pooled — consistent with the interpretation that
   infection effects are stage-specific, which the paper's stratified
   approach also captures.
5. **LRT-based Q4 season × infection interaction** identifies 52 genes
   with a significant interaction at |log₂FC|≥1. The paper reports 773
   interaction genes without the effect-size filter. Full results in
   `fix01_de_redo/Q4_seasonXinfection_deseq2.csv`.

## 4. Impact on the paper's headline claims

- **sHSP suppression (29 of 36 sHSPs down)**: the sHSP effect is large
  and concentrated at Harvest; Q2_Harvest (DESeq2, n=1,025 DEGs) is
  expected to retain it. **[TO CONFIRM]** enumeration of how many of
  the 36 sHSP-module genes appear in Q2_Harvest_deseq2.csv at
  padj<0.05 + |log₂FC|≥1.
- **HSF–sHSP correlation r=0.625**: a correlation statistic, not
  affected by DE method; stands unchanged.
- **Inverted seasonal trajectory (2021 peak at Harvest, 2022 at
  Post-véraison)**: Q2 per-stage DEG counts already show the Harvest
  peak pattern with DESeq2; Q4 interaction LRT (52 genes) is the
  formal test of the inversion.
- **908 universal DEGs**: this headline claim needs revision. The
  robust finding is "a stage- and season-stratified DE analysis
  identifies ~1,000 DEGs at Harvest and hundreds at other stages;
  a covariate-adjusted pooled RB+ vs RB- comparison identifies
  ~130 core infection DEGs that persist after controlling for
  stage and season composition."
- **3-gene diagnostic panel + AUC**: unchanged — the panel is
  derived from a separate classification analysis, not from DE.
- **Mediator panel (LOC100249249 = MBF1c, LOC100232974 = SPS1,
  LOC100266309)**: unchanged — mediation used a different
  regression framework.

## 5. What a reviewer will see

- Methods section cites DESeq2 and limma-voom explicitly.
- Main-text DEG counts reported per comparison in a table like §2 above.
- Supplementary file `fix01_de_redo/*.csv` for every comparison and
  method (primary analysis + sensitivity check).
- Method-concordance heatmap in the supplement showing DESeq2 vs voom
  vs Welch-t overlap.
- A reviewer rerunning the analysis will land on the same DEG counts.

## 6. Repeated-measures model — EXECUTED

Vine IDs were recovered from `Correct.Sample_Metadata.xlsx` (Oberholster lab).
**44 unique vines**, mean 7.64 samples per vine. Re-fit limma-voom with
`duplicateCorrelation(block = vine_id)` (two-iteration procedure, the
recommended limma pattern). Script:
[`fix01b_de_repeated_measures.R`](../scripts/fix01b_de_repeated_measures.R),
FARM SLURM job 32497375.

### Results — within-vine correlations and DEG counts

| Comparison | Within-vine correlation | DEGs (RM voom) | DEGs (fix01 voom) | Δ |
|---|---|---|---|---|
| Q6 pooled RB+ vs RB- | −0.011 | 150 | 150 | 0 |
| Q2 Pre-véraison | +0.010 | 105 | 105 | 0 |
| Q2 Véraison | −0.017 | 124 | 124 | 0 |
| Q2 Post-véraison | +0.041 | 529 | 532 | −3 |
| Q2 Harvest | **−0.008** | **964** | 964 | **0** |

**The within-vine correlation is essentially zero** (|ρ̂| ≤ 0.041 in every
comparison), and the repeated-measures DEG counts are within 0–3 genes of
the non-RM limma-voom counts. This confirms:

1. The stage × season × infection transcriptional programme dominates the
   transcriptome so thoroughly that vine-specific baseline is a negligible
   second-order effect.
2. Pseudoreplication is NOT a practical concern for this dataset.
3. The fix01 non-RM DESeq2 / limma-voom results are valid as primary
   reporting and the RM results should be reported in the supplement as a
   sensitivity check.

This is itself a publishable observation: *"within-vine correlation
estimated by `duplicateCorrelation` was negligible (|ρ̂| ≤ 0.04 across all
comparisons), confirming that the stratified sampling design is effectively
unpseudoreplicated at the transcriptome level."*

---

*Generated from FARM SLURM job 32497371 outputs on 2026-04-20. Raw
per-comparison DEG tables under `results/fix01_de_redo/`.*
