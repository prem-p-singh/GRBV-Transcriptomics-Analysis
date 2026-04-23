<p align="center">
  <img src="assets/banner.png" alt="GRBV Transcriptomics Analysis" width="100%">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Samples-327-1E2761?style=for-the-badge" alt="327 samples">
  <img src="https://img.shields.io/badge/Genes-14,604-B2182B?style=for-the-badge" alt="14604 genes">
  <img src="https://img.shields.io/badge/Seasons-2021%20%2B%202022-009E73?style=for-the-badge" alt="Two seasons">
  <img src="https://img.shields.io/badge/Pipeline-Snakemake%20%2B%20R-CC79A7?style=for-the-badge" alt="Snakemake + R">
  <img src="https://img.shields.io/badge/Status-Pre--submission%20(Nature%20Plants)-D55E00?style=for-the-badge" alt="Pre-submission">
</p>

---

> A 327-sample RNA-Seq study asking how Grapevine Red Blotch Virus (GRBV) rewires the *Vitis vinifera* (Merlot) transcriptome, how temperature modulates that rewiring, and which host genes carry the causal weight between the two.

Two consecutive seasons in the same vineyard turned this project into a natural experiment: 2021 was cool and typical; 2022 was hot, with a +7.6 °C mean excess at Harvest. The disease pattern flipped. Cool years peaked at Harvest, hot years peaked at Post-véraison, and once we looked hard enough to see why, the biology got specific.

The headline is not that GRBV is fighting the plant. It is that GRBV is quieting one specific arm of the plant's stress biology, the small heat-shock protein axis, while slipping past canonical immune defense without triggering it at all.

---

## Study at a glance

| | |
|---|---|
| Species | *Vitis vinifera* cv. Merlot, one commercial vineyard, 10+ blocks |
| Samples | 327 (264 RB+, 63 RB-) |
| Seasons | 2021 (cool), 2022 (hot, +7.6 °C at Harvest) |
| Maturity stages | Pre-véraison, Véraison, Post-véraison, Harvest |
| Infection cohorts | 7 onset years (2016 – 2022), n = 105 – 214 DEGs each |
| Sequencing | TAG-seq (UMI-indexed) · STAR · Pinot Noir reference ASM3070453v1 |
| Primary DE design | DESeq2 Wald, `~ Maturity_stage + Sampling_Season + Infection_Group`, BH-FDR < 0.05, \|log₂FC\| ≥ 1 |
| Figures | 10 Main composites + 47 Supplementary, PNG + SVG |

<p align="center">
  <img src="figures/main/Main03_inverted_trajectory.png" alt="Main03 — Inverted seasonal DEG trajectory" width="95%">
</p>

*Main Figure 3. The inverted seasonal DEG trajectory. 2021 peaks at Harvest (1,567 DEGs); 2022 peaks earlier at Post-véraison (772). Viral titer runs 3.6× higher in the hot year at Pre-véraison and 0.18× at Harvest: the biphasic titer pattern is the mechanistic hinge behind the inverted transcriptome.*

---

## Three threads, one story

### 1. Inverted seasonal trajectory (Q3, Q4)

2021 DEG counts peak at Harvest with 1,567 genes. 2022 peaks earlier at Post-véraison with 772, falling back to 759 by Harvest. Cross-season fold-change correlation rises with maturity: r = 0.34 at Pre-véraison, r = 0.67 at Post-véraison. The direction switch is not seasonal noise. It is mechanistically driven by biphasic viral titer: higher titer early in hot years, lower titer late, flipping the timing of disease expression.

### 2. Thermal defense silenced, canonical immunity untouched (sHSP axis)

29 of 36 small heat-shock proteins (sHSPs) are down-regulated in infected vines, the largest single-category signal in the 115-gene heat-responsive set. HSF ↔ sHSP coordination stays coherent within each infection group (r = 0.625, p = 7.2 × 10⁻³⁷) but the sHSP module score drops from +1.64 in healthy vines to +0.74 in infected, a 55% reduction. Meanwhile, every single one of 10 canonical plant-immunity categories (NLR, PRR, PR, SA, JA, ET, WRKY, callose, ROS, R-adapters) shows \|log₂FC\| < 0.08. GRBV is not triggering effector-triggered or pattern-triggered immunity. It is dampening thermal defense.

### 3. MBF1c, SPS1, and 37 other mediators (Q4 mechanism)

Causal-mediation sweep on the top 2,000 most variable genes returns 39 temperature → titer → gene mediators at joint BH-FDR < 0.05. The top hit, LOC100249249, turns out to be multiprotein-bridging factor 1c (MBF1c), a known HSF co-activator: joint p = 4.8 × 10⁻¹¹, 61.8% of the temperature effect mediated by viral titer. The second-tier robust mediator is LOC100232974, sucrose-phosphate synthase 1 (SPS1), the enzyme that makes sucrose: a direct molecular link from viral titer to the sugar-deficit phenotype.

---

## What changed under the reanalysis

Before submission, the project ran its own reviewer-style critical review (`docs/CRITICAL_REVIEW.md`). Thirteen methodological issues surfaced; all thirteen were addressed. Four of the corrections meaningfully moved published claims:

| Claim | Pre-fix | Post-fix |
|---|---|---|
| "Universal DEGs" | 908 from Welch-t, no covariate adjustment | 126 under DESeq2 with `~ Maturity_stage + Sampling_Season + Infection_Group`. Only 13.6% of the original set replicated. |
| ML biomarker AUC | 1.000 across 5 methods | 1.000 only when viral N761 transcripts are in the feature pool. Honest host-only LOSO AUC is 0.73 – 0.88. |
| Causal mediators | Zero significant at tv = 500 | 39 at tv = 2,000, with MBF1c as the cross-pipeline robust top hit. |
| QTL overlap | 27/27 inside 1 Mb of a QTL, no null | Same 27/27, now with a 10,000-permutation null: fold = 2.01, empirical p < 0.0002. |

The fix pass is documented in `docs/DE_REANALYSIS_REPORT.md`, `docs/SENSITIVITY_REPORT.md`, and `docs/METHODS_ADDENDUM.md`. The data layer came out of it richer (24 covariate-adjusted contrasts, including per-stage × per-season Q2x); the manuscript claims came out honest.

---

## Repository layout

```
.
├── figures/
│   ├── main/              Main01–Main10 composite figures, PNG + SVG
│   ├── supp/              S02–S55 supplementary figures, PNG + SVG
│   ├── scripts_R/         R/ggplot2 scripts for every figure
│   ├── GRBV_Main_Figures_Atlas.pdf     10-page main atlas
│   └── GRBV_Full_Figure_Atlas.pdf      65-page full atlas
├── farm_pipeline/
│   ├── Snakefile          13 rules, 24 DE contrasts, sensitivity sweeps
│   ├── scripts/
│   ├── results/           fix01, fix01c, fix03, fix10, fix11 + rules 00–10
│   ├── data/, resources/
│   └── docs/LITERATURE_SCHEMA.md
├── docs/
│   ├── MANUSCRIPT_DRAFT_v1.{md,pdf}
│   ├── INTEGRATED_ANALYSIS.{md,pdf}
│   ├── CRITICAL_REVIEW.{md,pdf}
│   ├── DE_REANALYSIS_REPORT.md
│   ├── METHODS_ADDENDUM.md
│   ├── HSFA2_VIGS_EXPERIMENT_PLAN.md
│   ├── SENSITIVITY_REPORT.{md,pdf}
│   ├── EXTENDED_RESULTS_REPORT.{md,pdf}
│   └── PLAIN_LANGUAGE_SUMMARY.md
├── Old-Analysis/          Pre-fix archive (55 old PNGs, 19 Python scripts,
│                          8 legacy deliverables, 11 superseded CSVs)
├── GRBV_Funding_Committee_Presentation.pptx   48-slide committee deck
├── counts_matrix.txt + sample_metadata.txt + 19 active-pipeline CSVs
├── org.Vvinifera.eg.db/, go-basic.obo
├── CLAUDE.md              Project memory / session context
└── README.md              This file
```

---

## Running it

### Figures

Every figure is an R/ggplot2 script in `figures/scripts_R/`. Each script is self-contained: it sources `theme_grbv.R` (the project palette + save helper) and writes a PNG + SVG pair to `figures/main/` or `figures/supp/`.

```bash
# One figure
Rscript figures/scripts_R/main03_inverted_trajectory.R

# A whole supplementary batch
Rscript figures/scripts_R/supp_de_overview.R          # S02 through S07

# Rebuild the atlases from the current PNG set
python figures/scripts_R/build_atlas.py               # 10-page main atlas
python figures/scripts_R/build_full_atlas.py          # 65-page full atlas
```

Figure exports default to **PNG (600 dpi via ragg)** and **SVG (via svglite)**. SVG opens cleanly in Illustrator, Inkscape, or Figma if you want to hand-tweak a panel. PDF is an opt-in: `save_grbv(..., formats = c("png","pdf"))`.

### FARM HPC pipeline

The extended analysis (DESeq2 stratified contrasts, mediation, WGCNA, GSEA, QTL permutation, NCBI annotation) runs as a Snakemake pipeline on UC Davis FARM.

```bash
cd farm_pipeline
bash setup_farm.sh         # install R packages + annotation resources
bash run_farm.sh           # submit to partition low, account publicgrp
bash pull.sh               # rsync results back to local figures/
```

### Funding-committee deck

Generated by `pptxgenjs` via a Node.js script. Output is 48 slides, 19 MB, with 10 Main + 6 Supp figures embedded and a Q-tag in the top-right of every result slide showing which research question each figure answers.

---

## Deliverables

| | |
|---|---|
| Manuscript draft | `docs/MANUSCRIPT_DRAFT_v1.pdf` (9.6 MB, 10 figures embedded) |
| Figure atlases | `figures/GRBV_Main_Figures_Atlas.pdf`, `GRBV_Full_Figure_Atlas.pdf` |
| Methods addendum | `docs/METHODS_ADDENDUM.md` (reviewer-facing) |
| Committee deck | `GRBV_Funding_Committee_Presentation.pptx` |
| Wet-lab protocol | `docs/HSFA2_VIGS_EXPERIMENT_PLAN.md` (7-month VIGS plan) |
| Integrated analysis | `docs/INTEGRATED_ANALYSIS.pdf` |

## What's next

1. qPCR validation of 10 key genes: MBF1c, SPS1, LOC100266309, 3 sHSPs, HSFA2, and 3 top host biomarkers. Six to eight weeks wet-lab.
2. GEO + SRA data deposition. Two to four weeks.
3. Nature Plants submission. Ready after (1) and (2).
4. VIGS agroinoculation of MBF1c and SPS1 in *N. benthamiana* and grapevine, per `docs/HSFA2_VIGS_EXPERIMENT_PLAN.md`. Seven months, elevates the bid to PNAS-competitive.
5. Tier-1 literature meta-analysis across 5 published Vitis stress studies (GLRaV-3, GPGV, powdery mildew, downy mildew, heat). Schema already scaffolded at `farm_pipeline/docs/LITERATURE_SCHEMA.md`.

---

## Acknowledgments

The founding PI of this study, **Dr. Anita Oberholster** (UC Davis Viticulture & Enology, 1980 – 2025), passed away in January 2025. Her lab built the field sampling design, ran the wet-lab work, and framed the year-one manuscript. Every result in this repository sits on that foundation.

Current mentorship from **Dr. David E. Block** (UC Davis V&E). GRBV virology co-authorship and field-site liaison from **Dr. Mysore R. Sudarshana** (USDA-ARS / UC Davis Plant Pathology). RNA-Seq pipeline execution by **Bradley N. Jenner** (UC Davis Bioinformatics Core, 2025-04-15 to 2025-05-12).

Funding: USDA-ARS grapevine pathology programs; California Department of Food and Agriculture (CDFA) specialty-crop block grants.

Compute: UC Davis FARM HPC, partition `low`, account `publicgrp`.

Annotation resources: NCBI, `org.Vvinifera.eg.db`, Gene Ontology Consortium, KEGG.

---

## Citing

> Singh, P.P., Block, D.E., Sudarshana, M.R., Oberholster, A. (in preparation). *Biphasic temperature drives viral titer dynamics and silences thermal defense in* Vitis vinifera *cv. Merlot under Grapevine Red Blotch Virus infection: a 327-sample transcriptomic study.*

Contact: [ppssingh@ucdavis.edu](mailto:ppssingh@ucdavis.edu) · ORCID [0000-0001-7921-9379](https://orcid.org/0000-0001-7921-9379) · [prempsingh.com](https://prempsingh.com)

## License

MIT. See [LICENSE](LICENSE).
