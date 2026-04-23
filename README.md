<p align="center">
  <img src="assets/banner.png" alt="GRBV Transcriptomics Analysis" width="100%">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Samples-327-1E2761?style=for-the-badge" alt="327 samples">
  <img src="https://img.shields.io/badge/Genes-14,604-B2182B?style=for-the-badge" alt="14604 genes">
  <img src="https://img.shields.io/badge/Seasons-2021%20%2B%202022-009E73?style=for-the-badge" alt="Two seasons">
  <img src="https://img.shields.io/badge/Pipeline-Snakemake%20%2B%20R-CC79A7?style=for-the-badge" alt="Snakemake + R">
  <img src="https://img.shields.io/badge/Stage-Pre--submission-D55E00?style=for-the-badge" alt="Pre-submission">
</p>

---

<div align="center">

*A 327-sample RNA-Seq study asking how Grapevine Red Blotch Virus rewires the* ***Vitis vinifera*** *(Merlot)*
<br/>
*transcriptome, how temperature modulates that rewiring, and which host genes carry*
<br/>
*the causal weight between the two.*

</div>

---

## The story in one diagram

```mermaid
flowchart TD
    classDef env fill:#E69F00,stroke:#333,stroke-width:1px,color:#fff
    classDef vir fill:#0072B2,stroke:#333,stroke-width:1px,color:#fff
    classDef host fill:#8E44AD,stroke:#333,stroke-width:1px,color:#fff
    classDef down fill:#2166AC,stroke:#333,stroke-width:1px,color:#fff
    classDef phe fill:#B2182B,stroke:#333,stroke-width:1px,color:#fff

    T["Hot season<br/>+7.6 degC at Harvest"]:::env
    V["Viral titer<br/>biphasic<br/>3.6x early, 0.18x late"]:::vir
    M["MBF1c<br/>HSF co-activator<br/>61.8 percent mediated"]:::host
    S["29 of 36 sHSPs<br/>suppressed<br/>module score -55 percent"]:::down
    I["Canonical immunity<br/>10 categories<br/>|log2FC| less than 0.08"]:::down
    P["Disease phenotype<br/>red blotch<br/>sugar deficit"]:::phe

    T -->|path a| V
    V -->|path b| M
    M --> S
    T -. direct effect .-> S
    V -. no trigger .-> I
    S --> P
    I -.->|silent| P

    linkStyle 0 stroke:#E69F00,stroke-width:2px
    linkStyle 1 stroke:#0072B2,stroke-width:2px
    linkStyle 2 stroke:#8E44AD,stroke-width:2px
```

<blockquote align="center">
<b>The headline is not that GRBV is fighting the plant.</b><br/>
<b>It is that GRBV is quieting one specific arm of the plant's stress biology,<br/>
the small heat-shock protein axis, while slipping past canonical immune defense<br/>
without triggering it at all.</b>
</blockquote>

---

## Study at a glance

| | |
|---|---|
| Species | *Vitis vinifera* cv. Merlot, one commercial vineyard, 10+ blocks |
| Samples | 327 (264 RB+, 63 RB-) |
| Seasons | 2021 (cool), 2022 (hot, +7.6 °C at Harvest) |
| Maturity stages | Pre-véraison, Véraison, Post-véraison, Harvest |
| Infection cohorts | 7 onset years (2016–2022), n = 105–214 DEGs each |
| Sequencing | TAG-seq (UMI-indexed) · STAR · Pinot Noir ASM3070453v1 |
| Primary DE design | DESeq2 Wald, `~ Maturity_stage + Sampling_Season + Infection_Group`, BH-FDR < 0.05, \|log₂FC\| ≥ 1 |
| Figures | 10 Main composites + 47 Supplementary, PNG + SVG |

<p align="center">
  <img src="figures/main/Main03_inverted_trajectory.png" alt="Main03 — Inverted seasonal DEG trajectory" width="95%">
</p>

<p align="center"><i>Main Figure 3. The inverted seasonal DEG trajectory. 2021 peaks at Harvest (1,567 DEGs); 2022 peaks earlier at Post-véraison (772). Viral titer runs 3.6x higher in the hot year at Pre-véraison and 0.18x at Harvest: the biphasic titer pattern is the mechanistic hinge.</i></p>

---

## Three threads, one story

<table>
<tr><td align="center" width="50"><b>1</b></td><td>

### Inverted seasonal trajectory *(Q3, Q4)*

2021 peaks at **Harvest with 1,567 DEGs**. 2022 peaks earlier at **Post-véraison with 772**, falling back to 759 by Harvest. Cross-season fold-change correlation rises with maturity: r = 0.34 at Pre-véraison, **r = 0.67 at Post-véraison**. The direction switch is mechanistically driven by biphasic viral titer: higher early in hot years, lower late, flipping the timing of disease expression.

</td></tr>
<tr><td align="center"><b>2</b></td><td>

### Thermal defense silenced, canonical immunity untouched

**29 of 36 small heat-shock proteins** (sHSPs) are down-regulated in infected vines. HSF ↔ sHSP coordination stays coherent (**r = 0.625, p = 7.2e-37**) but the sHSP module score drops from +1.64 in healthy vines to +0.74 in infected (**55 % reduction**). All 10 canonical plant-immunity categories (NLR, PRR, PR, SA, JA, ET, WRKY, callose, ROS, R-adapters) show \|log₂FC\| < 0.08. GRBV is not triggering effector-triggered or pattern-triggered immunity. It is dampening thermal defense.

</td></tr>
<tr><td align="center"><b>3</b></td><td>

### MBF1c, SPS1, and 37 other mediators

Causal-mediation sweep on the top 2,000 variable genes returns **39 temperature → titer → gene mediators** at joint BH-FDR < 0.05. The top hit, LOC100249249, turns out to be **multiprotein-bridging factor 1c (MBF1c)**, a known HSF co-activator: joint p = 4.8 × 10⁻¹¹, **61.8 % of the temperature effect mediated by viral titer**. Second-tier robust mediator is LOC100232974 = **sucrose-phosphate synthase 1 (SPS1)**. Paired harvest biochemistry (Figs S61-S63, n = 33) shows SPS1 correlates *negatively* with Brix within-season (2021 ρ = −0.73), so SPS1 is up-regulated in the infected, low-Brix berries. The mediator status stands; the role is compensatory, not rate-limiting.

</td></tr>
<tr><td align="center"><b>4</b></td><td>

### Harvest juice chemistry mirrors the inverted seasonal trajectory

Two-way Season × Infection ANOVA on n = 33 paired harvest samples confirms the transcriptomic inversion in the phenotype. **TA** is the cleanest mirror: no main effect of season (p = 0.69) or infection (p = 0.16), but the Season × Infection interaction reaches **η² = 0.77, p = 4.3 × 10⁻¹¹**. GRBV raises TA by **+0.62 g/L in 2021 cool** and lowers it by **−0.79 g/L in 2022 hot**; malic acid follows the same inversion (+398 vs −195 mg/L). Brix is suppressed in both seasons but **2.6× stronger in 2022** (−9.83 vs −3.77 °Bx). In 2021 a biphasic cohort-age recovery appears at harvest: Brix rises +0.10 °Bx / yr with years-since-infection (r = 0.76, p = 0.0045) while nitrogen markers decline.

</td></tr>
</table>

---

## What changed under the reanalysis

```mermaid
flowchart LR
    classDef old  fill:#fce8e6,stroke:#B2182B,color:#5a0000
    classDef fix  fill:#fff4d8,stroke:#E69F00,color:#5a3d00
    classDef new  fill:#e3f0fb,stroke:#2166AC,color:#0a2b52

    subgraph pre ["Pre-fix (Welch-t)"]
        P1["908 universal DEGs<br/>pooled, uncorrected"]:::old
        P2["AUC = 1.000<br/>no host/viral split"]:::old
        P3["0 mediators<br/>at tv = 500"]:::old
        P4["27/27 QTL overlap<br/>no null distribution"]:::old
    end

    subgraph fix ["Critical review: 13 fixes"]
        F1["Covariate-adjusted DESeq2<br/>Stage + Season + Infection"]:::fix
        F2["Leave-one-season-out CV<br/>host-only features"]:::fix
        F3["Raise tv to 2,000"]:::fix
        F4["10,000-permutation null<br/>chromosome-stratified"]:::fix
    end

    subgraph post ["Post-fix"]
        Q1["126 core infection DEGs<br/>13.6 percent of the original"]:::new
        Q2["Host LOSO AUC 0.73-0.88<br/>honest generalisation"]:::new
        Q3["39 mediators<br/>MBF1c top p=4.8e-11"]:::new
        Q4["fold = 2.01, p less than 0.0002<br/>defensible enrichment"]:::new
    end

    P1 --> F1 --> Q1
    P2 --> F2 --> Q2
    P3 --> F3 --> Q3
    P4 --> F4 --> Q4
```

<blockquote align="center">
<i>The data layer came out of the reanalysis richer (24 covariate-adjusted contrasts,</i>
<br/>
<i>including per-stage x per-season cells). The manuscript claims came out honest.</i>
</blockquote>

Full documentation: <a href="docs/CRITICAL_REVIEW.md">CRITICAL_REVIEW.md</a>, <a href="docs/DE_REANALYSIS_REPORT.md">DE_REANALYSIS_REPORT.md</a>, <a href="docs/SENSITIVITY_REPORT.md">SENSITIVITY_REPORT.md</a>, <a href="docs/METHODS_ADDENDUM.md">METHODS_ADDENDUM.md</a>.

---

## Analysis pipeline

```mermaid
flowchart LR
    classDef raw    fill:#f4f4f4,stroke:#64748B,color:#0F172A
    classDef pipe   fill:#d6e3f5,stroke:#2166AC,color:#0a2b52
    classDef out    fill:#e3dcf0,stroke:#8E44AD,color:#2a164a

    subgraph raw ["Raw data"]
        direction TB
        A1[TAG-seq reads<br/>UMI-indexed]:::raw
        A2[Sample metadata<br/>vine_id, stage, season]:::raw
        A3[Viral titer<br/>RT-qPCR + dPCR]:::raw
        A4[Field weather<br/>daily 2021, 2022]:::raw
    end

    subgraph pipe ["HPC pipeline (Snakemake on FARM)"]
        direction TB
        B1[STAR + featureCounts<br/>Pinot Noir ASM3070453v1]:::pipe
        B2[DESeq2<br/>24 stratified contrasts]:::pipe
        B3[Causal mediation<br/>tv = 2,000]:::pipe
        B4[WGCNA<br/>13 co-expression modules]:::pipe
        B5[GSEA<br/>375 sets x 8 conditions]:::pipe
        B6[QTL permutation<br/>10,000 chromosome-stratified]:::pipe
        B7[NCBI annotation<br/>82 of 86 LOC resolved]:::pipe
    end

    subgraph out ["Deliverables"]
        direction TB
        C1[Figure atlas<br/>10 Main + 47 Supp<br/>PNG + SVG]:::out
        C2[Manuscript draft<br/>docs/MANUSCRIPT_DRAFT_v1.pdf]:::out
        C3[Funding deck<br/>48 slides, 19 MB]:::out
        C4[Methods + sensitivity<br/>docs/]:::out
    end

    A1 --> B1 --> B2
    A2 --> B2
    A3 --> B3
    A4 --> B3
    B2 --> B3
    B2 --> B4
    B2 --> B5
    B2 --> B6
    B4 --> B7
    B2 --> C1
    B3 --> C1
    B4 --> C1
    B5 --> C1
    C1 --> C2
    C1 --> C3
    B2 --> C4
```

---

## Repository layout

<details>
<summary><b>Click to expand full tree</b></summary>

```
.
├── figures/
│   ├── main/               Main01-Main10 composite figures (PNG + SVG)
│   ├── supp/               S02-S55 supplementary figures (PNG + SVG)
│   ├── scripts_R/          R / ggplot2 scripts for every figure
│   │   └── theme_grbv.R    Shared theme + palette + save_grbv() helper
│   ├── GRBV_Main_Figures_Atlas.pdf       10-page main atlas
│   └── GRBV_Full_Figure_Atlas.pdf        65-page full atlas
├── farm_pipeline/
│   ├── Snakefile           13 rules, 24 DE contrasts, sensitivity sweeps
│   ├── scripts/
│   ├── results/            fix01, fix01c, fix03, fix10, fix11 + rules 00-10
│   ├── data/, resources/
│   └── docs/
│       └── LITERATURE_SCHEMA.md
├── docs/                   Project-level deliverables (14 files)
│   ├── MANUSCRIPT_DRAFT_v1.{md,pdf}
│   ├── INTEGRATED_ANALYSIS.{md,pdf}
│   ├── CRITICAL_REVIEW.{md,pdf}
│   ├── DE_REANALYSIS_REPORT.md
│   ├── METHODS_ADDENDUM.md
│   ├── HSFA2_VIGS_EXPERIMENT_PLAN.md
│   ├── SENSITIVITY_REPORT.{md,pdf}
│   ├── EXTENDED_RESULTS_REPORT.{md,pdf}
│   └── PLAIN_LANGUAGE_SUMMARY.md
├── Old-Analysis/           Pre-fix era archive
│   ├── images/             55 legacy PNGs from Python pipeline
│   ├── scripts_python/     19 pre-fix analysis scripts
│   ├── generate_scripts/   5 pre-fix PDF / PPTX builders
│   ├── reports/            8 legacy deliverables
│   ├── superseded_csvs/    11 CSVs with zero current-R-script references
│   └── farm_sensitivity_archive/        58 MB raw threshold sweeps
├── GRBV_Funding_Committee_Presentation.pptx     48-slide committee deck
├── counts_matrix.txt + sample_metadata.txt + 19 active-pipeline CSVs
├── org.Vvinifera.eg.db/, go-basic.obo
├── CLAUDE.md               Project memory / session context
└── README.md               This file
```

</details>

---

## Running it

### Figures

Every figure is an R / ggplot2 script in [`figures/scripts_R/`](figures/scripts_R/). Each script sources [`theme_grbv.R`](figures/scripts_R/theme_grbv.R) (project palette + save helper) and writes **PNG + SVG** to `figures/main/` or `figures/supp/`.

```bash
# One figure
Rscript figures/scripts_R/main03_inverted_trajectory.R

# A whole supplementary batch (S02 through S07)
Rscript figures/scripts_R/supp_de_overview.R

# Rebuild the atlases from the current PNG set
python figures/scripts_R/build_atlas.py
python figures/scripts_R/build_full_atlas.py
```

<details>
<summary><b>Output formats</b></summary>

- **PNG**: 600 dpi, via `ragg::agg_png`, Arial fonts reliable on macOS
- **SVG**: vector, via `svglite::svglite`, opens cleanly in Illustrator, Inkscape, Figma
- **PDF**: opt-in with `save_grbv(..., formats = c("png","pdf"))` via `cairo_pdf`
- **Atlases**: 2 PDFs assembled directly from PNGs (reportlab)

</details>

### FARM HPC pipeline

The extended analysis runs as a Snakemake pipeline on UC Davis FARM.

```bash
cd farm_pipeline
bash setup_farm.sh      # install R packages + annotation resources
bash run_farm.sh        # submit to partition low, account publicgrp
bash pull.sh            # rsync results back to local figures/
```

### Funding-committee deck

Generated by [pptxgenjs](https://gitbrent.github.io/PptxGenJS/) via a Node.js script. Output is 48 slides, 19 MB, with 10 Main + 6 Supp figures embedded and a Q-tag in the top-right of every result slide showing which research question each figure answers.

---

## Deliverables

<div align="center">

| Asset | Location | Size |
|---|---|---|
| Manuscript draft (with figures embedded) | [`docs/MANUSCRIPT_DRAFT_v1.pdf`](docs/MANUSCRIPT_DRAFT_v1.pdf) | 9.6 MB |
| Main figures atlas | [`figures/GRBV_Main_Figures_Atlas.pdf`](figures/GRBV_Main_Figures_Atlas.pdf) | 11 MB |
| Full atlas (main + supp) | [`figures/GRBV_Full_Figure_Atlas.pdf`](figures/GRBV_Full_Figure_Atlas.pdf) | 36 MB |
| Methods addendum | [`docs/METHODS_ADDENDUM.md`](docs/METHODS_ADDENDUM.md) | 16 KB |
| Committee deck | [`GRBV_Funding_Committee_Presentation.pptx`](GRBV_Funding_Committee_Presentation.pptx) | 19 MB |
| Wet-lab plan (VIGS) | [`docs/HSFA2_VIGS_EXPERIMENT_PLAN.md`](docs/HSFA2_VIGS_EXPERIMENT_PLAN.md) | 17 KB |
| Critical review + fix log | [`docs/CRITICAL_REVIEW.pdf`](docs/CRITICAL_REVIEW.pdf) | 23 KB |
| Integrated analysis | [`docs/INTEGRATED_ANALYSIS.pdf`](docs/INTEGRATED_ANALYSIS.pdf) | 28 KB |

</div>

---

## What's next

```mermaid
timeline
    title Submission roadmap
    section Validation
        May-Jun 2026 : qPCR on 10 key genes : MBF1c, SPS1, 3 sHSPs, HSFA2, 3 biomarkers
        Jun-Jul 2026 : GEO + SRA deposition : 14,604 genes x 327 samples
    section Submission
        Aug 2026 : Nature Plants submission : methods + data finalised
    section Extension
        Jul 2026 - Feb 2027 : VIGS on MBF1c + SPS1 : N. benthamiana then grapevine
        Aug 2026 : Tier-1 meta-analysis : 5 published Vitis stress studies
```

---

## Acknowledgments

<table>
<tr><td width="30%" align="right" valign="top"><b>In memoriam</b></td>
<td>
The founding PI, <b>Dr. Anita Oberholster</b> (UC Davis Viticulture and Enology, 1980 to 2025), passed away in January 2025. Her lab built the field-sampling design, ran the wet-lab, and framed the year-one manuscript. Every result here sits on that foundation.
</td></tr>
<tr><td align="right" valign="top"><b>Mentorship</b></td>
<td><b>Dr. David E. Block</b>, UC Davis V and E. Current PI of record; manuscript strategy; resource allocation on FARM.</td></tr>
<tr><td align="right" valign="top"><b>Virology co-authorship</b></td>
<td><b>Dr. Mysore R. Sudarshana</b>, USDA-ARS and UC Davis Plant Pathology. Field-site liaison; viral-titer wet-lab; co-author on two in-preparation manuscripts.</td></tr>
<tr><td align="right" valign="top"><b>RNA-Seq pipeline</b></td>
<td><b>Bradley N. Jenner</b>, UC Davis Bioinformatics Core (2025-04-15 to 2025-05-12). TAG-seq QC, HTStream, umi_tools dedup, STAR, featureCounts.</td></tr>
<tr><td align="right" valign="top"><b>Funding</b></td>
<td>USDA-ARS grapevine pathology programs; California Department of Food and Agriculture (CDFA) specialty-crop block grants.</td></tr>
<tr><td align="right" valign="top"><b>Compute</b></td>
<td>UC Davis FARM HPC, partition <code>low</code>, account <code>publicgrp</code>.</td></tr>
<tr><td align="right" valign="top"><b>Annotation</b></td>
<td>NCBI; <code>org.Vvinifera.eg.db</code>; Gene Ontology Consortium; KEGG.</td></tr>
</table>

---

## Citing

<blockquote>
Singh, P.P., Block, D.E., Sudarshana, M.R., Oberholster, A. (in preparation).<br/>
<i>Biphasic temperature drives viral titer dynamics and silences thermal defense in</i> Vitis vinifera <i>cv. Merlot under Grapevine Red Blotch Virus infection: a 327-sample transcriptomic study.</i>
</blockquote>

Contact: [ppssingh@ucdavis.edu](mailto:ppssingh@ucdavis.edu) · ORCID [0000-0001-7921-9379](https://orcid.org/0000-0001-7921-9379) · [prempsingh.com](https://prempsingh.com)

---

<p align="center"><sub>MIT License · see <a href="LICENSE">LICENSE</a> for details</sub></p>
