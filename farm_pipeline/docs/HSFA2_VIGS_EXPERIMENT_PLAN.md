# HSFA2 VIGS + temperature challenge — functional validation protocol

**Study:** Causal test of the HSF → sHSP → antiviral axis in GRBV–*V. vinifera* cv. Merlot.
**Author:** Prem Pratap Singh
**Version:** 1.0 (2026-04-19)

---

## 1. Objective

Test whether silencing *HSFA2* in *V. vinifera* cv. Merlot (i) reduces small heat-shock protein (sHSP) transcript accumulation, (ii) increases susceptibility to Grapevine Red Blotch Virus (GRBV) at baseline, and (iii) **abolishes the temperature-dependent reduction in viral titer** predicted by the transcriptomic mediation model (Temperature → HSF/sHSP axis → viral titer).

## 2. Primary hypotheses

| ID | Hypothesis | Predicted outcome if true |
|---|---|---|
| **H1** | HSFA2 drives sHSP expression in grapevine. | VIGS-HSFA2 plants show ≥ 50 % reduction in HSFA2 transcript and significant downregulation of ≥ 3 representative sHSPs vs TRV–empty controls. |
| **H2** | sHSPs suppress GRBV replication. | At baseline temperature, VIGS-HSFA2 plants carry significantly higher viral titer than TRV-empty controls at 8 and 12 weeks post-inoculation (wpi). |
| **H3** | Temperature's antiviral effect is mediated through the HSF–sHSP axis. | Control plants show lower titer at elevated temperature; VIGS-HSFA2 plants show no significant titer difference across temperatures — i.e., the interaction term *VIGS × temperature* is significant in a two-way ANOVA. |

The interaction test in **H3** is the single most discriminating analysis for PNAS-level causal language: it directly contrasts the two biological routes through which temperature could influence titer.

## 3. Background & rationale

The GRBV RNA-seq study (327 samples, 2 seasons, 4 maturity stages) identified three convergent findings that motivate this experiment:

1. **sHSP collapse**: 29 of 36 expressed sHSPs are downregulated in RB+ vines; module score drops ~55 % (+1.64 → +0.74).
2. **HSF–sHSP axis**: Pearson r = 0.625, p = 7.2 × 10⁻³⁷ across samples; expression of HSF family members covaries with the sHSP module.
3. **Temperature mediation** (sensitivity-corrected). A multi-threshold HIMA sweep (top-var = 3k, 5k, 7.5k, 10k, 14.6k) identified a panel of robust mediators rather than a single gene: `LOC100266309` (passes at 3/5 thresholds), `LOC100232974` (2/5), and `LOC100249249` (2/5). Independently, the main paper's causal-mediation analysis (Temperature → Titer → Gene) at tv = 2k and 5k identifies `LOC100249249` as the top-ranked mediator (p = 4.8 × 10⁻¹¹ and 1.2 × 10⁻¹⁰). **`LOC100249249` is the single gene to appear as a robust mediator in both independent causal pipelines at multiple thresholds, making it the primary functional-validation target.** See [SENSITIVITY_REPORT.md](SENSITIVITY_REPORT.md) for the full sweep.

Points 1–2 are correlative; point 3 is inferred from observational data. A perturbation experiment is required to establish causality — exactly the gap PNAS reviewers will flag.

## 4. Experimental design

### 4.1 Factors and groups

Factor A – construct (3 levels):
- `VIGS-HSFA2` (antisense/sense insert targeting grapevine HSFA2 CDS, pTRV2 backbone)
- `TRV-empty` (pTRV2 carrying a neutral stuffer; construct control)
- `Mock` (agro-infiltration buffer without T-DNA; procedure control)

Factor B – temperature regime (2 levels):
- `Control` — 22 °C day / 18 °C night, 16 h light
- `Elevated` — 32 °C day / 24 °C night, 16 h light (reflects 2022 Harvest-stage delta of +7.6 °C from the transcriptomic dataset)

Factor C – viral inoculation (2 levels): `GRBV` (inoculated) or `mock` (no virus) to establish baseline.

Full design is 3 × 2 × 2 = 12 cells. Cell composition:

| Construct | Temperature | Virus | n (biological replicates) |
|---|---|---|---|
| VIGS-HSFA2 | Control | GRBV | 8 |
| VIGS-HSFA2 | Control | mock | 4 |
| VIGS-HSFA2 | Elevated | GRBV | 8 |
| VIGS-HSFA2 | Elevated | mock | 4 |
| TRV-empty | Control | GRBV | 8 |
| TRV-empty | Control | mock | 4 |
| TRV-empty | Elevated | GRBV | 8 |
| TRV-empty | Elevated | mock | 4 |
| Mock | Control | GRBV | 4 |
| Mock | Control | mock | 2 |
| Mock | Elevated | GRBV | 4 |
| Mock | Elevated | mock | 2 |

**Total: n = 64 plants.** Each n is one independently propagated seedling. Three independent experimental rounds (≥ n = 24 per round) spaced ≥ 2 weeks apart to guard against batch effects.

### 4.2 Power calculation

Primary endpoint: log₁₀(GRBV titer) at 8 wpi.
Assumed SD of log₁₀(titer) from the observational dataset (263 RB+ samples) ≈ 0.42.
Minimum detectable effect: Δ log₁₀(titer) = 0.5 (≈ 3-fold titer change).
Two-sided *t*-test, α = 0.05, power = 0.80 → n = 6 per group. Using n = 8 per titer-comparison cell provides ~92 % power and tolerates one drop-out per cell.

## 5. VIGS construct design

### 5.1 Target

Grapevine *HSFA2*. Use the NCBI reference *V. vinifera* genome (GCF_030704535.1, already downloaded to `resources/genome/` on FARM) to identify the gene. Primary locus: the *V. vinifera* HSFA2 ortholog (expect a LOC-prefixed gene symbol). Verify by alignment to *Arabidopsis* AT2G26150 (HSFA2) coding sequence.

### 5.2 Insert

- **Length:** 300–450 bp fragment of the *HSFA2* CDS.
- **Position:** select a fragment in the first 60 % of the CDS, avoiding the conserved DNA-binding domain (bp 1–300 from ATG typically overlaps the DBD). A region 400–750 bp downstream of ATG is ideal.
- **Specificity check:** BLAST the candidate fragment against the full *V. vinifera* transcriptome; require ≥ 3 mismatches to any off-target in any continuous 21 nt window (`siRNA-Finder` / `VIGS tool`, Bio Rochester). Confirm no hits to GRBV sequences.
- **Orientation:** antisense preferred for TRV-based VIGS; include a sense control in a pilot run if HSFA2 is weakly silenced.
- **Cloning:** Gateway or restriction cloning into pTRV2; verify by Sanger sequencing.

### 5.3 Primer template (placeholder – confirm against the reference once locus ID is resolved)

```
HSFA2_VIGS_F  5'-CACC-[20–22 nt forward primer starting ~500 bp 3' of ATG]-3'
HSFA2_VIGS_R  5'-[20–22 nt reverse primer, insert length 300–450 bp]-3'
```

Design rules (lifted from the sgRNA design guide; rules that transfer to RT-primer / insert primer design):
- GC content 45–55 %
- Tm 58–62 °C, matched within 2 °C between forward and reverse
- No self-complementarity > 5 bp; no 3′ hairpins
- Avoid runs of > 4 identical nucleotides
- Verify primer specificity by *in silico* PCR against the reference genome

Once the exact *V. vinifera* HSFA2 locus is resolved from the NCBI GFF on FARM, I can fill in the actual sequences (next follow-up step).

## 6. Plant material, growth, and infection

- **Genotype:** *V. vinifera* cv. Merlot, micropropagated to guarantee uniformity; alternatively field cuttings rooted in perlite under 95 % humidity for 3 weeks.
- **Seedling stage at infiltration:** 4–5 true leaves, ~25 cm height.
- **Agroinfiltration:** pTRV1 + pTRV2-derivative co-infiltrated via needleless syringe into abaxial leaf surfaces. OD₆₀₀ = 0.8 per strain. Infiltrate 2 leaves per plant.
- **GRBV inoculation:** via viruliferous *Spissistilus festinus* (three-cornered alfalfa hopper) placed on plants for 7 days at 14 dpi with TRV — or, if the Yepes et al. infectious clone is available, agroinfiltrate the clone.
- **Temperature transfer:** 21 dpi (once VIGS silencing is established) → assign plants to temperature chambers.
- **Sampling:** leaf discs at 0 (pre-transfer), 4, 8, and 12 wpi (weeks post-GRBV inoculation). Flash-freeze in liquid N₂.

## 7. Readouts

| Endpoint | Method | Primers (placeholder) | Sampling |
|---|---|---|---|
| HSFA2 knockdown | RT-qPCR | *HSFA2*-qPCR-F/R; normalize to *VvUBQ*, *VvActin* | 21 dpi (pre-T shift), 4, 8 wpi |
| sHSP response | RT-qPCR panel of 5 representative sHSPs (choose top-5 downregulated from the transcriptomic dataset) | design from transcript sequences | Same time points |
| GRBV titer (primary) | qPCR of GRBV RNA-dependent RNA polymerase (target gene `LOC100268011` on the host side reflects viral RdRp1 activity; for viral DNA use a GRBV coat protein primer pair) | Published GRBV CP qPCR (e.g. Sudarshana et al. 2015) | 4, 8, 12 wpi |
| Mediator validation (primary) | RT-qPCR of `LOC100249249` — cross-pipeline robust mediator | design de novo | 4, 8 wpi |
| Mediator validation (secondary) | RT-qPCR of `LOC100266309` and `LOC100232974` — HIMA-robust panel | design de novo | 4, 8 wpi |
| Phenotype | Leaf red blotch score (0–4 ordinal scale), plant height | Visual, caliper | 4, 8, 12 wpi |
| sHSP protein | Immunoblot (anti-HSP17.6 commercial antibody, Agrisera AS07 255) | n/a | 8 wpi subset |

All qPCR runs in technical triplicate; calculate ΔΔCt against two reference genes; melting curves and no-RT controls mandatory.

## 8. Statistical analysis plan

Register all analyses **before data collection** (for pre-registration).

1. **H1 (HSFA2 knockdown):** Welch's *t*-test of ΔCt(HSFA2) between VIGS-HSFA2 and TRV-empty at 21 dpi. Report point estimate, 95 % CI, *t*, df, exact *p*.
2. **H1 (sHSP response):** Mixed-effects model `ΔCt ~ construct + (1|plant)` across the 5-gene panel; FDR-correct across panel.
3. **H2 (baseline susceptibility):** Welch's *t*-test of log₁₀(titer) at 8 wpi, VIGS-HSFA2 vs TRV-empty, control-temperature only. Pre-specified one-sided alternative (VIGS > TRV).
4. **H3 (interaction test):** Two-way ANOVA on log₁₀(titer) at 8 wpi: `titer ~ construct × temperature + experiment_round`. Primary inference is the F-test on the `construct × temperature` term. Effect size η² reported.
5. **Mediator validation:** Correlation of HSFA2 ΔCt and each of the three robust mediators (`LOC100249249` primary, `LOC100266309`, `LOC100232974`) ΔCt across all GRBV-inoculated plants (Pearson and Spearman). `LOC100249249` is the pre-specified primary readout because it is the one cross-pipeline robust mediator; the others are secondary.
6. **Effect sizes and CIs** reported alongside every *p*-value.
7. **Missing data:** pairwise deletion; sensitivity analysis with multiple imputation if > 10 % missingness on any endpoint.
8. **Multiplicity:** Bonferroni across the three primary hypotheses H1–H3 (adjusted α = 0.017).

## 9. Predicted outcomes — decision tree

| Outcome | H1 knockdown | H2 baseline titer | H3 interaction | Interpretation |
|---|---|---|---|---|
| **Expected** | ≥ 50 % HSFA2 reduction, sHSP panel down | VIGS > TRV at control temp | Significant interaction (VIGS flattens temperature effect) | Full support for HSF–sHSP antiviral axis; PNAS-level mechanistic claim justified |
| Partial | Knockdown present, sHSPs partial | Trend only, n.s. | Main effect of temp, no interaction | Rewrite to "HSF-sHSP axis contributes to but does not fully mediate the thermal antiviral effect" |
| Null H3 | Knockdown ok | VIGS = TRV | No interaction | Reject the mediation model; revisit alternative pathways (ROS, RNAi, metabolism) |
| Null H1 | No knockdown | — | — | VIGS failed; repeat with sense insert or alternate TRV construct |

## 10. Timeline (single round, all three rounds overlapped)

| Week | Milestone |
|---|---|
| 0 | Design + order VIGS primers (pending locus confirmation) |
| 1–3 | Clone TRV2-HSFA2, sequence-verify |
| 4 | Micropropagate Merlot seedlings (30 per round) |
| 7 | Agro-infiltrate TRV constructs |
| 9 | GRBV inoculation via *S. festinus* / agroinfection |
| 10 | Confirm VIGS silencing (RT-qPCR) → assign temperature chambers |
| 11–22 | Temperature incubation + weekly phenotyping |
| 14, 18, 22 | Titer sampling (4, 8, 12 wpi) |
| 23–26 | qPCR, immunoblot, data analysis, statistical tests |
| 27–30 | Cross-round reproducibility + writing |

**Total project duration: ~7 months** for one complete round with three temporal replicates overlapped.

## 11. Materials & reagents

| Category | Item | Source | Notes |
|---|---|---|---|
| Plant material | *V. vinifera* cv. Merlot | in-house / FPS UC Davis | Virus-tested stock |
| Vector | pTRV1, pTRV2 | Addgene #8446 / #8447 | Base VIGS system |
| Agrobacterium | GV3101 | Lab stock | — |
| Virus inoculum | GRBV-positive vines or infectious clone (Yepes 2018) | Collaborator | Titer-matched |
| Vector insect | *S. festinus* colony | Almeida lab UC Berkeley or in-house | Viruliferous for 7 d prior |
| qPCR reagents | PowerUp SYBR, nuclease-free water | Thermo | — |
| Antibody | Anti-HSP17.6 polyclonal | Agrisera AS07 255 | For subset immunoblot |
| Chambers | 2 × Conviron/Percival | Institutional core | Verify ± 0.5 °C, 16 h photoperiod |

## 12. Rigor & reproducibility (ARRIVE-adjacent principles for plant studies)

- **Randomisation:** random chamber position; plants coded, scorer blinded to construct identity during phenotyping.
- **Blinding:** qPCR plate layout randomised, sample IDs decoded only after Ct calling.
- **Independent rounds:** ≥ 3 rounds ≥ 2 weeks apart, using distinct seedling batches and Agrobacterium transformations.
- **Data availability:** raw Ct values, images, and analysis scripts deposited to the [GitHub repository](https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis) (companion branch for wet-lab data).
- **Pre-registration:** protocol frozen to the SHA of this commit before any round begins.

## 13. Risks and contingencies

| Risk | Probability | Mitigation |
|---|---|---|
| HSFA2 VIGS inefficient in Merlot | Medium | Screen 2–3 insert positions in parallel pilot; add sense+antisense hairpin (amiRNA) as backup |
| *S. festinus* colony unavailable | Medium | Fall back on GRBV agroinfectious clone (Yepes 2018, ATCC number depending) |
| Temperature chambers not at target ± 1 °C | Low | Pre-test with data loggers for ≥ 72 h before start |
| High biological variance saturates power | Medium | Pre-planned increase to n = 12 per titer cell if pilot variance > 0.5 |
| GRBV titer below qPCR LOD | Low | Include serial dilution standard curve; LOD ≈ 100 copies/μL |

## 14. Next concrete steps (from this protocol to bench start)

1. **Resolve the Merlot HSFA2 locus** from the NCBI GFF on FARM (`resources/genome/genome.gff`) — a small Rscript will suffice.
2. **Design insert + primers** against the resolved locus (siRNA-Finder for off-target screen).
3. **Order oligos** for cloning and qPCR.
4. **Seek IBC approval** for GRBV inoculation experiments if institutional policy requires.
5. **Pre-register** the protocol (OSF or institutional repository) with a frozen SHA of this document.

## 15. Analysis figure template (R / RGraphMaker style)

A two-panel figure drives the causal claim. Panel A shows titer as a function of construct × temperature (the interaction test). Panel B shows HSFA2 and sHSP knockdown validation. Accessible palette, grayscale-safe, error bars as 95 % CIs (per scientific-manuscript-writing guidance).

```r
# HSFA2_VIGS_figure.R — renders Figure 1 for the functional validation paper
library(ggplot2); library(patchwork); library(dplyr)

# Expected input:
#   qpcr  : data.frame(plant_id, construct, temperature, time_wpi, gene, dCt)
#   titer : data.frame(plant_id, construct, temperature, time_wpi, log10_titer)

accessible_palette <- c(VIGS_HSFA2 = "#E69F00", TRV_empty = "#0072B2", Mock = "#999999")

# Panel A: construct x temperature interaction on titer at 8 wpi
pA_df <- titer %>% filter(time_wpi == 8) %>%
  group_by(construct, temperature) %>%
  summarise(mean = mean(log10_titer, na.rm = TRUE),
            ci   = 1.96 * sd(log10_titer, na.rm = TRUE) / sqrt(sum(!is.na(log10_titer))),
            .groups = "drop")

pA <- ggplot(pA_df, aes(temperature, mean, color = construct, group = construct)) +
  geom_line(linewidth = 1.1) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = 0.1) +
  geom_point(size = 3) +
  scale_color_manual(values = accessible_palette, name = "Construct") +
  labs(x = "Temperature regime",
       y = expression(log[10]~"GRBV titer (copies mg"^{-1}~"leaf) at 8 wpi"),
       title = "A. Construct x temperature interaction on GRBV titer") +
  theme_classic(base_size = 11)

# Panel B: HSFA2 + sHSP panel knockdown (delta Ct relative to TRV-empty)
pB_df <- qpcr %>% filter(time_wpi == 4, gene %in% c("HSFA2", paste0("sHSP_", 1:5))) %>%
  group_by(construct, gene) %>%
  summarise(dCt_mean = mean(dCt, na.rm = TRUE),
            ci = 1.96 * sd(dCt, na.rm = TRUE) / sqrt(sum(!is.na(dCt))),
            .groups = "drop")

pB <- ggplot(pB_df, aes(gene, dCt_mean, fill = construct)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = dCt_mean - ci, ymax = dCt_mean + ci),
                position = position_dodge(0.8), width = 0.25) +
  scale_fill_manual(values = accessible_palette, name = "Construct") +
  labs(x = NULL, y = expression(Delta*"Ct (ref. genes)"),
       title = "B. HSFA2 and sHSP knockdown") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig1 <- (pA | pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("Figure1_HSFA2_VIGS_functional.pdf", fig1, width = 10, height = 4.5, device = cairo_pdf)
ggsave("Figure1_HSFA2_VIGS_functional.png", fig1, width = 10, height = 4.5, dpi = 300)
```

Save the script next to the data under `results/functional/` when the wet-lab round runs.

---

**Sign-off block (to be completed at pre-registration):**

| Role | Name | Date | Signature |
|---|---|---|---|
| PI | | | |
| Lab manager | | | |
| Biosafety officer (IBC) | | | |
