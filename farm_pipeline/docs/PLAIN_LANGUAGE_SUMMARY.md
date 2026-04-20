# The GRBV study in plain language

*A walkthrough of what we did, what we found, and why it matters — written for
anyone who isn't a grape-virology specialist.*

**Corrected version** (the earlier draft said "GRBV disables the thermostat";
the viral-titer data show that is too simple). For the rigorous, numeric
account see the companion [INTEGRATED_ANALYSIS.md](INTEGRATED_ANALYSIS.md).

## The problem

Grapevine Red Blotch Virus (GRBV) is a relatively new disease of wine grapes.
Infected vines produce fruit with less sugar, odd flavours, and delayed
ripening. The US grape industry loses tens of thousands of dollars per hectare
per year to it, and there is no cure — growers can only pull out infected
vines.

We wanted to know: what is the virus actually *doing* to the plant? If we can
understand the mechanism we can predict which years will be worse, breed more
resistant vines, or time harvests around the infection.

## What we did

We read every gene's activity level in 327 leaf samples of *Vitis vinifera* cv.
Merlot covering two climate-different seasons (2021 cool, 2022 about 7.6 °C
hotter at harvest), four ripening stages, and vines carrying the virus for
anywhere from one to six years, plus healthy controls. Then we asked three
kinds of question of the data: what changes (differential expression), how are
genes organised (network analysis), and which changes actually *cause* the
disease (causal inference).

## The big picture — in three mechanisms, not one

The virus causes disease by doing three things at the same time:

**1. Ignoring the plant's immune system.** We checked ten canonical defence
categories (pathogen receptors, resistance proteins, salicylic-acid signalling,
jasmonate signalling, ethylene signalling, WRKY TFs, PR proteins, callose
deposition, reactive-oxygen species, and R-protein adapters). None of them
showed meaningful up- or down-regulation in infected vines. The plant is
**not** mounting an immune response against GRBV.

**2. Selectively damping the plant's heat-shock chaperone system.** Small
heat-shock proteins (sHSPs) are a family of 36 expressed chaperones that
protect cellular machinery against stress. GRBV turns **29 of them off**. The
master regulator that normally switches sHSPs on, a transcription factor
called HSF, is also suppressed. Importantly, the **structure** of the
co-expression network is preserved — we tested all six gene-network modules
for preservation between healthy and infected tissue and every one of them
survived with very strong statistical preservation. So GRBV is not destroying
the wiring; it is tuning the output down.

**3. Rewiring the plant's sugar and metabolism regulators.** The host genes
that handle sugar accumulation and carbon metabolism end up in different
co-expression groups in infected tissue. The transcription factors whose
binding sites are enriched in these rewired groups are DOFs — the canonical
plant regulators of sugar transport. This is the direct regulatory link to the
grower-visible symptoms (low-sugar fruit, altered ripening). In our network
analyses, the virus's own genes sit among the most-rewired nodes — once
present, the virus is *in* the regulatory network, not outside it attacking.

## What the viral-titer data actually show

A reasonable question: *if the virus weakens heat defences, does heat kill the
virus?*

| Stage | 2021 cool season titer | 2022 hot season titer |
|---|---|---|
| Pre-véraison | 457,000 | **1,630,000** (3.6 × higher) |
| Véraison | 832,000 | 1,387,000 |
| Post-véraison | 755,000 | 729,000 |
| Harvest | **2,211,000** | 393,000 (5.6 × *lower*) |

Temperature doesn't act in one direction. It has **two opposing effects**:

- Early in the season, heat speeds up viral replication kinetically — hot
  years start with much more virus.
- Late in the season, under sustained heat, a residual antiviral pathway
  kicks in — hot years end with much less virus.

Net result: cool seasons see virus build up slowly to a harvest peak; hot
seasons see virus front-load and then crash. This is the "inverted seasonal
trajectory" that was the main paper's headline finding.

When we ran the same causal-inference method at five different gene-pool
sizes (from 3,000 to 14,604 genes) to check robustness, we found that a
**three-gene panel** (`LOC100266309`, `LOC100232974`, `LOC100249249`)
consistently comes up as the statistically detectable piece of the late-
season antiviral effect. **`LOC100249249` is the strongest candidate**: it
shows up as a significant mediator in *both* the gene-is-mediator analysis
(HIMA) and the opposite-direction titer-is-mediator analysis from the main
paper, at multiple gene-pool sizes. Earlier we had named a single gene
(`LOC132254135`) as "the" mediator — we've retired that specific claim
because it only appeared at one particular gene-pool size. The robust
finding is the three-gene panel, with `LOC100249249` as the best-supported
target for wet-lab follow-up.

## What dominates the transcriptome

**Maturity stage**, by a wide margin. In a three-way ANOVA, the main effect of
stage explains 43 % of gene-level variance — more than infection, season, or
any combination. When we clustered infected vines to look for hidden severity
subtypes, the two clusters that emerged were simply early vs late ripening,
with identical infection-year distributions and identical viral titer. There
is no latent "mild versus severe" infection subgroup in the data; maturity is
the landscape, and infection / season / temperature are weather on it.

## What happens in vines that have been infected for years

Ordering infected samples by years-since-infection reveals a **chronic drift
toward abiotic-stress signatures**, not an immune response. The top temporally
variable gene across infection age is RD22-C (Responsive to Desiccation 22),
an ABA / water-stress marker. Long-infected vines look progressively more
like they are dehydrated and heat-stressed. This is consistent with the
thermal-defence damping accumulating over years.

## Why it matters

Three practical consequences:

1. **For breeders.** The vulnerable node is the HSF master regulator, not a
   canonical immune gene. Screen cultivars for robust HSF expression under
   infection rather than for classical disease-resistance proteins.
2. **For growers.** Disease severity tracks seasonal temperature in a
   *biphasic*, predictable way: hot springs front-load the disease; hot
   late-season spells suppress viral load. Cool, even seasons are
   worst-case — they let the virus climb steadily to a harvest peak.
3. **For plant virology.** We have evidence that GRBV works by perturbing
   stress-circuit regulation rather than engaging with classical immunity.
   If this extends to other viruses, it's a conceptual reframe of how we
   think about plant–virus interactions.

## What we can predict (and haven't tested yet)

If the mechanism above is right, then silencing HSF in a healthy vine should
mimic part of what infection does — specifically, it should abolish the
late-season thermal antiviral effect. A companion protocol at
[HSFA2_VIGS_EXPERIMENT_PLAN.md](HSFA2_VIGS_EXPERIMENT_PLAN.md) lays out
exactly how to run this test. Either outcome is publishable: a positive
result confirms the HSF → sHSP mechanism; a negative result tells us the
late-season suppression runs through some other pathway we should pursue.

## Numbers for a grant or talk

| Metric | Number |
|---|---|
| Samples | 327 (264 RB+, 63 RB−) |
| Genes tested | 14,604 |
| Seasons compared | 2 (2021 cool, 2022 hot, +7.6 °C at harvest) |
| Maturity stages | 4 |
| Universal DEGs (robust across all analyses) | 908 |
| sHSPs suppressed | 29 of 36 expressed |
| HSF–sHSP correlation | r = 0.625, p = 7.2 × 10⁻³⁷ |
| Defense categories with |log₂FC| > 0.5 | 0 of 10 |
| Co-expression modules preserved | 6 of 6 |
| Pre-véraison titer, hot / cool ratio | 3.6 × |
| Harvest titer, hot / cool ratio | 0.18 × (i.e., 5.6 × lower in the hot year) |
| Genes mediating temperature → titer (robust across HIMA sweep + main-paper mediation sweep) | 3 candidates (`LOC100249249` convergent across both analyses) |
| Variance explained by maturity stage | η² = 0.43 |
| Directed TF → target regulatory edges built | 100,000 |
| DOF motifs enriched in sugar-metabolism modules | DOF3.6, DOF1.7 |
| Pipeline rules completed on FARM HPC | 13 of 13 |

For the full technical analysis with connections between every piece of
evidence, see [INTEGRATED_ANALYSIS.md](INTEGRATED_ANALYSIS.md).
