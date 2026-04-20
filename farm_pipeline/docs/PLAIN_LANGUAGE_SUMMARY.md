# The GRBV study in plain language

*A walkthrough of what we did, what we found, and why it matters — written for
anyone who isn't a grape-virology specialist.*

## The problem

Grapevine Red Blotch Virus (GRBV) is a relatively new disease of wine grapes.
Infected vines produce fruit with less sugar, odd flavours, and delayed ripening.
The US grape industry loses tens of thousands of dollars per hectare per year to
it, and there is no cure — growers can only pull out infected vines.

We wanted to know: what is the virus actually *doing* to the plant? If we can
understand the mechanism, we can predict which years will be worse, breed more
resistant vines, or time harvests around the infection.

## What we did

We read every gene's activity level in 327 grape leaf samples of cv. Merlot.
These samples covered:

- **Two climate-different growing seasons** (2021 cooler, 2022 much hotter —
  about 7.6 °C warmer at harvest)
- **Four berry ripening stages** (before colour change, at colour change, after,
  and at harvest)
- **Seven years of infection history** (some vines newly infected, some carrying
  the virus for six years)
- **Healthy control vines** matched to each condition

Then we asked three kinds of question of the data:

1. **Standard analysis** (our main paper): which genes change, by how much,
   and when?
2. **Network analysis** (the FARM pipeline): how do genes *organise* into
   coordinated groups, and how does infection re-arrange those groups?
3. **Causal analysis**: among all the changes, which ones actually *cause*
   the disease, versus just bystanders?

## The big picture finding

**The virus does not make the plant's immune system fight back. It disables
the plant's thermal defence system.**

Think of it this way. Plants have two lines of defence against viruses:

- An **immune response** (like our adaptive immunity — specific proteins made to
  recognise the threat)
- A **heat-shock response** (like a thermostat + repair kit — kicks in under any
  stress and includes chaperone proteins called HSPs that physically protect
  cellular machinery)

In GRBV-infected vines, the immune response is barely activated — we checked ten
different canonical immune-defence categories and none of them showed meaningful
up- or down-regulation. But the heat-shock response, specifically the small HSPs
(sHSPs), is **knocked down hard**: 29 of 36 of them are turned off. The master
switch that normally turns on the sHSPs (a transcription factor called HSF) is
also suppressed.

This is important because heat normally *helps* the plant fight the virus. Hot
weather triggers HSF, which turns on sHSPs, which suppress viral replication.
The virus has apparently evolved to block this line of defence. In a hot year,
the plant would normally win; with the virus's sabotage, it can't.

## Specific findings that support this story

- **Seasons matter differently.** In the cool 2021 season, disease signals
  peaked at harvest. In the hot 2022 season, they peaked earlier (at post-
  veraison). The virus's footprint *shifts* with temperature.

- **One gene mediates the whole temperature effect.** Using an objective causal-
  inference method on 3,000 candidate genes, we identified a single gene
  (`LOC132254135`) that carries 100 % of the detected temperature → viral-load
  effect. When it's up (during heat), viral load goes down. This is the single
  most specific, testable prediction in the paper.

- **The plant's gene-network skeleton stays intact.** We took the healthy-vine
  gene co-expression network and tested whether it survives in infected vines.
  All six modules pass with flying colours — the virus perturbs gene *levels*
  within existing groups but does not break the groups apart.

- **Two viral genes are the most rewired hubs.** In the differential co-
  expression analysis, two of the top five most-reconnected genes are the
  virus's own — a direct, independent signature of infection at the network
  level, useful as a quality control.

- **The sugar phenotype has a regulatory signature.** Gene modules that
  contain the sugar-handling genes are enriched for binding sites of DOF
  transcription factors — the canonical plant regulators of carbon
  metabolism. This ties the transcriptomic evidence to the grower-visible
  symptom (low-sugar fruit) through known biology.

- **No "severity subtype" exists in the data.** A common question is whether
  some infected vines have a mild form and others a severe form. We looked
  carefully with an unsupervised clustering method; the two clusters that
  emerged were simply early- versus late-season samples. Maturity is by far
  the biggest source of variation, not disease severity.

## What we can predict (and haven't yet tested)

If the above mechanism is right, then **silencing HSF in a healthy vine should
make it behave like an infected vine** — more susceptible to GRBV, and no
longer protected by hot weather. A companion protocol
([HSFA2_VIGS_EXPERIMENT_PLAN.md](HSFA2_VIGS_EXPERIMENT_PLAN.md)) lays out
exactly how to test this.

## Why it matters

Three practical consequences:

1. **For breeders.** The vulnerable node is HSF (the heat-stress master
   regulator), not a canonical immune gene. Breeding programmes should screen
   for varieties with more robust HSF expression under infection.
2. **For growers.** Disease severity will track seasonal temperature in a
   predictable way. Hot spells *after* post-veraison should reduce viral titer;
   cool late seasons will make the disease worse.
3. **For the field.** GRBV is the first plant virus where we have direct
   evidence (observational, to be confirmed functionally) that the infection
   mechanism is suppression of thermal defence rather than of classical
   immunity. That re-frames how we think about plant–virus interactions more
   generally.

## What's in the repository

- **Main paper pipeline** — 55 figures, 26 data tables, 20 analysis scripts
  covering DE, pathway enrichment, biomarkers, ML classification, GSEA,
  mediation, multifactor ANOVA.
- **Extended FARM HPC pipeline** — 13 rules producing the network-level,
  causal, and cis-regulatory analyses summarised above.
- **Functional protocol** — a wet-lab plan to confirm the mechanism using
  VIGS (virus-induced gene silencing) of the master heat-shock regulator.
- **Comparative pipeline** — harmonises your DEG signature against published
  GLRaV-3, GPGV, powdery-mildew and heat-stress signatures to ask how much of
  this is GRBV-specific versus generic grape stress.

## Numbers you can drop into a grant or talk

| Metric | Number |
|---|---|
| Samples | 327 (264 RB+, 63 RB−) |
| Genes tested | 14,604 |
| Seasons compared | 2 (2021, 2022) |
| Maturity stages | 4 |
| Universal DEGs (robust across all analyses) | 908 |
| sHSPs affected | 29 down / 36 expressed |
| HSF–sHSP correlation | r = 0.625, p = 7.2 × 10⁻³⁷ |
| Genes that mediate the temperature → titer link | 1 (`LOC132254135`) |
| Defense categories with |log₂FC| > 0.5 | 0 of 10 |
| Co-expression modules preserved in infected vines | 6 of 6 |
| Directed TF → target regulatory edges | 100,000 |
| Pipeline rules completed on FARM HPC | 13 of 13 |
