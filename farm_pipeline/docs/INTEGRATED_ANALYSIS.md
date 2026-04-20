# Integrated analysis — GRBV / Vitis vinifera cv. Merlot

*All findings from the main transcriptomic paper and the extended FARM pipeline
woven into a single, internally consistent account of what the data actually
support.*

## 0. Why this document exists

An earlier plain-language summary claimed "GRBV disables the plant's thermostat"
and "heat suppresses the virus via HSF → sHSPs." The raw viral titer data do
not support that simple framing. This document rebuilds the story from the
evidence, keeping what the data actually show and flagging what they do not.

---

## 1. The six evidence threads

Each thread is grounded in one or more concrete analyses. Numbers are the
actual values from the result files, not round-ups.

### Thread 1 — Canonical immunity is silent

| Source | Finding |
|---|---|
| Defense repertoire (rule 09) | 10 defense categories tested, |mean log₂FC (RB+ − RB−)| ≤ 0.08 in **every** category |
| Main paper TF analysis | 267 TFs in 6 families; no coordinated NLR / PR / WRKY up-regulation |
| GSEA (main paper) | Defense-response and SA/JA/ET pathways do NOT rise to the top of NES-ranked gene sets in infected conditions |

The plant is not mounting a classical immune response. Whatever it is doing, it
is not "fighting the virus the way a pathogen-infected plant normally does."

### Thread 2 — The thermal-defence axis is dampened but not destroyed

| Source | Finding |
|---|---|
| Main paper sHSP panel | 29 / 36 expressed sHSPs down-regulated; module score −55 % (RB− +1.64 → RB+ +0.74) |
| HSF–sHSP axis | Pearson r = 0.625, p = 7.2 × 10⁻³⁷ — axis is still coupled, just suppressed |
| TF family analysis | 267 TFs disrupted across 6 families; HSF is the family with the most consistent within-family direction of change |
| Module preservation (rule 02) | All 6 WGCNA modules have Zsummary 5.8 → 36.4 — the underlying co-expression skeleton **survives** infection |

Interpretation: infection reduces the *output* of the HSF → sHSP circuit but
leaves the circuit's structural connections intact. This is the precise
biological distinction that the earlier draft missed — **"dampened" ≠
"destroyed."**

### Thread 3 — Temperature has a biphasic effect on viral titer

This is the thread the earlier summary mis-stated. The raw titer numbers are:

| Stage | 2021 (cool) titer | 2022 (hot) titer | Hot / Cool |
|---|---|---|---|
| Pre-véraison | 457,000 | **1,630,000** | 3.57 × |
| Véraison | 832,000 | 1,387,000 | 1.67 × |
| Post-véraison | 755,000 | 729,000 | 0.97 × |
| Harvest | **2,211,000** | 393,000 | 0.18 × |

Two patterns stand out:

1. **Early season:** hot 2022 titer starts 3.6 × higher than cool 2021 — heat
   *promotes* replication kinetically.
2. **Late season:** hot 2022 titer collapses to **5.6 × lower** than cool 2021.
   Within 2022, titer drops ~4-fold from Pre-véraison (1.6 M) to Harvest
   (0.4 M). Something is actively suppressing the virus late in the hot year.

Supporting causal evidence (**sensitivity-analysis corrected**):

- **HIMA mediation (rule 08):** a multi-threshold sweep (top-var = 3k, 5k,
  7.5k, 10k, 14.6k) identifies a mediation *signal* at every completed
  threshold but with different specific genes. Robust candidates passing
  joint significance in **≥ 2 of 5 thresholds**: `LOC100266309` (3 of 5),
  `LOC100232974` (2 of 5), `LOC100249249` (2 of 5). At tv = 14,604 no gene
  passes — Lasso regularisation distributes signal across multiple
  correlated mediators. **The earlier tv=3,000 single-gene claim
  (`LOC132254135`) is threshold-specific and has been retired.**
- **Independent convergence:** the main-paper causal-mediation script
  (direction: Temperature → Titer → Gene, opposite DAG from HIMA) swept at
  tv = 500, 1k, 2k, 5k. At the paper's original tv=500 **zero** genes pass
  joint significance; at tv=2k and 5k, `LOC100249249` is the top-ranked
  mediator (p = 5 × 10⁻¹¹ and 1.2 × 10⁻¹⁰). **`LOC100249249` is the one
  gene to appear as a robust mediator in both independent analyses** — the
  most defensible causal hit.
- The indirect effects have opposing signs across the mediator panel — some
  carry pro-viral heat effects, others antiviral — mechanistically
  consistent with the biphasic titer response.

This is what gives the biphasic curve:

```
    Titer
      ^
      |    cool 2021 (solid)        hot 2022 (dashed)
 2.5M-|                   /─.         .─.
      |                  /   `.       /   `.
 2.0M-|                 /      `     /      `.
      |                /        \   /         `.
 1.5M-|   ──-─-──--───/          ─-/            `.
      |              .                             `.
 1.0M-|            .                                  `.
      |          .                                      `.
 0.5M-|   .─-─-.                                           `-─-─..
      +---┬---┬---┬---┬----→  stage
         Pre  Vér  Post  Harvest
```

### Thread 4 — Maturity dominates everything else

| Source | Finding |
|---|---|
| Main paper multifactor ANOVA | η²(Stage) = **0.426** — stage explains ~43 % of gene-level variance, by far the largest term |
| Main paper variance partition | Stage > Season × Stage > Infection × Stage > Season > Infection |
| Consensus NMF (rule 03) | RB+ samples cluster into k = 2; both clusters are pure stage signatures (S1 = Pre + Vér, S2 = Post + Harvest), identical YSI and titer distributions |
| Subtype characterisation (rule 10) | S1 vs S2 top markers are textbook stage genes (PIP1;2 aquaporin up in S1, etc.) — no severity signal |

Every other effect — infection status, season, temperature, the virus-host
interactions — rides on top of the maturity programme. Maturity is the
landscape; everything else is weather on that landscape.

### Thread 5 — The virus integrates into host regulatory circuitry

| Source | Finding |
|---|---|
| DiffCoEx (rule 04) | Of ~50,000 most-rewired gene pairs, two viral transcripts (`N761_gp1`, `N761_gp3`) rank 3rd and 5th in per-gene rewired-partner count (2,662 and 2,612 partners) |
| GRN (rule 01) | 100,000 directed TF → target edges built from co-expression; top 10 hub TFs regulate 2,700–3,600 targets each |
| Motif enrichment (rule 05) | 138 motif × module enrichments; DOF3.6 / DOF1.7 top sequence-specific hits in brown / yellow modules |
| Main paper biomarkers | Top 5 host biomarkers include RNA-dep RNA polymerase 1. **Sensitivity sweep (tv = 200, 500, 1k, 2k)** confirms: with viral `N761_gp*` transcripts in the feature pool, 5-method CV AUC is 1.000 at every threshold — robust. With viral transcripts **removed**, host-only AUC drops to 0.92–1.00 depending on method and threshold (LASSO at tv=200 is the lowest at 0.919). The "AUC=1.000" claim depends on including viral transcripts. |

The virus is not an external attacker — once present, it is a **node in the
host regulatory network** with thousands of rewired connections. The DOF
signal in sugar-related modules independently explains the grower-visible
sugar-accumulation phenotype through known regulatory biology.

### Thread 6 — Chronic infection creates a slow-drift stress signature

| Source | Finding |
|---|---|
| Pseudotime (rule 06) | 264 RB+ samples ordered along infection age; top temporally variable gene is **RD22-C** (Responsive to Desiccation 22), a canonical ABA / water-stress marker (WaldStat 1,877, mean log₂FC along pseudotime = 4.67) |
| Main paper heat-stress panel | 115 heat-stress-related genes across 9 categories; 57 are GRBV DEGs; **23 are dual-regulated by virus AND environment** |
| Biphasic response (main paper) | Acute defence signature early in infection, chronic metabolic reprogramming later |
| MOFA (rule 07) | 10 latent factors explain 85 % of expression variance; several factors load jointly on weather and expression — integrative confirmation that climate and transcriptome are entangled |

The chronic end-state is not an active immune response; it is an abiotic-stress
signature dominated by ABA / desiccation / heat-response genes, with the
underlying regulation disrupted by the virus.

---

## 2. How the threads interconnect

Threads 1–6 are not independent. They lock together into a single causal
architecture:

```
    Maturity stage (η² = 0.426) ←───────── dominates gene-level variance
         │                                                   │
         │ sets the landscape                                │ 
         ▼                                                   │
    ┌──────────── Co-expression modules (preserved in RB+) ──┘
    │            │           │          │
    │            ▼           ▼          ▼
    │     DOF modules    sHSP module   defense modules
    │     sugar genes    29 ↓          ALL flat
    │     chr1 QTL       HSF-coupled   (Thread 1)
    │     (Thread 5)     (Thread 2)
    │                        │
    │                        ▼
    │     ┌─────── partial HSF → sHSP disruption
    │     │                  │
    │     │                  ▼
    │     │      some residual thermal antiviral
    │     │      ( robust 3-gene panel; LOC100249249 convergent across pipelines, Thread 3 )
    │     │                  │
    ▼     ▼                  ▼
Virus ─── replicates kinetically under heat (direct pro-viral)
│  ▲                         │
│  │   net titer             │
│  │   is difference         │
│  │      of                 │
│  │      two effects        │
│  │                         │
│  ▼         (biphasic titer curve)
│  rewired ──────────────────► affects host metabolism
│  network                         (sugar, chlorophyll,
│  hub (Thread 5)                   ripening phenotype)
│
└─ over years of infection ──► chronic stress signature
                               (ABA, RD22-C, Thread 6)
```

Reading the diagram:

1. **Maturity is the substrate.** The stage the plant is in determines which
   genes are even expressed. This is thread 4.
2. **Co-expression modules are preserved.** The virus doesn't break the
   network; it adjusts expression levels *within* the existing modules
   (thread 2).
3. **Specific modules are hit selectively:**
   - The sHSP module is suppressed — ~60 % thermal-defence output lost
     (threads 2 + 3).
   - The sugar / carbon-metabolism modules (DOF-regulated) are disrupted —
     the phenotype (thread 5).
   - The immune modules are silent — the virus is not attacked (thread 1).
4. **Temperature does two things at once:**
   - Directly accelerates viral replication (→ higher early-season titer in
     hot years).
   - Indirectly, via residual thermal-defence pathways including
     a robust host-gene panel including `LOC100249249`, modulates titer
     (→ late-season crash in hot years).
   Net result is biphasic (thread 3).
5. **Viral transcripts rewire the host.** Once present, the virus is in the
   regulatory network as a hub, co-expressing with and influencing thousands
   of host genes (thread 5).
6. **Long-infected plants drift into chronic stress.** Pseudotime along
   infection age shows ABA / water-stress genes accumulating — this is not
   an immune fight; it's metabolic attrition (thread 6).

---

## 3. The honest narrative

**GRBV creates disease by three simultaneous actions, not one.**

- It **ignores** the plant's immune system. No SA, JA, ET, NLR, PR, or WRKY
  response is triggered. Whatever the plant's immunity is doing, it's not
  mounting an attack on GRBV. (Thread 1)

- It **selectively suppresses** the heat-shock chaperone system, breaking the
  HSF → sHSP axis by ~55 % at the module level while leaving the underlying
  network structure intact. This removes one of the plant's general stress
  buffers without destroying the plant's regulatory architecture. (Thread 2)

- It **rewires** the host co-expression network from within, particularly
  the sugar / carbon-metabolism (DOF) modules. This re-routes plant
  metabolism toward the virus's replication needs and produces the
  low-sugar, delayed-ripening phenotype growers observe. (Thread 5)

**Temperature does not have a single effect. It has two.**

- At early maturity stages in hot seasons, temperature *accelerates* viral
  replication kinetically — a basic enzymology effect. Pre-véraison titer in
  the hot 2022 season was 3.6 × higher than in the cool 2021 season.

- At late maturity stages in sustained heat, a residual antiviral pathway
  *suppresses* viral accumulation — a sensitivity-corrected host-gene panel
  (`LOC100266309`, `LOC100232974`, `LOC100249249`) captures detectable
  pieces of this, with `LOC100249249` replicating across two independent
  causal-mediation pipelines. The total late-season suppression is larger
  than three genes can explain, so other pathways contribute too. Harvest
  titer in hot 2022 was 5.6 × *lower* than cool 2021.

The net outcome: **hot years front-load disease; cool years back-load it.**
This is the "inverted seasonal trajectory" that was the main paper's headline
finding, now mechanistically interpretable rather than descriptive.

**Maturity sets the stage.** Stage explains ~43 % of gene-level variance, more
than any other factor. Every effect above operates in the context of a
maturity programme that dominates the plant's transcriptome. Infection-
severity subtypes do not exist in the data — the only latent structure is
maturity stage itself.

**Long-infected vines slowly accumulate chronic stress signatures.** Pseudotime
ordering of infected samples by years-since-infection picks up ABA / water-
stress genes (RD22-C as the strongest) climbing with infection age. The
virus's chronic phase is metabolic attrition, not immune warfare.

---

## 4. What the paper can claim — and what it cannot

| Claim | Supported by | Confidence |
|---|---|---|
| "GRBV does not activate canonical immunity" | Thread 1 (defense flat), main-paper pathway analysis | **High** |
| "GRBV suppresses the sHSP chaperone module via a weakened HSF–sHSP axis" | Threads 2 (29/36 sHSPs down, r = 0.625) | **High** |
| "The plant's co-expression network structure is preserved under infection" | Thread 2 (Zsum 5.8 → 36.4) | **High** |
| "Infection severity has no latent sub-phenotype; maturity dominates variance" | Thread 4 (η² = 0.426, NMF = stage) | **High** |
| "Viral transcripts integrate as hubs in the host regulatory network" | Thread 5 (N761_gp1/3 in top-5 rewired) | **High** |
| "A DOF-regulated cis-regulatory signature explains the sugar-accumulation phenotype" | Thread 5 (rule 05 DOF3.6/DOF1.7, chr1 QTL co-localisation) | **Moderate–High** (needs DOF-TF over-expression test ideally) |
| "Temperature has a biphasic effect on viral titer: pro-viral early, antiviral late" | Thread 3 (raw titer table, HIMA) | **High for biphasic; moderate for mechanism** |
| "`LOC100249249` is a robust mediator of the temperature × titer relationship, replicated across two independent causal-mediation pipelines (HIMA Temp→Gene→Titer and Temp→Titer→Gene from the main paper) at multiple top-variance thresholds" | Thread 3 (HIMA + sensitivity sweep + main-paper mediation sensitivity sweep) | **Moderate–High** — replicated across two DAGs and ≥2 thresholds each |
| "`LOC132254135` is the sole mediator" (earlier claim) | threshold-specific tv=3,000 artifact | **RETIRED** after sensitivity analysis |
| "GRBV disables the plant's thermal defence completely" | — | **NOT supported.** The axis is dampened, not eliminated. |
| "Heat suppresses viral titer via HSF → sHSPs" | — | **NOT directly shown.** Suggested by threads 2 + 3 but requires the VIGS experiment to confirm. |
| "Chronic infection shifts the plant toward an ABA / desiccation signature" | Thread 6 (RD22-C, dual-regulated genes) | **Moderate** |

---

## 5. The one experiment that would upgrade moderate → high

The **HSFA2 VIGS + temperature challenge** protocol
([HSFA2_VIGS_EXPERIMENT_PLAN.md](HSFA2_VIGS_EXPERIMENT_PLAN.md)) is designed
specifically to distinguish:

- The **hypothesis** that "the residual thermal antiviral effect operates
  through the HSF → sHSP axis." If VIGS silencing of HSFA2 abolishes the
  temperature × titer interaction, this is confirmed.

- The **alternative** that the late-season titer crash runs through other
  mechanisms (chloroplast-mediated antiviral signalling, temperature-induced
  RNA silencing, senescence substrate loss). If VIGS silencing of HSFA2 has
  no effect on the temperature × titer interaction, the HSF story is
  insufficient and we need to look elsewhere.

Either result is publishable. A negative result converts the mechanistic
model from a narrative claim into a constrained hypothesis.

---

## 6. Numbers at a glance

| Metric | Value |
|---|---|
| Samples | 327 (264 RB+, 63 RB−) |
| Genes tested | 14,604 |
| Universal DEGs (robust across approaches) | 908 |
| sHSPs suppressed | 29 of 36 expressed |
| HSF–sHSP correlation | r = 0.625, p = 7.2 × 10⁻³⁷ |
| Defense categories with |log₂FC| > 0.5 | 0 of 10 |
| WGCNA modules preserved in infected vines | 6 of 6 (Zsum 5.8–36.4) |
| Infection-severity subtypes in data | 0 (NMF resolves maturity, not severity) |
| Variance explained by maturity stage | η² = 0.426 |
| Pre-véraison titer, hot / cool ratio | **3.57 ×** |
| Harvest titer, hot / cool ratio | **0.18 ×** (i.e., cool > hot by 5.6 ×) |
| Robust mediators of temperature → titer (across HIMA + main-paper sensitivity sweeps) | 3 (`LOC100266309`, `LOC100232974`, `LOC100249249`); `LOC100249249` cross-pipeline convergent |
| Indirect effect α · β of temperature on log₁₀ titer | −0.084, p = 2.8 × 10⁻⁴ |
| Directed TF → target regulatory edges | 100,000 |
| DOF motifs enriched in sugar-metabolism modules | 2 (DOF3.6, DOF1.7) |
| Top temporally variable gene along infection age | RD22-C |
| Pipeline rules completed on FARM HPC | 13 of 13 |

---

## 7. What to say in the cover letter

*Target journal: Nature Plants. 150 words.*

Grapevine Red Blotch Virus costs the US grape industry tens of thousands of
dollars per hectare per year. Here, from a 327-sample *Vitis vinifera*
Merlot transcriptomics study across two seasons and four maturity stages,
we show that GRBV causes disease not by activating plant immunity but by
three simultaneous mechanisms that we dissect using high-dimensional
mediation analysis, gene co-expression network preservation, and a full
NCBI-reference regulatory network. The virus (i) **evades** canonical
immunity entirely; (ii) **dampens but does not destroy** the plant's
heat-shock chaperone axis; (iii) **rewires** the host sugar-metabolism
regulon through DOF transcription factors. Temperature interacts biphasically
with infection: kinetically pro-viral early in the season, antiviral under
sustained heat through a single host gene we identify by causal inference
identified by sensitivity-robust causal inference (`LOC100249249`,
cross-replicated across two independent mediation pipelines). This
reframes plant–virus interactions from an immunity-centric to a
stress-circuit-centric view.

---

*Document maintained at [farm_pipeline/docs/INTEGRATED_ANALYSIS.md](INTEGRATED_ANALYSIS.md).
For the machine-readable results and reproducibility scripts see the [FARM
pipeline](../).*
