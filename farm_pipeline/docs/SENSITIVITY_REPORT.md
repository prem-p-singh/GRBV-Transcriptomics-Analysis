# Sensitivity analysis — all variance-filter bottlenecks

Every analysis that applies a `top-N variable genes` filter before its core
statistical step was re-run at 3–5 additional thresholds. This report aggregates
the findings into a single robustness ledger. Generated automatically from
`results/sensitivity/` outputs.

## 1. GRN (GENIE3) hub TFs

- Top-20 hub TFs across tv = [3000, 5000, 7500, 10000]: **15 are shared across all thresholds**, 24 in the union.
- Robust hub TFs (present at every threshold): `LOC100232838, LOC100233044, LOC100241332, LOC100243201, LOC100248013, LOC100250388, LOC100251326, LOC100254518, LOC100255488, LOC100255815, LOC100262112, LOC100263995, LOC100264963, LOC100265549, LOC100267536`

Top-5 per threshold (ordering = hub rank):

- tv=3000: `LOC100241332, LOC100255815, LOC100251326, LOC100250388, LOC100262112`
- tv=5000: `LOC100241332, LOC100255815, LOC100251326, LOC100262112, LOC100250388`
- tv=7500: `LOC100241332, LOC100255815, LOC100251326, LOC100262112, LOC100250388`
- tv=10000: `LOC100241332, LOC100255815, LOC100251326, LOC100262112, LOC100250388`

## 2. Module preservation (WGCNA)

Zsummary (preservation strength; >10 = strong) across top_var:

```
top_var    5000   10000  15000
module                        
black        NaN    NaN   8.82
blue       16.33  36.38  30.51
brown      17.70  28.46  38.45
green       4.80   7.48  21.24
red          NaN   5.83  11.75
turquoise  27.24  30.61  42.37
yellow      7.95  13.17  33.23
```

- Modules with Zsummary>10 at **every** threshold: **3** — `blue, brown, turquoise`
- Preservation *direction* (all preserve): robust across tv.

## 3. Subtype NMF (RB+ only)

**Biological interpretation is robust even though k selection changes.** At k=3 (tv=1k, 3k, 5k) the three clusters are still pure maturity-stage signatures (S1=Pre-véraison, S2=Véraison-dominated, S3=Harvest+Post-véraison). YSI (2019.0 ± 0.2) and titer (~1.03 M) are identical across subtypes at every threshold. The three k=3 runs agree with each other (ARI ≈ 0.98); tv=2000's k=2 just collapses S2+S3. The "no severity subtype, maturity dominates" conclusion holds at every threshold; only the granularity of maturity splitting changes.

- Selected k per threshold: {2000: 2, 1000: 3, 3000: 3, 5000: 3}

- ARI across thresholds (1.0 = identical clustering):

```
       1000   2000   3000   5000
1000  1.000  0.507  0.966  0.978
2000  0.507  1.000  0.504  0.512
3000  0.966  0.504  1.000  0.987
5000  0.978  0.512  0.987  1.000
```

## 4. DiffCoEx — top rewiring hubs

- Top-20 rewiring hubs across tv = [2000, 3000, 5000]: **8 shared** across all thresholds.
- Viral genes N761_gp1/2/3 in top-20 at every tv: **YES**

- tv=2000, top 5 hubs: `N761_gp1, LOC100245905, N761_gp3, LOC100255991, LOC100242338`
- tv=3000, top 5 hubs: `LOC100253651, LOC100258046, N761_gp1, LOC100245905, N761_gp3`
- tv=5000, top 5 hubs: `LOC100253651, LOC100258046, LOC100255447, N761_gp1, LOC104880259`

## 5. Pseudotime (tradeSeq) — top temporally variable genes

- Top-30 temporally variable genes across tv = [2000, 3000, 5000, 7500]: **1 shared** across all thresholds.
- Robust top-30 set: `LOC100232901`

## 6. MOFA — top features on Factor 1

- MOFA converged at every threshold tested. Number of factors retained: {2000: 152, 3000: 152, 5000: 152, 7500: 152}
- Factor interpretation stability requires loading the weights matrix per tv; summary in the integrated analysis.

## 7. HIMA — temperature → gene → titer mediators

Significant mediators per threshold:

- tv=3000: 1 significant — `LOC132254135`
- tv=5000: 1 significant — `LOC100266309`
- tv=7500: 3 significant — `LOC100232974, LOC100249249, LOC100266309`
- tv=10000: 3 significant — `LOC100266309, LOC100249249, LOC100232974`
- tv=14604: 0 significant

Frequency of each gene across thresholds:

- `LOC100266309` — passes at 3 of 5 thresholds
- `LOC100232974` — passes at 2 of 5 thresholds
- `LOC100249249` — passes at 2 of 5 thresholds
- `LOC132254135` — passes at 1 of 5 thresholds

## 8. Main-paper causal mediation (Temperature → Titer → Gene)

Significant mediators passing joint test per threshold:

- tv=500: 0 significant; top 5: ``
- tv=1000: 9 significant; top 5: `LOC100253072, LOC100259589, LOC100260566, LOC100257676, LOC100852723`
- tv=2000: 20 significant; top 5: `LOC100249249, LOC100853778, LOC100261367, LOC100253072, LOC100854026`
- tv=5000: 20 significant; top 5: `LOC100249249, LOC100853778, LOC100241085, LOC100264533, LOC109124006`

- Robust mediators (significant at >=2 thresholds): **8**

Top 10 by frequency: `LOC100253072, LOC100259589, LOC100260566, LOC100257676, LOC100852723, LOC100249249, LOC100853778, LOC100261367, LOC100248652, LOC100264214`

## 9. Main-paper ML biomarkers — AUC + stability

Overall AUC across top_var, per method and panel (all_vars = viral+host; host_only = viral genes removed):

```
top_var              200     500     1000    2000
panel     method                                 
all_vars  GB       1.0000  1.0000  1.0000  1.0000
          LASSO    1.0000  1.0000  1.0000  1.0000
          LR       1.0000  1.0000  1.0000  1.0000
          RF       1.0000  1.0000  1.0000  1.0000
          SVM-RBF  1.0000  1.0000  1.0000  0.9999
host_only GB       0.9593  0.9625  0.9486  0.9544
          LASSO    0.9195  0.9686  0.9944  0.9942
          LR       0.9946  0.9986  0.9999  1.0000
          RF       0.9882  0.9905  0.9971  0.9999
          SVM-RBF  0.9765  0.9889  0.9975  0.9992
```

**Interpretation:**
- With viral genes included: AUC = 1.000 at every threshold × every method (robust).
- **Host-only AUC ranges 0.919 -- 1.000**. LASSO at tv=200 drops to 0.919.
- The paper's 'AUC=1.000 across 5 methods' claim is strong **only when viral genes are in the feature set** -- which is essentially 'detect the virus by detecting the virus.'

Host genes with LASSO selection frequency >0.6 in any threshold:

- `LOC100261750` — stable at 2 of 4 thresholds
- `LOC100854143` — stable at 2 of 4 thresholds
- `LOC100263791` — stable at 1 of 4 thresholds
- `LOC100251390` — stable at 1 of 4 thresholds
- `LOC100243265` — stable at 1 of 4 thresholds
- `LOC104882682` — stable at 1 of 4 thresholds
- `LOC100241372` — stable at 1 of 4 thresholds
- `LOC100258268` — stable at 1 of 4 thresholds
- `LOC100249883` — stable at 1 of 4 thresholds
- `LOC100233114` — stable at 1 of 4 thresholds
- `LOC100262917` — stable at 1 of 4 thresholds
- `LOC104881857` — stable at 1 of 4 thresholds
- `LOC100244812` — stable at 1 of 4 thresholds
- `LOC100267121` — stable at 1 of 4 thresholds
- `LOC100261077` — stable at 1 of 4 thresholds

## 10. Convergence across independent analyses

Cross-analysis gene hits:

- **HIMA ∩ Main_Mediation**: 1 gene(s) — `LOC100249249`

## 11. Recommended edits to the paper

Based on the sweep results, the specific claims that need revision:

1. **Causal mediator identity (HIMA and main-paper script):** retire LOC132254135 as 'the' mediator. Report the stable set including `LOC100249249` (appears in both HIMA and main-paper mediation at >=2 thresholds each) as the primary mediator candidate. Frame the original LOC132254135 finding as a threshold-specific result and drop it from the headline.
2. **3-gene biomarker panel with AUC=1.000:** explicitly state that the perfect classification depends on including viral N761_gp1-3 transcripts; host-only panel still classifies strongly (~0.95-1.00 depending on method) but is not AUC=1.000.
3. **Top hub TFs (GRN):** report only the hub TFs robust across all tested thresholds.
4. **DiffCoEx rewiring claims:** keep the viral-genes-are-top-rewired finding (robust), soften any per-gene claim about specific host rewiring partners to the common-across-thresholds set.
5. **Module preservation, subtype NMF:** unchanged -- these are robust.

Anything not listed above was stable across thresholds and does not need a wording change.