# Critical review of the GRBV analysis pipeline

*Simulated journal-reviewer critique. Every issue below is something a peer
reviewer at Nature Plants / Plant Cell / New Phytologist / PNAS is likely
to raise. Issues are graded by severity; remediation steps specify what
to rerun and against what published benchmark. **Read this before the
first manuscript submission.***

---

## Severity key

- 🔴 **Critical** — must fix before submission. A reviewer will reject
  or return for major revision.
- 🟡 **Major** — should fix before submission or have a defensible
  answer for rebuttal.
- 🟢 **Minor** — nice to address; strengthens the paper but unlikely
  to tank review.

---

## 1. Differential expression methodology 🔴

### What the pipeline does

- `analysis.py` and several downstream scripts run **Welch's *t*-test on
  `log2(CPM+1)`** with BH FDR and `|log₂FC| ≥ 1`.

### Field standard (what reviewers will expect)

Bulk RNA-seq DE in plant-virology papers published in Nature Plants,
Plant Cell, New Phytologist, and Plant Journal in the last 5 years
almost universally uses one of:

1. **DESeq2** (Love et al. 2014) — negative binomial GLM on raw counts
   with size-factor normalisation and `lfcShrink` (apeglm / ashr).
2. **edgeR** — quasi-likelihood F-test (`glmQLFit` / `glmQLFTest`) with
   TMM normalisation on raw counts.
3. **limma-voom** — mean-variance modelling on log-CPM, then empirical
   Bayes moderation with `eBayes`.

Welch's *t*-test on log-CPM **does not model the mean-variance trend of
count data**, so it under-powers low-count genes and over-rejects high-
count ones. It also does not shrink variance across genes, meaning the
per-gene SD is noisy with n ≈ 30–60 per group. Reviewers will ask why
a bulk-RNA-seq study with hundreds of samples used a generic 2-group
test rather than the field-standard negative-binomial / voom approach.

### Remediation

- **Rerun DE for Q1–Q6** with DESeq2 *and* limma-voom in parallel.
  Report overlap with the Welch-t results as a robustness check in
  the supplement.
- Expose the DE design matrix: `~ Sampling_Season + Maturity_stage +
  Infection_Group + Infection_Group:Maturity_stage` for Q2/Q4. Current
  analysis treats each pairwise comparison in isolation, losing power.
- Re-calculate the **908 universal DEGs** after the above.

**Expected outcome:** core findings (sHSP suppression, immune silence,
sugar rewiring, biphasic trajectory) will survive because the effect
sizes are large. But specific gene lists will shift. This is a 1–2 day
compute + write job.

---

## 2. Experimental-design reporting 🟡

### Issues

1. **No power calculation in the methods.** With 63 RB− controls vs
   264 RB+ infected, the group imbalance (4:1) is unusual. Reviewers
   will ask whether the control sample size was pre-specified.
2. **No allocation table.** How many samples per (season × stage ×
   infection-cohort) cell? A ≥ 12-row allocation table is standard.
3. **Repeated-measures structure not modelled.** Several samples come
   from the same vine across stages / seasons. The current DE model
   treats every sample as independent. For the Q2 stage-based
   comparisons and Q6 trajectory analysis, this is a genuine
   pseudoreplication concern.
4. **Single vineyard, single cultivar** (Merlot) — a field limitation,
   not a pipeline issue, but should be explicitly flagged in the
   Discussion with "limits generalisability" language.

### Remediation

- Add allocation table to supplement.
- Refit DE with `~ condition + (1 | vine_id)` or equivalent in `lmerSeq`
  / limma-voom with `duplicateCorrelation`. If vine identity isn't
  available, document why and use block bootstrapping in a sensitivity
  check.
- Explicit "limits of generalisation" paragraph in the Discussion.

---

## 3. Read-level QC and alignment not documented in the repo 🔴

### Issue

The Methods section a reviewer would expect covers:

- Read QC (FastQC, MultiQC)
- Adapter trimming (Trim Galore / fastp)
- Alignment / quantification (STAR + featureCounts, Salmon / Kallisto,
  or HISAT2 + StringTie)
- rRNA / chloroplast / GRBV decontamination
- Alignment rate, assigned-reads rate, mean read count per sample
- Multi-QC plots for the cohort

**None of this is in the repository or documented.** Reviewers will
demand it. Even if the work was outsourced (commercial quant), the
pipeline version, reference genome + annotation, and alignment
statistics must be reported.

### Remediation

- Recover or reconstruct the QC/alignment metadata. Generate a
  MultiQC report for the 327 libraries.
- State the reference genome (NCBI GCF_030704535.1 is what the
  extended pipeline uses — confirm this matches the original quant)
  and annotation version explicitly.
- Deposit raw FASTQ or at least counts + metadata to GEO / SRA.

---

## 4. Causal mediation analysis as published 🔴

### Issue

The main paper's `causal_mediation_analysis.py` filters to **top-500
variable genes** before the mediation regression (line 253: `nlargest(500)`).
The sensitivity sweep (SENSITIVITY_REPORT.md) showed **zero** genes pass
the joint-significance test at tv=500. Signal only emerges at tv=2k–5k,
with `LOC100249249` as the top mediator (p = 1.2 × 10⁻¹⁰).

**If Figs 49–50 in the existing report show specific mediator names, those
names may not be reproducible with the script as written.** A reviewer
who reruns the code will see zero significant mediators and will file a
reproducibility failure.

### Remediation

- Re-run `causal_mediation_analysis.py` at tv=2000 (or more). Regenerate
  Figs 49–50. Name the stable mediators (cross-referenced to the HIMA
  sweep in the extended pipeline): `LOC100249249` is the cross-pipeline
  convergent hit.
- Document in Methods: "Top-variable gene filter at 2,000 selected to
  balance multiple-testing burden against signal detection; sensitivity
  across 500–10,000 reported in Supplementary Fig. Sx."
- Include the full sensitivity table as a supplementary figure.

---

## 5. ML biomarker AUC=1.000 claim 🔴

### Issue

The "5 methods achieve AUC=1.000" headline is real for the stated
analysis — **but the feature pool includes viral transcripts
(`N761_gp1–3`) which are perfectly diagnostic of infection by
definition** (you cannot have a viral transcript without the virus).
Sensitivity sweep:

| Feature pool | AUC range across 5 methods × 4 thresholds |
|---|---|
| All variable genes (viral + host) | 1.000 in every cell |
| **Host-only** (viral genes removed) | **0.919 – 1.000** (LASSO at tv=200 is 0.919) |

At one sensible threshold (tv=200, LASSO, host-only): **AUC = 0.919**,
not 1.000. Reviewers who look at the feature-importance plot in Fig 46
will immediately see N761_gp1/gp3 dominating and ask exactly this
question.

### Remediation

- Reframe the AUC headline: "Host-transcriptome biomarkers alone achieve
  AUC = 0.92–1.00 across 5 methods; inclusion of viral transcripts pushes
  all methods to AUC = 1.000 as expected."
- Report both panels (all-vars and host-only) in Fig 46.
- For the 3-gene host-only panel, report the actual cross-validated AUC
  on host-only features, not the all-vars AUC with viral genes included.
- Consider leave-one-infection-cohort-out cross-validation as a stricter
  generalisation test.

---

## 6. No qPCR validation of DEGs 🔴

### Issue

Every RNA-seq paper in the target journals validates **at least 5–10 top
DEGs** by RT-qPCR on an independent sample aliquot or biological
replicate cohort. The repository contains no such validation.

### Remediation

- qPCR panel of **at least 10 key genes**:
  - 3 top host biomarkers (LOC100854143, LOC100268011, LOC100257913)
  - 3 sHSPs (to confirm suppression)
  - HSFA2 (master regulator)
  - 3 robust temperature-titer mediators (LOC100266309, LOC100232974,
    LOC100249249)
- Run on a **fresh sample aliquot**, not the same RNA pool used for
  RNA-seq. Ideally on vines from a subsequent season if possible.
- Expected scope: 6–8 weeks of lab time. If PI won't approve qPCR,
  a reviewer-safe fallback is to use a published external dataset
  (e.g. from the meta-analysis targets) to check concordance of the
  top-50 DEG directions. This is weaker but publishable.

---

## 7. GSEA and pathway enrichment framing 🟡

### Issue

- 270 significant gene set × condition pairs at FDR<0.05 with 500
  permutations — fine.
- But: what **universe** was used for Fisher ORA in
  `pathway_enrichment.py`? All 14,604 filtered genes, or all annotated
  *V. vinifera* genes? The choice changes p-values substantially.
- No multi-dataset normalisation: GSEA uses t-statistic-ranked genes,
  but if the DE method was Welch t (see #1), rankings will differ from
  DESeq2/edgeR rankings.

### Remediation

- State the GSEA ranking metric and universe explicitly in Methods.
- Rerun GSEA on DESeq2 / edgeR-based rankings (follows from #1).
- Confirm that KEGG:Ribosome and GO:translation initiation factor
  enrichment (your headline GSEA finding across 6/8 conditions)
  survives the re-ranking. It probably will — these are robust pathways
  — but confirm.

---

## 8. WGCNA methodology 🟡

### Issue

- 13 modules; soft power β chosen but **scale-free topology R² not
  reported**. Reviewers commonly ask for `R² > 0.8` at the chosen power.
- Module size distribution looks reasonable but no `mergeCutHeight` /
  `deepSplit` / `minModuleSize` reported in methods.
- **The MEs.csv on disk uses sample IDs that don't match counts_matrix.txt**
  (`16H1` vs `S01` — documented in CLAUDE.md). This is a legacy artifact.
  All extended-pipeline downstream analyses (preservation, subtype, etc.)
  recompute modules from counts. Document this explicitly so a reviewer
  doesn't get confused trying to cross-reference.

### Remediation

- Report scale-free topology R² at chosen β in the methods.
- Report all WGCNA parameters (deepSplit, minModuleSize, mergeCutHeight,
  networkType).
- Regenerate MEs.csv using the counts_matrix sample IDs so everything
  in the repo is self-consistent.

---

## 9. Multiple-testing across the whole paper 🟡

### Issue

Each of 6 research questions does its own BH FDR within that question.
If the reader is assessing family-wise evidence across 6 questions × up
to 8 conditions × 14,604 genes, global FDR is not controlled. A stringent
reviewer may ask for a hierarchical or stage-level adjustment.

### Remediation

- State clearly: "FDR controlled within each comparison." Usually
  sufficient, but flag it.
- For the headline claims (sHSP suppression, mediator panel), report
  fold-change + BH q-value together so readers can calibrate.

---

## 10. QTL co-localisation claim 🟡

### Issue

"60% biomarker-QTL co-localisation in 1 Mb windows" needs a null. What
fraction of 45 random gene sets of matched size would co-localise with
121 QTL intervals by chance? Without a permutation-based null, the 60%
could be anywhere from "highly enriched" to "indistinguishable from
random" depending on the density of QTL coverage on the Vitis genome.

### Remediation

- Permutation test: draw 1,000 random 45-gene sets from the same
  chromosomes as the biomarkers; compute the null distribution of
  1-Mb co-localisations; report empirical p and fold-enrichment.
- Expected: the enrichment will be real but possibly smaller than 60%
  in effect size. That's fine — just report honestly.

---

## 11. Gene annotations are mostly LOC identifiers 🟡

### Issue

Top biomarkers, mediators, and hub TFs are reported as NCBI `LOC`
identifiers. This is a pure reproducibility failure at review time —
reviewers can't interpret `LOC100249249` without looking it up. For
any published paper, every gene named in the Abstract / Figures should
have its curated annotation in the main-text caption.

### Remediation

- Pull the NCBI gene description for every headline gene (biomarkers,
  mediators, hub TFs, top pseudotime genes) and add to a supplementary
  gene-annotation table. `efetch gene -id <GeneID> -format xml` or the
  NCBI datasets CLI can batch this.
- Replace bare `LOC` references in main-text figures and tables with
  `LOC100249249 (heat-shock cognate 70-2-like)` style annotations
  (example — fill with actual gene names after the annotation pass).

---

## 12. Causal language in the manuscript 🟡

### Issue

"Causal mediation" and similar terms appear in the pipeline and docs.
Observational mediation with n ≈ 60 vines per season is **causal under
specific assumptions** (no unmeasured confounders, correct DAG, no
reverse causation). Plant-virology reviewers are mixed on whether to
accept "causal" without an intervention.

### Remediation

- Use "statistical mediation" or "objective-inference mediation" in
  the Methods and Results unless and until the VIGS experiment is run.
- In the Discussion, frame the mediator panel as "high-priority
  candidates for functional validation" rather than "confirmed causal
  mediators."
- If aiming at PNAS, the VIGS experiment upgrades the language to
  "causal" legitimately.

---

## 13. Missing external validation 🟡

### Issue

Everything is internal to one 327-sample cohort from Merlot at a single
vineyard. Even within the dataset, 5-fold CV (Fig 46) uses random sample
assignment, which is not a generalisation test — samples from the same
vine or the same year end up on both sides of the split.

### Remediation

- **Leave-one-cohort-out (LOCO) CV**: 7 infection cohorts × 2 seasons =
  up to 14 held-out groups. Train on 13, test on the 14th. Report AUC
  across folds. This is the actual "generalisation" test. If AUC drops
  substantially, the paper's diagnostic claim softens accordingly.
- **External validation with published Vitis stress datasets** — this
  is what rules 11–12 (Tier 1 meta-analysis, scaffolded) will deliver.
  Priority: run Tier 1 before submission.

---

## 14. Figure-level issues a reviewer will flag 🟢

Quick-scan issues likely to come up in review:

- **Color choice**: verify color-blind safe palettes for all heatmaps
  (viridis / RdBu_r is fine; tab10 / jet is not).
- **Axis ranges inconsistent** across volcano plots in Figs 20, 22, 26
  — match them or state why they differ.
- **Gene names illegible** at submission DPI. Recompute at final
  journal size before submission.
- **Figure 14 master pathway heatmap** is information-dense; consider
  splitting into two or moving to a supplement.
- **Network figure (Fig 53–55)** — reviewers often want the full edge
  list as a supplementary TSV. Available? Confirm.

---

## 15. The sensitivity sweep itself is a strength 🟢

Not an issue, but worth calling out: the fact that every threshold
bottleneck has now been swept is **above-average rigour** for a plant
RNA-seq paper. Include the sensitivity report as a supplementary file
and cite it in the Methods. This converts a weakness (threshold
sensitivity of specific gene claims) into a demonstration of
methodological care.

---

## Summary — pre-submission checklist

Order of priority:

1. 🔴 Refit DE with DESeq2 or limma-voom; regenerate all downstream
   gene lists (§1).
2. 🔴 Reproduce QC + alignment stats; produce MultiQC; deposit to GEO
   (§3).
3. 🔴 Rerun `causal_mediation_analysis.py` at tv ≥ 2,000; regenerate
   Figs 49–50 (§4).
4. 🔴 Reframe the AUC=1.000 claim with the host-only caveat (§5).
5. 🔴 Perform qPCR validation of ≥ 10 key genes (§6).
6. 🟡 Annotate every gene named in the paper by its NCBI product (§11).
7. 🟡 LOCO cross-validation for the biomarker panel (§13).
8. 🟡 QTL co-localisation permutation test (§10).
9. 🟡 Add repeated-measures term to the DE model (§2).
10. 🟡 Soften "causal" language to "statistical mediation" until VIGS
    is done (§12).
11. 🟡 Report WGCNA parameters and scale-free R² (§8).
12. 🟡 Tier-1 meta-analysis (rules 11–12) as external validation (§13).
13. 🟢 Figure typography / palette pass (§14).

**Time estimate to address everything above**: 4–6 weeks focused work
on the compute side, plus 6–8 weeks of lab time for qPCR (parallelisable
with analysis work).

**After these fixes, the paper is Nature Plants-competitive as-is. With
the HSFA2 VIGS experiment on top, it's PNAS-competitive.**

---

## Appendix: journal-specific red lines

**Nature Plants.** Expects state-of-the-art statistical methods (DESeq2 /
limma-voom mandatory in practice), qPCR validation, transparent data
availability, clean Methods. VIGS not strictly required if the causal
inference is rigorous.

**Plant Cell.** Similar but more tolerant of a purely bioinformatic
paper with solid qPCR. A functional experiment strengthens the case
but is not mandatory for landing.

**New Phytologist / Plant Journal.** More permissive on methods; qPCR
still expected. The current pipeline, after fixes 1–6 above, fits
well here.

**PNAS.** Effectively mandatory functional validation (VIGS, CRISPR,
heterologous expression, or biochemical assay) for plant-virology
papers. Without it, the paper will struggle at editorial triage.

**BMC Genomics / MPMI.** Forgiving on methodology; the paper is
submission-ready for these journals after fixes 1, 3, 4, 5 (the four
red items).

---

*This critical review was generated by comparing the existing analysis
pipeline against the methodological conventions of published plant
RNA-seq studies in the target journals over 2020–2025. It is
intentionally pessimistic (a reviewer's first pass) and should be
treated as a pre-submission stress test, not a final judgement.*
