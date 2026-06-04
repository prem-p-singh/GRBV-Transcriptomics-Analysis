"""
rumbaugh2023_concordance.py

Within-cultivar mechanism replication: Rumbaugh et al. 2023 GRBV / Merlot /
Paso Robles cell-wall + phenylpropanoid + PR + hexose-transporter gene set
vs the GRBV Merlot project's covariate-adjusted DESeq2 outputs.

Why this comparison matters
---------------------------
Rumbaugh 2023 is the only other published GRBV transcriptome in Merlot.
Same cultivar (Merlot clone 12), different rootstock (1103P), different
vineyard (Paso Robles), different sampling year (2019). Their published
Table 1 carries 66 genes covering cell-wall metabolism (39 genes),
pathogenesis-related protein synthesis (11 genes), phenylpropanoid
pathway (11 genes), and hexose transporters (5 genes), with stage-resolved
log2FC at pre-veraison, veraison, post-veraison, and harvest.

What this script does
---------------------
1. Reads Rumbaugh 2023 Table 1 from the literature CSV
   (`farm_pipeline/data/literature/Rumbaugh2023_GRBV_Merlot_PasoRobles.csv`).
2. Joins each gene to the project's Q2x per-stage-per-season DESeq2 contrasts,
   Q2 per-stage season-pooled contrasts, Q6 pooled overall, the five Q1b
   per-cohort contrasts, and the Q2 / Q6 repeated-measures variants.
3. Computes Spearman rank correlation, Pearson correlation, and sign
   agreement (overall and restricted to genes with our padj < 0.05) per
   functional category and per stage.
4. Writes the results as tidy CSVs at the repo root.

Sign-agreement threshold: padj < 0.05.

Outputs
-------
  rumbaugh2023_gene_concordance.csv      - per-gene Rumbaugh + project log2FC side-by-side
  rumbaugh2023_per_stage_concordance.csv - per-(category, stage) Spearman/sign agreement
  rumbaugh2023_q1b_concordance.csv       - per-(category, cohort) cohort-gradient
  rumbaugh2023_q2q6_concordance.csv      - per-(category, rollup) Q2/Q6 concordance
  rumbaugh2023_rm_robustness.csv         - same under vine random-intercept variant
  rumbaugh2023_run_manifest.csv          - reproducibility manifest

Rules.md compliance
-------------------
Sec 1 (no hallucination): every claim grounds in file output.
Sec 5 (reproducible analyses): constants pinned at top, sorted iteration,
  manifest with reconciliation counts.
Sec 7 (change protocol): no figure rendered, no atlas rebuild, no git ops.
Sec 11 (plain language): the prose in this docstring and the printed
  summaries below follows the plain-language convention for biologist
  readers.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr


# =============================================================================
# Constants (pinned)
# =============================================================================

PROJECT_ROOT   = Path("/Users/prempratapsingh/Data_Analysis")
LITERATURE_DIR = PROJECT_ROOT / "farm_pipeline" / "data" / "literature"
RESULTS_DIR    = PROJECT_ROOT / "farm_pipeline" / "results"
OUTPUT_DIR     = PROJECT_ROOT

RUMBAUGH_FILE  = LITERATURE_DIR / "Rumbaugh2023_GRBV_Merlot_PasoRobles.csv"
PADJ_THRESHOLD = 0.05

STAGES = ["Pre_veraison", "Veraison", "Post_veraison", "Harvest"]
RUMBAUGH_STAGE_COLS = {
    "Pre_veraison":  "log2fc_PV",
    "Veraison":      "log2fc_V",
    "Post_veraison": "log2fc_POV",
    "Harvest":       "log2fc_H",
}

# Project DESeq2 contrast file maps
Q2X_FILES = {
    (yr, stage): RESULTS_DIR / "de_extras" / f"Q2x_{yr}_{stage}_deseq2.csv"
    for yr in ("2021", "2022") for stage in STAGES
}
Q2_FILES = {
    stage: RESULTS_DIR / "de_redo" / f"Q2_{stage}_deseq2.csv"
    for stage in STAGES
}
Q6_FILE = RESULTS_DIR / "de_redo" / "Q6_pooled_deseq2.csv"
Q1B_COHORTS = ["cohort2016", "cohort2019", "cohort2020", "cohort2021", "cohort2022"]
Q1B_FILES = {
    c: RESULTS_DIR / "de_extras" / f"Q1b_{c}_vs_neg_deseq2.csv"
    for c in Q1B_COHORTS
}
Q2_RM_FILES = {
    stage: RESULTS_DIR / "de_repeated_measures" / f"Q2_{stage}_RM.csv"
    for stage in STAGES
}
Q6_RM_FILE = RESULTS_DIR / "de_repeated_measures" / "Q6_pooled_RM.csv"


# =============================================================================
# Loaders
# =============================================================================

def load_deseq2(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df[["gene", "log2FC", "padj"]].rename(columns={
        "gene": "loc_id", "log2FC": "log2fc_ours", "padj": "padj_ours"
    })


def load_rumbaugh():
    df = pd.read_csv(RUMBAUGH_FILE)
    df["loc_id"] = df["loc_id"].astype("string")
    return df


# =============================================================================
# Concordance metric
# =============================================================================

def concordance(rumbaugh_vec, ours_vec, padj_vec, padj_cut=PADJ_THRESHOLD):
    rv = pd.Series(rumbaugh_vec).astype(float)
    ov = pd.Series(ours_vec).astype(float)
    pv = pd.Series(padj_vec).astype(float) if padj_vec is not None else None
    mask = rv.notna() & ov.notna()
    n = int(mask.sum())
    out = {
        "n": n,
        "spearman_rho": np.nan,
        "pearson_r":    np.nan,
        "sign_agreement_all": np.nan,
        "n_sig_padj":   0,
        "sign_agreement_sig": np.nan,
    }
    if n >= 3:
        r = rv[mask].values
        o = ov[mask].values
        try:
            rho, _ = spearmanr(r, o)
            out["spearman_rho"] = float(rho)
        except Exception:
            pass
        try:
            pear, _ = pearsonr(r, o)
            out["pearson_r"] = float(pear)
        except Exception:
            pass
        out["sign_agreement_all"] = float((np.sign(r) == np.sign(o)).mean())
    if pv is not None:
        sig_mask = mask & pv.notna() & (pv < padj_cut)
        n_sig = int(sig_mask.sum())
        out["n_sig_padj"] = n_sig
        if n_sig >= 3:
            r_sig = rv[sig_mask].values
            o_sig = ov[sig_mask].values
            out["sign_agreement_sig"] = float((np.sign(r_sig) == np.sign(o_sig)).mean())
    return out


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 78)
    print("Rumbaugh 2023 (Merlot Paso Robles) -> GRBV Merlot project concordance")
    print("=" * 78)

    manifest = []

    # ---- Read Rumbaugh 2023 Table 1 -----------------------------------
    rumb = load_rumbaugh()
    print(f"\n[A] Loaded Rumbaugh 2023 Table 1: {len(rumb)} genes")
    n_mapped = rumb["loc_id"].notna().sum()
    print(f"    mapped to LOC: {n_mapped}  (unmapped: {len(rumb) - n_mapped})")
    manifest.append({"stage": "input", "label": "Rumbaugh2023_Table1",
                     "path": str(RUMBAUGH_FILE), "n_rows": len(rumb),
                     "n_mapped": int(n_mapped)})

    # Drop unmapped (no LOC ID)
    rumb_m = rumb.dropna(subset=["loc_id"]).copy()

    # Category counts
    print("\n    Category breakdown (after dropping unmapped):")
    for cat, sub in rumb_m.groupby("category"):
        print(f"      {cat:48s}: {len(sub):>3d} genes")

    # ---- Pivot to long ------------------------------------------------
    rumb_long = rumb_m.melt(
        id_vars=["vit_id", "loc_id", "description", "category"],
        value_vars=list(RUMBAUGH_STAGE_COLS.values()),
        var_name="rumb_stage_col",
        value_name="log2fc_rumbaugh"
    )
    inv_stage = {v: k for k, v in RUMBAUGH_STAGE_COLS.items()}
    rumb_long["stage"] = rumb_long["rumb_stage_col"].map(inv_stage)
    rumb_long = rumb_long.dropna(subset=["log2fc_rumbaugh"])
    print(f"\n    Rumbaugh long table (cells with significant log2FC): {len(rumb_long)} rows")

    # ---- Load all DESeq2 outputs --------------------------------------
    print("\n[B] Loading project DESeq2 outputs")
    q2x_frames = []
    for (yr, stage), path in sorted(Q2X_FILES.items()):
        f = load_deseq2(path); f["our_year"] = yr; f["stage"] = stage
        q2x_frames.append(f)
    q2x = pd.concat(q2x_frames, ignore_index=True)

    q2_frames = []
    for stage, path in sorted(Q2_FILES.items()):
        f = load_deseq2(path); f["stage"] = stage
        q2_frames.append(f)
    q2 = pd.concat(q2_frames, ignore_index=True)

    q6 = load_deseq2(Q6_FILE)

    q1b_frames = []
    for c, path in sorted(Q1B_FILES.items()):
        f = load_deseq2(path); f["cohort"] = c
        q1b_frames.append(f)
    q1b = pd.concat(q1b_frames, ignore_index=True)

    q2_rm_frames = []
    for stage, path in sorted(Q2_RM_FILES.items()):
        f = load_deseq2(path); f["stage"] = stage
        q2_rm_frames.append(f)
    q2_rm = pd.concat(q2_rm_frames, ignore_index=True)

    q6_rm = load_deseq2(Q6_RM_FILE)

    print(f"  Q2x  : {q2x['loc_id'].nunique()} genes x 8 cells")
    print(f"  Q2   : {q2['loc_id'].nunique()} genes x 4 stages")
    print(f"  Q6   : {q6['loc_id'].nunique()} genes")
    print(f"  Q1b  : {q1b['loc_id'].nunique()} genes x 5 cohorts")
    print(f"  Q2 RM: {q2_rm['loc_id'].nunique()} genes x 4 stages")
    print(f"  Q6 RM: {q6_rm['loc_id'].nunique()} genes")

    expressed = set(q6["loc_id"].dropna().unique())

    # ---- Phase C: per-gene wide table (rumbaugh side by side with ours) ----
    print("\n[C] Building per-gene wide table")

    # Pivot our Q2 (4 cells; season-pooled, most direct match to Rumbaugh's
    # single-season design) to wide
    q2_wide = (q2[["loc_id", "stage", "log2fc_ours", "padj_ours"]]
        .pivot_table(index="loc_id", columns="stage",
                     values=["log2fc_ours", "padj_ours"]))
    q2_wide.columns = [f"ours_Q2_{stage}_{metric}"
                       for metric, stage in q2_wide.columns]
    q2_wide = q2_wide.reset_index()

    # Pivot our Q2x (8 cells; year-stage) to wide
    q2x_wide = (q2x[["loc_id", "our_year", "stage", "log2fc_ours", "padj_ours"]]
        .pivot_table(index="loc_id",
                     columns=["our_year", "stage"],
                     values=["log2fc_ours", "padj_ours"]))
    q2x_wide.columns = [f"ours_Q2x_{year}_{stage}_{metric}"
                        for metric, year, stage in q2x_wide.columns]
    q2x_wide = q2x_wide.reset_index()

    wide = rumb_m.merge(q2_wide,  on="loc_id", how="left")
    wide = wide.merge(q2x_wide, on="loc_id", how="left")
    wide["in_expressed_background"] = wide["loc_id"].isin(expressed)
    print(f"  Per-gene wide table: {len(wide)} rows, "
          f"{int(wide['in_expressed_background'].sum())} in our expressed-gene background")

    # ---- Phase D: per-(category, stage) concordance vs Q2 (season-pooled) ---
    print("\n[D] Per-(category, stage) concordance against Q2 (season-pooled)")

    cat_stage_rows = []
    q2_keyed = q2[["loc_id", "stage", "log2fc_ours", "padj_ours"]]
    merged_q2 = rumb_long.merge(q2_keyed, on=["loc_id", "stage"], how="inner")
    for (cat, stage), grp in merged_q2.groupby(["category", "stage"], sort=True):
        m = concordance(grp["log2fc_rumbaugh"], grp["log2fc_ours"], grp["padj_ours"])
        m.update({"category": cat, "stage": stage,
                  "contrast": f"Q2_{stage}"})
        cat_stage_rows.append(m)

    # Per-(category, stage, our_year) concordance vs Q2x
    q2x_keyed = q2x[["loc_id", "our_year", "stage", "log2fc_ours", "padj_ours"]]
    merged_q2x = rumb_long.merge(q2x_keyed, on=["loc_id", "stage"], how="inner")
    for (cat, stage, yr), grp in merged_q2x.groupby(
            ["category", "stage", "our_year"], sort=True):
        m = concordance(grp["log2fc_rumbaugh"], grp["log2fc_ours"], grp["padj_ours"])
        m.update({"category": cat, "stage": stage,
                  "contrast": f"Q2x_{yr}_{stage}"})
        cat_stage_rows.append(m)

    per_stage = pd.DataFrame(cat_stage_rows)[
        ["category", "stage", "contrast", "n", "spearman_rho", "pearson_r",
         "sign_agreement_all", "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  Per-(category, stage, contrast) rows: {len(per_stage)}")

    # ---- Phase E: per-category overall (across all stages flattened) -------
    print("\n[E] Per-category overall concordance (flattened across stages)")
    overall_rows = []
    for cat, grp in merged_q2.groupby("category", sort=True):
        m = concordance(grp["log2fc_rumbaugh"], grp["log2fc_ours"], grp["padj_ours"])
        m.update({"category": cat, "contrast": "Q2_flattened",
                  "n_unique_genes": grp["loc_id"].nunique(),
                  "n_cells": grp.groupby("stage").ngroups})
        overall_rows.append(m)
    # And the all-network total
    m_all = concordance(merged_q2["log2fc_rumbaugh"],
                        merged_q2["log2fc_ours"],
                        merged_q2["padj_ours"])
    m_all.update({"category": "ALL",
                  "contrast": "Q2_flattened",
                  "n_unique_genes": merged_q2["loc_id"].nunique(),
                  "n_cells": merged_q2.groupby("stage").ngroups})
    overall_rows.append(m_all)
    overall = pd.DataFrame(overall_rows)[
        ["category", "contrast", "n_cells", "n_unique_genes",
         "n", "spearman_rho", "pearson_r", "sign_agreement_all",
         "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  Per-category overall rows: {len(overall)}")
    print(overall.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # ---- Phase F: Q1b cohort gradient using Rumbaugh harvest column -------
    print("\n[F] Q1b per-cohort concordance (Rumbaugh harvest log2FC used)")
    rumb_harvest = rumb_m[["category", "vit_id", "loc_id", "log2fc_H"]].rename(
        columns={"log2fc_H": "rumb_harvest_log2fc"}).dropna(subset=["rumb_harvest_log2fc"])

    q1b_rows = []
    for cohort, cohort_df in sorted(q1b.groupby("cohort")):
        joined = rumb_harvest.merge(
            cohort_df[["loc_id", "log2fc_ours", "padj_ours"]],
            on="loc_id", how="inner")
        for cat, grp in joined.groupby("category", sort=True):
            m = concordance(grp["rumb_harvest_log2fc"], grp["log2fc_ours"], grp["padj_ours"])
            m.update({"category": cat, "cohort": cohort,
                      "n_unique_genes": grp["loc_id"].nunique()})
            q1b_rows.append(m)
        # ALL category
        m_all = concordance(joined["rumb_harvest_log2fc"], joined["log2fc_ours"], joined["padj_ours"])
        m_all.update({"category": "ALL", "cohort": cohort,
                      "n_unique_genes": joined["loc_id"].nunique()})
        q1b_rows.append(m_all)
    q1b_concord = pd.DataFrame(q1b_rows)[
        ["category", "cohort", "n_unique_genes", "n", "spearman_rho",
         "pearson_r", "sign_agreement_all", "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  Q1b concordance rows: {len(q1b_concord)}")

    # ---- Phase G: Q6 pooled + Q2 per-stage rollup -------------------------
    print("\n[G] Q6 pooled + Q2 per-stage rollup")
    rollup_rows = []

    # Q6 pooled (collapse Rumbaugh to harvest)
    for cat, grp_r in rumb_harvest.groupby("category"):
        merged_q6 = grp_r.merge(q6[["loc_id", "log2fc_ours", "padj_ours"]],
                                on="loc_id", how="inner")
        m = concordance(merged_q6["rumb_harvest_log2fc"],
                        merged_q6["log2fc_ours"], merged_q6["padj_ours"])
        m.update({"category": cat, "contrast": "Q6_pooled",
                  "n_unique_genes": grp_r["loc_id"].nunique()})
        rollup_rows.append(m)
    merged_q6_all = rumb_harvest.merge(q6[["loc_id", "log2fc_ours", "padj_ours"]],
                                       on="loc_id", how="inner")
    m_all = concordance(merged_q6_all["rumb_harvest_log2fc"],
                        merged_q6_all["log2fc_ours"], merged_q6_all["padj_ours"])
    m_all.update({"category": "ALL", "contrast": "Q6_pooled",
                  "n_unique_genes": merged_q6_all["loc_id"].nunique()})
    rollup_rows.append(m_all)

    rollup = pd.DataFrame(rollup_rows)[
        ["category", "contrast", "n_unique_genes", "n", "spearman_rho",
         "pearson_r", "sign_agreement_all", "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  Q6 rollup rows: {len(rollup)}")

    # ---- Phase H: Repeated-measures robustness ----------------------------
    print("\n[H] Repeated-measures robustness")
    rm_rows = []
    rumb_long_for_rm = rumb_long.copy()
    merged_q2rm = rumb_long_for_rm.merge(
        q2_rm[["loc_id", "stage", "log2fc_ours", "padj_ours"]],
        on=["loc_id", "stage"], how="inner")
    for (cat, stage), grp in merged_q2rm.groupby(["category", "stage"], sort=True):
        m = concordance(grp["log2fc_rumbaugh"], grp["log2fc_ours"], grp["padj_ours"])
        m.update({"category": cat, "stage": stage,
                  "contrast": f"Q2_RM_{stage}"})
        rm_rows.append(m)
    # Flattened Q2 RM per category
    for cat, grp in merged_q2rm.groupby("category"):
        m = concordance(grp["log2fc_rumbaugh"], grp["log2fc_ours"], grp["padj_ours"])
        m.update({"category": cat, "stage": "ALL_STAGES",
                  "contrast": "Q2_RM_flattened"})
        rm_rows.append(m)
    # Q6 RM
    for cat, grp_r in rumb_harvest.groupby("category"):
        merged = grp_r.merge(q6_rm[["loc_id", "log2fc_ours", "padj_ours"]],
                             on="loc_id", how="inner")
        m = concordance(merged["rumb_harvest_log2fc"],
                        merged["log2fc_ours"], merged["padj_ours"])
        m.update({"category": cat, "stage": "Harvest_only",
                  "contrast": "Q6_pooled_RM"})
        rm_rows.append(m)
    rm = pd.DataFrame(rm_rows)[
        ["category", "stage", "contrast", "n", "spearman_rho", "pearson_r",
         "sign_agreement_all", "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  RM concordance rows: {len(rm)}")

    # ---- Phase I: Write outputs + manifest --------------------------------
    print("\n[I] Writing outputs to repo root")
    outputs = {
        "rumbaugh2023_gene_concordance.csv":       wide,
        "rumbaugh2023_per_stage_concordance.csv":  per_stage,
        "rumbaugh2023_per_category_summary.csv":   overall,
        "rumbaugh2023_q1b_concordance.csv":        q1b_concord,
        "rumbaugh2023_q2q6_concordance.csv":       rollup,
        "rumbaugh2023_rm_robustness.csv":          rm,
    }
    for fname, df in outputs.items():
        path = OUTPUT_DIR / fname
        df.to_csv(path, index=False)
        print(f"  wrote {fname}  ({df.shape[0]} rows, {df.shape[1]} cols)")
        manifest.append({"stage": "output", "label": fname,
                         "path": str(path), "n_rows": df.shape[0],
                         "n_mapped": df.shape[1]})

    # Reconciliation
    print("\n[I] Reconciliation per category:")
    for cat in sorted(rumb_m["category"].unique()):
        n_total = rumb_m[rumb_m["category"] == cat]["loc_id"].nunique()
        n_exp   = rumb_m[(rumb_m["category"] == cat)
                         & rumb_m["loc_id"].isin(expressed)]["loc_id"].nunique()
        print(f"  {cat:48s}: {n_total} mapped LOC, {n_exp} in expressed background")
        manifest.append({"stage": "reconcile", "label": f"category_{cat}",
                         "path": "", "n_rows": n_total, "n_mapped": n_exp})

    pd.DataFrame(manifest).to_csv(OUTPUT_DIR / "rumbaugh2023_run_manifest.csv",
                                  index=False)
    print(f"\n  wrote rumbaugh2023_run_manifest.csv ({len(manifest)} rows)")

    # ---- Phase J: Headline summary printed to stdout ---------------------
    print("\n" + "=" * 78)
    print("HEADLINE: Per-(category, stage) concordance against Q2 (season-pooled)")
    print("=" * 78)
    pivot_q2 = (per_stage[per_stage["contrast"].str.startswith("Q2_")
                          & ~per_stage["contrast"].str.startswith("Q2x")]
                .pivot_table(index="category", columns="stage",
                             values="spearman_rho", aggfunc="first"))
    print(pivot_q2.to_string(float_format=lambda v: f"{v:.3f}"))

    print("\n" + "=" * 78)
    print("HEADLINE: Q1b cohort gradient (Rumbaugh harvest log2FC used)")
    print("=" * 78)
    pivot_q1b = q1b_concord.pivot_table(
        index="category", columns="cohort",
        values="spearman_rho", aggfunc="first")
    print(pivot_q1b.to_string(float_format=lambda v: f"{v:.3f}"))

    print("\n" + "=" * 78)
    print("HEADLINE: Per-(category, stage, our_year) Q2x concordance")
    print("=" * 78)
    pivot_q2x = (per_stage[per_stage["contrast"].str.startswith("Q2x")]
                 .assign(our_year=lambda d: d["contrast"].str.extract(r"Q2x_(\d{4})_"))
                 .pivot_table(index="category", columns=["stage", "our_year"],
                              values="spearman_rho", aggfunc="first"))
    print(pivot_q2x.to_string(float_format=lambda v: f"{v:.3f}"))

    print("\nDone.")


if __name__ == "__main__":
    main()
