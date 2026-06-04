"""
rumbaugh2022_concordance.py

Cross-paper concordance between Rumbaugh et al. 2022 (CS multi-omics,
Int J Mol Sci 23, 13248) differential co-expression network gene sets
(Tables S1, S2, S3) and the GRBV Merlot project's covariate-adjusted
DESeq2 outputs.

Comparison framing
------------------
Rumbaugh 2022 provides per-gene log2FC at up to 14 cells per gene:
{Pre-Veraison, Veraison, Post-Veraison, Harvest}
x {110R, 420A} x {2016 cool, 2017 hot}
(Post-Veraison 2017 missing because of heat-wave sampling shift).

Cell mapping into the project (rootstocks kept separate per user
instruction):
  Rumbaugh 2016 (cool) <-> project 2021 (cool)
  Rumbaugh 2017 (hot)  <-> project 2022 (hot)

Three layers of DESeq2 contrasts feed this comparison:
  1. Q2x per-stage per-season (8 cells, primary per-cell match)
  2. Q2 per-stage season-pooled (4 cells, rollup sanity)
  3. Q6 pooled (1 value, overall sanity floor)
Plus:
  4. Q1b per-cohort (5 contrasts). Rumbaugh did not vary onset cohort,
     so for this layer we collapse Rumbaugh to a single per-gene
     log2FC per rootstock (harvest-mean across years) and compare that
     against each project cohort contrast separately. Tells us which
     cohort's transcriptome most resembles Rumbaugh's vines.
  5. Q2 / Q6 repeated-measures (vine random intercept) -- robustness
     against the project's non-independence finding.

Sign-agreement threshold: padj < 0.05 (in the project's DESeq2 output).

Outputs (repo root, alongside other cross-paper concordance CSVs):
  rumbaugh2022_network_members.csv      - wide per-gene table
  rumbaugh2022_concordance_per_cell.csv - per (network, rootstock, cell)
  rumbaugh2022_concordance_summary.csv  - per (network, rootstock)
  rumbaugh2022_q1b_concordance.csv      - per (network, rootstock, cohort)
  rumbaugh2022_q2q6_concordance.csv     - per (network, rootstock, rollup)
  rumbaugh2022_rm_robustness.csv        - same under RM variant
  rumbaugh2022_dcl2_trajectory.csv      - DCL2 single-gene side-by-side
  rumbaugh2022_cabp_trajectory.csv      - Ca-binding-protein single-gene
  rumbaugh2022_run_manifest.csv         - reproducibility manifest

Pipeline conventions followed (Rules.md):
  Sec 1 (no hallucination): every claim grounds in file output.
  Sec 5 (reproducible analyses): single parser, single loader,
    constants pinned at top, sorted iteration, run manifest
    with reconciliation counts.
  Sec 7 (change protocol): no figure rendered, no atlas rebuild,
    no git ops, no FARM dispatch.
"""

from pathlib import Path
import sys

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

PADJ_THRESHOLD = 0.05  # significance cut, per user

RUMBAUGH_TABLES = {
    "S1_CaBP_hub":  LITERATURE_DIR / "Rumbaugh et al 2022_Table S1.xlsx",
    "S2_DCL2_gain": LITERATURE_DIR / "Rumbaugh et al 2022_Table S2.xlsx",
    "S3_DCL2_loss": LITERATURE_DIR / "Rumbaugh et al 2022_Table S3.xlsx",
}

# Cell layout in Rumbaugh supplementary tables (col idx -> stage/rootstock/year).
# Verified by reading Tables S1/S2/S3 header bands.
RUMBAUGH_CELL_LAYOUT = [
    (4,  "Pre_veraison",  "110R", "2016"),
    (5,  "Pre_veraison",  "420A", "2016"),
    (6,  "Pre_veraison",  "110R", "2017"),
    (7,  "Pre_veraison",  "420A", "2017"),
    (8,  "Veraison",      "110R", "2016"),
    (9,  "Veraison",      "420A", "2016"),
    (10, "Veraison",      "110R", "2017"),
    (11, "Veraison",      "420A", "2017"),
    (12, "Post_veraison", "110R", "2016"),
    (13, "Post_veraison", "420A", "2016"),
    (14, "Harvest",       "110R", "2016"),
    (15, "Harvest",       "420A", "2016"),
    (16, "Harvest",       "110R", "2017"),
    (17, "Harvest",       "420A", "2017"),
]

# Rumbaugh sampling year -> matching project sampling year (climate-based)
RUMBAUGH_TO_OUR_YEAR = {"2016": "2021", "2017": "2022"}

STAGES = ["Pre_veraison", "Veraison", "Post_veraison", "Harvest"]

# DESeq2 contrast file maps
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
# Rumbaugh parser
# =============================================================================

def parse_rumbaugh_table(path: Path, network_label: str) -> pd.DataFrame:
    """Return long tidy frame: (network, vit_id, loc_id, description, pathway,
    rootstock, rumbaugh_year, our_year, stage, log2fc_rumbaugh)."""
    df = pd.read_excel(path, header=None)
    header_row = None
    for i in range(len(df)):
        v = df.iloc[i, 0]
        if isinstance(v, str) and v.strip() == "ensembl_gene_id":
            header_row = i
            break
    if header_row is None:
        raise ValueError(f"Header row 'ensembl_gene_id' not found in {path}")
    data_start = header_row + 2  # skip rootstock-year sub-header

    rows = []
    skipped_no_entrez = 0
    skipped_non_vit = 0
    for i in range(data_start, len(df)):
        row = df.iloc[i]
        vit = row.iloc[0]
        if not isinstance(vit, str) or not vit.startswith("VIT_"):
            skipped_non_vit += 1
            continue
        entrez = row.iloc[1]
        if pd.isna(entrez):
            skipped_no_entrez += 1
            continue
        try:
            loc_id = f"LOC{int(entrez)}"
        except (ValueError, TypeError):
            skipped_no_entrez += 1
            continue
        description = row.iloc[2] if pd.notna(row.iloc[2]) else ""
        pathway     = row.iloc[3] if pd.notna(row.iloc[3]) else ""
        for col_idx, stage, rootstock, year in RUMBAUGH_CELL_LAYOUT:
            val = row.iloc[col_idx]
            if pd.isna(val):
                continue
            rows.append({
                "network":         network_label,
                "vit_id":          vit,
                "loc_id":          loc_id,
                "description":     str(description),
                "pathway":         str(pathway),
                "rootstock":       rootstock,
                "rumbaugh_year":   year,
                "our_year":        RUMBAUGH_TO_OUR_YEAR[year],
                "stage":           stage,
                "log2fc_rumbaugh": float(val),
            })
    return (pd.DataFrame(rows), skipped_no_entrez, skipped_non_vit)


# =============================================================================
# DESeq2 loaders
# =============================================================================

def load_deseq2(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    out = df[["gene", "log2FC", "padj"]].rename(columns={
        "gene": "loc_id", "log2FC": "log2fc_ours", "padj": "padj_ours"
    })
    out["contrast"] = label
    return out


def load_all_q2x():
    frames = []
    for (yr, stage), path in sorted(Q2X_FILES.items()):
        f = load_deseq2(path, f"Q2x_{yr}_{stage}")
        f["our_year"] = yr
        f["stage"] = stage
        frames.append(f)
    return pd.concat(frames, ignore_index=True)


def load_all_q2():
    frames = []
    for stage, path in sorted(Q2_FILES.items()):
        f = load_deseq2(path, f"Q2_{stage}")
        f["stage"] = stage
        frames.append(f)
    return pd.concat(frames, ignore_index=True)


def load_q6():
    return load_deseq2(Q6_FILE, "Q6_pooled")


def load_all_q1b():
    frames = []
    for c, path in sorted(Q1B_FILES.items()):
        f = load_deseq2(path, f"Q1b_{c}")
        f["cohort"] = c
        frames.append(f)
    return pd.concat(frames, ignore_index=True)


def load_all_q2_rm():
    frames = []
    for stage, path in sorted(Q2_RM_FILES.items()):
        f = load_deseq2(path, f"Q2_RM_{stage}")
        f["stage"] = stage
        frames.append(f)
    return pd.concat(frames, ignore_index=True)


def load_q6_rm():
    return load_deseq2(Q6_RM_FILE, "Q6_pooled_RM")


# =============================================================================
# Concordance metric
# =============================================================================

def concordance(rumbaugh_vec, ours_vec, padj_vec, padj_cut=PADJ_THRESHOLD):
    """Spearman, Pearson, n, sign agreement (overall and padj-filtered)."""
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
# Phase A: load and parse
# =============================================================================

def main():
    print("=" * 78)
    print("Rumbaugh 2022 -> GRBV Merlot project concordance")
    print("=" * 78)

    manifest_rows = []

    # Parse Rumbaugh tables
    rumbaugh_long = []
    print("\n[A] Parsing Rumbaugh 2022 supplementary tables")
    for label, path in RUMBAUGH_TABLES.items():
        df, skip_no_entrez, skip_non_vit = parse_rumbaugh_table(path, label)
        n_genes = df["vit_id"].nunique()
        n_cells = len(df)
        print(f"  {label}: {n_genes} unique VIT IDs, "
              f"{n_cells} long-form (gene x cell) rows, "
              f"skipped {skip_no_entrez} (no entrez) + {skip_non_vit} (non-VIT)")
        manifest_rows.append({
            "stage": "input",
            "label": label,
            "path": str(path),
            "n_unique_genes_or_rows": n_genes,
            "n_total_rows": n_cells,
        })
        rumbaugh_long.append(df)
    rumbaugh_long = pd.concat(rumbaugh_long, ignore_index=True)
    print(f"  Combined Rumbaugh long table: {len(rumbaugh_long)} rows, "
          f"{rumbaugh_long['vit_id'].nunique()} unique VIT IDs across networks")

    # Load DESeq2 contrasts
    print("\n[A] Loading DESeq2 outputs")
    q2x = load_all_q2x()
    q2  = load_all_q2()
    q6  = load_q6()
    q1b = load_all_q1b()
    q2_rm = load_all_q2_rm()
    q6_rm = load_q6_rm()
    print(f"  Q2x (per-cell):      {q2x['loc_id'].nunique()} genes x "
          f"{q2x['contrast'].nunique()} contrasts = {len(q2x)} rows")
    print(f"  Q2  (per-stage):     {q2['loc_id'].nunique()} genes x "
          f"{q2['contrast'].nunique()} contrasts = {len(q2)} rows")
    print(f"  Q6  (overall):       {q6['loc_id'].nunique()} genes")
    print(f"  Q1b (per-cohort):    {q1b['loc_id'].nunique()} genes x "
          f"{q1b['cohort'].nunique()} cohorts = {len(q1b)} rows")
    print(f"  Q2 RM (per-stage):   {q2_rm['loc_id'].nunique()} genes x "
          f"{q2_rm['contrast'].nunique()} contrasts = {len(q2_rm)} rows")
    print(f"  Q6 RM (overall):     {q6_rm['loc_id'].nunique()} genes")

    expressed_set = set(q6["loc_id"].dropna().unique())
    print(f"\n  Background expressed-gene set (from Q6): {len(expressed_set)}")

    # =============================================================================
    # Phase B: per-cell concordance (Q2x, rootstocks separate)
    # =============================================================================

    print("\n[B] Per-cell concordance (Q2x, rootstocks separate)")

    # Join Rumbaugh long onto Q2x by (loc_id, our_year, stage)
    q2x_keyed = q2x[["loc_id", "our_year", "stage", "log2fc_ours", "padj_ours"]]
    merged = rumbaugh_long.merge(q2x_keyed, on=["loc_id", "our_year", "stage"], how="inner")
    print(f"  Joined rows: {len(merged)}")
    print(f"  Networks present after join: {sorted(merged['network'].unique())}")
    print(f"  Genes present after join: {merged['loc_id'].nunique()}")

    # Per-cell concordance: (network, rootstock, rumbaugh_year, our_year, stage)
    per_cell_rows = []
    for (net, rs, ryr, oyr, st), grp in merged.groupby(
            ["network", "rootstock", "rumbaugh_year", "our_year", "stage"], sort=True):
        m = concordance(grp["log2fc_rumbaugh"], grp["log2fc_ours"], grp["padj_ours"])
        m.update({"network": net, "rootstock": rs, "rumbaugh_year": ryr,
                  "our_year": oyr, "stage": st})
        per_cell_rows.append(m)
    per_cell = pd.DataFrame(per_cell_rows)[
        ["network", "rootstock", "rumbaugh_year", "our_year", "stage",
         "n", "spearman_rho", "pearson_r", "sign_agreement_all",
         "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  Per-cell concordance rows: {len(per_cell)}")

    # Per-(network, rootstock) flattened concordance across all cells
    summary_rows = []
    for (net, rs), grp in merged.groupby(["network", "rootstock"], sort=True):
        m = concordance(grp["log2fc_rumbaugh"], grp["log2fc_ours"], grp["padj_ours"])
        m.update({
            "network":      net,
            "rootstock":    rs,
            "n_cells":      grp.groupby(["rumbaugh_year", "stage"]).ngroups,
            "n_unique_genes": grp["loc_id"].nunique(),
        })
        summary_rows.append(m)
    summary = pd.DataFrame(summary_rows)[
        ["network", "rootstock", "n_cells", "n_unique_genes",
         "n", "spearman_rho", "pearson_r", "sign_agreement_all",
         "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  Flattened summary rows: {len(summary)}")

    # =============================================================================
    # Phase B2: wide network-members table
    # =============================================================================

    # Build wide layout: one row per (network, vit_id, loc_id, rootstock).
    # Columns: 7 Rumbaugh cells per rootstock (PV/V/PoV/H x 2016 + PV/V/H x 2017 -- 7 cells),
    # plus 8 of our cells (4 stages x 2 years).
    wide_meta = (rumbaugh_long
        .drop_duplicates(subset=["network", "vit_id", "loc_id", "rootstock"])
        [["network", "vit_id", "loc_id", "description", "pathway", "rootstock"]])

    # Pivot Rumbaugh log2fc to wide (one column per rumbaugh_year x stage)
    rumb_wide = (rumbaugh_long
        .assign(cell=lambda d: "rumb_" + d["rumbaugh_year"] + "_" + d["stage"])
        .pivot_table(
            index=["network", "vit_id", "loc_id", "rootstock"],
            columns="cell",
            values="log2fc_rumbaugh",
            aggfunc="first"
        ).reset_index())

    # Pivot our log2fc to wide (one column per our_year x stage). One value per
    # (loc_id, our_year, stage), independent of rootstock (our data has one).
    ours_wide = (q2x[["loc_id", "our_year", "stage", "log2fc_ours", "padj_ours"]]
        .assign(
            cell_log2=lambda d: "ours_" + d["our_year"] + "_" + d["stage"] + "_log2fc",
            cell_padj=lambda d: "ours_" + d["our_year"] + "_" + d["stage"] + "_padj",
        ))
    ours_log2 = ours_wide.pivot_table(index="loc_id", columns="cell_log2",
                                      values="log2fc_ours", aggfunc="first")
    ours_padj = ours_wide.pivot_table(index="loc_id", columns="cell_padj",
                                      values="padj_ours", aggfunc="first")
    ours_combined = ours_log2.merge(ours_padj, left_index=True, right_index=True,
                                    how="outer").reset_index()

    wide = wide_meta.merge(rumb_wide,
                           on=["network", "vit_id", "loc_id", "rootstock"],
                           how="left")
    wide = wide.merge(ours_combined, on="loc_id", how="left")
    wide["in_expressed_background"] = wide["loc_id"].isin(expressed_set)
    print(f"  Wide network-members table: {len(wide)} rows "
          f"(unique network x gene x rootstock combos)")

    # =============================================================================
    # Phase C: Q1b cohort concordance (Rumbaugh collapsed to one value per gene per rootstock)
    # =============================================================================

    print("\n[C] Q1b per-cohort concordance")

    # Collapse Rumbaugh to one log2FC per (gene, rootstock) using:
    #   (a) mean across all available cells (overall_mean)
    #   (b) mean across harvest cells only (harvest_mean) -- closest stable cell
    rumb_collapsed = (rumbaugh_long.groupby(
        ["network", "vit_id", "loc_id", "rootstock"])
        .agg(rumb_mean_all=("log2fc_rumbaugh", "mean"))
        .reset_index())
    harvest_only = (rumbaugh_long[rumbaugh_long["stage"] == "Harvest"]
        .groupby(["network", "vit_id", "loc_id", "rootstock"])
        .agg(rumb_mean_harvest=("log2fc_rumbaugh", "mean"))
        .reset_index())
    rumb_collapsed = rumb_collapsed.merge(
        harvest_only,
        on=["network", "vit_id", "loc_id", "rootstock"],
        how="left")

    q1b_rows = []
    for cohort, cohort_df in sorted(q1b.groupby("cohort")):
        joined = rumb_collapsed.merge(
            cohort_df[["loc_id", "log2fc_ours", "padj_ours"]],
            on="loc_id", how="inner")
        for (net, rs), grp in joined.groupby(["network", "rootstock"], sort=True):
            # Using rumb_mean_all
            m_all = concordance(grp["rumb_mean_all"], grp["log2fc_ours"], grp["padj_ours"])
            # Using rumb_mean_harvest
            m_hv  = concordance(grp["rumb_mean_harvest"], grp["log2fc_ours"], grp["padj_ours"])
            q1b_rows.append({
                "network":   net,
                "rootstock": rs,
                "cohort":    cohort,
                "n_genes":   grp["loc_id"].nunique(),
                "rho_meanAll": m_all["spearman_rho"],
                "signAgr_meanAll": m_all["sign_agreement_all"],
                "n_sig_meanAll":   m_all["n_sig_padj"],
                "signAgr_sig_meanAll": m_all["sign_agreement_sig"],
                "rho_harvest": m_hv["spearman_rho"],
                "signAgr_harvest": m_hv["sign_agreement_all"],
                "n_sig_harvest":   m_hv["n_sig_padj"],
                "signAgr_sig_harvest": m_hv["sign_agreement_sig"],
            })
    q1b_concord = pd.DataFrame(q1b_rows)
    print(f"  Q1b concordance rows: {len(q1b_concord)}")

    # =============================================================================
    # Phase D: Q2 (per-stage pooled) + Q6 (overall) rollups
    # =============================================================================

    print("\n[D] Q2 (per-stage pooled) and Q6 (overall) rollup concordance")

    # For Q2, collapse Rumbaugh to per-(gene, rootstock, stage) average across years
    rumb_q2 = (rumbaugh_long.groupby(["network", "vit_id", "loc_id", "rootstock", "stage"])
        .agg(rumb_mean=("log2fc_rumbaugh", "mean"))
        .reset_index())
    q2_rows = []
    for (st, rs, net), grp_r in rumb_q2.groupby(["stage", "rootstock", "network"]):
        ours_grp = q2[q2["stage"] == st][["loc_id", "log2fc_ours", "padj_ours"]]
        merged_st = grp_r.merge(ours_grp, on="loc_id", how="inner")
        m = concordance(merged_st["rumb_mean"], merged_st["log2fc_ours"], merged_st["padj_ours"])
        m.update({"network": net, "rootstock": rs, "contrast": f"Q2_{st}", "n_genes": len(merged_st)})
        q2_rows.append(m)

    # Q6: collapse Rumbaugh to per-(gene, rootstock) average across all cells
    q6_rows = []
    for (net, rs), grp_r in rumb_collapsed.groupby(["network", "rootstock"]):
        merged_q6 = grp_r.merge(q6[["loc_id", "log2fc_ours", "padj_ours"]],
                                on="loc_id", how="inner")
        m = concordance(merged_q6["rumb_mean_all"], merged_q6["log2fc_ours"], merged_q6["padj_ours"])
        m.update({"network": net, "rootstock": rs, "contrast": "Q6_pooled", "n_genes": len(merged_q6)})
        q6_rows.append(m)

    q2q6 = pd.DataFrame(q2_rows + q6_rows)[
        ["network", "rootstock", "contrast", "n_genes", "n", "spearman_rho",
         "pearson_r", "sign_agreement_all", "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  Q2 + Q6 concordance rows: {len(q2q6)}")

    # =============================================================================
    # Phase E: Repeated-measures robustness
    # =============================================================================

    print("\n[E] Repeated-measures (vine random intercept) robustness")
    rm_rows = []
    for (st, rs, net), grp_r in rumb_q2.groupby(["stage", "rootstock", "network"]):
        ours_grp = q2_rm[q2_rm["stage"] == st][["loc_id", "log2fc_ours", "padj_ours"]]
        merged_st = grp_r.merge(ours_grp, on="loc_id", how="inner")
        m = concordance(merged_st["rumb_mean"], merged_st["log2fc_ours"], merged_st["padj_ours"])
        m.update({"network": net, "rootstock": rs, "contrast": f"Q2_RM_{st}", "n_genes": len(merged_st)})
        rm_rows.append(m)
    for (net, rs), grp_r in rumb_collapsed.groupby(["network", "rootstock"]):
        merged_q6 = grp_r.merge(q6_rm[["loc_id", "log2fc_ours", "padj_ours"]],
                                on="loc_id", how="inner")
        m = concordance(merged_q6["rumb_mean_all"], merged_q6["log2fc_ours"], merged_q6["padj_ours"])
        m.update({"network": net, "rootstock": rs, "contrast": "Q6_pooled_RM", "n_genes": len(merged_q6)})
        rm_rows.append(m)
    rm = pd.DataFrame(rm_rows)[
        ["network", "rootstock", "contrast", "n_genes", "n", "spearman_rho",
         "pearson_r", "sign_agreement_all", "n_sig_padj", "sign_agreement_sig"]
    ]
    print(f"  RM concordance rows: {len(rm)}")

    # =============================================================================
    # Phase F: Single-gene trajectories (DCL2, calcium-binding protein)
    # =============================================================================

    print("\n[F] Single-gene trajectories")

    DCL2_VIT = "VIT_04s0023g00920"
    CABP_VIT = "VIT_14s0006g01400"

    def trajectory_table(vit_target):
        rumb = rumbaugh_long[rumbaugh_long["vit_id"] == vit_target]
        if rumb.empty:
            print(f"  WARNING: VIT {vit_target} not present in Rumbaugh tables")
            return pd.DataFrame()
        loc = rumb["loc_id"].iloc[0]
        ours = q2x[q2x["loc_id"] == loc][["our_year", "stage", "log2fc_ours", "padj_ours"]]
        if ours.empty:
            print(f"  WARNING: LOC {loc} not present in Q2x DESeq2 (no expressed-gene match)")
        out = rumb[["network", "vit_id", "loc_id", "description", "pathway",
                    "rootstock", "rumbaugh_year", "our_year", "stage",
                    "log2fc_rumbaugh"]].merge(
            ours.rename(columns={"log2fc_ours": "log2fc_ours_Q2x",
                                 "padj_ours":   "padj_ours_Q2x"}),
            on=["our_year", "stage"], how="left")
        return out

    dcl2_traj = trajectory_table(DCL2_VIT)
    cabp_traj = trajectory_table(CABP_VIT)
    print(f"  DCL2 trajectory rows: {len(dcl2_traj)}")
    print(f"  CaBP trajectory rows: {len(cabp_traj)}")

    # =============================================================================
    # Phase G: Write outputs + manifest
    # =============================================================================

    print("\n[G] Writing outputs to repo root")

    outputs = {
        "rumbaugh2022_network_members.csv":       wide,
        "rumbaugh2022_concordance_per_cell.csv":  per_cell,
        "rumbaugh2022_concordance_summary.csv":   summary,
        "rumbaugh2022_q1b_concordance.csv":       q1b_concord,
        "rumbaugh2022_q2q6_concordance.csv":      q2q6,
        "rumbaugh2022_rm_robustness.csv":         rm,
        "rumbaugh2022_dcl2_trajectory.csv":       dcl2_traj,
        "rumbaugh2022_cabp_trajectory.csv":       cabp_traj,
    }

    for fname, df in outputs.items():
        path = OUTPUT_DIR / fname
        df.to_csv(path, index=False)
        print(f"  wrote {fname}  ({df.shape[0]} rows, {df.shape[1]} cols)")
        manifest_rows.append({
            "stage": "output",
            "label": fname,
            "path": str(path),
            "n_unique_genes_or_rows": df.shape[0],
            "n_total_rows": df.shape[1],
        })

    # Reconciliation per network
    print("\n[G] Reconciliation per network:")
    for net in sorted(rumbaugh_long["network"].unique()):
        net_df  = rumbaugh_long[rumbaugh_long["network"] == net]
        n_total = net_df["loc_id"].nunique()
        n_exp   = net_df[net_df["loc_id"].isin(expressed_set)]["loc_id"].nunique()
        n_dropped = n_total - n_exp
        print(f"  {net}: {n_total} unique LOC IDs, "
              f"{n_exp} in expressed background, {n_dropped} not expressed")
        manifest_rows.append({
            "stage":                   "reconcile",
            "label":                   f"network_{net}",
            "path":                    "",
            "n_unique_genes_or_rows":  n_total,
            "n_total_rows":            n_exp,
        })

    manifest = pd.DataFrame(manifest_rows)
    mpath = OUTPUT_DIR / "rumbaugh2022_run_manifest.csv"
    manifest.to_csv(mpath, index=False)
    print(f"\n  wrote rumbaugh2022_run_manifest.csv ({len(manifest)} rows)")

    # =============================================================================
    # Headline summary printed to stdout
    # =============================================================================

    print("\n" + "=" * 78)
    print("HEADLINE: Q2x per-cell concordance (flattened across cells)")
    print("=" * 78)
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\n" + "=" * 78)
    print("HEADLINE: Q1b per-cohort concordance (Rumbaugh harvest-mean used)")
    print("=" * 78)
    pivot_q1b = q1b_concord.pivot_table(
        index=["network", "rootstock"], columns="cohort", values="rho_harvest", aggfunc="first")
    print(pivot_q1b.to_string(float_format=lambda v: f"{v:.3f}"))

    print("\n" + "=" * 78)
    print("HEADLINE: Q2 stage + Q6 pooled concordance")
    print("=" * 78)
    pivot_q2q6 = q2q6.pivot_table(
        index=["network", "rootstock"], columns="contrast",
        values="spearman_rho", aggfunc="first")
    print(pivot_q2q6.to_string(float_format=lambda v: f"{v:.3f}"))

    print("\n" + "=" * 78)
    print("HEADLINE: RM (vine random intercept) concordance vs DESeq2 rollup")
    print("=" * 78)
    pivot_rm = rm.pivot_table(
        index=["network", "rootstock"], columns="contrast",
        values="spearman_rho", aggfunc="first")
    print(pivot_rm.to_string(float_format=lambda v: f"{v:.3f}"))

    print("\nDone.")


if __name__ == "__main__":
    main()
