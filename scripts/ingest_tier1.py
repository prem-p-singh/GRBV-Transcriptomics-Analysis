"""
Ingest Tier 1 literature DEG tables into the schema expected by
farm_pipeline/scripts/11_harmonize_literature_degs.py.

Reads 6 heterogeneous supplementary files from
`farm_pipeline/data/literature/` and emits:

  - <stem>.csv  tidy (gene_id, log2FC, pvalue, padj [, baseMean, tested])
  - <stem>.yaml sidecar (stress, cultivar, tissue, reference, thresholds)

Naming convention (per LITERATURE_SCHEMA.md):
  <firstauthor><year>_<stress>_<cultivar>_<tissue>[_<stage>].csv

For papers with multiple contrasts (Blanco-Ulate, Song, Agudelo-Romero),
one tidy CSV is emitted per (cultivar, stage, tissue) slice so rule 12
can treat each as an independent comparison.

Run: python3 scripts/ingest_tier1.py
"""
from __future__ import annotations
from pathlib import Path
import re
import yaml
import numpy as np
import pandas as pd
from docx import Document

ROOT = Path("/Users/prempratapsingh/Data_Analysis")
LIT  = ROOT / "farm_pipeline/data/literature"
OUT  = LIT   # write CSV+YAML next to the sources


def write(name: str, df: pd.DataFrame, meta: dict) -> None:
    csv = OUT / f"{name}.csv"
    yml = OUT / f"{name}.yaml"
    # Coerce + reorder columns to schema expectations
    keep = ["gene_id", "log2FC", "pvalue", "padj"]
    for col in ("baseMean", "tested", "id_type"):
        if col in df.columns:
            keep.append(col)
    df_out = df[keep].copy()
    for col in ("log2FC", "pvalue", "padj", "baseMean"):
        if col in df_out.columns:
            df_out[col] = pd.to_numeric(df_out[col], errors="coerce")
    df_out = df_out.dropna(subset=["gene_id"])
    df_out["gene_id"] = df_out["gene_id"].astype(str).str.strip()
    df_out = df_out[df_out["gene_id"].ne("")].drop_duplicates("gene_id")
    df_out.to_csv(csv, index=False)
    with open(yml, "w") as f:
        yaml.safe_dump(meta, f, sort_keys=False)
    n_sig = int(((df_out["padj"] < 0.05) & (df_out["log2FC"].abs() >= 1)).sum()) \
            if df_out["padj"].notna().any() else int(len(df_out))
    print(f"  {name:64s}  {len(df_out):>6d} rows  (~{n_sig} sig @ padj<0.05|LFC|≥1)")


# =========================================================================
# 1. Blanco-Ulate et al. 2017 — GRBV × Zinfandel berries (erw506 Table S5)
# =========================================================================
def ingest_blanco_ulate():
    src = LIT / "erw506_suppl_Supplementary_Table_S5.xlsx"
    raw = pd.read_excel(src, skiprows=3)
    # Column 0 = VIT gene ID; columns 10..34 = 8 contrasts × (Base mean, Log2 FC, P value)
    # The paper states two vineyards (Healdsburg, Oakville) × 4 stages
    # (pre-véraison, véraison, post-véraison, harvest). Column order is H(4) O(4).
    vineyards = ["Healdsburg", "Oakville"]
    stages    = ["Pre-veraison", "Veraison", "Post-veraison", "Harvest"]
    raw["VIT"] = raw["Unnamed: 0"].astype(str).str.strip()
    # 8 contrasts: columns are (Base mean, Log2 FC GRBaV(+)/GRBaV(-), P value GRBaV(+)/GRBaV(-)) triplets
    fc_cols = [c for c in raw.columns if c.startswith("Log2 FC")]
    p_cols  = [c for c in raw.columns if c.startswith("P value")]
    bm_cols = [c for c in raw.columns if c.startswith("Base mean")]
    assert len(fc_cols) == 8, f"expected 8 FC cols, got {len(fc_cols)}"

    contrast_i = 0
    for vy in vineyards:
        for st in stages:
            df = pd.DataFrame({
                "gene_id":  raw["VIT"],
                "log2FC":   pd.to_numeric(raw[fc_cols[contrast_i]], errors="coerce"),
                "pvalue":   pd.to_numeric(raw[p_cols[contrast_i]],  errors="coerce"),
                "baseMean": pd.to_numeric(raw[bm_cols[contrast_i]], errors="coerce"),
            }).dropna(subset=["gene_id"])
            # Paper used q-value; recompute BH padj from pvalue for consistency
            from statsmodels.stats.multitest import multipletests
            mask = df["pvalue"].notna()
            padj = np.full(len(df), np.nan)
            if mask.sum():
                padj[mask.values] = multipletests(df.loc[mask, "pvalue"], method="fdr_bh")[1]
            df["padj"] = padj
            df["tested"] = 1
            df["id_type"] = "VIT"
            name = f"BlancoUlate2017_GRBV_Zinfandel_berry_{vy.lower()}_{st.lower()}"
            meta = dict(
                stress="GRBV",
                cultivar=f"Zinfandel ({vy})",
                tissue="berry",
                stage=st,
                reference="Blanco-Ulate et al. 2017, J Exp Bot 68(5):1225-1238, doi:10.1093/jxb/erw506",
                platform="RNA-Seq",
                deg_threshold_padj=0.05,
                deg_threshold_lfc=1.0,
                full_panel=True,
            )
            write(name, df, meta)
            contrast_i += 1


# =========================================================================
# 2. Song et al. 2022 — GLRaV-3 × Cabernet Franc leaves + berries (Table S4)
# =========================================================================
def ingest_song():
    src = LIT / "Song et al. RNA-Seq grapevine-GLRaV-3 STable S4 Identification of DEGs via RNA-Seq.xlsx"
    sheet_map = {
        "E-L31 Leaf DEGs":         ("leaf",  "E-L31"),
        "E-L35 Leaf DEGs":         ("leaf",  "E-L35"),
        "E-L38 Leaf DEGs (2018)":  ("leaf",  "E-L38"),
        "E-L31 Berry DEGs":        ("berry", "E-L31"),
        "E-L35 Berry DEGs":        ("berry", "E-L35"),
        "E-L38 Berry DEGs ":       ("berry", "E-L38"),
    }
    for sheet, (tissue, stage) in sheet_map.items():
        raw = pd.read_excel(src, sheet_name=sheet, skiprows=3)
        # Columns after skip: gene ID, baseMean, log2FoldChange, FDR, direction (may have trailing spaces)
        cols = {c.strip(): c for c in raw.columns}
        gid  = cols.get("gene ID")
        bm   = cols.get("baseMean")
        lfc  = cols.get("log2FoldChange")
        fdr  = cols.get("FDR")
        if not all((gid, lfc, fdr)):
            print(f"  [Song {sheet}] missing required column, skipping")
            continue
        df = pd.DataFrame({
            "gene_id": raw[gid].astype(str).str.strip(),
            "log2FC":  pd.to_numeric(raw[lfc], errors="coerce"),
            "padj":    pd.to_numeric(raw[fdr], errors="coerce"),
            "baseMean": pd.to_numeric(raw[bm], errors="coerce") if bm else np.nan,
        })
        df["pvalue"] = np.nan   # paper reports only adjusted FDR
        df["tested"] = 1
        df["id_type"] = "VIT"
        name = f"Song2022_GLRaV3_CabernetFranc_{tissue}_{stage.lower().replace('-', '')}"
        meta = dict(
            stress="GLRaV-3",
            cultivar="Cabernet Franc",
            tissue=tissue,
            stage=stage,
            reference="Song, Hanner, Meng 2022, Viruses 14(8):1831, doi:10.3390/v14081831",
            platform="RNA-Seq",
            deg_threshold_padj=0.05,
            deg_threshold_lfc=0.0,
            full_panel=False,
            note="Pre-filtered to FDR < 0.05 (paper threshold). No untested-gene universe in table.",
        )
        write(name, df, meta)


# =========================================================================
# 3. Agudelo-Romero et al. 2015 — Botrytis × Trincadeira berries
# =========================================================================
def ingest_agudelo_romero():
    src = LIT / "jexbot132548_file001.xlsx"
    raw = pd.read_excel(src, sheet_name="Table S4SignificantGenes2.0fold")
    gid = "12X V1 ID"
    contrasts = [
        ("EL33", "EL33Inf/EL33Mock significant FC", "adj. p-value Limma Infection in EL33"),
        ("EL35", "EL35Inf/EL35Mock significant FC", "adj. p-value Limma Infection in EL35"),
    ]
    for stage, fc_col, p_col in contrasts:
        df = pd.DataFrame({
            "gene_id": raw[gid].astype(str).str.strip(),
            "log2FC":  np.log2(pd.to_numeric(raw[fc_col], errors="coerce").abs().replace(0, np.nan))
                       * np.sign(pd.to_numeric(raw[fc_col], errors="coerce")),
            "padj":    pd.to_numeric(raw[p_col], errors="coerce"),
        })
        df["pvalue"] = np.nan
        df["tested"] = 1
        df["id_type"] = "VIT"
        # Drop rows with missing gene_id or both logFC and padj
        df = df[df["gene_id"].ne("") & df["gene_id"].ne("nan")]
        df = df[df["log2FC"].notna() | df["padj"].notna()]
        name = f"AgudeloRomero2015_Botrytis_Trincadeira_berry_{stage.lower()}"
        meta = dict(
            stress="Botrytis",
            cultivar="Trincadeira",
            tissue="berry",
            stage=stage,
            reference="Agudelo-Romero et al. 2015, J Exp Bot 66(7):1769-1785, doi:10.1093/jxb/eru517",
            platform="NimbleGen microarray",
            deg_threshold_padj=0.05,
            deg_threshold_lfc=1.0,
            full_panel=False,
            note="Pre-filtered to |FC| >= 2 + 5% FDR. FC column was raw fold-change; converted to log2(|FC|)*sign.",
        )
        write(name, df, meta)


# =========================================================================
# 4. Fung et al. 2008 — powdery mildew × Cabernet Sauvignon leaves
# =========================================================================
def ingest_fung():
    src = LIT / "108712Table_S1.xls"
    # Skip 3 title rows; layout has a "V. vinifera 'Cabernet sauvignon'" section divider row
    raw = pd.read_excel(src, sheet_name="s", engine="xlrd", skiprows=3)
    raw = raw.rename(columns={c: c.strip() for c in raw.columns})
    # Identify section: drop rows where Probeset ID is NaN or section-title text
    ps = raw["Probeset ID"].astype(str).str.strip()
    keep_mask = ps.str.match(r"^[0-9]+_at$")
    # Restrict to Cabernet Sauvignon block (the first contiguous block of probes)
    # The xls has CS first, then Norton — we split at the first blank line
    first_norton = ~keep_mask & ps.str.contains("Norton", na=False)
    if first_norton.any():
        split_i = first_norton.idxmax()
        raw = raw.loc[:split_i-1]
    raw = raw[keep_mask]
    # AT number column uses Arabidopsis homolog; keep as gene_id for cross-species mapping
    # Two P-value columns: 'P-value' for mock-vs-infected at each time, 'P-value.1' for post-hoc. Use 'P-value'.
    df = pd.DataFrame({
        "gene_id": raw["AT number"].astype(str).str.strip(),
        "log2FC":  np.nan,   # time-course with 6 points × 3 replicates; use 24h PM/mock mean ratio
        "pvalue":  pd.to_numeric(raw["P-value"], errors="coerce"),
    })
    # Compute log2FC at 24 hpi: mean across replicates in cols "24 h" (infected block)
    # Columns: 0h, 4h, 8h, 12h, 24h, 48h repeated 4 times (mock x 3 rep, infected x 3 rep — paper-specific layout).
    # Safely recover 24h PM / 24h mock ratio using columns indexed by name suffix.
    h24_cols = [c for c in raw.columns if re.match(r"^24 h(\.\d+)?$", str(c))]
    if len(h24_cols) >= 4:
        # First 3 are mock replicates, last 3 are PM replicates (convention from Fung 2008)
        # Layout confirmed: 0h,4h,..48h (mock rep1) then same (mock rep2) ... then (infected rep1)...
        # Simpler: take ratio of last-half mean / first-half mean.
        half = len(h24_cols) // 2
        mock = raw[h24_cols[:half]].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        inf  = raw[h24_cols[half:]].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        df["log2FC"] = np.log2((inf.replace(0, np.nan)) / (mock.replace(0, np.nan)))
    # BH padj
    from statsmodels.stats.multitest import multipletests
    mask = df["pvalue"].notna()
    padj = np.full(len(df), np.nan)
    if mask.sum():
        padj[mask.values] = multipletests(df.loc[mask, "pvalue"], method="fdr_bh")[1]
    df["padj"] = padj
    df["tested"] = 1
    df["id_type"] = "Symbol"  # AT numbers treated as symbols; won't map to LOC but preserved
    name = "Fung2008_PowderyMildew_CabernetSauvignon_leaf_24h"
    meta = dict(
        stress="PowderyMildew",
        cultivar="Cabernet Sauvignon",
        tissue="leaf",
        stage="24 hpi",
        reference="Fung et al. 2008, Plant Physiol 146(1):236-249, doi:10.1104/pp.107.108712",
        platform="Affymetrix Vitis GeneChip microarray",
        deg_threshold_padj=0.05,
        deg_threshold_lfc=1.0,
        full_panel=False,
        note=("626 PM-responsive probes (Cabernet Sauvignon only; Norton block dropped). gene_id = "
              "Arabidopsis homolog AT number. Direct VIT/LOC mapping not available; use Arabidopsis "
              "orthology tables for bridging in Tier 2."),
    )
    write(name, df, meta)


# =========================================================================
# 5. Polesani et al. 2010 — downy mildew × V. vinifera leaf disks (Addl 6)
# =========================================================================
def ingest_polesani():
    src = LIT / "12864_2009_2711_MOESM6_ESM.XLS"
    raw = pd.read_excel(src, sheet_name="modulated genes", engine="xlrd", skiprows=1)
    # Sheet has 4 side-by-side subtables. First subtable = "All modulated genes".
    # Columns 0..6 of raw (after skiprows=1) = Name, Target Description, Functional Category,
    #  Vr 12h fold change, Vr 24h fold change, Vv 12h fold change, Vv 24h fold change
    sub = raw.iloc[:, :7].copy()
    sub.columns = ["Name", "Target", "Category", "Vr12", "Vr24", "Vv12", "Vv24"]
    sub = sub.dropna(subset=["Name"])
    # Keep rows with a V. vinifera 24h fold change (our target contrast)
    sub["Vv24"] = pd.to_numeric(sub["Vv24"], errors="coerce")
    sub = sub[sub["Vv24"].notna()]
    df = pd.DataFrame({
        "gene_id": sub["Name"].astype(str).str.strip(),
        "log2FC":  np.log2(sub["Vv24"].abs().replace(0, np.nan)) * np.sign(sub["Vv24"]),
        "pvalue":  np.nan,
        "padj":    0.05,  # paper filter threshold; placeholder
    })
    df["tested"] = 1
    df["id_type"] = "Symbol"  # TIGR TC IDs -> no direct LOC mapping; keep as symbol
    name = "Polesani2010_DownyMildew_Vvinifera_leaf_24h"
    meta = dict(
        stress="DownyMildew",
        cultivar="Vitis vinifera (susceptible)",
        tissue="leaf",
        stage="24 hpi",
        reference="Polesani et al. 2010, BMC Genomics 11:117, doi:10.1186/1471-2164-11-117",
        platform="Combimatrix custom microarray",
        deg_threshold_padj=0.05,
        deg_threshold_lfc=1.0,
        full_panel=False,
        note=("gene_id = TIGR Tentative Consensus (TC#) ID. TC->VIT or TC->LOC mapping not available; "
              "will show 0% coverage under id_crosswalk.tsv. Sequence-based BLAST mapping needed for Tier 2."),
    )
    write(name, df, meta)


# =========================================================================
# 6. Liu et al. 2012 — heat × Cabernet Sauvignon leaves (Addl 7 DOCX)
# =========================================================================
def ingest_liu():
    src = LIT / "12870_2012_1124_MOESM7_ESM.docx"
    doc = Document(src)
    table = doc.tables[0]
    rows = []
    for r in table.rows[1:]:  # skip header
        cells = [c.text.strip() for c in r.cells]
        if len(cells) < 6:
            continue
        category, probe, accession, fc_heat, gene_short, gene_desc = cells[:6]
        if not probe or not re.match(r"^\d+_at$", probe):
            continue
        rows.append({
            "gene_id":   probe,
            "accession": accession,
            "gene_name": gene_short or gene_desc,
            "log2FC":    np.log2(abs(float(fc_heat))) * (1 if float(fc_heat) >= 0 else -1)
                          if fc_heat and fc_heat.replace("-","").replace(".","").isdigit() else np.nan,
        })
    df = pd.DataFrame(rows)
    df["pvalue"] = np.nan
    df["padj"]   = 0.05
    df["tested"] = 1
    df["id_type"] = "Symbol"
    name = "Liu2012_Heat_CabernetSauvignon_leaf"
    meta = dict(
        stress="Heat",
        cultivar="Cabernet Sauvignon",
        tissue="leaf",
        stage="heat stress",
        reference="Liu et al. 2012, BMC Plant Biol 12:174, doi:10.1186/1471-2229-12-174",
        platform="Affymetrix Grape Genome microarray",
        deg_threshold_padj=0.05,
        deg_threshold_lfc=1.0,
        full_panel=False,
        note=("gene_id = Affymetrix Vitis GeneChip probeset. Direct probe->LOC mapping unavailable; "
              "will show 0% coverage. Same chip as Fung 2008 -- one shared probe->VIT table would "
              "unlock both in Tier 2."),
    )
    write(name, df, meta)


# -------------------------------------------------------------------- main --
def main():
    print("[ingest_tier1] writing tidy CSVs + YAML sidecars into",
          LIT.relative_to(ROOT))
    print()
    print("Blanco-Ulate 2017 (GRBV, Zinfandel, berry x 8 contrasts):")
    ingest_blanco_ulate()
    print()
    print("Song 2022 (GLRaV-3, Cabernet Franc, leaf + berry x 3 stages):")
    ingest_song()
    print()
    print("Agudelo-Romero 2015 (Botrytis, Trincadeira, berry x 2 stages):")
    ingest_agudelo_romero()
    print()
    print("Fung 2008 (Powdery mildew, Cabernet Sauvignon, leaf):")
    ingest_fung()
    print()
    print("Polesani 2010 (Downy mildew, V. vinifera, leaf):")
    ingest_polesani()
    print()
    print("Liu 2012 (Heat, Cabernet Sauvignon, leaf):")
    ingest_liu()


if __name__ == "__main__":
    main()
