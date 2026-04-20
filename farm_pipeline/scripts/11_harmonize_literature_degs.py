#!/usr/bin/env python
"""11. Harmonize published Vitis DEG tables to NCBI LOC_ IDs.

Ingests user-supplied literature DEG CSVs from `data/literature/` (each
representing one comparison: stress vs control for a given pathogen,
cultivar, tissue). Maps heterogeneous identifiers (VIT_*, Vitvi*, LOC*,
gene symbols) to a common NCBI LOC_ ID space using a crosswalk table
(`resources/id_crosswalk.tsv`), and writes per-dataset harmonized
tables plus a coverage report.

Input schema for each literature CSV (required columns):
    gene_id      : identifier as published
    log2FC       : log2 fold change (stress vs control)
    pvalue       : raw p-value
    padj         : FDR-adjusted p-value
Optional:
    id_type      : VIT / Vvi / LOC / Symbol (inferred from prefix if missing)
    baseMean     : mean expression
    tested       : 1 if included in DE test universe, else 0 (defaults to 1)

YAML sidecar accepted next to each CSV (e.g. dataset.yaml) for metadata:
    stress: "GLRaV-3"   # GLRaV-3 | GPGV | PowderyMildew | DownyMildew | Heat | ...
    cultivar: "Cabernet Sauvignon"
    tissue: "leaf"
    reference: "Author et al. Year, DOI"
    deg_threshold_padj: 0.05
    deg_threshold_lfc: 1
"""
from __future__ import annotations
import argparse, os, sys, glob, re, json
from pathlib import Path
import pandas as pd
import yaml

ap = argparse.ArgumentParser()
ap.add_argument("--literature_dir", required=True)
ap.add_argument("--crosswalk",      required=True)
ap.add_argument("--out_dir",        required=True)
ap.add_argument("--out_coverage",   required=True)
args = ap.parse_args()

os.makedirs(args.out_dir, exist_ok=True)
os.makedirs(os.path.dirname(args.out_coverage), exist_ok=True)

REQUIRED_COLS = ["gene_id", "log2FC", "pvalue", "padj"]


def infer_id_type(gid: str) -> str:
    g = str(gid).strip()
    if g.startswith("VIT_") or g.startswith("VIT"):
        return "VIT"
    if g.startswith("Vitvi") or g.startswith("Vvi"):
        return "Vvi"
    if g.startswith("LOC"):
        return "LOC"
    return "Symbol"


def build_lookup(crosswalk_path: str) -> dict[str, dict[str, str]]:
    """Return {'VIT': {vit_id: loc_id}, 'Vvi': {...}, 'Symbol': {...}, 'LOC': {...}}."""
    cw = pd.read_csv(crosswalk_path, sep="\t")
    # Normalise column names
    cw.columns = [c.strip() for c in cw.columns]
    required = {"VIT_id", "Vvi_id", "LOC_id", "Symbol"}
    missing = required - set(cw.columns)
    if missing:
        sys.exit(f"[11] crosswalk missing columns: {missing}. Expected tab-separated "
                 f"file with VIT_id, Vvi_id, LOC_id, Symbol.")
    lookups: dict[str, dict[str, str]] = {"VIT": {}, "Vvi": {}, "Symbol": {}, "LOC": {}}
    for _, row in cw.iterrows():
        loc = str(row["LOC_id"]).strip()
        if not loc or loc == "nan":
            continue
        lookups["LOC"][loc] = loc
        for col, key in [("VIT_id", "VIT"), ("Vvi_id", "Vvi"), ("Symbol", "Symbol")]:
            v = row[col]
            if pd.isna(v):
                continue
            v = str(v).strip()
            if v:
                lookups[key][v] = loc
    return lookups


def load_dataset(csv_path: Path) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(csv_path)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{csv_path}: missing columns {missing}. Required: {REQUIRED_COLS}")
    # Coerce numerics
    for col in ("log2FC", "pvalue", "padj"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if "tested" not in df.columns:
        df["tested"] = 1
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df = df.dropna(subset=["gene_id"])

    meta_path = csv_path.with_suffix(".yaml")
    meta = {}
    if meta_path.exists():
        with open(meta_path) as f:
            meta = yaml.safe_load(f) or {}
    meta.setdefault("stress", csv_path.stem)
    meta.setdefault("cultivar", "unknown")
    meta.setdefault("tissue", "unknown")
    meta.setdefault("reference", "unknown")
    return df, meta


def harmonize(df: pd.DataFrame, lookups: dict[str, dict[str, str]]) -> pd.DataFrame:
    if "id_type" not in df.columns:
        df["id_type"] = df["gene_id"].map(infer_id_type)
    df["LOC_id"] = [
        lookups.get(t, {}).get(g)
        for g, t in zip(df["gene_id"], df["id_type"])
    ]
    return df


def main() -> None:
    lit_dir = Path(args.literature_dir)
    csvs = sorted(lit_dir.glob("*.csv"))
    if not csvs:
        print(f"[11] no literature CSVs in {lit_dir} - writing empty coverage report. "
              "Populate with published DEG tables (see docs/LITERATURE_SCHEMA.md) and rerun.",
              file=sys.stderr)
        pd.DataFrame(columns=["dataset","stress","cultivar","tissue","reference",
                              "n_total","n_mapped","coverage_pct",
                              "n_mapped_sig_padj05","harmonized_file"]
                     ).to_csv(args.out_coverage, index=False)
        return

    lookups = build_lookup(args.crosswalk)
    print(f"[11] crosswalk loaded: {len(lookups['LOC'])} LOC ids, "
          f"{len(lookups['VIT'])} VIT, {len(lookups['Vvi'])} Vvi, "
          f"{len(lookups['Symbol'])} symbols", file=sys.stderr)

    coverage_rows = []
    for csv in csvs:
        try:
            df, meta = load_dataset(csv)
        except ValueError as e:
            print(f"[11] SKIP {csv.name}: {e}", file=sys.stderr)
            continue
        df = harmonize(df, lookups)
        total = len(df)
        mapped = df.dropna(subset=["LOC_id"])
        # Write per-dataset harmonized output
        out_path = Path(args.out_dir) / f"{csv.stem}.harmonized.csv"
        out_cols = ["LOC_id", "gene_id", "id_type", "log2FC", "pvalue", "padj", "tested"]
        mapped[out_cols].to_csv(out_path, index=False)
        # Coverage row
        coverage_rows.append({
            "dataset": csv.stem,
            "stress": meta["stress"],
            "cultivar": meta["cultivar"],
            "tissue": meta["tissue"],
            "reference": meta["reference"],
            "n_total": total,
            "n_mapped": len(mapped),
            "coverage_pct": round(100 * len(mapped) / total, 2) if total else 0.0,
            "n_mapped_sig_padj05": int((mapped["padj"] < 0.05).sum()),
            "harmonized_file": str(out_path.relative_to(out_path.parents[2])),
        })
        print(f"[11] {csv.stem}: {len(mapped)}/{total} mapped "
              f"({100*len(mapped)/max(total,1):.1f}%), "
              f"{(mapped['padj']<0.05).sum()} DEGs at padj<0.05", file=sys.stderr)

    cov = pd.DataFrame(coverage_rows).sort_values("dataset")
    cov.to_csv(args.out_coverage, index=False)
    print(f"[11] wrote coverage report: {args.out_coverage}", file=sys.stderr)


if __name__ == "__main__":
    main()
