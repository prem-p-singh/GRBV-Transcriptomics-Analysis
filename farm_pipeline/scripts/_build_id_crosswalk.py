#!/usr/bin/env python
"""Build a VIT_ / Vvi_ / Symbol -> NCBI LOC_ ID crosswalk from the NCBI GFF.

Parses the downloaded V. vinifera reference GFF and extracts every gene
record's NCBI LOC symbol plus any Dbxref= entries (e.g. VGNC, IGGP IDs,
GeneID mappings). Many NCBI grape records include the old IGGP locus tag
as `old_locus_tag=VIT_XXsXXXXXg00000`; this script captures those.

For Vitvi* (12X.v2 / PN40024 annotation) and gene symbols, we rely on
either:
  1. a user-supplied supplementary `ortholog_map.tsv` from Ensembl
     BioMart with columns `VIT_id`, `Vvi_id`, `Symbol`, `LOC_id`, which
     this script will merge in, OR
  2. pattern-based heuristics (last resort).

Output: resources/id_crosswalk.tsv with columns
    VIT_id  Vvi_id  LOC_id  Symbol  source
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--gff",         required=True)
ap.add_argument("--ortholog_map", default=None,
                help="Optional user-supplied TSV with VIT_id/Vvi_id/LOC_id/Symbol columns")
ap.add_argument("--out",         required=True)
args = ap.parse_args()

rows = []
seen_loc = set()
n_with_oldlt = 0
print(f"[crosswalk] parsing {args.gff}", file=sys.stderr)

with open(args.gff) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        attr = p[8]
        name_m = re.search(r"gene=([^;]+)", attr) or re.search(r"Name=([^;]+)", attr)
        if not name_m:
            continue
        loc = name_m.group(1)
        seen_loc.add(loc)
        symbol = None
        if "gene_synonym=" in attr:
            symbol = re.search(r"gene_synonym=([^;]+)", attr).group(1).split(",")[0]
        # Old IGGP locus tag (often stored here)
        vit = None
        for mt in re.finditer(r"old_locus_tag=([^;]+)", attr):
            for tag in mt.group(1).split(","):
                if tag.startswith("VIT_") or tag.startswith("VIT"):
                    vit = tag.strip()
                    break
        if vit:
            n_with_oldlt += 1
        rows.append({
            "VIT_id": vit,
            "Vvi_id": None,
            "LOC_id": loc,
            "Symbol": symbol,
            "source": "ncbi_gff",
        })

cw = pd.DataFrame(rows)
print(f"[crosswalk] genes parsed: {len(cw)}; with VIT old_locus_tag: {n_with_oldlt}",
      file=sys.stderr)

# Merge user-supplied ortholog map if available
if args.ortholog_map and Path(args.ortholog_map).exists():
    print(f"[crosswalk] merging user ortholog map {args.ortholog_map}", file=sys.stderr)
    ext = pd.read_csv(args.ortholog_map, sep="\t")
    for col in ("VIT_id", "Vvi_id", "LOC_id", "Symbol"):
        if col not in ext.columns:
            ext[col] = None
    ext["source"] = "user_orthologs"
    # Prefer existing GFF rows; append new from ortholog map keyed on LOC_id
    combined = pd.concat([cw, ext[["VIT_id", "Vvi_id", "LOC_id", "Symbol", "source"]]],
                         ignore_index=True)
    # Collapse per LOC_id — take first non-null in each column
    def collapse(g):
        return pd.Series({
            "VIT_id":  g["VIT_id"].dropna().iloc[0] if g["VIT_id"].dropna().size else None,
            "Vvi_id":  g["Vvi_id"].dropna().iloc[0] if g["Vvi_id"].dropna().size else None,
            "LOC_id":  g.name,
            "Symbol":  g["Symbol"].dropna().iloc[0] if g["Symbol"].dropna().size else None,
            "source":  "+".join(sorted(set(g["source"].dropna().astype(str)))) or "unknown",
        })
    cw = (combined.dropna(subset=["LOC_id"])
                   .groupby("LOC_id", sort=False).apply(collapse, include_groups=False)
                   .reset_index(drop=True))

cw = cw[["VIT_id", "Vvi_id", "LOC_id", "Symbol", "source"]]
cw.to_csv(args.out, sep="\t", index=False)
print(f"[crosswalk] wrote {len(cw)} rows to {args.out}; "
      f"coverage: VIT {cw['VIT_id'].notna().sum()}, "
      f"Vvi {cw['Vvi_id'].notna().sum()}, "
      f"Symbol {cw['Symbol'].notna().sum()}", file=sys.stderr)
