#!/usr/bin/env python
"""Fix #10: permutation null for biomarker-QTL co-localisation.

Uses the NCBI V. vinifera GFF (genome.gff) to build a genome-wide universe
of ~45,000 positioned genes. Biomarkers and QTL intervals come from the
main paper's GRBV_biomarker_QTL_crossref.csv. The permutation draws random
gene sets of the same size from the same chromosomes and counts how many
would co-localise within a 1 Mb window.
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--crossref", default="/Users/prempratapsingh/Data_Analysis/GRBV_biomarker_QTL_crossref.csv")
ap.add_argument("--gff",      default="/Users/prempratapsingh/Data_Analysis/farm_pipeline/resources/genome/genome.gff")
ap.add_argument("--window_mb", type=float, default=1.0)
ap.add_argument("--n_perm",    type=int,   default=5000)
ap.add_argument("--out", default="/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/fix10_qtl/qtl_permutation_null.csv")
args = ap.parse_args()
Path(args.out).parent.mkdir(parents=True, exist_ok=True)

cx = pd.read_csv(args.crossref)
bm = cx[["gene_id","gene_chr","gene_pos_Mb"]].dropna().drop_duplicates("gene_id").reset_index(drop=True)
bm.columns = ["gene_id","chromosome","start_Mb"]
bm["chromosome"] = bm["chromosome"].astype(str)
print(f"[fix10] {len(bm)} unique biomarkers with positions", file=sys.stderr)

qtl = cx[["qtl_locus","gene_chr","qtl_start_Mb","qtl_end_Mb"]].dropna().drop_duplicates(["qtl_locus","gene_chr"])
qtl["gene_chr"] = qtl["gene_chr"].astype(str)
print(f"[fix10] {qtl['qtl_locus'].nunique()} unique QTL intervals on {qtl['gene_chr'].nunique()} chromosomes", file=sys.stderr)

# Build universe from GFF
if not Path(args.gff).exists():
    sys.exit(f"[fix10] GFF missing: {args.gff}. Pull from FARM (rsync farm:~/grbv_extended/resources/genome/).")

# Map NCBI chromosome accessions (NC_*.*) to integer chromosomes used in the paper
univ_rows = []
chrom_pat = re.compile(r"gene=([^;]+)")
# First pass: identify chromosome numbering. The main paper uses integers 1..19.
# NCBI chromosome sequence IDs for V. vinifera GCF_030704535.1 follow NC_081832..NC_081850
# Parse chromosome accessions from ##sequence-region lines and map them in order.
chrom_ids = []
with open(args.gff) as f:
    for ln in f:
        if ln.startswith("##sequence-region"):
            parts = ln.split()
            if len(parts) >= 2 and parts[1].startswith("NC_"):
                chrom_ids.append(parts[1])
# Heuristic: first N NC_* accessions are chromosomes 1..N
chrom_map = {cid: str(i + 1) for i, cid in enumerate(chrom_ids[:19])}
print(f"[fix10] chromosome map: {len(chrom_map)} NC_ -> 1..N", file=sys.stderr)

with open(args.gff) as f:
    for ln in f:
        if ln.startswith("#"): continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene": continue
        chrom = chrom_map.get(p[0])
        if chrom is None: continue
        m = chrom_pat.search(p[8])
        if not m: continue
        name = m.group(1)
        start_mb = float(p[3]) / 1e6
        univ_rows.append((name, chrom, start_mb))
universe = pd.DataFrame(univ_rows, columns=["gene_id","chromosome","start_Mb"]).drop_duplicates("gene_id")
print(f"[fix10] universe: {len(universe)} positioned genes across chromosomes {sorted(universe['chromosome'].unique(), key=lambda x: int(x))}", file=sys.stderr)

bm_chr_counts = bm["chromosome"].value_counts().to_dict()

def overlaps(sample_df):
    n = 0
    for _, r in sample_df.iterrows():
        chrom = str(r["chromosome"])
        qs = qtl[qtl["gene_chr"] == chrom]
        if qs.empty: continue
        if ((r["start_Mb"] >= qs["qtl_start_Mb"] - args.window_mb) &
            (r["start_Mb"] <= qs["qtl_end_Mb"]   + args.window_mb)).any():
            n += 1
    return n

obs = overlaps(bm)
print(f"[fix10] observed: {obs}/{len(bm)} ({100*obs/len(bm):.1f}%) biomarkers within {args.window_mb} Mb of a QTL", file=sys.stderr)

rng = np.random.RandomState(42)
null = np.zeros(args.n_perm, dtype=int)
for i in range(args.n_perm):
    draws = []
    for chrom, k in bm_chr_counts.items():
        pool = universe[universe["chromosome"] == str(chrom)]
        if len(pool) < k: k = len(pool)
        if k == 0: continue
        draws.append(pool.sample(n=k, random_state=rng.randint(1, 2**31 - 1)))
    s = pd.concat(draws) if draws else universe.head(0)
    null[i] = overlaps(s)

obs_frac = obs / len(bm)
null_frac = null / max(len(bm), 1)
p_emp = (null >= obs).mean()
fold_enrich = obs / max(null.mean(), 0.5)

print(f"[fix10] null mean = {null.mean():.2f}  sd = {null.std():.2f}", file=sys.stderr)
print(f"[fix10] observed fraction {obs_frac:.3f} vs null mean {null_frac.mean():.3f}", file=sys.stderr)
print(f"[fix10] empirical one-sided p = {p_emp:.4f}  |  fold enrichment = {fold_enrich:.2f}x", file=sys.stderr)

pd.DataFrame([{
    "n_biomarkers_with_position": len(bm),
    "n_universe_genes": len(universe),
    "n_qtl_intervals": qtl["qtl_locus"].nunique(),
    "window_mb": args.window_mb,
    "n_permutations": args.n_perm,
    "observed_overlaps": int(obs),
    "observed_fraction": obs_frac,
    "null_mean_overlaps": float(null.mean()),
    "null_mean_fraction": float(null_frac.mean()),
    "null_sd": float(null.std()),
    "null_97_5_pct": float(np.percentile(null, 97.5)),
    "null_2_5_pct":  float(np.percentile(null, 2.5)),
    "empirical_p_one_sided": float(p_emp),
    "fold_enrichment": fold_enrich,
}]).to_csv(args.out, index=False)
print(f"[fix10] wrote {args.out}", file=sys.stderr)
