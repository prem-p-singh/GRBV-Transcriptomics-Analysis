#!/usr/bin/env python
"""05b. Motif enrichment per WGCNA module using FIMO.

For each module (Module_genes.xlsx sheet): subset promoter FASTA,
run FIMO vs JASPAR plants, count motif hits, compute Fisher enrichment
against the full-promoter background.
"""
import argparse, os, subprocess, tempfile, sys
import pandas as pd
from collections import defaultdict
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

ap = argparse.ArgumentParser()
ap.add_argument("--promoters", required=True)
ap.add_argument("--motifs",    required=True)
ap.add_argument("--modules",   required=True)
ap.add_argument("--out",       required=True)
ap.add_argument("--qvalue",    type=float, default=0.05)
ap.add_argument("--threads",   type=int, default=4)
args = ap.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

# Read all promoters
prom = {}
with open(args.promoters) as f:
    name = None
    for line in f:
        if line.startswith(">"):
            name = line[1:].strip().split()[0]; prom[name] = []
        else:
            prom[name].append(line.strip())
prom = {k: "".join(v) for k, v in prom.items()}
all_genes = set(prom)
print(f"[05b] loaded {len(prom)} promoters", file=sys.stderr)

def run_fimo(fa_path):
    """Return set of (motif_id, gene) hits passing q-value filter."""
    with tempfile.TemporaryDirectory() as tmp:
        cmd = ["fimo", "--thresh", "1e-4", "--max-stored-scores", "1000000",
               "--oc", tmp, args.motifs, fa_path]
        subprocess.run(cmd, check=True)
        hits = set()
        fimo_tsv = os.path.join(tmp, "fimo.tsv")
        if not os.path.exists(fimo_tsv):
            return hits
        with open(fimo_tsv) as f:
            header = None
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                if header is None:
                    header = line.rstrip("\n").split("\t"); continue
                fields = dict(zip(header, line.rstrip("\n").split("\t")))
                try:
                    q = float(fields.get("q-value", "1"))
                except ValueError:
                    continue
                if q <= args.qvalue:
                    hits.add((fields["motif_id"], fields["sequence_name"]))
        return hits

# Global FIMO scan (background)
print("[05b] running FIMO on full promoter set (background) …", file=sys.stderr)
bg_hits = run_fimo(args.promoters)
bg_by_motif = defaultdict(set)
for m, g in bg_hits:
    bg_by_motif[m].add(g)
N = len(all_genes)

# Per module: subset FASTA, rerun FIMO, Fisher-test each motif
xlsx = pd.ExcelFile(args.modules)
rows = []
for sheet in xlsx.sheet_names:
    df = pd.read_excel(xlsx, sheet_name=sheet)
    # assume first column is gene ID
    genes = set(df.iloc[:, 0].dropna().astype(str).tolist()) & all_genes
    if len(genes) < 10:
        continue
    K = len(genes)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tf:
        for g in genes:
            tf.write(f">{g}\n{prom[g]}\n")
        sub_fa = tf.name
    try:
        mod_hits = run_fimo(sub_fa)
    finally:
        os.unlink(sub_fa)
    mod_by_motif = defaultdict(set)
    for m, g in mod_hits:
        mod_by_motif[m].add(g)
    for motif, bg_g in bg_by_motif.items():
        a = len(mod_by_motif.get(motif, set()) & genes)
        b = K - a
        c = len(bg_g) - a
        d = (N - K) - c
        if a < 3 or c < 0 or d < 0:
            continue
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"module": sheet, "motif": motif,
                     "module_hits": a, "module_size": K,
                     "bg_hits": len(bg_g), "bg_size": N, "p": p})
    print(f"[05b] module={sheet} genes={K} motifs_scored={len(bg_by_motif)}", file=sys.stderr)

out = pd.DataFrame(rows)
if len(out):
    out["q"] = multipletests(out["p"], method="fdr_bh")[1]
    out = out.sort_values(["module", "q"])
out.to_csv(args.out, index=False)
print(f"[05b] wrote {len(out)} rows", file=sys.stderr)
