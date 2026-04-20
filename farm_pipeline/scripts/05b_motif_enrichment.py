#!/usr/bin/env python
"""05b. Motif enrichment per WGCNA module using AME (MEME suite).

For each module (a sheet in Module_genes.xlsx): AME scans the module
promoters (foreground) against non-module promoters (background) and
reports per-motif enrichment directly — far faster than FIMO scans.
"""
import argparse, os, subprocess, tempfile, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--promoters", required=True)
ap.add_argument("--motifs",    required=True)
ap.add_argument("--modules",   required=True)
ap.add_argument("--out",       required=True)
ap.add_argument("--qvalue",    type=float, default=0.05)
ap.add_argument("--threads",   type=int, default=4)
args = ap.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

# Load promoters into a dict
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

def write_fa(genes, path):
    with open(path, "w") as f:
        for g in genes:
            if g in prom:
                f.write(f">{g}\n{prom[g]}\n")

def run_ame(fg_path, bg_path, tmp):
    # --text streams TSV to stdout, bypasses HTML template bug in conda build
    cmd = ["ame", "--text",
           "--control", bg_path,
           "--evalue-report-threshold", "10",
           fg_path, args.motifs]
    proc = subprocess.run(cmd, check=False, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, text=True)
    rows, header = [], None
    for line in proc.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if header is None:
            header = parts; continue
        rows.append(dict(zip(header, parts)))
    return pd.DataFrame(rows)

import random
random.seed(42)
BG_SAMPLE_SIZE = 3000

xlsx = pd.ExcelFile(args.modules)
all_rows = []
for sheet in xlsx.sheet_names:
    df = pd.read_excel(xlsx, sheet_name=sheet)
    module_genes = set(df.iloc[:, 0].dropna().astype(str)) & all_genes
    bg_pool = all_genes - module_genes
    if len(module_genes) < 10 or len(bg_pool) < 50:
        continue
    # Random-subsample background to keep AME runtime bounded
    bg_genes = set(random.sample(sorted(bg_pool), min(BG_SAMPLE_SIZE, len(bg_pool))))
    with tempfile.TemporaryDirectory() as tmp:
        fg = os.path.join(tmp, "fg.fa"); bg = os.path.join(tmp, "bg.fa")
        write_fa(module_genes, fg); write_fa(bg_genes, bg)
        out_dir = os.path.join(tmp, "ame_out")
        res = run_ame(fg, bg, out_dir)
    if len(res):
        res.insert(0, "module", sheet)
        res["module_size"] = len(module_genes)
        all_rows.append(res)
    print(f"[05b] module={sheet} genes={len(module_genes)} enriched_motifs={len(res)}",
          file=sys.stderr)

if all_rows:
    out = pd.concat(all_rows, ignore_index=True)
    # Keep useful cols if AME produced them
    keep = [c for c in ["module","module_size","motif_ID","motif_alt_ID","motif_DB",
                         "rank","E-value","adj_p-value","consensus","TP","FP",
                         "pos","neg","tp_percent","fp_percent"] if c in out.columns]
    out = out[keep] if keep else out
else:
    out = pd.DataFrame()
out.to_csv(args.out, index=False)
print(f"[05b] wrote {len(out)} rows", file=sys.stderr)
