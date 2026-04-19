#!/usr/bin/env python
"""09. Defense gene repertoire profiling.

Categorises V. vinifera genes into defense classes via GFF `product` / `description`
keywords, then summarises expression and differential regulation (RB+ vs RB-) per class.
"""
import argparse, os, re, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True); ap.add_argument("--meta", required=True)
ap.add_argument("--gff",    required=True)
ap.add_argument("--out_summary", required=True); ap.add_argument("--out_heatmap", required=True)
ap.add_argument("--out_table",   required=True)
args = ap.parse_args()
os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)

# -------- Categorisation rules --------
RULES = [
    ("NLR (NBS-LRR)",      [r"\bNBS-LRR\b", r"disease resistance protein RPP", r"disease resistance protein RGA",
                            r"TMV resistance", r"CC-NBS-LRR", r"TIR-NBS-LRR"]),
    ("PRR (RLK/LRR-RLK)",  [r"receptor-like (kinase|protein)", r"LRR receptor", r"FLS2",
                            r"EFR", r"CERK"]),
    ("PR proteins",        [r"\bpathogenesis-related\b", r"\bPR-?\d+\b", r"thaumatin-like",
                            r"chitinase", r"\bbeta-1,3-glucanase\b", r"osmotin"]),
    ("WRKY TF",            [r"\bWRKY\b"]),
    ("SA biosynthesis/signal", [r"salicylic acid", r"isochorismate synthase", r"\bICS\b",
                                r"\bNPR1\b", r"\bTGA\b"]),
    ("JA biosynthesis/signal", [r"jasmon", r"lipoxygenase", r"allene oxide synthase",
                                r"\bOPR3\b", r"\bCOI1\b", r"\bMYC2\b", r"\bJAZ\b"]),
    ("ET biosynthesis/signal", [r"ethylene", r"\bACS\b", r"\bACO\b", r"\bEIN\b", r"\bETR\b",
                                r"\bERF\b", r"\bAP2\b"]),
    ("Callose / cell-wall defense", [r"callose", r"\bPMR\b", r"cellulose synthase-like"]),
    ("ROS / antioxidant",  [r"glutathione", r"peroxidase", r"catalase", r"superoxide dismutase",
                            r"thioredoxin"]),
    ("R-protein adapters", [r"\bRIN4\b", r"\bEDS1\b", r"\bPAD4\b", r"\bSGT1\b", r"\bRAR1\b"]),
]

# -------- Parse GFF for gene products --------
gene_desc = {}
with open(args.gff) as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] not in ("gene", "mRNA", "CDS"): continue
        attr = p[8]
        name = re.search(r'gene=([^;]+)', attr) or re.search(r'Name=([^;]+)', attr)
        prod = re.search(r'product=([^;]+)', attr)
        desc = re.search(r'description=([^;]+)', attr)
        if not name: continue
        gn = name.group(1)
        text = (prod.group(1) if prod else "") + " " + (desc.group(1) if desc else "")
        text = text.strip()
        if gn not in gene_desc or (len(text) > len(gene_desc[gn])):
            gene_desc[gn] = text

print(f"[09] parsed {len(gene_desc)} gene annotations", file=sys.stderr)

# Classify
def classify(text):
    cats = []
    for cat, pats in RULES:
        if any(re.search(p, text, re.IGNORECASE) for p in pats):
            cats.append(cat)
    return cats

classif = []
for gn, text in gene_desc.items():
    for c in classify(text):
        classif.append({"gene": gn, "category": c, "description": text[:160]})
classif = pd.DataFrame(classif)
classif.to_csv(args.out_table, index=False)

# -------- Expression summaries --------
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta   = pd.read_csv(args.meta,   sep="\t")
lib = counts.sum()
lcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)

rb_neg = meta.loc[meta["Infection_Group"] == "RB(-)", "Sample"].tolist()
rb_pos = meta.loc[meta["Infection_Group"] == "RB(+)", "Sample"].tolist()
rb_neg = [s for s in rb_neg if s in lcpm.columns]
rb_pos = [s for s in rb_pos if s in lcpm.columns]

rows = []
per_cat_mean = {}
for cat, grp in classif.groupby("category"):
    genes = [g for g in grp["gene"].unique() if g in lcpm.index]
    if len(genes) < 3:
        continue
    sub = lcpm.loc[genes]
    mean_neg = sub[rb_neg].mean(axis=1)
    mean_pos = sub[rb_pos].mean(axis=1)
    log2fc = mean_pos - mean_neg
    t = stats.ttest_ind(sub[rb_pos].values, sub[rb_neg].values, axis=1, equal_var=False)
    p = np.nanmedian(t.pvalue)
    rows.append({"category": cat, "n_genes": len(genes),
                 "mean_log2FC_pos_vs_neg": float(log2fc.mean()),
                 "median_gene_pvalue": float(p),
                 "frac_up (log2FC>0.5)": float((log2fc > 0.5).mean()),
                 "frac_down (log2FC<-0.5)": float((log2fc < -0.5).mean())})
    per_cat_mean[cat] = pd.concat([mean_neg.rename("RB-"),
                                   mean_pos.rename("RB+"),
                                   log2fc.rename("log2FC")], axis=1)
summary = pd.DataFrame(rows).sort_values("mean_log2FC_pos_vs_neg")
summary.to_csv(args.out_summary, index=False)

# Heatmap: per-category mean log2FC across maturity stages x season
combos = meta.groupby(["Sampling_Season","Maturity_stage"])["Sample"].apply(list).to_dict()
mat = []
labels = []
for cat, grp in classif.groupby("category"):
    genes = [g for g in grp["gene"].unique() if g in lcpm.index]
    if len(genes) < 3: continue
    vec = []
    cols = []
    for (season, stage), samples in combos.items():
        s_pos = [s for s in samples if s in rb_pos]
        s_neg = [s for s in samples if s in rb_neg]
        if len(s_pos) < 3 or len(s_neg) < 3:
            vec.append(np.nan); cols.append(f"{season}|{stage}"); continue
        vec.append((lcpm.loc[genes, s_pos].mean(axis=1) -
                    lcpm.loc[genes, s_neg].mean(axis=1)).mean())
        cols.append(f"{season}|{stage}")
    mat.append(vec); labels.append(cat)
heat = pd.DataFrame(mat, index=labels, columns=cols)
fig, ax = plt.subplots(figsize=(10, 6))
sns.heatmap(heat, center=0, cmap="RdBu_r", annot=True, fmt=".2f", ax=ax,
            cbar_kws={"label":"mean log2FC (RB+ − RB−)"})
ax.set_title("Defense category × condition"); fig.tight_layout()
fig.savefig(args.out_heatmap, dpi=200); plt.close(fig)
print("Done.")
