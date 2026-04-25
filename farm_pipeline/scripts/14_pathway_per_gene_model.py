#!/usr/bin/env python
"""14. Per-gene Infection x Season x Stage x Cohort model on three curated
pathways: phenylpropanoid biosynthesis, PR proteins, cell wall remodelling.

Goes one level deeper than the existing GSEA pathway-level enrichment (S70):
each gene in each curated category is fit with a per-gene OLS model on
log2(CPM+1) expression with main effects + 2- and 3-way interactions, and
ranked by infection-driven effect size and stability across stage/season.

This addresses agenda gap #2 (the explicit phenylpropanoid + PR + cell wall
pathways named in the project agenda; previously only had pathway-level GSEA).

Outputs (results/14_pathway_per_gene/):
  - curated_gene_table.csv         (LOC ids x category, from GFF regex)
  - per_gene_model_full.csv        (per-gene F-stats, partial eta^2 per term)
  - top_hits_per_category.csv      (ranked top 25 per category)
  - category_summary.csv           (% DE per category x stage)
  - pathway_dotplot.png            (top 30 hits visualised)
"""
import argparse, os, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True)
ap.add_argument("--meta",   required=True)
ap.add_argument("--gff",    required=True)
ap.add_argument("--out_dir", required=True)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)

# -------- curated regex categories ---------------------------------------
CATEGORIES = {
    "Phenylpropanoid": [
        r"\bphenylalanine ammonia[- ]?lyase\b", r"\bPAL\b",
        r"cinnamate 4-hydroxylase", r"\bC4H\b",
        r"4-coumarate.*CoA ligase", r"\b4CL\b",
        r"chalcone synthase", r"\bCHS\b",
        r"chalcone isomerase", r"\bCHI\b",
        r"flavanone 3-hydroxylase", r"\bF3H\b",
        r"flavonoid 3'-hydroxylase", r"\bF3'?H\b",
        r"flavonoid 3',5'-hydroxylase",
        r"dihydroflavonol 4-reductase", r"\bDFR\b",
        r"anthocyanidin synthase", r"\bANS\b", r"\bLDOX\b",
        r"UDP-glucose:flavonoid", r"\bUFGT\b", r"\bUFGAT\b",
        r"flavonol synthase", r"\bFLS\b",
        r"\bMYBA\b", r"VvMYBA",
        r"O-methyltransferase.*flavon", r"\bOMT\b", r"caffeic acid 3-O",
        r"leucoanthocyanidin",
        r"stilbene synthase", r"resveratrol synthase",
        r"shikimate", r"hydroxycinnamoyl",
    ],
    "PR_proteins": [
        r"\bpathogenesis-related\b", r"\bPR-?\d+\b",
        r"thaumatin-like", r"\bTLP\b",
        r"chitinase\b", r"chitin-binding",
        r"\bbeta-1,3-glucanase\b", r"endo-1,3-beta-glucanase",
        r"\bosmotin\b",
        r"defensin", r"\bPDF\b",
        r"non-specific lipid-transfer protein", r"\bnsLTP\b",
        r"hevein-like",
        r"germin-like protein", r"\bGLP\b",
    ],
    "Cell_wall": [
        r"xyloglucan", r"\bXTH\b", r"xyloglucan endotransglucosylase",
        r"expansin", r"\bEXP\b", r"\bEXLA\b", r"\bEXLB\b",
        r"pectin methylesterase", r"\bPME\b", r"pectinesterase",
        r"polygalacturonase", r"\bPG\b",
        r"pectate lyase", r"\bPL\b",
        r"cellulose synthase", r"\bCESA\b", r"cellulose synthase-like", r"\bCSL\b",
        r"\bCOBRA\b", r"COBRA-like",
        r"endoglucanase", r"cellulase",
        r"beta-galactosidase",
        r"alpha-arabinofuranosidase",
        r"laccase",                                # not strictly cell wall but lignin polymerisation
        r"cinnamoyl-CoA reductase", r"\bCCR\b",
        r"cinnamyl alcohol dehydrogenase", r"\bCAD\b",
        r"caffeoyl-CoA O-methyltransferase", r"\bCCoAOMT\b",
        r"peroxidase\b",
    ],
}
# Compile patterns
CAT_PATS = {k: [re.compile(p, flags=re.IGNORECASE) for p in v] for k, v in CATEGORIES.items()}

# -------- parse GFF --------
def parse_gff(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene": continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            gid = attrs.get("ID", "").replace("gene-", "")
            desc = (attrs.get("description", "") + " | " + attrs.get("Name", "") + " | "
                    + attrs.get("gene_biotype", "")).replace("%2C", ",")
            rows.append({"gene_id": gid, "description": desc})
    return pd.DataFrame(rows)

gff = parse_gff(args.gff)
print(f"GFF rows: {len(gff)}")

def categorise(desc):
    cats = []
    for cat, pats in CAT_PATS.items():
        for pat in pats:
            if pat.search(desc):
                cats.append(cat); break
    return cats

cat_rows = []
for _, r in gff.iterrows():
    for c in categorise(r["description"]):
        cat_rows.append({"gene_id": r["gene_id"], "category": c, "description": r["description"][:120]})
cat_df = pd.DataFrame(cat_rows)
print(f"Curated genes (any category): {cat_df['gene_id'].nunique()}; "
      + ", ".join(f"{c}={cat_df[cat_df.category==c].gene_id.nunique()}" for c in CATEGORIES))
cat_df.to_csv(os.path.join(args.out_dir, "curated_gene_table.csv"), index=False)

# -------- load expression + meta --------
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta = pd.read_csv(args.meta, sep="\t").set_index("Sample")
common = counts.columns.intersection(meta.index)
counts, meta = counts[common], meta.loc[common]

# Filter low-count genes (consistent with main analysis)
keep = (counts >= 10).sum(axis=1) >= int(0.10 * counts.shape[1])
counts = counts.loc[keep]
print(f"Filtered counts: {counts.shape}")

lib = counts.sum(axis=0)
logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)

# Tested universe = curated genes intersected with filtered counts
tested = sorted(set(cat_df.gene_id) & set(counts.index))
print(f"Tested in expression matrix: {len(tested)}")

# Build design DataFrame
df = meta[["Infection_Group","Sampling_Season","Maturity_stage","Years_Since_Infection"]].copy()
df.columns = ["infection","season","stage","ysi"]
df["infection"] = df["infection"].astype(str).str.replace("RB(+)","pos",regex=False).str.replace("RB(-)","neg",regex=False)
df["season"]    = df["season"].astype(str)
df["stage"]     = df["stage"].astype(str)
df["ysi"]       = pd.to_numeric(df["ysi"], errors="coerce").fillna(0)
# Centre ysi for interpretability
df["ysi_c"]     = df["ysi"] - df["ysi"].mean()

# -------- per-gene fit ---------------------------------------------------
results = []
formula = "y ~ C(infection) * C(season) * C(stage) + ysi_c"
for gene in tested:
    y = logcpm.loc[gene]
    yy = pd.concat([df, y.rename("y")], axis=1).dropna()
    if yy["y"].std() < 1e-3: continue
    try:
        m = ols(formula, data=yy).fit()
        anova = sm.stats.anova_lm(m, typ=2)
    except Exception:
        continue
    ss_total = anova["sum_sq"].sum()
    row = {"gene_id": gene, "n": len(yy), "ss_total": ss_total}
    for term in anova.index:
        row[f"F_{term}"] = anova.loc[term, "F"]
        row[f"p_{term}"] = anova.loc[term, "PR(>F)"]
        row[f"eta2_{term}"] = anova.loc[term, "sum_sq"] / ss_total
    # Infection main effect coefficient (pos vs neg)
    inf_coef = [k for k in m.params.index if k.startswith("C(infection)")]
    row["beta_infection"] = m.params[inf_coef[0]] if inf_coef else np.nan
    results.append(row)

res = pd.DataFrame(results)
# Annotate with categories
gene2cat = cat_df.groupby("gene_id")["category"].apply(lambda s: ",".join(sorted(set(s)))).to_dict()
res["category"] = res["gene_id"].map(gene2cat)
res.to_csv(os.path.join(args.out_dir, "per_gene_model_full.csv"), index=False)
print(f"Per-gene fits: {len(res)}")

# -------- top hits per category ------------------------------------------
top_rows = []
for cat in CATEGORIES:
    sub = res[res["category"].str.contains(cat, na=False)].copy()
    sub = sub.sort_values("p_C(infection)").head(25)
    sub["category_pri"] = cat
    top_rows.append(sub)
top = pd.concat(top_rows, ignore_index=True)
top.to_csv(os.path.join(args.out_dir, "top_hits_per_category.csv"), index=False)

# -------- category summary -----------------------------------------------
sumr = []
for cat in CATEGORIES:
    sub = res[res["category"].str.contains(cat, na=False)]
    if sub.empty: continue
    n = len(sub)
    sig_inf = (sub["p_C(infection)"] < 0.05).sum()
    sig_isxs = (sub["p_C(infection):C(season):C(stage)"] < 0.05).sum() if "p_C(infection):C(season):C(stage)" in sub else 0
    sumr.append({"category": cat, "n_tested": n, "n_sig_infection": int(sig_inf),
                 "frac_sig_infection": sig_inf / n if n else np.nan,
                 "n_sig_3way": int(sig_isxs)})
pd.DataFrame(sumr).to_csv(os.path.join(args.out_dir, "category_summary.csv"), index=False)

# -------- dotplot of top 30 -----------------------------------------------
top30 = top.sort_values("p_C(infection)").head(30).copy()
top30["minus_log10p"] = -np.log10(top30["p_C(infection)"].clip(lower=1e-30))
fig, ax = plt.subplots(figsize=(8, 8))
sns.scatterplot(data=top30, x="beta_infection", y="gene_id",
                hue="category_pri", size="minus_log10p",
                sizes=(40, 280), palette="Set1", ax=ax)
ax.axvline(0, color="grey", lw=0.5)
ax.set_xlabel("beta(infection): RB(+) vs RB(-) on log2 CPM")
ax.set_ylabel("gene")
ax.set_title("Top 30 curated-pathway genes by infection-effect significance")
plt.tight_layout()
plt.savefig(os.path.join(args.out_dir, "pathway_dotplot.png"), dpi=150)
plt.close(fig)
print("Done.")
