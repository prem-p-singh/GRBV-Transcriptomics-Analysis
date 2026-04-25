#!/usr/bin/env python
"""15. Carbon translocation / phloem-loading module analysis.

Tests Martinez-Luscher 2018's hypothesis that GRBV blocks hexose
translocation rather than primary carbon assimilation.  Builds a curated
sugar-transport module from GFF (SWEETs, SUTs/SUC, STPs/HXT, CWINVs,
vacuolar invertases, sucrose synthases, hexokinase, plus the already-
flagged SPS1).

Per-gene per-season Brix correlation, module eigengene per sample, module-
trait correlations against (titer, season, stage, cohort).  Reports whether
the sugar-transport module recapitulates the biphasic titer-Brix pattern.

Outputs (results/15_carbon_translocation/):
  - sugar_module_genes.csv          (curated gene list with subclass)
  - sugar_module_eigengene.csv      (per-sample module score, restricted to RB+)
  - sugar_module_traits.csv         (module-trait correlations)
  - sugar_brix_per_gene.csv         (per-gene Brix Spearman, partial-rho controlling for season)
  - sugar_module_trajectory.png     (module score x stage x season)
"""
import argparse, os, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from openpyxl import load_workbook
from sklearn.decomposition import PCA

ap = argparse.ArgumentParser()
ap.add_argument("--counts",   required=True)
ap.add_argument("--meta",     required=True)
ap.add_argument("--gff",      required=True)
ap.add_argument("--titer",    required=True)
ap.add_argument("--biochem",  required=True)
ap.add_argument("--out_dir",  required=True)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)

# -------- curated regex --------------------------------------------------
SUBCLASS = {
    "SWEET":           [r"\bSWEET\b", r"sugars will eventually be exported"],
    "SUT_SUC":         [r"sucrose transporter", r"\bSUT\b", r"\bSUC\b", r"sucrose-proton symporter"],
    "STP_HXT":         [r"hexose transporter", r"\bSTP\b", r"\bHXT\b", r"sugar transport protein",
                        r"glucose transporter"],
    "Cell_wall_INV":   [r"cell wall.*invertase", r"\bCWINV\b", r"\bCWIN\b",
                        r"beta-fructofuranosidase", r"acid.*invertase"],
    "Vacuolar_INV":    [r"vacuolar.*invertase", r"\bVIN\b", r"\bVI\b.*invertase"],
    "Sucrose_synth":   [r"sucrose synthase", r"\bSuSy\b", r"\bSUS\b"],
    "SPS":             [r"sucrose-phosphate synthase", r"\bSPS\b"],
    "Sucrose_phosphat":[r"sucrose-phosphate phosphatase", r"\bSPP\b"],
    "Hexokinase":      [r"\bhexokinase\b", r"\bHXK\b"],
    "Trehalose":       [r"trehalose-6-phosphate", r"\bTPS\b", r"\bTPP\b",
                        r"trehalase"],
    "Phloem_marker":   [r"phloem", r"sieve element", r"companion cell"],
}
PATS = {k: [re.compile(p, flags=re.IGNORECASE) for p in v] for k, v in SUBCLASS.items()}

# -------- parse GFF ------------------------------------------------------
def parse_gff(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene": continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            gid = attrs.get("ID", "").replace("gene-", "")
            desc = (attrs.get("description", "") + " | " + attrs.get("Name", "")).replace("%2C", ",")
            rows.append({"gene_id": gid, "description": desc})
    return pd.DataFrame(rows)

gff = parse_gff(args.gff)
sub_rows = []
for _, r in gff.iterrows():
    for sub, pats in PATS.items():
        for p in pats:
            if p.search(r["description"]):
                sub_rows.append({"gene_id": r["gene_id"], "subclass": sub,
                                 "description": r["description"][:120]})
                break
sub_df = pd.DataFrame(sub_rows).drop_duplicates(subset=["gene_id"])
print(f"Sugar-transport curated genes: {len(sub_df)}")
print(sub_df.groupby("subclass").size())
sub_df.to_csv(os.path.join(args.out_dir, "sugar_module_genes.csv"), index=False)

# -------- expression -----------------------------------------------------
counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta = pd.read_csv(args.meta, sep="\t").set_index("Sample")
common = counts.columns.intersection(meta.index)
counts, meta = counts[common], meta.loc[common]
meta.index.name = "Sample"
keep = (counts >= 10).sum(axis=1) >= int(0.10 * counts.shape[1])
counts = counts.loc[keep]
lib = counts.sum(axis=0)
logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
module_genes = sorted(set(sub_df.gene_id) & set(counts.index))
print(f"Sugar-module genes in expression matrix: {len(module_genes)}")

# -------- module eigengene (PC1 of standardised module expression) -------
X = logcpm.loc[module_genes].T  # samples x genes
Xs = (X - X.mean()) / X.std(ddof=0).replace(0, 1)
pca = PCA(n_components=2).fit(Xs.fillna(0))
eig = pca.transform(Xs.fillna(0))
# Sign convention: positive = high-sugar-transport
if pca.components_[0].mean() < 0: eig[:, 0] *= -1
me = pd.DataFrame({"Sample": Xs.index, "ME_sugar": eig[:, 0],
                   "ME_sugar_pc2": eig[:, 1]})
meta_cols_needed = ["Infection_Group","Sampling_Season","Maturity_stage",
                    "Years_Since_Infection"]
if "biorep_id" in meta.columns: meta_cols_needed.append("biorep_id")
me = me.merge(meta[meta_cols_needed],
              left_on="Sample", right_index=True, how="left")
me["Sampling_Season"] = me["Sampling_Season"].astype(str)
me.to_csv(os.path.join(args.out_dir, "sugar_module_eigengene.csv"), index=False)
print(f"Eigengene var explained: PC1={pca.explained_variance_ratio_[0]:.2%}, "
      f"PC2={pca.explained_variance_ratio_[1]:.2%}")

# -------- module-trait correlations --------------------------------------
trait_rows = []
me_rb = me[me["Infection_Group"]=="RB(+)"].copy()
me_rb["stage_rank"] = me_rb["Maturity_stage"].map({"Pre_veraison":1,"Veraison":2,
                                                    "Post_veraison":3,"Harvest":4})
me_rb["season_num"] = (me_rb["Sampling_Season"]=="2022").astype(int)
for trait in ["stage_rank","season_num","Years_Since_Infection"]:
    s = me_rb[["ME_sugar", trait]].dropna()
    if len(s) < 10: continue
    rho, p = stats.spearmanr(s["ME_sugar"], s[trait])
    trait_rows.append({"trait": trait, "n": len(s), "rho": rho, "p": p})
pd.DataFrame(trait_rows).to_csv(os.path.join(args.out_dir, "sugar_module_traits.csv"), index=False)

# -------- per-gene Brix correlation (Harvest, RB+ + RB-) -----------------
# Use tidied biochemistry_harvest.csv (parser lives in repo-root
# biochemistry_harvest_prep.py).
bc = pd.read_csv(args.biochem)
bc = bc.rename(columns={"Biorep":"biorep_id"})
bc["Sampling_Season"] = bc["Sampling_Season"].astype(str)
bc["Years_Since_Infection"] = pd.to_numeric(bc["Infection_Year"], errors="coerce").astype("Int64")
bc["biorep_id"] = pd.to_numeric(bc["biorep_id"], errors="coerce").astype("Int64")

mh = meta[meta["Maturity_stage"]=="Harvest"].reset_index().copy()
mh["Sampling_Season"] = mh["Sampling_Season"].astype(str)
if "biorep_id" not in mh.columns:
    raise SystemExit("meta lacks biorep_id (need sample_metadata_v2.txt)")
mh["biorep_id"] = pd.to_numeric(mh["biorep_id"], errors="coerce").astype("Int64")
mh["Years_Since_Infection"] = pd.to_numeric(mh["Years_Since_Infection"], errors="coerce").astype("Int64")
mh = mh.merge(bc[["Sampling_Season","Years_Since_Infection","biorep_id","Infection_Group","Brix"]],
              on=["Sampling_Season","Years_Since_Infection","biorep_id","Infection_Group"], how="inner")
print(f"Harvest+biochem merged: {len(mh)}")

per_gene = []
for g in module_genes:
    expr = logcpm.loc[g, mh["Sample"]].values
    sub = pd.DataFrame({"expr": expr, "Brix": mh["Brix"].values,
                        "season": mh["Sampling_Season"].values}).dropna()
    if len(sub) < 8: continue
    rho_all, p_all = stats.spearmanr(sub["expr"], sub["Brix"])
    rows = {"gene_id": g, "n": len(sub),
            "rho_pooled": rho_all, "p_pooled": p_all}
    for s in ["2021","2022"]:
        ssub = sub[sub.season==s]
        if len(ssub) >= 5:
            r, p = stats.spearmanr(ssub["expr"], ssub["Brix"])
            rows[f"rho_{s}"] = r; rows[f"p_{s}"] = p
    per_gene.append(rows)
pg = pd.DataFrame(per_gene)
pg = pg.merge(sub_df[["gene_id","subclass","description"]], on="gene_id", how="left")
pg.to_csv(os.path.join(args.out_dir, "sugar_brix_per_gene.csv"), index=False)

# -------- trajectory plot -------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4))
order = ["Pre_veraison","Veraison","Post_veraison","Harvest"]
agg = me.groupby(["Sampling_Season","Maturity_stage","Infection_Group"])["ME_sugar"].mean().reset_index()
agg["stage_idx"] = agg["Maturity_stage"].map({s:i for i,s in enumerate(order)})
sns.lineplot(data=agg, x="stage_idx", y="ME_sugar",
             hue="Infection_Group", style="Sampling_Season", markers=True, ax=ax)
ax.set_xticks(range(4)); ax.set_xticklabels(order, rotation=20)
ax.axhline(0, color="grey", lw=0.5)
ax.set_ylabel("Sugar-transport module eigengene (PC1)")
ax.set_title("Carbon-translocation module across stage x season x infection")
plt.tight_layout()
plt.savefig(os.path.join(args.out_dir, "sugar_module_trajectory.png"), dpi=150)
plt.close(fig)
print("Done.")
