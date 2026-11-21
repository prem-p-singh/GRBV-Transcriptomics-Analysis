#!/usr/bin/env python3
"""T14. Every number in the manuscript must sit inside a phrase rebuilt from its
source CSV or script. Numbers left uncovered, missing phrases, retired terms and
word-limit breaches all fail. Run from the repository root."""
import csv, gzip, pathlib, re, sys
import os
ROOT = pathlib.Path(__file__).resolve().parents[2]
A = ROOT / "Current_Analysis/simplified_reanalysis_2026-11"
D = A / "data"
M = pathlib.Path(os.environ.get("GRBV_MS_DIR", str(ROOT / "Current_Analysis/manuscripts/PNAS_2026-11_simplified")))
rows = lambda f: list(csv.DictReader(open(D / f, encoding="utf-8")))
s3 = lambda x, n=3: f"{float(x):.{n}f}".replace("-", "−")
sg = lambda x, n=3: f"{float(x):+.{n}f}".replace("-", "−")
fdr = lambda x: f"{float(x):.4f}".rstrip("0")
com = lambda x: f"{int(float(x)):,}"

inf = {r["cell"]: r for r in rows("dev_displacement_inference.csv")}
diff = rows("dev_displacement_season_difference.csv")[0]
val = {(r["validation"], r["score"]): r for r in rows("dev_score_validation.csv")}
summ = {(r["season"], r["stage"]): r for r in rows("adjusted_vs_unadjusted_summary.csv")}
agree = {(r["season"], r["stage"]): r for r in rows("adjustment_method_agreement.csv")}
cls = {(r["group"], r["season"]): r for r in rows("headline_group_classes.csv") if r["stage"] == "Harvest"}
meth = {(r["group"], r["season"]): r for r in rows("headline_group_method_agreement.csv") if r["stage"] == "Harvest"}
nl = {(r["group"], r["statistic"], r["season"]): r for r in rows("matched_null_percentiles.csv") if r["stage"] == "Harvest"}
tight = {(r["group"], r["statistic"], r["season"]): r for r in rows("matched_null_percentiles_tight.csv") if r["stage"] == "Harvest"}
perm = {r["test"]: r for r in rows("joint_test_permutation_null.csv")}
draws = rows("joint_test_permutation_draws.csv")
reff = {(r["model"], r["test"]): r for r in rows("dream_reffect_comparison.csv")}
ov = [r for r in rows("interaction_overlap.csv") if r["viral"].upper() == "FALSE"][0]
ext = {r["stage"]: r for r in rows("external_stage_concordance.csv")}
named = {(r["symbol"], r["stage"], r["season"]): r for r in rows("named_gene_effects.csv")}
atlas = {(r["name"], r["stage"], r["season"]): r for r in rows("atlas_branch_cell_summary.csv")}
chem = rows("chemistry_summary.csv")
audit = rows("sample_audit.csv")
sett = rows("settest_old_vs_new.csv")
zen = rows("zenith_unadjusted.csv")
wi = {(r["outcome"], r["term"]): r for r in rows("within_infected_gradients.csv")}
vc = rows("variance_components.csv")
ti = rows("titer_vs_viralcpm.csv")[0]
ms = {(r["cell"]): r for r in rows("dev_displacement.csv")}
qt = {r["group"]: r for r in rows("matched_null_quality_tight.csv")}
extb = {r["set_id"]: r for r in rows("external_branch_directions_SI.csv")}
bounds = {(r["score"], r["infection"], r["season"], r["stage"]): r for r in rows("dev_score_bounds.csv")}

neg = [float(r["viral_cpm"]) for r in audit if r["infection"] == "RB(-)"]
pos = [float(r["viral_cpm"]) for r in audit if r["infection"] == "RB(+)"]
vines = {r["vine"]: r["infection"] for r in audit}
nv_pos, nv_neg = sum(v == "RB(+)" for v in vines.values()), sum(v == "RB(-)" for v in vines.values())
HS, H1, H2 = "Harvest", "Harvest__2021", "Harvest__2022"
g = lambda k, y, c: cls[(k, y)][c]
SH, CY, CH, PS = "sHSP", "Cytosolic_translation", "Chloroplast_translation", "Photosystem_light_harvesting"
adj, off = "mean_adjusted_effect", "score_offset"
pv = lambda d, k, st, y: f"{float(d[(k, st, y)]['p_two_sided']):.3f}"
assert fdr(inf[H1]["perm_fdr"]) == "0.0016"
ci = lambda c: f"{s3(inf[c]['nested_ci_low'])} to {sg(inf[c]['nested_ci_high'])}"
ci2 = lambda c: f"{s3(inf[c]['nested_ci_low'])} to {s3(inf[c]['nested_ci_high'])}"
brix = lambda y: sorted(float(r["Brix_deficit"]) for r in chem if r["season"] == y and r["group"] != "RB(-)")
b21, b22 = brix("2021"), brix("2022")
n_same = sum(r["same_direction"].upper() == "TRUE" for r in sett)
n_sig = sum(float(r["FDR_all"]) < .05 for r in zen)
sh_p = {r["season"]: r for r in sett if r["set_id"] == "BRANCH:50" and r["stage_code"] == "Harvest"}
at = lambda n, c, y: atlas[(n, HS, y)][c]
mx = max(draws, key=lambda r: int(r["infection_stage"]))
OLDB = bounds[("old_score", "neg", "2021", HS)], bounds[("old_score", "neg", "2022", HS)]
limit = sum(int(r["exactly_0"]) + int(r["exactly_3"]) for r in bounds.values() if r["score"] == "new_score")
assert limit == 0, "a library sits at a limit"
assert int(ov["n_three_way"]) - int(ov["n_both"]) == int(ov["n_three_way_only"])

PHRASES = [
 # Abstract
 f"327 berry RNA-sequencing libraries from {len(vines)} Merlot vines ({nv_pos} infected, {nv_neg} healthy)",
 f"by {s3(inf[H1]['estimate'])} stage units in 2021 (95% CI, {ci2(H1)}) and {s3(inf[H2]['estimate'])} in 2022 ({ci(H2)})",
 f"{g(SH,'2021','beyond_delay')} of {g(SH,'2021','n_genes')} small heat-shock protein genes in 2021 and {g(SH,'2022','beyond_delay')} of {g(SH,'2022','n_genes')} in 2022",
 f"at {com(ov['n_infection_stage'])} and {com(ov['n_three_way'])} host genes",
 # Results 1
 f"327 libraries from {len(vines)} vines, {nv_pos} infected and {nv_neg} healthy",
 f"at most {max(neg):.2f} viral counts per million and every infected library at least {min(pos):.1f}",
 f"held {OLDB[0]['exactly_3']} of {OLDB[0]['n']} and {OLDB[1]['exactly_3']} of {OLDB[1]['n']} healthy harvest libraries",
 f"for {val[('held_out_vine','new_score')]['nearest_correct']} of {val[('held_out_vine','new_score')]['n']} healthy libraries, with a root mean squared error (RMSE) of {s3(val[('held_out_vine','new_score')]['RMSE'])}",
 f"for {val[('held_out_season','new_score')]['nearest_correct']} of {val[('held_out_season','new_score')]['n']}, with an RMSE of {s3(val[('held_out_season','new_score')]['RMSE'])}",
 f"was {s3(inf[H1]['estimate'])} stage units in 2021 (95% CI, {ci2(H1)}; FDR {fdr(inf[H1]['perm_fdr'])}) and {s3(inf[H2]['estimate'])} in 2022 ({ci(H2)}; FDR {fdr(inf[H2]['perm_fdr'])})",
 f"({s3(inf['Post_veraison__2021']['estimate'])}; {ci('Post_veraison__2021')}; FDR {fdr(inf['Post_veraison__2021']['perm_fdr'])})",
 f"(difference {s3(diff['estimate'])}; 95% CI, {s3(diff['nested_ci_low'])} to {s3(diff['nested_ci_high'])}; permutation P = {float(diff['perm_p']):.4f})",
 f"Brix deficits of {b21[0]:.2f} to {b21[-1]:.2f} in 2021 and {b22[0]:.2f} to {b22[-1]:.2f} in 2022",
 # Results 2
 f"In 2021, {com(summ[('2021',HS)]['n_unadjusted'])} genes were significant before adjustment; {com(summ[('2021',HS)]['n_explained_by_delay'])} of these were no longer significant after it, whereas {com(summ[('2021',HS)]['n_beyond_delay'])} genes",
 f"In 2022 the matching counts were {com(summ[('2022',HS)]['n_unadjusted'])}, {com(summ[('2022',HS)]['n_explained_by_delay'])}, and {com(summ[('2022',HS)]['n_beyond_delay'])}",
 f"it found {com(agree[('2021',HS)]['n_covariate'])} genes and {com(agree[('2021',HS)]['n_both'])} were shared",
 f"{int(g(SH,'2021','beyond_delay'))+int(g(SH,'2021','explained_by_delay'))} of {g(SH,'2021','n_genes')} tested genes were lower at harvest in each season, with mean log2 fold changes of {s3(g(SH,'2021','mean_log2FC_unadjusted'))} and {s3(g(SH,'2022','mean_log2FC_unadjusted'))}",
 f"the means were {s3(g(SH,'2021','mean_log2FC_adjusted'))} and {s3(g(SH,'2022','mean_log2FC_adjusted'))}, and {g(SH,'2021','beyond_delay')} and {g(SH,'2022','beyond_delay')} of {g(SH,'2021','n_genes')} genes",
 f"(P = {pv(nl,SH,adj,'2021')} in each season)",
 f"Cytosolic means were {s3(g(CY,'2021','mean_log2FC_unadjusted'))} and {s3(g(CY,'2022','mean_log2FC_unadjusted'))} before adjustment and {s3(g(CY,'2021','mean_log2FC_adjusted'))} and {s3(g(CY,'2022','mean_log2FC_adjusted'))} after it",
 f"Chloroplast means were {sg(g(CH,'2021','mean_log2FC_unadjusted'))} and {sg(g(CH,'2022','mean_log2FC_unadjusted'))} before and {sg(g(CH,'2021','mean_log2FC_adjusted'))} and {sg(g(CH,'2022','mean_log2FC_adjusted'))} after",
 f"both were clear in 2022 (P = {pv(nl,CY,adj,'2022')})",
 f"borderline in 2021 (P = {pv(nl,CY,adj,'2021')})",
 f"not significant in 2021 (P = {pv(nl,CH,adj,'2021')})",
 f"at {s3(nl[(PS,off,'2021')]['real_value'])} stage units in 2021 and {s3(nl[(PS,off,'2022')]['real_value'])} in 2022",
 f"fell outside the matched null (P = {pv(nl,PS,off,'2021')})",
 f"fell from {sg(g(PS,'2021','mean_log2FC_unadjusted'))} to {sg(g(PS,'2021','mean_log2FC_adjusted'))} in 2021 and from {sg(g(PS,'2022','mean_log2FC_unadjusted'))} to {sg(g(PS,'2022','mean_log2FC_adjusted'))} in 2022",
 f"only {g(PS,'2021','beyond_delay')} and {g(PS,'2022','beyond_delay')} of {g(PS,'2021','n_genes')} genes",
 f"{n_sig} of {com(len(zen))} set-by-cell tests reached FDR 0.05, and all {n_sig}",
 f"in {100*n_same/len(sett):.1f}% of {len(sett)} branch-by-cell comparisons",
 f"(P = {float(sh_p['2021']['zenith_unadj_p']):.3f} and {float(sh_p['2022']['zenith_unadj_p']):.4f} before adjustment)",
 f"({s3(at('Stilbene synthesis','mean_log2FC_unadjusted','2021'))} and {s3(at('Stilbene synthesis','mean_log2FC_unadjusted','2022'))})",
 f"({sg(at('Stilbene synthesis','mean_log2FC_adjusted','2021'))} and {s3(at('Stilbene synthesis','mean_log2FC_adjusted','2022'))})",
 f"({sg(named[('DCL2','Veraison','2021')]['log2FC_unadj'])} and {sg(named[('DCL2','Veraison','2022')]['log2FC_unadj'])}; {sg(named[('RDR1','Veraison','2021')]['log2FC_unadj'])} and {sg(named[('RDR1','Veraison','2022')]['log2FC_unadj'])})",
 # Results 3
 f"averaged over seasons, identified {com(ov['n_infection_stage'])} host genes at FDR 0.05, and the infection-by-stage-by-season test identified {com(ov['n_three_way'])}. Of the latter, {com(ov['n_three_way_only'])}",
 f"in {len(draws)} full model refits",
 f"The null median was {int(float(perm['infection_stage']['null_median']))} for each test, and the 95th percentile was {float(perm['infection_stage']['null_p95']):.2f}, {int(float(perm['infection_stage_season']['null_p95']))}, and {int(float(perm['only_season_specific']['null_p95']))} genes",
 f"The largest null draw, {com(mx['infection_stage'])} genes, kept {mx['n_real_controls_as_pseudo']} of the {nv_neg} true controls",
 f"({com(reff[('vine_season_row','infection_stage')]['n_fdr05'])} and {com(reff[('vine_season_row','infection_stage_season')]['n_fdr05'])} genes)",
 # Results 4
 f"rose from {s3(ext['Pre_veraison']['spearman_rho'])} before veraison to {s3(ext['Veraison']['spearman_rho'])} at veraison, {s3(ext['Post_veraison']['spearman_rho'])} after veraison, and {s3(ext[HS]['spearman_rho'])} at harvest, and {float(ext[HS]['sign_agreement_pct']):.1f}% of {com(ext[HS]['n_genes'])} matched genes",
 f"There are {nv_neg} control vines",
]
# Fixed design and method values, each checked against the script that sets it.
SCRIPTS = {"common": A/"00_common.R", "t3": A/"03_displacement_inference.R", "t5": A/"05_matched_null.R", "t6": A/"06_zenith_sets.R"}
def src(key, pat):
    assert re.search(pat, SCRIPTS[key].read_text()), f"{key}: {pat}"
src("common", r"N_MARKERS <- 1000L"); src("common", r"MIN_DYNAMIC <- 0\.5"); src("common", r"MIN_SCALE <- 0\.25")
src("t3", r"1999L"); src("t3", r"4999L"); src("t5", r"B <- 1999L"); src("t5", r"N_NEAREST <- 25L"); src("t6", r"== 1369"); src("t6", r"n_branch == 59")
assert len(draws) == 200 and len({r["draw"] for r in draws}) == 200
rownums = sorted({int(r["row_number"]) for r in audit}); assert (rownums[0], rownums[-1]) == (19, 22)
cohorts = sorted({r["cohort"] for r in audit if r["cohort"] != "-"}); assert cohorts == ["2016", "2019", "2020", "2021", "2022"]
seasons = {}
for r in audit: seasons.setdefault(r["vine"], set()).add(r["season"])
two, one = sum(len(v) == 2 for v in seasons.values()), sum(len(v) == 1 for v in seasons.values())
assert (two, one) == (40, 4)
PHRASES += [
 f"{two} vines were sampled in both seasons and {one} only in 2022", "in rows 19 to 22", "labeled 2016, 2019, 2020, 2021, and 2022",
 "Of 28,168 features, 14,604 had at least 10 counts in at least 32 libraries", "Nine metadata rows",
 "over 16 infection, stage, and season cells, testing eight", "8 of the 40 two-season vines were assigned as pseudo-controls in each of 200 draws",
 "up to 1,000 markers", "centroid range at least 0.5; residual SD floor 0.25", "1,999 nested bootstrap draws", "Permutation used 4,999 draws, with 8 of the 40 two-season vines",
 "in each of 1,999 draws", "from the 25 nearest genes", "on 1,369 sets of 5 to 500 genes across eight cells, including 59 branches", "The 70 branch definitions",
 "giving 11 pools", "excluded genes in 70 annotation-defined", "Filtering retained 14,604 features, including 14,599 host and 5 viral features",
 "held 4 of 8", "stopped scores at 0 and 3", "27 and 30 of 34 genes", "the 10% threshold", "(1 + count of null values", "(B + 1)", "1,999 matched",
]
# Legends, and values typed into the figure script.
fig = (A / "13_figures.R").read_text(encoding="utf-8")
for pat in [r"\[1:30\]", r"n_genes >= 10", r"LIM <- 1\.5", r"1,999 matched sets", r"200 permuted refits"]:
    assert re.search(pat, fig), pat
for typed in (f"{g(SH,'2021','beyond_delay')} and {g(SH,'2022','beyond_delay')} of {g(SH,'2021','n_genes')} genes",
              f"{g(PS,'2021','beyond_delay')} and {g(PS,'2022','beyond_delay')} of {g(PS,'2021','n_genes')} genes remain"):
    assert typed in fig, typed
PHRASES += [
 f"repeated sampling of {len(vines)} vines over two seasons", f"in 2021 and 2022 from {len(vines)} vines",
 "Bars are 95% intervals from 1,999 nested vine bootstrap draws",
 f"for {g(SH,'2021','n_genes')} small heat-shock protein, {g(CY,'2021','n_genes')} cytosolic translation, {g(CH,'2021','n_genes')} chloroplast translation, and {g(PS,'2021','n_genes')} photosystem and light-harvesting genes",
 "compared with 1,999 matched random gene sets",
 f"Counts from {len(draws)} refits",
 "for the 30 branches of at least 10 genes", "Values beyond ±1.5 are clipped",
]
tx = (M / "manuscript.md").read_text(encoding="utf-8")
head, legends = tx.split("## References")[0], tx.split("## Figure legends")[1]
body = head + legends
body = re.sub(r"\(\d+(?:, ?\d+)*\)", " ", body) if False else body
missing = [p for p in PHRASES if p not in body]
masked = body
for p in sorted(set(PHRASES), key=len, reverse=True):
    masked = masked.replace(p, " " * len(p))
strip = [r"\[AUTHOR[^\]]*\]", r"!\[\]\([^)]*\)(\{[^}]*\})?", r"`[^`]*`", r"LOC\d+", r"HSP18\.2", r"HSFA2e", r"DCL2", r"RDR1", r"log2|log10", r"ASM3070453v1", r"GSE\d+",
         r"Figs?\. \d[A-C]?(?: (?:A|and|B))*", r"\*\*Fig\. \d\.\*\*", r"Tables? S?\d", r"Datasets? S\d", r"\(\d+(?:, \d+)*\)", r"\(3; ", r"<!--.*?-->"]
for pat in strip:
    masked = re.sub(pat, " ", masked)
tokens = re.findall(r"[−+-]?\d[\d,]*(?:\.\d+)?%?", masked)
# Season labels and plain design terms left in running text.
ALLOW = {"2021": "season", "2022": "season", "95%": "interval level", "0.05": "FDR threshold", "0": "lower score bound", "3": "upper score bound", "4": "stage count in legend or 4 onlysee", }
ALLOW.pop("4")
left = [t for t in tokens if t.strip("+-−,") not in ALLOW and t.rstrip(",") not in ALLOW]
fails = []
print(f"Phrases built from sources: {len(PHRASES)}; missing from manuscript: {len(missing)}")
for p in missing: print("  MISSING:", p); fails.append(p)
print(f"Numbers left uncovered: {len(left)}")
for t in left:
    i = masked.find(t); print("  UNCOVERED:", t, "|", re.sub(r"\s+", " ", masked[max(0, i-50):i+40])); fails.append(t)

# SI Appendix numbers.
si = (M / "SI_Appendix.md").read_text(encoding="utf-8")
b50 = extb["BRANCH:50"]
SI = [f"({s3(b50['mean_external'])} vs {s3(b50['mean_current'])}, season-averaged Merlot, {b50['n']} matched genes)",
      f"(Spearman {float(ti['spearman_rho']):.3f}, n = {ti['n']})",
      f"({s3(wi[('developmental_offset','linear trend per cohort year')]['estimate'])} per cohort year, FDR {float(wi[('developmental_offset','linear trend per cohort year')]['fdr']):.3f})",
      f"was not significant (FDR {float(wi[('developmental_offset','slope per log10')]['fdr']):.3f})",
      f"(FDR at least {min(float(wi[('shsp_adjusted_residual',t)]['fdr']) for t in ('slope per log10','joint F, 3 df','linear trend per cohort year')):.2f})",
      f"harvest differences {s3(ms[H1]['spline_score'])} and {s3(ms[H2]['spline_score'])}",
      f"with extrapolation: {s3(ms[H1]['seasonal_extrap'])} and {s3(ms[H2]['seasonal_extrap'])}",
      f"({float(qt[PS]['null_mean_logcpm']):.2f} vs {float(qt[PS]['real_mean_logcpm']):.2f} log-CPM)",
      f"| {s3(inf[H1]['old_estimate'])} ({s3(inf[H1]['old_ci_low'])} to {s3(inf[H1]['old_ci_high'])}) | {s3(inf[H1]['estimate'])} ({ci2(H1)}) |",
      f"| {s3(inf[H2]['old_estimate'])} ({s3(inf[H2]['old_ci_low'])} to {s3(inf[H2]['old_ci_high'])}) | {s3(inf[H2]['estimate'])} ({ci(H2)}) |",
      f"| {nl[(PS,off,'2021')]['real_value'] and ''}"[:0] + f"{s3(nl[(PS,off,'2021')]['real_value'])}, {s3(nl[(PS,off,'2022')]['real_value'])} |"]
for p in SI:
    if p not in si: print("  SI MISSING:", p); fails.append(p)
print(f"SI phrases checked: {len(SI)}")
med = lambda k: sorted(float(r[k]) for r in vc)[len(vc)//2]
assert f"{med('stage'):.3f}" == "0.420", med("stage")

# Terms and limits (T14 items 2 and 3).
BANNED = ["908", "AUC", "MBF1c", "three phases", "temperature caused", "temperature", "climate", "memory", "scar", "tolerance gradient", "thermostat", "60%", "0.6 °Bx", "early burden"]
main = tx.split("## References")[0] + legends + si
for w in BANNED:
    pat = r"\b" + re.escape(w) + r"\b"
    n = len(re.findall(pat, main, re.I))
    if n: print("  BANNED TERM:", w, n); fails.append(w)
print("Retired terms checked:", len(BANNED))
if "—" in main: fails.append("em dash"); print("  EM DASH")
def words(a, b):
    t = tx[tx.index(a) + len(a):tx.index(b)]
    return len(re.sub(r"\[AUTHOR[^\]]*\]", "", t).split())
wc = {"Abstract": (words("## Abstract", "## Significance"), 250), "Significance": (words("## Significance Statement", "Keywords:"), 120),
      "Methods": (words("## Materials and Methods", "<!-- END_MAIN_TEXT -->"), 950)}
for k, (n, lim) in wc.items():
    print(f"{k}: {n} words (limit {lim})")
    if n > lim: fails.append(k)
print("RESULT:", "FAIL" if fails else "PASS", f"({len(fails)} problems)")
sys.exit(1 if fails else 0)
