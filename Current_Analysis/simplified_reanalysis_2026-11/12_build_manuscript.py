#!/usr/bin/env python3
"""Fills Tables 1-3 from this revision's CSVs and builds the .docx with pandoc.
Pandoc takes its styles from the previous PNAS manuscript, so the look matches."""
import csv, pathlib, re, subprocess, sys, zipfile
ROOT = pathlib.Path(__file__).resolve().parents[2]
A = ROOT / "Current_Analysis/simplified_reanalysis_2026-11/data"
M = ROOT / "Current_Analysis/manuscripts/PNAS_2026-11_simplified"
STYLE = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else (ROOT / "Current_Analysis/prior_era_inputs/Era7_PNAS_Manuscript.docx")
REF = M / "style_reference.docx"
# Style-only copy: pandoc would otherwise carry over the old file's figure images.
with zipfile.ZipFile(STYLE) as src, zipfile.ZipFile(REF, "w", zipfile.ZIP_DEFLATED) as dst:
    for info in src.infolist():
        if info.filename.startswith("word/media/"):
            continue
        data = src.read(info.filename)
        if info.filename.endswith(".rels"):
            data = re.sub(rb'<Relationship [^>]*media/[^>]*/>', b"", data)
        dst.writestr(info, data)
STAGES = ["Pre_veraison", "Veraison", "Post_veraison", "Harvest"]
LABEL = dict(zip(STAGES, ["Pre-veraison", "Veraison", "Post-veraison", "Harvest"]))
rows = lambda f: list(csv.DictReader(open(A / f)))
num = lambda x, d=3: f"{float(x):+.{d}f}".replace("-", "−")
fdr = lambda x: f"{float(x):.4f}"

c = rows("sample_cell_counts.csv")
n = lambda i, y, s: next(r for r in c if r["infection"] == i and r["season"] == y and r["stage"] == s)["libraries"]
t1 = ["Table 1. RNA-sequencing libraries in each stage and season", "", "| Season | Stage | Healthy libraries | Infected libraries |", "| --- | --- | --- | --- |"]
for y in ("2021", "2022"):
    t1 += [f"| {y} | {LABEL[s]} | {n('RB(-)', y, s)} | {n('RB(+)', y, s)} |" for s in STAGES]
neg = sum(int(r["libraries"]) for r in c if r["infection"] == "RB(-)"); pos = sum(int(r["libraries"]) for r in c if r["infection"] == "RB(+)")
assert (neg, pos) == (63, 264)
t1 += [f"| Both | All | {neg} | {pos} |", "", "Libraries are repeated observations from 8 healthy and 36 infected vines."]

d = {r["cell"]: r for r in rows("dev_displacement_inference.csv")}
t2 = ["Table 2. Infected minus healthy developmental score", "", "| Season | Stage | Difference (stage units) | 95% CI | Permutation FDR |", "| --- | --- | --- | --- | --- |"]
for y in ("2021", "2022"):
    for s in STAGES:
        r = d[f"{s}__{y}"]
        t2.append(f"| {y} | {LABEL[s]} | {num(r['estimate'])} | {num(r['nested_ci_low'])} to {num(r['nested_ci_high'])} | {fdr(r['perm_fdr'])} |")
t2 += ["", "Intervals are from 1,999 nested vine bootstrap draws that rebuilt the reference. FDRs are from 4,999 vine-label permutations."]

g = {(r["group"], r["season"]): r for r in rows("headline_group_classes.csv") if r["stage"] == "Harvest"}
b = {(r["group"], r["season"]): r for r in rows("headline_group_method_agreement.csv") if r["stage"] == "Harvest"}
nl = {(r["group"], r["statistic"], r["season"]): r for r in rows("matched_null_percentiles.csv") if r["stage"] == "Harvest"}
names = [("sHSP", "Small heat-shock proteins"), ("Cytosolic_translation", "Cytosolic translation"), ("Chloroplast_translation", "Chloroplast translation"), ("Photosystem_light_harvesting", "Photosystem and light harvesting")]
t3 = ["Table 3. Headline gene groups at harvest", "", "| Group (genes) | Season | Mean log2FC, unadjusted | Mean log2FC, adjusted | Genes beyond delay | Adjusted effect, null P | Score offset | Offset, null P |", "| --- | --- | --- | --- | --- | --- | --- | --- |"]
for k, lab in names:
    for y in ("2021", "2022"):
        r, m, o = g[(k, y)], nl[(k, "mean_adjusted_effect", y)], nl[(k, "score_offset", y)]
        t3.append(f"| {lab} ({r['n_genes']}) | {y} | {num(r['mean_log2FC_unadjusted'])} | {num(r['mean_log2FC_adjusted'])} | {r['beyond_delay']} | {float(m['p_two_sided']):.3f} | {num(o['real_value'])} | {float(o['p_two_sided']):.3f} |")
t3 += ["", "Null P values are two-sided, from 1,999 matched random sets."]

src = (M / "manuscript.md").read_text()
assert "<!-- TABLE1 -->" in src
out = src.replace("<!-- TABLE1 -->", "\n".join(t1)).replace("<!-- TABLE2 -->", "\n".join(t2))
out = out.replace("## Discussion", "\n".join(t3) + "\n\n## Discussion", 1).replace("<!-- MAIN_TEXT -->\n", "").replace("<!-- END_MAIN_TEXT -->\n", "")
ch = rows("chemistry_summary.csv")
ts1 = ["| Season | Pool | Brix (mean ± SD) | pH (mean ± SD) | Brix deficit |", "| --- | --- | --- | --- | --- |"]
for r in sorted(ch, key=lambda r: (r["season"], r["group"] != "RB(-)", r["group"])):
    deficit = "" if r["group"] == "RB(-)" else f"{float(r['Brix_deficit']):.2f}"
    ts1.append(f"| {r['season']} | {r['group']} | {float(r['Brix']):.2f} ± {float(r['Brix_sd']):.2f} | {float(r['pH']):.3f} ± {float(r['pH_sd']):.3f} | {deficit} |")
si = (M / "SI_Appendix.md").read_text().replace("<!-- TABLES1 -->", "\n".join(ts1))
(M / "SI_Appendix_with_tables.md").write_text(si)
subprocess.run(["pandoc", str(M / "SI_Appendix_with_tables.md"), "-o", str(M / "PNAS_SI_Appendix_simplified.docx"), f"--reference-doc={REF}"], check=True)
full = M / "manuscript_with_tables.md"
full.write_text(out)
subprocess.run(["pandoc", str(full), "-f", "markdown-implicit_figures", "--resource-path=.", "-o", str(M / "PNAS_Manuscript_simplified.docx"), f"--reference-doc={REF}"], check=True)
print("built")
