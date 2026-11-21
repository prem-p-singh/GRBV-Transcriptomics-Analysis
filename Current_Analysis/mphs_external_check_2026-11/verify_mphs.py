"""Re-read every CSV from this check and confirm counts, draws, and that each decimal in MPHS_SUMMARY.md is in a CSV."""
import csv, glob, re, pathlib
here = pathlib.Path(__file__).parent
rd = lambda p: list(csv.DictReader(open(p, newline="")))

proj = rd(here / "data/mphs_all_projection.csv")
assert len(proj) == 327 and len({r["Sample"] for r in proj}) == 327, "library count"
assert sum(r["infection"] == "RB(-)" for r in proj) == 63, "healthy count"
assert len({r["vine"] for r in proj}) == 44, "vine count"
assert "2103" in {r["Sample"] for r in proj}, "sample 2103 kept as text"
assert len(rd(here / "data/mphs_healthy_projection.csv")) == 63
for f in ["mphs_delay_inference.csv", "mphs_independent_genes_delay.csv"]:
    for r in rd(here / "data" / f):
        assert int(r["n_boot"]) == 4999 and int(r["n_perm"]) == 4999, f
for r in rd(here / "data/mphs_season_difference.csv"):
    assert int(r["n_boot"]) == 4999

vals = []
for p in glob.glob(str(here / "data/*.csv")) + glob.glob(str(here / "../simplified_reanalysis_2026-11/data/dev_displacement_inference.csv")) + glob.glob(str(here / "../simplified_reanalysis_2026-11/data/headline_group_classes.csv")):
    for r in rd(p):
        for v in r.values():
            try: vals.append(float(v))
            except (TypeError, ValueError): pass
summary = (here / "logs/MPHS_SUMMARY.md").read_text()
summary = summary.split("## Results", 1)[1]
nums = re.findall(r"(?<![\w.])[-+]?\d+\.\d+", summary)
missing = []
for n in nums:
    dp = len(n.split(".")[1]); x = float(n)
    if not any(round(v, dp) == x or round(abs(v), dp) == abs(x) for v in vals):
        missing.append(n)
g = list(csv.DictReader(open(here / "data/mphs_adjusted_headline_groups.csv")))
assert len(g) == 8 and {r["n_genes"] for r in g} == {"34", "44", "95", "20"}, "step 5 groups"
print(f"libraries 327, healthy 63, vines 44, draws 4999: OK")
print(f"decimal numbers in summary: {len(nums)}; not found in any CSV: {missing or 'none'}")
assert not missing
print("VERIFY OK")
