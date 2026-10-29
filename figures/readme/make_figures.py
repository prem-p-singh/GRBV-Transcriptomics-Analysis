"""Write the README figures as SVG sources into figures/readme/src/.

Values come from the PNAS draft (27 Oct 2026), its Table 2 and docs/ANALYSIS_ERAS.md.
Each source is then opened in Illustrator and exported to figures/readme/<name>.svg
with text converted to outlines, so GitHub renders it without the fonts.

Run from the repo root: python figures/readme/make_figures.py
"""
from pathlib import Path

SRC = Path(__file__).parent / "src"
FONT = "Helvetica, Arial, sans-serif"
INK, MUTED, GRID, CARD, EDGE = "#231F22", "#6F6670", "#E9E5E7", "#FFFFFF", "#E2DCDF"
HEALTHY, INFECTED, GOLD, SLATE = "#2F8A82", "#A23B4F", "#B8892E", "#5B6B8C"
STAGES = ["Pre-veraison", "Veraison", "Post-veraison", "Harvest"]


def t(x, y, s, size=14, fill=INK, weight="400", anchor="start"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
            f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{s}</text>')


def line(x1, y1, x2, y2, stroke=GRID, w=1.5, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{w}"{d}/>'


def dot(x, y, r, fill, stroke="none"):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="2"/>'


def rect(x, y, w, h, fill, rx=8, stroke="none"):
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}" stroke="{stroke}"/>'


def arrow(x1, y, x2, color, w=3):
    """Horizontal arrow from x1 to x2 with a filled head."""
    s = -1 if x2 < x1 else 1
    head = f'<path d="M{x2:.1f} {y:.1f} l{-10 * s} -6 l0 12 z" fill="{color}"/>'
    return line(x1, y, x2 - 8 * s, y, color, w) + head


def card(w, h, title, subtitle, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
            + rect(1, 1, w - 2, h - 2, CARD, 14, EDGE)
            + t(40, 52, title, 24, INK, "700") + t(40, 80, subtitle, 15, MUTED) + body + "</svg>\n")


def hero():
    # Developmental axis: 0 = pre-veraison ... 3 = harvest, healthy-only reference.
    x0, x1 = 230, 1060
    sx = lambda u: x0 + (x1 - x0) * u / 3
    rows = [("2021", 175, 0.804, 0.339), ("2022", 285, 0.148, 0.396)]
    body = []
    for i, s in enumerate(STAGES):
        body.append(line(sx(i), 120, sx(i), 330, GRID, 1.2))
        body.append(t(sx(i), 360, s, 13, MUTED, anchor="middle"))
    body.append(t(x0, 388, "Developmental position on the healthy-only reference (stage units)", 12, MUTED))
    for season, y, common, extra in rows:
        inf = 3 - common
        body += [t(60, y + 5, season, 18, INK, "700"),
                 line(x0, y, sx(3), y, "#D9D2D6", 3),
                 arrow(sx(3), y - 14, sx(inf), INFECTED, 3),
                 t((sx(3) + sx(inf)) / 2, y - 24, f"common delay &#8722;{common:.3f}", 12, INFECTED, "700",
                   "middle"),
                 arrow(sx(inf), y + 16, sx(inf - extra), GOLD, 3),
                 t(sx(inf - extra) - 8, y + 36, f"photosystem extra &#8722;{extra:.3f}", 12, GOLD, "700", "end"),
                 dot(sx(3), y, 9, HEALTHY), dot(sx(inf), y, 9, INFECTED)]
    legend = [dot(70, 430, 7, HEALTHY), t(84, 435, "Healthy at harvest", 13),
              dot(250, 430, 7, INFECTED), t(264, 435, "Infected at harvest (common shift)", 13),
              line(510, 430, 540, 430, GOLD, 3), t(550, 435, "Additional photosystem shift", 13)]
    note = [rect(40, 458, 1120, 40, "#F6F1F3", 8),
            t(60, 484, "Small heat-shock genes do not shift in time: 31 of 34 are lower at harvest in each season "
                       "(an expression deficit, not a later clock).", 13, INK)]
    return card(1200, 520, "Ripening, out of step",
                "Infected Merlot berries sit behind healthy ones, and photosystem genes fall further behind still",
                "".join(body + legend + note))


def design():
    body = []
    for i in range(44):  # 36 infected + 8 healthy vines
        col, row = i % 11, i // 11
        body.append(dot(60 + col * 30, 150 + row * 32, 10, INFECTED if i < 36 else HEALTHY))
    body += [t(60, 300, "44 vines", 20, INK, "700"),
             t(60, 324, "36 GRBV-positive, 8 GRBV-negative", 14, MUTED),
             t(60, 346, "Infection onset cohorts 2016 to 2022", 14, MUTED),
             arrow(400, 230, 470, MUTED, 2)]
    for j, s in enumerate(STAGES):
        for k, season in enumerate(["2021", "2022"]):
            body.append(rect(500 + j * 92, 140 + k * 70, 84, 56, "#F6F1F3", 8))
            body.append(t(542 + j * 92, 164 + k * 70, s, 11, INK, "700", "middle"))
            body.append(t(542 + j * 92, 182 + k * 70, season, 11, MUTED, anchor="middle"))
    body += [t(500, 300, "4 stages &#215; 2 seasons", 20, INK, "700"),
             t(500, 324, "At most one library per vine per cell", 14, MUTED),
             arrow(880, 230, 940, MUTED, 2),
             t(960, 205, "327", 48, INFECTED, "700"),
             t(960, 232, "RNA-seq libraries", 16, INK, "700"),
             t(960, 256, "264 infected, 63 healthy", 14, MUTED),
             t(960, 300, "11", 30, GOLD, "700"),
             t(1000, 300, "harvest chemistry pools", 14, INK),
             t(960, 324, "one pool per season and group,", 13, MUTED),
             t(960, 344, "measured in triplicate", 13, MUTED)]
    return card(1200, 390, "Study design",
                "Repeated sampling of the same vines: libraries are not independent plants", "".join(body))


def pipeline():
    steps = [("1", "Sequence and count", "TAG-seq, STAR,|featureCounts", "external"),
             ("2", "Filter and audit", "14,604 features;|vine-level units", "local"),
             ("3", "Model expression", "DESeq2, voom,|dream mixed models", "partly in repo"),
             ("4", "Networks and omics", "Snakemake rules|00 to 19 on FARM", "in repo"),
             ("5", "Clock and atlas", "Healthy-only clock,|70 branches", "local"),
             ("6", "Validate and draw", "3 published studies,|R figures", "in repo")]
    badge = {"in repo": HEALTHY, "partly in repo": GOLD, "local": SLATE, "external": MUTED}
    body = []
    for i, (n, title, sub, where) in enumerate(steps):
        x = 40 + i * 188
        body += [rect(x, 120, 172, 150, "#F6F1F3", 10),
                 dot(x + 30, 150, 16, INFECTED), t(x + 30, 156, n, 16, CARD, "700", "middle"),
                 t(x + 16, 196, title, 15, INK, "700")]
        for k, part in enumerate(sub.split("|")):
            body.append(t(x + 16, 220 + k * 18, part, 12, MUTED))
        body += [rect(x + 16, 286, 140, 26, badge[where], 13), t(x + 86, 304, where, 12, CARD, "700", "middle")]
        if i < 5:
            body.append(arrow(x + 174, 195, x + 186, MUTED, 2))
    return card(1200, 340, "How the analysis runs",
                "Six stages from reads to results; the badge shows where each stage lives", "".join(body))


def effects():
    rows = [("Small heat-shock proteins", "harvest", -2.294, -1.581),
            ("Stilbene synthesis", "harvest", -0.70, -1.70),
            ("Cytosolic translation", "harvest", -0.205, -0.302),
            ("Flavonoid hydroxylation", "harvest", 0.29, 0.82),
            ("Chloroplast translation", "harvest", 0.400, 0.440),
            ("Sucrose synthesis and cycling", "harvest", 0.67, 0.38),
            ("DCL2", "veraison", 1.610, 0.933),
            ("RDR1 (probable)", "veraison", 4.254, 3.045)]
    x0, x1, lo, hi = 420, 1140, -3, 5
    sx = lambda v: x0 + (x1 - x0) * (v - lo) / (hi - lo)
    body = []
    for v in range(lo, hi + 1):
        body.append(line(sx(v), 110, sx(v), 470, "#C9C1C5" if v == 0 else GRID, 2 if v == 0 else 1))
        body.append(t(sx(v), 492, f"{v:+d}" if v else "0", 12, MUTED, anchor="middle"))
    body.append(t((x0 + x1) / 2, 516, "Mean log2 fold change, infected vs healthy", 13, MUTED, anchor="middle"))
    for i, (name, stage, a, b) in enumerate(rows):
        y = 130 + i * 44
        body += [t(40, y + 5, name, 15, INK, "700"), t(300, y + 5, stage, 13, MUTED),
                 line(sx(a), y, sx(b), y, "#CFC6CB", 3), dot(sx(a), y, 8, INFECTED), dot(sx(b), y, 8, GOLD)]
    body += [dot(820, 70, 7, INFECTED), t(834, 75, "2021", 13), dot(900, 70, 7, GOLD), t(914, 75, "2022", 13)]
    return card(1200, 540, "Programs move in different directions",
                "Same sign in both seasons, different sizes", "".join(body))


def chemistry():
    data = {"2021": (26.73, [23.37, 22.70, 22.97, 22.83], "3.37 to 4.03"),
            "2022": (34.43, [25.50, 24.47, 23.50, 25.23, 24.30], "8.93 to 10.93")}
    y0, y1, lo, hi = 400, 110, 20, 36
    sy = lambda v: y0 - (y0 - y1) * (v - lo) / (hi - lo)
    body = []
    for v in range(lo, hi + 1, 4):
        body += [line(120, sy(v), 1100, sy(v), GRID, 1), t(108, sy(v) + 4, str(v), 12, MUTED, anchor="end")]
    body.append(t(40, 104, "Brix", 13, MUTED, "700"))
    for i, (season, (healthy, infected, deficit)) in enumerate(data.items()):
        cx = 360 + i * 480
        # The axis starts at 20 Brix, so the healthy pool is a reference line, not a bar.
        body += [line(cx - 140, sy(healthy), cx + 140, sy(healthy), HEALTHY, 4),
                 t(cx - 150, sy(healthy) + 5, f"{healthy:.2f}", 13, HEALTHY, "700", "end")]
        for k, v in enumerate(infected):
            xk = cx - 60 + k * 30
            body += [line(xk, sy(healthy), xk, sy(v), "#E3C9CF", 2), dot(xk, sy(v), 8, INFECTED)]
        body += [t(cx, 428, season, 16, INK, "700", "middle"),
                 t(cx + 100, (sy(healthy) + sy(min(infected))) / 2, f"deficit {deficit}", 13, INFECTED, "700")]
    body += [line(850, 47, 878, 47, HEALTHY, 4), t(886, 52, "Healthy pool", 13),
             dot(1000, 47, 7, INFECTED), t(1012, 52, "Infected cohort pools", 13)]
    return card(1200, 450, "Lower soluble solids in every infected pool",
                "Pooled harvest samples, triplicate measurements; descriptive only", "".join(body))


def eras():
    items = [("1", "Mar", "Welch t-test", "908 DEGs withdrawn"),
             ("2", "Apr", "Networks, multi-omics", "immunity flat"),
             ("3", "Apr", "Threshold sweeps", "mediator retired"),
             ("4", "Apr", "DESeq2 and voom", "908 to 126 genes"),
             ("5", "Apr to Jun", "External validation", "3 framings dropped"),
             ("6", "Aug", "Mixed models, claim gates", "every claim judged"),
             ("7", "Sep", "Clock and pathway atlas", "timing vs level")]
    body = [line(80, 170, 1120, 170, "#D9D2D6", 4)]
    for i, (n, when, what, out) in enumerate(items):
        x = 80 + i * (1040 / 6)
        body += [dot(x, 170, 18, INFECTED if i < 4 else HEALTHY), t(x, 176, n, 15, CARD, "700", "middle"),
                 t(x, 132, when, 12, MUTED, anchor="middle"),
                 t(x, 214, what, 12, INK, "700", "middle"), t(x, 234, out, 12, MUTED, anchor="middle")]
    body += [dot(80, 290, 7, INFECTED), t(94, 295, "Per-gene and network eras, mostly revised", 13),
             dot(560, 290, 7, HEALTHY), t(574, 295, "Validation and repeated-measures eras, current", 13)]
    return card(1200, 320, "Seven analysis eras, 2026",
                "Each era asked a better-posed question; counts from different eras are not comparable",
                "".join(body))


FIGURES = {"hero": hero, "design": design, "pipeline": pipeline, "effects": effects,
           "chemistry": chemistry, "eras": eras}

if __name__ == "__main__":
    SRC.mkdir(parents=True, exist_ok=True)
    for name, fn in FIGURES.items():
        (SRC / f"{name}.svg").write_text(fn())
    print(f"wrote {len(FIGURES)} sources to {SRC}")
