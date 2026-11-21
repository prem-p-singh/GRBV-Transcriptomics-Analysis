#!/usr/bin/env python3
"""Generate clean, illustrated visual SVG figures for README.md in Manuscript/figures/readme/.

Uses fresh filenames to bust GitHub's SVG image cache and centers the hero graphic on
stage-by-stage GRBV disease pathogenesis and disrupted host pathways in ripening berries.
"""
from pathlib import Path

OUT = Path(__file__).resolve().parent
FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"
INK, MUTED, GRID, CARD, EDGE = "#1E1A20", "#645D66", "#EAE6E8", "#FFFFFF", "#E2DCDF"
TEAL, BURG, GOLD, SLATE, GREEN = "#1F8A84", "#9E2A47", "#B8862D", "#54657E", "#4D8B55"
TEAL_BG, BURG_BG, GOLD_BG, SLATE_BG, GREEN_BG = "#EEF7F6", "#FBF1F4", "#FCF7EC", "#F2F4F8", "#EFF7F0"


def t(x, y, s, size=14, fill=INK, weight="400", anchor="start"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
            f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{s}</text>')


def line(x1, y1, x2, y2, stroke=GRID, w=1.5, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{w}"{d}/>'


def poly(pts, stroke, w=2.8, fill="none"):
    p = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    return f'<polyline points="{p}" fill="{fill}" stroke="{stroke}" stroke-width="{w}" stroke-linecap="round" stroke-linejoin="round"/>'


def dot(x, y, r, fill, stroke="none", sw=1.8):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'


def rect(x, y, w, h, fill, rx=12, stroke="none", sw=1.2):
    s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke != "none" else ""
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}"{s}/>'


def pill(x, y, w, h, label, bg, fg, size=11.5, weight="700"):
    return rect(x, y, w, h, bg, h / 2) + t(x + w / 2, y + h * 0.68, label, size, fg, weight, "middle")


def arrow(x1, y, x2, color, w=2.4):
    s = -1 if x2 < x1 else 1
    head = f'<path d="M{x2:.1f} {y:.1f} l{-8.5 * s} -5 l0 10 z" fill="{color}"/>'
    return line(x1, y, x2 - 6.5 * s, y, color, w) + head


def v_arrow(x, y1, y2, color, w=2.2):
    s = -1 if y2 < y1 else 1
    head = f'<path d="M{x:.1f} {y2:.1f} l-5 {-8.5 * s} l10 0 z" fill="{color}"/>'
    return line(x, y1, x, y2 - 6.5 * s, color, w) + head


def grape_cluster(cx, cy, colors):
    """Draw a stylized 10-berry grape cluster with stem at (cx, cy)."""
    stem = (f'<path d="M{cx:.1f} {cy - 24:.1f} l0 10 M{cx - 8:.1f} {cy - 22:.1f} q8 4 16 0" '
            f'stroke="#5E6B54" stroke-width="2.2" fill="none" stroke-linecap="round"/>')
    offsets = [
        (-12, -8), (0, -8), (12, -8),
        (-17, 2), (-6, 2), (6, 2), (17, 2),
        (-11, 12), (1, 12), (11, 12),
        (-5, 22), (5, 22),
    ]
    berries = []
    for idx, (dx, dy) in enumerate(offsets):
        col = colors[idx % len(colors)]
        berries.append(dot(cx + dx, cy + dy, 6.8, col, CARD, 1.2))
    return stem + "".join(berries)


def leaf_virus_icon(cx, cy):
    """Stylized red-blotched grape leaf + circular ssDNA GRBV icon."""
    leaf = (f'<path d="M{cx:.1f} {cy - 18:.1f} c14 -8 26 6 20 20 c-4 10 -14 16 -20 20 '
            f'c-6 -4 -16 -10 -20 -20 c-6 -14 6 -28 20 -20 z" fill="#6A996E" stroke="#3F6A44" stroke-width="1.5"/>')
    blotch1 = dot(cx - 7, cy - 3, 5.5, BURG)
    blotch2 = dot(cx + 8, cy + 3, 6.5, BURG)
    blotch3 = dot(cx - 2, cy + 10, 4.5, "#B83B59")
    vein = line(cx, cy - 15, cx, cy + 18, "#E8F2E8", 1.3)
    return leaf + blotch1 + blotch2 + blotch3 + vein


def card(w, h, title, subtitle, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">'
            + rect(1, 1, w - 2, h - 2, CARD, 16, EDGE)
            + t(38, 45, title, 22, INK, "700")
            + t(38, 70, subtitle, 13.8, MUTED)
            + body + "</svg>\n")


def graphical_abstract():
    body = []
    # Top Phloem / Virus Banner connecting Source Leaf to Developing Berries
    body += [
        rect(38, 90, 1124, 62, "#F8F4F6", 12, "#E6D8DE"),
        leaf_virus_icon(78, 120),
        t(114, 116, "Phloem-Limited GRBV Infection", 14.5, BURG, "700"),
        t(114, 136, " Traps leaf carbon &#8594; lowers must sugar (&#8722;3.4 to &#8722;10.9 &#176;Bx) &amp; delays ripening (&#8722;0.81 stage lag in 2021)", 12.5, INK),
        pill(918, 108, 226, 26, "681 Host Genes Coupled to Viral RNA", BURG, CARD, 11),
    ]

    # Four Ripening Stage Cards with Illustrated Grape Clusters
    stages = [
        ("Pre-veraison", "Metabolic Priming", ["#63A36B", "#52935A"], GREEN, GREEN_BG,
         "+ RDR1 &amp; NAC62 defense", "+ RSL1 ubiquitin ligase", "&#8722; Nitrate &amp; TCA enzymes", "349 dysregulated genes"),
        ("Veraison", "Antiviral Silencing Peak", ["#63A36B", "#9E2A47", "#73AD78", "#B84A66"], TEAL, TEAL_BG,
         "+ RDR1 (+4.4) &amp; DCL2 (+1.6)", "+ NTM1-like 9 NAC (+4.3)", "160 NLRs kept off (no death)", "74 dysregulated genes"),
        ("Post-veraison", "Viral C2 Switch", ["#9E2A47", "#7D2038", "#B33958"], GOLD, GOLD_BG,
         "+ Viral C2 suppressor (+1.9)", "DCL2 silencing subsides", "&#8722; XTH7 wall, proteasome, GST", "957 dysregulated genes"),
        ("Harvest", "Proteostasis Collapse", ["#6E172E", "#8A203B", "#591124"], BURG, BURG_BG,
         "&#8722; 30/34 sHSP chaperones", "&#8722; HSFA2e &amp; HSP18.2 (&#8722;3.9)", "+ Metacaspase-9 (+3.6)", "1,069 dysregulated genes"),
    ]
    w_box, h_box, gap = 260, 268, 28
    x0, y0 = 38, 170
    for i, (st, phase, cols, fg, bg, l1, l2, l3, count_lbl) in enumerate(stages):
        x = x0 + i * (w_box + gap)
        body += [
            rect(x, y0, w_box, h_box, CARD, 14, fg if i == 3 else EDGE, 2.0 if i == 3 else 1.3),
            rect(x + 1, y0 + 1, w_box - 2, 42, bg, 13),
            t(x + w_box / 2, y0 + 27, st, 15, fg, "700", "middle"),
            grape_cluster(x + w_box / 2, y0 + 82, cols),
            t(x + w_box / 2, y0 + 132, phase, 13.8, INK, "700", "middle"),
            line(x + 20, y0 + 144, x + w_box - 20, y0 + 144, GRID, 1.2),
            t(x + 18, y0 + 168, l1, 12.5, BURG if l1.startswith("+") else TEAL, "600"),
            t(x + 18, y0 + 191, l2, 12.5, INK),
            t(x + 18, y0 + 214, l3, 12.5, TEAL if l3.startswith("&#8722;") else MUTED, "600"),
            pill(x + 34, y0 + 230, w_box - 68, 24, count_lbl, bg, fg, 11.2),
        ]
        if i < 3:
            body.append(arrow(x + w_box + 4, y0 + h_box / 2, x + w_box + gap - 4, SLATE, 2.3))

    return card(
        1200, 460,
        "GRBV Disease Pathogenesis Across Grape Berry Ripening",
        "How viral C2 activation and host pathway disruption progress from early antiviral silencing to harvest chaperone collapse",
        "".join(body),
    )


def study_overview():
    body = []
    body.append(rect(38, 94, 310, 224, "#FAF9FA", 12, EDGE))
    for i in range(44):
        col, row = i % 11, i // 11
        body.append(dot(64 + col * 25.5, 124 + row * 24, 7.8, BURG if i < 36 else TEAL))
    body += [
        t(58, 238, "44 Merlot Vines", 18.5, INK, "700"),
        t(58, 260, "Clone 12 on 1103P · Paso Robles, CA", 12.2, MUTED),
        dot(64, 282, 5.8, BURG), t(76, 286, "36 GRBV(+) vines (2016–2022)", 12.2, INK, "600"),
        dot(64, 304, 5.8, TEAL), t(76, 308, "8 healthy controls (vine unit)", 12.2, INK, "600"),
        arrow(356, 206, 396, SLATE, 2.4),
    ]

    body.append(rect(404, 94, 404, 224, "#FAF9FA", 12, EDGE))
    dates = {
        ("2021", 0): "2 Aug", ("2021", 1): "16 Aug", ("2021", 2): "8 Sep", ("2021", 3): "21 Sep",
        ("2022", 0): "26 Jul", ("2022", 1): "9 Aug", ("2022", 2): "24 Aug", ("2022", 3): "13 Sep",
    }
    for j, s in enumerate(["Pre-Ver", "Veraison", "Post-Ver", "Harvest"]):
        for k, season in enumerate(["2021", "2022"]):
            rx, ry = 422 + j * 93, 112 + k * 68
            bg = TEAL_BG if k == 0 else BURG_BG
            fg = TEAL if k == 0 else BURG
            body += [
                rect(rx, ry, 85, 56, bg, 8, EDGE),
                t(rx + 42.5, ry + 20, s, 11.2, INK, "700", "middle"),
                t(rx + 42.5, ry + 36, season, 11.2, fg, "700", "middle"),
                t(rx + 42.5, ry + 50, dates[(season, j)], 10.2, MUTED, anchor="middle"),
            ]
    body += [
        t(422, 276, "Same-Day Leaf &amp; Berry Collection", 16, INK, "700"),
        t(422, 298, "Deseeded skin + pulp (Qiagen Midi) &amp; petioles", 12.2, MUTED),
        arrow(816, 206, 856, SLATE, 2.4),
    ]

    body += [
        rect(864, 94, 298, 68, BURG_BG, 10, "#E9CFD7"),
        t(882, 134, "327", 29, BURG, "700"),
        t(948, 122, "Berry 3' TAG-seq", 13.8, INK, "700"),
        t(948, 142, "Element AVITI · UMI · 14,604 genes", 11.5, MUTED),

        rect(864, 172, 298, 68, TEAL_BG, 10, "#CFE5E3"),
        t(882, 212, "267", 29, TEAL, "700"),
        t(948, 200, "Leaf-Petiole dPCR", 13.8, INK, "700"),
        t(948, 220, "261 paired · copies/&#181;L DNA", 11.5, MUTED),

        rect(864, 250, 298, 68, GOLD_BG, 10, "#E8DAB8"),
        t(892, 290, "11", 29, GOLD, "700"),
        t(948, 278, "Harvest Must Pools", 13.8, INK, "700"),
        t(948, 298, "Brix, TA, pH, malate (triplicate)", 11.5, MUTED),
    ]
    return card(
        1200, 340,
        "Field Cohort &amp; Paired Sampling",
        "Repeated longitudinal measurements on the same 44 Merlot vines across two growing seasons",
        "".join(body),
    )


def host_pathways():
    rows = [
        ("RDR1 antiviral polymerase", "Veraison", 4.37, 3.10),
        ("NTM1-like 9 NAC factor*", "Veraison", 4.30, 2.54),
        ("RSL1 ABA-receptor ligase", "Veraison", 4.23, 2.16),
        ("Metacaspase-9 protease", "Harvest", 3.57, 2.81),
        ("DCL2 Dicer endonuclease", "Veraison", 1.61, 0.93),
        ("Cytosolic ribosomes (95)", "Harvest", -0.09, -0.27),
        ("Small HSP chaperones (34)", "Harvest", -1.85, -1.58),
        ("XTH7 wall xyloglucanase*", "Post-veraison", -2.43, -1.88),
        ("HSP18.2 class I chaperone", "Harvest", -3.87, -2.35),
    ]
    x0, x1, lo, hi = 430, 1140, -4.5, 5.0
    sx = lambda v: x0 + (x1 - x0) * (v - lo) / (hi - lo)
    body = [
        rect(x0, 98, sx(0) - x0, 382, "#F4F9F9", 0),
        rect(sx(0), 98, x1 - sx(0), 382, "#FDF7F9", 0),
        t((x0 + sx(0)) / 2, 118, "Repressed by GRBV beyond delay", 11.5, TEAL, "700", "middle"),
        t((sx(0) + x1) / 2, 118, "Induced by GRBV beyond delay", 11.5, BURG, "700", "middle"),
    ]
    for v in range(-4, 6):
        body.append(line(sx(v), 128, sx(v), 480, "#B8B0B5" if v == 0 else GRID, 2 if v == 0 else 1))
        body.append(t(sx(v), 500, f"{v:+d}" if v else "0", 11.5, MUTED, anchor="middle"))
    body.append(t((x0 + x1) / 2, 522, "Log2 fold change (infected vs healthy; *coupled within vines to berry viral RNA)", 12.2, MUTED, anchor="middle"))

    for i, (name, stage, a, b) in enumerate(rows):
        y = 148 + i * 38
        body += [
            t(38, y + 5, name, 13.6, INK, "700"),
            pill(292, y - 10, 106, 22, stage, SLATE_BG, SLATE, 10.8, "600"),
            line(sx(a), y, sx(b), y, "#C7BFC4", 3),
            dot(sx(a), y, 7, TEAL),
            dot(sx(b), y, 7, BURG),
        ]
    body += [
        dot(875, 45, 6.5, TEAL), t(888, 50, "2021", 13, INK, "600"),
        dot(955, 45, 6.5, BURG), t(968, 50, "2022", 13, INK, "600"),
    ]
    return card(
        1200, 542,
        "Major Host Factors &amp; Pathways Disrupted by GRBV",
        "Maturation-independent berry responses replicated across both seasons (dream + MPhS; P = 0.001 for sHSPs)",
        "".join(body),
    )


def viral_host_dynamics():
    body = []
    # Panel A: Petiole dPCR
    body.append(rect(38, 94, 360, 266, "#FAF9FA", 12, EDGE))
    body.append(t(56, 120, "A · Leaf-Petiole Viral DNA", 14.8, INK, "700"))
    body.append(t(56, 139, "Opposite seasonal slopes (P = 4.4e-13)", 11.8, MUTED))
    lx0, lx1, ly0, ly1 = 92, 366, 300, 166
    sx_l = lambda i: lx0 + (lx1 - lx0) * i / 3.0
    sy_l = lambda v: ly0 - (ly0 - ly1) * (v - 5.3) / (6.3 - 5.3)
    p21, p22 = [5.58, 5.81, 5.74, 6.21], [5.90, 5.50, 5.62, 5.37]
    for v in [5.4, 5.8, 6.2]:
        body += [line(lx0, sy_l(v), lx1, sy_l(v), GRID, 1), t(lx0 - 8, sy_l(v) + 4, f"{v:.1f}", 11, MUTED, anchor="end")]
    for i, s in enumerate(["Pre", "Ver", "Post", "Har"]):
        body.append(t(sx_l(i), 322, s, 11.5, MUTED, anchor="middle"))
    body += [poly([(sx_l(i), sy_l(v)) for i, v in enumerate(p21)], TEAL, 3),
             poly([(sx_l(i), sy_l(v)) for i, v in enumerate(p22)], BURG, 3)]
    for i in range(4):
        body += [dot(sx_l(i), sy_l(p21[i]), 5, TEAL), dot(sx_l(i), sy_l(p22[i]), 5, BURG)]
    body += [
        t(sx_l(3) - 6, sy_l(p21[3]) - 9, "2021 (+4.2x)", 11.5, TEAL, "700", "end"),
        t(sx_l(3) - 6, sy_l(p22[3]) - 10, "2022 (&#8722;3.3x)", 11.5, BURG, "700", "end"),
        t(218, 346, "log10 copies/&#181;L extracted DNA", 11.2, MUTED, anchor="middle"),
    ]

    # Panel B: Clock Lag vs Sugar Deficit
    body.append(rect(418, 94, 360, 266, "#FAF9FA", 12, EDGE))
    body.append(t(436, 120, "B · Delay vs. Must Sugar Loss", 14.8, INK, "700"))
    body.append(t(436, 139, "Healthy-vine clock (MPhS r = 0.986) &amp; Brix", 11.8, MUTED))
    body += [
        t(442, 180, "Harvest Maturation Delay (stage units)", 12, INK, "700"),
        t(442, 204, "2021", 12, TEAL, "700"), rect(484, 192, 210, 16, TEAL, 4), t(702, 205, "&#8722;0.811", 12, TEAL, "700"),
        t(442, 228, "2022", 12, BURG, "700"), rect(484, 216, 52, 16, BURG, 4), t(544, 229, "&#8722;0.202", 12, BURG, "700"),
        line(442, 248, 756, 248, GRID, 1.2),
        t(442, 272, "Harvest Sugar Deficit (&#176;Bx below healthy)", 12, INK, "700"),
        t(442, 296, "2021", 12, TEAL, "700"), rect(484, 284, 78, 16, TEAL, 4), t(570, 297, "&#8722;3.4 to &#8722;4.0", 12, TEAL, "700"),
        t(442, 320, "2022", 12, BURG, "700"), rect(484, 308, 205, 16, BURG, 4), t(696, 321, "&#8722;8.9 to &#8722;10.9", 12, BURG, "700"),
        t(598, 346, "Early foliar load &#8594; severe Brix deficit", 11.2, BURG, "600", "middle"),
    ]

    # Panel C: Berry Viral ORFs
    body.append(rect(798, 94, 364, 266, "#FAF9FA", 12, EDGE))
    body.append(t(816, 120, "C · Berry Viral ORF Shift", 14.8, INK, "700"))
    body.append(t(816, 139, "Complementary C2 rises as V1–V3 decline", 11.8, MUTED))
    rx0, rx1, ry0, ry1 = 850, 1134, 300, 168
    sx_r = lambda i: rx0 + (rx1 - rx0) * i / 3.0
    sy_r = lambda v: ry0 - (ry0 - ry1) * (v - (-3.3)) / (1.8 - (-3.3))
    c2, v2, v1, c1 = [0.0, 0.68, 1.44, 0.62], [0.0, -1.56, -2.14, -3.04], [0.0, -0.88, -1.32, -2.25], [0.0, -1.01, -1.03, -1.38]
    for v in [-3, -1.5, 0, 1.5]:
        body += [line(rx0, sy_r(v), rx1, sy_r(v), "#C9C1C5" if v == 0 else GRID, 1.5 if v == 0 else 1),
                 t(rx0 - 8, sy_r(v) + 4, f"{v:+.1f}" if v else "0", 10.5, MUTED, anchor="end")]
    for i, s in enumerate(["Pre", "Ver", "Post", "Har"]):
        body.append(t(sx_r(i), 322, s, 11.5, MUTED, anchor="middle"))
    body += [poly([(sx_r(i), sy_r(v)) for i, v in enumerate(v2)], "#7A8990", 2),
             poly([(sx_r(i), sy_r(v)) for i, v in enumerate(v1)], SLATE, 2),
             poly([(sx_r(i), sy_r(v)) for i, v in enumerate(c1)], GOLD, 2),
             poly([(sx_r(i), sy_r(v)) for i, v in enumerate(c2)], BURG, 3.4)]
    for i in range(4):
        body.append(dot(sx_r(i), sy_r(c2[i]), 5.5, BURG))
    body += [
        t(sx_r(2), sy_r(c2[2]) - 10, "C2 (+1.92, +0.95)", 11.5, BURG, "700", "middle"),
        t(sx_r(3) - 4, sy_r(v2[3]) - 7, "V1–V3, C1 drop", 11, "#526369", "700", "end"),
        t(980, 346, "681 host genes track berry RNA vs 0 petiole", 11.2, BURG, "600", "middle"),
    ]
    return card(
        1200, 382,
        "Viral Dynamics &amp; Ripening Symptoms",
        "Vegetative viral burden governs maturation delay (P = 0.0049), while berry C2 enrichment drives local host remodeling",
        "".join(body),
    )


def analysis_evolution():
    row1 = [
        ("Era 1", "Naive Screening", "Welch t-test &amp; WGCNA", "Exposed stage confounding", SLATE, SLATE_BG),
        ("Era 2", "Network &amp; Omics", "FARM Snakemake (00–12)", "160 berry NLRs stay quiet", SLATE, SLATE_BG),
        ("Era 3", "Stress-Testing", "Cutoff &amp; ML sweeps", "Retired viral AUC = 1.000", SLATE, SLATE_BG),
        ("Era 4", "Stage Covariates", "DESeq2 &amp; limma-voom", "908 pooled DEGs &#8594; 126", GOLD, GOLD_BG),
    ]
    row2 = [
        ("Era 5", "Orthology Check", "Pinot–Merlot MCScan", "87.5% 1:1 orthologs verified", TEAL, TEAL_BG),
        ("Era 6", "Mixed Models", "dream + claim gates", "44 vines as biological unit", TEAL, TEAL_BG),
        ("Era 7", "Ripening Clock", "Healthy-only staging", "Post-veraison viral C2 rise", TEAL, TEAL_BG),
        ("Era 8", "Pathology Synthesis", "MPhS + 267 Petiole dPCR", "Active pipeline in repo", BURG, BURG_BG),
    ]
    body = []
    w_box, h_box, gap = 256, 110, 32
    x_start = 38

    for i, (tag, title, method, outcome, fg, bg) in enumerate(row1):
        x = x_start + i * (w_box + gap)
        y = 100
        body += [
            rect(x, y, w_box, h_box, bg, 12, EDGE),
            pill(x + 14, y + 12, 60, 22, tag, fg, CARD, 11.2),
            t(x + 84, y + 28, title, 13.8, INK, "700"),
            t(x + 14, y + 58, method, 12.2, MUTED),
            t(x + 14, y + 84, outcome, 12.2, fg, "700"),
        ]
        if i < 3:
            body.append(arrow(x + w_box + 4, y + h_box / 2, x + w_box + gap - 4, SLATE, 2.2))

    body += [
        line(x_start + 3 * (w_box + gap) + w_box / 2, 214, x_start + 3 * (w_box + gap) + w_box / 2, 232, GOLD, 2.2),
        line(x_start + w_box / 2, 232, x_start + 3 * (w_box + gap) + w_box / 2, 232, GOLD, 2.2, "5,4"),
        v_arrow(x_start + w_box / 2, 232, 252, TEAL, 2.2),
    ]

    for i, (tag, title, method, outcome, fg, bg) in enumerate(row2):
        x = x_start + i * (w_box + gap)
        y = 256
        body += [
            rect(x, y, w_box, h_box, bg, 12, BURG if i == 3 else EDGE, 2.2 if i == 3 else 1.2),
            pill(x + 14, y + 12, 60, 22, tag, fg, CARD, 11.2),
            t(x + 84, y + 28, title, 13.8, INK, "700"),
            t(x + 14, y + 58, method, 12.2, MUTED),
            t(x + 14, y + 84, outcome, 12.2, fg, "700"),
        ]
        if i < 3:
            body.append(arrow(x + w_box + 4, y + h_box / 2, x + w_box + gap - 4, TEAL if i < 2 else BURG, 2.2))

    return card(
        1200, 390,
        "Analytical Evolution (Eras 1 &#8594; 8)",
        "How confounding artifacts were retired on the way to the final repeated-measures pathology pipeline",
        "".join(body),
    )


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("*.svg"):
        old.unlink()
    for fn in (graphical_abstract, study_overview, host_pathways, viral_host_dynamics, analysis_evolution):
        out_path = OUT / f"{fn.__name__}.svg"
        out_path.write_text(fn(), encoding="utf-8")
        print(f"Wrote {out_path.name}")
