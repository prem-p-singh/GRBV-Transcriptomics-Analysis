#!/usr/bin/env python
"""Assemble the 10 Main composite PNGs into a single navigable PDF atlas.

Since the figures are now written only as PNG + SVG (no PDF), the atlas
is rebuilt directly from the PNG rasters using reportlab. Each page fits
one Main figure to page width, preserving aspect ratio.
"""
from pathlib import Path
from PIL import Image
from reportlab.lib.pagesizes import LETTER
from reportlab.pdfgen import canvas

FIG_DIR = Path("/Users/prempratapsingh/Data_Analysis/figures/main")
OUT     = Path("/Users/prempratapsingh/Data_Analysis/figures/GRBV_Main_Figures_Atlas.pdf")

MAINS = [
    "Main01_overview", "Main02_de_landscape", "Main03_inverted_trajectory",
    "Main04_heat_stress", "Main05_mediation",  "Main06_biomarker",
    "Main07_wgcna_tf",   "Main08_chr_qtl",     "Main09_gsea",
    "Main10_biology_cascade",
]

MARGIN_IN = 0.4 * 72    # 0.4 inch margin
page_w, page_h = LETTER

c = canvas.Canvas(str(OUT), pagesize=LETTER)
for stem in MAINS:
    png = FIG_DIR / f"{stem}.png"
    if not png.exists():
        raise SystemExit(f"missing: {png}")
    img = Image.open(png)
    w, h = img.size
    usable_w = page_w - 2 * MARGIN_IN
    usable_h = page_h - 2 * MARGIN_IN - 20   # leave room for caption
    scale = min(usable_w / w, usable_h / h)
    draw_w = w * scale
    draw_h = h * scale
    x = (page_w - draw_w) / 2
    y = (page_h - draw_h) / 2 + 10
    c.drawImage(str(png), x, y, draw_w, draw_h,
                preserveAspectRatio=True, anchor="c")
    c.setFont("Helvetica", 9)
    c.setFillGray(0.4)
    c.drawString(MARGIN_IN, MARGIN_IN, stem.replace("_", " "))
    c.drawRightString(page_w - MARGIN_IN, MARGIN_IN,
                      f"Page {MAINS.index(stem) + 1} / {len(MAINS)}")
    c.showPage()
c.save()
print(f"wrote {OUT} ({len(MAINS)} composites)")
