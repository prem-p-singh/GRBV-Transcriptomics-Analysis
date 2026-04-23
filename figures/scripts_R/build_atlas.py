#!/usr/bin/env python
"""Concatenate the 10 Main composite PDFs into one figure atlas."""
from pathlib import Path
from pypdf import PdfWriter

HERE = Path(__file__).resolve().parent
FIG_DIR = Path("/Users/prempratapsingh/Data_Analysis/figures/main")
OUT = Path("/Users/prempratapsingh/Data_Analysis/figures/GRBV_Main_Figures_Atlas.pdf")

mains = [
    "Main01_overview", "Main02_de_landscape", "Main03_inverted_trajectory",
    "Main04_heat_stress", "Main05_mediation",  "Main06_biomarker",
    "Main07_wgcna_tf",   "Main08_chr_qtl",     "Main09_gsea",
    "Main10_biology_cascade",
]

writer = PdfWriter()
n_pages = 0
for stem in mains:
    p = FIG_DIR / f"{stem}.pdf"
    if not p.exists():
        raise SystemExit(f"missing: {p}")
    writer.append(str(p))
    n_pages += 1

with open(OUT, "wb") as fh:
    writer.write(fh)
print(f"wrote {OUT} ({n_pages} composites)")
