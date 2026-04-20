#!/usr/bin/env python3
"""Render SENSITIVITY_REPORT.md to PDF using reportlab (same engine as INTEGRATED)."""
import re, sys
from pathlib import Path
sys.path.insert(0, "/Users/prempratapsingh/Data_Analysis/farm_pipeline/scripts")
# Reuse the renderer from build_integrated_pdf.py
import importlib.util
spec = importlib.util.spec_from_file_location(
    "rl_renderer",
    "/Users/prempratapsingh/Data_Analysis/farm_pipeline/scripts/build_integrated_pdf.py",
)
# That script's output path is hardcoded — we re-implement minimally here.

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak

ROOT = Path("/Users/prempratapsingh/Data_Analysis/farm_pipeline")
SRC  = ROOT / "docs/SENSITIVITY_REPORT.md"
OUT  = ROOT / "docs/SENSITIVITY_REPORT.pdf"
md = SRC.read_text(encoding="utf-8")
md = re.sub(r"^# .*?\n\n.*?ledger\. Generated automatically from\n`results/sensitivity/` outputs\.\n\n", "", md, count=1, flags=re.DOTALL)

def md_inline(s):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"(?<!\w)\*(?!\s)(.+?)(?<!\s)\*(?!\w)", r"<i>\1</i>", s)
    s = re.sub(r"`([^`]+)`", r'<font face="Courier" color="#222244">\1</font>', s)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2" color="#1a4d8c">\1</a>', s)
    return s

body = ParagraphStyle("body", fontName="Times-Roman", fontSize=10, leading=13,
                     alignment=4, spaceAfter=3, spaceBefore=2)
h1 = ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=17, leading=20,
                    textColor=colors.HexColor("#1a3960"), spaceBefore=14, spaceAfter=6)
h2 = ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=13, leading=16,
                    textColor=colors.HexColor("#25497a"), spaceBefore=12, spaceAfter=5)
h3 = ParagraphStyle("h3", fontName="Helvetica-Bold", fontSize=11.5, leading=14,
                    textColor=colors.HexColor("#3a577a"), spaceBefore=9, spaceAfter=3)
li = ParagraphStyle("li", parent=body, leftIndent=14, bulletIndent=4, spaceAfter=1, spaceBefore=1)
code = ParagraphStyle("code", parent=body, fontName="Courier", fontSize=8.5, leading=11,
                      backColor=colors.HexColor("#f4f5f7"), borderPadding=3,
                      leftIndent=8, rightIndent=8, spaceBefore=4, spaceAfter=4, alignment=0)
tbl_style = TableStyle([
    ("BACKGROUND",(0,0),(-1,0), colors.HexColor("#eef2f7")),
    ("TEXTCOLOR",(0,0),(-1,0), colors.HexColor("#1a3960")),
    ("FONTNAME",(0,0),(-1,0), "Helvetica-Bold"),
    ("FONTNAME",(0,1),(-1,-1), "Helvetica"),
    ("FONTSIZE",(0,0),(-1,-1), 8.5),
    ("VALIGN",(0,0),(-1,-1), "TOP"),
    ("GRID",(0,0),(-1,-1), 0.25, colors.HexColor("#c7d1de")),
    ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white, colors.HexColor("#fafbfc")]),
    ("LEFTPADDING",(0,0),(-1,-1), 4),("RIGHTPADDING",(0,0),(-1,-1), 4),
    ("TOPPADDING",(0,0),(-1,-1), 3),("BOTTOMPADDING",(0,0),(-1,-1), 3),
])

flow = [
    Spacer(1, 18*mm),
    Paragraph("Sensitivity analysis — all variance-filter bottlenecks",
              ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=20, leading=24,
                             textColor=colors.HexColor("#122a4d"), alignment=1, spaceAfter=4)),
    Paragraph("GRBV / <i>Vitis vinifera</i> cv. Merlot", ParagraphStyle(
        "sub", fontName="Times-Italic", fontSize=13, leading=16,
        textColor=colors.HexColor("#3a557d"), alignment=1, spaceAfter=12)),
    Paragraph("Every analysis that applies a top-N variable-gene filter was re-run at 3–5 additional "
              "thresholds. Results aggregated into a single robustness ledger with recommended paper edits.",
              ParagraphStyle("abstract", fontName="Times-Roman", fontSize=10.5, leading=14,
                             alignment=1, textColor=colors.HexColor("#333"), spaceAfter=8)),
    Paragraph("<b>Author:</b> Prem Pratap Singh &nbsp;•&nbsp; <b>Date:</b> 2026-04-20",
              ParagraphStyle("meta", fontName="Helvetica", fontSize=9,
                             alignment=1, textColor=colors.HexColor("#555"))),
    PageBreak(),
]

def parse_table(lines):
    rows = []
    for ln in lines:
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        rows.append(cells)
    rows = [r for r in rows if not all(set(c) <= set("-: ") for c in r)]
    return rows

lines = md.splitlines()
i = 0; para = []
def flush():
    global para
    if para:
        t = " ".join(para).strip()
        if t:
            flow.append(Paragraph(md_inline(t), body))
        para = []

while i < len(lines):
    ln = lines[i]
    m = re.match(r"^(#{1,3})\s+(.+)$", ln)
    if m:
        flush()
        flow.append(Paragraph(md_inline(m.group(2)), [h1, h2, h3][len(m.group(1))-1]))
        i += 1; continue
    if ln.startswith("```"):
        flush(); i += 1
        buf = []
        while i < len(lines) and not lines[i].startswith("```"):
            buf.append(lines[i]); i += 1
        i += 1
        esc = "\n".join(buf).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        flow.append(Paragraph("<br/>".join(esc.split("\n")), code))
        continue
    if ln.startswith("|") and i+1 < len(lines) and re.match(r"^\|[\s|:\-]+\|$", lines[i+1].strip()):
        flush()
        tbl = [ln]; j = i+1
        while j < len(lines) and lines[j].startswith("|"):
            tbl.append(lines[j]); j += 1
        rows = parse_table(tbl)
        rl_rows = [[Paragraph(md_inline(c), body) for c in r] for r in rows]
        ncols = max(len(r) for r in rl_rows)
        tbl_obj = Table(rl_rows, colWidths=[175*mm/ncols]*ncols, repeatRows=1)
        tbl_obj.setStyle(tbl_style)
        flow.append(Spacer(1, 2*mm)); flow.append(tbl_obj); flow.append(Spacer(1, 3*mm))
        i = j; continue
    if re.match(r"^\s*[-*]\s+", ln):
        flush()
        flow.append(Paragraph("&bull;&nbsp;" + md_inline(re.sub(r"^\s*[-*]\s+", "", ln)), li))
        i += 1; continue
    m_ol = re.match(r"^(\s*)(\d+)\.\s+(.+)$", ln)
    if m_ol:
        flush()
        flow.append(Paragraph(f"{m_ol.group(2)}.&nbsp; " + md_inline(m_ol.group(3)), li))
        i += 1; continue
    if not ln.strip():
        flush(); i += 1; continue
    para.append(ln.strip()); i += 1
flush()

def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica-Oblique", 8)
    canvas.setFillColorRGB(0.45,0.45,0.45)
    if doc.page > 1:
        canvas.drawString(18*mm, A4[1]-14*mm, "GRBV sensitivity analysis")
    canvas.drawCentredString(A4[0]/2, 10*mm, f"Page {doc.page}")
    canvas.restoreState()

doc = SimpleDocTemplate(str(OUT), pagesize=A4,
                        leftMargin=18*mm, rightMargin=18*mm,
                        topMargin=18*mm, bottomMargin=15*mm,
                        title="GRBV sensitivity analysis")
doc.build(flow, onLaterPages=footer, onFirstPage=footer)
print(f"Wrote {OUT}  ({OUT.stat().st_size/1024:.1f} KB)")
