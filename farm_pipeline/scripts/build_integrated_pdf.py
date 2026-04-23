#!/usr/bin/env python3
"""Render INTEGRATED_ANALYSIS.md to PDF using reportlab Platypus.

Tables are rendered as true reportlab Table objects with column-width
management; paragraphs use Paragraph with basic markup translation.
"""
from __future__ import annotations
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether
)

ROOT = Path("/Users/prempratapsingh/Data_Analysis")
SRC  = ROOT / "docs/INTEGRATED_ANALYSIS.md"
OUT  = ROOT / "docs/INTEGRATED_ANALYSIS.pdf"

md = SRC.read_text(encoding="utf-8")
# Strip the H1 + italic lead (we use a title page instead)
md = re.sub(r"^# .*?\n\n\*.*?thread\.\*\n\n", "", md, count=1, flags=re.DOTALL)
# Drop horizontal rules
md = re.sub(r"^\s*---\s*$", "", md, flags=re.MULTILINE)

# --- Markdown token stream ---------------------------------------------------

def md_inline_to_rl(s: str) -> str:
    """Convert basic inline Markdown to reportlab Paragraph markup."""
    # Escape angle brackets
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    # bold
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    # italic (single * or _)
    s = re.sub(r"(?<!\w)\*(?!\s)(.+?)(?<!\s)\*(?!\w)", r"<i>\1</i>", s)
    s = re.sub(r"(?<!\w)_(.+?)_(?!\w)", r"<i>\1</i>", s)
    # inline code
    s = re.sub(r"`([^`]+)`", r'<font face="Courier" color="#222244">\1</font>', s)
    # markdown links [text](url) -> anchor (reportlab supports <a href>)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2" color="#1a4d8c">\1</a>', s)
    # common unicode substitutions
    repl = {"—": "—", "–": "–", "α": "α", "β": "β", "η": "η", "μ": "μ",
            "²": "²", "³": "³", "Δ": "Δ", "°": "°", "×": "×", "·": "·",
            "≤": "≤", "≥": "≥", "→": "→", "←": "←"}
    for k, v in repl.items():
        s = s.replace(k, v)
    return s


def parse_table(lines: list[str]):
    """Accept a list of table lines (pipe-separated) and return rows."""
    rows = []
    for ln in lines:
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        rows.append(cells)
    # Drop the separator row (the one with ---)
    rows = [r for r in rows if not all(set(c) <= set("-: ") for c in r)]
    return rows


# --- reportlab styles -------------------------------------------------------

styles = getSampleStyleSheet()
body = ParagraphStyle(
    "body", parent=styles["BodyText"], fontName="Times-Roman", fontSize=10,
    leading=13, alignment=4, spaceAfter=3, spaceBefore=2,
)
h1 = ParagraphStyle(
    "h1", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=17,
    leading=20, textColor=colors.HexColor("#1a3960"), spaceBefore=14, spaceAfter=6,
)
h2 = ParagraphStyle(
    "h2", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=13,
    leading=16, textColor=colors.HexColor("#25497a"), spaceBefore=12, spaceAfter=5,
    borderPadding=(0, 0, 2, 0),
)
h3 = ParagraphStyle(
    "h3", parent=styles["Heading3"], fontName="Helvetica-Bold", fontSize=11.5,
    leading=14, textColor=colors.HexColor("#3a577a"), spaceBefore=9, spaceAfter=3,
)
li = ParagraphStyle(
    "li", parent=body, leftIndent=14, bulletIndent=4, spaceAfter=1, spaceBefore=1,
)
code = ParagraphStyle(
    "code", parent=body, fontName="Courier", fontSize=8.5, leading=11,
    backColor=colors.HexColor("#f4f5f7"), borderPadding=3, textColor=colors.HexColor("#222"),
    leftIndent=8, rightIndent=8, spaceBefore=4, spaceAfter=4, alignment=0,
)
title_style = ParagraphStyle(
    "title", fontName="Helvetica-Bold", fontSize=22, leading=26,
    textColor=colors.HexColor("#122a4d"), alignment=1, spaceAfter=4,
)
subtitle_style = ParagraphStyle(
    "subtitle", fontName="Times-Italic", fontSize=13, leading=16,
    textColor=colors.HexColor("#3a557d"), alignment=1, spaceAfter=12,
)
meta_style = ParagraphStyle(
    "meta", fontName="Helvetica", fontSize=9, leading=12,
    textColor=colors.HexColor("#444"), alignment=1,
)

table_base_style = TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef2f7")),
    ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#1a3960")),
    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
    ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
    ("FONTSIZE", (0, 0), (-1, -1), 8.5),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#c7d1de")),
    ("ROWBACKGROUNDS", (0, 1), (-1, -1),
     [colors.white, colors.HexColor("#fafbfc")]),
    ("LEFTPADDING", (0, 0), (-1, -1), 4),
    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ("TOPPADDING", (0, 0), (-1, -1), 3),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
])


# --- Render --------------------------------------------------------------

flow: list = []

# Title block on first page
flow.append(Spacer(1, 18 * mm))
flow.append(Paragraph("Integrated analysis", title_style))
flow.append(Paragraph("GRBV / <i>Vitis vinifera</i> cv. Merlot", subtitle_style))
flow.append(Paragraph(
    "All findings from the main transcriptomic paper and the extended "
    "FARM pipeline woven into a single, internally consistent account "
    "of what the data actually support.", meta_style))
flow.append(Spacer(1, 8 * mm))
flow.append(Paragraph(
    "<b>Author:</b> Prem Pratap Singh &nbsp; &bull; &nbsp; "
    "<b>Date:</b> 2026-04-20 &nbsp; &bull; &nbsp; "
    '<b>Repository:</b> <a href="https://github.com/prem-p-singh/GRBV-Transcriptomics-Analysis" '
    'color="#1a4d8c">github.com/prem-p-singh/GRBV-Transcriptomics-Analysis</a>',
    meta_style))
flow.append(PageBreak())

# Walk the markdown line by line
lines = md.splitlines()
i = 0
para_buf: list[str] = []

def flush_para():
    global para_buf
    if para_buf:
        text = " ".join(para_buf).strip()
        if text:
            flow.append(Paragraph(md_inline_to_rl(text), body))
        para_buf = []

while i < len(lines):
    ln = lines[i]
    # Headings
    m = re.match(r"^(#{1,3})\s+(.+)$", ln)
    if m:
        flush_para()
        level = len(m.group(1))
        txt = md_inline_to_rl(m.group(2))
        flow.append(Paragraph(txt, [h1, h2, h3][level - 1]))
        i += 1; continue
    # Fenced code block
    if ln.startswith("```"):
        flush_para()
        i += 1
        buf = []
        while i < len(lines) and not lines[i].startswith("```"):
            buf.append(lines[i]); i += 1
        i += 1
        code_text = "\n".join(buf)
        # reportlab Preformatted-style via a single Paragraph with Courier
        esc = code_text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        flow.append(Paragraph("<br/>".join(esc.split("\n")), code))
        continue
    # Table detection: pipe-separated line followed by a separator row
    if ln.startswith("|") and i + 1 < len(lines) and re.match(r"^\|[\s|:\-]+\|$", lines[i + 1].strip()):
        flush_para()
        tbl = [ln]
        j = i + 1
        while j < len(lines) and lines[j].startswith("|"):
            tbl.append(lines[j]); j += 1
        rows = parse_table(tbl)
        # Convert cells to Paragraphs to let text wrap
        rl_rows = [[Paragraph(md_inline_to_rl(c), body) for c in r] for r in rows]
        ncols = max(len(r) for r in rl_rows)
        usable_width = 175 * mm  # A4 portrait minus margins
        col_widths = [usable_width / ncols] * ncols
        tbl_obj = Table(rl_rows, colWidths=col_widths, repeatRows=1)
        tbl_obj.setStyle(table_base_style)
        flow.append(Spacer(1, 2 * mm))
        flow.append(tbl_obj)
        flow.append(Spacer(1, 3 * mm))
        i = j; continue
    # Bulleted list
    m_ul = re.match(r"^(\s*)[\-*]\s+(.+)$", ln)
    if m_ul:
        flush_para()
        flow.append(Paragraph(
            "&bull;&nbsp;" + md_inline_to_rl(m_ul.group(2)), li))
        i += 1; continue
    # Numbered list
    m_ol = re.match(r"^(\s*)(\d+)\.\s+(.+)$", ln)
    if m_ol:
        flush_para()
        flow.append(Paragraph(
            f"{m_ol.group(2)}.&nbsp; " + md_inline_to_rl(m_ol.group(3)), li))
        i += 1; continue
    # Blank line -> paragraph break
    if not ln.strip():
        flush_para(); i += 1; continue
    # Default: accumulate into paragraph buffer
    para_buf.append(ln.strip()); i += 1

flush_para()


# --- Build the document ---------------------------------------------------

def _header_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica-Oblique", 8)
    canvas.setFillColorRGB(0.45, 0.45, 0.45)
    if doc.page > 1:
        canvas.drawString(18 * mm, A4[1] - 14 * mm, "GRBV integrated analysis")
    canvas.drawCentredString(A4[0] / 2, 10 * mm, f"Page {doc.page}")
    canvas.restoreState()


doc = SimpleDocTemplate(
    str(OUT), pagesize=A4,
    leftMargin=18 * mm, rightMargin=18 * mm,
    topMargin=18 * mm, bottomMargin=15 * mm,
    title="GRBV integrated analysis",
    author="Prem Pratap Singh",
)
doc.build(flow, onLaterPages=_header_footer, onFirstPage=_header_footer)

print(f"Wrote {OUT}  ({OUT.stat().st_size/1024:.1f} KB)")
