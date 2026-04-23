"""Build a PDF of MANUSCRIPT_DRAFT_v1.md with all 10 Main figures embedded.

Output:
    /Users/prempratapsingh/Data_Analysis/docs/MANUSCRIPT_DRAFT_v1.pdf

Figure placement:
    After each Results subsection that references a figure, the corresponding
    Main0X_*.png is embedded and shown at page width. Supplementary Main figures
    not cited inline (Main09 GSEA, Main10 biology cascade) appear at the end.
"""

from __future__ import annotations
import os
import re
from pathlib import Path
import matplotlib
from fpdf import FPDF

ROOT   = Path("/Users/prempratapsingh/Data_Analysis")
MD     = ROOT / "docs" / "MANUSCRIPT_DRAFT_v1.md"
FIGDIR = ROOT / "figures" / "main"
OUT    = ROOT / "docs" / "MANUSCRIPT_DRAFT_v1.pdf"

# Unicode-capable font (DejaVu Sans ships with matplotlib)
_MPL_FONTS = Path(os.path.dirname(matplotlib.__file__)) / "mpl-data" / "fonts" / "ttf"
FONT_REGULAR = _MPL_FONTS / "DejaVuSans.ttf"
FONT_BOLD    = _MPL_FONTS / "DejaVuSans-Bold.ttf"
FONT_ITALIC  = _MPL_FONTS / "DejaVuSans-Oblique.ttf"
FONT_BI      = _MPL_FONTS / "DejaVuSans-BoldOblique.ttf"
FAMILY = "DejaVu"

# In-text figure insertion: after which Results ### heading to embed which Main figure.
# Keys are the exact ### subheading text (after "### ").
RESULTS_FIG_MAP = {
    "Global variance structure motivates covariate-adjusted differential expression":
        ("Main01_overview.png", "Figure 1"),
    "A compact pooled infection signature amplifies through fruit maturity":
        ("Main02_de_landscape.png", "Figure 2"),
    "The HSF to sHSP thermal-defence axis is dampened but topologically preserved":
        ("Main04_heat_stress.png", "Figure 4"),
    "Co-expression network preservation and DOF transcription-factor motif enrichment":
        ("Main07_wgcna_tf.png", "Figure 7"),
    "Biphasic temperature-titer relationship and multiple-mediator analysis":
        ("Main05_mediation.png", "Figure 5"),
    "Host-only biomarker panel validates with leave-one-season-out cross-validation":
        ("Main06_biomarker.png", "Figure 6"),
}

# Extra figures shown at the end under a "Supporting main figures" section.
TAIL_FIGS = [
    ("Main03_inverted_trajectory.png", "Figure 3",
     "Inverted seasonal trajectory of differential expression and biphasic viral "
     "titer dynamics across stages."),
    ("Main08_chr_qtl.png", "Figure 8",
     "Genomic localisation of host biomarkers and biomarker-to-QTL overlap versus "
     "a chromosome-stratified permutation null."),
    ("Main09_gsea.png", "Figure 9",
     "GSEA across stage and season conditions; NES heatmap, suppressed pathway "
     "ranking, and cross-season NES agreement per stage."),
    ("Main10_biology_cascade.png", "Figure 10",
     "Integrated biology cascade summarising the five-step GRBV disruption model: "
     "establishment, silent immunity, thermal dampening, sugar-pathway rewiring, "
     "and disease phenotype."),
]

# ---------------------------------------------------------------------------
# Parse markdown into structured blocks
# ---------------------------------------------------------------------------

def parse_blocks(md_text: str) -> list[tuple[str, str]]:
    """Return a flat list of (block_kind, content) tuples.

    Kinds: 'h1', 'h2', 'h3', 'para', 'hr', 'table', 'em_line'.
    Tables are stored as a single string with newlines preserved.
    """
    blocks: list[tuple[str, str]] = []
    lines = md_text.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        stripped = line.rstrip()
        if not stripped:
            i += 1
            continue
        if stripped == "---":
            blocks.append(("hr", ""))
            i += 1
            continue
        if stripped.startswith("# "):
            blocks.append(("h1", stripped[2:].strip()))
            i += 1
            continue
        if stripped.startswith("## "):
            blocks.append(("h2", stripped[3:].strip()))
            i += 1
            continue
        if stripped.startswith("### "):
            blocks.append(("h3", stripped[4:].strip()))
            i += 1
            continue
        # Italic-line (single-line emphasised blocks used for subtitles)
        if stripped.startswith("*") and stripped.endswith("*") and stripped.count("*") == 2:
            blocks.append(("em_line", stripped.strip("*").strip()))
            i += 1
            continue
        # Tables: pipe-delimited, consecutive lines starting with |
        if stripped.startswith("|"):
            tbl = [stripped]
            i += 1
            while i < n and lines[i].lstrip().startswith("|"):
                tbl.append(lines[i].rstrip())
                i += 1
            blocks.append(("table", "\n".join(tbl)))
            continue
        # Ordered list item: collect consecutive items as one para block with
        # leading markers preserved.
        if re.match(r"^\s*\d+\.\s", stripped):
            items = [stripped]
            i += 1
            while i < n and (re.match(r"^\s*\d+\.\s", lines[i].lstrip())
                             or (lines[i].startswith("   ") and lines[i].strip())):
                items.append(lines[i].rstrip())
                i += 1
            blocks.append(("para", "\n".join(items)))
            continue
        # Default: accumulate consecutive non-blank, non-special lines as one paragraph
        para = [stripped]
        i += 1
        while i < n:
            nxt = lines[i].rstrip()
            if not nxt:
                break
            if nxt.startswith("#") or nxt.startswith("|") or nxt == "---":
                break
            para.append(nxt)
            i += 1
        blocks.append(("para", " ".join(para)))
    return blocks


# ---------------------------------------------------------------------------
# Inline markdown simplification for fpdf2
# ---------------------------------------------------------------------------

def md_inline_to_text(s: str) -> str:
    """Strip inline markdown that fpdf2 won't honour; keep unicode intact.

    DejaVu Sans supports Greek letters, en/em dashes, arrows, multiplication
    sign, combining diacritics, etc., so no character substitutions needed.
    """
    s = re.sub(r"\*\*(.+?)\*\*", r"\1", s)     # **bold**
    s = re.sub(r"\*(.+?)\*", r"\1", s)         # *italic*
    s = re.sub(r"`([^`]+)`", r"\1", s)         # `code`
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", s)  # [link](url)
    return s


# ---------------------------------------------------------------------------
# Build the PDF
# ---------------------------------------------------------------------------

class ManuscriptPDF(FPDF):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_font(FAMILY, "",  str(FONT_REGULAR))
        self.add_font(FAMILY, "B", str(FONT_BOLD))
        self.add_font(FAMILY, "I", str(FONT_ITALIC))
        self.add_font(FAMILY, "BI", str(FONT_BI))

    def header(self):
        if self.page_no() == 1:
            return
        self.set_font(FAMILY, "I", 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 6, "GRBV manuscript draft v1", align="L")
        self.cell(0, 6, f"Page {self.page_no()}", align="R")
        self.ln(8)
        self.set_text_color(0, 0, 0)


def render_table(pdf: ManuscriptPDF, md_table: str):
    rows = [r for r in md_table.splitlines() if r.strip()]
    if len(rows) < 2:
        return
    parsed = []
    for r in rows:
        cells = [c.strip() for c in r.strip().strip("|").split("|")]
        parsed.append(cells)
    # Drop the header-separator row (all dashes / pipes)
    parsed = [row for row in parsed
              if not all(re.match(r"^:?-+:?$", c) for c in row if c)]
    if not parsed:
        return
    ncols = max(len(r) for r in parsed)
    avail = pdf.w - pdf.l_margin - pdf.r_margin
    col_w = avail / ncols
    pdf.set_font(FAMILY, "B", 9)
    # header
    for c in parsed[0]:
        pdf.cell(col_w, 6, md_inline_to_text(c), border=1)
    pdf.ln()
    pdf.set_font(FAMILY, "", 9)
    for row in parsed[1:]:
        # pad row to ncols
        row = row + [""] * (ncols - len(row))
        # compute row height by longest cell
        heights = []
        for c in row:
            txt = md_inline_to_text(c)
            lines = pdf.multi_cell(col_w, 5, txt, split_only=True)
            heights.append(len(lines) * 5)
        rh = max(heights) if heights else 5
        # check page break
        if pdf.get_y() + rh > pdf.h - pdf.b_margin:
            pdf.add_page()
        x0, y0 = pdf.get_x(), pdf.get_y()
        for i, c in enumerate(row):
            txt = md_inline_to_text(c)
            pdf.set_xy(x0 + i * col_w, y0)
            pdf.multi_cell(col_w, 5, txt, border=1)
        pdf.set_xy(x0, y0 + rh)
    pdf.ln(2)


def insert_figure(pdf: ManuscriptPDF, img_path: Path, fig_label: str,
                  caption: str | None = None):
    """Insert an image sized to the text column, with a bold figure label."""
    if not img_path.exists():
        pdf.set_font(FAMILY, "I", 9)
        pdf.set_text_color(180, 30, 30)
        pdf.multi_cell(0, 5, f"[MISSING FIGURE: {img_path.name}]")
        pdf.set_text_color(0, 0, 0)
        return

    # Compute image width
    avail = pdf.w - pdf.l_margin - pdf.r_margin
    img_w = avail

    # Compute image display height. fpdf2 exposes image intrinsic size via
    # a reader; we rely on simple aspect handling by letting fpdf scale.
    # We set h=0 to preserve aspect.
    # First check if we need a new page for the image plus caption.
    est_height = 100  # conservative estimate; fpdf2 will render natural size up to max
    if pdf.get_y() + est_height > pdf.h - pdf.b_margin:
        pdf.add_page()

    pdf.set_font(FAMILY, "B", 10)
    pdf.ln(3)
    pdf.cell(0, 6, fig_label, ln=1)
    pdf.image(str(img_path), x=pdf.l_margin, w=img_w)
    if caption:
        pdf.set_font(FAMILY, "", 8.5)
        pdf.ln(1)
        pdf.multi_cell(0, 4.5, md_inline_to_text(caption))
    pdf.ln(3)


def render_para(pdf: ManuscriptPDF, text: str):
    pdf.set_font(FAMILY, "", 10.5)
    pdf.multi_cell(0, 5.6, md_inline_to_text(text))
    pdf.ln(1.5)


def render_em_line(pdf: ManuscriptPDF, text: str):
    pdf.set_font(FAMILY, "I", 9.5)
    pdf.set_text_color(80, 80, 80)
    pdf.multi_cell(0, 5, md_inline_to_text(text))
    pdf.set_text_color(0, 0, 0)
    pdf.ln(1)


def render_h1(pdf: ManuscriptPDF, text: str):
    pdf.set_font(FAMILY, "B", 15)
    pdf.multi_cell(0, 8, md_inline_to_text(text))
    pdf.ln(2)


def render_h2(pdf: ManuscriptPDF, text: str):
    pdf.ln(2)
    pdf.set_font(FAMILY, "B", 13)
    pdf.multi_cell(0, 7, md_inline_to_text(text))
    pdf.ln(1)


def render_h3(pdf: ManuscriptPDF, text: str):
    pdf.ln(1)
    pdf.set_font(FAMILY, "B", 11)
    pdf.multi_cell(0, 6, md_inline_to_text(text))
    pdf.ln(0.5)


def build():
    md_text = MD.read_text()
    blocks = parse_blocks(md_text)

    pdf = ManuscriptPDF(unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.set_left_margin(18)
    pdf.set_right_margin(18)
    pdf.set_top_margin(15)
    pdf.add_page()

    # Walk the blocks and emit PDF, inserting figures after matching ### blocks.
    inside_results = False
    for kind, content in blocks:
        if kind == "h1":
            render_h1(pdf, content)
        elif kind == "h2":
            inside_results = (content == "Results")
            render_h2(pdf, content)
        elif kind == "h3":
            render_h3(pdf, content)
            if inside_results and content in RESULTS_FIG_MAP:
                # Peek ahead: we want to render the following paragraphs first,
                # then the figure. We'll defer figure insertion until we see the
                # next h3 or h2 or hr. Simplest: stash the figure and emit after
                # the next blank/hr. Implementation: use a lookup to insert
                # *after* this subsection's paragraphs. We achieve that by
                # consuming subsequent blocks until the next heading/hr.
                pass
        elif kind == "hr":
            pdf.ln(2)
            y = pdf.get_y()
            pdf.set_draw_color(200, 200, 200)
            pdf.line(pdf.l_margin, y, pdf.w - pdf.r_margin, y)
            pdf.set_draw_color(0, 0, 0)
            pdf.ln(2)
        elif kind == "em_line":
            render_em_line(pdf, content)
        elif kind == "table":
            render_table(pdf, content)
        elif kind == "para":
            render_para(pdf, content)

    # Simplified strategy wasn't sufficient for in-place figure insertion —
    # the first pass above built the text-only body. We now rebuild with a
    # state machine that defers figure insertion until the end of each mapped
    # subsection (marked by the next h3, h2, or hr).
    pass  # placeholder; real work done below


def build_v2():
    md_text = MD.read_text()
    blocks = parse_blocks(md_text)

    pdf = ManuscriptPDF(unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.set_left_margin(18)
    pdf.set_right_margin(18)
    pdf.set_top_margin(15)
    pdf.add_page()

    inside_results = False
    pending_fig: tuple[str, str] | None = None

    def flush_pending_figure():
        nonlocal pending_fig
        if pending_fig is None:
            return
        fname, label = pending_fig
        insert_figure(pdf, FIGDIR / fname, label)
        pending_fig = None

    for kind, content in blocks:
        if kind == "h1":
            flush_pending_figure()
            render_h1(pdf, content)
        elif kind == "h2":
            flush_pending_figure()
            inside_results = (content == "Results")
            render_h2(pdf, content)
        elif kind == "h3":
            flush_pending_figure()
            render_h3(pdf, content)
            if inside_results and content in RESULTS_FIG_MAP:
                pending_fig = RESULTS_FIG_MAP[content]
        elif kind == "hr":
            flush_pending_figure()
            pdf.ln(2)
            y = pdf.get_y()
            pdf.set_draw_color(200, 200, 200)
            pdf.line(pdf.l_margin, y, pdf.w - pdf.r_margin, y)
            pdf.set_draw_color(0, 0, 0)
            pdf.ln(2)
        elif kind == "em_line":
            render_em_line(pdf, content)
        elif kind == "table":
            render_table(pdf, content)
        elif kind == "para":
            render_para(pdf, content)

    # Flush any trailing figure
    flush_pending_figure()

    # Supporting main figures at the end
    pdf.add_page()
    render_h2(pdf, "Supporting main figures")
    for fname, label, caption in TAIL_FIGS:
        insert_figure(pdf, FIGDIR / fname, label, caption)

    pdf.output(str(OUT))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    build_v2()
