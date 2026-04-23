#!/usr/bin/env python3
"""Render EXTENDED_RESULTS_REPORT.md to a publication-style PDF using fpdf2."""
import re, os, html as htmllib
from pathlib import Path
import markdown
from fpdf import FPDF
from fpdf.enums import XPos, YPos

ROOT = Path("/Users/prempratapsingh/Data_Analysis")
SRC  = ROOT / "docs/EXTENDED_RESULTS_REPORT.md"
OUT  = ROOT / "docs/EXTENDED_RESULTS_REPORT.pdf"
RES  = ROOT / "results"

md_text = SRC.read_text(encoding="utf-8")
md_body = re.sub(r"^# .*?\n\n\*\*Study:.*?pipeline\)\n\n---\n\n", "", md_text, count=1, flags=re.DOTALL)

# Inject figure refs as markdown image tags where relevant
figure_map = [
    ("### 3.3 Rule 02",  "02_module_preservation/preservation_plot.png",
     "Module preservation (Zsummary vs module size); dashed lines mark Z=2 and Z=10 reference cut-offs."),
    ("### 3.4 Rule 03",  "03_subtype_nmf/cophenetic_plot.png",
     "Consensus NMF model-selection: cophenetic coefficient and Brunet dispersion across k."),
    ("### 3.4 Rule 03",  "03_subtype_nmf/consensus_matrix.png",
     "Consensus co-clustering matrix at the selected k=2."),
    ("### 3.7 Rule 06",  "06_pseudotime/trajectory_plot.png",
     "RB+ pseudotime trajectory on PC1-PC2, colored by slingshot pseudotime."),
    ("### 3.8 Rule 07",  "07_mofa/variance_explained.png",
     "MOFA+ per-view R^2 across 10 latent factors."),
    ("### 3.10 Rule 09", "09_defense/defense_heatmap.png",
     "Defense-category mean log2FC across season x stage combinations."),
    ("### 3.11 Rule 10", "10_subtype_characterize/subtype_overview.png",
     "Clinical overview of the two NMF subtypes (stage composition, season composition, YSI)."),
]

# fpdf2 write_html supports a bounded set of tags; build our own lightweight renderer.
# Easiest: parse MD -> HTML -> custom walk. But fpdf2 has basic write_html; give that a try
# with images replaced by figures inserted inline.

# Preprocess: strip characters fpdf's core fonts can't render (alpha/beta/Greek, smart quotes, em-dashes)
def sanitize(s):
    repl = {
        "—": "-", "–": "-", "−": "-",
        "’": "'", "‘": "'", "“": '"', "”": '"',
        "α": "alpha", "β": "beta", "η": "eta", "μ": "mu", "Δ": "Delta", "·": ".",
        "°": "deg", "²": "^2", "₂": "2", "↔": "<->", "→": "->", "←": "<-",
        "•": "-", "×": "x", "≤": "<=", "≥": ">=", "≠": "!=",
        "✓": "[ok]", "✗": "[x]", "—": "-",
        "\xa0": " ",
    }
    for k, v in repl.items():
        s = s.replace(k, v)
    # Drop any non-latin1 character remaining
    s = s.encode("latin-1", "replace").decode("latin-1")
    return s

# Parse MD -> HTML
html_body = markdown.markdown(
    md_body,
    extensions=["tables", "fenced_code", "sane_lists", "attr_list"],
    output_format="html5",
)

# --- Light HTML processing: insert figures right after anchor sections ---
anchor_patterns = []
for anchor, rel_path, caption in figure_map:
    img = RES / rel_path
    if not img.exists():
        continue
    section_id = anchor.replace("### ", "").strip()
    # Normalize anchor to match rendered heading (section_id contains text like "3.3 Rule 02 - ...")
    # We'll match on the raw heading text prefix.
    heading_key = section_id.split(" - ")[0]
    fig_html = (
        f'<figure><img src="{img}" width="520" /><br/>'
        f'<i>{sanitize(caption)}</i></figure>'
    )
    anchor_patterns.append((heading_key, fig_html))

for key, fig in anchor_patterns:
    # Find first <h3> whose text starts with this key and insert figure after its closing </h3>
    pattern = rf"(<h3[^>]*>\s*{re.escape(key)}[^<]*</h3>)"
    html_body = re.sub(pattern, r"\1" + fig, html_body, count=1)

# Sanitize text (leave HTML tags untouched by working tag-by-tag)
def sanitize_html(s):
    # Split by tags; sanitize text nodes only
    parts = re.split(r"(<[^>]+>)", s)
    out = []
    for p in parts:
        if p.startswith("<") and p.endswith(">"):
            out.append(p)
        else:
            out.append(sanitize(p))
    return "".join(out)

html_body = sanitize_html(html_body)
# fpdf2 write_html doesn't allow <strong>/<em> inside <td>; strip nested tags inside table cells.
html_body = re.sub(
    r"(<td[^>]*>)(.*?)(</td>)",
    lambda m: m.group(1) + re.sub(r"</?(strong|em|b|i|code)[^>]*>", "", m.group(2)) + m.group(3),
    html_body, flags=re.DOTALL,
)
html_body = re.sub(
    r"(<th[^>]*>)(.*?)(</th>)",
    lambda m: m.group(1) + re.sub(r"</?(strong|em|b|i|code)[^>]*>", "", m.group(2)) + m.group(3),
    html_body, flags=re.DOTALL,
)

# --- Build PDF ---
class Report(FPDF):
    def __init__(self):
        super().__init__(orientation="P", unit="mm", format="A4")
        self.set_auto_page_break(auto=True, margin=20)
        self.set_margins(left=18, top=20, right=18)

    def header(self):
        if self.page_no() > 1:
            self.set_font("Helvetica", "I", 8)
            self.set_text_color(120, 120, 120)
            self.cell(0, 5, "GRBV extended analysis - FARM pipeline",
                      new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="L")
            self.set_draw_color(180, 190, 200)
            self.line(18, 27, 192, 27)
            self.ln(3)

    def footer(self):
        self.set_y(-13)
        self.set_font("Helvetica", "I", 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 8, f"Page {self.page_no()} / {{nb}}", align="C")

    def title_page(self):
        self.add_page()
        self.ln(30)
        self.set_font("Helvetica", "B", 20)
        self.set_text_color(26, 57, 96)
        self.multi_cell(0, 11,
            "Extended Analysis Report\nFARM HPC Pipeline (Rules 00-10)", align="C")
        self.ln(6)
        self.set_font("Helvetica", "I", 13)
        self.set_text_color(58, 85, 125)
        self.multi_cell(0, 8,
            "GRBV Transcriptomics, Vitis vinifera cv. Merlot", align="C")
        self.ln(3)
        self.set_draw_color(26, 57, 96)
        self.set_line_width(0.4)
        self.line(60, self.get_y(), 150, self.get_y())
        self.ln(12)
        self.set_font("Helvetica", "", 11)
        self.set_text_color(60, 60, 60)
        self.multi_cell(0, 6,
            "Complementary analyses executed on top of the primary differential-expression study:\n"
            "gene regulatory network inference, WGCNA module preservation, consensus NMF subtype\n"
            "discovery, DGCA differential co-expression, pseudotime trajectory, MOFA+ multi-block\n"
            "integration, HIMA high-dimensional mediation, defense repertoire profiling,\n"
            "and subtype characterization.", align="C")
        self.ln(12)
        self.set_font("Helvetica", "B", 10)
        self.set_text_color(26, 57, 96)
        self.cell(0, 6, "Key numbers", align="C",
                  new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_font("Helvetica", "", 10)
        self.set_text_color(60, 60, 60)
        summary_lines = [
            "327 samples (264 RB+, 63 RB-)  |  14,604 retained genes  |  2 seasons  |  4 maturity stages",
            "12 of 13 rules completed  |  22 output files  |  ~17 MB of results",
            "GRN: 100,000 directed TF-target edges  |  Module preservation: all 6 modules Zsum >= 5.8",
            "HIMA: LOC132254135 identified as sole temperature-titer mediator (p = 2.8e-4)",
            "Defense repertoire: all 10 categories flat (|log2FC| < 0.08)",
        ]
        for line in summary_lines:
            self.cell(0, 5.5, sanitize(line), align="C",
                      new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(20)
        self.set_font("Helvetica", "B", 10)
        self.set_text_color(26, 57, 96)
        self.cell(0, 6, "Author, date, repository",
                  new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
        self.set_font("Helvetica", "", 10)
        self.set_text_color(60, 60, 60)
        self.cell(0, 6, "Prem Pratap Singh",
                  new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
        self.cell(0, 6, "2026-04-19",
                  new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
        self.cell(0, 6, "github.com/prem-p-singh/GRBV-Transcriptomics-Analysis",
                  new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")

pdf = Report()
pdf.add_page()                              # Dummy so footer {nb} numbering is set up
pdf.set_font("Helvetica", "", 10)
pdf = Report()
pdf.title_page()

# --- Main body ---
pdf.add_page()
pdf.set_font("Helvetica", "", 10)
# Tell fpdf to treat our CSS loosely
pdf.write_html(
    html_body,
    tag_styles={
        "h1": {"b": True, "size": 16, "color": (26, 57, 96)},
        "h2": {"b": True, "size": 13, "color": (37, 73, 122)},
        "h3": {"b": True, "size": 11, "color": (58, 87, 122)},
        "h4": {"b": True, "size": 10, "color": (70, 70, 90)},
    } if False else None,   # fpdf2 version-sensitive; keep default if this fails
)

pdf.output(str(OUT))
print(f"Wrote {OUT}  ({OUT.stat().st_size/1024:.1f} KB, {pdf.page_no()} pages)")
