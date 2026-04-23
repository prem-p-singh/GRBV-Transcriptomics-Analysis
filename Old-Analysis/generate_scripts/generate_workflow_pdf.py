#!/usr/bin/env python3
"""
Generate a professional PDF from the Paper-Writing Workflow document.
Uses fpdf2 with clean typography and structured layout.
"""

from fpdf import FPDF
import re

# ============================================================
# CUSTOM PDF CLASS
# ============================================================
def safe_text(text):
    """Replace Unicode characters with ASCII-safe equivalents."""
    replacements = {
        ' --': '--',   # em-dash
        '\u2013': '-',    # en-dash
        '\u2022': '-',    # bullet
        '\u2610': '[ ]',  # unchecked box
        '\u2611': '[x]',  # checked box
        '\u2192': '->',   # right arrow
        '\u2019': "'",    # right single quote
        '\u2018': "'",    # left single quote
        '\u201c': '"',    # left double quote
        '\u201d': '"',    # right double quote
        '\u2026': '...',  # ellipsis
        '\u00d7': 'x',    # multiplication sign
        '\u2265': '>=',   # >=
        '\u2264': '<=',   # <=
        '\u00b1': '+/-',  # plus-minus
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


class WorkflowPDF(FPDF):
    def __init__(self):
        super().__init__(orientation='P', unit='mm', format='letter')
        self.set_auto_page_break(auto=True, margin=20)
        # Colors
        self.dark_blue = (27, 58, 92)
        self.accent_blue = (46, 134, 171)
        self.accent_red = (214, 39, 40)
        self.dark_gray = (51, 51, 51)
        self.med_gray = (102, 102, 102)
        self.light_gray = (240, 240, 240)
        self.white = (255, 255, 255)
        self.table_header_bg = (27, 58, 92)
        self.table_alt_bg = (245, 248, 252)

    def header(self):
        if self.page_no() > 1:
            self.set_font('Helvetica', 'I', 8)
            self.set_text_color(*self.med_gray)
            self.cell(0, 6, 'Paper-Writing Workflow: GRBV Transcriptomics Study', align='L')
            self.cell(0, 6, f'Page {self.page_no()}', align='R', new_x='LMARGIN', new_y='NEXT')
            self.set_draw_color(*self.accent_blue)
            self.set_line_width(0.3)
            self.line(self.l_margin, self.get_y(), self.w - self.r_margin, self.get_y())
            self.ln(4)

    def footer(self):
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 7)
        self.set_text_color(*self.med_gray)
        self.cell(0, 10, 'GRBV Transcriptomics Analysis Pipeline', align='C')

    # ------ CONTENT HELPERS ------

    def add_title_page(self):
        self.add_page()
        self.ln(60)
        # Title
        self.set_font('Helvetica', 'B', 28)
        self.set_text_color(*self.dark_blue)
        self.multi_cell(0, 12, 'Paper-Writing Workflow', align='C')
        self.ln(3)
        self.set_font('Helvetica', '', 16)
        self.set_text_color(*self.accent_blue)
        self.multi_cell(0, 8, 'GRBV Transcriptomics Study', align='C')
        self.ln(8)
        # Divider
        self.set_draw_color(*self.accent_blue)
        self.set_line_width(0.8)
        cx = self.w / 2
        self.line(cx - 40, self.get_y(), cx + 40, self.get_y())
        self.ln(10)
        # Subtitle
        self.set_font('Helvetica', '', 11)
        self.set_text_color(*self.dark_gray)
        self.multi_cell(0, 6,
            safe_text('A comprehensive guide for transforming 35 figures, 13 data tables,\n'
            'and 10 analysis scripts into a publication-ready manuscript\n'
            'and conference presentation.'),
            align='C')
        self.ln(25)
        # Quick ref box
        self.set_fill_color(*self.light_gray)
        self.set_draw_color(*self.accent_blue)
        self.set_line_width(0.4)
        box_x = 35
        box_w = self.w - 70
        box_y = self.get_y()
        box_h = 48
        self.rect(box_x, box_y, box_w, box_h, style='DF')
        self.set_xy(box_x + 5, box_y + 4)
        self.set_font('Helvetica', 'B', 11)
        self.set_text_color(*self.dark_blue)
        self.cell(box_w - 10, 6, 'Quick Reference', align='C', new_x='LMARGIN', new_y='NEXT')
        self.set_xy(box_x + 10, box_y + 13)
        self.set_font('Helvetica', '', 10)
        self.set_text_color(*self.dark_gray)
        refs = [
            ('Dataset', '327 samples, 14,604 genes, 2 seasons, 4 maturity stages'),
            ('Figures ready', '35 publication-quality PNGs'),
            ('CSV tables', '13 supplementary data files'),
            ('Draft report', 'GRBV_Transcriptomics_Report.md (500 lines)'),
            ('PPT', 'GRBV_Transcriptomics_Presentation.pptx (auto-generated)'),
        ]
        for label, value in refs:
            self.set_x(box_x + 10)
            self.set_font('Helvetica', 'B', 10)
            self.cell(35, 6, label + ':', new_x='END')
            self.set_font('Helvetica', '', 10)
            self.cell(0, 6, '  ' + value, new_x='LMARGIN', new_y='NEXT')

    def add_phase_header(self, phase_num, title, timeline):
        """Add a prominent phase header."""
        self.ln(4)
        if self.get_y() > 240:
            self.add_page()
        # Phase banner
        self.set_fill_color(*self.dark_blue)
        y = self.get_y()
        self.rect(self.l_margin, y, self.w - self.l_margin - self.r_margin, 12, style='F')
        self.set_xy(self.l_margin + 3, y + 1.5)
        self.set_font('Helvetica', 'B', 14)
        self.set_text_color(*self.white)
        self.cell(0, 9, f'Phase {phase_num}: {title}', new_x='LMARGIN')
        # Timeline tag
        self.set_xy(self.w - self.r_margin - 42, y + 2)
        self.set_font('Helvetica', 'I', 10)
        self.set_text_color(200, 220, 240)
        self.cell(40, 8, timeline, align='R')
        self.set_xy(self.l_margin, y + 14)
        self.set_text_color(*self.dark_gray)
        self.ln(2)

    def add_section_header(self, text):
        """Add a subsection header."""
        if self.get_y() > 250:
            self.add_page()
        self.ln(3)
        self.set_font('Helvetica', 'B', 12)
        self.set_text_color(*self.accent_blue)
        self.cell(0, 7, text, new_x='LMARGIN', new_y='NEXT')
        self.set_draw_color(*self.accent_blue)
        self.set_line_width(0.2)
        self.line(self.l_margin, self.get_y(), self.l_margin + 60, self.get_y())
        self.ln(3)
        self.set_text_color(*self.dark_gray)

    def add_body_text(self, text, bold=False, italic=False):
        """Add body text paragraph."""
        style = ''
        if bold: style += 'B'
        if italic: style += 'I'
        self.set_font('Helvetica', style, 10)
        self.set_text_color(*self.dark_gray)
        self.multi_cell(0, 5, safe_text(text))
        self.ln(1)

    def add_blockquote(self, text):
        """Add an indented quote block."""
        text = safe_text(text)
        self.set_fill_color(245, 248, 252)
        self.set_draw_color(*self.accent_blue)
        x = self.l_margin + 5
        w = self.w - self.l_margin - self.r_margin - 10
        y0 = self.get_y()
        self.set_x(x + 4)
        self.set_font('Helvetica', 'BI', 10)
        self.set_text_color(*self.dark_blue)
        self.multi_cell(w - 8, 5.5, text)
        y1 = self.get_y()
        # Background and left bar
        self.set_fill_color(245, 248, 252)
        self.rect(x, y0 - 1, w, y1 - y0 + 3, style='F')
        self.set_line_width(0.8)
        self.line(x, y0 - 1, x, y1 + 2)
        # Re-draw text on top
        self.set_xy(x + 5, y0)
        self.set_font('Helvetica', 'BI', 10)
        self.set_text_color(*self.dark_blue)
        self.multi_cell(w - 10, 5.5, text)
        self.ln(3)

    def add_bullet(self, text, indent=0, bullet_char='-'):
        """Add a bullet point."""
        self.set_font('Helvetica', '', 10)
        self.set_text_color(*self.dark_gray)
        x = self.l_margin + 3 + indent * 5
        self.set_x(x)
        self.cell(5, 5, safe_text(bullet_char))
        self.multi_cell(self.w - x - self.r_margin - 5, 5, safe_text(text))
        self.ln(0.5)

    def add_numbered_item(self, num, text):
        """Add a numbered list item."""
        self.set_font('Helvetica', 'B', 10)
        self.set_text_color(*self.accent_blue)
        x = self.l_margin + 3
        self.set_x(x)
        self.cell(8, 5, f'{num}.')
        self.set_font('Helvetica', '', 10)
        self.set_text_color(*self.dark_gray)
        self.multi_cell(self.w - x - self.r_margin - 8, 5, safe_text(text))
        self.ln(0.5)

    def add_checkbox(self, text, checked=False):
        """Add a checkbox item."""
        self.set_font('Helvetica', '', 10)
        self.set_text_color(*self.dark_gray)
        x = self.l_margin + 5
        self.set_x(x)
        mark = '[x]' if checked else '[ ]'
        self.cell(8, 5, mark)
        self.multi_cell(self.w - x - self.r_margin - 8, 5, safe_text(text))
        self.ln(0.5)

    def add_table(self, headers, rows, col_widths=None):
        """Add a formatted table."""
        if self.get_y() > 230:
            self.add_page()
        n_cols = len(headers)
        available = self.w - self.l_margin - self.r_margin
        if col_widths is None:
            col_widths = [available / n_cols] * n_cols
        else:
            # Scale to fit
            total = sum(col_widths)
            col_widths = [w * available / total for w in col_widths]

        row_height = 6
        # Header
        self.set_fill_color(*self.table_header_bg)
        self.set_text_color(*self.white)
        self.set_font('Helvetica', 'B', 8)
        for i, h in enumerate(headers):
            self.cell(col_widths[i], row_height + 1, safe_text(h), border=1, fill=True, align='C')
        self.ln()

        # Rows
        self.set_text_color(*self.dark_gray)
        self.set_font('Helvetica', '', 8)
        for r_idx, row in enumerate(rows):
            if self.get_y() > 258:
                self.add_page()
                # Re-draw header
                self.set_fill_color(*self.table_header_bg)
                self.set_text_color(*self.white)
                self.set_font('Helvetica', 'B', 8)
                for i, h in enumerate(headers):
                    self.cell(col_widths[i], row_height + 1, safe_text(h), border=1, fill=True, align='C')
                self.ln()
                self.set_text_color(*self.dark_gray)
                self.set_font('Helvetica', '', 8)

            if r_idx % 2 == 1:
                self.set_fill_color(*self.table_alt_bg)
                fill = True
            else:
                fill = False
            for i, cell_text in enumerate(row):
                al = 'L' if i == 0 or len(cell_text) > 10 else 'C'
                self.cell(col_widths[i], row_height, ' ' + safe_text(cell_text), border=1, fill=fill, align=al)
            self.ln()
        self.ln(3)

    def add_code_block(self, text):
        """Add a monospaced code/data block."""
        if self.get_y() > 220:
            self.add_page()
        self.set_fill_color(248, 248, 248)
        self.set_draw_color(200, 200, 200)
        x = self.l_margin + 2
        w = self.w - self.l_margin - self.r_margin - 4
        y0 = self.get_y()
        lines = safe_text(text).strip().split('\n')
        line_h = 4.2
        block_h = len(lines) * line_h + 6
        if y0 + block_h > 260:
            self.add_page()
            y0 = self.get_y()
        self.rect(x, y0, w, block_h, style='DF')
        self.set_xy(x + 3, y0 + 3)
        self.set_font('Courier', '', 7.5)
        self.set_text_color(60, 60, 60)
        for line in lines:
            self.set_x(x + 3)
            self.cell(w - 6, line_h, line[:120], new_x='LMARGIN', new_y='NEXT')
        self.set_y(y0 + block_h + 2)
        self.ln(2)

    def add_divider(self):
        """Add a thin horizontal rule."""
        self.ln(3)
        self.set_draw_color(*self.med_gray)
        self.set_line_width(0.15)
        self.line(self.l_margin, self.get_y(), self.w - self.r_margin, self.get_y())
        self.ln(3)


# ============================================================
# BUILD THE PDF
# ============================================================
pdf = WorkflowPDF()
pdf.set_margins(18, 18, 18)

# --- TITLE PAGE ---
pdf.add_title_page()

# --- PHASE 1 ---
pdf.add_page()
pdf.add_phase_header(1, 'Organize Your Story Arc', 'Day 1')

pdf.add_section_header('1.1  Define the Central Narrative')
pdf.add_body_text('Your data supports this story:')
pdf.add_blockquote(
    'GRBV hijacks host transcription in a development- and environment-dependent manner, '
    'suppressing heat-stress defenses, inverting seasonal response kinetics, and targeting genes '
    'that co-localize with known disease-resistance and berry-quality QTLs  -- with implications '
    'for diagnostics and breeding.'
)

pdf.add_section_header('1.2  Recommended Paper Structure')
pdf.add_table(
    ['Section', 'Key Content', 'Figures', 'Tables'],
    [
        ['Abstract', 'Multi-dimensional RNA-Seq; 6 research Qs; inverted trajectory; biomarkers', ' --', ' --'],
        ['Introduction', 'GRBV biology, RBD symptoms, knowledge gaps, objectives', 'Fig 23', ' --'],
        ['Methods', 'RNA-Seq pipeline, DE analysis, pathway, biomarker, QTL, heat-stress', ' --', 'Table 1'],
        ['Results 1: Global', 'PCA, variance hierarchy (stage > season > infection)', 'Fig 1', ' --'],
        ['Results 2: Duration (Q1)', 'Biphasic pattern, early vs late responders', 'Fig 2, 8', 'Table S1'],
        ['Results 3: Stage (Q2)', 'Stage-dependent response, Post-veraison/Harvest peaks', 'Fig 3, 9', 'Table S2'],
        ['Results 4: Season (Q3-Q5)', 'KEY: Inverted trajectory, environment modulation', 'Fig 4-6, 10-12, 15-19', 'Table S3-S5'],
        ['Results 5: Trajectory (Q6)', 'Biphasic: acute defense -> chronic metabolism', 'Fig 7, 13', 'Table S6'],
        ['Results 6: Pathways', 'Phenylpropanoid, plant-pathogen, ABC transporters', 'Fig 14', 'Table S7'],
        ['Results 7: Heat-Stress', '115 heat genes, sHSP suppression, virus-heat split', 'Fig 34-35', 'Table S4-S6'],
        ['Results 8: Env-Titer', 'Temperature-titer-DEG correlations', 'Fig 30-33', 'Table S9'],
        ['Results 9: Biomarkers', 'Top 50 biomarkers, 3-gene panel (AUC=1.0)', 'Fig 20-22', 'Table S1'],
        ['Results 10: GWAS/QTL', '60% co-localization, resistance clusters, sugar QTL', 'Fig 24-25', 'Table S2'],
        ['Discussion', 'Developmental modulation, seasonal inversion, breeding', 'Fig 26-29', ' --'],
        ['Conclusions', '6 key conclusions', ' --', ' --'],
    ],
    col_widths=[28, 55, 22, 15]
)

pdf.add_section_header('1.3  Suggested Journal Targets')
pdf.add_table(
    ['Journal', 'Impact', 'Fit', 'Notes'],
    [
        ['Plant Cell', '~17', 'Excellent', 'Multi-dimensional systems biology, mechanistic insight'],
        ['New Phytologist', '~10', 'Excellent', 'Plant-pathogen interaction, environmental modulation'],
        ['Plant Physiology', '~8', 'Very Good', 'Comprehensive transcriptomics, heat-stress biology'],
        ['Mol. Plant-Microbe Interact.', '~5', 'Very Good', 'Virus-host interaction focus'],
        ['Frontiers in Plant Science', '~5', 'Good', 'Broad plant biology, open access'],
        ['Phytopathology', '~3', 'Good', 'Plant disease focus'],
    ],
    col_widths=[32, 10, 12, 55]
)

# --- PHASE 2 ---
pdf.add_phase_header(2, 'Figure Selection & Refinement', 'Day 2-3')

pdf.add_section_header('2.1  Main Text Figures (recommended 6-8)')
pdf.add_table(
    ['Priority', 'Figure', 'File', 'Rationale'],
    [
        ['1', 'Analysis Pipeline', '23_analysis_schematic.png', 'Orients the reader'],
        ['2', 'Global PCA', '01_global_pca_overview.png', 'Shows variance hierarchy'],
        ['3', 'Seasonal Maturity DEGs', '15_seasonal_maturity_DEG_comparison.png', 'KEY: inverted trajectory'],
        ['4', 'Pathway Landscape', '14_master_pathway_landscape.png', 'Systems-level view'],
        ['5', 'Heat vs Virus', '35_heat_vs_virus_disentanglement.png', 'Novel biology'],
        ['6', 'Biomarker Overview', '20_biomarker_overview.png', 'Translational impact'],
        ['7', 'QTL Chromosome Map', '24_biomarker_QTL_chromosome_map.png', 'Breeding relevance'],
        ['8', 'Summary Dashboard', '29_comprehensive_summary_dashboard.png', 'Integrative conclusion'],
    ],
    col_widths=[10, 22, 40, 30]
)

pdf.add_section_header('2.2  Supplementary Figures')
pdf.add_body_text('All other figures (02-07, 08-13, 16-19, 21-22, 25-28, 30-34) go as supplementary material.')

pdf.add_section_header('2.3  Figure Refinement Checklist')
checklist = [
    'Resolution >= 300 DPI (current: 150 DPI  -- re-render at 300)',
    'Font size >= 8pt in all panels',
    'Consistent color scheme across figures',
    'Panel labels (A, B, C, D) clearly visible',
    'Axis labels fully spelled out (no abbreviations)',
    'Legend positioned without overlapping data',
    'Color-blind friendly palette (check with Coblis simulator)',
]
for item in checklist:
    pdf.add_checkbox(item)

pdf.ln(2)
pdf.add_body_text('To re-render at 300 DPI: In each .py script, change \'figure.dpi\': 150 to \'figure.dpi\': 300.', bold=True)

# --- PHASE 3 ---
pdf.add_phase_header(3, 'Write the Manuscript', 'Day 3-10')

pdf.add_section_header('3.1  Writing Order (most efficient)')
writing_order = [
    'Methods  -- straightforward, sets the foundation (use Report Section 2)',
    'Results  -- present findings in logical order (use Report Section 3 as draft)',
    'Discussion  -- interpret results, compare to literature',
    'Introduction  -- frame the problem knowing what you found',
    'Abstract  -- summarize after everything else is written',
    'Title  -- finalize last',
]
for i, item in enumerate(writing_order, 1):
    pdf.add_numbered_item(i, item)

pdf.add_section_header('3.2  Key Numbers to Cite')
pdf.add_code_block("""Dataset:         327 samples, 14,604 genes (from 28,168 raw)
Infection groups: RB(+) infected vs RB(-) uninfected controls
Seasons:         2021 (cooler) and 2022 (hotter, +7.6C at Harvest)
Stages:          Pre-veraison, Veraison, Post-veraison, Harvest
Cohorts:         7 infection-year cohorts (2016-2022)

DE thresholds:   padj < 0.05, |log2FC| >= 1
Enrichment:      Fisher's exact, BH FDR correction

Key DEG counts:
  - 2021 Harvest: 1,192 DEGs (peak in 2021)
  - 2022 Post-veraison: 565 DEGs (peak in 2022)
  - Cross-season correlation: Pearson r = 0.33-0.69

Heat-stress:
  - 115 heat-stress genes expressed (9 categories)
  - 57 heat genes are GRBV DEGs
  - 29/36 sHSPs downregulated at 2021 Harvest in RB(+)
  - 23 heat genes with BOTH virus and environment effects
  - HSF-sHSP module correlation: r=0.625, p=7.2e-37

Biomarkers:
  - 12 universal biomarkers (robust across all conditions)
  - 3-gene panel AUC = 1.000
  - 60% of biomarkers near known QTLs (within 1 Mb)
  - 18 biomarkers in disease-resistance clusters
  - Top host biomarker: F-box/LRR protein (AUC = 0.880)

Environmental:
  - 2022 Harvest was +7.6C warmer than 2021
  - sHSP module score: +1.64 (2022 Harvest RB-) vs +0.46 (2021)
  - GRBV reduces sHSP induction by ~50% (RB+ vs RB-)""")

pdf.add_section_header('3.3  Manuscript Draft Sources')
pdf.add_table(
    ['Section', 'Source File', 'Lines'],
    [
        ['Introduction', 'GRBV_Transcriptomics_Report.md', '29-41'],
        ['Methods', 'GRBV_Transcriptomics_Report.md', '44-94'],
        ['Results (Q1-Q6)', 'GRBV_Transcriptomics_Report.md', '98-394'],
        ['Discussion', 'GRBV_Transcriptomics_Report.md', '396-425'],
        ['Conclusions', 'GRBV_Transcriptomics_Report.md', '435-448'],
        ['Heat stress results', 'heat_stress_signatures.py output', '(run script)'],
        ['Env-titer results', 'environment_titer_integration.py output', '(run script)'],
    ],
    col_widths=[22, 45, 15]
)

pdf.add_section_header('3.4  Literature to Cite')
refs = [
    'GRBV original characterization (Al Rwahnih et al., 2013)',
    'GRBV genome and classification (Varsani et al., 2017)',
    'GRBV transcriptomics (Blanco-Ulate et al., 2017; Rumbaugh et al., 2021)',
    'GRBV impacts on berry metabolism (Girardello et al., 2020)',
    'Grapevine heat-stress HSPs (Liu et al., 2012; Wang et al., 2014)',
    'Grapevine QTL/GWAS resistance loci (Zyprian et al., 2016; Possamai et al., 2020)',
    'Biomarker RT-qPCR diagnostics in plants (various)',
]
for r in refs:
    pdf.add_bullet(r)

# --- PHASE 4 ---
pdf.add_phase_header(4, 'Supplementary Materials', 'Day 8-12')

pdf.add_section_header('4.1  Supplementary Tables')
pdf.add_table(
    ['Table', 'Source CSV', 'Description'],
    [
        ['S1', 'GRBV_biomarkers_top50.csv', 'Complete biomarker list with scores'],
        ['S2', 'GRBV_biomarker_QTL_crossref.csv', 'QTL cross-reference details'],
        ['S3', 'universal_DEGs_across_approaches.csv', '908 robust DEGs'],
        ['S4', 'heat_stress_DEGs_all_conditions.csv', 'Heat-gene DE results'],
        ['S5', 'heat_gene_temperature_correlations.csv', 'Temperature correlations'],
        ['S6', 'heat_DEGs_recurrence.csv', 'Recurrent heat DEGs'],
        ['S7', 'Q1_pooled_vs_stratified.csv', 'Analytical robustness (Q1)'],
        ['S8', 'Q2_pooled_vs_stratified.csv', 'Analytical robustness (Q2)'],
        ['S9', 'integrated_environment_titer_DEG.csv', 'Env-titer-DEG integration'],
    ],
    col_widths=[8, 42, 42]
)

pdf.add_section_header('4.2  Data Availability Statement')
pdf.add_code_block("""Raw RNA-Seq data have been deposited in [NCBI GEO/SRA under accession GSEXXXXX].
All analysis scripts, processed count matrices, and supplementary data files are
available at [GitHub repository URL]. The Vitis vinifera functional annotation
database (org.Vvinifera.eg.db) was obtained from Bioconductor.""")

# --- PHASE 5 ---
pdf.add_phase_header(5, 'Review & Polish', 'Day 12-15')

pdf.add_section_header('5.1  Internal Review Checklist')
review_items = [
    'All figure references in text match actual figure numbers',
    'All p-values and fold-changes verified against script outputs',
    'Methods section matches what scripts actually do (check thresholds)',
    'Supplementary table column headers explained',
    'No orphan references (every citation in text is in bibliography)',
    'Author contributions drafted',
    'Acknowledgments and funding sources listed',
    'Conflict of interest statement prepared',
]
for item in review_items:
    pdf.add_checkbox(item)

pdf.add_section_header('5.2  Common Reviewer Concerns')
pdf.add_table(
    ['Concern', 'Pre-emptive Response'],
    [
        ['Why Welch\'s t-test?', 'Large sample (n=327) makes parametric robust; show supplementary DESeq2 comparison'],
        ['Multiple testing burden?', 'Each question uses independent BH correction; 908 universal DEGs validated'],
        ['Batch effects between seasons?', 'Addressed by Q5 (baseline environmental in RB- controls); no PCA batch clustering'],
        ['Small sub-group sample sizes?', 'Minimum n=3 per group; stratified analysis confirms pooled results'],
        ['Biomarker validation?', 'RT-qPCR as future work; 5-fold stratified cross-validation already performed'],
    ],
    col_widths=[30, 75]
)

# --- PHASE 6 ---
pdf.add_phase_header(6, 'Generate Presentation', 'Automated')

pdf.add_body_text('Run the following command to generate the slide deck:')
pdf.add_code_block('python3 generate_presentation.py')
pdf.add_body_text(
    'This produces GRBV_Transcriptomics_Presentation.pptx with all 35 figures organized '
    'into a 44-slide conference-ready deck with section dividers, key-number panels, '
    'bullet-point annotations, and speaker notes.'
)

# --- FILE INVENTORY ---
pdf.add_divider()
pdf.add_section_header('File Inventory')

pdf.add_body_text('Analysis Scripts (10)', bold=True)
scripts = [
    'analysis.py  -- Core DE analysis (Q1-Q6), PCA, Figs 1-7',
    'pathway_enrichment.py  -- GO/KEGG enrichment, Figs 8-14',
    'seasonal_maturity_progression.py  -- Season x Stage, Figs 15-19',
    'biomarker_discovery.py  -- Biomarker scoring & classification, Figs 20-22',
    'gwas_qtl_crossref.py  -- QTL cross-reference, Figs 24-25',
    'stratified_vs_pooled_comparison.py  -- Robustness check, Figs 26-29',
    'environment_titer_integration.py  -- Weather + titer, Figs 30-33',
    'heat_stress_signatures.py  -- Heat-stress signatures, Figs 34-35',
    'analysis_schematic.py  -- Pipeline diagram, Fig 23',
    'generate_report_pdf.py  -- PDF report generator',
]
for i, s in enumerate(scripts, 1):
    pdf.add_numbered_item(i, s)

pdf.ln(2)
pdf.add_body_text('Key Output Files', bold=True)
outputs = [
    'GRBV_Transcriptomics_Report.md  -- Full draft report (500 lines)',
    'GRBV_Transcriptomics_Presentation.pptx  -- 44-slide conference deck',
    '35 PNG figures (01-35)  -- Publication-quality plots',
    '13 CSV files  -- Supplementary data tables',
]
for o in outputs:
    pdf.add_bullet(o)

# --- SAVE ---
output_path = '/Users/prempratapsingh/Data_Analysis/PAPER_WRITING_WORKFLOW.pdf'
pdf.output(output_path)
print(f'PDF saved: {output_path}')
print(f'Pages: {pdf.page_no()}')
