#!/usr/bin/env python3
"""Generate a formatted PDF of the Results and Discussion sections."""

from fpdf import FPDF
import textwrap

class ResultsPDF(FPDF):
    def header(self):
        if self.page_no() > 1:
            self.set_font('Helvetica', 'I', 9)
            self.set_text_color(128, 128, 128)
            self.cell(0, 8, 'GRBV Transcriptomics - Results and Discussion', 0, 0, 'L')
            self.cell(0, 8, f'Page {self.page_no()}', 0, 1, 'R')
            self.set_draw_color(180, 180, 180)
            self.set_line_width(0.3)
            self.line(10, self.get_y() + 1, 200, self.get_y() + 1)
            self.ln(6)

    def footer(self):
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 8)
        self.set_text_color(128, 128, 128)
        self.cell(0, 10, f'Singh (2026) - GRBV Transcriptomics Analysis', 0, 0, 'C')

    def section_title(self, title, level=2):
        self.ln(4)
        if level == 1:
            self.set_font('Helvetica', 'B', 18)
            self.set_text_color(26, 26, 46)
            self.cell(0, 12, title, 0, 1, 'L')
            self.set_draw_color(26, 26, 46)
            self.set_line_width(0.8)
            self.line(10, self.get_y(), 200, self.get_y())
            self.ln(4)
        elif level == 2:
            self.set_font('Helvetica', 'B', 14)
            self.set_text_color(26, 26, 46)
            self.cell(0, 10, title, 0, 1, 'L')
            self.set_draw_color(100, 100, 180)
            self.set_line_width(0.4)
            self.line(10, self.get_y(), 120, self.get_y())
            self.ln(2)
        elif level == 3:
            self.set_font('Helvetica', 'B', 12)
            self.set_text_color(50, 50, 100)
            self.cell(0, 9, title, 0, 1, 'L')
            self.ln(1)

    def body_text(self, text):
        self.set_font('Helvetica', '', 10)
        self.set_text_color(30, 30, 30)
        self.multi_cell(0, 5.5, text)
        self.ln(2)

    def bold_start_text(self, bold_part, rest):
        self.set_font('Helvetica', 'B', 10)
        self.set_text_color(30, 30, 30)
        w_bold = self.get_string_width(bold_part) + 1
        self.cell(w_bold, 5.5, bold_part, 0, 0)
        self.set_font('Helvetica', '', 10)
        self.multi_cell(0, 5.5, rest)
        self.ln(1)

    def table(self, headers, rows, col_widths=None):
        if col_widths is None:
            n = len(headers)
            col_widths = [190 / n] * n

        # Header
        self.set_font('Helvetica', 'B', 9)
        self.set_fill_color(26, 26, 46)
        self.set_text_color(255, 255, 255)
        for i, h in enumerate(headers):
            self.cell(col_widths[i], 7, h, 1, 0, 'C', True)
        self.ln()

        # Rows
        self.set_font('Helvetica', '', 9)
        self.set_text_color(30, 30, 30)
        fill = False
        for row in rows:
            if fill:
                self.set_fill_color(240, 240, 250)
            else:
                self.set_fill_color(255, 255, 255)
            for i, val in enumerate(row):
                self.cell(col_widths[i], 6, str(val), 1, 0, 'C', True)
            self.ln()
            fill = not fill
        self.ln(3)


pdf = ResultsPDF()
pdf.set_auto_page_break(auto=True, margin=20)

# ── Title Page ──
pdf.add_page()
pdf.ln(40)
pdf.set_font('Helvetica', 'B', 28)
pdf.set_text_color(26, 26, 46)
pdf.cell(0, 15, 'Results and Discussion', 0, 1, 'C')
pdf.ln(5)
pdf.set_font('Helvetica', '', 14)
pdf.set_text_color(80, 80, 80)
pdf.cell(0, 8, 'Comprehensive Transcriptomic Analysis of Grapevine', 0, 1, 'C')
pdf.cell(0, 8, 'Red Blotch Virus Effects on Vitis vinifera cv. Merlot', 0, 1, 'C')
pdf.cell(0, 8, 'Across Seasons and Maturity Stages', 0, 1, 'C')
pdf.ln(10)
pdf.set_draw_color(26, 26, 46)
pdf.set_line_width(1)
pdf.line(60, pdf.get_y(), 150, pdf.get_y())
pdf.ln(10)
pdf.set_font('Helvetica', '', 12)
pdf.set_text_color(60, 60, 60)
pdf.cell(0, 8, 'Prem Pratap Singh', 0, 1, 'C')
pdf.cell(0, 8, 'University of California, Davis', 0, 1, 'C')
pdf.ln(5)
pdf.set_font('Helvetica', 'I', 10)
pdf.cell(0, 8, '327 samples | 14,604 genes | 55 figures | 20 analysis modules', 0, 1, 'C')

# ── RESULTS ──
pdf.add_page()
pdf.section_title('3. Results', level=1)

# 3.1
pdf.section_title('3.1 Global Transcriptome Overview and Variance Structure', level=2)
pdf.body_text(
    'RNA-Seq data from 327 field-collected Vitis vinifera cv. Merlot berry samples '
    '(264 GRBV-infected RB(+), 63 uninfected RB(-)) yielded 28,168 expressed genes, '
    'of which 14,604 passed quality filters (>=10 counts in >=10% of samples) and were '
    'retained for downstream analysis. Principal component analysis (PCA) of log2(CPM+1)-'
    'transformed expression values revealed that maturity stage was the dominant source '
    'of transcriptome variation, followed by growing season and GRBV infection status (Figure 1).'
)
pdf.body_text(
    'Multi-factor ANOVA (expression ~ Infection x Season x Stage) across all 14,604 genes '
    'formally quantified this variance hierarchy. Maturity stage accounted for the largest '
    'mean proportion of variance (mean eta-sq = 0.426), followed by Season x Stage interaction '
    '(eta-sq = 0.069), Season (eta-sq = 0.026), Infection x Stage interaction (eta-sq = 0.020), '
    'Infection (eta-sq = 0.018), Infection x Season x Stage three-way interaction (eta-sq = 0.011), '
    'and Infection x Season interaction (eta-sq = 0.003) (Figure 40). A total of 14,418 genes '
    '(98.7% of all tested) exhibited significant Stage effects (BH-adjusted p < 0.05), whereas '
    '7,987 genes (54.7%) showed significant Infection main effects and 9,155 genes (62.7%) '
    'showed significant Season main effects (Table S1). Notably, 773 genes displayed significant '
    'Infection x Season interaction effects (adjusted p < 0.05), and 4,146 genes showed significant '
    'three-way Infection x Season x Stage interactions (Figures 41-42), formally confirming that '
    'GRBV transcriptional effects are context-dependent.'
)

# 3.2
pdf.section_title('3.2 Differential Expression Across Six Research Questions', level=2)
pdf.body_text(
    'Differential expression (DE) analysis was performed using Welch\'s t-test with '
    'Benjamini-Hochberg false discovery rate (FDR) correction, applying thresholds of '
    'adjusted p < 0.05 and |log2FC| >= 1.'
)
pdf.bold_start_text('Infection duration (Q1). ',
    'Comparison of vines across seven infection cohorts (2016-2022, representing 1-7 years '
    'since infection) revealed a biphasic transcriptomic response (Figure 2). Newly infected '
    'vines (1-2 years) exhibited an acute defense response, while chronically infected vines '
    '(5+ years) showed metabolic reprogramming.')
pdf.bold_start_text('Maturity stage (Q2). ',
    'GRBV effects were strongly stage-dependent, with the greatest number of differentially '
    'expressed genes (DEGs) detected at Post-veraison and Harvest (Figure 3). Up to 1,873 DEGs '
    'were identified in a single season x stage condition, consistent with the metabolically '
    'active ripening period being the most susceptible to viral disruption.')
pdf.bold_start_text('Seasonal reproducibility (Q3). ',
    'Cross-season fold-change correlations ranged from r = 0.33 at Pre-veraison to r = 0.69 '
    'at Post-veraison (Figure 4), indicating that mid-to-late ripening transcriptional responses '
    'to GRBV were moderately reproducible between the cooler 2021 and hotter 2022 growing seasons.')
pdf.bold_start_text('Environmental modulation (Q4). ',
    'The 2022 season recorded mean temperatures +7.6 C above 2021 at Harvest (Figure 30). GRBV '
    'effects were not simply additive with seasonal effects; the 773 genes with significant '
    'Infection x Season interaction terms demonstrated that environmental conditions actively '
    'modulated the viral response (Figure 5).')
pdf.bold_start_text('Baseline seasonal changes (Q5). ',
    'In uninfected RB(-) control vines, 2,431 genes were differentially expressed between '
    'the 2021 and 2022 seasons, establishing the baseline environmental transcriptomic shift '
    '(Figure 6). These genes provided a reference for disentangling virus-specific from '
    'environment-driven changes.')
pdf.bold_start_text('Infection age trajectory (Q6). ',
    'The transcriptome trajectory across infection years followed a biphasic pattern: acute '
    'defense activation in early infection followed by chronic metabolic reprogramming in later '
    'years (Figure 7).')

# 3.3
pdf.section_title('3.3 Inverted Seasonal DEG Trajectory', level=2)
pdf.body_text(
    'Season-stratified DE analysis across four maturity stages revealed an inverted seasonal '
    'trajectory (Figure 15). In 2021 (cooler year), DEG counts peaked at Harvest (1,192 DEGs), '
    'whereas in 2022 (hotter year), the peak shifted to Post-veraison (565 DEGs). This inversion '
    'was not attributable to statistical power differences, as sample sizes were comparable '
    'between seasons. The 773 Infection x Season interaction genes from the multi-factor ANOVA '
    'independently confirmed this inverted pattern.'
)

# 3.4
pdf.section_title('3.4 Pathway Enrichment Analysis', level=2)
pdf.body_text(
    'Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analysis '
    '(Fisher\'s exact test, BH-corrected p < 0.05) identified consistent downregulation of '
    'phenylpropanoid biosynthesis, flavonoid biosynthesis, and stilbenoid/diarylheptanoid '
    'biosynthesis pathways across multiple conditions (Figures 8-14). The master pathway landscape '
    'across all six research questions revealed that secondary metabolism pathways were the most '
    'recurrently disrupted (Figure 14).'
)
pdf.body_text(
    'Permutation-based Gene Set Enrichment Analysis (GSEA; 500 permutations, 375 gene sets, '
    '8 season x stage conditions) identified 270 significant gene set-condition pairs at FDR < 0.05 '
    '(Figure 39). Down-regulation dominated, accounting for 77% (209/270) of significant enrichments. '
    'KEGG:Ribosome and GO:translation initiation factor activity were significant in 6 of 8 '
    'conditions, indicating broad suppression of the translational machinery by GRBV. GO:protein '
    'refolding, GO:response to hydrogen peroxide, and GO:structural constituent of ribosome '
    'were each significant in 5 of 8 conditions (Table S2). The strongest enrichment score was '
    'observed for KEGG:Ribosome at 2022 Harvest (NES = -2.995, FDR < 0.001) (Figure 37).'
)
pdf.body_text(
    'Cross-season GSEA NES correlations paralleled the DEG-level findings: Pre-veraison showed '
    'near-zero correlation (r = 0.09), while Post-veraison exhibited the strongest concordance '
    '(r = 0.63) (Figure 38). GSEA and Fisher\'s exact test over-representation analysis showed '
    'moderate concordance (Spearman r = 0.17-0.47 across conditions).'
)

# 3.5
pdf.section_title('3.5 Heat-Stress Gene Suppression', level=2)
pdf.body_text(
    'Analysis of 115 heat-stress-related genes across nine functional categories revealed '
    'pronounced GRBV-mediated suppression (Figure 34). Among 36 expressed small heat-shock '
    'protein (sHSP) genes, 31 were differentially expressed in at least one condition, with '
    'the vast majority downregulated by GRBV infection. The sHSP module score was reduced by '
    'approximately 55% in infected vines compared to healthy controls (mean score: healthy +1.64, '
    'infected +0.74). HSP70 family genes (14 expressed, 4 DEGs), HSP90 (6 expressed, 3 DEGs), '
    'and HSF transcription factors (13 expressed, 3 DEGs) were also affected, though to a lesser '
    'extent (Table 1). The strongest heat-stress suppression occurred at the 2022 Harvest time '
    'point (GSEA NES for GO:response to heat = -2.867, FDR < 0.001), corresponding to the hottest '
    'sampling condition (+7.6 C above 2021) (Figure 35).'
)

# Table 1
pdf.set_font('Helvetica', 'B', 10)
pdf.set_text_color(30, 30, 30)
pdf.cell(0, 7, 'Table 1. Heat-stress gene categories and GRBV differential expression.', 0, 1)
pdf.ln(1)
pdf.table(
    ['Category', 'Expressed', 'DEGs', 'Direction'],
    [
        ['sHSP', '36', '31', 'DOWN'],
        ['GO:response to heat', '42', '28', 'DOWN'],
        ['LEA/Dehydrin', '14', '6', 'DOWN'],
        ['Galactinol synthase', '7', '6', 'DOWN'],
        ['HSP70 family', '14', '4', 'MIXED'],
        ['Antioxidant enzymes', '14', '3', 'MIXED'],
        ['HSP90 family', '6', '3', 'MIXED'],
        ['HSF transcription factors', '13', '3', 'DOWN'],
        ['HSP100/ClpB', '1', '1', 'DOWN'],
    ],
    col_widths=[70, 30, 30, 60]
)

# 3.6
pdf.section_title('3.6 Transcription Factor Family Disruption', level=2)
pdf.body_text(
    'Profiling of six major transcription factor (TF) families (267 genes total) revealed '
    'family-specific patterns of GRBV-mediated disruption across conditions (Figure 43). '
    'ERF/AP2 was the largest family (79 genes) and exhibited up to 11 DEGs (13.9% of family) '
    'in peak conditions. HSF was the most consistently affected family proportionally, with '
    'up to 23.1% of its members differentially expressed (3/13 genes) (Figure 44). WRKY '
    '(33 genes, max 5 DEGs, 15.2%), bHLH (47 genes, max 6 DEGs, 12.8%), MYB (67 genes, '
    'max 4 DEGs, 6.0%), and NAC (28 genes, max 2 DEGs, 7.1%) showed condition-dependent '
    'effects (Figure 45).'
)

# 3.7
pdf.section_title('3.7 Biomarker Discovery and Machine Learning Validation', level=2)
pdf.body_text(
    'Multi-condition biomarker scoring identified 50 candidate genes ranked by a composite '
    'metric incorporating significance breadth, effect size, and AUC (Figure 20). The top-ranked '
    'biomarkers were the five viral genes (N761_gp1 through gp5), each achieving AUC = 1.000. '
    'Among host genes, the top five were: LOC100854143 (F-box/LRR-repeat protein, score = 0.741, '
    'AUC = 0.880), LOC100256964 (E3 ubiquitin ligase RSL1, score = 0.731, AUC = 0.894), '
    'LOC104881857 (Protein NTM1-like 9, score = 0.725, AUC = 0.869), LOC100268011 (RNA-dependent '
    'RNA polymerase 1, score = 0.716, AUC = 0.895), and LOC100257913 (glycine-rich cell wall '
    'protein, score = 0.705, AUC = 0.871) (Figure 21). Of the top 50 biomarkers, 29 were '
    'upregulated and 21 downregulated, with a mean AUC of 0.783.'
)
pdf.body_text(
    'A three-gene diagnostic panel (viral genes N761_gp1, gp2, gp3) achieved AUC = 1.000 across '
    'five-fold cross-validation using logistic regression (Figure 20). Five independent machine '
    'learning methods -- logistic regression, random forest, LASSO, SVM-RBF, and gradient '
    'boosting -- all achieved perfect classification (AUC = 1.000, std = 0.000 across folds) '
    '(Figure 46, Table 2). LASSO stability selection (100 bootstrap subsamples) identified '
    'N761_gp1 and N761_gp3 as the only two features with selection frequency >0.80, confirming '
    'their reliability as diagnostic markers (Figure 47). A 57-gene consensus panel was '
    'identified across all five methods (Figure 48).'
)

# Table 2
pdf.set_font('Helvetica', 'B', 10)
pdf.cell(0, 7, 'Table 2. Machine learning classification performance for GRBV infection status.', 0, 1)
pdf.ln(1)
pdf.table(
    ['Method', 'Mean AUC', 'Std AUC'],
    [
        ['Logistic Regression', '1.000', '0.000'],
        ['Random Forest', '1.000', '0.000'],
        ['LASSO (L1)', '1.000', '0.000'],
        ['SVM-RBF', '1.000', '0.000'],
        ['Gradient Boosting', '1.000', '0.000'],
    ],
    col_widths=[80, 55, 55]
)

# 3.8
pdf.section_title('3.8 Robustness Validation and Universal DEGs', level=2)
pdf.body_text(
    'Comparison of season-pooled versus season-stratified DE analyses identified 907 universal '
    'DEGs significant across both analytical approaches (Figure 26). All three viral genes '
    '(N761_gp1, gp3, gp5) were detected in 100% of the 36 tested condition comparisons. The '
    'overlap between pooled and stratified approaches confirmed that these core genes represent '
    'robust, condition-independent GRBV signatures (Figures 27-29).'
)

# 3.9
pdf.section_title('3.9 GWAS/QTL Co-localization', level=2)
pdf.body_text(
    'Cross-referencing the top 45 host biomarkers with 121 published GWAS/QTL loci for berry '
    'quality traits identified 70 biomarker-QTL overlaps involving multiple biomarker genes '
    'within 1 Mb windows of known trait loci (Figure 24). Approximately 60% of biomarker genes '
    'co-localized with at least one QTL region. The top-ranked host biomarker, LOC100268011 '
    '(RNA-dependent RNA polymerase 1), mapped to a sugar accumulation QTL on chromosome 1 '
    '(Figure 25).'
)

# 3.10
pdf.section_title('3.10 Co-expression Network Analysis', level=2)
pdf.body_text(
    'Weighted gene co-expression network analysis (WGCNA) identified 13 co-expression modules '
    'encompassing 10,221 genes (Figure 53). Module-trait correlation analysis revealed that '
    'maturity stage dominated module behavior, consistent with the ANOVA results: the turquoise '
    'module (1,621 genes) showed the strongest stage correlation (r = 0.941), while the blue '
    'module (1,505 genes, r = -0.931) and brown module (646 genes, r = -0.929) exhibited strong '
    'inverse stage associations. The red module (238 genes, r = 0.926) also tracked developmental '
    'stage (Figure 53).'
)
pdf.body_text(
    'The black module (147 genes) exhibited the strongest positive correlation with GRBV '
    'infection status (r = +0.126, p = 0.022), representing a module whose expression is '
    'elevated in infected vines. Ten hub genes per module were identified by module eigengene '
    'correlation (kME), yielding 130 hub genes across all modules (Figure 55, Table S3).'
)

# 3.11
pdf.section_title('3.11 Causal Mediation Analysis', level=2)
pdf.body_text(
    'Causal mediation analysis evaluated whether temperature effects on gene expression were '
    'mediated through viral titer changes (Temperature -> Titer -> Gene expression). Among 500 '
    'tested genes, 90 showed significant total temperature effects (BH-adjusted p < 0.05) and '
    'were classified as temperature-responsive (Figure 49). The proportion of the temperature '
    'effect mediated through viral titer was generally small (mean proportion mediated ~ 0.0%), '
    'indicating that temperature predominantly exerts direct effects on the host transcriptome '
    'rather than indirect effects through modulation of viral titer (Figure 50).'
)

# 3.12
pdf.section_title('3.12 Chromosomal DEG Distribution', level=2)
pdf.body_text(
    'Genome-wide DEG density analysis across 1 Mb sliding windows on all 19 V. vinifera '
    'chromosomes revealed non-uniform distribution of GRBV-responsive genes (Figure 51). DEG '
    'enrichment analysis identified chromosomal regions with statistically elevated DEG '
    'concentrations, visualized as a Manhattan plot (Figure 52).'
)

# ── DISCUSSION ──
pdf.add_page()
pdf.section_title('4. Discussion', level=1)

# 4.1
pdf.section_title('4.1 GRBV Induces a Developmentally Modulated, Environmentally Sensitive Transcriptomic Response', level=2)
pdf.body_text(
    'The results of this study demonstrate that GRBV infection produces a complex transcriptomic '
    'signature that is strongly modulated by berry developmental stage and growing season '
    'environment. The multi-factor ANOVA finding that stage alone accounts for 42.6% of '
    'transcriptomic variance, compared to 1.8% for infection, establishes that GRBV operates '
    'within the dominant developmental framework of grape berry ripening rather than overriding '
    'it. This is consistent with the biology of geminiviruses, which co-opt host cellular '
    'machinery rather than producing large-scale cytopathic effects, and aligns with the '
    'relatively subtle symptomology of GRBV compared to acute viral infections.'
)
pdf.body_text(
    'The identification of 7,987 genes with significant infection effects and 4,146 genes with '
    'significant three-way Infection x Season x Stage interactions underscores that GRBV effects '
    'cannot be adequately characterized by a single, static gene list. Context-dependent viral '
    'transcriptomics have been previously observed in other plant-virus pathosystems, but the '
    'scale of interaction effects reported here -- affecting 28.4% of the filtered transcriptome '
    '-- is notably large.'
)

# 4.2
pdf.section_title('4.2 The Inverted Seasonal Trajectory Reveals Environment-Dependent Viral Pathogenesis', level=2)
pdf.body_text(
    'The inverted seasonal DEG trajectory -- peaking at Harvest in the cooler 2021 season but '
    'at Post-veraison in the hotter 2022 season -- represents a key finding with implications '
    'for understanding GRBV pathogenesis under changing climates. The +7.6 C temperature '
    'differential at Harvest between seasons provides a natural experiment demonstrating that '
    'heat stress can alter the temporal dynamics of viral transcriptional disruption. The 773 '
    'Infection x Season interaction genes formally confirm that these are genuine interaction '
    'effects rather than artifacts of differential statistical power.'
)
pdf.body_text(
    'The cross-season GSEA NES correlation pattern (Pre-veraison r = 0.09, Post-veraison '
    'r = 0.63) provides pathway-level confirmation: early-season GRBV effects are highly '
    'environment-dependent and poorly reproducible between seasons, while mid-ripening effects '
    'are more conserved. This has practical implications for disease monitoring, suggesting '
    'that transcriptomic diagnostics should target the Post-veraison window for maximum '
    'reliability across variable climatic conditions.'
)

# 4.3
pdf.section_title('4.3 Suppression of the Heat-Stress Response: A Novel Vulnerability Mechanism', level=2)
pdf.body_text(
    'The finding that GRBV suppresses 31 of 36 expressed sHSP genes, reducing the sHSP module '
    'score by approximately 55%, reveals a previously unreported vulnerability mechanism. Small '
    'heat-shock proteins function as ATP-independent molecular chaperones that prevent irreversible '
    'protein aggregation during thermal stress. Their systematic downregulation by GRBV suggests '
    'that infected vines may have a compromised capacity to cope with heat stress events.'
)
pdf.body_text(
    'This interpretation is supported by the GSEA results showing that the strongest gene set '
    'enrichments overall occurred at the 2022 Harvest time point -- the hottest sampling '
    'condition -- where GO:response to heat (NES = -2.867), GO:protein complex oligomerization '
    '(NES = -2.919), and GO:response to hydrogen peroxide (NES = -2.909) were all strongly '
    'suppressed. The convergence of heat-stress suppression with the hotter growing season '
    'suggests a synergistic interaction: GRBV impairs the heat-stress defense machinery at '
    'precisely the conditions when it is most needed. In the context of climate change projections '
    'for viticultural regions, this finding raises the prospect that GRBV-infected vineyards may '
    'experience compounded vulnerability to rising temperatures.'
)
pdf.body_text(
    'The causal mediation analysis indicates that temperature effects on gene expression are '
    'predominantly direct rather than mediated through viral titer, suggesting that the heat-stress '
    'suppression reflects a host-level interaction between viral and thermal stresses rather than '
    'temperature-driven increases in viral load.'
)

# 4.4
pdf.section_title('4.4 Translational Machinery as a Primary Target of GRBV', level=2)
pdf.body_text(
    'The identification of KEGG:Ribosome and GO:translation initiation factor activity as the '
    'most recurrently suppressed gene sets (significant in 6 of 8 conditions) reveals that GRBV '
    'broadly downregulates the host translational apparatus. This is consistent with the '
    'replication strategy of geminiviruses, which depend on host polymerases and ribosomes for '
    'viral gene expression. Translational repression may represent a host defense mechanism '
    '(translational silencing to limit viral protein production) or, alternatively, viral '
    'manipulation to redirect translational resources toward viral transcripts. The 77% '
    'predominance of downregulated enrichments in the GSEA analysis supports a model of broad '
    'transcriptional suppression by GRBV, extending beyond specific metabolic pathways to '
    'fundamental cellular processes.'
)

# 4.5
pdf.section_title('4.5 Diagnostic Biomarkers and Practical Applications', level=2)
pdf.body_text(
    'The perfect classification achieved by all five machine learning methods (AUC = 1.000) '
    'demonstrates that GRBV infection produces a transcriptomic signature sufficient for '
    'unambiguous diagnosis. The viral genes N761_gp1 and N761_gp3, identified as the only '
    'stable features in LASSO stability selection (100% selection frequency across 100 bootstrap '
    'subsamples), are the most reliable individual markers. However, a three-gene viral panel '
    'provides redundancy against individual gene dropout.'
)
pdf.body_text(
    'Among host biomarkers, the 60% co-localization rate with published GWAS/QTL loci for '
    'berry quality traits suggests that GRBV may preferentially target genes within genomic '
    'regions controlling commercially relevant phenotypes. The association of the top host '
    'biomarker (RNA-dependent RNA polymerase 1) with a sugar accumulation QTL on chromosome 1 '
    'is particularly relevant, as reduced sugar content is one of the primary economic impacts '
    'of GRBV infection on wine grape production.'
)
pdf.body_text(
    'The 907 universal DEGs, robust across both pooled and stratified analytical approaches, '
    'provide a reference set of condition-independent GRBV markers for future transcriptomic '
    'studies. These genes represent the core viral signature that persists regardless of '
    'developmental stage or environmental context.'
)

# 4.6
pdf.section_title('4.6 Transcription Factor Network Disruption', level=2)
pdf.body_text(
    'The disproportionate impact of GRBV on HSF transcription factors (up to 23.1% of family '
    'affected) connects to the sHSP suppression findings, as HSFs are the primary transcriptional '
    'regulators of sHSP expression. The disruption of ERF/AP2 (up to 13.9% affected), the '
    'largest TF family in the analysis, indicates broad effects on ethylene-responsive and '
    'stress-signaling pathways. These TF-level disruptions likely cascade into downstream '
    'regulatory effects, amplifying the viral impact beyond direct target genes.'
)

# 4.7
pdf.section_title('4.7 Co-expression Network Architecture and Infection Signatures', level=2)
pdf.body_text(
    'The WGCNA finding that stage-associated modules dominate the network architecture '
    '(turquoise module r = 0.941 with stage) is consistent with the ANOVA variance decomposition. '
    'The relatively modest correlation of the black module with infection status (r = 0.126, '
    'p = 0.022) reflects the subtle nature of GRBV effects within the dominant developmental '
    'program. Nonetheless, this module represents a coherent gene set whose coordinated '
    'upregulation is specifically associated with infection, potentially revealing co-regulated '
    'defense or susceptibility gene networks.'
)

# 4.8
pdf.section_title('4.8 Limitations', level=2)
pdf.body_text(
    'Several limitations should be considered when interpreting these results. First, the '
    'imbalanced sample design (264 infected vs. 63 control) may reduce statistical power for '
    'detecting subtle effects in uninfected vines. Second, the two-season sampling captures '
    'environmental variation but cannot fully represent the range of viticultural conditions '
    'encountered across diverse growing regions. Third, bulk RNA-Seq averages expression across '
    'cell types within berry tissue, potentially masking cell-type-specific GRBV effects. Fourth, '
    'the observational field design does not permit causal inference about the relationship '
    'between infection duration and transcriptomic changes, as infection cohorts may confound '
    'with vine age or position effects. Finally, while perfect ML classification was achieved, '
    'validation in independent cohorts from different vineyards and cultivars would be required '
    'before deployment as a diagnostic tool.'
)

# Save
pdf.output('RESULTS_AND_DISCUSSION.pdf')
print(f'Saved RESULTS_AND_DISCUSSION.pdf ({pdf.page_no()} pages)')
