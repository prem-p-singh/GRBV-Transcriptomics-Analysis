"""
Add 7 agenda-gap-extension slides to GRBV_Funding_Committee_Presentation.pptx,
inserted after the existing "New dimensions" block (slide 64) and before
SECTION 7 / Synthesis.

Each slide pairs one of the new supplementary figures (S73-S79) with the
headline finding from the corresponding FARM rule (13-19, run 2026-04-24).

Reuses the slide-cloning + 2-column layout pattern from
update_ppt_new_dimensions.py, with the post-consistency-fix typography:
  * Title: Cambria 16 pt explicit, white #FFFFFF (navy header band).
  * Tag:   Cambria  9 pt explicit, white #FFFFFF.
  * Body:  Cambria 10 pt explicit, navy #0F172A (matches dense-slide convention).
  * Sidebar root label: "6 / Gap-fill" (a new sub-category to distinguish
    these from the four existing post-viral PLUS sub-categories).

Run: python3 scripts/update_ppt_agenda_gaps.py
"""
from copy import deepcopy
from pathlib import Path

from pptx import Presentation
from pptx.util import Pt
from pptx.dml.color import RGBColor
from pptx.oxml.ns import qn

ROOT = Path("/Users/prempratapsingh/Data_Analysis")
PPT  = ROOT / "GRBV_Funding_Committee_Presentation.pptx"
FIGS = ROOT / "figures" / "supp"

TEMPLATE_INDEX = 27   # "Direct vs indirect temperature effects" 2-col layout
SECTION_LABEL  = "PLUS: agenda gap"
SECTION_ROOT   = "6 / Gap-fill"

NAVY  = RGBColor(0x0F, 0x17, 0x2A)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

NEW_SLIDES = [
    {
        "title": "PLUS (viral ORF): per-ORF prediction beyond qPCR titer + module coupling",
        "image": "S73_viral_orf_extensions.png",
        "body": [
            "Per-ORF expression as predictor of harvest juice phenotype on n = 27 RB(+) Harvest samples, fit with OLS + LRT against a titer-only baseline (Brix, TA, Malic).",
            "V2 (N761_gp1, silencing suppressor) adds significant beyond-titer signal for Malic (LRT p = 1.2e-3, beta = +271 mg/L per log2 CPM unit) and Brix (p = 0.033, beta = -0.37 degBx).  V1 (N761_gp2, coat) similar for Malic (p = 4.4e-4); V3 and C1 also reach significance for Malic.",
            "Per-ORF x 13-WGCNA-module Spearman correlation on 264 RB+ samples: 78 cells, BH-adjusted per ORF.  Operationalises the slide 51 framing that mRNA-level viral expression is decoupled from qPCR-DNA titer.",
            "Reading: viral mRNA load explains malic-acid variance that titer alone misses.  Direct answer to project Objective 2 ('impact of viral titer on disease outcomes'), and a candidate covariate for the next Brix/TA models.",
            "Caveat: N761_gp6 (C3 replication enhancer) is near-zero in our libraries (slide 49); its rows in the panel are uninformative as expected.",
        ],
    },
    {
        "title": "PLUS (per-gene path): phenylpropanoid + PR + cell-wall dissection",
        "image": "S74_pathway_per_gene.png",
        "body": [
            "Per-gene Infection x Season x Stage x Cohort OLS on three curated pathways called out in the project agenda.  Goes one level deeper than the existing pathway-level GSEA (S70).",
            "Tested universe (curated by GFF regex): Phenylpropanoid 51, PR proteins 40, Cell-wall 161 -- 252 genes total.  Filtered against the 14,604-gene background.",
            "Fraction with significant infection main effect (p < 0.05): Phenylpropanoid 31 % (16/51), Cell-wall 26 % (41/161), PR 20 % (8/40).  Phenylpropanoid is the most-disrupted of the three.",
            "Three-way interaction (Inf x Season x Stage) significant in 45-51 % of genes across all three pathways -- consistent with the inverted seasonal trajectory affecting these specific ripening axes.",
            "Reading: GRBV does not flatten these pathways uniformly; specific genes shift while others stay flat.  This is the granularity reviewers will look for, given the agenda explicitly names these three pathways.",
        ],
    },
    {
        "title": "PLUS (carbon translocation): SWEET / SUT / STP / CWINV module test",
        "image": "S75_carbon_translocation.png",
        "body": [
            "Tests the Martinez-Luscher 2018 (project ref [4]) hypothesis that GRBV blocks hexose translocation rather than primary carbon assimilation.  77 sugar-transport genes curated by GFF regex; 33 expressed in our 14,604-gene matrix.",
            "Subclasses present: SWEETs 0, SUTs 5, STPs/HXTs 30, cell-wall invertases 6, vacuolar invertases 0, sucrose synthases 6, SPS 4, hexokinase 5, trehalose 1, phloem markers 20.",
            "Module eigengene (PC1, 37.3 % VE) tracks Maturity_stage strongly (rho = 0.94, p = 2e-128).  No main season effect (rho = 0.01, p = 0.83) and no main cohort-age effect (rho = -0.07, p = 0.23).",
            "Reading: the sugar-transport panel is dominantly stage-driven, NOT virally rewired at the panel-level.  SPS1's mediator status (rule 08) and per-gene Brix correlations remain individually relevant; the panel-level Martinez-Luscher framing is not supported here without further filtering.",
            "Useful as a NULL counter-finding to anchor the SPS1 reframe ('compensatory upregulation, not rate-limiting step') against panel context.",
        ],
    },
    {
        "title": "PLUS (vine heterogeneity): tolerant vs susceptible Merlot vines within one vineyard",
        "image": "S76_vine_heterogeneity.png",
        "body": [
            "44 unique vines x up to 8 samples each (4 stages x 2 seasons) using vine_id from sample_metadata_v2.txt.  Per-vine response score = mean |log2FC| from cohort+season+stage matched RB(-) centroid on the top 2,000 variable genes.",
            "KMeans 3-cluster splits 44 vines into 13 tolerant / 18 intermediate / 13 susceptible.  All within one Merlot vineyard, controlling for cultivar genetic background.",
            "Phenotypic translation: tolerant vines reach Harvest Brix of 24.22 +/- 1.09 (n = 43 samples); susceptible vines 23.64 +/- 0.87 (n = 42).  ~0.6 Brix delta within the same vineyard year-equivalent.",
            "Implication: GRBV response is NOT genotype-uniform within Merlot.  Either rootstock, microenvironment, or transcriptional state heterogeneity matters.  This is publishable on its own -- a direct breeding lead if any of the tolerant vines have heritable differences.",
            "Caveat: 13 vines per cluster is small; the tolerance gradient is robust at the population level but per-vine assignments need replication.",
        ],
    },
    {
        "title": "PLUS (TF activity): regulon-based scoring re-ranks TF disruption",
        "image": "S77_tf_activity.png",
        "body": [
            "Per-sample TF activity scores derived from each TF's GENIE3 regulon (top 200 weighted targets per TF; rule 01 output).  Weighted-sum z-score on log2 CPM; Welch-t RB(+) vs RB(-) on the resulting activity vector per TF.",
            "236 TFs scored after filtering for >=5 expressed targets.  Activity-level DE re-ranks the family vs the expression-level DE that drives slide 36-37.",
            "Top hits by activity DE: LOC100262652 (uncharacterised, padj = 1.4e-64), LOC100246521 (padj = 7.6e-26), MYB30 (DOWN in RB+, padj = 8.9e-10), LOC100255447 (padj = 3.5e-7).  These are NOT the same hits as the expression-DE TF panel.",
            "MYB30 down-regulation in infected samples is consistent with classical PR-protein and HR-suppression pathways and complements the slide 27 'canonical immunity stays silent' framing.",
            "Caveat: decoupleR was unavailable in the FARM Rlibs; results use the lighter weighted-sum fallback.  If decoupleR is added, re-run for primary results -- the rank ordering should be similar but ulm gives more defensible activity calls.",
        ],
    },
    {
        "title": "PLUS (segmented): cohort age is a continuous gradient, not phase-transitions",
        "image": "S78_segmented_cohort.png",
        "body": [
            "Tests whether the deck's '3 temporal phases of chronic infection' (acute Yr 1 / peak Yrs 3-5 / chronic Yrs 6-7) reflects discrete phase transitions or a continuous gradient.  Fits segmented regression with k = 0/1/2 break-points on PC1, PC2, and 13 module eigengenes against years-since-infection on n = 264 RB(+) samples.",
            "Result: best_k = 0 (linear) for ALL 15 responses by BIC.  Delta BIC (k=1 minus k=0) is positive in every case, meaning the simpler linear model is preferred even after the break-point fits.",
            "Notable null finding: the deck's three-phase framing is a descriptive simplification, not a statistically supported phase architecture.  Reframe to 'cohort-age trajectory' (continuous) with 'short / mid / long-term cohort bins' as labels for clarity.",
            "Implication: mechanistic models do not need to invoke phase transitions; a smooth drift model is sufficient.  This simplifies the discussion and removes a target reviewers could push back on.",
            "Useful for reviewer-response: 'We tested for break-points formally and prefer the linear model; phase labels in the figures are descriptive only.'",
        ],
    },
    {
        "title": "PLUS (Bayesian): vine-random-intercept posteriors temper the MBF1c claim",
        "image": "S79_bayesian_posterior.png",
        "body": [
            "brms hierarchical model y ~ Infection + Season + Stage + (1 | vine_id) on 6 high-interest genes (top mediators + top biomarkers + top DOWN biomarker).  4 chains x 2000 iterations; weakly-informative N(0,1) priors on coefficients.  All converged (rhat <= 1.005, n_eff_min >= 746).",
            "Strong posteriors with 95 % CIs entirely outside ROPE [-0.1, +0.1]: RdRp1 (LOC100268011) +1.78 [+1.54, +2.02]; F-box LRR (LOC100854143) +2.01 [+1.76, +2.25]; glycine-rich (LOC100257913) -1.33 [-1.65, -1.01]; LOC100266309 (third robust mediator) +0.80 [+0.47, +1.11].",
            "Weaker posterior: SPS1 (LOC100232974) +0.52 [+0.09, +0.95], 2.3 % in ROPE -- credibly positive but small.",
            "Notable softening: MBF1c (LOC100249249) posterior is -0.25 [-0.81, +0.33] with 19 % of the posterior inside ROPE.  Straddles zero.  This is meaningfully more uncertain than the frequentist mediation sweep suggested.",
            "Reading: the headline mediator's effect needs the qPCR / VIGS validation -- the random-intercept Bayesian model already pushes back on the strength of the call.  Worth one cautious sentence in the Methods/Discussion before submission.",
        ],
    },
]


def set_text(tf, text, *, size=12, bold=False, color=None, name="Cambria"):
    tf.clear()
    p = tf.paragraphs[0]
    r = p.add_run(); r.text = text
    r.font.size = Pt(size); r.font.bold = bold; r.font.name = name
    if color:
        r.font.color.rgb = color if isinstance(color, RGBColor) else RGBColor(*color)


def find_shape(slide, name):
    for sh in slide.shapes:
        if sh.name == name:
            return sh
    return None


def clone_slide(prs, template_idx):
    src = prs.slides[template_idx]
    new = prs.slides.add_slide(src.slide_layout)
    for ph in list(new.placeholders):
        ph._element.getparent().remove(ph._element)
    skip = {qn("p:nvGrpSpPr"), qn("p:grpSpPr")}
    for child in src.shapes._spTree:
        if child.tag in skip:
            continue
        new.shapes._spTree.append(deepcopy(child))
    return new


def replace_picture(slide, new_image_path):
    for sh in slide.shapes:
        if sh.shape_type == 13:
            left, top, w, h = sh.left, sh.top, sh.width, sh.height
            sp = sh._element
            sp.getparent().remove(sp)
            slide.shapes.add_picture(str(new_image_path), left, top, width=w, height=h)
            return
    raise RuntimeError("no picture shape found in template")


def replace_body_text(slide, bullets):
    """Find largest non-structural Text box and refill with Cambria 10pt navy."""
    body = None
    for sh in slide.shapes:
        if sh.has_text_frame and sh.name.startswith("Text") and sh.name not in (
                "Text 1", "Text 3", "Text 6", "Text 7"):
            if body is None or sh.width * sh.height > body.width * body.height:
                body = sh
    if body is None:
        raise RuntimeError("no body text frame found")
    tf = body.text_frame
    tf.clear(); tf.word_wrap = True
    for i, line in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(4)
        r = p.add_run(); r.text = line
        r.font.size = Pt(10); r.font.name = "Cambria"
        r.font.color.rgb = NAVY


def update_title(slide, title):
    t = find_shape(slide, "Text 1")
    if t is not None:
        set_text(t.text_frame, title, size=16, bold=True, color=WHITE, name="Cambria")


def update_section_label(slide, label, root_label):
    lbl = find_shape(slide, "Text 3")
    if lbl is not None:
        set_text(lbl.text_frame, label, size=9, bold=True, color=WHITE, name="Cambria")
    root = find_shape(slide, "Text 6")
    if root is not None:
        set_text(root.text_frame, root_label, size=8, name="Cambria")


def reorder_slide(prs, from_idx, to_idx):
    xml_slides = prs.slides._sldIdLst
    slides = list(xml_slides)
    moving = slides[from_idx]
    xml_slides.remove(moving)
    xml_slides.insert(to_idx, moving)


def find_slide_by_title(prs, title_fragment):
    for i, s in enumerate(prs.slides):
        t = find_shape(s, "Text 1")
        if t is not None and title_fragment.lower() in t.text_frame.text.lower():
            return i
    return None


def append_reanalysis_bullet(prs):
    """Add an 'Agenda-gap pass added 7 analyses' bullet to slide 68."""
    idx = find_slide_by_title(prs, "What the reanalysis corrected")
    if idx is None:
        print("  warn: could not find 'What the reanalysis corrected' slide")
        return
    slide = prs.slides[idx]
    body = None
    for sh in slide.shapes:
        if sh.has_text_frame and "Universal DEGs" in sh.text_frame.text:
            body = sh; break
    if body is None:
        print("  warn: could not locate body text on reanalysis-corrected slide")
        return
    tf = body.text_frame
    p = tf.add_paragraph(); p.space_after = Pt(4)
    r = p.add_run()
    r.text = ("Agenda-gap extension pass (FARM rules 13-19, run 2026-04-24): "
              "viral ORFs add Brix/Malic prediction beyond titer; per-gene phenylpropanoid + PR + cell-wall "
              "(31 / 20 / 26 % infection-significant); 77-gene sugar-transport module is stage-driven not virus-driven; "
              "44 vines split 13 tolerant / 13 susceptible with 0.6 Bx delta; regulon-based TF activity re-ranks TFs "
              "(MYB30 down padj=8.9e-10); cohort age is a continuous gradient (best_k=0 for all 15 responses); "
              "Bayesian random-intercept posteriors confirm RdRp1, F-box, glycine-rich, LOC100266309 strongly, "
              "but soften MBF1c (95 % CI [-0.81, +0.33], 19 % in ROPE).  Figs S73-S79.")
    r.font.size = Pt(10); r.font.name = "Cambria"; r.font.color.rgb = NAVY
    print(f"  added agenda-gaps bullet to slide {idx + 1}")


def main():
    prs = Presentation(str(PPT))
    total_before = len(prs.slides)
    print(f"Opened presentation with {total_before} slides.")

    # 1. Extend "What the reanalysis corrected" BEFORE shuffling indices
    append_reanalysis_bullet(prs)

    # 2. Find the last "New dimension" slide; insert after it (i.e., before SECTION 7)
    anchor_idx = find_slide_by_title(prs, "New dimension 3")
    if anchor_idx is None:
        # fall back: anchor before SECTION 7 Synthesis
        anchor_idx = find_slide_by_title(prs, "Integrated disease cascade")
        if anchor_idx is not None:
            anchor_idx -= 1
    if anchor_idx is None:
        print("  warn: could not find anchor; appending at end")
        anchor_idx = len(prs.slides) - 1
    print(f"  inserting agenda-gap slides after slide {anchor_idx + 1}")

    added_idx = []
    for entry in NEW_SLIDES:
        new = clone_slide(prs, TEMPLATE_INDEX)
        update_title(new, entry["title"])
        update_section_label(new, SECTION_LABEL, SECTION_ROOT)
        img = FIGS / entry["image"]
        if not img.exists():
            raise FileNotFoundError(img)
        replace_picture(new, img)
        replace_body_text(new, entry["body"])
        added_idx.append(len(prs.slides) - 1)
        print(f"  cloned template -> '{entry['title'][:60]}...'")

    # Reorder cloned slides into position
    target_pos = [anchor_idx + 1 + i for i in range(len(NEW_SLIDES))]
    for i, _ in enumerate(added_idx):
        cur = len(prs.slides) - (len(added_idx) - i)
        reorder_slide(prs, cur, target_pos[i])

    # 3. Refresh slide-number footers (both Text 5 and Text 7 variants)
    total_after = len(prs.slides)
    for i, s in enumerate(prs.slides, start=1):
        for shape_name in ("Text 5", "Text 7"):
            num = find_shape(s, shape_name)
            if num is not None and "/" in num.text_frame.text:
                # preserve existing font-size by reading first run
                cur_size = None
                for p in num.text_frame.paragraphs:
                    for r in p.runs:
                        if r.font.size:
                            cur_size = r.font.size.pt; break
                    if cur_size: break
                cur_size = cur_size or 8
                set_text(num.text_frame, f"{i} / {total_after}",
                         size=cur_size, name="Cambria",
                         color=RGBColor(0x64, 0x74, 0x8B) if shape_name == "Text 5" else None)

    prs.save(str(PPT))
    print(f"Saved {PPT} with {total_after} slides (+{total_after - total_before}).")


if __name__ == "__main__":
    main()
