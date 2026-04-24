"""
Add 5 slides for the new analyses that lack dedicated PPT coverage:
  S29 - DE method concordance heatmap (Welch-t / DESeq2 / voom / RM-voom)
  S68 - Season x stage enrichment with season-flip catalog
  S69 - Per-cohort pathway hyperenrichment (2016 chronic vs 2022 acute)
  S70 - Literature-curated GRBV pathway shifts (Sudarshana, Blanco-Ulate, etc.)
  S72 - Per-gene variance ternary (Singh 2025 bioRxiv framework adapted)

Each slide goes into the natural section of the deck:
  S29 after 'Statistical design' (Section 2)
  S68 after 'GSEA NES heatmap, full gene-set panel'
  S69 before 'New dimension 2: chronic infection has three temporal phases'
  S70 after 'Master pathway landscape'
  S72 after 'Variance decomposition, per-gene detail'

Run: python3 scripts/add_ppt_new_analyses.py
"""
from copy import deepcopy
from pathlib import Path

from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE
from pptx.util import Pt
from pptx.dml.color import RGBColor
from pptx.oxml.ns import qn

ROOT = Path("/Users/prempratapsingh/Data_Analysis")
PPT  = ROOT / "GRBV_Funding_Committee_Presentation.pptx"

TEMPLATE_INDEX = 29           # 2-col layout template (still at this index)
SECTION_ROOT_METHODS = "2 / Methods"
SECTION_ROOT_RESULTS = "6 / Results"

NEW_SLIDES = [
    {
        "after": "Statistical design",
        "title": "PLUS (DE rigor): 4-method concordance across every contrast",
        "section": "PLUS: method rigor",
        "root": SECTION_ROOT_METHODS,
        "image": "supp/S29_comprehensive_summary_dashboard.png",
        "body": [
            "Primary DE pipeline uses DESeq2 with covariate adjustment (~ Maturity_stage + Sampling_Season + Infection_Group). limma-voom + repeated-measures limma-voom run as sensitivity checks. Welch-t on log2(CPM+1) is the originally-published comparator.",
            "Q6 pooled RB+ vs RB-: DESeq2 = 126 DEGs; voom = 150; RM-voom = 150; Welch-t reported ~908. Jaccard(DESeq2, voom) = 0.62. The published 908-gene 'universal' set was inflated by stage and season confounding that covariate adjustment removes.",
            "Q2 stage-stratified contrasts all converge across the 3 covariate-adjusted methods: Pre-veraison 105/111/105, Veraison 124/130/124, Post-veraison 532/541/529, Harvest 964/1025/964. Q2 Harvest is the single largest infection signal.",
            "Repeated-measures voom (duplicateCorrelation on vine_id, sensitivity check only) shows within-vine correlation is negligible (|rho-hat| <= 0.04 everywhere). Pseudoreplication is not a practical concern; de_deseq2_voom.R is the primary result.",
            "Methods-level upshot: the paper's headline counts shift (~130 covariate-adjusted core Q6 DEGs, ~1,025 at Harvest), but directional story is preserved under all three modern methods.",
        ],
    },
    {
        "after": "Variance decomposition, per-gene detail",
        "title": "PLUS (variance ternary): gene-level decomposition (Singh 2025 framework)",
        "section": "PLUS: variance ternary",
        "root": SECTION_ROOT_RESULTS,
        "image": "supp/S72_variance_decomposition_ternary.png",
        "body": [
            "Conceptual replica of Singh et al. 2025 bioRxiv Fig 5 (Botrytis host x isolate variance partitioning), adapted to our single-species design with axes: Developmental = stage eta^2, Viral = infection + interactions with infection, Environmental = season + season:stage interactions.",
            "Panel A (hex-bin count on log10 scale): 14,604 genes with >2% total variance explained. Density peaks at the Dev apex - most genes are stage-driven, consistent with stage eta^2 = 0.426.",
            "Panel B (biological categories overlaid): GRBV 55-gene signature pulls toward the Viral apex; sHSPs, HSFs, and WRKY TFs cluster along the Dev-Env edge rather than the Viral apex.",
            "Panel C (per-category box plots): GRBV signature carries viral variance share 0.33 (vs baseline 0.12) - 2.8x the all-genes background. sHSPs carry environmental share 0.33 vs baseline 0.19 with viral share 0.11 ~ baseline 0.12 - viral contribution indistinguishable from noise.",
            "Interpretation: sHSPs are environment-driven, not virus-driven. The observed sHSP 'suppression' reflects failure to mount the normal heat-induced programme, not active GRBV-driven downregulation. Quantitative support for the passive-failure reframe.",
        ],
    },
    {
        "after": "GSEA NES heatmap, full gene-set panel",
        "title": "PLUS (enrichment, season): season x stage GSEA + season-flip catalog",
        "section": "PLUS: enrichment shifts",
        "root": SECTION_ROOT_RESULTS,
        "image": "supp/S68_enrichment_season_stage.png",
        "body": [
            "Three-panel companion to the existing NES heatmap (slide 39), using the same 375-gene-set GSEA output but visualised in three ways that reveal structure not in the original plot.",
            "Panel A (top 30 most variable pathways): rows clustered on 8-condition NES profile; two clean groups emerge - 'suppressed everywhere' (ribosome, translation) vs 'season-dependent' (stress, defense blocks that flip sign between 2021 and 2022).",
            "Panel B (season-flip catalog): 20 pathways where NES reverses sign between the two seasons at the same maturity stage. This is the pathway-level measurable equivalent of the inverted seasonal trajectory.",
            "Panel C (trajectory lines, 8 key families): ribosome and translation drop across stages in both seasons but more severely in 2022. Heat response drops to NES = -2.5 at 2022 Harvest. Sugar, flavonoid, photosynthesis all show the 2021-vs-2022 inversion; HSP/chaperone and ribosome stay suppressed across both seasons.",
            "Reading: GRBV effects split into 'season-independent core suppression' (HSP, ribosome) and 'season-dependent ripening disruption' (sugar, flavonoid, photosynthesis).",
        ],
    },
    {
        "after": "Master pathway landscape",
        "title": "PLUS (enrichment, literature): GRBV-specific pathway shifts (10 curated categories)",
        "section": "PLUS: literature-curated",
        "root": SECTION_ROOT_RESULTS,
        "image": "supp/S70_literature_curated_pathways.png",
        "body": [
            "Filtered the 375-pathway GSEA panel to 90 gene sets matching 10 categories reported in the GRBV / Vitis-virus literature (Sudarshana 2015, Blanco-Ulate 2017, Yepes 2018, Chitarra 2018, Agudelo-Romero 2015, Song 2022). Companion CSV at repo root has the top-4 per category.",
            "Ribosome / translation: strongly DOWN in every cell (-1.55 to -1.84) - the single most consistent GRBV signature in our data.",
            "Heat shock / chaperone: DOWN, intensifying through stages, maxes at -1.94 NES in 2022 Harvest.",
            "Sugar / starch + photosynthesis + phenylpropanoid / flavonoid: all show the 2021-vs-2022 inversion (positive late in 2022, flat or negative in 2021). These are the season-dependent ripening-disruption pathways.",
            "Defense / immunity: transient UP at 2021 Pre-ver + Ver (+0.70, +1.09 NES) then collapses to baseline after veraison. Hormone signaling: no strong season-stage trajectory. Sits alongside the 10-defense-category silence claim on slide 24.",
            "Reframing: GRBV effects split cleanly into 'core suppression' (HSP, ribosome) and 'ripening disruption' (sugar, flavonoid, photosynthesis).",
        ],
    },
    {
        "after": "Tier 1 log2FC concordance",
        "title": "PLUS (enrichment, cohort): per-cohort pathway hyperenrichment reveals chronic phase",
        "section": "PLUS: cohort enrichment",
        "root": SECTION_ROOT_RESULTS,
        "image": "supp/S69_enrichment_cohort.png",
        "body": [
            "Fisher-exact enrichment of per-cohort DEG lists (Q1b DESeq2) against 13 curated gene-set families: 9 heat-stress categories, 3 TF families, and the GRBV 55-gene signature. Universe = all tested genes per cohort.",
            "Panel A (cohort x family heatmap): 2016 cohort (year 6-7) reaches 65x for sHSP (17 overlap), 138x for HSP100/ClpB (1 overlap), 49x for GO:response to heat (15 overlap), 46x for HSP90. 2022 cohort (year 1) only reaches 15x sHSP, 11x HSP90, 5x HSF.",
            "GRBV 55-gene signature is enriched 52-100x in every cohort - the persistent transcriptomic fingerprint across all three temporal phases.",
            "Panel B (per-cohort DEG counts): monotonic from 105 in 2016 to 214 in 2022. The newest cohort has the most DEGs but the weakest enrichment.",
            "Inverse relationship: enrichment intensity rises as DEG breadth falls. Older infections focus narrowly on heat / chaperone pathways; newer infections scatter broadly across the transcriptome without hitting any pathway hard. Direct support for the three-phase temporal architecture (slide 57).",
            "No prior GRBV study has had the 7-year cohort span to see this. Pairs with slide 57 as the enrichment counterpart to the PC-space response-strength view.",
        ],
    },
]


def set_text(tf, text, size=12, bold=False, color=None):
    tf.clear()
    p = tf.paragraphs[0]
    r = p.add_run(); r.text = text
    r.font.size = Pt(size); r.font.bold = bold
    if color:
        r.font.color.rgb = RGBColor(*color)


def find_shape(slide, name):
    for sh in slide.shapes:
        if sh.name == name:
            return sh
    return None


def find_slide_by_title(prs, fragment):
    frag = fragment.lower()
    for i, s in enumerate(prs.slides):
        for sh in s.shapes:
            if sh.has_text_frame and sh.text_frame.text:
                if frag in sh.text_frame.text.lower():
                    return i
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
        if sh.shape_type == MSO_SHAPE_TYPE.PICTURE:
            left, top, w, h = sh.left, sh.top, sh.width, sh.height
            sp = sh._element
            sp.getparent().remove(sp)
            slide.shapes.add_picture(str(new_image_path), left, top, width=w, height=h)
            return True
    return False


def replace_body_text(slide, bullets):
    body = None
    for sh in slide.shapes:
        if sh.has_text_frame and sh.name.startswith("Text") and sh.name not in (
                "Text 1", "Text 3", "Text 6", "Text 7"):
            if body is None or sh.width * sh.height > body.width * body.height:
                body = sh
    if body is None:
        return False
    tf = body.text_frame
    tf.clear(); tf.word_wrap = True
    for i, line in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(4)
        r = p.add_run(); r.text = line
        r.font.size = Pt(10); r.font.name = "Arial"
    return True


def update_title(slide, title):
    t = find_shape(slide, "Text 1")
    if t is not None:
        set_text(t.text_frame, title, size=16, bold=True)


def update_section_label(slide, label, root_label):
    lbl = find_shape(slide, "Text 3")
    if lbl is not None:
        set_text(lbl.text_frame, label, size=9, bold=True, color=(255, 255, 255))
    root = find_shape(slide, "Text 6")
    if root is not None:
        set_text(root.text_frame, root_label, size=8)


def reorder_slide(prs, from_idx, to_idx):
    xml_slides = prs.slides._sldIdLst
    slides = list(xml_slides)
    moving = slides[from_idx]
    xml_slides.remove(moving)
    xml_slides.insert(to_idx, moving)


def refresh_slide_numbers(prs):
    total = len(prs.slides)
    for i, s in enumerate(prs.slides, start=1):
        num = find_shape(s, "Text 7")
        if num is not None and "/" in num.text_frame.text:
            set_text(num.text_frame, f"{i} / {total}", size=8)


def main():
    prs = Presentation(str(PPT))
    before = len(prs.slides)
    print(f"Opened presentation with {before} slides.")

    # Insert from BOTTOM of deck upward so earlier inserts don't shift
    # later anchor positions.  Compute all positions first.
    positions = []
    for entry in NEW_SLIDES:
        idx = find_slide_by_title(prs, entry["after"])
        if idx is None:
            print(f"  WARN: anchor not found: {entry['after']}")
            continue
        positions.append((idx + 1, entry))

    # Sort by target position descending
    positions.sort(key=lambda p: -p[0])

    for target_pos, entry in positions:
        new = clone_slide(prs, TEMPLATE_INDEX)
        update_title(new, entry["title"])
        update_section_label(new, entry["section"], entry["root"])
        img = ROOT / "figures" / entry["image"]
        if not img.exists():
            print(f"  MISSING image: {img}")
            xml = prs.slides._sldIdLst
            xml.remove(list(xml)[-1])
            continue
        replace_picture(new, img)
        replace_body_text(new, entry["body"])
        cur = len(prs.slides) - 1
        reorder_slide(prs, cur, target_pos)
        print(f"  inserted at position {target_pos+1}: {entry['title'][:65]}")

    refresh_slide_numbers(prs)
    prs.save(str(PPT))
    after = len(prs.slides)
    print(f"\nSaved {PPT} with {after} slides (+{after - before}).")


if __name__ == "__main__":
    main()
