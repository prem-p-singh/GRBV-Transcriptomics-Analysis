"""
Add slides for the three new scientific dimensions unlocked by S66-S72:
  1. Passive-failure sHSP mechanism (S66 + S67 comparison -> passive failure)
  2. Three temporal phases of chronic infection (S69 + S71)
  3. Disease cost is temperature-amplified 2.6x at harvest phenotype (S66 panel D + S67 panel C)

The 4 new slides are inserted after the existing Tier 1 block (the last PLUS
block before SECTION 7 / Synthesis).  The "What the reanalysis corrected"
slide is extended with a bullet covering the new dimensions.

Reuses the helper functions from update_ppt_tier1.py.

Run: python3 scripts/update_ppt_new_dimensions.py
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
SECTION_LABEL  = "PLUS: new dimensions"
SECTION_ROOT   = "6 / Results"

NEW_SLIDES = [
    {
        "title": "New dimension 1: sHSP effect is passive failure, not active suppression",
        "image": "S66_disease_progression_cube.png",
        "body": [
            "Disease progression cube for RB(+) vines across cohort (2016-2022 onset) x sampling season (2021, 2022) x maturity stage (Pre-ver, Ver, Post-ver, Harvest), stacked as 4 layers: viral titer, sHSP module z-score, GRBV 55-gene signature, harvest juice Brix.",
            "sHSP panel (B) at the 2022 Harvest cell shows a mean RB(+) z-score of -0.67 -- the expected sHSP suppression finding.",
            "But the companion figure S67 shows uninfected RB(-) vines at the same 2022 Harvest cell reach +1.64 z -- a strong natural heat-driven upregulation.  The GRBV effect is therefore a 2.3-unit DEFICIT from healthy baseline, not an active downregulation from neutral.",
            "Gene-level variance decomposition S72 confirms: sHSP genes carry environmental variance share 0.33 (vs 0.19 baseline) but viral variance share only 0.11, indistinguishable from the 0.12 baseline across 14,604 genes.  GRBV does not drive sHSP transcription; it blocks the environmental programme from running.",
            "Reframing: wet-lab target shifts from 'what represses sHSPs' to 'what blocks HSF / MBF1c transactivation under heat' -- directly strengthens the planned HSFA2 VIGS experiment.",
        ],
    },
    {
        "title": "New dimension 1 (cont.): RB(-) baseline is the comparator that makes passive failure visible",
        "image": "S67_rbminus_baseline_cube.png",
        "body": [
            "RB(-) baseline cube: same grid as S66 but uninfected controls, so cohort axis collapses to a single control group per season x stage cell.",
            "Panel A (sHSP module z-score across stage x season): RB(-) vines reach +1.64 z at 2022 Harvest -- the single strongest cell in the entire control grid.  This is the healthy heat-stress response operating normally.",
            "Panel B (GRBV 55-gene signature z-score): RB(-) vines mostly sit at -1.0 to -1.4 as expected -- the GRBV signature is specific to infection, not a general plant stress response.",
            "Panel C (harvest juice Brix): RB(-) 26.7 Bx in cool 2021, 34.4 Bx in hot 2022 -- normal vines ripen sugar normally in both seasons.",
            "Direct comparison of S66 panel B and S67 panel A at the same cell (2022 Harvest) establishes the passive-failure interpretation: healthy vines switch on sHSPs under heat, infected vines do not.  The 2.3-unit deficit is the most reliable GRBV phenotype in the dataset.",
        ],
    },
    {
        "title": "New dimension 2: chronic infection has three temporal phases",
        "image": "S71_transcriptional_response_strength.png",
        "body": [
            "PCA-based response-strength metric: for each RB(+) sample, Euclidean distance in PC1-PC2 space (72.9 % of variance on top-2,000 variable genes) from the centroid of matched-condition (same season x stage) RB(-) controls.  Isolates infection drift from the dominant stage and season variation.",
            "Acute phase (2022 cohort, year 1): 214 DEGs -- the largest count -- but the weakest pathway enrichment (<=5x fold for most curated families, S69).  Broad and shallow.",
            "Peak phase (2020 cohort, years 3-5): PC-space distance reaches 21.9 at 2021 Harvest -- the top-ranked cell across all 40 cohort x season x stage cells.  Sharpest transcriptome drift from healthy baseline.",
            "Chronic phase (2016 cohort, years 6-7): only 105 DEGs but hyper-focused enrichment: 65x sHSP, 138x HSP100/ClpB, 49x GO:response to heat.  Narrow and intense.",
            "GRBV 55-gene signature is the one family enriched 52-100x in every cohort -- the shared fingerprint that persists through all three phases.  No prior GRBV study has had the cohort span to see this temporal architecture.",
        ],
    },
    {
        "title": "New dimension 3: disease cost is temperature-amplified 2.6x at harvest",
        "image": "S66_disease_progression_cube.png",
        "body": [
            "Harvest juice Brix from the biochemistry batch, cross-referenced between RB(+) and RB(-) samples within the same season.",
            "Cool 2021 season: RB(-) control vines reach 26.7 Bx; RB(+) infected vines reach 22.8-23.4 Bx across 5 infection-onset cohorts.  Disease cost ~ 4 Bx.",
            "Hot 2022 season: RB(-) control vines reach 34.4 Bx; RB(+) infected vines reach 23.5-25.0 Bx across 5 cohorts.  Disease cost ~ 10 Bx.",
            "Hot-year yield cost is 2.6x larger than cool-year cost -- a direct crop-value translation of the transcriptional 'inverted seasonal trajectory' finding.",
            "Combined with the passive-failure mechanism: the hot-year environment is where healthy vines mount the biggest thermal defence, and therefore the hot year is where infected vines fail the hardest.  The temperature-amplified yield penalty is the phenotypic signature of that mechanism.",
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
        r.font.size = Pt(10); r.font.name = "Arial"


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


def find_slide_by_title(prs, title_fragment):
    for i, s in enumerate(prs.slides):
        t = find_shape(s, "Text 1")
        if t is not None and title_fragment.lower() in t.text_frame.text.lower():
            return i
    return None


def append_reanalysis_bullet(prs):
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
    r.text = ("Three new dimensions from the cohort x season x stage cube: "
              "(1) sHSP effect reframed as passive failure to activate under heat "
              "(S66 -0.67 z vs S67 +1.64 z at 2022 Harvest; S72 viral share 0.11 ~ baseline 0.12); "
              "(2) three temporal phases of chronic infection visible across 7-year cohort span "
              "(acute 2022 / peak 2020 / chronic 2016; S69 + S71); "
              "(3) hot-year disease cost is 2.6x larger than cool-year cost "
              "(10 Bx vs 4 Bx; S66 panel D, S67 panel C).")
    r.font.size = Pt(10); r.font.name = "Arial"
    print(f"  added new-dimensions bullet to slide {idx + 1}")


def main():
    prs = Presentation(str(PPT))
    total_before = len(prs.slides)
    print(f"Opened presentation with {total_before} slides.")

    # 1. Extend "What the reanalysis corrected" BEFORE shuffling indices
    append_reanalysis_bullet(prs)

    # 2. Find the Tier 1 concordance slide (end of Tier 1 block); insert after
    anchor_idx = find_slide_by_title(prs, "Tier 1 log2FC concordance")
    if anchor_idx is None:
        anchor_idx = find_slide_by_title(prs, "SPS1 reframed")
    if anchor_idx is None:
        print("  warn: could not find anchor slide; appending at end")
        anchor_idx = len(prs.slides) - 1
    print(f"  inserting new-dimension slides after slide {anchor_idx + 1}")

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
        print(f"  cloned template -> new slide '{entry['title'][:60]}...'")

    target_pos = [anchor_idx + 1 + i for i in range(len(NEW_SLIDES))]
    for i, _ in enumerate(added_idx):
        cur = len(prs.slides) - (len(added_idx) - i)
        reorder_slide(prs, cur, target_pos[i])

    # 3. Refresh slide-number footers
    total_after = len(prs.slides)
    for i, s in enumerate(prs.slides, start=1):
        num = find_shape(s, "Text 7")
        if num is not None and "/" in num.text_frame.text:
            set_text(num.text_frame, f"{i} / {total_after}", size=8)

    prs.save(str(PPT))
    print(f"Saved {PPT} with {total_after} slides (+{total_after - total_before}).")


if __name__ == "__main__":
    main()
