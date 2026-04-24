"""
Refresh the PPTX to reflect all figure upgrades and additions.

1. Replace stale embedded images on 13 slides with the current figures
   on disk under figures/main/ and figures/supp/.
2. Add a new slide 2 (after the title) pointing to the Summary_study_design
   graphical abstract so the deck opens with the unified project overview.
3. Add a new slide 5 (after 'GRBV disease biology') pointing to the
   MoA_GRBV_mechanism figure so the mechanism-of-action schematic sits in
   the biology-setup section rather than only in the cascade slide.

Run: python3 scripts/refresh_ppt_figures.py
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
FIGS_MAIN = ROOT / "figures" / "main"
FIGS_SUPP = ROOT / "figures" / "supp"

# Map from slide title fragment to figure file.  The updater looks up the
# slide by title fragment so index shifts don't break the match.
FIGURE_MAP = [
    ("Q1, Q2, Q6: DE landscape",                "main/Main02_de_landscape.png"),
    ("Q3, Q4, Q5: inverted seasonal",           "main/Main03_inverted_trajectory.png"),
    ("PLUS (heat-stress): sHSP suppression",    "main/Main04_heat_stress.png"),
    ("Q4 mechanism: causal mediation",          "main/Main05_mediation.png"),
    ("Direct vs indirect temperature effects",  "supp/S50_direct_vs_indirect_effects.png"),
    ("PLUS (biomarker): ML panel",              "main/Main06_biomarker.png"),
    ("Q2 (network view): WGCNA modules",        "main/Main07_wgcna_tf.png"),
    ("PLUS (GWAS / QTL): chromosome hotspots",  "main/Main08_chr_qtl.png"),
    ("PLUS (biochemistry): juice chemistry",    "supp/S63_biochem_season_infection.png"),
    ("Cohort trajectory at harvest",            "supp/S62_biochem_cohort_trajectory.png"),
    ("SPS1 reframed",                           "supp/S61_biochem_validation.png"),
    ("Tier 1 log2FC concordance",               "supp/S65_tier1_fc_concordance.png"),
    ("Synthesis: integrated molecular cascade", "main/Main10_biology_cascade.png"),
]

# Retitle + reframe some slides whose narrative the new analyses have shifted
TITLE_REFRAMES = {
    "PLUS (heat-stress): sHSP suppression and HSF coordination":
        "PLUS (heat-stress): HSF to sHSP axis dampened (passive failure -- see slides 54-55)",
    "Q4 mechanism: causal mediation via viral titer":
        "Q4 mechanism: Temp -> Titer -> Gene mediation (MBF1c, SPS1, LOC100266309)",
}

TEMPLATE_INDEX = 27   # 2-col layout template used for new slides
SECTION_LABEL_OVERVIEW = "OVERVIEW"
SECTION_LABEL_BIOLOGY  = "SECTION 1 / Biology"
SECTION_ROOT = "0 / Graphical"

NEW_FRONTMATTER_SLIDES = [
    {
        # Insert at position 2 (right after title slide 0, keep Executive 1)
        "after_title_fragment": "Executive summary",
        "title":   "Study at a glance",
        "section": SECTION_LABEL_OVERVIEW,
        "root":    "0 / Overview",
        "image":   "main/Summary_study_design.png",
        "body": [
            "Vitis vinifera cv. Merlot, Oakville CA, 2021 cool season + 2022 hot season, 4 maturity stages, 7 infection-onset cohorts (2016-2022), 327 RNA-Seq libraries.",
            "Workflow: field samples -> TAG-Seq library prep (B. N. Jenner, UCD) -> covariate-adjusted DESeq2 + limma-voom -> integrated multi-omics analysis (13 FARM pipeline rules) -> disease phenotype.",
            "Integrated findings hub (station 4): twin geminate virion + sHSP chaperone ring + DOF-TF motif + sugar hexagon + cracked shield (immunity dampened) + red-blotched leaf.",
            "Headline outcomes: inverted seasonal trajectory, sHSP axis dampened, MBF1c + SPS1 mediators, GRBV signature replicates in Zinfandel (Blanco-Ulate 2017).",
        ],
    },
    {
        # Insert after the "GRBV disease biology" slide
        "after_title_fragment": "GRBV disease biology",
        "title":   "Mechanism of action schematic (vector -> phloem -> host cell)",
        "section": SECTION_LABEL_BIOLOGY,
        "root":    "1 / Biology",
        "image":   "main/MoA_GRBV_mechanism.png",
        "body": [
            "Three-zone scientific scene summarising the known GRBV pathway at the tissue, cellular, and molecular levels.",
            "Zone A (left): three-cornered alfalfa treehopper feeds on a healthy Merlot vine and delivers GRBV to the vascular tissue. Transmission contact shown as a small red dot at the feeding site.",
            "Zone B (middle): twinned geminate virion (two fused icosahedral capsids, GRBV family signature) enters phloem sieve elements.",
            "Zone C (right): rolling-circle replication of the circular ssDNA genome in the host cell nucleus. Downstream effects shown as dashed arrows: dampened immune shield, suppressed sHSP chaperone ring, rewired DNA + transcription factor with sugar glucose, and emerging red-blotched leaf as phenotype.",
            "Silent narrative: the full chain from vector to leaf symptom in one frame, matching the mechanism-cascade of slide 58.",
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


def find_slide_by_title(prs, title_fragment):
    frag = title_fragment.lower()
    for i, s in enumerate(prs.slides):
        for sh in s.shapes:
            if sh.has_text_frame and sh.text_frame.text:
                if frag in sh.text_frame.text.lower():
                    return i
    return None


def get_slide_title_shape(slide):
    t = find_shape(slide, "Text 1")
    if t is not None:
        return t
    for sh in slide.shapes:
        if sh.has_text_frame and sh.text_frame.text:
            return sh
    return None


def replace_picture(slide, new_image_path):
    for sh in slide.shapes:
        if sh.shape_type == MSO_SHAPE_TYPE.PICTURE:
            left, top, w, h = sh.left, sh.top, sh.width, sh.height
            sp = sh._element
            sp.getparent().remove(sp)
            slide.shapes.add_picture(str(new_image_path), left, top, width=w, height=h)
            return True
    return False


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


def retitle_slides(prs):
    n = 0
    for old, new in TITLE_REFRAMES.items():
        idx = find_slide_by_title(prs, old)
        if idx is None:
            print(f"  retitle skipped (not found): {old[:60]}")
            continue
        t = get_slide_title_shape(prs.slides[idx])
        if t is not None:
            set_text(t.text_frame, new, size=16, bold=True)
            n += 1
            print(f"  retitled slide {idx+1}: {new[:70]}")
    print(f"  total retitles: {n}")


def refresh_figures(prs):
    n_ok = 0
    n_missing = 0
    for fragment, rel_path in FIGURE_MAP:
        img = ROOT / "figures" / rel_path
        if not img.exists():
            print(f"  MISSING on disk: {img}")
            n_missing += 1
            continue
        idx = find_slide_by_title(prs, fragment)
        if idx is None:
            print(f"  slide not found for: {fragment}")
            n_missing += 1
            continue
        if replace_picture(prs.slides[idx], img):
            print(f"  refreshed slide {idx+1}: {rel_path}")
            n_ok += 1
        else:
            print(f"  no picture shape on slide {idx+1} (add manually)")
            n_missing += 1
    print(f"  refreshed {n_ok} figures, {n_missing} unresolved")


def insert_new_slides(prs):
    # We insert in descending position order so earlier inserts don't shift
    # later insertion points.
    for entry in sorted(NEW_FRONTMATTER_SLIDES,
                        key=lambda e: -(find_slide_by_title(prs, e["after_title_fragment"]) or 0)):
        anchor_idx = find_slide_by_title(prs, entry["after_title_fragment"])
        if anchor_idx is None:
            print(f"  skip new slide (anchor not found): {entry['title'][:60]}")
            continue
        new = clone_slide(prs, TEMPLATE_INDEX)
        update_title(new, entry["title"])
        update_section_label(new, entry["section"], entry["root"])
        img = ROOT / "figures" / entry["image"]
        if not img.exists():
            print(f"  skip new slide (image missing): {img}")
            # still remove the appended clone to avoid dangling
            xml = prs.slides._sldIdLst
            xml.remove(list(xml)[-1])
            continue
        replace_picture(new, img)
        replace_body_text(new, entry["body"])
        target_pos = anchor_idx + 1
        cur = len(prs.slides) - 1
        reorder_slide(prs, cur, target_pos)
        print(f"  inserted new slide at position {target_pos+1}: {entry['title'][:60]}")


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

    print("\n[1/4] Refreshing stale embedded figures...")
    refresh_figures(prs)

    print("\n[2/4] Reframing slide titles for new narrative...")
    retitle_slides(prs)

    print("\n[3/4] Inserting Summary + MoA front-matter slides...")
    insert_new_slides(prs)

    print("\n[4/4] Refreshing slide-number footers...")
    refresh_slide_numbers(prs)

    prs.save(str(PPT))
    print(f"\nSaved {PPT} with {len(prs.slides)} slides (+{len(prs.slides) - before}).")


if __name__ == "__main__":
    main()
