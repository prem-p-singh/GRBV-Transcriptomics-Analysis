"""
Replace the single 'PLUS (per-gene, full)' slide (S80) with two dedicated slides:
  - S80A: Phenylpropanoid + PR proteins (91 genes, 2 stacked panels)
  - S80B: Cell wall (161 genes, sub-categorised, 5 functional sub-categories)

Run: python3 scripts/update_ppt_split_S80.py
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

TEMPLATE_INDEX = 27
SECTION_LABEL  = "PLUS: agenda gap"
SECTION_ROOT   = "6 / Gap-fill"

NAVY  = RGBColor(0x0F, 0x17, 0x2A)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
GREY  = RGBColor(0x64, 0x74, 0x8B)

OLD_SLIDE_TITLE_FRAGMENT = "PLUS (per-gene, full)"

NEW_SLIDES = [
    {
        "title": "PLUS (per-gene, S80A): Phenylpropanoid + PR proteins log2 FC across season x cohort x stage",
        "image": "S80A_phenylpropanoid_PR_heatmap.png",
        "body": [
            "91 genes total (Phenylpropanoid 51 + PR proteins 40), shown in two stacked panels.  Companion to S74 summary view.  Cell value = mean log2(CPM+1) RB(+) at (Sampling_Season, Maturity_stage, Years-Since-Infection cohort) - mean log2(CPM+1) RB(-) at matched (Sampling_Season, Maturity_stage).",
            "Layout fixes: facet_nested season > stage so strip labels do not overlap; 2-digit cohort years on x-axis ('16, '19, '20, '21, '22) at 90 deg; gene labels on y at 6 pt for full readability over 14-inch height.",
            "Row order: hierarchical clustering (Euclidean / Ward.D2) per panel on the full 8-stage-season x 5-cohort log2 FC matrix.  Genes with similar response trajectories cluster together.",
            "Color scale: diverging blue / white / red, capped at +/-3 log2 units via scale_fill_gradient2(limits = c(-3, 3), oob = squish) so extreme cells do not blow out midrange contrast.",
            "Cohort convention: '16 = 5-6 yrs since infection at 2021 sampling, 9-10 yrs at 2022 sampling.  '22 = 0-1 yr at sampling.",
        ],
    },
    {
        "title": "PLUS (per-gene, S80B): Cell-wall genes log2 FC, sub-categorised functionally",
        "image": "S80B_cell_wall_heatmap.png",
        "body": [
            "161 cell-wall genes split into 5 functional sub-categories: Pectin (PME / PG / pectinesterase / pectate lyase), Hemicellulose (XTH / xyloglucan / beta-galactosidase), Cellulose (CESA / CSL / endoglucanase / COBRA), Lignin (laccase / CCR / CAD / CCoAOMT / peroxidase), Expansin.  Sub-category assignment by GFF description regex.",
            "Sub-categories rendered as facet rows via facet_nested, separating functional groups for biological interpretation.  Each sub-category cluster ordered independently by Ward.D2 hierarchical clustering on Euclidean distance.",
            "Gene labelling threshold: only the 13 genes with |log2FC| >= 1 in at least 4 of 8 (season x stage) outer cells receive y-axis labels.  Other 148 genes shown as anonymous tile rows for the global pattern -- prevents the 161 LOC IDs from overlapping into illegibility.",
            "Layout: 11 x 18 inch figure at 600 DPI.  Cohort years '16/'19/'20/'21/'22 on x-axis, season > stage facets on top, sub-category facets on left.",
            "Reading: most cell-wall genes show muted infection effects (pale tiles); the labelled subset highlights the genes that consistently shift across multiple stage-season cells.  Use this for gene-by-gene reviewer queries on cell-wall remodelling.",
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


def find_slide_by_title(prs, title_fragment):
    for i, s in enumerate(prs.slides):
        t = find_shape(s, "Text 1")
        if t is not None and title_fragment.lower() in t.text_frame.text.lower():
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
        r.font.size = Pt(10); r.font.name = "Cambria"
        r.font.color.rgb = NAVY


def update_title(slide, title):
    t = find_shape(slide, "Text 1")
    if t is not None:
        set_text(t.text_frame, title, size=16, bold=True, color=WHITE)


def update_section_label(slide, label, root_label):
    lbl = find_shape(slide, "Text 3")
    if lbl is not None:
        set_text(lbl.text_frame, label, size=9, bold=True, color=WHITE)
    root = find_shape(slide, "Text 6")
    if root is not None:
        set_text(root.text_frame, root_label, size=8)


def reorder_slide(prs, from_idx, to_idx):
    xml_slides = prs.slides._sldIdLst
    slides = list(xml_slides)
    moving = slides[from_idx]
    xml_slides.remove(moving)
    xml_slides.insert(to_idx, moving)


def remove_slide(prs, idx):
    xml_slides = prs.slides._sldIdLst
    slides = list(xml_slides)
    xml_slides.remove(slides[idx])


def main():
    prs = Presentation(str(PPT))
    total_before = len(prs.slides)
    print(f"Opened presentation with {total_before} slides.")

    # 1. Find and remove the old single S80 slide
    old_idx = find_slide_by_title(prs, OLD_SLIDE_TITLE_FRAGMENT)
    if old_idx is None:
        raise SystemExit(f"Could not find slide containing '{OLD_SLIDE_TITLE_FRAGMENT}'")
    print(f"  removing old slide {old_idx + 1}")
    remove_slide(prs, old_idx)

    # 2. Insert two new slides at old_idx (after slide 66, the per-gene-path summary)
    insertions = []
    for entry in NEW_SLIDES:
        new = clone_slide(prs, TEMPLATE_INDEX)
        update_title(new, entry["title"])
        update_section_label(new, SECTION_LABEL, SECTION_ROOT)
        img = FIGS / entry["image"]
        if not img.exists():
            raise FileNotFoundError(img)
        replace_picture(new, img)
        replace_body_text(new, entry["body"])
        insertions.append(len(prs.slides) - 1)
        print(f"  cloned -> '{entry['title'][:60]}...'")

    # Reorder cloned slides (which were appended at end) to old_idx and old_idx+1
    target_pos = [old_idx + i for i in range(len(NEW_SLIDES))]
    for i, _ in enumerate(insertions):
        cur = len(prs.slides) - (len(insertions) - i)
        reorder_slide(prs, cur, target_pos[i])
        print(f"  moved appended slide -> position {target_pos[i] + 1}")

    # 3. Refresh footers
    total_after = len(prs.slides)
    for i, s in enumerate(prs.slides, start=1):
        for shape_name in ("Text 5", "Text 7"):
            num = find_shape(s, shape_name)
            if num is not None and "/" in num.text_frame.text:
                cur_size = None
                for p in num.text_frame.paragraphs:
                    for r in p.runs:
                        if r.font.size:
                            cur_size = r.font.size.pt; break
                    if cur_size: break
                cur_size = cur_size or 8
                set_text(num.text_frame, f"{i} / {total_after}",
                         size=cur_size, name="Cambria",
                         color=GREY if shape_name == "Text 5" else None)

    # 4. Update slide-4 roadmap if section ranges shifted
    slide4 = prs.slides[3]
    s7_idx = find_slide_by_title(prs, "Integrated disease cascade")
    s8_idx = find_slide_by_title(prs, "Reanalysis, deliverables")
    if s7_idx is not None and s8_idx is not None:
        s7_label = f"slides {s7_idx + 1} – {s8_idx}"
        s8_label = f"slides {s8_idx + 1} – {total_after}"
        for shape_name, label in (("Text 33", s7_label), ("Text 37", s8_label)):
            sh = find_shape(slide4, shape_name)
            if sh is not None:
                set_text(sh.text_frame, label, size=13, bold=False, name="Cambria", color=NAVY)
                print(f"  slide 4 {shape_name} -> {label!r}")

    prs.save(str(PPT))
    print(f"Saved {PPT} with {total_after} slides ({total_before} - 1 + 2 = {total_after}).")


if __name__ == "__main__":
    main()
