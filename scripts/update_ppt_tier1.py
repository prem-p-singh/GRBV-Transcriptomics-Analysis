"""
Add two Tier 1 meta-analysis slides to GRBV_Funding_Committee_Presentation.pptx
after the SPS1-reframed slide (biochem batch), and extend the
"What the reanalysis corrected" slide with a Tier 1 bullet.

Reuses the slide-cloning and template pattern from update_ppt_biochem.py.

Run: python3 scripts/update_ppt_tier1.py
"""
from copy import deepcopy
from pathlib import Path
from io import BytesIO

from pptx import Presentation
from pptx.util import Pt, Emu
from pptx.dml.color import RGBColor
from pptx.oxml.ns import qn

ROOT = Path("/Users/prempratapsingh/Data_Analysis")
PPT  = ROOT / "GRBV_Funding_Committee_Presentation.pptx"
FIGS = ROOT / "figures" / "supp"

TEMPLATE_INDEX = 27           # slide 28 "Direct vs indirect temperature effects" -- 2-col layout
SECTION_LABEL  = "PLUS: Tier 1 meta"
SECTION_ROOT   = "6 / Results"

TIER1_SLIDES = [
    {
        "title": "PLUS (Tier 1): GRBV signature vs 6 published Vitis stress datasets",
        "image": "S64_tier1_overlap_summary.png",
        "body": [
            "List-level meta-analysis on 6 verified papers: Blanco-Ulate 2017 GRBV/Zinfandel, Song 2022 GLRaV-3/Cabernet Franc, Fung 2008 powdery mildew, Polesani 2010 downy mildew, Agudelo-Romero 2015 Botrytis, Liu 2012 heat.",
            "ID-space crosswalk (VIT 12X <-> NCBI LOC, 24,817 rows) gives 68-91 % mapping for the 3 RNA-Seq / VIT papers. The 3 probe/TC-based papers map at 0-0.5 % and are carried forward as Tier 2 (raw-count re-analysis) candidates.",
            "Signature size (panel a): Botrytis dominates at 1670 DEGs; GRBV this-study 55 DEGs; GLRaV-3 131; Blanco-Ulate GRBV/Zinfandel 43. The three probe-based papers contribute 0-2 mapped DEGs.",
            "GRBV signature partition (panel b): 35 of 55 DEGs are GRBV-only (64 %); 13 overlap with Botrytis (berry-stress core); 3 overlap exclusively with Blanco-Ulate GRBV/Zinfandel; 7 are shared with two or more comparators; 1 overlaps GLRaV-3.",
            "Fisher enrichment (panel c): GRBV (Merlot) vs Blanco-Ulate GRBV (Zinfandel) shows log2(OR) = 5.5 with BH q = 2.1e-5. Botrytis comparator reaches log2(OR) = 1.5 at q = 1.3e-2. GLRaV-3 and the three probe-based papers show no significant enrichment.",
        ],
    },
    {
        "title": "Tier 1 log2FC concordance: GRBV replicates across cultivars, diverges from GLRaV-3",
        "image": "S65_tier1_fc_concordance.png",
        "body": [
            "Spearman rho of log2 fold-change on shared tested genes across all pairs of Tier 1 datasets. Square heatmap; NA cells indicate zero ID-space overlap.",
            "Cross-cultivar GRBV replication: GRBV (Merlot, this study) vs Blanco-Ulate 2017 GRBV (Zinfandel, Oakville harvest) reaches rho = 0.47 on 11,886 shared tested genes. First cross-cultivar validation of the inverted-trajectory / sHSP / MBF1c findings.",
            "GRBV vs GLRaV-3: rho ~ 0 with Song 2022 at every developmental stage. Two viruses produce distinct host-transcriptome signatures; the Merlot GRBV response is not a generic Vitis viral response.",
            "Blanco-Ulate GRBV vs Song GLRaV-3: rho = -0.21. Mildly anti-correlated, strengthens the distinctness of the two viruses.",
            "AgudeloRomero Botrytis is only loosely coupled (rho = 0.00-0.03 with the two GRBV signatures), confirming that the Botrytis overlap in the intersection counts is set-level only, not effect-size-level.",
            "Tier 2 path forward: raw-count re-analysis from GEO/SRA for all 6 papers under a unified DESeq2 pipeline would unlock the 3 zero-coverage papers and give effect-size-level meta-analysis across all 6. 1-2 weeks of work on FARM.",
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
    r.text = ("Tier 1 cross-study meta-analysis added: GRBV signature compared "
              "against 6 published Vitis stress papers (verified citations). "
              "Cross-cultivar replication confirmed (GRBV Merlot vs GRBV Zinfandel, "
              "Spearman rho = 0.47, Fisher q = 2.1e-5). GRBV vs GLRaV-3 rho = 0, "
              "supporting stress-specific host signature. Figs S64 and S65.")
    r.font.size = Pt(10); r.font.name = "Arial"
    print(f"  added Tier 1 correction bullet to slide {idx + 1}")


def main():
    prs = Presentation(str(PPT))
    total_before = len(prs.slides)
    print(f"Opened presentation with {total_before} slides.")

    # 1. Extend "What the reanalysis corrected" BEFORE we shuffle indices
    append_reanalysis_bullet(prs)

    # 2. Locate the SPS1-reframed slide (end of biochem batch); insert after
    sps1_idx = find_slide_by_title(prs, "SPS1 reframed")
    if sps1_idx is None:
        print("  warn: could not find 'SPS1 reframed' slide; appending at end")
        sps1_idx = len(prs.slides) - 1
    print(f"  inserting Tier 1 slides after slide {sps1_idx + 1}")

    # 3. Append the two new slides then reorder them to positions sps1_idx+1, +2
    added_idx = []
    for entry in TIER1_SLIDES:
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

    # Reorder the two freshly-appended slides to immediately follow sps1_idx
    target_pos = [sps1_idx + 1, sps1_idx + 2]
    for i, new_abs_idx in enumerate(added_idx):
        cur = len(prs.slides) - (len(added_idx) - i)
        reorder_slide(prs, cur, target_pos[i])

    # 4. Refresh slide-number footers
    total_after = len(prs.slides)
    for i, s in enumerate(prs.slides, start=1):
        num = find_shape(s, "Text 7")
        if num is not None and "/" in num.text_frame.text:
            set_text(num.text_frame, f"{i} / {total_after}", size=8)

    prs.save(str(PPT))
    print(f"Saved {PPT} with {total_after} slides (+{total_after - total_before}).")


if __name__ == "__main__":
    main()
