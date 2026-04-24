"""
Consistency + style fixes applied to GRBV_Funding_Committee_Presentation.pptx
on 2026-04-24.

Idempotent replay script (skips edits already applied). Fixes three bug classes
surfaced by a systematic review of the 71-slide deck:

  1. Pagination footers: 12 slides carried stale `N / 55` footers from the
     earlier 55-slide deck. Corrected to `N / 71`.
  2. Slide 4 roadmap: all 8 section ranges pointed to old slide numbers.
     Rewritten to match current section-header locations (5-10 / 11-15 /
     16-25 / 26-34 / 35-47 / 48-54 / 65-66 / 67-71).
  3. Sidebar section labels: 13 slides had wrong or ambiguous `N / Topic`
     tags. Fixed to distinguish post-viral PLUS content (Biochem / Tier 1 /
     Cohort / Dimensions) from the Section 6 viral slides.

Narrative + cross-reference fixes:

  4. Slide 7: "mechanism-cascade of slide 58" -> "slide 66" (Main10 location).
  5. Slide 10: all 7 Q-to-slide pointers updated to current slide numbers.
  6. Slide 43: "NES heatmap (slide 39)" -> "slide 42".
  7. Slide 14: "Twenty-four stratified contrasts" -> "Twenty-two" (arithmetic
     fix: 3+5+4+8+2=22).
  8. Slide 64: figure callout clarified to name panel (S66 panel D; S67 panel C).
  9. Slide 69: deliverable count updated to 12 Main + 64 supplementary figures
     (was 10 + 55, stale after the figure-upgrade pass).
 10. Slide 29: bridge sentence appended acknowledging slide 61's
     variance-decomposition view, resolving a contradiction between the
     "infection-dominant" framing (slide 29) and the passive-failure reframe
     (slides 61-64).

Typography + color unification:

 11. Font family: 89 body-text runs were Arial 10pt on later-added PLUS slides;
     normalised to Cambria 10pt + color #0F172A to match the canonical body
     style. Size kept at 10pt because 16 of 18 boxes would overflow at 13pt.
 12. Remaining 71 theme-inherited (font=None) runs set to explicit Cambria
     to prevent rendering under the theme's sans-serif major font.
 13. Slide 2 body sublabels ("covariate-adjusted" / "core infection DEGs")
     12pt -> 13pt to match canonical body size.
 14. 18 content-slide titles (slides 3, 7, 12, 19, 27, 30, 43, 47, 55-64) were
     rendering BLACK on the navy header band because color=None defaulted to
     black. Set explicitly to #FFFFFF white.
 15. Slide 4 roadmap first range ("slides 5 - 10") had gray #64748B while the
     other 7 were navy #0F172A. Normalised to navy.

Structural items NOT changed (need design input):
  - Section-header colors: 8 sections use 5 colors with 2 reuses
    (S1/S6 orange, S3/S4/S5/S7 blue). Intentional grouping or drift unknown.
  - Two footer-shape variants coexist (Text 5 at 10pt Cambria #64748B vs
    Text 7 at 8pt theme-default).
  - Body-size split: slides added during the fix/figure-upgrade pass use
    Cambria 10pt; canonical body is 13pt. Unification requires resizing
    text boxes, not just the font.

Run: python3 scripts/update_ppt_consistency_fixes.py
"""
from copy import copy
from pathlib import Path

from pptx import Presentation
from pptx.util import Pt
from pptx.dml.color import RGBColor

ROOT = Path("/Users/prempratapsingh/Data_Analysis")
PPT  = ROOT / "GRBV_Funding_Committee_Presentation.pptx"

NAVY  = RGBColor(0x0F, 0x17, 0x2A)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

prs = Presentation(PPT)


# ---------- helpers ------------------------------------------------------

def set_shape_text(slide, shape_name, new_text, *, slide_num=None):
    """Replace shape's first run text; keep its font. Clear extra runs/paragraphs."""
    for shape in slide.shapes:
        if shape.name == shape_name and shape.has_text_frame:
            tf = shape.text_frame
            p0 = tf.paragraphs[0]
            if not p0.runs:
                raise RuntimeError(f"slide {slide_num} {shape_name}: no runs")
            p0.runs[0].text = new_text
            for r in list(p0.runs[1:]):
                r._r.getparent().remove(r._r)
            for p in list(tf.paragraphs[1:]):
                p._p.getparent().remove(p._p)
            return True
    return False


def replace_line_in_shape(slide, shape_name, old_substr, new_substr, *, slide_num=None):
    """Find run containing old_substr; replace in place. Idempotent if already changed."""
    for shape in slide.shapes:
        if shape.name == shape_name and shape.has_text_frame:
            for para in shape.text_frame.paragraphs:
                for run in para.runs:
                    if new_substr in run.text:
                        return False  # already applied
                    if old_substr in run.text:
                        run.text = run.text.replace(old_substr, new_substr)
                        return True
            for para in shape.text_frame.paragraphs:
                full = "".join(r.text for r in para.runs)
                if new_substr in full:
                    return False
                if old_substr in full:
                    para.runs[0].text = full.replace(old_substr, new_substr)
                    for r in list(para.runs[1:]):
                        r._r.getparent().remove(r._r)
                    return True
            return False
    return False


def append_paragraph(slide, shape_name, new_text, marker=None, *, slide_num=None):
    """Clone last paragraph's formatting, append new paragraph. Skip if marker already present."""
    for shape in slide.shapes:
        if shape.name == shape_name and shape.has_text_frame:
            tf = shape.text_frame
            full = tf.text
            if marker and marker in full:
                return False
            last_p = tf.paragraphs[-1]
            new_p = copy(last_p._p)
            ns = {"a": "http://schemas.openxmlformats.org/drawingml/2006/main"}
            runs = new_p.findall("a:r", ns)
            for r in runs[1:]:
                new_p.remove(r)
            for tag in ("a:br", "a:fld"):
                for el in new_p.findall(tag, ns):
                    new_p.remove(el)
            first_run = new_p.find("a:r", ns)
            t = first_run.find("a:t", ns)
            t.text = new_text
            last_p._p.addnext(new_p)
            return True
    return False


# ---------- 1. Pagination footers ----------------------------------------
footer_fixes = {2: "2 / 71", 4: "4 / 71", 6: "6 / 71", 8: "8 / 71",
                9: "9 / 71", 10: "10 / 71", 13: "13 / 71", 14: "14 / 71",
                54: "54 / 71", 68: "68 / 71", 69: "69 / 71", 70: "70 / 71"}
n = 0
for sn, new_footer in footer_fixes.items():
    if set_shape_text(prs.slides[sn-1], "Text 5", new_footer, slide_num=sn):
        n += 1
print(f"[1] Pagination footers rewritten: {n}/{len(footer_fixes)}")


# ---------- 2. Slide 4 roadmap ranges ------------------------------------
roadmap = [("Text 9",  "slides 5 – 10"),   ("Text 13", "slides 11 – 15"),
           ("Text 17", "slides 16 – 25"),  ("Text 21", "slides 26 – 34"),
           ("Text 25", "slides 35 – 47"),  ("Text 29", "slides 48 – 54"),
           ("Text 33", "slides 65 – 66"),  ("Text 37", "slides 67 – 71")]
n = 0
for shape_name, new in roadmap:
    if set_shape_text(prs.slides[3], shape_name, new, slide_num=4):
        n += 1
print(f"[2] Roadmap ranges on slide 4: {n}/{len(roadmap)}")


# ---------- 3. Sidebar section labels ------------------------------------
sidebars = [(19, "4 / Results"),    (43, "5 / Results"),     (47, "5 / Results"),
            (55, "6 / Biochem"),    (56, "6 / Biochem"),     (57, "6 / Biochem"),
            (58, "6 / Tier 1"),     (59, "6 / Tier 1"),      (60, "6 / Cohort"),
            (61, "6 / Dimensions"), (62, "6 / Dimensions"),
            (63, "6 / Dimensions"), (64, "6 / Dimensions")]
n = 0
for sn, new in sidebars:
    if set_shape_text(prs.slides[sn-1], "Text 6", new, slide_num=sn):
        n += 1
print(f"[3] Sidebar labels: {n}/{len(sidebars)}")


# ---------- 4-10. Narrative / cross-reference / arithmetic ---------------
changes = [
    (7,  "Text 9", "matching the mechanism-cascade of slide 58.",
                   "matching the mechanism-cascade of slide 66."),
    (43, "Text 9", "existing NES heatmap (slide 39)",
                   "existing NES heatmap (slide 42)"),
    (14, "Text 6", "Twenty-four stratified contrasts computed: 3 infection-duration bins",
                   "Twenty-two stratified contrasts computed: 3 infection-duration bins"),
    (64, "Text 9", "Harvest juice Brix from the biochemistry batch, cross-referenced",
                   "Harvest juice Brix (S66 panel D; S67 panel C) cross-referenced"),
    (69, "Text 6", "Completed: 10 Main composite figures + 55 supplementary figures",
                   "Completed: 12 Main composite figures + 64 supplementary figures"),
]
n = 0
for sn, sh, old, new in changes:
    if replace_line_in_shape(prs.slides[sn-1], sh, old, new, slide_num=sn):
        n += 1
print(f"[4-10] Narrative/cross-ref/arithmetic edits: {n}/{len(changes)}")

# Slide 10 Q-routing
q_routing = [("Text 11", "slides 20, 21, 44"),
             ("Text 17", "slides 20, 22, 36, 45"),
             ("Text 23", "slides 23, 24, 25, 41, 42"),
             ("Text 29", "slides 23, 30, 31, 32"),
             ("Text 35", "slide 23"),
             ("Text 41", "slides 20, 46, 66"),
             ("Text 47", "slides 27–30, 33, 34, 39")]
n = 0
for shape_name, new in q_routing:
    if set_shape_text(prs.slides[9], shape_name, new, slide_num=10):
        n += 1
print(f"       Slide 10 Q-pointers: {n}/{len(q_routing)}")


# ---------- 10. Slide 29 bridge sentence ---------------------------------
bridge = ("Bridge: slide 61 shows sHSP viral-variance share (0.11) matches the 0.12 "
          "all-genes baseline, so 'infection-dominant' here reflects the mean-FC deficit, "
          "not active virus-driven transcription.")
added = append_paragraph(prs.slides[28], "Text 9", bridge, marker="Bridge: slide 61", slide_num=29)
print(f"[10] Slide 29 bridge sentence: {'added' if added else 'already present'}")


# ---------- 11-12. Font normalisation (Arial + None -> Cambria) ----------
ARIAL_BODY_SLIDES = [3, 7, 12, 19, 30, 43, 47, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 68]
arial_fixed = 0
for sn in ARIAL_BODY_SLIDES:
    for sh in prs.slides[sn-1].shapes:
        if sh.name == "Text 9" and sh.has_text_frame:
            for p in sh.text_frame.paragraphs:
                for r in p.runs:
                    if r.font.name == "Arial":
                        r.font.name = "Cambria"
                        r.font.color.rgb = NAVY
                        arial_fixed += 1

# Any remaining Arial anywhere
for slide in prs.slides:
    for sh in slide.shapes:
        if sh.has_text_frame:
            for p in sh.text_frame.paragraphs:
                for r in p.runs:
                    if r.font.name == "Arial":
                        r.font.name = "Cambria"
                        r.font.color.rgb = NAVY
                        arial_fixed += 1

# None-font runs -> Cambria explicit
none_fixed = 0
for slide in prs.slides:
    for sh in slide.shapes:
        if sh.has_text_frame:
            for p in sh.text_frame.paragraphs:
                for r in p.runs:
                    if r.font.name is None and r.text.strip():
                        r.font.name = "Cambria"
                        none_fixed += 1
print(f"[11-12] Arial -> Cambria: {arial_fixed}; None -> Cambria: {none_fixed}")


# ---------- 13. Slide 2 body sublabels 12pt -> 13pt ----------------------
n = 0
for sh in prs.slides[1].shapes:
    if sh.name == "Text 9" and sh.has_text_frame:
        for p in sh.text_frame.paragraphs:
            for r in p.runs:
                if r.font.size and r.font.size.pt == 12.0:
                    r.font.size = Pt(13)
                    n += 1
print(f"[13] Slide 2 body 12pt -> 13pt: {n} runs")


# ---------- 14. Black-on-navy titles -> white ----------------------------
title_fixed = 0
for slide in prs.slides:
    for sh in slide.shapes:
        if sh.name == "Text 1" and sh.has_text_frame:
            for p in sh.text_frame.paragraphs:
                for r in p.runs:
                    if not (r.font.size and r.font.size.pt == 16.0): continue
                    has_col = False
                    try:
                        if r.font.color.type == 1 and r.font.color.rgb:
                            has_col = True
                    except: pass
                    if not has_col and r.text.strip():
                        r.font.color.rgb = WHITE
                        title_fixed += 1
print(f"[14] Black-on-navy titles -> white: {title_fixed}")


# ---------- 15. Slide 4 roadmap gray range -> navy -----------------------
n = 0
for sh in prs.slides[3].shapes:
    if sh.name == "Text 9" and sh.has_text_frame:
        for p in sh.text_frame.paragraphs:
            for r in p.runs:
                try:
                    if (r.font.color.type == 1 and r.font.color.rgb
                            and f"#{r.font.color.rgb}" == "#64748B"):
                        r.font.color.rgb = NAVY
                        n += 1
                except: pass
print(f"[15] Slide 4 roadmap gray -> navy: {n}")


prs.save(PPT)
print(f"\nSaved -> {PPT}")
