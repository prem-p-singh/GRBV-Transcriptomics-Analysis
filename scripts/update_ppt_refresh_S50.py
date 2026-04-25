"""
Refresh slide 32 to use the regenerated S50 PNG (Panel A readability pass +
Panel B 1,441 label fix).  Replaces the picture placeholder in-place; keeps
the body bullets unchanged.

Run: python3 scripts/update_ppt_refresh_S50.py
"""
from pathlib import Path
from pptx import Presentation
from pptx.util import Pt
from pptx.dml.color import RGBColor

ROOT = Path("/Users/prempratapsingh/Data_Analysis")
PPT  = ROOT / "GRBV_Funding_Committee_Presentation.pptx"
NEW_IMG = ROOT / "figures" / "supp" / "S50_direct_vs_indirect_effects.png"

ANCHOR = "Direct vs indirect temperature effects"


def find_shape(slide, name):
    for sh in slide.shapes:
        if sh.name == name:
            return sh
    return None


def main():
    if not NEW_IMG.exists():
        raise SystemExit(f"missing {NEW_IMG}")
    prs = Presentation(str(PPT))
    target = None
    for i, s in enumerate(prs.slides):
        t = find_shape(s, "Text 1")
        if t is not None and ANCHOR.lower() in t.text_frame.text.lower():
            target = (i, s); break
    if target is None:
        raise SystemExit(f"could not find anchor slide '{ANCHOR}'")
    idx, slide = target
    print(f"refreshing slide {idx + 1} ({ANCHOR})")

    pic = None
    for sh in slide.shapes:
        if sh.shape_type == 13:
            pic = sh; break
    if pic is None:
        raise SystemExit("no picture shape on target slide")

    left, top, w, h = pic.left, pic.top, pic.width, pic.height
    pic._element.getparent().remove(pic._element)
    slide.shapes.add_picture(str(NEW_IMG), left, top, width=w, height=h)
    print(f"  replaced picture (left={left}, top={top}, w={w}, h={h})")

    # Append a single-line note that the image was refreshed
    body = None
    for sh in slide.shapes:
        if sh.has_text_frame and sh.name.startswith("Text") and sh.name not in (
                "Text 1", "Text 3", "Text 4", "Text 5", "Text 6", "Text 7"):
            if body is None or sh.width * sh.height > body.width * body.height:
                body = sh
    if body is not None:
        tf = body.text_frame
        already = "readability pass" in tf.text.lower()
        if not already:
            p = tf.add_paragraph(); p.space_after = Pt(4)
            r = p.add_run()
            r.text = ("Refreshed 2026-04-24: legend moved out of plot, "
                      "z-order reordered (titer-mediated points layered on "
                      "top), white-filled labels, quadrant annotations "
                      "('Concordant'/'Discordant'), patchwork title bug "
                      "fixed, and the 1,441 'Not significant' label in panel "
                      "B now renders fully (was clipped to '1,44').")
            r.font.size = Pt(10); r.font.name = "Cambria"
            r.font.color.rgb = RGBColor(0x0F, 0x17, 0x2A)
            print("  appended readability-pass bullet")
        else:
            print("  bullet already present; skipping")

    prs.save(str(PPT))
    print("Saved.")


if __name__ == "__main__":
    main()
