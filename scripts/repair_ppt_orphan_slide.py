"""
Repair GRBV_Funding_Committee_Presentation.pptx by removing the orphan
slide67.xml + its rels + the unreferenced rId in presentation.xml.rels +
the unused PartName in [Content_Types].xml.

Root cause: an earlier scripted slide-replace operation
(update_ppt_split_S80.py) removed slide 67 from the slide list but left
the underlying slide67.xml file and its relationship entry in the zip.
PowerPoint flagged the inconsistency on open with 'PowerPoint found a
problem with content...'.

This script:
  1. Verifies the orphan via three independent checks (no sldId refers to
     slide67.xml, only one rId points to it, the rId is not in the slide
     list).
  2. Reads every entry from the original zip.
  3. Skips ppt/slides/slide67.xml + ppt/slides/_rels/slide67.xml.rels.
  4. Strips the orphan <Relationship Id="rIdXX" Target="slides/slide67.xml" />
     line from ppt/_rels/presentation.xml.rels.
  5. Strips the <Override PartName="/ppt/slides/slide67.xml" /> line from
     [Content_Types].xml.
  6. Writes a clean zip to the same path.

A backup copy is taken before modification.
Run: python3 scripts/repair_ppt_orphan_slide.py
"""
from pathlib import Path
import re
import shutil
import zipfile
import sys
from datetime import datetime

PPT = Path("/Users/prempratapsingh/Data_Analysis/GRBV_Funding_Committee_Presentation.pptx")
ORPHAN_FILE = "ppt/slides/slide67.xml"
ORPHAN_RELS = "ppt/slides/_rels/slide67.xml.rels"


def diagnose(path):
    """Verify the orphan is genuinely orphaned. Returns the orphan rId."""
    with zipfile.ZipFile(path) as z:
        pres = z.read("ppt/presentation.xml").decode()
        rels = z.read("ppt/_rels/presentation.xml.rels").decode()

    sld_rids_used = set(re.findall(r'<p:sldId\s+id="\d+"\s+r:id="(rId\d+)"', pres))
    rid_target = dict(re.findall(
        r'Id="(rId\d+)"\s+Type="[^"]+slide"\s+Target="([^"]+)"', rels))

    pointing_to_67 = [rid for rid, tgt in rid_target.items()
                      if tgt == "slides/slide67.xml"]
    if len(pointing_to_67) != 1:
        sys.exit(f"abort: expected exactly 1 rId pointing to slide67.xml, "
                 f"got {len(pointing_to_67)}: {pointing_to_67}")
    orphan_rid = pointing_to_67[0]
    if orphan_rid in sld_rids_used:
        sys.exit(f"abort: rId {orphan_rid} IS referenced in the slide list "
                 f"-- not safe to remove.")
    print(f"  diagnostic: rId {orphan_rid} is the orphan (points to "
          f"slide67.xml but absent from sldIdLst)")
    return orphan_rid


def repair(path, orphan_rid):
    backup = path.with_suffix(
        f".pre_orphan_repair_{datetime.now():%Y%m%d_%H%M%S}.pptx")
    shutil.copy2(path, backup)
    print(f"  backup -> {backup}")

    tmp_out = path.with_suffix(".pptx.tmp")
    skipped_files = 0
    rels_modified = ct_modified = False

    with zipfile.ZipFile(path, "r") as zin, \
         zipfile.ZipFile(tmp_out, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename in (ORPHAN_FILE, ORPHAN_RELS):
                skipped_files += 1
                continue
            data = zin.read(item.filename)

            if item.filename == "ppt/_rels/presentation.xml.rels":
                txt = data.decode()
                # Remove the orphan <Relationship .../> element
                pattern = (rf'<Relationship\s+Id="{re.escape(orphan_rid)}"'
                           rf'\s+Type="[^"]+slide"'
                           rf'\s+Target="slides/slide67\.xml"\s*/>')
                new_txt, n = re.subn(pattern, "", txt)
                if n != 1:
                    sys.exit(f"abort: could not match orphan relationship "
                             f"line (matched {n} times)")
                data = new_txt.encode()
                rels_modified = True

            elif item.filename == "[Content_Types].xml":
                txt = data.decode()
                pattern = (r'<Override\s+PartName="/ppt/slides/slide67\.xml"'
                           r'\s+ContentType="[^"]+"\s*/>')
                new_txt, n = re.subn(pattern, "", txt)
                if n != 1:
                    sys.exit(f"abort: could not match Content_Types override "
                             f"for slide67.xml (matched {n} times)")
                data = new_txt.encode()
                ct_modified = True

            zout.writestr(item, data)

    if not (rels_modified and ct_modified and skipped_files == 2):
        sys.exit(f"abort: incomplete repair (rels_modified={rels_modified}, "
                 f"ct_modified={ct_modified}, files_removed={skipped_files})")

    tmp_out.replace(path)
    print(f"  repair complete: removed {skipped_files} files + 1 rel + 1 CT entry")


def verify(path):
    with zipfile.ZipFile(path) as z:
        files = set(z.namelist())
        assert ORPHAN_FILE not in files, "slide67.xml still present"
        assert ORPHAN_RELS not in files, "slide67.xml.rels still present"
        rels = z.read("ppt/_rels/presentation.xml.rels").decode()
        assert "slide67.xml" not in rels, \
            "presentation.xml.rels still mentions slide67.xml"
        ct = z.read("[Content_Types].xml").decode()
        assert "slide67.xml" not in ct, \
            "[Content_Types].xml still mentions slide67.xml"
        # Slide-rId count = sldId count
        pres = z.read("ppt/presentation.xml").decode()
        n_sld = len(re.findall(r'<p:sldId\s+id="\d+"\s+r:id="rId\d+"', pres))
        n_rid = len(re.findall(
            r'Id="rId\d+"\s+Type="[^"]+slide"\s+Target="', rels))
        assert n_sld == n_rid, f"counts still mismatch ({n_sld} vs {n_rid})"
        print(f"  verify: deck has {n_sld} slides, {n_rid} slide rIds (match)")


def main():
    print(f"Repairing {PPT}")
    orphan_rid = diagnose(PPT)
    repair(PPT, orphan_rid)
    verify(PPT)
    print("Done.")


if __name__ == "__main__":
    main()
