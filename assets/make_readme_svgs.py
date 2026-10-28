"""Write the README header graphics (banner, section bars, step bars) to assets/.

Run from the repo root: python assets/make_readme_svgs.py
"""
from pathlib import Path

OUT = Path(__file__).parent
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"
TEAL, BURG, GOLD, MUTED, TEXT, SUB = "#3E9C93", "#B04A5A", "#C9A14A", "#8C7F86", "#F3EEF0", "#CBBDC3"
BG = ('<linearGradient id="bg" x1="0" y1="0" x2="1" y2="1">'
      '<stop offset="0%" stop-color="#1E1A1D"/><stop offset="55%" stop-color="#2A2127"/>'
      '<stop offset="100%" stop-color="#5E1E31"/></linearGradient>')
STAGES = ["Pre-veraison", "Veraison", "Post-veraison", "Harvest"]


def svg(w, h, label, body, extra_defs=""):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            f'role="img" aria-label="{label}"><defs>{BG}{extra_defs}</defs>'
            f'<rect width="{w}" height="{h}" rx="12" fill="url(#bg)"/>{body}</svg>\n')


def text(x, y, s, size, fill, weight="400", anchor="start", spacing=None):
    ls = f' letter-spacing="{spacing}"' if spacing else ""
    return (f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
            f'fill="{fill}" text-anchor="{anchor}"{ls}>{s}</text>')


def berries(x0, y, lag, color, label):
    """Four stage berries; `lag` shifts the row left to show developmental displacement."""
    out = [text(x0 - 90, y + 5, label, 11, MUTED, anchor="end")]
    for i in range(4):
        cx = x0 + i * 62 - lag * 62
        r = 9 + i * 3
        out.append(f'<circle cx="{cx:.0f}" cy="{y}" r="{r}" fill="{color}" opacity="{0.45 + i * 0.17:.2f}"/>')
    return "".join(out)


def banner():
    glow = ('<radialGradient id="glow" cx="0.85" cy="0.15" r="0.7"><stop offset="0%" stop-color="#3E9C93" '
            'stop-opacity=".22"/><stop offset="100%" stop-color="#3E9C93" stop-opacity="0"/></radialGradient>')
    axis = "".join(text(900 + i * 62, 300, s, 10, MUTED, anchor="middle") for i, s in enumerate(STAGES))
    body = (
        '<rect width="1200" height="360" rx="12" fill="url(#glow)"/>'
        f'<rect x="56" y="70" width="6" height="96" rx="3" fill="{TEAL}"/>'
        + text(80, 92, "GRBV TRANSCRIPTOMICS", 15, TEAL, "700", spacing="2.4")
        + text(80, 140, "Ripening, out of step", 44, TEXT, "700")
        + text(80, 176, "Grapevine red blotch disease and the timing of berry ripening programs", 19, SUB)
        + text(80, 250, "327 libraries · 44 Merlot vines · 4 stages · 2 seasons", 16, SUB)
        + text(80, 280, "Repeated-measures models · healthy-only developmental clock · pathway atlas", 14, MUTED)
        + berries(900, 140, 0, TEAL, "Healthy")
        + berries(900, 210, 0.8, BURG, "Infected")
        + f'<path d="M790 250 H1110" stroke="#4A3B42" stroke-width="2"/>' + axis
        + text(950, 330, "infected berries sit earlier on the ripening axis", 11, MUTED, anchor="middle")
    )
    return svg(1200, 360, "GRBV transcriptomics: grapevine red blotch disease and the timing of berry ripening "
                          "programs, 327 libraries from 44 Merlot vines", body, glow)


def section(title, subtitle):
    dots = "".join(
        f'<circle cx="{1010 + i * 34}" cy="60" r="{6 if i else 9}" fill="{c if not i else c}" '
        f'opacity="{1 if not i else 0.55}"/>' for i, c in enumerate([TEAL, BURG, GOLD, TEAL, BURG]))
    body = (f'<rect x="28" y="30" width="6" height="60" rx="3" fill="{TEAL}"/>'
            + text(52, 52, title.upper(), 15, TEAL, "700", spacing="2.4")
            + text(52, 84, subtitle, 22, TEXT)
            + f'<path d="M1010 60 H1146" stroke="#4A3B42" stroke-width="2"/>' + dots)
    return svg(1200, 120, f"{title}: {subtitle}", body)


def step(n, total, title, subtitle):
    colors = [TEAL, BURG, GOLD, "#6C7DA6"]
    dots = "".join(
        f'<circle cx="{1044 + i * 34}" cy="42" r="{9 if i + 1 == n else 6}" '
        f'fill="{colors[i] if i + 1 == n else "none"}" stroke="{colors[i]}" stroke-width="2" '
        f'opacity="{1 if i + 1 == n else 0.55}"/>' for i in range(total))
    body = (f'<circle cx="52" cy="42" r="24" fill="#150F12" stroke="{TEAL}" stroke-width="3"/>'
            + text(52, 50, n, 24, TEAL, "700", anchor="middle")
            + text(96, 36, title, 22, TEXT, "700")
            + text(96, 61, subtitle, 16, SUB)
            + f'<path d="M1044 42 H{1044 + (total - 1) * 34}" stroke="#4A3B42" stroke-width="2"/>' + dots)
    return svg(1200, 84, f"Step {n} of {total}: {title}", body)


SECTIONS = {
    "findings": ("What the study found", "Common delay plus program-specific timing and level changes"),
    "workflow": ("How the analysis runs", "From TAG-seq counts to claim-gated results"),
    "quickstart": ("Quick start", "Clone, set up the environments, add the data, run"),
    "results": ("Find and interpret results", "Where each output lands and how to read it"),
    "history": ("How the analysis got here", "Seven analysis eras and the claims they withdrew"),
    "limits": ("Known limits", "What this study can and cannot say"),
    "license": ("License and citation", "MIT source, citable analysis"),
}
STEPS = [
    ("Get the code", "Clone the repository and choose a working folder"),
    ("Set up the environments", "Python for local scripts, Conda and R on FARM"),
    ("Add the input data", "Count matrix, sample metadata, titer and weather tables"),
    ("Run the analysis", "Snakemake on FARM, then the R figure scripts"),
]

if __name__ == "__main__":
    (OUT / "readme-banner.svg").write_text(banner())
    for key, (title, sub) in SECTIONS.items():
        (OUT / f"section-{key}.svg").write_text(section(title, sub))
    for i, (title, sub) in enumerate(STEPS, 1):
        (OUT / f"step-{i}.svg").write_text(step(i, len(STEPS), title, sub))
    print(f"wrote {1 + len(SECTIONS) + len(STEPS)} SVGs to {OUT}")
