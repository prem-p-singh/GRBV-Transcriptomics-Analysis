"""Publication style for every figure in the GRBV manuscript.

Single source of truth for fonts (Arial, strict), colour palette, figure
dimensions, and panel-label helpers. Every figure script imports this
module first thing and calls `apply()` before any plt call.
"""
from __future__ import annotations

import warnings
from typing import Iterable

import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

# ---------------------------------------------------------------------
# Sizes (millimetres, converted to inches inside helpers)
# ---------------------------------------------------------------------
MM_PER_INCH = 25.4

# Standard publication widths
W_SINGLE_MM  = 89.0    # single-column
W_ONEHALF_MM = 120.0   # 1.5-column
W_DOUBLE_MM  = 183.0   # double-column / full page width

# Standard heights
H_SHORT_MM  = 60.0
H_MEDIUM_MM = 100.0
H_TALL_MM   = 160.0
H_PAGE_MM   = 240.0

def fig_size(width_mm: float, height_mm: float) -> tuple[float, float]:
    return width_mm / MM_PER_INCH, height_mm / MM_PER_INCH

# Convenience size constants
SIZE_SINGLE   = fig_size(W_SINGLE_MM,  H_MEDIUM_MM)
SIZE_DOUBLE_S = fig_size(W_DOUBLE_MM,  H_MEDIUM_MM)
SIZE_DOUBLE_M = fig_size(W_DOUBLE_MM,  H_TALL_MM)
SIZE_PAGE     = fig_size(W_DOUBLE_MM,  H_PAGE_MM)

# ---------------------------------------------------------------------
# Colour palette (consistent across every figure)
# ---------------------------------------------------------------------
COLORS = {
    "RB_pos":       "#c0392b",
    "RB_neg":       "#27ae60",
    "season_2021":  "#3498db",
    "season_2022":  "#e67e22",
    "stage_pre":    "#74b9ff",
    "stage_ver":    "#0984e3",
    "stage_post":   "#6c5ce7",
    "stage_harv":   "#2d3436",
    "up":           "#c0392b",
    "down":         "#2980b9",
    "neutral":      "#95a5a6",
    "highlight":    "#8e44ad",
    "warm":         "#d35400",
    "cool":         "#2980b9",
}

STAGE_ORDER  = ["Pre_veraison", "Veraison", "Post_veraison", "Harvest"]
STAGE_LABEL  = {"Pre_veraison": "Pre-véraison",
                "Veraison":      "Véraison",
                "Post_veraison": "Post-véraison",
                "Harvest":       "Harvest"}
STAGE_COLOR  = dict(zip(STAGE_ORDER,
                        [COLORS["stage_pre"], COLORS["stage_ver"],
                         COLORS["stage_post"], COLORS["stage_harv"]]))
SEASON_COLOR = {2021: COLORS["season_2021"], 2022: COLORS["season_2022"]}
INFECT_COLOR = {"pos": COLORS["RB_pos"], "neg": COLORS["RB_neg"],
                "RB(+)": COLORS["RB_pos"], "RB(-)": COLORS["RB_neg"]}

# ---------------------------------------------------------------------
# Apply rcParams (call once at the top of every figure script)
# ---------------------------------------------------------------------
def apply() -> None:
    """Set matplotlib rcParams for publication output. Idempotent."""
    arial_present = any("arial" in f.name.lower() for f in fm.fontManager.ttflist)
    if not arial_present:
        warnings.warn("Arial not found; falling back to DejaVu Sans.",
                      RuntimeWarning, stacklevel=2)

    matplotlib.rcParams.update({
        "font.family":      "sans-serif",
        "font.sans-serif":  ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size":        10.0,
        "axes.titlesize":   12.0,
        "axes.titleweight": "bold",
        "axes.labelsize":   10.0,
        "axes.labelweight": "normal",
        "xtick.labelsize":   9.0,
        "ytick.labelsize":   9.0,
        "legend.fontsize":   8.5,
        "figure.titlesize": 13.0,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.linewidth":    0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size":  3.0,
        "ytick.major.size":  3.0,
        "axes.grid":         False,
        "grid.alpha":        0.25,
        "grid.linestyle":    "-",
        "grid.linewidth":    0.5,
        "savefig.dpi":       300,
        "figure.dpi":        150,
        "figure.facecolor":  "white",
        "axes.facecolor":    "white",
        "savefig.facecolor": "white",
        "savefig.bbox":      "tight",
        "savefig.pad_inches": 0.05,
        "pdf.fonttype":      42,
        "ps.fonttype":       42,
        "svg.fonttype":      "none",
    })

def panel_label(ax: Axes, label: str,
                x: float = -0.12, y: float = 1.06,
                fontsize: float = 13.0) -> None:
    """Add an A/B/C panel label in the upper-left of an axes."""
    ax.text(x, y, label,
            transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold",
            va="bottom", ha="left")

def cleanup_axes(ax: Axes,
                 keep_top: bool = False, keep_right: bool = False) -> None:
    """Hide top/right spines (default) and tighten ticks."""
    ax.spines["top"].set_visible(keep_top)
    ax.spines["right"].set_visible(keep_right)
    ax.tick_params(axis="both", which="both", direction="out", length=3, width=0.8)

def save(fig: Figure, path_no_ext: str,
         formats: Iterable[str] = ("png", "pdf")) -> list[str]:
    """Save a figure under one stem in multiple formats. Returns paths written."""
    out = []
    for ext in formats:
        p = f"{path_no_ext}.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        out.append(p)
    return out

if __name__ == "__main__":
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    apply()
    arial = [f.name for f in fm.fontManager.ttflist if "arial" in f.name.lower()]
    print(f"[pub_style] Arial fonts available: {len(arial)}")
    print(f"[pub_style] resolved family: {plt.rcParams['font.sans-serif'][0]}")
    fig, ax = plt.subplots(figsize=SIZE_SINGLE)
    x = np.linspace(0, 10, 100); ax.plot(x, np.sin(x), color=COLORS["RB_pos"])
    ax.set_title("pub_style smoke test"); ax.set_xlabel("x"); ax.set_ylabel("sin(x)")
    cleanup_axes(ax); panel_label(ax, "A")
    out_dir = os.path.dirname(os.path.abspath(__file__))
    save(fig, os.path.join(out_dir, "_pub_style_smoketest"))
    print(f"[pub_style] wrote {out_dir}/_pub_style_smoketest.png/.pdf")
