# theme_grbv.R — shared ggplot2 theme + palette for the GRBV manuscript.
# Sourced at the top of every figure script.
#
# 2026-04-22 revision:
#   * Re-worked palette for strict one-colour-per-semantic-role usage.
#   * Red/blue reserved for up/down regulation (DEG direction).
#   * RB(+)/RB(-) moved to Okabe-Ito vermilion / sky-blue pair.
#   * Season moved to bluish-green (2021, cool) / magenta (2022, hot).
#   * Stage is a viridis 4-step ordinal gradient (not blue, not crimson).
#   * Mediation palette is amber / blue / grey (no red-green).
#   * New named palettes: pal_direction, pal_infection, pal_season,
#     pal_stage, pal_tf_family, pal_module, pal_fc_diverging, pal_mediation.
#
# All legacy names (GRBV_COLORS, GRBV_STAGE_COLOR, GRBV_SEASON_COLOR,
# GRBV_INFECTION_COLOR, theme_grbv, save_grbv) remain available so existing
# scripts keep working; only their hex values have shifted.

suppressPackageStartupMessages({
  library(ggplot2); library(scales); library(grid)
})

# Force UTF-8 so Unicode glyphs (é, ×) survive ggplot rendering on macOS
invisible(suppressWarnings(Sys.setlocale("LC_ALL", "en_US.UTF-8")))

# -----------------------------------------------------------------
# Dimensions (mm -> inches conversion handled by ggsave)
# -----------------------------------------------------------------
GRBV <- list(
  W_SINGLE_MM  =  89,
  W_ONEHALF_MM = 120,
  W_DOUBLE_MM  = 183,
  H_SHORT_MM   =  60,
  H_MEDIUM_MM  = 100,
  H_TALL_MM    = 160,
  H_PAGE_MM    = 240
)

# -----------------------------------------------------------------
# Palette — one hue family per semantic role
# -----------------------------------------------------------------
# Role table (keep in sync with manuscript):
#   Direction   red / blue          up / down regulation
#   Infection   vermilion / sky     RB(+) / RB(-)
#   Season      bluish-green / mag. 2021 / 2022
#   Stage       viridis 4-step      Pre -> Harvest (ordinal)
#   Mediation   blue / amber / grey titer / temperature / ns
#   Highlight   purple              accents only
#   TF family   Okabe-Ito 6-cat     HSF/NAC/bHLH/WRKY/MYB/ERF-AP2
# -----------------------------------------------------------------

GRBV_COLORS <- c(
  # Infection status (RB+ / RB-)           — Okabe-Ito vermilion / sky-blue
  RB_pos      = "#D55E00",
  RB_neg      = "#56B4E9",
  # Season                                  — Okabe-Ito bluish-green / magenta
  season_2021 = "#009E73",
  season_2022 = "#CC79A7",
  # Stage (ordinal, viridis 4-step)         — Pre -> Harvest
  stage_pre   = "#3B528B",
  stage_ver   = "#21908C",
  stage_post  = "#5DC863",
  stage_harv  = "#FDE725",
  # DE regulation direction                 — RdBu anchors
  up          = "#B2182B",
  down        = "#2166AC",
  neutral     = "#95A5A6",
  highlight   = "#8E44AD",
  # Mediation categorical                   — blue / amber / grey
  med_titer   = "#0072B2",
  med_temp    = "#E69F00",
  # Legacy aliases (kept so "warm" / "cool" still resolve)
  warm        = "#E69F00",
  cool        = "#2166AC"
)

# --- Semantic palettes (canonical; prefer these in new code) ---

pal_direction <- c(up   = unname(GRBV_COLORS["up"]),
                   down = unname(GRBV_COLORS["down"]),
                   ns   = unname(GRBV_COLORS["neutral"]))

pal_infection <- c(`RB(+)` = unname(GRBV_COLORS["RB_pos"]),
                   `RB(-)` = unname(GRBV_COLORS["RB_neg"]))

pal_season    <- c(`2021` = unname(GRBV_COLORS["season_2021"]),
                   `2022` = unname(GRBV_COLORS["season_2022"]))

pal_stage     <- c(`Pre-véraison`  = unname(GRBV_COLORS["stage_pre"]),
                   `Véraison`      = unname(GRBV_COLORS["stage_ver"]),
                   `Post-véraison` = unname(GRBV_COLORS["stage_post"]),
                   `Harvest`       = unname(GRBV_COLORS["stage_harv"]))

pal_tf_family <- c(HSF       = "#E69F00",
                   NAC       = "#0072B2",
                   bHLH      = "#009E73",
                   WRKY      = "#F0E442",
                   MYB       = "#CC79A7",
                   `ERF/AP2` = "#D55E00")

pal_module <- c(black       = "#2D3436", blue       = "#3498DB",
                brown       = "#A0522D", green      = "#2ECC71",
                turquoise   = "#17A2B8", yellow     = "#F1C40F",
                salmon      = "#FA8072", pink       = "#F48FB1",
                red         = "#E74C3C", purple     = "#8E44AD",
                magenta     = "#C2185B", greenyellow= "#B8E994",
                tan         = "#D4A373", grey       = "#95A5A6")

pal_fc_diverging <- c(low = unname(GRBV_COLORS["down"]),
                      mid = "#F7F7F7",
                      high= unname(GRBV_COLORS["up"]))

pal_mediation <- c(`Titer-mediated`   = unname(GRBV_COLORS["med_titer"]),
                   `Temperature-only` = unname(GRBV_COLORS["med_temp"]),
                   `Not significant`  = unname(GRBV_COLORS["neutral"]))

# --- Legacy stage/season/infection aliases (unchanged names) ---
GRBV_STAGES      <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
GRBV_STAGE_LABEL <- c(
  Pre_veraison  = "Pre-véraison",
  Veraison      = "Véraison",
  Post_veraison = "Post-véraison",
  Harvest       = "Harvest"
)
GRBV_STAGE_COLOR <- unname(GRBV_COLORS[c("stage_pre","stage_ver","stage_post","stage_harv")])
names(GRBV_STAGE_COLOR) <- GRBV_STAGES
GRBV_SEASON_COLOR <- c("2021" = unname(GRBV_COLORS["season_2021"]),
                       "2022" = unname(GRBV_COLORS["season_2022"]))
GRBV_INFECTION_COLOR <- c("pos" = unname(GRBV_COLORS["RB_pos"]),
                          "neg" = unname(GRBV_COLORS["RB_neg"]),
                          "RB(+)" = unname(GRBV_COLORS["RB_pos"]),
                          "RB(-)" = unname(GRBV_COLORS["RB_neg"]))

# -----------------------------------------------------------------
# Theme
# -----------------------------------------------------------------
theme_grbv <- function(base_size = 10, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text             = element_text(family = base_family, colour = "black"),
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                      hjust = 0, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size, colour = "grey25",
                                      margin = margin(b = 3)),
      plot.caption     = element_text(size = base_size - 2, colour = "grey35",
                                      hjust = 0),
      plot.tag         = element_text(face = "bold", size = base_size + 3),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      axis.line        = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 1.5),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      legend.key.height = unit(9, "pt"),
      legend.key.width  = unit(11, "pt"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size, margin = margin(b = 2)),
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      plot.margin      = margin(4, 6, 4, 6)
    )
}

# Tight publication variant (for Nature Plants column widths).  Base 8 pt.
theme_grbv_pub <- function(base_size = 8, base_family = "Arial") {
  theme_grbv(base_size = base_size, base_family = base_family) +
    theme(plot.margin = margin(3, 4, 3, 3, "mm"))
}

# -----------------------------------------------------------------
# Save helper: PNG (raster, 600 DPI) + SVG (vector, editable in Illustrator /
# Inkscape / Figma).  Default stack is PNG + SVG; pass formats = c("png","pdf")
# to keep the legacy Cairo PDF export.
# -----------------------------------------------------------------
save_grbv <- function(plot, stem, width_mm, height_mm,
                      formats = c("png","svg")) {
  for (fmt in formats) {
    out <- paste0(stem, ".", fmt)
    if (fmt == "pdf") {
      ggsave(out, plot = plot, width = width_mm, height = height_mm,
             units = "mm", device = cairo_pdf, limitsize = FALSE)
    } else if (fmt == "svg") {
      if (!requireNamespace("svglite", quietly = TRUE))
        stop("svglite is required for SVG export; install.packages('svglite')")
      ggsave(out, plot = plot, width = width_mm, height = height_mm,
             units = "mm", bg = "white", limitsize = FALSE,
             device = svglite::svglite)
    } else {   # png
      ggsave(out, plot = plot, width = width_mm, height = height_mm,
             units = "mm", dpi = 600, bg = "white", limitsize = FALSE,
             device = ragg::agg_png)
    }
    message(sprintf("wrote %s", out))
  }
}

# Explicit alias matching the skill's pub_save() naming
pub_save <- save_grbv
