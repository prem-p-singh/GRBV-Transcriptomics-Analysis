# Project summary / graphical abstract.
#
# A horizontal 5-roundel workflow infographic telling the full GRBV study at
# a glance, matching the same icon family used in Main10 and the MoA figure.
#
# Roundels (left -> right):
#   1  Field samples        : grapevine (leaf + berry cluster); sun + cloud
#                             as season context indicator
#   2  Sample prep          : test tube + gel
#   3  Sequencing / DE      : DNA-TF icon + a small bar-chart drawn in ggplot
#   4  Integrated findings  : hexagonal HUB larger than others, with 6 small
#                             glyphs around its rim (twin virion, sHSP ring,
#                             DNA-TF, sugar hexagon, cracked shield, blotched
#                             leaf); the central tile shows the Merlot berry
#   5  Disease phenotype    : red-blotched leaf + berry cluster; amber tint
#
# Roundels connected by thick curved arrows.
#
# Output: figures/main/Summary_study_design.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(grid); library(png); library(patchwork)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ICON_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/assets/icons"
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "Summary_study_design")

load_icon <- function(name) {
  f <- file.path(ICON_DIR, paste0(name, ".png"))
  if (!file.exists(f)) stop("icon missing: ", f)
  rasterGrob(readPNG(f), interpolate = TRUE)
}

ic <- list(
  leaf     = load_icon("leaf_healthy"),
  leaf_rb  = load_icon("leaf_redblotch"),
  berry    = load_icon("berry_cluster"),
  sun_cl   = load_icon("sun_cloud_pair"),
  tube_gel = load_icon("test_tube_gel"),
  dna_tf   = load_icon("dna_tf"),
  virion   = load_icon("virion_geminivirus"),
  shsp     = load_icon("shsp_oligomer"),
  sugar    = load_icon("sugar_glucose"),
  shield   = load_icon("shield_dampened"),
  nucleus  = load_icon("nucleus")
)

# ============================================================================
# Canvas and roundel centres (coord system 0..100 x, 0..30 y)
# ============================================================================
xlim <- c(0, 100); ylim <- c(0, 30)

roundels <- data.frame(
  id   = 1:5,
  cx   = c(10, 30, 50, 75, 95),
  cy   = c(17, 17, 17, 17, 17),
  r    = c( 6,  6,  6, 10,  6),                  # HUB (id 4) is larger
  fill = c(unname(GRBV_COLORS["season_2021"]),   # pale green
           "#cfdfc6",
           unname(GRBV_COLORS["med_titer"]),     # data blue
           "#F6EEDC",                            # cream HUB
           unname(GRBV_COLORS["warm"])),         # amber
  alpha_ = c(0.22, 0.22, 0.22, 0.30, 0.22),
  title = c("Field samples",
            "Sample prep",
            "Sequencing & DE",
            "Integrated findings",
            "Disease phenotype")
)

# Circle path generator
circle_df <- function(cx, cy, r, n = 180) {
  t <- seq(0, 2*pi, length.out = n)
  data.frame(x = cx + r*cos(t), y = cy + r*sin(t))
}

# Draw filled translucent discs first, then the outline in solid colour
p <- ggplot() +
  coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
  annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = 30,
           fill = "#FBF7EE", colour = NA) +
  theme_void() +
  theme(plot.margin = margin(2, 4, 2, 4))

for (i in seq_len(nrow(roundels))) {
  d <- circle_df(roundels$cx[i], roundels$cy[i], roundels$r[i])
  p <- p +
    annotate("polygon", x = d$x, y = d$y,
             fill = roundels$fill[i], alpha = roundels$alpha_[i],
             colour = NA) +
    annotate("path", x = d$x, y = d$y,
             colour = alpha(roundels$fill[i], 0.9), linewidth = 0.9)
}

# Connector arrows between consecutive roundels
arrow_xy <- function(i) {
  x1 <- roundels$cx[i]   + roundels$r[i]
  x2 <- roundels$cx[i+1] - roundels$r[i+1]
  y  <- roundels$cy[i]
  list(x = x1 + 0.3, xend = x2 - 0.3, y = y, yend = y)
}
for (i in 1:4) {
  a <- arrow_xy(i)
  p <- p +
    geom_segment(aes(x = !!a$x, y = !!a$y, xend = !!a$xend, yend = !!a$yend),
                 linewidth = 0.9, colour = unname(GRBV_COLORS["med_titer"]),
                 arrow = arrow(length = unit(8, "pt"), type = "closed"))
}

# ============================================================================
# ROUNDEL 1 — Field samples
#   healthy leaf top, berry cluster bottom-left, sun+cloud small top-right
# ============================================================================
c1 <- roundels[1, ]
p <- p +
  annotation_custom(ic$leaf,
                    xmin = c1$cx - 3.5, xmax = c1$cx + 0.5,
                    ymin = c1$cy + 0.2, ymax = c1$cy + 4.6) +
  annotation_custom(ic$berry,
                    xmin = c1$cx - 2.8, xmax = c1$cx + 1.2,
                    ymin = c1$cy - 4.4, ymax = c1$cy - 0.3) +
  annotation_custom(ic$sun_cl,
                    xmin = c1$cx + 0.8, xmax = c1$cx + 4.6,
                    ymin = c1$cy + 1.2, ymax = c1$cy + 3.6)

# ============================================================================
# ROUNDEL 2 — Sample prep : test tube + gel (single icon)
# ============================================================================
c2 <- roundels[2, ]
p <- p +
  annotation_custom(ic$tube_gel,
                    xmin = c2$cx - 4.5, xmax = c2$cx + 4.5,
                    ymin = c2$cy - 4.5, ymax = c2$cy + 4.5)

# ============================================================================
# ROUNDEL 3 — Sequencing & DE : DNA+TF + small bar chart
# ============================================================================
c3 <- roundels[3, ]
p <- p +
  annotation_custom(ic$dna_tf,
                    xmin = c3$cx - 5.0, xmax = c3$cx + 5.0,
                    ymin = c3$cy + 0.2, ymax = c3$cy + 3.8)

# Descending 3-bar chart inside lower half of roundel 3 (drawn in ggplot)
bar_df <- data.frame(
  x = c(c3$cx - 2.7, c3$cx - 0.4, c3$cx + 1.9),
  h = c(4.0, 2.6, 1.5)
)
for (i in seq_len(nrow(bar_df))) {
  p <- p +
    annotate("rect",
             xmin = bar_df$x[i] - 0.9, xmax = bar_df$x[i] + 0.9,
             ymin = c3$cy - 4.0, ymax = c3$cy - 4.0 + bar_df$h[i],
             fill = "white", colour = "grey15", linewidth = 0.5)
}

# ============================================================================
# ROUNDEL 4 — Integrated findings (HUB, larger)
#   Nucleus at center, 6 glyphs around at 60-degree spokes
# ============================================================================
c4 <- roundels[4, ]
# Central nucleus
p <- p +
  annotation_custom(ic$nucleus,
                    xmin = c4$cx - 2.2, xmax = c4$cx + 2.2,
                    ymin = c4$cy - 2.2, ymax = c4$cy + 2.2)

# 6 glyphs around the rim at angles starting from top, clockwise
hub_glyphs <- list(
  list(icon = ic$virion,   angle =  90),  # top
  list(icon = ic$shsp,     angle =  30),  # upper right
  list(icon = ic$dna_tf,   angle = -30),  # lower right
  list(icon = ic$sugar,    angle = -90),  # bottom
  list(icon = ic$shield,   angle = 210),  # lower left
  list(icon = ic$leaf_rb,  angle = 150)   # upper left
)
R_hub   <- c4$r * 0.72
halfsz <- 1.6
for (g in hub_glyphs) {
  th <- g$angle * pi / 180
  cx <- c4$cx + R_hub * cos(th)
  cy <- c4$cy + R_hub * sin(th)
  p <- p + annotation_custom(
    g$icon,
    xmin = cx - halfsz, xmax = cx + halfsz,
    ymin = cy - halfsz, ymax = cy + halfsz
  )
}

# ============================================================================
# ROUNDEL 5 — Disease phenotype
#   red-blotched leaf + berry (phenotype)
# ============================================================================
c5 <- roundels[5, ]
p <- p +
  annotation_custom(ic$leaf_rb,
                    xmin = c5$cx - 3.5, xmax = c5$cx + 0.5,
                    ymin = c5$cy + 0.2, ymax = c5$cy + 4.6) +
  annotation_custom(ic$berry,
                    xmin = c5$cx - 2.0, xmax = c5$cx + 2.2,
                    ymin = c5$cy - 4.4, ymax = c5$cy - 0.3)

# ============================================================================
# Roundel labels below each circle
# ============================================================================
lab_df <- roundels %>% mutate(y = cy - (r + 2))
p <- p +
  geom_text(data = lab_df,
            aes(x = cx, y = y, label = title),
            family = "Arial", fontface = "bold", size = 3.3,
            colour = "grey15")

# Roundel sub-captions (key numbers per stage)
sub_df <- data.frame(
  cx = roundels$cx,
  y  = roundels$cy - (roundels$r + 3.5),
  label = c(
    "327 samples \u00b7 4 stages\n2 seasons \u00b7 7 cohorts",
    "TAG-Seq RNA\nBradley N. Jenner (UCD)",
    "DESeq2 + voom\ncovariate-adjusted",
    "13 FARM pipeline rules\nTier 1 meta-analysis",
    "Red-blotch leaves;\nsugar-deficit berries"
  )
)
p <- p +
  geom_text(data = sub_df,
            aes(x = cx, y = y, label = label),
            family = "Arial", size = 2.3, colour = "grey30",
            lineheight = 1.05)

# ============================================================================
# COMPOSE with title + footer
# ============================================================================
title_block <- ggplot() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  annotate("text", x = 0.5, y = 0.74,
           label = "GRBV transcriptomics: study at a glance",
           family = "Arial", size = 5.6, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 0.5, y = 0.26,
           label = paste("Vitis vinifera cv. Merlot \u00b7 Oakville, CA \u00b7",
                         "2021 (cool) and 2022 (hot) seasons \u00b7",
                         "TAG-Seq \u2192 DESeq2 \u2192 multi-omics integration"),
           family = "Arial", size = 3.3, colour = "grey25", hjust = 0.5) +
  theme_void() + theme(plot.margin = margin(4, 4, 2, 4))

caption_block <- ggplot() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  annotate("text", x = 0.5, y = 0.5,
           label = paste("Headline findings: inverted seasonal trajectory \u00b7",
                         "sHSP axis dampened \u00b7 MBF1c + SPS1 mediators \u00b7",
                         "GRBV signature replicates in Zinfandel (Blanco-Ulate 2017)"),
           family = "Arial", size = 2.7, colour = "grey35", hjust = 0.5) +
  theme_void() + theme(plot.margin = margin(2, 4, 4, 4))

final <- title_block / p / caption_block +
  plot_layout(heights = c(0.55, 4.8, 0.35))

save_grbv(final, OUT_STEM,
          width_mm = 260, height_mm = 125)

cat("[summary] figure written to ", OUT_STEM, ".png / .svg\n", sep = "")
