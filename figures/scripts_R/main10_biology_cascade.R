# Main 10 v2 — GRBV molecular cascade, upgraded to a true biology schematic.
#
# Redesign brief (SciAgent-Skills review):
#   * Original Main10 read as a text flowchart (5 stacked boxes).
#   * v2 keeps the 5-step narrative but now has:
#       - A tissue-level banner across the top (leaf -> virion -> phloem ->
#         nucleus), so the reader immediately sees where GRBV acts.
#       - Each of the 5 cascade boxes carries an inline icon on the left
#         conveying the mechanism at a glance (virion / shield / sHSP /
#         DNA-TF+sugar / red-blotched leaf).
#
# Icon assets are Adobe-Firefly PNGs in figures/assets/icons/.
# No magick dependency: uses png::readPNG + grid::rasterGrob via
# ggplot2::annotation_custom.
#
# Output: Main10_biology_cascade.png + .svg in figures/main/.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(grid); library(png); library(patchwork)
})

HERE       <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ICON_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/assets/icons"
OUT_DIR    <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM   <- file.path(OUT_DIR, "Main10_biology_cascade")

load_icon <- function(name) {
  f <- file.path(ICON_DIR, paste0(name, ".png"))
  if (!file.exists(f)) stop("icon missing: ", f)
  rasterGrob(readPNG(f), interpolate = TRUE)
}

icons <- list(
  leaf_healthy      = load_icon("leaf_healthy"),
  virion            = load_icon("virion_geminivirus"),
  phloem            = load_icon("phloem_sieve"),
  nucleus           = load_icon("nucleus"),
  shield            = load_icon("shield_dampened"),
  shsp              = load_icon("shsp_oligomer"),
  dna_tf            = load_icon("dna_tf"),
  sugar             = load_icon("sugar_glucose"),
  leaf_redblotch    = load_icon("leaf_redblotch")
)

# ============================================================================
# PART A — TISSUE BANNER (leaf -> virion -> phloem -> nucleus)
# ============================================================================
# Use a 0..10 x, 0..3 y coord system for easy sizing
banner_nodes <- data.frame(
  cx = c(1.25, 3.75, 6.25, 8.75),
  cy = 1.7,
  icon  = c("leaf_healthy", "virion", "phloem", "nucleus"),
  label = c("Infected leaf", "Geminivirus virion",
            "Phloem sieve element", "Host cell nucleus")
)
ic_halfw <- 0.80; ic_halfh <- 0.95

arrow_segs <- data.frame(
  x  = banner_nodes$cx[-4] + ic_halfw + 0.05,
  xe = banner_nodes$cx[-1] - ic_halfw - 0.05,
  y  = banner_nodes$cy[-4], ye = banner_nodes$cy[-1]
)

p_banner <- ggplot() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 3), expand = FALSE)

# Insert each icon via annotation_custom (data coords)
for (i in seq_len(nrow(banner_nodes))) {
  p_banner <- p_banner + annotation_custom(
    icons[[banner_nodes$icon[i]]],
    xmin = banner_nodes$cx[i] - ic_halfw,
    xmax = banner_nodes$cx[i] + ic_halfw,
    ymin = banner_nodes$cy[i] - ic_halfh,
    ymax = banner_nodes$cy[i] + ic_halfh
  )
}

p_banner <- p_banner +
  geom_segment(data = arrow_segs,
               aes(x = x, y = y, xend = xe, yend = ye),
               linewidth = 0.7, color = "grey25",
               arrow = arrow(length = unit(6, "pt"), type = "closed")) +
  geom_text(data = banner_nodes,
            aes(x = cx, y = 0.30, label = label),
            family = "Arial", size = 3.0, colour = "grey20") +
  annotate("text", x = 0.1, y = 2.85,
           label = "A  Tissue-level view: where GRBV acts",
           family = "Arial", size = 3.6, fontface = "bold",
           hjust = 0, colour = "black") +
  theme_void() +
  theme(plot.margin = margin(2, 4, 2, 4))

# ============================================================================
# PART B — 5-STEP MECHANISM CASCADE WITH INLINE ICONS
# ============================================================================
box_h <- 1.60
y_centers <- seq(10.8 - box_h/2, 1.2 + box_h/2, length.out = 5)
station_fill <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
station_text <- c("white",   "white",   "white",   "grey10",  "grey10")

nodes <- data.frame(
  y        = y_centers,
  fill     = station_fill,
  text_col = station_text,
  icon     = c("virion", "shield", "shsp", "dna_tf", "leaf_redblotch"),
  title = c(
    "1. GRBV infection establishes",
    "2. Canonical plant immunity dampened",
    "3. HSF \u2194 sHSP thermal-defence axis suppressed",
    "4. Sugar-metabolism transcription rewired",
    "5. Disease phenotype: red blotch + sugar deficit"),
  body = c(
    "Geminivirus replicates in phloem;\nN761_gp1\u20133 dominate the Q6 volcano (Main 2);\nbiphasic titer dynamics across stages (Main 3, panel C)",
    "All 10 immune categories show |log2FC| < 0.08\n(NLR, PRR, PR, SA, JA, ET, WRKY, callose, ROS, R-adapters);\nWGCNA modules preserved across infected vines (Zsum 5.8\u201336.4)",
    "29 / 36 sHSPs down-regulated in infected vines (Main 4, panel A);\nHSF \u2194 sHSP coordinated, r = 0.625 (Main 4, panel B);\nMBF1c (LOC100249249) is the cross-pipeline robust\ntemperature \u2192 titer mediator (Main 5)",
    "DOF3.6 / DOF1.7 motifs enriched in brown + yellow modules (Main 7);\nBrix suppressed 4\u20139.8 \u00b0Bx in infected berries, with TA / Malic\ninverted between 2021 cool and 2022 hot seasons (S63);\ntop host biomarker maps to sugar-accumulation QTL on chr 1 (Main 8)",
    "Visible red blotching;  reduced \u00b0Brix and yield;\n\u22651,500 stage-specific DEGs at Harvest in cool seasons (Main 3, panel A);\nbiomarker panel co-localises with known grape QTL\n(Main 8, panel C;  2.0\u00d7 enrichment over a 5,000-permutation null)"
  ),
  stringsAsFactors = FALSE
)

# Icon column on the left; text column on the right
icon_x_left <- 0.35; icon_x_right <- 2.15
text_x_left <- 2.35; text_x_right <- 9.60

box <- data.frame(
  xmin = 0.2, xmax = 9.8,
  ymin = nodes$y - box_h/2, ymax = nodes$y + box_h/2,
  fill = nodes$fill
)

arrows <- data.frame(
  x = 5, xend = 5,
  y = head(nodes$y, -1) - box_h/2,
  yend = tail(nodes$y, -1) + box_h/2
)

p_cascade <- ggplot() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 12), expand = FALSE) +
  geom_rect(data = box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = box$fill, color = "black", linewidth = 0.45, alpha = 0.95) +
  # White tile behind each icon so colourful icons read cleanly over dark fills
  geom_rect(data = nodes,
            aes(xmin = icon_x_left, xmax = icon_x_right,
                ymin = y - box_h/2 + 0.15, ymax = y + box_h/2 - 0.15),
            fill = "white", color = NA)

# Embed each row icon
for (i in seq_len(nrow(nodes))) {
  icon_name <- nodes$icon[i]
  # Slight shrink to sit inside the white tile with a margin
  p_cascade <- p_cascade + annotation_custom(
    icons[[icon_name]],
    xmin = icon_x_left + 0.07,
    xmax = icon_x_right - 0.07,
    ymin = nodes$y[i] - box_h/2 + 0.22,
    ymax = nodes$y[i] + box_h/2 - 0.22
  )
}

# Box 4 (sugar rewiring): overlay a small sugar glucose icon in the
# bottom-right corner of the icon tile so the glyph reads "TF + sugar"
p_cascade <- p_cascade + annotation_custom(
  icons$sugar,
  xmin = icon_x_right - 0.65, xmax = icon_x_right - 0.07,
  ymin = nodes$y[4] - box_h/2 + 0.18,
  ymax = nodes$y[4] - box_h/2 + 0.78
)

p_cascade <- p_cascade +
  geom_text(data = nodes,
            aes(x = text_x_left, y = y + box_h/2 - 0.28, label = title),
            family = "Arial", color = nodes$text_col, fontface = "bold",
            hjust = 0, vjust = 1, size = 3.5) +
  geom_text(data = nodes,
            aes(x = text_x_left, y = y + box_h/2 - 0.58, label = body),
            family = "Arial", color = nodes$text_col, hjust = 0,
            vjust = 1, size = 2.6, lineheight = 1.08) +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.7, color = "grey25",
               arrow = arrow(length = unit(7, "pt"), type = "closed")) +
  annotate("text", x = 0.1, y = 11.85,
           label = "B  Mechanism cascade",
           family = "Arial", size = 3.6, fontface = "bold",
           hjust = 0, colour = "black") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(2, 4, 2, 4))

# ============================================================================
# COMPOSE
# ============================================================================
title_plot <- ggplot() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  annotate("text", x = 0.5, y = 0.75, label = "GRBV molecular disruption cascade",
           family = "Arial", size = 5.2, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 0.5, y = 0.28,
           label = paste("Viral establishment \u2192 silent immunity",
                         "\u2192 thermal-defence dampening \u2192",
                         "sugar-pathway rewiring \u2192 disease phenotype"),
           family = "Arial", size = 3.2, colour = "grey25", hjust = 0.5) +
  theme_void() + theme(plot.margin = margin(4, 4, 2, 4))

caption_plot <- ggplot() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  annotate("text", x = 0.5, y = 0.5,
           label = paste("Robust temp \u2192 titer mediators (multi-pipeline):",
                         "MBF1c (LOC100249249), SPS1 (LOC100232974),",
                         "uncharacterised LOC100266309"),
           family = "Arial", size = 2.8, colour = "grey35", hjust = 0.5) +
  theme_void() + theme(plot.margin = margin(2, 4, 4, 4))

final <- title_plot / p_banner / p_cascade / caption_plot +
  plot_layout(heights = c(0.45, 1.40, 5.00, 0.30))

save_grbv(final, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 260)

cat("[main10_v2] cascade diagram with icons complete\n")
