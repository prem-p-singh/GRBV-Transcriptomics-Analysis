# Mechanism-of-action figure for GRBV.
#
# Composition: single horizontal 3-zone scene
#   LEFT  : vineyard / field context
#           leaf_healthy, leaf_redblotch, berry_cluster, leafhopper_vector
#   MIDDLE: phloem entry
#           zoom-frame circle containing phloem_tissue_crosssection +
#           virion_geminivirus entering
#   RIGHT : host cell / molecular mechanism
#           zoom-frame circle containing nucleus + rolling_circle_replication;
#           around the nucleus: shield_dampened (faded),
#           shsp_oligomer (faded), dna_tf with small sugar_glucose annotation
#
# Icons live in figures/assets/icons/ (all Firefly PNGs with
# transparent or near-transparent backgrounds).
#
# Output: figures/main/MoA_GRBV_mechanism.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(grid); library(png); library(patchwork)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ICON_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/assets/icons"
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "MoA_GRBV_mechanism")

load_icon <- function(name) {
  f <- file.path(ICON_DIR, paste0(name, ".png"))
  if (!file.exists(f)) stop("icon missing: ", f)
  rasterGrob(readPNG(f), interpolate = TRUE)
}

icons <- list(
  leaf      = load_icon("leaf_healthy"),
  leaf_rb   = load_icon("leaf_redblotch"),
  berry     = load_icon("berry_cluster"),
  hopper    = load_icon("leafhopper_vector"),
  phloem    = load_icon("phloem_tissue_crosssection"),
  virion    = load_icon("virion_geminivirus"),
  nucleus   = load_icon("nucleus"),
  rcr       = load_icon("rolling_circle_replication"),
  shield    = load_icon("shield_dampened"),
  shsp      = load_icon("shsp_oligomer"),
  dna_tf    = load_icon("dna_tf"),
  sugar     = load_icon("sugar_glucose")
)

# ============================================================================
# BASE CANVAS
# ============================================================================
# Coord system 0..30 x 0..10 y, landscape. Three zones ~10 units wide each.
xlim <- c(0, 30); ylim <- c(0, 10)

zone_divider <- data.frame(x = c(10, 20), xend = c(10, 20),
                           y = 0.3, yend = 9.3)
zone_labels <- data.frame(
  x  = c(5, 15, 25),
  y  = 9.6,
  label = c("A  Vector delivers virus to the vine",
            "B  Geminate virion enters the phloem",
            "C  Host cell is reprogrammed")
)

# Zoom-frame circles for zones B and C
zoom_circ <- data.frame(
  cx = c(15, 25), cy = c(4.5, 4.5), r = 3.5
)

base <- ggplot() +
  coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
  # Soft green background band covering full canvas
  annotate("rect", xmin = 0, xmax = 30, ymin = 0, ymax = 10,
           fill = "#F4F7F1", colour = NA) +
  # Thin vertical dividers between zones
  geom_segment(data = zone_divider,
               aes(x = x, y = y, xend = xend, yend = yend),
               colour = "grey70", linewidth = 0.3, linetype = "12") +
  # Zone headers
  geom_text(data = zone_labels,
            aes(x = x, y = y, label = label),
            family = "Arial", fontface = "bold", size = 3.5,
            colour = "grey15") +
  theme_void()

# ============================================================================
# ZONE A (LEFT) - vineyard / vector context
# ============================================================================
# Correction 2026-04-23: panel A now shows ONLY the pre-symptomatic vine
# (healthy leaves + unripe grapes + treehopper vector feeding on the plant).
# Red-blotch phenotype belongs in panel C (outcome), not here.
#
# Layout:
#   top-left        : healthy leaf (main)
#   top-right       : second healthy leaf (smaller, angled) so the vine looks
#                     like a plant rather than a single leaf
#   middle-right    : treehopper vector, enlarged and positioned as if feeding
#                     on the leaf edge (mouthpart pointing at the leaf)
#   bottom-center   : grape cluster (green/unripe appearance comes from icon)
p <- base +
  # Primary healthy leaf, upper-left of zone A
  annotation_custom(icons$leaf,
                    xmin = 0.3, xmax = 4.8, ymin = 4.6, ymax = 8.8) +
  # Secondary small healthy leaf, upper-right, implies "vine" not "leaf"
  annotation_custom(icons$leaf,
                    xmin = 5.6, xmax = 8.2, ymin = 6.0, ymax = 8.5) +
  # Berry cluster, bottom-left
  annotation_custom(icons$berry,
                    xmin = 0.7, xmax = 4.6, ymin = 0.8, ymax = 4.6) +
  # Leafhopper vector, enlarged, positioned to the RIGHT of the main leaf
  # with its pointed mouthpart touching the leaf edge (feeding posture)
  annotation_custom(icons$hopper,
                    xmin = 4.6, xmax = 9.2, ymin = 2.4, ymax = 6.2) +
  # Small red stylized "virus transmission" dot at the feeding contact point
  annotate("point", x = 5.0, y = 5.2, size = 1.5,
           colour = unname(GRBV_COLORS["RB_pos"]), shape = 16) +
  # Thin feeding line from hopper mouthpart to leaf (highlights the
  # piercing contact, NOT a curved "from-to" arrow)
  geom_segment(aes(x = 5.2, y = 5.3, xend = 4.7, yend = 5.6),
               linewidth = 0.4, colour = "grey25",
               linetype = "solid",
               arrow = arrow(length = unit(3, "pt"), type = "closed"))

# ============================================================================
# ZONE B (MIDDLE) - phloem + virion entering
# ============================================================================
# Draw the zoom-frame circle first so icons sit on top of the pale fill
zone_B_circ <- zoom_circ[1, ]
p <- p +
  annotate("path",
           x = zone_B_circ$cx + zone_B_circ$r * cos(seq(0, 2*pi, length.out = 120)),
           y = zone_B_circ$cy + zone_B_circ$r * sin(seq(0, 2*pi, length.out = 120)),
           colour = "grey40", linewidth = 0.45) +
  annotate("point", x = zone_B_circ$cx, y = zone_B_circ$cy,
           size = 0, alpha = 0) +
  # Phloem tissue cross-section, inside the zoom circle, slightly to the right
  annotation_custom(icons$phloem,
                    xmin = 14.2, xmax = 17.2, ymin = 2.6, ymax = 6.6) +
  # Virion entering from the left edge of the circle
  annotation_custom(icons$virion,
                    xmin = 11.8, xmax = 14.0, ymin = 4.0, ymax = 5.4) +
  # Arrow: virion enters phloem
  geom_segment(aes(x = 14.0, y = 4.8, xend = 14.8, yend = 4.8),
               linewidth = 0.5, colour = "grey20",
               arrow = arrow(length = unit(5, "pt"), type = "closed"))

# Small connector between zone A (leafhopper mouthpart) and zone B (virion)
p <- p +
  geom_curve(aes(x = 9.3, y = 3.2, xend = 12.0, yend = 4.6),
             curvature = 0.15, linewidth = 0.45, colour = "grey45",
             linetype = "22",
             arrow = arrow(length = unit(4, "pt"), type = "closed"))

# ============================================================================
# ZONE C (RIGHT) - host cell / molecular mechanism
# ============================================================================
zone_C_circ <- zoom_circ[2, ]
p <- p +
  annotate("path",
           x = zone_C_circ$cx + zone_C_circ$r * cos(seq(0, 2*pi, length.out = 120)),
           y = zone_C_circ$cy + zone_C_circ$r * sin(seq(0, 2*pi, length.out = 120)),
           colour = "grey40", linewidth = 0.45) +
  # Nucleus center of zone C
  annotation_custom(icons$nucleus,
                    xmin = 23.5, xmax = 26.5, ymin = 4.2, ymax = 7.0) +
  # Rolling-circle replication, overlapping nucleus to show it happens inside
  annotation_custom(icons$rcr,
                    xmin = 23.9, xmax = 26.1, ymin = 4.5, ymax = 6.5) +
  # Shield (dampened) top-left of nucleus
  annotation_custom(icons$shield,
                    xmin = 21.2, xmax = 23.1, ymin = 5.6, ymax = 7.3) +
  # sHSP oligomer bottom-left of nucleus
  annotation_custom(icons$shsp,
                    xmin = 21.2, xmax = 23.0, ymin = 2.3, ymax = 4.0) +
  # DNA + TF bottom-right of nucleus
  annotation_custom(icons$dna_tf,
                    xmin = 25.8, xmax = 28.2, ymin = 2.2, ymax = 3.8) +
  # Sugar glucose top-right of nucleus
  annotation_custom(icons$sugar,
                    xmin = 26.4, xmax = 27.9, ymin = 5.5, ymax = 7.0) +
  # Disease OUTCOME: small red-blotched leaf in the lower-right of the
  # zoom circle, representing the phenotypic consequence of nuclear
  # reprogramming. Kept small so it reads as a downstream output.
  annotation_custom(icons$leaf_rb,
                    xmin = 25.8, xmax = 28.4, ymin = 1.2, ymax = 3.0)

# DCL2 antiviral-silencing thread (added 2026-05-19): a labelled callout in
# the empty space below the nucleus, with a downward curved arrow from the
# bottom of the nucleus to the top of the callout. Represents the
# veraison-specific antiviral RNA-silencing response confirmed by
# cross-cultivar replication against Rumbaugh 2022 (S96) and Rumbaugh 2023
# (S97).
p <- p +
  annotate("rect", xmin = 23.0, xmax = 25.7, ymin = 1.5, ymax = 2.8,
           fill = "#E8F0E4", colour = unname(GRBV_COLORS["RB_neg"]),
           linewidth = 0.45, alpha = 0.85) +
  annotate("text", x = 24.35, y = 2.45,
           label = "DCL2 ↑ at veraison",
           family = "Arial", fontface = "bold", size = 2.6,
           colour = "grey15") +
  annotate("text", x = 24.35, y = 1.85,
           label = "antiviral RNA silencing",
           family = "Arial", size = 2.2, colour = "grey25")

# Small arrows from nucleus outward to the five mechanism glyphs (dashed
# to convey "downstream effect", not direct physical link)
p <- p +
  geom_curve(aes(x = 24.3, y = 5.8, xend = 22.7, yend = 6.3),   # nucleus -> shield
             curvature = 0.10, linewidth = 0.35,
             colour = "grey45", linetype = "13",
             arrow = arrow(length = unit(3, "pt"))) +
  geom_curve(aes(x = 24.3, y = 5.3, xend = 22.6, yend = 3.5),   # nucleus -> sHSP
             curvature = -0.10, linewidth = 0.35,
             colour = "grey45", linetype = "13",
             arrow = arrow(length = unit(3, "pt"))) +
  geom_curve(aes(x = 25.7, y = 5.3, xend = 26.8, yend = 3.5),   # nucleus -> DNA-TF
             curvature = 0.12, linewidth = 0.35,
             colour = "grey45", linetype = "13",
             arrow = arrow(length = unit(3, "pt"))) +
  geom_curve(aes(x = 25.8, y = 5.9, xend = 27.0, yend = 6.1),   # nucleus -> sugar
             curvature = -0.08, linewidth = 0.35,
             colour = "grey45", linetype = "13",
             arrow = arrow(length = unit(3, "pt"))) +
  geom_curve(aes(x = 25.0, y = 4.3, xend = 24.45, yend = 2.85),  # nucleus -> DCL2
             curvature = 0.06, linewidth = 0.4,
             colour = unname(GRBV_COLORS["RB_neg"]), linetype = "13",
             arrow = arrow(length = unit(3, "pt")))

# Small connector between zone B circle and zone C circle
p <- p +
  geom_curve(aes(x = 18.2, y = 4.6, xend = 21.5, yend = 4.8),
             curvature = 0.18, linewidth = 0.45, colour = "grey45",
             linetype = "22",
             arrow = arrow(length = unit(4, "pt"), type = "closed"))

# ============================================================================
# CAPTIONS UNDER EACH ZONE
# ============================================================================
# Zone captions rendered as a separate short row beneath the scene (so they
# don't collide with icons in zones A and C).  Built as its own ggplot on a
# matching 0..30 x-axis so text positions line up with the scene zones above.
caption_plot <- ggplot() +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 1), expand = FALSE) +
  annotate("text", x = 5, y = 0.55,
           label = "Three-cornered alfalfa treehopper feeds\non a healthy Merlot vine and delivers GRBV",
           family = "Arial", size = 2.8, colour = "grey25", lineheight = 1.05) +
  annotate("text", x = 15, y = 0.55,
           label = "Twinned geminate virion enters\nphloem sieve elements",
           family = "Arial", size = 2.8, colour = "grey25", lineheight = 1.05) +
  annotate("text", x = 25, y = 0.55,
           label = "Rolling-circle replication in the nucleus\nreprograms the host transcriptome,\nproducing red-blotch phenotype",
           family = "Arial", size = 2.8, colour = "grey25", lineheight = 1.05) +
  theme_void() + theme(plot.margin = margin(0, 6, 2, 6))

# ============================================================================
# TITLE + COMPOSE
# ============================================================================
title_block <- ggplot() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  annotate("text", x = 0.5, y = 0.72,
           label = "Mechanism of action of Grapevine Red Blotch Virus",
           family = "Arial", size = 5.4, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 0.5, y = 0.25,
           label = paste("Insect vector \u2192 phloem entry \u2192",
                         "nuclear replication \u2192 rewired host transcriptome \u2192",
                         "red-blotch phenotype"),
           family = "Arial", size = 3.2, colour = "grey25", hjust = 0.5) +
  theme_void() + theme(plot.margin = margin(4, 4, 2, 4))

caption_block <- ggplot() +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  annotate("text", x = 0.5, y = 0.5,
           label = paste("Downstream mechanism hits: DCL2 (LOC100258190) antiviral RNA-silencing peaks at veraison",
                         "and replicates cross-cultivar (Rumbaugh 2022, 2023);",
                         "MBF1c (LOC100249249) candidate mediator; SPS1 (LOC100232974) sugar-metabolism marker;",
                         "DOF3.6 / DOF1.7 motif enrichment; 9 of 10 defence categories remain silent."),
           family = "Arial", size = 2.7, colour = "grey35", hjust = 0.5) +
  theme_void() + theme(plot.margin = margin(2, 4, 4, 4))

final <- title_block / p / caption_plot / caption_block +
  plot_layout(heights = c(0.55, 4.8, 0.75, 0.35))

save_grbv(final, OUT_STEM,
          width_mm = 260, height_mm = 180)

cat("[MoA] figure written to ", OUT_STEM, ".png / .svg\n", sep = "")
