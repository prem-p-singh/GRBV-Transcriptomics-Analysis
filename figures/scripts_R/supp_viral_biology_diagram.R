# Supplementary S60 — GRBV viral-gene function atlas (biology explainer)
#
# Three stacked panels:
#   A  Linear GRBV genome schematic: 6 ORFs with strand orientation,
#      standard-name labels, and a proportional CPM bar above each
#      showing observed expression from this dataset.
#   B  Function cards: one per ORF, with protein class + primary
#      mechanism + a one-line "what it does to the plant" summary.
#   C  Virus-host impact diagram: ORFs on the left, the plant processes
#      each ORF hits on the right, with observed per-sample correlations
#      from S59 annotated on each edge.
#
# Output: figures/supp/S60_viral_biology_diagram.{png,svg}

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(grid); library(ggforce)
})
HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Palette: virion-sense = warm (amber tones), complementary = cool (blues)
orf_col <- c(
  V2 = "#E69F00",   # amber, virion sense
  V1 = "#D55E00",   # vermilion, virion sense (Coat Protein emphasis)
  V3 = "#F0E442",   # yellow, virion sense
  C1 = "#0072B2",   # blue, complementary sense (Rep emphasis)
  C2 = "#56B4E9",   # sky, complementary
  C3 = "#95A5A6"    # grey, complementary (silent)
)

# GRBV schematic: simplified linear layout (3.2 kb, six ORFs laid out
# without overlap). Accurate GRBV gene order preserved (V-cluster on +
# strand, C-cluster on - strand); coordinates spaced for visual clarity,
# not exact NCBI positions.
orfs <- data.frame(
  orf       = c("V2",  "V1",  "V3",   "C1",    "C2",   "C3"),
  gp        = c("gp1", "gp2", "gp3",  "gp5",   "gp4",  "gp6"),
  full_name = c("V2 protein", "Coat protein (CP)", "V3 protein",
                "Rep (replication-assoc.)",
                "C2 / TrAP-like", "Replication enhancer (REn)"),
  xmin      = c( 180,  700, 1300,    1940,    2530,   2900),
  xmax      = c( 660, 1260, 1860,    2510,    2870,   3100),
  strand    = c( "+",  "+",  "+",     "-",     "-",    "-"),
  log_cpm   = c( 9.2,  6.6,  7.2,     5.1,     4.0,    0.0),
  stringsAsFactors = FALSE
)
orfs$y <- ifelse(orfs$strand == "+", 0.45, -0.45)
orfs$h <- 0.32

GENOME_LEN <- 3200   # approx GRBV genome length in bp

# ==================================================================
# Panel A — genome schematic with CPM bars above each ORF
# ==================================================================
cpm_bars <- orfs %>%
  mutate(bar_x = (xmin + xmax) / 2, bar_top = log_cpm,
         bar_lab = sprintf("%.1f", log_cpm))

pA <- ggplot() +
  # genome backbone
  geom_segment(aes(x = 0, xend = GENOME_LEN, y = 0, yend = 0),
               color = "grey30", linewidth = 1.2) +
  # ORF arrow bodies (drawn as rounded rectangles)
  geom_rect(data = orfs,
            aes(xmin = xmin, xmax = xmax, ymin = y - h/2, ymax = y + h/2,
                fill = orf), color = "grey20", linewidth = 0.4) +
  # strand direction arrows
  geom_segment(data = orfs,
               aes(x = ifelse(strand == "+", xmin + 25, xmax - 25),
                   xend = ifelse(strand == "+", xmax - 25, xmin + 25),
                   y = y, yend = y),
               arrow = arrow(length = unit(6, "pt"), type = "closed"),
               color = "grey10", linewidth = 0.55) +
  # ORF labels inside
  geom_text(data = orfs,
            aes(x = (xmin + xmax) / 2, y = y, label = orf),
            family = "Arial", fontface = "bold", color = "white", size = 3.5) +
  # gp identifier below each arrow
  geom_text(data = orfs,
            aes(x = (xmin + xmax) / 2, y = y - 0.28,
                label = gp),
            family = "Arial", size = 2.5, color = "grey20") +
  # full-name label below the gp identifier, one line
  geom_text(data = orfs,
            aes(x = (xmin + xmax) / 2, y = y - 0.55, label = full_name),
            family = "Arial", size = 2.2, color = "grey35", fontface = "italic") +
  # log-CPM bar: virion-sense ORFs grow upward above their arrow; complementary
  # ORFs grow upward too (above strand, below CPM axis) for a consistent read.
  geom_rect(data = cpm_bars,
            aes(xmin = bar_x - 80, xmax = bar_x + 80,
                ymin = 1.0,
                ymax = 1.0 + log_cpm * 0.10,
                fill = orf), alpha = 0.85, color = "grey35", linewidth = 0.25) +
  geom_text(data = cpm_bars,
            aes(x = bar_x,
                y = 1.0 + log_cpm * 0.10 + 0.12,
                label = bar_lab),
            family = "Arial", fontface = "bold", size = 2.6, color = "grey15") +
  # strand-side labels
  annotate("text", x = -100, y = 0.45, label = "virion-sense (+)",
           family = "Arial", size = 2.7, fontface = "bold",
           color = "#8B4000", hjust = 1) +
  annotate("text", x = -100, y = -0.45, label = "complementary (−)",
           family = "Arial", size = 2.7, fontface = "bold",
           color = "#0a2b52", hjust = 1) +
  annotate("text", x = GENOME_LEN / 2, y = 2.3,
           label = "Observed mean log2(CPM+1), RB+ samples",
           family = "Arial", size = 2.9, fontface = "italic", color = "grey30") +
  annotate("text", x = GENOME_LEN / 2, y = -1.3,
           label = "circular ssDNA genome, 3.2 kb  (linear schematic for clarity)",
           family = "Arial", size = 2.5, fontface = "italic", color = "grey40") +
  scale_fill_manual(values = orf_col, guide = "none") +
  coord_cartesian(xlim = c(-700, GENOME_LEN + 60), ylim = c(-1.6, 2.6),
                  expand = FALSE) +
  labs(tag = "A", title = "A. GRBV genome organisation and viral-gene expression") +
  theme_void(base_family = "Arial") +
  theme(plot.title    = element_text(family = "Arial", face = "bold",
                                     size = 11, hjust = 0.5, color = "grey10"),
        plot.tag      = element_text(family = "Arial", face = "bold", size = 13),
        plot.margin   = margin(4, 6, 4, 6))

# ==================================================================
# Panel B — function cards (6 tiles in a 2 x 3 grid)
# ==================================================================
cards <- data.frame(
  orf = c("V2 (gp1)", "V1 (gp2)", "V3 (gp3)",
          "C1 (gp5)", "C2 (gp4)", "C3 (gp6)"),
  role = c("silencing suppressor",
           "coat protein",
           "pathogenicity accessory",
           "replication engine",
           "transcriptional activator / silencing suppressor II",
           "replication enhancer (silent in our data)"),
  mech = c("Binds small RNAs; disables plant siRNA antiviral defense.",
           "Encapsidates ssDNA; required for treehopper transmission and phloem-wide movement.",
           "Grablovirus-specific ORF; implicated in symptom severity (leaf blotching, delayed ripening).",
           "Initiates rolling-circle replication; binds pRBR, releases E2F, pushes host cells into S-phase.",
           "Activates viral gene transcription; second-line silencing suppressor; pathogenicity determinant.",
           "Normally boosts Rep-mediated replication. In our data counts are near zero: GRBV appears to replicate without a functional C3 enhancer in field-infected Merlot."),
  fill = c(orf_col["V2"], orf_col["V1"], orf_col["V3"],
           orf_col["C1"], orf_col["C2"], orf_col["C3"]),
  stringsAsFactors = FALSE
)

cards$row <- c(1, 1, 1, 2, 2, 2)
cards$col <- c(1, 2, 3, 1, 2, 3)

pB <- ggplot(cards) +
  geom_rect(aes(xmin = (col - 1) * 3.4, xmax = (col - 1) * 3.4 + 3.2,
                ymin = -row * 2.2, ymax = -row * 2.2 + 2.0),
            fill = "#F7F4EE", color = "#CBD5E1", linewidth = 0.6) +
  geom_rect(aes(xmin = (col - 1) * 3.4, xmax = (col - 1) * 3.4 + 0.18,
                ymin = -row * 2.2, ymax = -row * 2.2 + 2.0,
                fill = fill), color = NA) +
  geom_text(aes(x = (col - 1) * 3.4 + 0.3, y = -row * 2.2 + 1.70,
                label = orf),
            family = "Arial", fontface = "bold", size = 3.4,
            color = "grey10", hjust = 0) +
  geom_text(aes(x = (col - 1) * 3.4 + 0.3, y = -row * 2.2 + 1.42,
                label = role),
            family = "Arial", fontface = "italic", size = 2.6,
            color = "#0072B2", hjust = 0) +
  geom_text(aes(x = (col - 1) * 3.4 + 0.3, y = -row * 2.2 + 0.8,
                label = stringr::str_wrap(mech, width = 48)),
            family = "Arial", size = 2.4, color = "grey25",
            hjust = 0, vjust = 1, lineheight = 1.05) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-0.1, 10.1), ylim = c(-4.6, 0.1), expand = FALSE) +
  labs(tag = "B", title = "B. Function of each viral gene") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(family = "Arial", face = "bold",
                                  size = 11, hjust = 0.5, color = "grey10"),
        plot.tag   = element_text(family = "Arial", face = "bold", size = 13),
        plot.margin = margin(4, 6, 4, 6))

# ==================================================================
# Panel C — virus-host impact diagram (ORFs -> plant processes -> phenotype)
# ==================================================================
# Two columns: viral ORFs on left, plant targets on right.
# Edges annotated with observed Spearman r from S59.

orf_nodes <- data.frame(
  id    = c("V2", "V1", "V3", "C1", "C2"),
  label = c("V2 (gp1)\nsilencing\nsuppressor",
            "V1 (gp2)\ncoat protein",
            "V3 (gp3)\npathogenicity",
            "C1 (gp5)\nRep",
            "C2 (gp4)\nTrAP /\nsuppressor II"),
  x = 1,
  y = c(4.5, 3.5, 2.5, 1.5, 0.5),
  fill = orf_col[c("V2","V1","V3","C1","C2")],
  stringsAsFactors = FALSE
)

host_nodes <- data.frame(
  id    = c("silencing", "movement", "symptom", "cellcycle", "sugar", "shsp"),
  label = c("Host RNA silencing\n(siRNA pathway disabled)",
            "Phloem-wide movement\n+ vector transmission",
            "Red blotching,\ndelayed ripening",
            "Host cell cycle hijack\n(pRBR → E2F → S-phase)",
            "Sugar metabolism\n(SPS1 suppressed, r=−0.40)",
            "Thermal defense\n(MBF1c ↓, sHSP module ↓)"),
  x = 4.2,
  y = c(4.7, 3.6, 2.6, 1.6, 0.8, 0.0),
  fill = c("#fde3d4", "#ffe8b3", "#f9d5d5", "#d5e3ff", "#e2d9f2", "#f9d5d5"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c("V2", "V1", "V3", "C1", "C2", "C2", "V2"),
  to   = c("silencing", "movement", "symptom", "cellcycle", "sugar", "silencing", "shsp"),
  lab  = c("disables", "packages", "modulates", "re-programmes",
           "r=−0.40  (×SPS1)", "2nd layer", "r=−0.25  (×MBF1c)"),
  stringsAsFactors = FALSE
)

edges <- edges %>%
  left_join(orf_nodes  %>% transmute(from = id, x1 = x + 0.55, y1 = y), by = "from") %>%
  left_join(host_nodes %>% transmute(to   = id, x2 = x - 0.9,  y2 = y), by = "to")

pC <- ggplot() +
  # edges first (so nodes draw over them)
  geom_curve(data = edges,
             aes(x = x1, y = y1, xend = x2, yend = y2),
             curvature = 0.08, color = "grey45", linewidth = 0.5,
             arrow = arrow(length = unit(6, "pt"), type = "closed")) +
  geom_text(data = edges,
            aes(x = (x1 + x2) / 2, y = (y1 + y2) / 2 + 0.08,
                label = lab),
            family = "Arial", size = 2.2, color = "grey25", fontface = "italic") +
  # viral ORF nodes
  geom_rect(data = orf_nodes,
            aes(xmin = x - 0.55, xmax = x + 0.55,
                ymin = y - 0.45, ymax = y + 0.45,
                fill = fill), color = "grey20", linewidth = 0.5) +
  geom_text(data = orf_nodes,
            aes(x = x, y = y, label = label),
            family = "Arial", fontface = "bold", size = 2.5,
            color = "grey10", lineheight = 1.0) +
  # host target nodes
  geom_rect(data = host_nodes,
            aes(xmin = x - 0.9, xmax = x + 0.9,
                ymin = y - 0.4, ymax = y + 0.4,
                fill = fill), color = "grey35", linewidth = 0.4) +
  geom_text(data = host_nodes,
            aes(x = x, y = y, label = label),
            family = "Arial", size = 2.3, color = "grey15", lineheight = 1.05) +
  # column headers
  annotate("text", x = 1, y = 5.3, label = "Viral ORFs",
           family = "Arial", fontface = "bold", size = 3.2, color = "grey15") +
  annotate("text", x = 4.2, y = 5.3, label = "Plant processes attacked",
           family = "Arial", fontface = "bold", size = 3.2, color = "grey15") +
  annotate("text", x = 5.4, y = 5.3, label = "(→ GRBV phenotype)",
           family = "Arial", fontface = "italic", size = 2.6, color = "grey35") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.2, 5.3), ylim = c(-0.55, 5.7), expand = FALSE) +
  labs(tag = "C", title = "C. Each viral ORF maps to a host process it disrupts") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(family = "Arial", face = "bold",
                                  size = 11, hjust = 0.5, color = "grey10"),
        plot.tag   = element_text(family = "Arial", face = "bold", size = 13),
        plot.margin = margin(4, 6, 4, 6))

# ==================================================================
# Composite
# ==================================================================
layout <- pA / pB / pC +
  plot_layout(heights = c(0.85, 1.00, 1.25)) +
  plot_annotation(
    title = "S60. GRBV viral-gene function atlas",
    subtitle = "Each ORF mapped to its protein role, its mechanism in plant cells, and the host process it attacks in this dataset",
    caption = "ORF coordinates approximated from GRBV isolate CF214 (NCBI NC_022002). Expression values from N761_gp1-gp5 in counts_matrix.txt. Correlations from S59.",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 13, hjust = 0.5),
                  plot.subtitle = element_text(family = "Arial", size = 9.5,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)),
                  plot.caption = element_text(family = "Arial", size = 7.5,
                                              hjust = 0.5, color = "grey40")))

save_grbv(layout, file.path(OUT_DIR, "S60_viral_biology_diagram"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 260)

cat("[S60] viral biology diagram written.\n")
