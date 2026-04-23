# Main Composite 10 — GRBV molecular cascade summary diagram.
#
# A single-panel publication cartoon showing the disease cascade from viral
# infection through host transcriptome reprogramming to phenotype. Each
# stage of the cascade is annotated with the corresponding evidence figure.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(grid)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "Main10_biology_cascade")

# Five cascade stages, evenly spaced with taller boxes so 4-line bodies fit.
# Box fills follow a viridis 5-step ordinal (early -> late in cascade),
# matching the manuscript's stage palette; hue conveys position, not meaning.
# Light text is used on dark stations (1-3), dark text on light stations (4-5).
box_h <- 1.45
y_centers <- seq(10.5 - box_h/2, 1.0 + box_h/2, length.out = 5)
station_fill <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
station_text <- c("white",   "white",   "white",   "grey10",  "grey10")
nodes <- data.frame(
  y = y_centers,
  fill = station_fill,
  text_col = station_text,
  title = c("1. GRBV infection establishes",
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
  )
)

box <- data.frame(
  xmin = 0.4, xmax = 9.6,
  ymin = nodes$y - box_h/2, ymax = nodes$y + box_h/2,
  fill = nodes$fill
)

# Arrows connecting boxes
arrows <- data.frame(
  x = 5, xend = 5,
  y = head(nodes$y, -1) - box_h/2,
  yend = tail(nodes$y, -1) + box_h/2
)

p <- ggplot() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 11), expand = FALSE) +
  geom_rect(data = box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = box$fill, color = "black", linewidth = 0.5, alpha = 0.95) +
  geom_text(data = nodes,
            aes(x = 0.7, y = y + box_h/2 - 0.22, label = title),
            family = "Arial", color = nodes$text_col, fontface = "bold",
            hjust = 0, vjust = 1, size = 3.6) +
  geom_text(data = nodes,
            aes(x = 0.7, y = y + box_h/2 - 0.52, label = body),
            family = "Arial", color = nodes$text_col, hjust = 0,
            vjust = 1, size = 2.7, lineheight = 1.08) +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.7, color = "grey25",
               arrow = arrow(length = unit(7, "pt"), type = "closed")) +
  labs(title = "GRBV molecular disruption cascade",
       subtitle = "Viral establishment \u2192 silent immunity \u2192 thermal-defence dampening\n\u2192 sugar-pathway rewiring \u2192 disease phenotype",
       caption = "Robust temp \u2192 titer mediators (multi-pipeline): MBF1c (LOC100249249), SPS1 (LOC100232974), uncharacterised LOC100266309") +
  theme_void(base_family = "Arial") +
  theme(plot.title    = element_text(family = "Arial", face = "bold",
                                     size = 14, hjust = 0.5,
                                     margin = margin(t = 4, b = 2)),
        plot.subtitle = element_text(family = "Arial", size = 9,
                                     hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        plot.caption  = element_text(family = "Arial", size = 8,
                                     hjust = 0.5, color = "grey35",
                                     margin = margin(t = 6)),
        plot.margin   = margin(8, 8, 8, 8))

save_grbv(p, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 220)

cat("[main10] cascade diagram complete\n")
