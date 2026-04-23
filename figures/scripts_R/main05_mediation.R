# Main Composite 5 — Causal mediation: Temperature -> Viral titer -> Gene
# (R/ggplot2 version, replaces main05_mediation.py).
#
# Panels:
#   A  Headline mediation path diagram for MBF1c (LOC100249249)
#   B  Direct vs indirect effect scatter for 2,000 variable genes
#   C  Top-10 mediators ranked by joint significance
#   D  Category counts
#
# Inputs:
#   /Users/prempratapsingh/Data_Analysis/mediation_results.csv  (fix03 tv=2000)
#   /Users/prempratapsingh/Data_Analysis/farm_pipeline/results/fix11_annotations/gene_annotations.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
  library(patchwork); library(ggrepel); library(scales)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
MED <- file.path(DATA_ROOT, "mediation_results.csv")
ANN <- file.path(DATA_ROOT, "farm_pipeline/results/fix11_annotations/gene_annotations.csv")
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "Main05_mediation")

med <- read_csv(MED, show_col_types = FALSE)
ann <- read_csv(ANN, show_col_types = FALSE) %>%
  select(gene_id, product) %>% distinct()
name_map <- setNames(ann$product, ann$gene_id)

short_label <- c(
  LOC100249249 = "MBF1c",
  LOC100853778 = "PSK6",
  LOC100261367 = "PP2-B15",
  LOC100253072 = "DREB1E",
  LOC100247528 = "GRX-C9",
  LOC100232974 = "SPS1"
)

med <- med %>%
  mutate(
    category = case_when(
      sig_joint                   ~ "Titer-mediated",
      total_padj < 0.05           ~ "Temperature-only",
      TRUE                        ~ "Not significant"
    ),
    category = factor(category,
                      levels = c("Titer-mediated","Temperature-only","Not significant"))
  )

n_total <- nrow(med)
sig <- med %>% filter(sig_joint) %>% arrange(path_b_padj)
n_sig       <- nrow(sig)
n_temp_only <- sum(med$category == "Temperature-only")
n_neither   <- sum(med$category == "Not significant")

cat_color <- pal_mediation

# -----------------------------------------------------------------
# Panel A — MBF1c headline path diagram
# -----------------------------------------------------------------
mbf <- sig[1, ]

# Coordinates for boxes and arrows (data-coords on a 0-10 x 0-10 canvas)
nodes <- data.frame(
  x = c(1.9, 5.0, 8.1),
  y = c(4.75, 8.15, 4.75),
  w = c(3.4, 3.2, 3.4),
  h = c(1.6, 1.6, 1.6),
  # Temperature = amber; Viral titer = blue (mediator); MBF1c = purple (outcome)
  fill = c(unname(GRBV_COLORS["med_temp"]),
           unname(GRBV_COLORS["med_titer"]),
           unname(GRBV_COLORS["highlight"])),
  label_top = c("Temperature", "Viral titer", "MBF1c"),
  label_sub = c("",            "",            "(LOC100249249)")
)

arrows <- data.frame(
  x    = c(3.45, 6.55, 3.6),
  y    = c(5.55, 7.75, 4.75),
  xend = c(3.95, 6.05, 6.4),
  yend = c(7.75, 5.55, 4.75),
  type = c("a", "b", "cprime")
)
arrow_labels <- data.frame(
  x = c(3.20, 6.80),
  y = c(6.7, 6.7),
  label = c("a", "b"),
  size  = c(4.2, 4.2),
  color = c("black", "black")
)

stats_text <- sprintf(
  paste(
    "Total c = %+.3f   (p = %s)",
    "Direct c' = %+.3f",
    "Indirect a \u00d7 b = %+.3f",
    "joint padj = %s",
    "%% mediated = %+.1f%%",
    sep = "\n"),
  mbf$total_c, scales::scientific(mbf$total_p, digits = 2),
  mbf$direct_cprime,
  mbf$indirect_ab,
  scales::scientific(mbf$path_b_padj, digits = 2),
  mbf$prop_mediated * 100
)

pA <- ggplot() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10), expand = FALSE) +
  geom_rect(data = nodes,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2),
            fill = nodes$fill, color = "black", linewidth = 0.4) +
  geom_text(data = nodes,
            aes(x = x, y = y + 0.18, label = label_top),
            family = "Arial", color = "white", fontface = "bold", size = 3.4) +
  geom_text(data = nodes,
            aes(x = x, y = y - 0.34, label = label_sub),
            family = "Arial", color = "white", size = 2.1) +
  geom_segment(data = arrows %>% filter(type %in% c("a","b")),
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "black", linewidth = 0.55,
               arrow = arrow(length = unit(7, "pt"), type = "closed")) +
  geom_segment(data = arrows %>% filter(type == "cprime"),
               aes(x = x, y = y, xend = xend, yend = yend),
               color = unname(GRBV_COLORS["RB_pos"]),
               linewidth = 0.6, linetype = "dashed",
               arrow = arrow(length = unit(7, "pt"), type = "closed")) +
  geom_text(data = arrow_labels,
            aes(x = x, y = y, label = label),
            family = "Arial", fontface = "bold",
            color = arrow_labels$color, size = arrow_labels$size) +
  annotate("text", x = 5.0, y = 2.15, label = stats_text,
           family = "Arial", size = 2.7, hjust = 0.5, vjust = 0.5,
           lineheight = 1.15) +
  annotate("text", x = 5.0, y = 0.55,
           label = "multiprotein-bridging factor 1c",
           family = "Arial", size = 2.7, fontface = "italic", color = "grey25") +
  labs(tag = "A") +
  theme_void(base_family = "Arial") +
  theme(plot.tag = element_text(face = "bold", size = 13),
        plot.tag.position = c(0.01, 0.98),
        plot.margin = margin(2, 2, 2, 2))

# -----------------------------------------------------------------
# Panel B — direct vs indirect scatter
# -----------------------------------------------------------------
top_label <- sig %>% slice_head(n = 5) %>%
  mutate(label = ifelse(gene %in% names(short_label),
                        short_label[gene], gene))

pB <- ggplot(med, aes(x = direct_cprime, y = indirect_ab, color = category)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.55, size = 0.9) +
  geom_text_repel(data = top_label,
                  aes(label = label), color = "black",
                  family = "Arial", fontface = "bold", size = 2.6,
                  box.padding = 0.4, max.overlaps = Inf,
                  segment.size = 0.25, segment.color = "grey40") +
  scale_color_manual(values = cat_color,
                     labels = function(x) {
                       n <- sapply(x, function(c) sum(med$category == c))
                       sprintf("%s (%d)", x, n)
                     }) +
  labs(x = "Direct effect c' (Temp \u2192 Gene | Titer)",
       y = "Indirect effect a \u00d7 b",
       color = NULL, tag = "B") +
  theme_grbv() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7.5),
        legend.key.height = unit(7, "pt"),
        legend.margin = margin(t = 0, b = 0))

# -----------------------------------------------------------------
# Panel C — top 10 mediators ranked by joint padj
# -----------------------------------------------------------------
top10 <- sig %>% slice_head(n = 10) %>%
  mutate(short = ifelse(gene %in% names(short_label),
                        short_label[gene], gene),
         product = name_map[gene],
         product = ifelse(is.na(product), "", substr(product, 1, 38)),
         label = paste0(short, "\n", product),
         label = factor(label, levels = rev(label)),
         neg_log = -log10(pmax(path_b_padj, 1e-300)))

pC <- ggplot(top10, aes(x = neg_log, y = label)) +
  geom_col(fill = unname(GRBV_COLORS["med_titer"]), alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", neg_log)),
            hjust = -0.15, family = "Arial", size = 2.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression(-log[10]*"(joint padj)"), y = NULL, tag = "C") +
  theme_grbv() +
  theme(axis.text.y = element_text(size = 7.0, lineheight = 0.9))

# -----------------------------------------------------------------
# Panel D — category bar chart
# -----------------------------------------------------------------
counts <- med %>% count(category) %>%
  mutate(category = factor(category,
                           levels = c("Titer-mediated","Temperature-only","Not significant")))
pD <- ggplot(counts, aes(x = n, y = category, fill = category)) +
  geom_col(alpha = 0.92, width = 0.7, color = "white") +
  geom_text(aes(label = scales::comma(n)), hjust = -0.15,
            family = "Arial", size = 2.8) +
  scale_fill_manual(values = cat_color, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20)),
                     labels = scales::label_comma()) +
  scale_y_discrete(limits = rev) +
  labs(x = "Number of genes", y = NULL, tag = "D",
       caption = "BH FDR<0.05 on both Temp\u2192Gene and Titer\u2192Gene paths;  |a\u00d7b|>0.05;  top 2,000 most variable genes") +
  theme_grbv() +
  theme(plot.caption = element_text(size = 6.5, hjust = 1, color = "grey35"))

# -----------------------------------------------------------------
# Composite layout
# -----------------------------------------------------------------
layout <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1, 1), heights = c(1.0, 1.0)) +
  plot_annotation(
    title = "Causal mediation: temperature \u2192 viral titer \u2192 gene expression",
    theme = theme(plot.title = element_text(family = "Arial",
                                            face = "bold", size = 13,
                                            hjust = 0.5,
                                            margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 175)

cat(sprintf("\n[main05] joint-significant: %d/%d; temp-only: %d; none: %d\n",
            n_sig, n_total, n_temp_only, n_neither))
