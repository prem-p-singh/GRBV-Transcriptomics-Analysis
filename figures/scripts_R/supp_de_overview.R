# Supplementary batch: regenerate per-question DE overview figures (S2 - S7)
# using DESeq2 results from de_deseq2_voom and de_deseq2_extras. Each figure is a 4-panel
# composite: DEG counts bar, volcano plot, log2FC histogram, top-20 DEG table.
#
# Output: farm_pipeline/figures/supp/S{02..07}_*.{png,pdf}

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(ggrepel); library(scales); library(data.table)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results"
ANN_F   <- file.path(RESULTS, "annotations/gene_annotations.csv")
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ <- 0.05; LFC <- 1.0

# --------- load ---------
all_de <- rbindlist(lapply(c(
  list.files(file.path(RESULTS, "de_redo"),
             pattern = "_deseq2\\.csv$", full.names = TRUE),
  list.files(file.path(RESULTS, "de_extras"),
             pattern = "_deseq2\\.csv$", full.names = TRUE)
), function(f) { d <- fread(f); d$file <- basename(f); d }),
fill = TRUE) %>% as.data.frame()

ann <- fread(ANN_F) %>% as.data.frame()
name_map <- setNames(ann$product, ann$gene_id)

# --------- helper: one DE-overview figure for a set of comparisons ---
#
# Inputs
#   comparisons : named list of comparison -> display label
#   primary     : the "primary" comparison to use for volcano + top-20 panel
#   title       : top title for the figure
#
render_overview <- function(comparisons, primary, title, out_stem) {

  # Panel A — DEG counts per sub-comparison
  cnt <- do.call(rbind, lapply(names(comparisons), function(nm) {
    d <- all_de %>% filter(comparison == nm)
    if (!nrow(d)) return(NULL)
    ok <- !is.na(d$padj) & d$padj < PADJ & abs(d$log2FC) >= LFC
    data.frame(label = comparisons[[nm]],
               up = sum(ok & d$log2FC >  0),
               down = sum(ok & d$log2FC <  0))
  })) %>%
    mutate(label = factor(label, levels = rev(comparisons)))
  cnt_long <- cnt %>% pivot_longer(c(up, down), names_to = "dir", values_to = "n") %>%
    mutate(n_signed = ifelse(dir == "up", n, -n))
  totals <- cnt %>% mutate(total = up + down) %>% select(label, total)

  # Panel A — horizontal DEG counts with totals placed at the longer-side
  # bar tip (Main01-style) instead of on the zero line.
  tot_labels <- cnt %>%
    mutate(total = up + down,
           # place the "n = X" label on the side with more DEGs
           x_end = ifelse(up >= down,  up, -down),
           # modest outside offset, sign depends on side
           x_lab = x_end + ifelse(up >= down,  0.03,  -0.03) * max(cnt$up + cnt$down),
           hjust = ifelse(up >= down, 0, 1))

  pA <- ggplot(cnt_long, aes(x = n_signed, y = label, fill = dir)) +
    geom_col(width = 0.70, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey25", linewidth = 0.4) +
    geom_text(data = tot_labels,
              aes(x = x_lab, y = label, label = sprintf("n = %s", scales::comma(total)),
                  hjust = hjust),
              inherit.aes = FALSE, family = "Arial", size = 2.5,
              fontface = "bold", color = "grey15") +
    scale_x_continuous(labels = function(x) scales::comma(abs(x)),
                       name = "DEGs (DOWN \u2190 0 \u2192 UP)",
                       expand = expansion(mult = c(0.22, 0.22))) +
    scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                                 down = unname(GRBV_COLORS["down"])),
                      labels = c(up = "Up", down = "Down"), name = NULL) +
    labs(y = NULL, tag = "A", title = "DEG counts") +
    theme_grbv() +
    theme(legend.position = "top",
          legend.margin = margin(b = -3),
          plot.title = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(face = "bold", size = 8))

  # Panel B — volcano for the primary comparison (clean labels: top 6 only)
  prim <- all_de %>% filter(comparison == primary) %>%
    filter(!is.na(padj), !is.na(log2FC))
  prim$status <- with(prim,
                      ifelse(padj < PADJ & log2FC >=  LFC, "Up",
                      ifelse(padj < PADJ & log2FC <= -LFC, "Down", "ns")))
  prim$negLog <- -log10(pmax(prim$padj, 1e-300))
  top_lab <- prim %>% filter(status != "ns") %>%
    arrange(padj) %>% slice_head(n = 6)

  pB <- ggplot(prim, aes(x = log2FC, y = negLog, color = status)) +
    # significance bands as soft background shading
    geom_rect(data = data.frame(x = 1), inherit.aes = FALSE,
              aes(xmin = -LFC, xmax = LFC, ymin = -Inf, ymax = Inf),
              fill = "grey95", alpha = 0.6) +
    geom_hline(yintercept = -log10(PADJ), linetype = "dashed",
               color = "grey55", linewidth = 0.35) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed",
               color = "grey55", linewidth = 0.35) +
    geom_point(alpha = 0.6, size = 0.75) +
    scale_color_manual(values = c(Up   = unname(GRBV_COLORS["up"]),
                                  Down = unname(GRBV_COLORS["down"]),
                                  ns   = unname(GRBV_COLORS["neutral"])),
                       name = NULL,
                       labels = function(x) {
                         n <- sapply(x, function(s) sum(prim$status == s))
                         sprintf("%s (%s)", x, scales::comma(n))
                       }) +
    geom_text_repel(data = top_lab, aes(label = gene),
                    family = "Arial", size = 2.4, fontface = "bold",
                    color = "black", max.overlaps = Inf,
                    box.padding = 0.45, point.padding = 0.15,
                    segment.size = 0.2, segment.color = "grey45",
                    min.segment.length = 0) +
    labs(x = expression(log[2]*" fold change"),
         y = expression(-log[10]*"(padj)"),
         tag = "B", title = sprintf("Volcano: %s", comparisons[[primary]])) +
    theme_grbv() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          legend.position = c(0.02, 0.97),
          legend.justification = c(0, 1),
          legend.text = element_text(size = 7),
          legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
          legend.key.height = unit(8, "pt"))

  # Panel C — log2FC histogram + density overlay (Main01-panel-C polish)
  sig <- prim %>% filter(status != "ns")
  pC <- ggplot(sig, aes(x = log2FC)) +
    geom_histogram(aes(fill = status, y = after_stat(count)),
                   binwidth = 0.25, color = "white", alpha = 0.88) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = c(-LFC, LFC), color = "grey60",
               linetype = "dashed", linewidth = 0.3) +
    scale_fill_manual(values = c(Up = unname(GRBV_COLORS["up"]),
                                 Down = unname(GRBV_COLORS["down"])),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = expression(log[2]*" fold change"),
         y = "Significant genes", tag = "C",
         title = "FC distribution") +
    theme_grbv() +
    theme(plot.title = element_text(size = 10, face = "bold"))

  # Panel D — top 20 DEGs table
  top20 <- prim %>% filter(status != "ns") %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(padj) %>% slice_head(n = 20) %>%
    mutate(product = ifelse(gene %in% names(name_map),
                            substr(name_map[gene], 1, 26), ""),
           label = paste0(gene, " \u2014 ", product),
           label = factor(label, levels = rev(label)))

  # Always place the p-value OUTSIDE the bar on the right, regardless of sign,
  # so text never collides with the y-axis label gutter.
  pval_x <- pmax(top20$log2FC, 0) + 0.35

  pD <- ggplot(top20, aes(x = log2FC, y = label, fill = log2FC)) +
    geom_col(alpha = 0.92, width = 0.72, color = "white") +
    geom_text(aes(x = pval_x, label = sprintf("p=%.1e", padj)),
              hjust = 0, family = "Arial", size = 2.1, color = "grey20") +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
    scale_fill_gradient2(low = unname(GRBV_COLORS["down"]), mid = "white",
                         high = unname(GRBV_COLORS["up"]), midpoint = 0,
                         guide = "none", oob = scales::squish,
                         limits = c(-6, 6)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.40))) +
    labs(x = "log2 fold change", y = NULL, tag = "D",
         title = "Top-20 DEGs by padj") +
    theme_grbv() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 6.5))

  # Top row: equal widths. Bottom row: give panel D ~60% of horizontal space
  # so long gene-product labels + right-side p-values fit.
  top_row    <- (pA | pB) + plot_layout(widths = c(1, 1))
  bottom_row <- (pC | pD) + plot_layout(widths = c(0.75, 1.25))
  layout <- top_row / bottom_row +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(
      title = title,
      subtitle = "DESeq2 (Wald), padj < 0.05, |log2FC| \u2265 1",
      theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                              size = 12, hjust = 0.5,
                                              margin = margin(b = 2)),
                    plot.subtitle = element_text(family = "Arial", size = 8.5,
                                                 hjust = 0.5, color = "grey25",
                                                 margin = margin(b = 6)))
    )

  save_grbv(layout, file.path(OUT_DIR, out_stem),
            width_mm = GRBV$W_DOUBLE_MM, height_mm = 200)
}

# ------------ S2 Q1 infection years ----------
render_overview(
  comparisons = list(
    Q1a_recent_le2yr_vs_neg = "Recent (\u22642 yr)",
    Q1a_mid_3to4yr_vs_neg   = "Mid (3-4 yr)",
    Q1a_old_ge5yr_vs_neg    = "Old (\u22655 yr)",
    Q1b_cohort2016_vs_neg   = "Cohort 2016",
    Q1b_cohort2019_vs_neg   = "Cohort 2019",
    Q1b_cohort2020_vs_neg   = "Cohort 2020",
    Q1b_cohort2021_vs_neg   = "Cohort 2021",
    Q1b_cohort2022_vs_neg   = "Cohort 2022"
  ),
  primary = "Q1a_old_ge5yr_vs_neg",
  title = "S2. Q1 \u2014 Infection duration (years since infection)",
  out_stem = "S02_Q1_infection_years"
)

# ------------ S3 Q2 maturity stages ----------
render_overview(
  comparisons = list(
    Q2_Pre_veraison  = "Pre-v\u00e9raison",
    Q2_Veraison      = "V\u00e9raison",
    Q2_Post_veraison = "Post-v\u00e9raison",
    Q2_Harvest       = "Harvest"
  ),
  primary = "Q2_Harvest",
  title = "S3. Q2 \u2014 Maturity stages (RB+ vs RB\u2212)",
  out_stem = "S03_Q2_maturity_stages"
)

# ------------ S4 Q3 seasonal effects in RB+ ----------
render_overview(
  comparisons = list(Q3_seasonRBpos = "2022 vs 2021 (RB+)"),
  primary = "Q3_seasonRBpos",
  title = "S4. Q3 \u2014 Seasonal effects in infected vines",
  out_stem = "S04_Q3_seasonal_effects"
)

# ------------ S5 Q4 environmental modulation ----------
render_overview(
  comparisons = list(
    Q4_seasonXinfection      = "Q4 LRT (Sea \u00d7 Inf)",
    Q4d_seasonXinfection_dir = "Q4 interaction direction"
  ),
  primary = "Q4d_seasonXinfection_dir",
  title = "S5. Q4 \u2014 Environmental modulation of GRBV effects",
  out_stem = "S05_Q4_environmental_effects"
)

# ------------ S6 Q5 baseline seasonal change in RB- ----------
render_overview(
  comparisons = list(Q5_seasonRBneg = "2022 vs 2021 (RB\u2212)"),
  primary = "Q5_seasonRBneg",
  title = "S6. Q5 \u2014 Baseline seasonal change in uninfected vines",
  out_stem = "S06_Q5_uninfected_environment"
)

# ------------ S7 Q6 pooled transcriptome trajectory ----------
render_overview(
  comparisons = list(Q6_pooled = "RB+ vs RB\u2212 (pooled)"),
  primary = "Q6_pooled",
  title = "S7. Q6 \u2014 Pooled transcriptome response to GRBV",
  out_stem = "S07_Q6_transcriptome_trajectory"
)

cat("[supp S2-S7] done.\n")
