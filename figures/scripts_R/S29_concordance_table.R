# S29 v2 — Method-concordance table (replaces truncated Headline-findings panel).
#
# Redesign brief (SciAgent-Skills review):
#   * Original S29 panel B was free-form text on a plot and got truncated
#     ("MBF1c (LOC100249249) = cross-pipeline medi...").
#   * v2 replaces panel B with a real method-concordance table showing
#     Welch-t (paper) vs DESeq2 vs voom vs RM-voom DEG counts per contrast.
#     Numbers sourced from docs/DE_REANALYSIS_REPORT.md and CLAUDE.md (fix 1b).
#   * Panel A (core-count bars) kept but shrunk to 35 % width.
#
# Output: S29_comprehensive_summary_dashboard.png + .svg in figures/supp/.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(scales); library(data.table)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------- Panel A: core counts bar (trimmed from original) --------------
sum_df <- data.frame(
  metric = c("Samples", "Genes tested",
             "Q6 DEGs (DESeq2)", "Q2 DEGs (union)",
             "Q3 season DEGs", "Q4 interaction DEGs",
             "Biomarker panel", "WGCNA modules"),
  value  = c(327, 14604,
             126, 1329,
             121, 52,
             45, 13)
) %>% mutate(metric = factor(metric, levels = rev(metric)))

pS29a <- ggplot(sum_df, aes(x = value, y = metric)) +
  geom_col(aes(fill = metric), alpha = 0.9, width = 0.72, color = "white") +
  geom_text(aes(label = scales::comma(value)),
            hjust = -0.1, family = "Arial", size = 2.8, color = "grey15") +
  scale_fill_manual(values = rep(c(unname(GRBV_COLORS["RB_pos"]),
                                    unname(GRBV_COLORS["season_2021"]),
                                    unname(GRBV_COLORS["warm"])),
                                  length.out = nrow(sum_df)),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35)),
                     labels = scales::label_comma()) +
  labs(x = "Count", y = NULL, tag = "A",
       title = "Core counts") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 8),
        plot.margin = margin(6, 10, 6, 6))

# -------------- Panel B: DE method-concordance table --------------
# Numbers from docs/DE_REANALYSIS_REPORT.md (Welch / DESeq2 / voom) and
# CLAUDE.md fix-1b (RM-voom).  Empty cells shown as "-" (not tested).
concord <- tribble(
  ~contrast,                 ~welch,    ~deseq2, ~voom, ~rmvoom, ~jaccard,
  "Q6 pooled RB+ vs RB-",    "~908",    "126",   "150", "150",   0.62,
  "Q2 Pre-veraison",         "-",       "111",   "105", "105",    NA,
  "Q2 Veraison",             "-",       "130",   "124", "124",    NA,
  "Q2 Post-veraison",        "-",       "541",   "532", "529",    NA,
  "Q2 Harvest",              "-",       "1,025", "964", "964",    NA,
  "Q3 season (RB+)",         "-",       "121",   "153", "-",      NA,
  "Q4 season x infection",   "~773",    "52",    "-",   "-",      NA,
  "Q5 season (RB-)",         "-",       "182",   "228", "-",      NA
)

# Wide -> long for a heatmap-style rendering
long <- concord %>%
  pivot_longer(cols = c(welch, deseq2, voom, rmvoom),
               names_to = "method", values_to = "count") %>%
  mutate(
    method = factor(method,
                    levels = c("welch", "deseq2", "voom", "rmvoom"),
                    labels = c("Welch-t\n(paper)", "DESeq2", "voom", "RM-voom")),
    contrast = factor(contrast, levels = rev(concord$contrast)),
    count_num = suppressWarnings(as.numeric(gsub("[,~-]", "", count))),
    fill_alpha = ifelse(is.na(count_num), 0, pmin(1, count_num / 1100))
  )

# Fill: grey for "-", warm for numeric (scaled by value)
pS29b <- ggplot(long, aes(x = method, y = contrast)) +
  geom_tile(aes(fill = count_num), color = "white", linewidth = 0.4) +
  geom_text(aes(label = count),
            family = "Arial", size = 3.0, color = "grey10") +
  scale_fill_gradient(low = "#F7EDE1", high = unname(GRBV_COLORS["warm"]),
                      na.value = "grey92", guide = "none") +
  # Jaccard shown as a footnote-style annotation inside the Q6 row instead of
  # a sixth column (kept the figure width tight)
  annotate("text", x = 2, y = 8, label = "Jaccard(DESeq2, voom) = 0.62",
           family = "Arial", size = 2.4, colour = "grey30",
           fontface = "italic", hjust = 0.5, vjust = -1.6) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(0.5, 8.5), clip = "off") +
  labs(x = NULL, y = NULL, tag = "B",
       title = "DE method concordance",
       subtitle = "padj < 0.05 and |log\u2082FC| \u2265 1 throughout") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        axis.text.x   = element_text(size = 8, hjust = 0.5),
        axis.text.y   = element_text(size = 8, hjust = 1),
        axis.ticks    = element_blank(),
        axis.line     = element_blank(),
        panel.grid    = element_blank(),
        plot.margin   = margin(6, 8, 6, 6))

# -------------- Compose --------------
p29 <- pS29a + pS29b + plot_layout(widths = c(0.38, 0.62))
p29 <- p29 + plot_annotation(
  title = "S29. Integrated dashboard: core counts and DE-method concordance",
  subtitle = "Paper's Welch-t counts shown next to covariate-adjusted DESeq2 / voom / repeated-measures-voom re-analysis (fix #1 / #1b)",
  theme = theme(plot.title    = element_text(family = "Arial", face = "bold",
                                             size = 13, hjust = 0.5,
                                             margin = margin(b = 2)),
                plot.subtitle = element_text(family = "Arial", size = 9,
                                             hjust = 0.5, colour = "grey25",
                                             margin = margin(b = 6))))

save_grbv(p29, file.path(OUT_DIR, "S29_comprehensive_summary_dashboard"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 130)

cat("[S29_v2] concordance dashboard complete\n")
