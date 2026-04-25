# S73 — Viral ORF extensions: phenotype prediction + module correlation.
#
# Backs slides 49-54 with a cleaner publication panel.  Two parts:
#   A: per-ORF effect on Brix/TA/Malic AFTER controlling for qPCR titer
#      (LRT p_orf_added).  Asks "does the viral mRNA add information beyond
#      what titer already captures?"
#   B: per-ORF x 13 WGCNA module Spearman correlation (RB+ only) — which
#      host module each viral gene tracks.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/13_viral_orf_ext"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S73_viral_orf_extensions")

# ---- panel A: phenotype prediction LRT --------------------------------------
phen <- read_csv(file.path(RESULTS, "orf_phenotype_prediction.csv"), show_col_types = FALSE)

phen <- phen %>%
  filter(!is.na(LRT_p_orf_added)) %>%
  mutate(
    minus_log10p = -log10(pmax(LRT_p_orf_added, 1e-30)),
    sig_label = case_when(
      LRT_p_orf_added < 0.001 ~ "***",
      LRT_p_orf_added < 0.01  ~ "**",
      LRT_p_orf_added < 0.05  ~ "*",
      TRUE ~ ""),
    delta_R2 = R2_titer_plus_orf - R2_titer_only,
    ORF = factor(ORF, levels = c("N761_gp1","N761_gp2","N761_gp3",
                                 "N761_gp4","N761_gp5","N761_gp6"))
  )

pA <- ggplot(phen, aes(x = ORF, y = phenotype, fill = delta_R2)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), vjust = 0.7, size = 3.6, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.05, 0.20),
                       oob = squish, name = expression(Delta*R^2)) +
  labs(title = "(A) ORF beyond qPCR titer (LRT)",
       subtitle = "Stars = LRT p_orf_added (* < 0.05, ** < 0.01, *** < 0.001)",
       x = NULL, y = NULL) +
  theme_grbv(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---- panel B: ORF x module heatmap -----------------------------------------
mod <- read_csv(file.path(RESULTS, "orf_module_correlation.csv"), show_col_types = FALSE)
mod_wide <- mod %>%
  mutate(
    sig_label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE ~ ""),
    ORF = factor(ORF, levels = c("N761_gp1","N761_gp2","N761_gp3",
                                 "N761_gp4","N761_gp5","N761_gp6"))
  )

pB <- ggplot(mod_wide, aes(x = ORF, y = module, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_label), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-0.6, 0.6), oob = squish,
                       name = expression(rho)) +
  labs(title = "(B) Per-ORF x WGCNA module correlation (RB+ only)",
       subtitle = "Spearman; n=264; padj BH per ORF",
       x = NULL, y = NULL) +
  theme_grbv(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---- compose -----------------------------------------------------------------
fig <- (pA / pB) + plot_layout(heights = c(1, 1.5)) +
  plot_annotation(title = "S73 — Viral ORF extensions: phenotype + module dissection",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 200, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 200, units = "mm")
cat("[S73] saved\n")
