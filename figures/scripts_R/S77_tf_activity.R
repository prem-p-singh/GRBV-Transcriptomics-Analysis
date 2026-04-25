# S77 — TF activity inference (regulon-based) vs expression-only TF DE.
# Per-sample TF activity scores derived from each TF's GENIE3 regulon
# (top 200 targets per TF, weighted-sum z-score on log2 CPM).  Welch-t
# RB(+) vs RB(-) on activity ranks TFs by regulatory disruption rather
# than expression change.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/17_tf_activity"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S77_tf_activity")

de <- read_csv(file.path(RESULTS, "tf_activity_de.csv"), show_col_types = FALSE) %>%
  filter(!is.na(p_value)) %>%
  mutate(minus_log10_padj = -log10(pmax(padj, 1e-300)),
         direction = ifelse(diff > 0, "Higher in RB+", "Lower in RB+"))

# ---- panel A: top 25 TFs by activity DE ------------------------------------
top <- de %>% arrange(p_value) %>% slice_head(n = 25) %>%
  mutate(tf = factor(tf, levels = rev(tf)))

pA <- ggplot(top, aes(x = diff, y = tf, fill = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("padj=%.1e", padj)),
            x = 0.02, hjust = 0, size = 2.4, color = "grey30") +
  scale_fill_manual(values = c("Higher in RB+" = "#B2182B",
                                "Lower in RB+" = "#2166AC"),
                    name = NULL) +
  labs(title = "(A) Top 25 TFs by activity disruption",
       subtitle = "Welch-t on regulon-derived activity scores; BH-adjusted",
       x = "Activity diff (RB+ - RB-)", y = NULL) +
  theme_grbv(base_size = 9) +
  theme(legend.position = "top")

# ---- panel B: volcano on all TFs -------------------------------------------
pB <- ggplot(de, aes(x = diff, y = minus_log10_padj, color = direction)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.65, size = 1.5) +
  scale_color_manual(values = c("Higher in RB+" = "#B2182B",
                                 "Lower in RB+" = "#2166AC"),
                     guide = "none") +
  labs(title = sprintf("(B) Volcano: all %d TFs tested", nrow(de)),
       subtitle = "Dashed grey: BH FDR = 0.05",
       x = "Activity diff (RB+ - RB-)", y = "-log10 padj") +
  theme_grbv(base_size = 9)

fig <- (pA + pB) + plot_layout(widths = c(1.4, 1)) +
  plot_annotation(title = "S77 — Regulon-based TF activity (RB+ vs RB-)",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 170, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 170, units = "mm")
cat("[S77] saved\n")
