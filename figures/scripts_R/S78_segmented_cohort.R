# S78 — Segmented regression on cohort age.
# Tests whether the deck's "three temporal phases" (acute Yr 1 / peak Yrs 3-5 /
# chronic Yrs 6-7) reflect continuous gradient vs phase transitions.  Fits
# segmented regression with k=0/1/2 break-points on PC1, PC2, and 13 module
# eigengenes against years-since-infection on n=264 RB(+) samples.
#
# Result: best_k=0 (linear) for ALL 15 responses by BIC.  The phase framing
# is descriptive simplification, not statistically supported phase transitions.
# Notable null finding.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/18_segmented"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S78_segmented_cohort")

s <- read_csv(file.path(RESULTS, "segmented_summary.csv"), show_col_types = FALSE) %>%
  rowwise() %>%
  mutate(BIC_min = min(c_across(starts_with("BIC_")), na.rm = TRUE),
         dBIC_k1 = BIC_k1 - BIC_k0,
         dBIC_k2 = ifelse(is.na(BIC_k2), NA_real_, BIC_k2 - BIC_k0)) %>%
  ungroup()

# ---- panel A: dBIC for k=1 vs k=0 (positive = k=0 wins) --------------------
sA <- s %>%
  select(response, dBIC_k1) %>%
  arrange(desc(dBIC_k1)) %>%
  mutate(response = factor(response, levels = response))

pA <- ggplot(sA, aes(x = dBIC_k1, y = response,
                      fill = ifelse(dBIC_k1 > 0, "linear preferred", "1-break preferred"))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("linear preferred" = "#2CA02C",
                                "1-break preferred" = "#B2182B"),
                    name = NULL) +
  labs(title = "(A) Delta BIC: 1-break minus linear",
       subtitle = "Positive = linear (k=0) preferred; n=264 RB+ samples",
       x = "BIC(k=1) - BIC(k=0)", y = NULL) +
  theme_grbv(base_size = 9) +
  theme(legend.position = "top")

# ---- panel B: best_k counts -------------------------------------------------
sB <- s %>% count(best_k, name = "n_responses") %>% mutate(best_k = factor(best_k))

pB <- ggplot(sB, aes(x = best_k, y = n_responses, fill = best_k)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n_responses), vjust = -0.4, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("0" = "#2CA02C", "1" = "#FFA500", "2" = "#B2182B")) +
  scale_y_continuous(limits = c(0, max(sB$n_responses) * 1.2),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "(B) Best-fitting break-point count by BIC",
       subtitle = "PC1, PC2, + 13 module eigengenes",
       x = "Number of break-points (k)", y = "# responses") +
  theme_grbv(base_size = 9)

fig <- (pA + pB) + plot_layout(widths = c(1.4, 1)) +
  plot_annotation(title = "S78 — Segmented regression on cohort age (null result: continuous gradient)",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 150, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 150, units = "mm")
cat("[S78] saved\n")
