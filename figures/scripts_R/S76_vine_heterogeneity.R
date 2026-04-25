# S76 — Vine-level heterogeneity within genotype-matched cohort cells.
# Within a single Merlot vineyard, do all 44 vines respond uniformly to GRBV?
# Per-vine response score = mean |log2FC| from cohort-matched RB(-) centroid
# on the top 2,000 variable genes; KMeans 3-cluster splits vines into
# tolerant / intermediate / susceptible.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/16_vine_heterogeneity"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S76_vine_heterogeneity")

cluster_pal <- c(tolerant = "#2CA02C",
                 intermediate = "#999999",
                 susceptible = "#B2182B")
cluster_levels <- c("tolerant", "intermediate", "susceptible")

# ---- panel A: per-vine response distribution -------------------------------
vines <- read_csv(file.path(RESULTS, "per_vine_response_scores.csv"), show_col_types = FALSE) %>%
  mutate(cluster = factor(cluster, levels = cluster_levels))

pA <- ggplot(vines, aes(x = mean_response, fill = cluster)) +
  geom_histogram(bins = 14, color = "white", linewidth = 0.4) +
  scale_fill_manual(values = cluster_pal, name = NULL) +
  labs(title = sprintf("(A) Per-vine infection response (n=%d vines)", nrow(vines)),
       subtitle = "Mean |log2FC vs cohort RB- centroid| on top 2000 variable genes",
       x = "Mean |log2FC|", y = "# vines") +
  theme_grbv(base_size = 9) +
  theme(legend.position = "top")

# ---- panel B: tolerant vs susceptible Brix ---------------------------------
brix <- read_csv(file.path(RESULTS, "tolerant_vs_susceptible_traits.csv"),
                 show_col_types = FALSE) %>%
  mutate(cluster = factor(cluster, levels = cluster_levels))

pB <- ggplot(brix, aes(x = cluster, y = mean, fill = cluster)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean - std/sqrt(count), ymax = mean + std/sqrt(count)),
                width = 0.18, linewidth = 0.6) +
  geom_text(aes(label = sprintf("n=%d\n%.2f Bx", count, mean)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = cluster_pal, guide = "none") +
  scale_y_continuous(limits = c(0, 27), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "(B) Harvest Brix by tolerance class",
       subtitle = sprintf("Tolerant - Susceptible delta = %.2f Bx",
                          brix$mean[brix$cluster=="tolerant"] -
                          brix$mean[brix$cluster=="susceptible"]),
       x = NULL, y = "Mean Brix +/- SE") +
  theme_grbv(base_size = 9)

# ---- panel C: cluster-stage matrix ----------------------------------------
vines_n <- vines %>%
  select(cluster, n_samples) %>%
  group_by(cluster) %>%
  summarise(median_n_samples = median(n_samples),
            n_vines = n(),
            .groups = "drop")

pC <- ggplot(vines_n, aes(x = cluster, y = n_vines, fill = cluster)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n_vines), vjust = -0.4, size = 3.6, fontface = "bold") +
  scale_fill_manual(values = cluster_pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "(C) Cluster sizes",
       x = NULL, y = "# vines") +
  theme_grbv(base_size = 9)

fig <- (pA / (pB + pC)) + plot_layout(heights = c(0.9, 1)) +
  plot_annotation(title = "S76 — Vine-level heterogeneity (tolerant vs susceptible)",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 180, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 180, units = "mm")
cat("[S76] saved\n")
