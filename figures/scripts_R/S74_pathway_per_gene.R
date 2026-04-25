# S74 — Per-gene Phenylpropanoid + PR + Cell-wall Inf x Season x Stage x Cohort
# model.  Goes one level deeper than the existing GSEA pathway-level enrichment
# (S70) by ranking individual genes within each curated pathway.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/14_pathway_per_gene"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S74_pathway_per_gene")

CAT_LEVELS <- c("Phenylpropanoid", "PR_proteins", "Cell_wall")
CAT_LABELS <- c("Phenylpropanoid", "PR proteins", "Cell wall")
CAT_COLORS <- c("#1F77B4", "#D62728", "#2CA02C")

# ---- panel A: category-level summary ----------------------------------------
sumr <- read_csv(file.path(RESULTS, "category_summary.csv"), show_col_types = FALSE) %>%
  mutate(category_label = factor(category, levels = CAT_LEVELS, labels = CAT_LABELS),
         pct_inf = frac_sig_infection * 100,
         pct_3way = n_sig_3way / n_tested * 100)

sumr_long <- sumr %>%
  select(category_label, pct_inf, pct_3way) %>%
  pivot_longer(c(pct_inf, pct_3way), names_to = "term", values_to = "pct") %>%
  mutate(term = recode(term,
                       pct_inf  = "Infection main",
                       pct_3way = "Inf x Season x Stage"))

pA <- ggplot(sumr_long, aes(x = category_label, y = pct, fill = term)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_dodge(0.8), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Infection main" = "#B2182B",
                               "Inf x Season x Stage" = "#2166AC"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "(A) Per-category significant gene fraction (p<0.05)",
       subtitle = "n_tested per category in axis annotation",
       x = NULL, y = "% significant") +
  theme_grbv(base_size = 9) +
  theme(legend.position = "top")

# ---- panel B: top 15 hits per category --------------------------------------
top <- read_csv(file.path(RESULTS, "top_hits_per_category.csv"), show_col_types = FALSE) %>%
  group_by(category_pri) %>%
  arrange(`p_C(infection)`) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  mutate(category_label = factor(category_pri, levels = CAT_LEVELS, labels = CAT_LABELS),
         minus_log10p = -log10(pmax(`p_C(infection)`, 1e-30)),
         dir = ifelse(beta_infection > 0, "UP in RB+", "DOWN in RB+"))

pB <- ggplot(top, aes(x = beta_infection, y = reorder(gene_id, beta_infection),
                       fill = category_label, size = minus_log10p)) +
  geom_point(alpha = 0.85, shape = 21, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = setNames(CAT_COLORS, CAT_LABELS), name = NULL) +
  scale_size_continuous(range = c(2.5, 8), name = "-log10 p_inf") +
  labs(title = "(B) Top 15 per-category hits by infection-effect significance",
       x = "beta(infection): RB(+) - RB(-) on log2 CPM",
       y = NULL) +
  theme_grbv(base_size = 8) +
  theme(legend.position = "right")

fig <- (pA / pB) + plot_layout(heights = c(0.55, 1.45)) +
  plot_annotation(title = "S74 — Phenylpropanoid + PR + Cell-wall per-gene model",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 240, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 240, units = "mm")
cat("[S74] saved\n")
