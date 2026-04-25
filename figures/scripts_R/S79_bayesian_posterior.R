# S79 — Bayesian hierarchical posteriors (vine random intercept).  Six high-
# interest genes fit with brms ~ Infection + Season + Stage + (1 | vine_id),
# 4 chains x 2000 iter, weakly-informative priors.  Forest plot with 95% CIs
# and ROPE band [-0.1, +0.1].

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results/19_bayesian"
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S79_bayesian_posterior")

GENE_LABELS <- c(
  "LOC100268011" = "RdRp1 (top biomarker)",
  "LOC100854143" = "F-box LRR (top biomarker)",
  "LOC100257913" = "Glycine-rich (top DOWN)",
  "LOC100266309" = "LOC100266309 (mediator)",
  "LOC100232974" = "SPS1 (mediator)",
  "LOC100249249" = "MBF1c (top mediator)"
)

post <- read_csv(file.path(RESULTS, "bayes_posterior_summary.csv"), show_col_types = FALSE) %>%
  mutate(label = GENE_LABELS[gene],
         label = factor(label, levels = rev(GENE_LABELS)),
         straddles_zero = ci_lo < 0 & ci_hi > 0,
         in_rope_str = sprintf("p_ROPE = %.0f%%", p_in_rope * 100),
         color_class = case_when(
           p_in_rope > 0.10 ~ "Substantial ROPE overlap",
           ci_lo > 0 ~ "Positive (95% CI > 0)",
           ci_hi < 0 ~ "Negative (95% CI < 0)",
           TRUE ~ "Crosses zero"))

# ---- panel A: forest plot --------------------------------------------------
pA <- ggplot(post, aes(x = posterior_mean, y = label, color = color_class)) +
  annotate("rect", xmin = -0.1, xmax = 0.1, ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f [%.2f, %.2f]", posterior_mean, ci_lo, ci_hi)),
            x = 2.35, hjust = 0, size = 2.7, color = "grey25") +
  geom_text(aes(label = in_rope_str), x = 3.40, hjust = 0, size = 2.7, color = "grey25") +
  scale_color_manual(values = c("Positive (95% CI > 0)" = "#B2182B",
                                "Negative (95% CI < 0)" = "#2166AC",
                                "Crosses zero" = "#999999",
                                "Substantial ROPE overlap" = "#E69F00"),
                     name = NULL) +
  scale_x_continuous(limits = c(-2, 4.5), breaks = seq(-2, 3, 1)) +
  labs(title = "S79 — Bayesian posterior intervals: infection effect on log2 CPM",
       subtitle = "brms with vine_id random intercept; ROPE = [-0.1, +0.1]; 4 chains x 2000 iter",
       x = "Posterior mean (95% credible interval)", y = NULL) +
  theme_grbv(base_size = 9) +
  theme(legend.position = "top",
        plot.subtitle = element_text(color = "grey40", size = 8))

ggsave(paste0(OUT_STEM, ".png"), pA,
       width = GRBV$W_DOUBLE_MM, height = 110, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), pA,
       width = GRBV$W_DOUBLE_MM, height = 110, units = "mm")
cat("[S79] saved\n")
