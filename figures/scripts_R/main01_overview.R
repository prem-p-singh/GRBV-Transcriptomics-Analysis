# Main Composite 1 — Study overview + variance decomposition.
#
# Panels:
#   A  Sample composition tile chart (stage x season x infection)
#   B  Global PCA on top-2,000 variable log-CPM (PC1 vs PC2),
#      coloured by stage, shaped by season
#   C  ANOVA variance decomposition: distribution of per-gene eta-sq
#      across the seven factors (boxplots ordered by mean)
#   D  Top-20 variance-explained genes per factor (mean eta-sq bar)
#
# Inputs:
#   /Users/prempratapsingh/Data_Analysis/counts_matrix.txt
#   /Users/prempratapsingh/Data_Analysis/sample_metadata.txt
#   /Users/prempratapsingh/Data_Analysis/multifactor_anova_results.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "Main01_overview")

# -----------------------------------------------------------------
# Load
# -----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
.cnt_dt <- fread(file.path(DATA_ROOT, "counts_matrix.txt"),
                 sep = "\t", header = TRUE, data.table = FALSE)
rownames(.cnt_dt) <- .cnt_dt[, 1]; .cnt_dt[, 1] <- NULL
counts <- as.matrix(.cnt_dt); storage.mode(counts) <- "integer"; rm(.cnt_dt)
meta   <- read_tsv(file.path(DATA_ROOT, "sample_metadata.txt"),
                   show_col_types = FALSE) %>% as.data.frame()
rownames(meta) <- meta$Sample
meta <- meta[colnames(counts), ]

meta$Maturity_stage  <- factor(meta$Maturity_stage, levels = GRBV_STAGES)
meta$Sampling_Season <- factor(meta$Sampling_Season)
meta$Infection_Group <- factor(meta$Infection_Group, levels = c("RB(-)","RB(+)"))

# Filter genes (>=10 in >=10% samples)
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]

# log-CPM
lib <- colSums(counts)
lcpm <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)

# Top-2000 variable genes
vars <- apply(lcpm, 1, var)
top2k <- names(sort(vars, decreasing = TRUE)[1:2000])

# PCA (samples in rows)
pcafit <- prcomp(t(lcpm[top2k, ]), scale. = TRUE)
ve <- (pcafit$sdev ^ 2) / sum(pcafit$sdev ^ 2)
pca_df <- data.frame(
  PC1 = pcafit$x[, 1],
  PC2 = pcafit$x[, 2],
  Sample = rownames(pcafit$x)
) %>% left_join(meta, by = "Sample")

# ANOVA
anova_df <- read_csv(file.path(DATA_ROOT, "multifactor_anova_results.csv"),
                     show_col_types = FALSE)
eta_long <- anova_df %>%
  select(starts_with("eta2_")) %>%
  pivot_longer(everything(), names_to = "factor", values_to = "eta2",
               names_prefix = "eta2_") %>%
  filter(!is.na(eta2)) %>%
  mutate(
    factor = recode(factor,
      "Inf"         = "Infection",
      "Sea"         = "Season",
      "Stg"         = "Stage",
      "Inf:Sea"     = "Inf \u00d7 Sea",
      "Inf:Stg"     = "Inf \u00d7 Stg",
      "Sea:Stg"     = "Sea \u00d7 Stg",
      "Inf:Sea:Stg" = "Inf \u00d7 Sea \u00d7 Stg"
    )
  )

mean_eta <- eta_long %>%
  group_by(factor) %>% summarise(mean_eta2 = mean(eta2), .groups = "drop") %>%
  arrange(desc(mean_eta2))
order_factors <- mean_eta$factor
eta_long$factor <- factor(eta_long$factor, levels = order_factors)
mean_eta$factor <- factor(mean_eta$factor, levels = order_factors)

# -----------------------------------------------------------------
# Panel A — sample composition tile chart
# -----------------------------------------------------------------
comp <- meta %>%
  count(Sampling_Season, Maturity_stage, Infection_Group, name = "n") %>%
  complete(Sampling_Season, Maturity_stage, Infection_Group, fill = list(n = 0))

pA <- ggplot(comp,
             aes(x = Maturity_stage, y = Infection_Group, fill = n)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n), family = "Arial", size = 2.7,
            color = ifelse(comp$n > max(comp$n) * 0.55, "white", "black")) +
  facet_wrap(~ Sampling_Season, ncol = 1, strip.position = "top") +
  scale_fill_gradient(low = "#fef9e7", high = unname(GRBV_COLORS["RB_pos"]),
                      guide = "none") +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  labs(x = NULL, y = NULL, tag = "A",
       title = "Sample composition") +
  theme_grbv() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 8),
        plot.title  = element_text(size = 10, face = "bold"),
        strip.text  = element_text(size = 9, face = "bold"))

# -----------------------------------------------------------------
# Panel B — global PCA
# -----------------------------------------------------------------
pB <- ggplot(pca_df,
             aes(x = PC1, y = PC2,
                 color = Maturity_stage, shape = Sampling_Season)) +
  geom_point(alpha = 0.85, size = 1.4, stroke = 0.4) +
  scale_color_manual(values = GRBV_STAGE_COLOR, labels = GRBV_STAGE_LABEL,
                     name = "Maturity stage") +
  scale_shape_manual(values = c("2021" = 16, "2022" = 17), name = "Season") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1] * 100),
       y = sprintf("PC2 (%.1f%%)", ve[2] * 100),
       tag = "B",
       title = "Global PCA (top-2,000 variable genes)") +
  theme_grbv() +
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = unit(2, "pt"),
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 7.5, face = "bold"),
        legend.key.height = unit(7, "pt"),
        plot.title = element_text(size = 10, face = "bold"))

# -----------------------------------------------------------------
# Panel C — variance decomposition boxplot
# -----------------------------------------------------------------
factor_color <- c(
  "Stage"             = unname(GRBV_COLORS["stage_post"]),
  "Sea \u00d7 Stg"    = unname(GRBV_COLORS["season_2022"]),
  "Inf \u00d7 Stg"    = unname(GRBV_COLORS["RB_pos"]),
  "Season"            = unname(GRBV_COLORS["season_2021"]),
  "Inf \u00d7 Sea \u00d7 Stg" = unname(GRBV_COLORS["highlight"]),
  "Infection"         = unname(GRBV_COLORS["RB_neg"]),
  "Inf \u00d7 Sea"    = unname(GRBV_COLORS["warm"])
)
pal_C <- factor_color[order_factors]

pC <- ggplot(eta_long,
             aes(x = factor, y = eta2, fill = factor)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.25,
               color = "grey25", linewidth = 0.3, width = 0.65) +
  scale_fill_manual(values = pal_C, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = expression(paste("Per-gene ", eta^2)),
       tag = "C",
       title = "Variance decomposition (n = 14,604 genes)") +
  theme_grbv() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5),
        plot.title = element_text(size = 10, face = "bold"))

# -----------------------------------------------------------------
# Panel D — mean eta-sq per factor
# -----------------------------------------------------------------
pD <- ggplot(mean_eta,
             aes(x = mean_eta2, y = factor, fill = factor)) +
  geom_col(alpha = 0.92, width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("%.3f", mean_eta2)),
            hjust = -0.18, family = "Arial", size = 2.7) +
  scale_fill_manual(values = pal_C, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_y_discrete(limits = rev) +
  labs(x = expression(paste("Mean per-gene ", eta^2)),
       y = NULL, tag = "D",
       title = "Mean variance share") +
  theme_grbv() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold"))

# -----------------------------------------------------------------
# Composite layout
# -----------------------------------------------------------------
layout <- (pA + pB + plot_layout(widths = c(1.0, 1.4))) /
  (pC + pD + plot_layout(widths = c(1.4, 1.0))) +
  plot_layout(heights = c(1.0, 0.95)) +
  plot_annotation(
    title = "Study overview and global variance structure",
    subtitle = sprintf("327 samples (264 RB+, 63 RB\u2212), 14,604 expressed genes; %s",
                       "stage explains the most transcriptome variance"),
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 13, hjust = 0.5,
                                            margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 175)

cat(sprintf("[main01] PC1 = %.1f%%, PC2 = %.1f%%; mean eta-sq Stage = %.3f\n",
            ve[1] * 100, ve[2] * 100,
            mean_eta$mean_eta2[mean_eta$factor == "Stage"]))
