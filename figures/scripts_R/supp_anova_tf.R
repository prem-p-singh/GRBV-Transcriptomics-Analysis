# Supplementary: ANOVA variance decomposition + TF family detail —
# S40..S45.
#
#   S40  Variance decomposition: per-factor eta-sq distribution + mean
#   S41  Interaction-term gene counts (padj<0.05) per interaction
#   S42  Top Infection x Season interaction gene profiles
#   S43  TF family heatmap (%DE per condition)
#   S44  TF family balloon plot (mean |log2FC| x %DE per condition)
#   S45  Top TF gene expression boxplots per family

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(ggrepel)
})
HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))
DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- file.path(DATA_ROOT, "figures", "supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stages <- GRBV_STAGES

# -----------------------------------------------------------------
# ANOVA
# -----------------------------------------------------------------
anv <- fread(file.path(DATA_ROOT, "multifactor_anova_results.csv")) %>% as.data.frame()

eta_long <- anv %>%
  select(gene, starts_with("eta2_")) %>%
  pivot_longer(-gene, names_to = "factor", values_to = "eta2",
               names_prefix = "eta2_") %>%
  mutate(factor = recode(factor,
                         "Inf" = "Infection",
                         "Sea" = "Season",
                         "Stg" = "Stage",
                         "Inf:Sea" = "Inf x Sea",
                         "Inf:Stg" = "Inf x Stg",
                         "Sea:Stg" = "Sea x Stg",
                         "Inf:Sea:Stg" = "Inf x Sea x Stg"))
mean_eta <- eta_long %>% group_by(factor) %>%
  summarise(mean_eta2 = mean(eta2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_eta2))
order_fac <- mean_eta$factor
eta_long$factor <- factor(eta_long$factor, levels = order_fac)

pal_fac <- c(
  "Stage" = unname(GRBV_COLORS["stage_post"]),
  "Sea x Stg" = unname(GRBV_COLORS["season_2022"]),
  "Inf x Stg" = unname(GRBV_COLORS["RB_pos"]),
  "Season" = unname(GRBV_COLORS["season_2021"]),
  "Inf x Sea x Stg" = unname(GRBV_COLORS["highlight"]),
  "Infection" = unname(GRBV_COLORS["RB_neg"]),
  "Inf x Sea" = unname(GRBV_COLORS["warm"])
)

# S40 — per-factor eta2 boxplot + mean bar
pA40 <- ggplot(eta_long, aes(x = factor, y = eta2, fill = factor)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.25,
               color = "grey25", linewidth = 0.3, width = 0.65) +
  scale_fill_manual(values = pal_fac[order_fac], guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = expression(paste("Per-gene ", eta^2)),
       title = "A. Per-gene eta-squared distribution") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5))

pB40 <- ggplot(mean_eta %>% mutate(factor = factor(factor, levels = order_fac)),
               aes(x = mean_eta2, y = factor, fill = factor)) +
  geom_col(alpha = 0.92, width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("%.3f", mean_eta2)),
            hjust = -0.18, family = "Arial", size = 2.6) +
  scale_fill_manual(values = pal_fac[order_fac], guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_y_discrete(limits = rev) +
  labs(x = expression(paste("Mean per-gene ", eta^2)), y = NULL,
       title = "B. Mean variance share") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"))

S40 <- pA40 | pB40 +
  plot_annotation(title = "Supplementary S40. Variance decomposition (ANOVA on 14,604 genes)",
                  theme = theme(plot.title = element_text(family = "Arial",
                                                          face = "bold", size = 12,
                                                          hjust = 0.5)))
save_grbv(S40, file.path(OUT_DIR, "S40_variance_decomposition"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 120)

# -----------------------------------------------------------------
# S41 — interaction-term gene counts (padj<0.05) per interaction
# -----------------------------------------------------------------
PADJ <- 0.05
counts_df <- data.frame(
  interaction = c("Inf", "Sea", "Stg",
                  "Inf x Sea", "Inf x Stg", "Sea x Stg", "Inf x Sea x Stg"),
  n = c(sum(anv$padj_Inf < PADJ, na.rm = TRUE),
        sum(anv$padj_Sea < PADJ, na.rm = TRUE),
        sum(anv$padj_Stg < PADJ, na.rm = TRUE),
        sum(anv$`padj_Inf:Sea` < PADJ, na.rm = TRUE),
        sum(anv$`padj_Inf:Stg` < PADJ, na.rm = TRUE),
        sum(anv$`padj_Sea:Stg` < PADJ, na.rm = TRUE),
        sum(anv$`padj_Inf:Sea:Stg` < PADJ, na.rm = TRUE))
) %>% arrange(desc(n)) %>%
  mutate(interaction = factor(interaction, levels = rev(interaction)))

S41 <- ggplot(counts_df, aes(x = n, y = interaction, fill = interaction)) +
  geom_col(alpha = 0.92, width = 0.7, color = "white") +
  geom_text(aes(label = scales::comma(n)), hjust = -0.15,
            family = "Arial", size = 2.8) +
  scale_fill_manual(values = pal_fac[as.character(counts_df$interaction)],
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = scales::label_comma()) +
  labs(x = "Genes with significant term (padj < 0.05)", y = NULL,
       title = "Supplementary S41. Interaction-term gene counts") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"))
save_grbv(S41, file.path(OUT_DIR, "S41_interaction_gene_sets"),
          width_mm = GRBV$W_ONEHALF_MM, height_mm = 100)

# -----------------------------------------------------------------
# S42 — top Inf x Sea interaction gene profiles (mean lcpm per season x infection)
# -----------------------------------------------------------------
top_int <- anv %>% arrange(`padj_Inf:Sea`) %>% head(12) %>% pull(gene)

counts <- fread(file.path(DATA_ROOT, "counts_matrix.txt"),
                sep = "\t", header = TRUE, data.table = FALSE)
rownames(counts) <- counts[,1]; counts[,1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts); counts <- counts[keep, ]
lib <- colSums(counts); lcpm <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)
meta <- read_tsv(file.path(DATA_ROOT, "sample_metadata.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
rownames(meta) <- meta$Sample; meta <- meta[colnames(lcpm), ]

top_in <- intersect(top_int, rownames(lcpm))
prof <- as.data.frame(lcpm[top_in, ]) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "Sample", values_to = "lcpm") %>%
  left_join(meta, by = "Sample") %>%
  mutate(gene = factor(gene, levels = top_in),
         Sampling_Season = factor(Sampling_Season))

S42 <- ggplot(prof,
              aes(x = Sampling_Season, y = lcpm,
                  fill = Infection_Group)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3,
               linewidth = 0.25, width = 0.75) +
  facet_wrap(~ gene, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = GRBV_INFECTION_COLOR, name = NULL) +
  labs(x = "Season", y = "log2(CPM + 1)",
       title = "Supplementary S42. Top-12 Infection x Season interaction genes") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 7.5, face = "bold"),
        legend.position = "top")
save_grbv(S42, file.path(OUT_DIR, "S42_interaction_profiles"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 160)

# -----------------------------------------------------------------
# TF family
# -----------------------------------------------------------------
tf <- fread(file.path(DATA_ROOT, "tf_family_DE_summary.csv")) %>% as.data.frame() %>%
  separate(condition, into = c("season","sa","sb"), sep = "_",
           remove = FALSE, fill = "right") %>%
  mutate(stage = ifelse(is.na(sb), sa, paste(sa, sb, sep = "_")),
         stage = factor(stage, levels = stages),
         season = factor(season, levels = c("2021","2022")),
         tf_family = factor(tf_family,
                            levels = c("ERF/AP2","MYB","WRKY","bHLH","NAC","HSF")))

# S43 — TF family %DE heatmap
S43 <- ggplot(tf, aes(x = stage, y = tf_family, fill = percent_de)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", percent_de),
                color = percent_de > 12),
            family = "Arial", size = 2.4) +
  facet_wrap(~ season, ncol = 2) +
  scale_fill_gradient(low = "#fffaf0", high = unname(GRBV_COLORS["RB_pos"]),
                      name = "% DE") +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                     guide = "none") +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  labs(x = NULL, y = NULL,
       title = "Supplementary S43. TF family DEG percentages per stage x season") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
        strip.text = element_text(size = 9, face = "bold"),
        legend.key.height = unit(18, "pt"))
save_grbv(S43, file.path(OUT_DIR, "S43_tf_family_heatmap"),
          width_mm = GRBV$W_ONEHALF_MM, height_mm = 110)

# S44 — balloon plot (size = %DE, color = mean log2FC)
S44 <- ggplot(tf, aes(x = stage, y = tf_family)) +
  geom_point(aes(size = percent_de, fill = mean_log2FC),
             shape = 21, color = "grey15", stroke = 0.3) +
  facet_wrap(~ season, ncol = 2) +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white", high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, name = "Mean log2FC") +
  scale_size_continuous(name = "% DE", range = c(0.5, 7)) +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  labs(x = NULL, y = NULL,
       title = "Supplementary S44. TF family enrichment balloon plot") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
        strip.text = element_text(size = 9, face = "bold"),
        legend.key.height = unit(14, "pt"))
save_grbv(S44, file.path(OUT_DIR, "S44_tf_enrichment_balloon"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 110)

# S45 — top TF gene expression boxplots per family
tfg_f <- file.path(DATA_ROOT, "tf_family_gene_details.csv")
if (file.exists(tfg_f)) {
  tfg <- fread(tfg_f) %>% as.data.frame()
  # Pick top-2 by |log2FC| per family
  tfg_top <- tfg %>% group_by(tf_family) %>%
    slice_max(abs(log2FC), n = 2, with_ties = FALSE) %>%
    ungroup() %>% distinct(gene, tf_family, gene_name, log2FC)
  pick_genes <- intersect(tfg_top$gene, rownames(lcpm))
  if (length(pick_genes) > 0) {
    tfg_top <- tfg_top %>% filter(gene %in% pick_genes)
    expr45 <- as.data.frame(lcpm[pick_genes, , drop = FALSE]) %>%
      tibble::rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "Sample", values_to = "lcpm") %>%
      left_join(meta, by = "Sample") %>%
      left_join(tfg_top %>% select(gene, tf_family, gene_name),
                by = "gene") %>%
      mutate(facet = paste0(tf_family, ": ", gene))
    S45 <- ggplot(expr45, aes(x = Infection_Group, y = lcpm,
                              fill = Infection_Group)) +
      geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3,
                   linewidth = 0.25, width = 0.65) +
      facet_wrap(~ facet, ncol = 4, scales = "free_y") +
      scale_fill_manual(values = GRBV_INFECTION_COLOR, guide = "none") +
      labs(x = NULL, y = "log2(CPM + 1)",
           title = "Supplementary S45. Top-2 DE TF gene expression per family") +
      theme_grbv() +
      theme(plot.title = element_text(size = 11, face = "bold"),
            strip.text = element_text(size = 7, face = "bold"),
            axis.text.x = element_text(size = 7))
    save_grbv(S45, file.path(OUT_DIR, "S45_tf_expression_boxplots"),
              width_mm = GRBV$W_DOUBLE_MM, height_mm = 150)
  }
}

cat("[supp_anova_tf] wrote S40, S41, S42, S43, S44, S45\n")
