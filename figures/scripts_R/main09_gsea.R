# Main Composite 9 — GSEA across stage x season conditions.
#
# Panels:
#   A  NES heatmap of top 30 most consistently regulated gene sets across the
#      8 stage x season conditions
#   B  Up vs Down regulated gene-set count per condition (bar chart)
#   C  Top-15 most-suppressed gene sets (mean NES across all conditions)
#   D  Cross-season NES correlation per stage (mirrors inverted trajectory)
#
# Inputs:
#   /Users/prempratapsingh/Data_Analysis/gsea_results_all_conditions.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM  <- file.path(OUT_DIR, "Main09_gsea")

gsea <- fread(file.path(DATA_ROOT, "gsea_results_all_conditions.csv")) %>%
  as.data.frame() %>%
  mutate(condition = gsub("_", " ", condition)) %>%
  separate(condition, into = c("season","stage_a","stage_b"),
           sep = " ", remove = FALSE, fill = "right") %>%
  mutate(stage = ifelse(is.na(stage_b), stage_a,
                        paste(stage_a, stage_b, sep = "_")),
         stage = factor(stage, levels = GRBV_STAGES),
         season = factor(season, levels = c("2021","2022")),
         sig = !is.na(fdr) & fdr < 0.05)

# -----------------------------------------------------------------
# Panel A — NES heatmap (top 30 by max |NES| across conditions)
# -----------------------------------------------------------------
top_sets <- gsea %>% filter(sig) %>%
  group_by(gene_set) %>%
  summarise(max_abs = max(abs(nes), na.rm = TRUE),
            n_cond = n_distinct(condition[sig]),
            .groups = "drop") %>%
  arrange(desc(n_cond), desc(max_abs)) %>%
  head(20)

heat <- gsea %>% filter(gene_set %in% top_sets$gene_set) %>%
  mutate(condition_full = paste(season, GRBV_STAGE_LABEL[as.character(stage)],
                                sep = " "),
         condition_full = factor(condition_full,
                                 levels = unique(condition_full[order(season, stage)])),
         gene_set = factor(gene_set, levels = rev(top_sets$gene_set)))

pA <- ggplot(heat, aes(x = condition_full, y = gene_set, fill = nes)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(data = subset(heat, sig),
            aes(label = ifelse(abs(nes) > 1.5, "*", "")),
            family = "Arial", size = 2.6,
            color = ifelse(subset(heat, sig)$nes > 0, "white", "white")) +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, name = "NES",
                       limits = c(-3, 3), oob = scales::squish) +
  labs(x = NULL, y = NULL, tag = "A",
       title = "Top recurrently-perturbed gene sets (* = FDR<0.05 with |NES|>1.5)") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 9.5, face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5),
        legend.key.height = unit(20, "pt"))

# -----------------------------------------------------------------
# Panel B — Up vs Down gene-set counts per condition
# -----------------------------------------------------------------
cond_counts <- gsea %>% filter(sig) %>%
  mutate(direction = ifelse(nes > 0, "up", "down"),
         condition_full = paste(season, GRBV_STAGE_LABEL[as.character(stage)],
                                sep = " ")) %>%
  count(season, stage, condition_full, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(condition_full = factor(condition_full,
                                 levels = unique(condition_full[order(season, stage)]))) %>%
  pivot_longer(c(up, down), names_to = "direction", values_to = "n")

pB <- ggplot(cond_counts, aes(x = condition_full, y = n, fill = direction)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7,
           color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                               down = unname(GRBV_COLORS["down"])),
                    labels = c(up = "Enriched (NES>0)",
                               down = "Suppressed (NES<0)"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Significant gene sets (FDR<0.05)",
       tag = "B",
       title = "Significant sets per condition") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 9.5, face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
        legend.position = "top",
        legend.text = element_text(size = 7.5))

# -----------------------------------------------------------------
# Panel C — top 15 most-suppressed gene sets (mean NES)
# -----------------------------------------------------------------
suppressed <- gsea %>%
  group_by(gene_set) %>%
  summarise(mean_nes = mean(nes, na.rm = TRUE),
            n_sig = sum(sig), .groups = "drop") %>%
  filter(n_sig >= 4) %>%
  arrange(mean_nes) %>%
  head(15) %>%
  mutate(label = factor(gene_set, levels = rev(gene_set)))

pC <- ggplot(suppressed, aes(x = mean_nes, y = label, fill = mean_nes)) +
  geom_col(alpha = 0.92, width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("n=%d", n_sig)),
            hjust = 1.18, family = "Arial", size = 2.4, color = "white") +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.05))) +
  labs(x = "Mean NES (across 8 conditions)", y = NULL, tag = "C",
       title = "Most suppressed gene sets") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 9.5, face = "bold"),
        axis.text.y = element_text(size = 6.5))

# -----------------------------------------------------------------
# Panel D — cross-season NES correlation per stage
# -----------------------------------------------------------------
nes_wide <- gsea %>% select(gene_set, season, stage, nes) %>%
  pivot_wider(names_from = season, values_from = nes, names_prefix = "y") %>%
  filter(!is.na(y2021), !is.na(y2022))

cor_per_stage <- nes_wide %>% group_by(stage) %>%
  summarise(r = cor(y2021, y2022, use = "pairwise.complete.obs"),
            n = n(), .groups = "drop") %>%
  mutate(stage = factor(stage, levels = GRBV_STAGES),
         label = sprintf("%s\nr = %.2f  (n = %d)",
                         GRBV_STAGE_LABEL[as.character(stage)], r, n))

nes_wide$stage <- factor(nes_wide$stage,
                         levels = cor_per_stage$stage,
                         labels = cor_per_stage$label)

pD <- ggplot(nes_wide, aes(x = y2021, y = y2022)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "grey50",
              linetype = "dashed", linewidth = 0.3) +
  geom_point(alpha = 0.4, size = 0.7, color = unname(GRBV_COLORS["RB_pos"])) +
  facet_wrap(~ stage, ncol = 4) +
  coord_fixed(xlim = c(-3, 3), ylim = c(-3, 3)) +
  labs(x = "NES (2021)", y = "NES (2022)", tag = "D",
       title = "Cross-season NES agreement (mirrors inverted DEG trajectory)") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 9.5, face = "bold"),
        strip.text  = element_text(size = 7.5, face = "bold", lineheight = 1.0))

# -----------------------------------------------------------------
# Composite layout
# -----------------------------------------------------------------
layout <- pA / (pB | pC) / pD +
  plot_layout(heights = c(1.3, 1.0, 0.85)) +
  plot_annotation(
    title = "Gene-set enrichment (GSEA) across stage \u00d7 season conditions",
    subtitle = "375 gene sets \u00d7 8 conditions (500 permutations); 3,000 tests, 270 significant\n77% of significant calls are suppression (NES<0)",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 13, hjust = 0.5,
                                            margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 8.5,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)))
  )

# Widen composite so: (i) pA can fit 8 rotated condition labels; (ii) pB
# + pC at half-width each still give 110+ mm per panel, enough for pB's 8
# stage x season labels and pC's long gene-set names; (iii) pD 1x4 strip
# text "Pre-veraison\nr = ..." fits without clipping.
save_grbv(layout, OUT_STEM,
          width_mm = 230, height_mm = 280)

n_sig <- sum(gsea$sig)
n_up  <- sum(gsea$sig & gsea$nes > 0)
n_dn  <- sum(gsea$sig & gsea$nes < 0)
cat(sprintf("[main09] %d total set\u00d7cond tests, %d significant (%.1f%% suppressed)\n",
            nrow(gsea), n_sig, 100 * n_dn / n_sig))
cat("[main09] cross-season correlation per stage:\n"); print(cor_per_stage)
