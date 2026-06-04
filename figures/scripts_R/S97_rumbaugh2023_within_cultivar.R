# S97 — Within-cultivar replication against Rumbaugh et al. 2023 Merlot (Paso Robles).
#
# Rumbaugh 2023 is the only other published GRBV transcriptome in Merlot:
# same cultivar (Merlot clone 12), different rootstock (1103P), different
# vineyard (Paso Robles), different year (2019), focused on cell-wall +
# phenylpropanoid + pathogenesis-related + hexose-transporter biology
# (66 genes total, of which 64 map to LOC IDs and 51 are in our expressed-
# gene background).
#
# Four panels:
#
#   A. Per-gene log2FC concordance scatter (Rumbaugh harvest vs our Q6 pooled).
#      Each point is one Rumbaugh gene, colour by category. The within-cultivar
#      replication is the diagonal trend.
#
#   B. Per-category overall Spearman concordance bars across all stages
#      (Rumbaugh log2FC vs our Q2 season-pooled), with sign agreement at our
#      padj < 0.05 annotated. Cell wall metabolism is the strongest category.
#
#   C. Per-(category, stage) Spearman heatmap. Shows where in development
#      each category replicates strongest.
#
#   D. Q1b cohort gradient. The Rumbaugh 2023 concordance is broadly flat
#      across our 5 onset-cohort contrasts, unlike Rumbaugh 2022 where it
#      peaked at cohort2020. This suggests their Paso Robles vineyard's
#      transcriptome is cohort-age-portable for the cell-wall-centric
#      gene set they profiled.
#
# Source CSVs (produced by farm_pipeline/scripts/rumbaugh2023_concordance.py):
#   rumbaugh2023_gene_concordance.csv
#   rumbaugh2023_per_category_summary.csv
#   rumbaugh2023_per_stage_concordance.csv
#   rumbaugh2023_q1b_concordance.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(stringr)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT    <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR <- file.path(ROOT, "figures/supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S97_rumbaugh2023_within_cultivar")

# Friendly category labels
CAT_LABS <- c(
  "Cell wall metabolism"                   = "Cell wall metabolism",
  "Phenylpropanoid pathway"                = "Phenylpropanoid pathway",
  "Pathogenesis-related protein synthesis" = "PR proteins",
  "Hexose transporters"                    = "Hexose transporters"
)
cat_order <- c("Cell wall metabolism", "Phenylpropanoid pathway",
               "PR proteins", "Hexose transporters")

# Distinct, colour-blind-friendly palette per category
cat_pal <- c(
  "Cell wall metabolism"  = "#1F8A70",   # bluish-green
  "Phenylpropanoid pathway" = "#B5174E", # deep red
  "PR proteins"            = "#2E70B5",  # blue
  "Hexose transporters"    = "#E69F00"   # orange
)

STAGE_ORDER <- c("Pre_veraison", "Veraison", "Post_veraison", "Harvest")
STAGE_LABS  <- c(Pre_veraison = "Pre-V", Veraison = "Veraison",
                 Post_veraison = "Post-V", Harvest = "Harvest")

# ---------------------------------------------------------------------------
# Panel A — per-gene log2FC concordance scatter (Rumbaugh harvest vs our Q6)
# ---------------------------------------------------------------------------
gene_df <- read_csv(file.path(ROOT, "rumbaugh2023_gene_concordance.csv"),
                    show_col_types = FALSE)

# Use Rumbaugh harvest log2fc on x; project Q6 pooled log2fc on y
# Q6 pooled column does not exist in the wide table; reconstruct from Q2 Harvest
# as a close proxy when Q6 was not pivoted, but the gene_concordance file does
# carry the harvest log2fc column.
# Looking at column names: "ours_Q2_Harvest_log2fc_ours" is the Q2 harvest cell.

scatter <- gene_df %>%
  rename(rumb_harvest = log2fc_H,
         ours_harvest = `ours_Q2_Harvest_log2fc_ours`,
         ours_harvest_padj = `ours_Q2_Harvest_padj_ours`) %>%
  filter(!is.na(rumb_harvest), !is.na(ours_harvest)) %>%
  mutate(category = factor(category, levels = cat_order),
         is_sig = !is.na(ours_harvest_padj) & ours_harvest_padj < 0.05)

# Limits for symmetric coord_fixed
lim_x <- max(abs(scatter$rumb_harvest), na.rm = TRUE)
lim_y <- max(abs(scatter$ours_harvest), na.rm = TRUE)
lim   <- max(lim_x, lim_y) * 1.05

# Concordance stats annotation
spear <- cor(scatter$rumb_harvest, scatter$ours_harvest, method = "spearman")
n     <- nrow(scatter)
sign_pct <- mean(sign(scatter$rumb_harvest) == sign(scatter$ours_harvest)) * 100

pA <- ggplot(scatter, aes(x = rumb_harvest, y = ours_harvest, colour = category)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50", linewidth = 0.4) +
  geom_point(aes(shape = is_sig, size = is_sig), alpha = 0.85) +
  scale_colour_manual(values = cat_pal, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "padj < 0.05", `FALSE` = "n.s."),
                     name = "Our Q2 Harvest") +
  scale_size_manual(values = c(`TRUE` = 2.6, `FALSE` = 1.6), guide = "none") +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "(A) Per-gene log2FC concordance at Harvest, Rumbaugh 2023 (Paso Robles) vs ours (Oakville)",
       subtitle = sprintf(
         "Same cultivar (Merlot), different vineyard and year. n = %d genes, Spearman rho = %.2f, %.0f%% sign agreement overall.",
         n, spear, sign_pct),
       x = "Rumbaugh 2023 log2 fold-change at Harvest",
       y = "Our Q2 Harvest log2 fold-change") +
  theme_grbv(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

# ---------------------------------------------------------------------------
# Panel B — per-category overall concordance bars
# ---------------------------------------------------------------------------
per_cat <- read_csv(file.path(ROOT, "rumbaugh2023_per_category_summary.csv"),
                    show_col_types = FALSE)

panelB_df <- per_cat %>%
  filter(category != "ALL", category != "Hexose transporters") %>%
  mutate(label = CAT_LABS[category],
         label = factor(label, levels = cat_order),
         annot = sprintf(
           "rho = %.2f\nn_genes = %d\n%.0f%% sign agree\n(n_sig = %d)",
           spearman_rho, n_unique_genes,
           100 * sign_agreement_sig, n_sig_padj))

pB <- ggplot(panelB_df, aes(x = label, y = spearman_rho, fill = label)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = annot), vjust = -0.18, size = 2.8, lineheight = 0.85,
            colour = "grey15") +
  scale_fill_manual(values = cat_pal, guide = "none") +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "(B) Cell wall biology is the strongest replication, with 100% sign agreement at significant genes",
       subtitle = "Spearman rho on flattened per-gene per-stage log2FC vector; significance computed in our Q2 (season-pooled) DESeq2 output",
       x = "Functional category",
       y = "Spearman rho (flattened across stages)") +
  theme_grbv(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# ---------------------------------------------------------------------------
# Panel C — per-(category, stage) Spearman heatmap
# ---------------------------------------------------------------------------
per_stage <- read_csv(file.path(ROOT, "rumbaugh2023_per_stage_concordance.csv"),
                      show_col_types = FALSE)

heat_df <- per_stage %>%
  filter(contrast %in% paste0("Q2_", STAGE_ORDER),
         category != "Hexose transporters") %>%
  mutate(category = factor(CAT_LABS[category], levels = cat_order),
         stage = factor(stage, levels = STAGE_ORDER),
         label = ifelse(is.na(spearman_rho),
                        ifelse(n < 3, sprintf("n=%d", n), "NA"),
                        sprintf("%.2f\n(n=%d)", spearman_rho, n)))

pC <- ggplot(heat_df, aes(x = stage, y = category, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.9, lineheight = 0.85, colour = "grey15") +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#1B7837",
                       midpoint = 0, limits = c(-1, 1),
                       na.value = "grey90", name = "Spearman\nrho") +
  scale_x_discrete(labels = STAGE_LABS) +
  labs(title = "(C) Per-category concordance is strongest at Harvest, where most Rumbaugh data sits",
       subtitle = "Spearman rho between Rumbaugh 2023 per-stage log2FC and our Q2 season-pooled log2FC; grey cells = fewer than 3 genes with data in that stage",
       x = "Ripening stage", y = NULL) +
  theme_grbv(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

# ---------------------------------------------------------------------------
# Panel D — Q1b cohort gradient
# ---------------------------------------------------------------------------
q1b <- read_csv(file.path(ROOT, "rumbaugh2023_q1b_concordance.csv"),
                show_col_types = FALSE)

q1b_plot <- q1b %>%
  filter(category != "Hexose transporters") %>%
  mutate(category = factor(if_else(category == "ALL", "All categories",
                                   CAT_LABS[category]),
                           levels = c("All categories", cat_order)),
         cohort_year = as.integer(str_extract(cohort, "\\d{4}")),
         years_label = case_when(
           cohort_year == 2016 ~ "2016\n(5-6 y)",
           cohort_year == 2019 ~ "2019\n(2-3 y)",
           cohort_year == 2020 ~ "2020\n(1-2 y)",
           cohort_year == 2021 ~ "2021\n(0-1 y)",
           cohort_year == 2022 ~ "2022\n(year 1)"),
         years_label = factor(years_label, levels = c(
           "2016\n(5-6 y)", "2019\n(2-3 y)", "2020\n(1-2 y)",
           "2021\n(0-1 y)", "2022\n(year 1)")))

cat_pal_with_all <- c(cat_pal, "All categories" = "grey25")

pD <- ggplot(q1b_plot, aes(x = years_label, y = spearman_rho, colour = category,
                            group = category)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(aes(linewidth = ifelse(category == "All categories", 1.1, 0.8)),
            show.legend = FALSE) +
  geom_point(aes(size = ifelse(category == "All categories", 2.8, 2.0)),
             show.legend = FALSE) +
  geom_line(linewidth = 0.0) +
  geom_point(size = 0) +
  scale_colour_manual(values = cat_pal_with_all, name = NULL) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2)) +
  labs(title = "(D) Within-cultivar agreement is broadly cohort-age-portable",
       subtitle = "Spearman rho between Rumbaugh 2023 harvest log2FC and each of our infection-onset cohort contrasts. Cell wall and phenylpropanoid categories agree across the full age gradient, in contrast to Rumbaugh 2022 (Cabernet Sauvignon) where the agreement peaked at our 1-2 year cohort.",
       x = "Our infection-onset cohort (years post-infection at 2021 sampling)",
       y = "Spearman rho (Rumbaugh harvest vs Q1b cohort)") +
  theme_grbv(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
fig <- (pA / pB / pC / pD) +
  plot_layout(heights = c(1.15, 1.0, 0.85, 1.0)) +
  plot_annotation(
    title    = "S97 - Within-cultivar replication against Rumbaugh 2023 Merlot (Paso Robles)",
    subtitle = paste0(
      "66 genes published in Rumbaugh 2023 Table 1 (cell wall, phenylpropanoid, PR proteins, hexose transporters); ",
      "64 map to LOC, 51 expressed in our data. Same cultivar (Merlot), different rootstock (1103P), ",
      "different vineyard (Paso Robles), different year (2019). Overall Spearman rho = 0.72 across 51 genes x stages, ",
      "95 percent sign agreement among our padj < 0.05 genes. Cell-wall biology replicates most strongly (rho = 0.69, 100 percent sign agreement)."),
    theme = theme(plot.title    = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey25",
                                                margin = margin(b = 6)))
  )

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 360, units = "mm", dpi = 600)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 360, units = "mm")

cat(sprintf(
  "[S97] wrote PNG + SVG to %s\n  Panels: A=harvest scatter, B=per-category bars, C=stage heatmap, D=cohort gradient\n",
  OUT_STEM))
