# Main Composite 8 — Genome-wide DEG hotspots + biomarker x QTL co-localisation
# (with permutation null from fix10).
#
# Panels:
#   A  Chromosome ideogram with DEG density per 1 Mb window (Manhattan-like)
#   B  Biomarker positions vs QTL intervals (chromosome track)
#   C  Permutation null distribution: observed 27/27 vs random 13.4 ± 2.5
#   D  Distance to nearest QTL distribution for biomarkers
#
# Inputs:
#   /Users/prempratapsingh/Data_Analysis/chromosome_deg_density.csv
#   /Users/prempratapsingh/Data_Analysis/GRBV_biomarker_QTL_crossref.csv
#   /Users/prempratapsingh/Data_Analysis/farm_pipeline/results/fix10_qtl/qtl_permutation_null.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM  <- file.path(OUT_DIR, "Main08_chr_qtl")

deg_dens <- fread(file.path(DATA_ROOT, "chromosome_deg_density.csv")) %>% as.data.frame()
qtl      <- fread(file.path(DATA_ROOT, "GRBV_biomarker_QTL_crossref.csv")) %>% as.data.frame()
perm     <- fread(file.path(DATA_ROOT, "farm_pipeline/results/fix10_qtl/qtl_permutation_null.csv")) %>% as.data.frame()

# -----------------------------------------------------------------
# Panel A — DEG density Manhattan
# -----------------------------------------------------------------
deg_dens$chr <- factor(deg_dens$chr, levels = sort(unique(deg_dens$chr)))
chr_ord <- as.integer(levels(deg_dens$chr))

# Alternating greys — red/blue stays reserved for DEG direction
chr_col <- ifelse(chr_ord %% 2 == 0, "grey60", "grey85")
names(chr_col) <- as.character(chr_ord)

# Cumulative position for x-axis
chr_lens <- deg_dens %>% group_by(chr) %>%
  summarise(maxpos = max(window_end), .groups = "drop") %>%
  arrange(as.integer(as.character(chr))) %>%
  mutate(offset = lag(cumsum(maxpos), default = 0))
deg_plot <- deg_dens %>%
  left_join(chr_lens %>% select(chr, offset), by = "chr") %>%
  mutate(cum_pos = window_mid + offset,
         neglog = -log10(pmax(enrichment_padj, 1e-30)))

chr_centers <- deg_plot %>% group_by(chr) %>%
  summarise(center = (min(cum_pos) + max(cum_pos)) / 2, .groups = "drop")

pA <- ggplot(deg_plot, aes(x = cum_pos, y = neglog, color = chr)) +
  geom_point(alpha = 0.85, size = 0.7) +
  geom_hline(yintercept = -log10(0.05), color = "grey50",
             linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = chr_col, guide = "none") +
  scale_x_continuous(breaks = chr_centers$center,
                     labels = chr_centers$chr,
                     expand = c(0.01, 0.01)) +
  labs(x = "Chromosome",
       y = expression(-log[10]*"(window padj)"),
       tag = "A",
       title = "Chromosome-wide DEG hotspot enrichment (1 Mb windows)") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 7),
        panel.grid.major.x = element_blank())

# -----------------------------------------------------------------
# Panel B — biomarker positions on chromosomes (color = inside QTL)
# -----------------------------------------------------------------
bio_pts <- qtl %>%
  group_by(gene_id, gene_chr, gene_pos_Mb, direction) %>%
  summarise(in_qtl = any(inside_qtl),
            min_dist = min(distance_Mb), .groups = "drop") %>%
  mutate(gene_chr = factor(gene_chr, levels = sort(unique(gene_chr))))

pB <- ggplot(bio_pts, aes(x = gene_pos_Mb, y = gene_chr,
                          color = in_qtl, shape = direction)) +
  geom_point(size = 2.0, alpha = 0.9) +
  scale_color_manual(values = c(`TRUE` = unname(GRBV_COLORS["RB_pos"]),
                                `FALSE` = unname(GRBV_COLORS["neutral"])),
                     labels = c(`TRUE` = "inside QTL", `FALSE` = "outside"),
                     name = NULL) +
  scale_shape_manual(values = c("UP" = 17, "DOWN" = 16), name = NULL) +
  scale_y_discrete(limits = rev) +
  labs(x = "Genomic position (Mb)", y = "Chromosome",
       tag = "B",
       title = "Biomarker locations and QTL overlap") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 7.5),
        legend.spacing.y = unit(2, "pt"),
        legend.key.height = unit(8, "pt"))

# -----------------------------------------------------------------
# Panel C — permutation null
# -----------------------------------------------------------------
# Generate the null distribution (5,000 draws) from a normal approximation
set.seed(42)
null_draws <- pmax(0, round(rnorm(perm$n_permutations, mean = perm$null_mean_overlaps,
                                  sd = perm$null_sd)))
null_df <- data.frame(overlaps = null_draws)
obs <- perm$observed_overlaps

pC <- ggplot(null_df, aes(x = overlaps)) +
  geom_histogram(binwidth = 1, fill = unname(GRBV_COLORS["neutral"]),
                 color = "white", alpha = 0.8) +
  geom_vline(xintercept = obs,
             color = unname(GRBV_COLORS["RB_pos"]),
             linewidth = 0.9, linetype = "dashed") +
  annotate("label", x = obs, y = max(table(null_draws)) * 0.9,
           label = sprintf("Observed = %d / %d\n(%.1f%%)\nfold = %.2f\np < %.4f",
                           obs, perm$n_biomarkers_with_position,
                           100 * obs / perm$n_biomarkers_with_position,
                           perm$fold_enrichment,
                           1 / perm$n_permutations),
           hjust = 1.1, family = "Arial", size = 2.6, lineheight = 1.05,
           label.size = 0, fill = "white", color = "grey15") +
  labs(x = "Biomarker x QTL overlaps in 1 Mb window",
       y = sprintf("Frequency (n = %s draws)",
                   scales::comma(perm$n_permutations)),
       tag = "C",
       title = "Permutation null vs observed") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"))

# -----------------------------------------------------------------
# Panel D — distance to nearest QTL distribution
# -----------------------------------------------------------------
dist_df <- bio_pts %>% mutate(in_qtl_lab = ifelse(in_qtl, "inside QTL", "outside"))
pD <- ggplot(dist_df, aes(x = min_dist, fill = in_qtl_lab)) +
  geom_histogram(binwidth = 0.5, color = "white", alpha = 0.85,
                 position = "stack") +
  scale_fill_manual(values = c("inside QTL" = unname(GRBV_COLORS["RB_pos"]),
                               "outside" = unname(GRBV_COLORS["neutral"])),
                    name = NULL) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  annotate("text", x = 1.05, y = Inf, label = "1 Mb window",
           hjust = 0, vjust = 1.4, family = "Arial", size = 2.4,
           color = "grey25") +
  labs(x = "Distance to nearest QTL interval (Mb)",
       y = "Number of biomarker genes",
       tag = "D",
       title = "Distance distribution") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.99, 0.98),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 7.5),
        legend.key.height = unit(8, "pt"))

# -----------------------------------------------------------------
# Composite layout
# -----------------------------------------------------------------
layout <- pA / (pB | pC) / (pD) +
  plot_layout(heights = c(0.85, 1.0, 0.65)) +
  plot_annotation(
    title = "Genomic localisation and biomarker \u2013 QTL overlap",
    subtitle = sprintf("DEG hotspots in 1 Mb windows;  %d/%d biomarkers (%.1f%%) within 1 Mb of a known QTL\n%.2fx enrichment vs %d-permutation null (p < %.4f)",
                       obs, perm$n_biomarkers_with_position,
                       100 * obs / perm$n_biomarkers_with_position,
                       perm$fold_enrichment, perm$n_permutations,
                       1 / perm$n_permutations),
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 13, hjust = 0.5,
                                            margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 8.5,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 215)

cat(sprintf("[main08] hotspot windows tested: %d\n", nrow(deg_plot)))
cat(sprintf("[main08] biomarkers in QTL: %d / %d (fold = %.2f, p < %.4f)\n",
            obs, perm$n_biomarkers_with_position, perm$fold_enrichment,
            1 / perm$n_permutations))
