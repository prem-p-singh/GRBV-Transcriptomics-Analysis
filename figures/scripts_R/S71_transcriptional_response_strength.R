# S71 — Transcriptional response strength from PCA across cohort x season x stage
#
# Metric:
#   For each infected RB(+) sample, we define "response strength" as the
#   Euclidean distance in PC1-PC2 space from the centroid of RB(-) samples
#   matched on season and maturity stage.  This normalises out the normal
#   season x stage variation (which dominates transcriptome variance,
#   eta^2=0.426) and isolates how far each infected sample has drifted
#   from its healthy baseline.
#
#   Input matrix: log2(CPM+1) of top 2,000 most variable genes (same
#   as Main01).  PCA = prcomp with centered, scaled input.
#
# Panels:
#   A. PCA biplot  - PC1 vs PC2, samples colored by Infection_Group,
#                    faceted by season, healthy samples as grey anchors
#   B. Response strength heatmap - cohort (y) x stage (x), faceted by
#                                  season; fill = mean PC-space distance
#                                  from matched RB(-) centroid
#   C. Ranking bar plot - top 20 cohort x season x stage cells sorted by
#                         response strength; colour = cohort
#
# Output: figures/supp/S71_transcriptional_response_strength.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S71_transcriptional_response_strength")

STAGES    <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STAGE_LAB <- c(Pre_veraison = "Pre-v\u00e9r", Veraison = "V\u00e9r",
               Post_veraison = "Post-v\u00e9r", Harvest = "Harvest")
COHORTS   <- c("2016","2019","2020","2021","2022")
SEASONS   <- c("2021","2022")

# Cohort palette: viridis 5-step (years since infection ordinal)
cohort_pal <- setNames(c("#440154","#3B528B","#21908C","#5DC863","#FDE725"),
                      COHORTS)

# ============================================================================
# LOAD
# ============================================================================
meta <- read_tsv(file.path(ROOT, "sample_metadata.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
meta$Maturity_stage  <- factor(meta$Maturity_stage, levels = STAGES)
meta$Sampling_Season <- as.character(meta$Sampling_Season)
meta$Cohort          <- as.character(meta$Years_Since_Infection)

counts <- fread(file.path(ROOT, "counts_matrix.txt"),
                sep = "\t", header = TRUE, data.table = FALSE)
rownames(counts) <- counts[, 1]; counts[, 1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"
keep   <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]
lib    <- colSums(counts)
lcpm   <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)

# Align metadata to expression column order
meta <- meta[match(colnames(lcpm), meta$Sample), ]

# Top-2000 variable genes (match Main01's "global PCA" panel)
row_vars <- apply(lcpm, 1, var)
top_idx  <- order(row_vars, decreasing = TRUE)[seq_len(min(2000, length(row_vars)))]
X <- lcpm[top_idx, ]

# PCA: rows = samples => transpose
pca <- prcomp(t(X), center = TRUE, scale. = TRUE)
pc_scores <- as.data.frame(pca$x[, 1:2])
names(pc_scores) <- c("PC1","PC2")
pc_scores$Sample <- rownames(pc_scores)
pc_scores <- pc_scores %>% left_join(meta, by = "Sample")

pc_var <- summary(pca)$importance[2, 1:2] * 100
cat(sprintf("[S71] PCA variance: PC1 = %.1f%%, PC2 = %.1f%%\n",
            pc_var[1], pc_var[2]))

# ============================================================================
# RESPONSE STRENGTH: distance from matched (Season, Stage) RB(-) centroid
# ============================================================================
rbn_centroids <- pc_scores %>% filter(Infection_Group == "RB(-)") %>%
  group_by(Sampling_Season, Maturity_stage) %>%
  summarise(cx = mean(PC1), cy = mean(PC2), .groups = "drop")

pc_scores <- pc_scores %>%
  left_join(rbn_centroids, by = c("Sampling_Season","Maturity_stage")) %>%
  mutate(response = sqrt((PC1 - cx)^2 + (PC2 - cy)^2))

# Aggregate RB(+) to cohort x season x stage cells
cell_resp <- pc_scores %>%
  filter(Infection_Group == "RB(+)",
         Cohort %in% COHORTS,
         Sampling_Season %in% SEASONS) %>%
  group_by(Cohort, Sampling_Season, Maturity_stage) %>%
  summarise(mean_resp = mean(response, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  rename(Season = Sampling_Season)

# Expand to full grid
full_grid <- expand_grid(Cohort = COHORTS, Season = SEASONS,
                         Maturity_stage = factor(STAGES, levels = STAGES))
cell_resp <- full_grid %>%
  left_join(cell_resp, by = c("Cohort","Season","Maturity_stage")) %>%
  mutate(Cohort = factor(Cohort, levels = rev(COHORTS)))

# ============================================================================
# Panel A: PCA biplot, faceted by season, colored by Infection_Group,
# RB(+) further shaded by Cohort
# ============================================================================
RBNEG <- "RB(-)"
pc_plot <- pc_scores %>%
  filter(Sampling_Season %in% SEASONS) %>%
  mutate(group = ifelse(Infection_Group == "RB(-)", RBNEG, Cohort),
         group = factor(group, levels = c(RBNEG, COHORTS)))

grp_pal <- c(setNames(unname(GRBV_COLORS["neutral"]), RBNEG),
             cohort_pal)

pA <- ggplot(pc_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = group, shape = Maturity_stage),
             size = 1.5, alpha = 0.8) +
  # Mark the season-stage RB(-) centroid as a black cross so the
  # "distance from baseline" metric is visually explicit
  geom_point(data = rbn_centroids,
             aes(x = cx, y = cy),
             inherit.aes = FALSE,
             shape = 4, size = 2.5, stroke = 1.1, colour = "black") +
  facet_wrap(~ Sampling_Season, ncol = 2) +
  scale_colour_manual(values = grp_pal, name = "Cohort", drop = FALSE) +
  scale_shape_manual(values = c(Pre_veraison = 16, Veraison = 17,
                                Post_veraison = 15, Harvest = 18),
                     labels = STAGE_LAB, name = "Stage") +
  labs(x = sprintf("PC1 (%.1f%%)", pc_var[1]),
       y = sprintf("PC2 (%.1f%%)", pc_var[2]),
       title = "A  PCA on top 2,000 variable genes",
       subtitle = "black \u00d7 = RB(\u2212) centroid for each season \u00d7 stage; response strength = distance to matched \u00d7") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.text    = element_text(face = "bold", size = 9),
        legend.position = "right",
        legend.text   = element_text(size = 7.5),
        legend.key.height = unit(9, "pt"))

# ============================================================================
# Panel B: Response strength heatmap (cohort x stage x season)
# ============================================================================
cell_resp$cell_label <- ifelse(!is.na(cell_resp$mean_resp),
                               sprintf("%.1f\n(n=%d)",
                                       cell_resp$mean_resp, cell_resp$n),
                               "")
resp_lim <- max(cell_resp$mean_resp, na.rm = TRUE)

pB <- ggplot(cell_resp, aes(x = Maturity_stage, y = Cohort)) +
  geom_tile(fill = "#E8E8E8", colour = "white", linewidth = 0.4) +
  geom_tile(aes(fill = mean_resp), colour = "white", linewidth = 0.4,
            data = cell_resp %>% filter(!is.na(mean_resp))) +
  geom_text(aes(label = cell_label),
            family = "Arial", size = 2.3, colour = "grey15",
            lineheight = 0.95) +
  scale_fill_gradient(low = "#F4EDE2",
                      high = unname(GRBV_COLORS["RB_pos"]),
                      name = "Response\nstrength\n(PC distance)",
                      na.value = NA, limits = c(0, resp_lim)) +
  facet_wrap(~ Season, ncol = 2) +
  scale_x_discrete(labels = STAGE_LAB) +
  labs(x = NULL, y = "Infection cohort",
       title = "B  Mean transcriptional response strength per cell",
       subtitle = "Euclidean distance in PC1\u2013PC2 space from matched season \u00d7 stage RB(\u2212) centroid") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.text    = element_text(face = "bold", size = 9),
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 8),
        axis.text.y   = element_text(size = 8),
        legend.key.width  = unit(8, "pt"),
        legend.key.height = unit(14, "pt"))

# ============================================================================
# Panel C: Top-ranked cells (horizontal bar plot)
# ============================================================================
rank_df <- cell_resp %>% filter(!is.na(mean_resp)) %>%
  arrange(desc(mean_resp)) %>%
  mutate(cell_label = sprintf("%s \u00b7 %s \u00b7 %s",
                              Season, STAGE_LAB[as.character(Maturity_stage)],
                              as.character(Cohort)),
         cell_label = factor(cell_label, levels = rev(cell_label)))

# Show top 20 for legibility
top20 <- rank_df %>% slice_head(n = 20)
top20$Cohort_val <- factor(as.character(top20$Cohort), levels = COHORTS)

pC <- ggplot(top20, aes(x = mean_resp, y = cell_label, fill = Cohort_val)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", mean_resp)),
            hjust = -0.18, family = "Arial", size = 2.4, colour = "grey15") +
  scale_fill_manual(values = cohort_pal, name = "Cohort") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Response strength (PC distance)",
       y = NULL,
       title = "C  Top 20 cells by transcriptional response strength",
       subtitle = "the strongest GRBV-driven transcriptome drifts across the 3-axis cube") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.y   = element_text(size = 7.5),
        legend.position = "top",
        legend.text   = element_text(size = 8),
        legend.key.height = unit(8, "pt"))

# ============================================================================
# Compose
# ============================================================================
layout <- pA / pB / pC +
  plot_layout(heights = c(1.15, 1.0, 1.15)) +
  plot_annotation(
    title = "Fig. S71 \u2014 Transcriptional response strength across cohort \u00d7 season \u00d7 stage",
    subtitle = "PCA-based distance metric isolates infection drift from normal season/stage variation.",
    theme = theme(plot.title    = element_text(family = "Arial", face = "bold",
                                               size = 12, hjust = 0.5,
                                               margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, colour = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 290)

cat(sprintf("[S71] top cell: %s %s %s = %.2f\n",
            rank_df$Season[1],
            STAGE_LAB[as.character(rank_df$Maturity_stage[1])],
            as.character(rank_df$Cohort[1]),
            rank_df$mean_resp[1]))
