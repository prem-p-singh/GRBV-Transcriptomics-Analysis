# S67 — RB(-) baseline companion to S66.
#
# Answers "what happens in the uninfected controls across the same grid?",
# providing a null comparator for S66.  RB(-) samples have no infection
# cohort dimension (all controls), so the grid collapses to:
#   stage x season   (for expression layers)
#   season           (for Brix; Harvest stage only)
#
# Three layers:
#   A. sHSP module z-score   (stage x season, RB- only)
#   B. GRBV signature score  (stage x season, RB- only)
#   C. Harvest juice Brix    (season, RB- only, single row)
#
# Output: figures/supp/S67_rbminus_baseline_cube.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S67_rbminus_baseline_cube")

STAGES   <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STAGE_LAB <- c(Pre_veraison = "Pre-v\u00e9r", Veraison = "V\u00e9r",
               Post_veraison = "Post-v\u00e9r", Harvest = "Harvest")
SEASONS  <- c("2021","2022")

# ============================================================================
# LOAD + compute sample-level scores (same recipe as S66)
# ============================================================================
meta <- read_tsv(file.path(ROOT, "sample_metadata.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
meta$Maturity_stage <- factor(meta$Maturity_stage, levels = STAGES)
meta$Sampling_Season <- as.character(meta$Sampling_Season)

counts <- fread(file.path(ROOT, "counts_matrix.txt"),
                sep = "\t", header = TRUE, data.table = FALSE)
rownames(counts) <- counts[, 1]; counts[, 1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"
keep   <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]
lib    <- colSums(counts)
lcpm   <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)
meta_e <- meta[match(colnames(lcpm), meta$Sample), ]

heat_de <- fread(file.path(ROOT, "heat_stress_DEGs_all_conditions.csv")) %>%
  as.data.frame()
shsp_genes <- intersect(unique(heat_de$Gene[grepl("sHSP", heat_de$Category)]),
                        rownames(lcpm))
sig <- read_csv(file.path(ROOT, "tier1_grbv_signature.csv"),
                show_col_types = FALSE) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FC) >= 1) %>%
  filter(LOC_id %in% rownames(lcpm))

zscore_rows <- function(M) {
  t(apply(M, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
  }))
}
shsp_z <- zscore_rows(lcpm[shsp_genes, , drop = FALSE])
meta_e$shsp_score <- colMeans(shsp_z, na.rm = TRUE)

sig_z   <- zscore_rows(lcpm[sig$LOC_id, , drop = FALSE])
sign_w  <- sign(sig$log2FC)
meta_e$sig_score <- colSums(sig_z * sign_w) / length(sign_w)

# ============================================================================
# RB(-) ONLY aggregates
# ============================================================================
meta_rbn <- meta_e %>% filter(Infection_Group == "RB(-)",
                              Sampling_Season %in% SEASONS)

cell_shsp <- meta_rbn %>%
  group_by(Sampling_Season, Maturity_stage) %>%
  summarise(value = mean(shsp_score, na.rm = TRUE), n = n(), .groups = "drop") %>%
  rename(Season = Sampling_Season)

cell_sig <- meta_rbn %>%
  group_by(Sampling_Season, Maturity_stage) %>%
  summarise(value = mean(sig_score, na.rm = TRUE), n = n(), .groups = "drop") %>%
  rename(Season = Sampling_Season)

full_grid <- expand_grid(Season = SEASONS,
                         Maturity_stage = factor(STAGES, levels = STAGES))
fill_grid <- function(df) {
  full_grid %>% left_join(df, by = c("Season","Maturity_stage"))
}
cell_shsp <- fill_grid(cell_shsp)
cell_sig  <- fill_grid(cell_sig)

# Harvest Brix, RB(-), aggregated per season
biochem <- read_csv(file.path(ROOT, "biochemistry_harvest.csv"),
                    show_col_types = FALSE) %>% as.data.frame()
brix_rbn <- biochem %>% filter(Infection_Group == "RB(-)") %>%
  mutate(Season = as.character(Sampling_Season)) %>%
  group_by(Season) %>%
  summarise(value = mean(Brix, na.rm = TRUE), n = n(), .groups = "drop")

# ============================================================================
# RENDER
# ============================================================================
heatmap_panel <- function(df, title, subtitle, fill_name, tag = "") {
  lim <- max(abs(df$value), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  df$cell_label <- ifelse(!is.na(df$value) & !is.na(df$n) & df$n > 0,
                          sprintf("%+.2f\n(n=%d)", df$value, df$n), "")
  ggplot(df, aes(x = Maturity_stage, y = Season)) +
    geom_tile(fill = "#E8E8E8", colour = "white", linewidth = 0.4) +
    geom_tile(aes(fill = value), colour = "white", linewidth = 0.4,
              data = df %>% filter(!is.na(value))) +
    geom_text(aes(label = cell_label),
              family = "Arial", size = 2.4, colour = "grey15",
              lineheight = 0.95) +
    scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                         mid = "white",
                         high = unname(GRBV_COLORS["up"]),
                         midpoint = 0, limits = c(-lim, lim),
                         na.value = NA, name = fill_name) +
    scale_x_discrete(labels = STAGE_LAB) +
    labs(title = title, subtitle = subtitle, tag = tag,
         x = NULL, y = "Season") +
    theme_grbv() +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          axis.text.x   = element_text(angle = 25, hjust = 1, size = 8),
          axis.text.y   = element_text(size = 9),
          legend.key.width  = unit(8, "pt"),
          legend.key.height = unit(14, "pt"))
}

pA <- heatmap_panel(
  cell_shsp,
  title = "A  sHSP module score (RB\u2212 controls)",
  subtitle = "mean z-score across 36 sHSP genes; baseline host thermal response",
  fill_name = "z-score"
)

pB <- heatmap_panel(
  cell_sig,
  title = "B  GRBV signature score (RB\u2212 controls)",
  subtitle = sprintf("mean sign-weighted z-score across %d GRBV DEGs in uninfected plants",
                     nrow(sig)),
  fill_name = "z-score"
)

# Harvest Brix for RB(-): a compact row of two tiles (one per season)
brix_rbn$cell_label <- sprintf("%.1f \u00b0Bx\n(n=%d)", brix_rbn$value, brix_rbn$n)
pC <- ggplot(brix_rbn, aes(x = Season, y = 1)) +
  geom_tile(aes(fill = value), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_label),
            family = "Arial", size = 2.7, colour = "grey15",
            lineheight = 0.95) +
  scale_fill_gradient(low = "#F4F2D6",
                      high = unname(GRBV_COLORS["stage_harv"]),
                      name = "\u00b0Brix") +
  scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
  labs(title = "C  Harvest juice \u00b0Brix (RB\u2212 controls)",
       subtitle = "mean per season across uninfected controls",
       x = NULL, y = "Harvest") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.text.x   = element_text(size = 9),
        legend.key.width  = unit(8, "pt"),
        legend.key.height = unit(14, "pt"))

layout <- pA / pB / pC +
  plot_layout(heights = c(1, 1, 0.5)) +
  plot_annotation(
    title = "Fig. S67 \u2014 RB(\u2212) baseline: same grid, uninfected controls",
    subtitle = "Null comparator for S66. No cohort axis: controls form a single group per season \u00d7 stage. Compare directly to S66 B, C, D.",
    theme = theme(plot.title    = element_text(family = "Arial", face = "bold",
                                               size = 12, hjust = 0.5,
                                               margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, colour = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 160)

cat(sprintf("[S67] rendered: RB\u2212 grid with %d stage\u00d7season cells, Brix: %d seasons\n",
            nrow(cell_shsp), nrow(brix_rbn)))
