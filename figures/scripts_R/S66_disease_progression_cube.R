# S66 — Disease progression cube: cohort x stage x season, 3 layers
#
# Answers "how does disease change across year, season, and maturity stage?"
# by assembling a 3D atlas where each of three biological layers is rendered
# as a cohort x stage heatmap, faceted by sampling season.
#
# Layers:
#   A. Viral titer        : mean per cell on log10 scale (Virus titer.txt)
#   B. sHSP module score  : mean z-score across sHSP genes, RB(+) samples
#   C. GRBV signature     : mean sign-weighted z-score across significant
#                           GRBV DEGs (tier1_grbv_signature.csv), RB(+) samples
#
# Empty cells (no samples) drawn as hatched grey so readers see data coverage.
#
# Output: figures/supp/S66_disease_progression_cube.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S66_disease_progression_cube")

STAGES  <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STAGE_LAB <- c(Pre_veraison = "Pre-v\u00e9r", Veraison = "V\u00e9r",
               Post_veraison = "Post-v\u00e9r", Harvest = "Harvest")
COHORTS <- c("2016","2019","2020","2021","2022")
SEASONS <- c("2021","2022")

# ============================================================================
# LOAD
# ============================================================================
meta <- read_tsv(file.path(ROOT, "sample_metadata.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
meta$Maturity_stage <- factor(meta$Maturity_stage, levels = STAGES)
meta$Sampling_Season <- as.character(meta$Sampling_Season)
# Years_Since_Infection stores the cohort ONSET year for RB(+); "0" for RB(-)
meta$Cohort <- as.character(meta$Years_Since_Infection)

titer_raw <- read_tsv(file.path(ROOT, "Virus titer.txt"),
                      show_col_types = FALSE) %>% as.data.frame()
titer_long <- titer_raw %>%
  rename(Pre_veraison = `Pre_véraison`, Veraison = `Véraison`,
         Post_veraison = `Post_véraison`) %>%
  pivot_longer(cols = c(Pre_veraison, Veraison, Post_veraison, Harvest),
               names_to = "Maturity_stage", values_to = "titer") %>%
  filter(!is.na(titer)) %>%
  mutate(Season = as.character(Season),
         Cohort = as.character(Year),
         Maturity_stage = factor(Maturity_stage, levels = STAGES))

# Counts matrix (for expression-derived layers)
counts <- fread(file.path(ROOT, "counts_matrix.txt"),
                sep = "\t", header = TRUE, data.table = FALSE)
rownames(counts) <- counts[, 1]; counts[, 1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"
keep   <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]
lib    <- colSums(counts)
lcpm   <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)
meta_e <- meta[match(colnames(lcpm), meta$Sample), ]

# sHSP gene list (from Main04's heat gene panel)
heat_de <- fread(file.path(ROOT, "heat_stress_DEGs_all_conditions.csv")) %>%
  as.data.frame()
shsp_genes <- intersect(unique(heat_de$Gene[grepl("sHSP", heat_de$Category)]),
                        rownames(lcpm))

# GRBV signature (Tier 1 table): significant DEGs at padj<0.05, |log2FC|>=1
sig <- read_csv(file.path(ROOT, "tier1_grbv_signature.csv"),
                show_col_types = FALSE) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FC) >= 1) %>%
  filter(LOC_id %in% rownames(lcpm))
cat(sprintf("[S66] sHSP genes: %d | GRBV signature size: %d\n",
            length(shsp_genes), nrow(sig)))

# ============================================================================
# PER-SAMPLE SCORES
# ============================================================================
zscore_rows <- function(M) {
  t(apply(M, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
  }))
}

# sHSP module z (mean z-score across sHSP genes, per sample)
shsp_z <- zscore_rows(lcpm[shsp_genes, , drop = FALSE])
meta_e$shsp_score <- colMeans(shsp_z, na.rm = TRUE)

# Signature score: sign-weighted mean z (up genes +1, down genes -1)
sig_z <- zscore_rows(lcpm[sig$LOC_id, , drop = FALSE])
sign_w <- sign(sig$log2FC)
meta_e$sig_score <- colSums(sig_z * sign_w) / length(sign_w)

# ============================================================================
# AGGREGATE PER CELL (cohort x season x stage)
# ============================================================================
#  Layer A (titer): mean across bioreps in Virus titer table; RB(+) by design
cell_titer <- titer_long %>%
  filter(Cohort %in% COHORTS, Season %in% SEASONS) %>%
  group_by(Cohort, Season, Maturity_stage) %>%
  summarise(value = mean(titer, na.rm = TRUE),
            n = n(), .groups = "drop")

# Layer B + C: per-sample scores aggregated to cell level (RB+ only)
meta_rbp <- meta_e %>% filter(Infection_Group == "RB(+)",
                              Cohort %in% COHORTS,
                              Sampling_Season %in% SEASONS)

cell_shsp <- meta_rbp %>%
  group_by(Cohort, Sampling_Season, Maturity_stage) %>%
  summarise(value = mean(shsp_score, na.rm = TRUE), n = n(), .groups = "drop") %>%
  rename(Season = Sampling_Season)

cell_sig <- meta_rbp %>%
  group_by(Cohort, Sampling_Season, Maturity_stage) %>%
  summarise(value = mean(sig_score, na.rm = TRUE), n = n(), .groups = "drop") %>%
  rename(Season = Sampling_Season)

# Layer D: harvest juice Brix (RB+ only) aggregated by cohort x season
biochem <- read_csv(file.path(ROOT, "biochemistry_harvest.csv"),
                    show_col_types = FALSE) %>% as.data.frame()
cell_brix <- biochem %>%
  filter(Infection_Group == "RB(+)") %>%
  mutate(Cohort = as.character(Infection_Year),
         Season = as.character(Sampling_Season)) %>%
  filter(Cohort %in% COHORTS, Season %in% SEASONS) %>%
  group_by(Cohort, Season) %>%
  summarise(value = mean(Brix, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(Maturity_stage = factor("Harvest", levels = STAGES))

# Expand to the full cohort x season x stage grid; missing cells become NA
full_grid <- expand_grid(Cohort = COHORTS, Season = SEASONS,
                         Maturity_stage = factor(STAGES, levels = STAGES))

fill_grid <- function(df) {
  full_grid %>%
    left_join(df, by = c("Cohort","Season","Maturity_stage")) %>%
    mutate(Cohort = factor(Cohort, levels = rev(COHORTS)))
}
cell_titer <- fill_grid(cell_titer)
cell_shsp  <- fill_grid(cell_shsp)
cell_sig   <- fill_grid(cell_sig)

# Brix has no stage dim, so build a cohort x season grid only
brix_grid <- expand_grid(Cohort = COHORTS, Season = SEASONS) %>%
  left_join(cell_brix %>% select(-Maturity_stage),
            by = c("Cohort","Season")) %>%
  mutate(Cohort = factor(Cohort, levels = rev(COHORTS)))

# ============================================================================
# RENDER
# ============================================================================
base_heatmap <- function(df, title, subtitle, fill_name, fill_trans = "identity",
                         palette = "divergent", tag = "") {
  df$value_display <- df$value
  if (fill_trans == "log10") {
    df$value_display <- log10(pmax(df$value, 1))
  }

  # Pre-build per-cell label text as a plain character vector (avoids
  # ifelse-scalar-vs-vector bug where all cells end up with the same text).
  si_fmt <- scales::label_number(scale_cut = scales::cut_short_scale(),
                                 accuracy = 1)
  df$cell_label <- ""
  have <- !is.na(df$n) & df$n > 0 & !is.na(df$value)
  if (fill_trans == "log10") {
    df$cell_label[have] <- si_fmt(df$value[have])
  } else {
    df$cell_label[have] <- sprintf("%+.2f", df$value[have])
  }

  if (palette == "divergent") {
    lim <- max(abs(df$value_display), na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) lim <- 1
    fill_layer <- scale_fill_gradient2(
      low = unname(GRBV_COLORS["down"]),
      mid = "white",
      high = unname(GRBV_COLORS["up"]),
      midpoint = 0, limits = c(-lim, lim), na.value = NA,
      name = fill_name
    )
  } else {                                                # sequential_warm
    fill_layer <- scale_fill_gradient(
      low = "#F4E4CE",
      high = unname(GRBV_COLORS["warm"]),
      na.value = NA, name = fill_name
    )
  }

  ggplot(df, aes(x = Maturity_stage, y = Cohort)) +
    # Grey under every cell so NA cells read as "no data"
    geom_tile(fill = "#E8E8E8", colour = "white", linewidth = 0.4) +
    geom_tile(aes(fill = value_display), colour = "white", linewidth = 0.4,
              data = df %>% filter(!is.na(value_display))) +
    geom_text(aes(label = cell_label),
              family = "Arial", size = 2.1, colour = "grey15") +
    fill_layer +
    facet_wrap(~ Season, ncol = 2) +
    scale_x_discrete(labels = STAGE_LAB) +
    scale_y_discrete() +
    labs(title = title, subtitle = subtitle, tag = tag,
         x = NULL, y = "Infection cohort") +
    theme_grbv() +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          strip.text    = element_text(face = "bold", size = 9),
          axis.text.x   = element_text(angle = 25, hjust = 1, size = 8),
          legend.key.width  = unit(8, "pt"),
          legend.key.height = unit(14, "pt"),
          panel.spacing.x   = unit(6, "pt"))
}

pA <- base_heatmap(
  cell_titer,
  title = "A  Viral titer (mean across bioreps)",
  subtitle = "log\u2081\u2080(copies); blanks = no samples in that cell",
  fill_name = "log\u2081\u2080(titer)",
  fill_trans = "log10",
  palette = "sequential_warm"
)

pB <- base_heatmap(
  cell_shsp,
  title = "B  sHSP module score (host thermal defence)",
  subtitle = "mean z-score across 36 sHSP genes; negative = suppressed",
  fill_name = "z-score",
  palette = "divergent"
)

pC <- base_heatmap(
  cell_sig,
  title = "C  GRBV infection signature score (host reprogramming)",
  subtitle = sprintf("mean sign-weighted z-score across %d GRBV DEGs (padj<0.05, |log\u2082FC|\u22651)",
                     nrow(sig)),
  fill_name = "z-score",
  palette = "divergent"
)

# --------- Panel D: Harvest juice Brix (cohort x season, single cell row) ----
brix_grid$cell_label <- ifelse(is.na(brix_grid$value), "",
                               sprintf("%.1f \u00b0Bx", brix_grid$value))
brix_lim <- range(brix_grid$value, na.rm = TRUE)
pD <- ggplot(brix_grid, aes(x = 1, y = Cohort)) +
  geom_tile(fill = "#E8E8E8", colour = "white", linewidth = 0.4) +
  geom_tile(aes(fill = value), colour = "white", linewidth = 0.4,
            data = brix_grid %>% filter(!is.na(value))) +
  geom_text(aes(label = cell_label),
            family = "Arial", size = 2.4, colour = "grey15") +
  scale_fill_gradient(low = "#F4F2D6",
                      high = unname(GRBV_COLORS["stage_harv"]),
                      name = "\u00b0Brix", na.value = NA,
                      limits = brix_lim) +
  facet_wrap(~ Season, ncol = 2) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  labs(x = "Harvest", y = "Infection cohort",
       title = "D  Harvest juice \u00b0Brix (RB+ samples)",
       subtitle = "mean per cohort \u00d7 season; GRBV lowers \u00b0Brix more in hot 2022 than cool 2021") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.text    = element_text(face = "bold", size = 9),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        legend.key.width  = unit(8, "pt"),
        legend.key.height = unit(14, "pt"),
        panel.spacing.x   = unit(6, "pt"))

layout <- pA / pB / pC / pD +
  plot_layout(heights = c(1.05, 1.05, 1.05, 0.65)) +
  plot_annotation(
    title = "Fig. S66 \u2014 Disease progression across cohort \u00d7 season \u00d7 maturity stage (RB+)",
    subtitle = "Four biological layers: viral load (A), host thermal defence (B), host transcriptome signature (C), harvest juice \u00b0Brix (D).",
    theme = theme(plot.title    = element_text(family = "Arial", face = "bold",
                                               size = 12, hjust = 0.5,
                                               margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, colour = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 295)

cat(sprintf("[S66] rendered, cells: %d per layer; cohorts %s\n",
            nrow(cell_titer), paste(COHORTS, collapse = ", ")))
