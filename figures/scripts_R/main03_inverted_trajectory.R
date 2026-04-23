# Main Composite 3 — Inverted seasonal trajectory + biphasic viral titer.
#
# Panels:
#   A  Stage x season DEG counts (Q2x); 2021 peak at Harvest, 2022 at Post-vér.
#   B  Cross-season log2FC correlation per stage (scatter, r in title).
#   C  Mean viral titer per stage x season; biphasic temperature effect.
#
# Inputs:
#   farm_pipeline/results/de_extras/Q2x_*_deseq2.csv
#   /Users/prempratapsingh/Data_Analysis/Virus titer.txt

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT      <- "/Users/prempratapsingh/Data_Analysis"
RESULTS   <- file.path(ROOT, "farm_pipeline/results/de_extras")
TITER_F   <- file.path(ROOT, "Virus titer.txt")
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM  <- file.path(OUT_DIR, "Main03_inverted_trajectory")

PADJ <- 0.05; LFC <- 1.0
seasons <- c("2021","2022"); stages <- GRBV_STAGES

# -----------------------------------------------------------------
# Load all Q2x_*_deseq2 results into one data frame
# -----------------------------------------------------------------
load_q2x <- function() {
  rows <- list()
  for (sn in seasons) for (st in stages) {
    f <- file.path(RESULTS, sprintf("Q2x_%s_%s_deseq2.csv", sn, st))
    d <- fread(f) %>% as.data.frame()
    d$Sampling_Season <- sn; d$Maturity_stage <- st
    rows[[length(rows)+1]] <- d
  }
  rbindlist(rows, fill = TRUE) %>% as.data.frame()
}
q2x <- load_q2x()

deg_counts <- q2x %>%
  mutate(sig = !is.na(padj) & padj < PADJ & abs(log2FC) >= LFC,
         dir = ifelse(log2FC > 0, "up", "down")) %>%
  filter(sig) %>%
  count(Sampling_Season, Maturity_stage, dir) %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0) %>%
  mutate(total = up + down,
         Maturity_stage = factor(Maturity_stage, levels = stages))

# -----------------------------------------------------------------
# Panel A — DEG counts per stage x season (inverted trajectory)
# -----------------------------------------------------------------
deg_long <- deg_counts %>% select(Sampling_Season, Maturity_stage, up, down) %>%
  pivot_longer(c(up, down), names_to = "dir", values_to = "n")
totals <- deg_counts %>% select(Sampling_Season, Maturity_stage, total)

pA <- ggplot(deg_long, aes(x = Maturity_stage, y = n, fill = dir)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7,
           color = "white", linewidth = 0.2) +
  geom_text(data = totals,
            aes(x = Maturity_stage, y = total, label = scales::comma(total)),
            inherit.aes = FALSE, vjust = -0.4,
            family = "Arial", size = 2.7) +
  facet_wrap(~ Sampling_Season, ncol = 2, strip.position = "top") +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                               down = unname(GRBV_COLORS["down"])),
                    labels = c(up = "Up", down = "Down"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = scales::label_comma()) +
  labs(x = NULL, y = "DEGs (padj<0.05, |log2FC|\u22651)",
       tag = "A",
       title = "Inverted seasonal DEG trajectory") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
        strip.text  = element_text(size = 9, face = "bold"),
        legend.position = c(0.99, 0.98),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 7.5),
        legend.key.height = unit(7, "pt"))

# -----------------------------------------------------------------
# Panel B — cross-season log2FC correlation per stage
# -----------------------------------------------------------------
fc_wide <- q2x %>%
  select(gene, Sampling_Season, Maturity_stage, log2FC) %>%
  pivot_wider(names_from = Sampling_Season, values_from = log2FC,
              names_prefix = "fc_") %>%
  filter(!is.na(fc_2021), !is.na(fc_2022))

cor_per_stage <- fc_wide %>%
  group_by(Maturity_stage) %>%
  summarise(r = cor(fc_2021, fc_2022, use = "pairwise.complete.obs"),
            n = n(), .groups = "drop") %>%
  mutate(Maturity_stage = factor(Maturity_stage, levels = stages),
         label = sprintf("%s\nr = %.2f  (n = %s)",
                         GRBV_STAGE_LABEL[as.character(Maturity_stage)],
                         r, scales::comma(n)))

fc_wide_lab <- fc_wide %>%
  mutate(Maturity_stage = factor(Maturity_stage, levels = stages))
fc_wide_lab$Maturity_stage <- factor(fc_wide_lab$Maturity_stage,
                                     levels = cor_per_stage$Maturity_stage,
                                     labels = cor_per_stage$label)

pB <- ggplot(fc_wide_lab, aes(x = fc_2021, y = fc_2022)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "grey50",
              linetype = "dashed", linewidth = 0.3) +
  geom_point(alpha = 0.30, size = 0.45,
             color = unname(GRBV_COLORS["RB_pos"])) +
  facet_wrap(~ Maturity_stage, ncol = 4) +
  coord_fixed(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(x = "log2FC (2021)", y = "log2FC (2022)",
       tag = "B",
       title = "Cross-season fold-change agreement") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        strip.text  = element_text(size = 8, face = "bold", lineheight = 1.0))

# -----------------------------------------------------------------
# Panel C — mean viral titer per stage x season
# -----------------------------------------------------------------
titer_raw <- read.delim(TITER_F, header = TRUE, sep = "\t",
                        check.names = FALSE)
# Long-format
titer_long <- titer_raw %>%
  rename(Pre_veraison = `Pre_véraison`, Veraison = `Véraison`,
         Post_veraison = `Post_véraison`, Harvest = Harvest) %>%
  pivot_longer(c(Pre_veraison, Veraison, Post_veraison, Harvest),
               names_to = "Maturity_stage", values_to = "titer") %>%
  filter(!is.na(titer)) %>%
  mutate(Maturity_stage = factor(Maturity_stage, levels = stages),
         Sampling_Season = factor(Season))

titer_summary <- titer_long %>%
  group_by(Sampling_Season, Maturity_stage) %>%
  summarise(mean = mean(titer), sd = sd(titer), n = n(),
            se = sd / sqrt(n), .groups = "drop")

pC <- ggplot(titer_summary, aes(x = Maturity_stage, y = mean / 1e6,
                                 color = Sampling_Season,
                                 group = Sampling_Season,
                                 shape = Sampling_Season)) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = (mean - se) / 1e6, ymax = (mean + se) / 1e6),
                width = 0.10, linewidth = 0.5) +
  geom_point(size = 2.6, fill = "white", stroke = 1.2) +
  scale_color_manual(values = GRBV_SEASON_COLOR, name = "Season") +
  scale_shape_manual(values = c("2021" = 21, "2022" = 24), name = "Season") +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  labs(x = NULL,
       y = "Mean viral titer (\u00d7 10\u2076 copies)",
       tag = "C",
       title = "Biphasic viral titer dynamics") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5),
        legend.position = c(0.99, 0.98),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 7.5),
        legend.key.height = unit(8, "pt"))

# -----------------------------------------------------------------
# Composite layout
# -----------------------------------------------------------------
layout <- pA / pB / pC +
  plot_layout(heights = c(1.0, 1.05, 0.95)) +
  plot_annotation(
    title = "Inverted seasonal trajectory: cool 2021 peaks late, hot 2022 peaks early",
    subtitle = "Per-stage \u00d7 per-season DESeq2 contrasts (RB+ vs RB\u2212);  raw titer averaged across infection-onset cohorts",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 13, hjust = 0.5,
                                            margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 8.5,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 215)

cat("[main03] DEG counts per stage x season:\n"); print(deg_counts)
cat("\n[main03] cross-season FC correlation per stage:\n"); print(cor_per_stage)
