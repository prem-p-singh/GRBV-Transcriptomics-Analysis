# Main Composite 2 — Differential expression landscape across Q1-Q6
# under covariate-adjusted DESeq2 (fix01 + fix01c).
#
# Panels:
#   A  DEG counts per comparison (UP vs DOWN), all 6 research questions
#   B  Volcano of Q6 pooled RB+ vs RB-
#   C  Per-cohort-year DEG counts (Q1b biphasic pattern)
#   D  Per-stage DEG counts (Q2; 4 stages, both seasons combined)
#
# Inputs:
#   farm_pipeline/results/fix01_de_redo/
#   farm_pipeline/results/fix01c_de_extras/
#   farm_pipeline/results/fix11_annotations/gene_annotations.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(ggrepel); library(scales); library(data.table)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results"
ANN_F    <- file.path(ROOT, "fix11_annotations/gene_annotations.csv")
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "Main02_de_landscape")

PADJ <- 0.05; LFC <- 1.0

read_de <- function(path) {
  d <- fread(path); d <- as.data.frame(d)
  d$path <- basename(path)
  d
}

de_files <- c(
  list.files(file.path(ROOT, "fix01_de_redo"),
             pattern = "_deseq2\\.csv$", full.names = TRUE),
  list.files(file.path(ROOT, "fix01c_de_extras"),
             pattern = "_deseq2\\.csv$", full.names = TRUE)
)
all_de <- rbindlist(lapply(de_files, function(f) {
  d <- fread(f); d$path <- basename(f); d
}), fill = TRUE) %>% as.data.frame()

count_degs <- function(df, name) {
  d <- df[df$comparison == name, ]
  ok <- !is.na(d$padj) & d$padj < PADJ & abs(d$log2FC) >= LFC
  c(up = sum(ok & d$log2FC >  0),
    down = sum(ok & d$log2FC <  0))
}

contrasts_in_order <- list(
  list(label = "Q1a recent (\u22642 yr)",       comp = "Q1a_recent_le2yr_vs_neg",     group = "Q1"),
  list(label = "Q1a mid (3-4 yr)",              comp = "Q1a_mid_3to4yr_vs_neg",       group = "Q1"),
  list(label = "Q1a old (\u22655 yr)",          comp = "Q1a_old_ge5yr_vs_neg",        group = "Q1"),
  list(label = "Q2 Pre-véraison",               comp = "Q2_Pre_veraison",             group = "Q2"),
  list(label = "Q2 Véraison",                   comp = "Q2_Veraison",                 group = "Q2"),
  list(label = "Q2 Post-véraison",              comp = "Q2_Post_veraison",            group = "Q2"),
  list(label = "Q2 Harvest",                    comp = "Q2_Harvest",                  group = "Q2"),
  list(label = "Q3 Season (RB+)",               comp = "Q3_seasonRBpos",              group = "Q3"),
  list(label = "Q4 Season \u00d7 Infection",    comp = "Q4d_seasonXinfection_dir",    group = "Q4"),
  list(label = "Q5 Season (RB\u2212)",          comp = "Q5_seasonRBneg",              group = "Q5"),
  list(label = "Q6 Pooled RB+ vs RB\u2212",     comp = "Q6_pooled",                   group = "Q6")
)

panel_A_df <- do.call(rbind, lapply(contrasts_in_order, function(x) {
  cnt <- count_degs(all_de, x$comp)
  data.frame(label = x$label, group = x$group,
             up = unname(cnt["up"]), down = unname(cnt["down"]))
}))
panel_A_long <- panel_A_df %>%
  mutate(label = factor(label, levels = rev(panel_A_df$label))) %>%
  pivot_longer(c(up, down), names_to = "dir", values_to = "n") %>%
  mutate(n_signed = ifelse(dir == "up", n, -n))

q_color <- c(Q1 = "#3498db", Q2 = "#e67e22", Q3 = "#16a085",
             Q4 = "#8e44ad", Q5 = "#7f8c8d", Q6 = "#c0392b")
panel_A_df$grpcolor <- q_color[panel_A_df$group]

# -----------------------------------------------------------------
# Panel A — diverging bar of UP vs DOWN per contrast
# -----------------------------------------------------------------
total <- panel_A_long %>% group_by(label) %>%
  summarise(total = sum(n),
            up_end = sum(n[dir == "up"]),
            down_end = -sum(n[dir == "down"]),
            .groups = "drop")

x_extent <- max(abs(total$up_end), abs(total$down_end))

pA <- ggplot(panel_A_long, aes(x = n_signed, y = label, fill = dir)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  # UP counts to the right end of positive bars
  geom_text(data = total %>% filter(up_end > 0),
            aes(x = up_end, y = label,
                label = scales::comma(up_end)),
            inherit.aes = FALSE,
            family = "Arial", size = 2.5, hjust = -0.15,
            color = "grey15") +
  # DOWN counts to the left end of negative bars
  geom_text(data = total %>% filter(down_end < 0),
            aes(x = down_end, y = label,
                label = scales::comma(abs(down_end))),
            inherit.aes = FALSE,
            family = "Arial", size = 2.5, hjust = 1.15,
            color = "grey15") +
  scale_x_continuous(limits = c(-x_extent * 1.25, x_extent * 1.25),
                     labels = function(x) scales::comma(abs(x)),
                     name = "DEGs (DOWN \u2190 0 \u2192 UP)") +
  scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                               down = unname(GRBV_COLORS["down"])),
                    labels = c(up = "Up", down = "Down"),
                    name = NULL) +
  labs(y = NULL, tag = "A",
       title = "DEG counts per comparison") +
  theme_grbv() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7.5),
        legend.text = element_text(size = 8))

# -----------------------------------------------------------------
# Panel B — Q6 pooled volcano
# -----------------------------------------------------------------
ann <- fread(ANN_F) %>% as.data.frame()
name_map <- setNames(ann$product, ann$gene_id)

q6 <- all_de[all_de$comparison == "Q6_pooled", ]
q6 <- q6[!is.na(q6$padj) & !is.na(q6$log2FC), ]
q6$status <- with(q6,
                  ifelse(padj < PADJ & log2FC >=  LFC, "Up",
                  ifelse(padj < PADJ & log2FC <= -LFC, "Down", "ns")))
q6$negLog <- -log10(pmax(q6$padj, 1e-300))
top_lab <- q6 %>% filter(status != "ns") %>%
  arrange(padj) %>% slice_head(n = 12) %>%
  mutate(short = ifelse(gene %in% names(name_map),
                        substr(name_map[gene], 1, 24), gene),
         label = ifelse(grepl("^LOC", gene), gene, gene))

pB <- ggplot(q6, aes(x = log2FC, y = negLog, color = status)) +
  geom_point(alpha = 0.55, size = 0.7) +
  scale_color_manual(values = c(Up   = unname(GRBV_COLORS["up"]),
                                Down = unname(GRBV_COLORS["down"]),
                                ns   = unname(GRBV_COLORS["neutral"])),
                     name = NULL,
                     labels = function(x) {
                       n <- sapply(x, function(s) sum(q6$status == s))
                       sprintf("%s (%s)", x, scales::comma(n))
                     }) +
  geom_text_repel(data = top_lab, aes(label = label),
                  family = "Arial", size = 2.3, color = "black",
                  box.padding = 0.3, max.overlaps = 14,
                  segment.size = 0.2, segment.color = "grey45") +
  geom_vline(xintercept = c(-LFC, LFC), color = "grey60",
             linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(PADJ), color = "grey60",
             linetype = "dashed", linewidth = 0.3) +
  labs(x = "log2 fold change (RB+ / RB-)",
       y = expression(-log[10]*"(padj)"),
       tag = "B",
       title = "Q6 pooled volcano (DESeq2)") +
  theme_grbv() +
  theme(plot.title = element_text(size = 9.5, face = "bold"),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.text  = element_text(size = 7),
        legend.key.height = unit(7, "pt"))

# -----------------------------------------------------------------
# Panel C — per-cohort-year DEGs (Q1b)
# -----------------------------------------------------------------
cohorts <- c("2016","2019","2020","2021","2022")
panel_C_df <- do.call(rbind, lapply(cohorts, function(yr) {
  comp <- sprintf("Q1b_cohort%s_vs_neg", yr)
  cnt <- count_degs(all_de, comp)
  data.frame(cohort = yr, up = unname(cnt["up"]), down = unname(cnt["down"]))
})) %>%
  mutate(total = up + down,
         years_since = as.integer(2022) - as.integer(cohort))
panel_C_long <- panel_C_df %>%
  pivot_longer(c(up, down), names_to = "dir", values_to = "n")

pC <- ggplot(panel_C_long, aes(x = cohort, y = n, fill = dir)) +
  geom_col(position = position_stack(reverse = TRUE),
           width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(data = panel_C_df, aes(x = cohort, y = total, label = scales::comma(total)),
            inherit.aes = FALSE, family = "Arial", size = 2.6, vjust = -0.3) +
  scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                               down = unname(GRBV_COLORS["down"])),
                    labels = c(up = "Up", down = "Down"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Infection-onset cohort year",
       y = "DEGs (padj<0.05, |log2FC|\u22651)",
       tag = "C",
       title = "Per-cohort DE response (Q1b)") +
  theme_grbv() +
  theme(plot.title = element_text(size = 9.5, face = "bold"),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.text  = element_text(size = 7),
        legend.key.height = unit(7, "pt"))

# -----------------------------------------------------------------
# Panel D — per-stage DEGs (Q2)
# -----------------------------------------------------------------
stages <- GRBV_STAGES
panel_D_df <- do.call(rbind, lapply(stages, function(st) {
  comp <- sprintf("Q2_%s", st)
  cnt <- count_degs(all_de, comp)
  data.frame(stage = st, up = unname(cnt["up"]), down = unname(cnt["down"]))
})) %>%
  mutate(stage = factor(stage, levels = stages),
         total = up + down)
panel_D_long <- panel_D_df %>%
  pivot_longer(c(up, down), names_to = "dir", values_to = "n")

pD <- ggplot(panel_D_long, aes(x = stage, y = n, fill = dir)) +
  geom_col(position = position_stack(reverse = TRUE),
           width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(data = panel_D_df, aes(x = stage, y = total, label = scales::comma(total)),
            inherit.aes = FALSE, family = "Arial", size = 2.6, vjust = -0.3) +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                               down = unname(GRBV_COLORS["down"])),
                    labels = c(up = "Up", down = "Down"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "DEGs",
       tag = "D",
       title = "Per-stage DE response (Q2)") +
  theme_grbv() +
  theme(plot.title = element_text(size = 9.5, face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
        legend.position = "none")

# -----------------------------------------------------------------
# Composite layout
# -----------------------------------------------------------------
layout <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Differential expression landscape (DESeq2 reanalysis)",
    subtitle = "All contrasts refit with covariate-adjusted DESeq2;  significance: padj < 0.05 and |log2FC| \u2265 1",
    theme = theme(plot.title    = element_text(family = "Arial", face = "bold",
                                               size = 13, hjust = 0.5,
                                               margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 8.5,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 200)

cat(sprintf("[main02] panel A rows = %d\n", nrow(panel_A_df)))
cat("[main02] Q6 pooled UP=", sum(q6$status == "Up"),
    "  DOWN=", sum(q6$status == "Down"), "\n")
