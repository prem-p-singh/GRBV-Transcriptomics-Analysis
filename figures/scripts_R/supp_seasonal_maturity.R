# Supplementary batch: seasonal-maturity series (S15, S16, S17, S19)
# refit on fix01c Q2x contrasts (per-stage x per-season RB+ vs RB-).

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(ggrepel)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results"
ANN_F   <- file.path(RESULTS, "fix11_annotations/gene_annotations.csv")
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ <- 0.05; LFC <- 1.0
stages  <- GRBV_STAGES
seasons <- c("2021","2022")

# ------ load Q2x ------
load_q2x <- function() {
  rows <- list()
  for (sn in seasons) for (st in stages) {
    f <- file.path(RESULTS, "fix01c_de_extras",
                   sprintf("Q2x_%s_%s_deseq2.csv", sn, st))
    d <- fread(f) %>% as.data.frame()
    d$Sampling_Season <- sn; d$Maturity_stage <- st
    rows[[length(rows)+1]] <- d
  }
  rbindlist(rows, fill = TRUE) %>% as.data.frame()
}
q2x <- load_q2x()
ann <- fread(ANN_F) %>% as.data.frame()
name_map <- setNames(ann$product, ann$gene_id)

# ------ shared: DEG counts ------
deg <- q2x %>%
  mutate(sig = !is.na(padj) & padj < PADJ & abs(log2FC) >= LFC,
         dir = ifelse(log2FC > 0, "up", "down")) %>%
  filter(sig)
deg_counts <- deg %>% count(Sampling_Season, Maturity_stage, dir) %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0) %>%
  mutate(total = up + down,
         Maturity_stage = factor(Maturity_stage, levels = stages))

# -----------------------------------------------------------------
# S15 — KEY inverted seasonal trajectory detail (already condensed in Main 3)
# -----------------------------------------------------------------
p_a <- ggplot(deg_counts %>%
                 pivot_longer(c(up, down), names_to = "dir", values_to = "n"),
              aes(x = Maturity_stage, y = n, fill = dir)) +
  geom_col(position = position_stack(reverse = TRUE),
           width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(data = deg_counts,
            aes(x = Maturity_stage, y = total, label = scales::comma(total)),
            inherit.aes = FALSE, vjust = -0.4, family = "Arial", size = 2.7) +
  facet_wrap(~ Sampling_Season, ncol = 2) +
  scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                               down = unname(GRBV_COLORS["down"])),
                    labels = c(up = "Up", down = "Down"), name = NULL) +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::label_comma()) +
  labs(x = NULL, y = "DEGs (padj<0.05, |log2FC|\u22651)",
       title = "S15. Inverted seasonal DEG trajectory",
       subtitle = "Per-stage \u00d7 per-season DESeq2 contrasts (RB+ vs RB\u2212);\ncool 2021 peaks late at Harvest, hot 2022 peaks early at Post-v\u00e9raison") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6), lineheight = 1.1),
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 8),
        strip.text    = element_text(face = "bold", size = 10),
        legend.position = c(0.99, 0.98),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 8))

save_grbv(p_a, file.path(OUT_DIR, "S15_seasonal_maturity_DEG_comparison"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 120)

# -----------------------------------------------------------------
# S16 — stage x season DEG directionality balloon plot
# -----------------------------------------------------------------
# One bubble per (stage x season) cell.  Size = total DEGs,
# fill = net (up - down) on a diverging RdBu scale through white at 0,
# label = total count + "% up".  Mirrors Main04 panel A's visual language.
cells <- deg %>% count(Sampling_Season, Maturity_stage, dir) %>%
  mutate(Maturity_stage = factor(Maturity_stage, levels = stages)) %>%
  tidyr::pivot_wider(names_from = dir, values_from = n, values_fill = 0) %>%
  mutate(total  = up + down,
         net    = up - down,
         pct_up = ifelse(total == 0, 0, up / total),
         label  = sprintf("%s\n(%.0f%% up)",
                          scales::comma(total), 100 * pct_up))

net_max <- max(abs(cells$net), na.rm = TRUE)

p_16 <- ggplot(cells,
               aes(x = Maturity_stage, y = Sampling_Season,
                   size = total, fill = net)) +
  geom_point(shape = 21, color = "grey25", stroke = 0.4, alpha = 0.95) +
  geom_text(aes(label = label),
            size = 2.3, family = "Arial", color = "grey15",
            vjust = -1.3, lineheight = 0.9) +
  scale_fill_gradient2(low  = unname(GRBV_COLORS["down"]),
                       mid  = "#F7F7F7",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0,
                       limits = c(-net_max, net_max),
                       name = "Net (up \u2212 down)") +
  scale_size_area(max_size = 20, name = "Total DEGs",
                  breaks = c(100, 300, 600, 1000)) +
  scale_x_discrete(labels = GRBV_STAGE_LABEL, expand = expansion(add = 0.6)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(x = NULL, y = "Sampling season",
       title = "S16. Stage \u00d7 season DEG directionality",
       subtitle = "Bubble size = total DEGs;  fill = net (up \u2212 down);  label = total and % up") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 8),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.height = unit(10, "pt"))

save_grbv(p_16, file.path(OUT_DIR, "S16_seasonal_maturity_pathway_enrichment"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 110)

# -----------------------------------------------------------------
# S17 — per-season stage landscape: DEG density and overlap heatmap
# -----------------------------------------------------------------
dir_mat <- q2x %>%
  mutate(sig = !is.na(padj) & padj < PADJ & abs(log2FC) >= LFC) %>%
  filter(sig) %>%
  transmute(gene, Sampling_Season, Maturity_stage,
            up = log2FC > 0) %>%
  distinct()

# Overlap between stages within each season
overlap_fn <- function(season) {
  d <- dir_mat %>% filter(Sampling_Season == season)
  out <- expand.grid(a = stages, b = stages) %>%
    mutate(n_shared = mapply(function(x, y) {
      ga <- d$gene[d$Maturity_stage == x]
      gb <- d$gene[d$Maturity_stage == y]
      length(intersect(ga, gb))
    }, a, b))
  out$season <- season
  out
}
ov <- bind_rows(overlap_fn("2021"), overlap_fn("2022")) %>%
  mutate(a = factor(a, levels = stages),
         b = factor(b, levels = rev(stages)))

p_17 <- ggplot(ov, aes(x = a, y = b, fill = n_shared)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = scales::comma(n_shared)),
            family = "Arial", size = 2.3,
            color = ifelse(ov$n_shared > quantile(ov$n_shared, 0.75), "white", "grey15")) +
  facet_wrap(~ season, ncol = 2) +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  scale_y_discrete(labels = GRBV_STAGE_LABEL) +
  scale_fill_gradient(low = "#fff5e1", high = unname(GRBV_COLORS["RB_pos"]),
                      name = "Shared DEGs") +
  labs(x = NULL, y = NULL,
       title = "S17. Pairwise DEG overlap between stages (per season)",
       subtitle = "Upper + lower triangles are symmetric;  diagonal = stage DEG count") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 8),
        axis.text.y   = element_text(size = 8),
        strip.text    = element_text(face = "bold", size = 10),
        legend.key.height = unit(15, "pt"))

save_grbv(p_17, file.path(OUT_DIR, "S17_seasonal_maturity_pathway_landscape"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 130)

# -----------------------------------------------------------------
# S19 — top-gene trajectories across stages, colored by season
# -----------------------------------------------------------------
# pick 10 most-recurrent DEGs (sig in >=4 stage×season cells)
rec <- q2x %>%
  mutate(sig = !is.na(padj) & padj < PADJ & abs(log2FC) >= LFC) %>%
  filter(sig) %>%
  count(gene, name = "n_cond") %>%
  arrange(desc(n_cond)) %>% head(12)

top_prof <- q2x %>% filter(gene %in% rec$gene) %>%
  mutate(Maturity_stage = factor(Maturity_stage, levels = stages),
         Sampling_Season = factor(Sampling_Season),
         product = ifelse(gene %in% names(name_map),
                          substr(name_map[gene], 1, 26), ""),
         panel_label = paste0(gene, "\n", product))

p_19 <- ggplot(top_prof,
               aes(x = Maturity_stage, y = log2FC,
                   color = Sampling_Season, group = Sampling_Season)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  facet_wrap(~ panel_label, ncol = 4, scales = "free_y") +
  scale_color_manual(values = GRBV_SEASON_COLOR, name = "Season") +
  scale_x_discrete(labels = c(Pre_veraison = "Pre", Veraison = "V\u00e9r",
                              Post_veraison = "Post", Harvest = "H")) +
  labs(x = NULL, y = "log2 fold change",
       title = "S19. Fold-change trajectories of the most recurrent DEGs",
       subtitle = "Top 12 genes significant in \u22654 of 8 stage \u00d7 season conditions") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        strip.text    = element_text(size = 7.5, face = "bold",
                                     lineheight = 1.0),
        axis.text.x   = element_text(size = 7),
        legend.position = "top")

save_grbv(p_19, file.path(OUT_DIR, "S19_seasonal_maturity_pathway_detail"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 200)

cat("[supp S15-S19] done.\n")
