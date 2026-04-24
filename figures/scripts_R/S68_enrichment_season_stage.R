# S68 — Pathway enrichment shifts across season x stage (RB+ vs RB-).
#
# Third diagnostic layer alongside S66 (expression cube) and S67 (baseline):
# how does pathway enrichment change as a function of stage and season?
#
# Panels:
#   A  Heatmap of top-30 most variable pathways x 8 conditions (NES)
#   B  Season-flip pathways: sign(NES_2021) != sign(NES_2022) at same stage
#   C  NES trajectory lines for 8 key pathway families
#
# Source: gsea_results_all_conditions.csv (3,000 pathway x condition rows,
# 375 gene sets x 8 stage x season conditions, 500 permutations).
#
# Output: figures/supp/S68_enrichment_season_stage.png + .svg

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(stringr)
})

HERE     <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR  <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S68_enrichment_season_stage")

STAGES    <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STAGE_LAB <- c(Pre_veraison = "Pre-v\u00e9r", Veraison = "V\u00e9r",
               Post_veraison = "Post-v\u00e9r", Harvest = "Harvest")

# Canonical condition order (matches the season x stage grid in S66/S67)
COND_ORDER <- c("2021_Pre_veraison","2021_Veraison",
                "2021_Post_veraison","2021_Harvest",
                "2022_Pre_veraison","2022_Veraison",
                "2022_Post_veraison","2022_Harvest")

# ============================================================================
# LOAD + preprocess
# ============================================================================
gs <- read_csv(file.path(ROOT, "gsea_results_all_conditions.csv"),
               show_col_types = FALSE) %>%
  mutate(Season = str_extract(condition, "^\\d{4}"),
         Stage  = str_remove(condition, "^\\d{4}_"),
         Stage  = factor(Stage, levels = STAGES),
         condition = factor(condition, levels = COND_ORDER))

# Short pathway labels: strip DB prefix for axis readability; keep DB as a tag
gs <- gs %>%
  mutate(db       = str_extract(gene_set, "^[A-Z]+"),
         gs_short = str_remove(gene_set, "^[A-Z]+:"),
         gs_short = ifelse(nchar(gs_short) > 34,
                           paste0(substr(gs_short, 1, 33), "\u2026"), gs_short))

# ============================================================================
# Panel A — top-30 most variable pathways heatmap
# ============================================================================
top_var <- gs %>%
  group_by(gene_set) %>%
  summarise(var_nes = max(abs(nes), na.rm = TRUE),
            any_sig = any(fdr < 0.05, na.rm = TRUE),
            .groups = "drop") %>%
  filter(any_sig) %>%
  arrange(desc(var_nes)) %>%
  slice_head(n = 30)

matA <- gs %>% filter(gene_set %in% top_var$gene_set) %>%
  select(gene_set, gs_short, db, condition, nes) %>%
  pivot_wider(names_from = condition, values_from = nes,
              id_cols = c(gene_set, gs_short, db)) %>%
  as.data.frame()

# Row cluster on the 8-condition NES profile
ord_idx <- hclust(dist(as.matrix(matA[, COND_ORDER])))$order
row_order <- matA$gs_short[ord_idx]

matA_long <- gs %>% filter(gene_set %in% top_var$gene_set) %>%
  mutate(gs_short = factor(gs_short, levels = row_order))

pA <- ggplot(matA_long, aes(x = condition, y = gs_short, fill = nes)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")),
            family = "Arial", size = 2.3, colour = "grey15",
            vjust = 0.75) +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish, name = "NES") +
  scale_x_discrete(labels = function(x) {
    xs <- str_split_fixed(x, "_", 2)
    paste0(xs[,1], "\n", STAGE_LAB[xs[,2]])
  }) +
  labs(x = NULL, y = NULL,
       title = "A  Top 30 most-variable pathways (FDR<0.05 somewhere)",
       subtitle = "* marks FDR<0.05 in that condition; rows clustered on 8-condition NES profile") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.x   = element_text(size = 7, lineheight = 1.0),
        axis.text.y   = element_text(size = 7),
        legend.key.width  = unit(8, "pt"),
        legend.key.height = unit(14, "pt"))

# ============================================================================
# Panel B — season-flip pathways
#   pairs where sign(NES_2021) differs from sign(NES_2022) at the same stage
#   AND |NES| > 1.0 in at least one season (non-trivial flip)
# ============================================================================
wide <- gs %>% select(gene_set, gs_short, Stage, Season, nes, fdr) %>%
  pivot_wider(names_from = Season, values_from = c(nes, fdr),
              id_cols = c(gene_set, gs_short, Stage))

flip <- wide %>%
  filter(!is.na(nes_2021), !is.na(nes_2022),
         sign(nes_2021) != sign(nes_2022),
         pmax(abs(nes_2021), abs(nes_2022)) >= 1.0,
         (fdr_2021 < 0.05 | fdr_2022 < 0.05)) %>%
  mutate(delta = nes_2022 - nes_2021) %>%
  arrange(desc(abs(delta))) %>%
  slice_head(n = 20)

flip_long <- flip %>%
  select(gs_short, Stage, nes_2021, nes_2022) %>%
  pivot_longer(c(nes_2021, nes_2022),
               names_to = "Season", values_to = "nes") %>%
  mutate(Season = sub("nes_", "", Season)) %>%
  mutate(label = sprintf("%s  \u00b7  %s", STAGE_LAB[as.character(Stage)], gs_short))

row_order_B <- flip %>% arrange(delta) %>%
  mutate(label = sprintf("%s  \u00b7  %s", STAGE_LAB[as.character(Stage)], gs_short)) %>%
  pull(label)
flip_long$label <- factor(flip_long$label, levels = row_order_B)

pB <- ggplot(flip_long, aes(x = nes, y = label, colour = Season)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(group = label), colour = "grey75", linewidth = 0.5) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c("2021" = unname(GRBV_COLORS["season_2021"]),
                                 "2022" = unname(GRBV_COLORS["season_2022"])),
                      name = "Season") +
  labs(x = "NES (sign flips between 2021 and 2022)", y = NULL,
       title = "B  Season-flip pathways",
       subtitle = "pathways where NES reverses sign between 2021 and 2022 at the same stage (top 20 by |\u0394NES|)") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.y   = element_text(size = 6.5),
        legend.position = "top")

# ============================================================================
# Panel C — NES trajectory for 8 key pathway families
# ============================================================================
key_patterns <- c(
  "ribosome"        = "[Rr]ibosom",
  "chaperone"       = "protein folding|unfolded protein|chaperone",
  "heat response"   = "response to heat|heat shock",
  "translation"     = "translation initiation|[Tt]ranslation",
  "photosynthesis"  = "[Pp]hotosynthesi|light harvesting|photosystem",
  "sugar metab."    = "[Gg]lycolysis|[Ss]ucrose|[Ss]tarch",
  "defense"         = "defense|immune|pathogenesis",
  "2nd metab."      = "[Ff]lavonoid|[Pp]henylpropan|[Tt]erpen"
)

traj_rows <- lapply(names(key_patterns), function(k) {
  rx <- key_patterns[[k]]
  hits <- gs %>%
    filter(str_detect(gene_set, rx)) %>%
    group_by(Season, Stage) %>%
    summarise(nes = mean(nes, na.rm = TRUE),
              n_sets = n_distinct(gene_set), .groups = "drop") %>%
    mutate(family = k)
  hits
})
traj <- bind_rows(traj_rows) %>%
  mutate(family = factor(family, levels = names(key_patterns)),
         Stage  = factor(Stage, levels = STAGES))

pC <- ggplot(traj, aes(x = Stage, y = nes,
                       colour = Season, group = Season)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("2021" = unname(GRBV_COLORS["season_2021"]),
                                 "2022" = unname(GRBV_COLORS["season_2022"])),
                      name = "Season") +
  scale_x_discrete(labels = STAGE_LAB) +
  facet_wrap(~ family, ncol = 4, scales = "free_y") +
  labs(x = NULL, y = "Mean NES across matched gene sets",
       title = "C  Pathway trajectory across stages, per season",
       subtitle = "mean NES across gene sets matched by family keyword") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.text    = element_text(face = "bold", size = 8),
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 7),
        axis.text.y   = element_text(size = 7),
        legend.position = "top",
        panel.spacing = unit(6, "pt"))

# ============================================================================
# COMPOSE
# ============================================================================
layout <- pA / pB / pC +
  plot_layout(heights = c(1.6, 1.15, 1.25)) +
  plot_annotation(
    title = "Fig. S68 \u2014 Pathway enrichment across season \u00d7 maturity stage",
    subtitle = "GSEA (NES) across 8 conditions; A top variable sets, B season-flip sets, C trajectory of 8 key pathway families.",
    theme = theme(plot.title    = element_text(family = "Arial", face = "bold",
                                               size = 12, hjust = 0.5,
                                               margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, colour = "grey25",
                                               margin = margin(b = 6)))
  )

save_grbv(layout, OUT_STEM,
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 310)

cat(sprintf("[S68] top-variable: %d pathways | flips: %d | trajectories: %d families\n",
            nrow(top_var), nrow(flip), length(key_patterns)))
