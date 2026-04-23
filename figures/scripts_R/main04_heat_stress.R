# Main Composite 4 — Heat-stress disruption + sHSP suppression (REDESIGN).
#
# Redesign (rgraphmaker): single-column 3-panel stack with a teal/rose
# palette local to this figure.
#
# Panels:
#   A  Balloon plot: 9 heat-gene categories x 8 stage/season conditions.
#      Size encodes total DEGs in that cell (|up|+|down|), color encodes
#      net direction (up - down).
#   B  HSF vs sHSP module z-score scatter FACETED by maturity stage
#      (Pre-v, Vér, Post-v, Harvest). Per-infection-group linear fits,
#      per-facet r values annotated.
#   C  Expression heatmap of the top-20 recurrent heat-DEGs
#      (n_conditions >= 4): rows = genes, columns = 8 conditions,
#      tile fill = log2FC. Rows ordered by hierarchical clustering.
#
# Inputs:
#   heat_stress_DEGs_all_conditions.csv
#   heat_DEGs_recurrence.csv
#   counts_matrix.txt + sample_metadata.txt (for HSF/sHSP module scores)
#
# Output: figures/main/Main04_heat_stress.{png,pdf}
# Engine: ggplot2 + patchwork; rendered via ragg (PNG 600 dpi) + cairo_pdf.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(ggrepel)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM  <- file.path(OUT_DIR, "Main04_heat_stress")

PADJ   <- 0.05; LFC <- 1.0
stages <- GRBV_STAGES

# -----------------------------------------------------------------
# Use canonical semantic palette (no local override).
# Direction = red/blue (RdBu); Infection = vermilion/sky (Okabe-Ito).
# -----------------------------------------------------------------
PAL <- list(
  down    = unname(GRBV_COLORS["down"]),
  up      = unname(GRBV_COLORS["up"]),
  neutral = unname(GRBV_COLORS["neutral"]),
  RB_pos  = unname(GRBV_COLORS["RB_pos"]),
  RB_neg  = unname(GRBV_COLORS["RB_neg"])
)

# -----------------------------------------------------------------
# Load
# -----------------------------------------------------------------
heat_de <- fread(file.path(DATA_ROOT, "heat_stress_DEGs_all_conditions.csv")) %>%
  as.data.frame() %>%
  mutate(Stage  = factor(Stage, levels = stages),
         Season = factor(Season),
         dir    = ifelse(log2FC > 0, "up", "down"),
         is_DEG = is_DEG == "True" | is_DEG == TRUE)

heat_rec <- fread(file.path(DATA_ROOT, "heat_DEGs_recurrence.csv")) %>%
  as.data.frame()

# Category ordering by total distinct genes per category (most-recurrent first)
cat_order <- heat_de %>% group_by(Category) %>%
  summarise(n = n_distinct(Gene), .groups = "drop") %>%
  arrange(desc(n)) %>% pull(Category)

# ==================================================================
# Panel A — BALLOON PLOT (9 categories x 8 conditions)
# ==================================================================
cells <- heat_de %>% filter(is_DEG) %>%
  count(Category, Stage, Season, dir) %>%
  group_by(Category, Stage, Season) %>%
  summarise(up    = sum(n[dir == "up"]),
            down  = sum(n[dir == "down"]),
            .groups = "drop") %>%
  mutate(total = up + down,
         net   = up - down) %>%
  complete(Category = cat_order,
           Stage    = stages,
           Season   = c("2021","2022"),
           fill     = list(up = 0, down = 0, total = 0, net = 0)) %>%
  mutate(Category = factor(Category, levels = rev(cat_order)),
         Stage    = factor(Stage, levels = stages),
         Season   = factor(Season, levels = c("2021","2022")),
         cond     = paste(Season, GRBV_STAGE_LABEL[as.character(Stage)],
                          sep = " "),
         cond     = factor(cond,
                           levels = unique(cond[order(Season, Stage)])))

pA <- ggplot(cells, aes(x = cond, y = Category)) +
  geom_point(data = subset(cells, total > 0),
             aes(size = total, fill = net),
             shape = 21, color = "grey25", stroke = 0.4) +
  geom_point(data = subset(cells, total == 0),
             shape = 4, size = 1.0, color = "grey75") +
  scale_fill_gradient2(low = PAL$down, mid = "white", high = PAL$up,
                       midpoint = 0, name = "Net (up \u2212 down)",
                       limits = c(-30, 30), oob = scales::squish) +
  scale_size_area(max_size = 9, name = "Total DEGs",
                  breaks = c(5, 10, 20, 30, 40)) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(x = NULL, y = NULL, tag = "A",
       title = "Heat-stress gene disruption per category \u00d7 condition",
       subtitle = "Bubble size = total DEGs;  fill = net direction;  \u2715 = no DEG in that cell") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, color = "grey35"),
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 8),
        axis.text.y   = element_text(size = 8.5),
        legend.box    = "horizontal",
        legend.key.height = unit(10, "pt"),
        legend.spacing.x  = unit(8, "pt"))

# ==================================================================
# Panel B — HSF vs sHSP scatter, FACETED BY STAGE
# ==================================================================
counts <- fread(file.path(DATA_ROOT, "counts_matrix.txt"),
                sep = "\t", header = TRUE, data.table = FALSE)
rownames(counts) <- counts[,1]; counts[,1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"
keep   <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]
lib    <- colSums(counts)
lcpm   <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)

meta <- read_tsv(file.path(DATA_ROOT, "sample_metadata.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
rownames(meta) <- meta$Sample
meta <- meta[colnames(lcpm), ]
meta$Maturity_stage <- factor(meta$Maturity_stage, levels = stages)

shsp_genes <- intersect(unique(heat_de$Gene[grepl("sHSP", heat_de$Category)]),
                        rownames(lcpm))
hsf_genes  <- intersect(unique(heat_de$Gene[grepl("HSF",  heat_de$Category)]),
                        rownames(lcpm))

z_row <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - m) / s
}
shsp_z <- t(apply(lcpm[shsp_genes, ], 1, z_row))
hsf_z  <- t(apply(lcpm[hsf_genes,  ], 1, z_row))

mod <- data.frame(
  Sample     = colnames(lcpm),
  HSF_score  = colMeans(hsf_z),
  sHSP_score = colMeans(shsp_z)
) %>% left_join(meta, by = "Sample") %>%
  mutate(Infection_Group = factor(Infection_Group, levels = c("RB(-)","RB(+)")))

# Per-facet, per-group r values — placed INSIDE each panel (top-left corner)
# so the narrow 1x4 strip only has to hold the stage name.
cor_annot <- mod %>% group_by(Maturity_stage, Infection_Group) %>%
  summarise(r = cor(HSF_score, sHSP_score, use = "pairwise.complete.obs"),
            .groups = "drop") %>%
  mutate(label = sprintf("%s  r = %.2f",
                         ifelse(Infection_Group == "RB(-)",
                                "RB\u2212", "RB+"),
                         r),
         # RB- on top line, RB+ below
         vjust_pos = ifelse(Infection_Group == "RB(-)", 1.4, 2.9))

# Strip labeller: full stage names (panel B is widened below to accommodate)
stage_names <- unique(as.character(mod$Maturity_stage))
stage_strip_label <- setNames(GRBV_STAGE_LABEL[stage_names], stage_names)

infection_cols <- c("RB(-)" = PAL$RB_neg, "RB(+)" = PAL$RB_pos)

pB <- ggplot(mod, aes(x = HSF_score, y = sHSP_score,
                      color = Infection_Group)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.65, size = 1.3) +
  geom_smooth(aes(group = Infection_Group), method = "lm",
              se = FALSE, linewidth = 0.8) +
  geom_text(data = cor_annot,
            aes(x = -Inf, y = Inf, label = label,
                color = Infection_Group, vjust = vjust_pos),
            hjust = -0.08, size = 2.5, fontface = "bold",
            show.legend = FALSE, inherit.aes = FALSE) +
  facet_wrap(~ Maturity_stage, ncol = 4,
             labeller = as_labeller(stage_strip_label)) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1, 2), expand = expansion(mult = c(0.05, 0.18))) +
  scale_color_manual(values = infection_cols, name = NULL) +
  labs(x = "HSF module z-score",
       y = "sHSP module z-score",
       tag = "B",
       title = "HSF \u2194 sHSP coordination, split by maturity stage") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 11, face = "bold"),
        strip.text  = element_text(size = 8.5, face = "bold",
                                   margin = margin(t = 3, b = 3)),
        legend.position = "top",
        legend.text  = element_text(size = 9),
        panel.spacing.x = unit(8, "pt"),
        axis.text    = element_text(size = 8))

# ==================================================================
# Panel C — EXPRESSION HEATMAP of top-20 recurrent heat-DEGs
# ==================================================================
top20 <- heat_rec %>%
  filter(n_conditions >= 4) %>%
  arrange(desc(n_conditions), mean_FC) %>%
  distinct(Gene, .keep_all = TRUE) %>%   # dedup across Category
  slice_head(n = 20)

# Mean log2FC per (gene, season, stage) from the per-condition DE table
fc_mat <- heat_de %>%
  filter(Gene %in% top20$Gene) %>%
  group_by(Gene, Season, Stage) %>%
  summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
  mutate(cond = paste(Season, GRBV_STAGE_LABEL[as.character(Stage)],
                      sep = " "),
         cond = factor(cond,
                       levels = unique(cond[order(Season, Stage)])))

# Hierarchical clustering of rows on the log2FC matrix
m <- fc_mat %>%
  select(Gene, cond, log2FC) %>%
  pivot_wider(names_from = cond, values_from = log2FC, values_fill = 0) %>%
  as.data.frame()
rownames(m) <- m$Gene; m$Gene <- NULL
gene_order <- rownames(m)[hclust(dist(as.matrix(m)))$order]

fc_mat$Gene <- factor(fc_mat$Gene, levels = rev(gene_order))
fc_mat <- fc_mat %>%
  left_join(top20 %>% distinct(Gene, Gene_Name), by = "Gene") %>%
  mutate(row_label = sprintf("%s \u2014 %s", Gene,
                             substr(Gene_Name, 1, 28)),
         row_label = factor(row_label,
                            levels = unique(row_label[order(Gene)])))

# Keep row_label and Gene aligned so clustering order is respected
row_lookup <- fc_mat %>% distinct(Gene, row_label) %>%
  mutate(Gene = factor(Gene, levels = rev(gene_order))) %>%
  arrange(Gene)
fc_mat$row_label <- factor(fc_mat$row_label, levels = row_lookup$row_label)

fc_limit <- max(abs(fc_mat$log2FC), na.rm = TRUE)
fc_limit <- min(fc_limit, 5)

pC <- ggplot(fc_mat, aes(x = cond, y = row_label, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", log2FC),
                color = abs(log2FC) > fc_limit * 0.55),
            family = "Arial", size = 1.9, show.legend = FALSE) +
  scale_fill_gradient2(low = PAL$down, mid = "white", high = PAL$up,
                       midpoint = 0,
                       limits = c(-fc_limit, fc_limit),
                       oob = scales::squish,
                       name = "log2FC") +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, tag = "C",
       title = "Top-20 recurrent heat-DEGs (\u22654 of 8 conditions)",
       subtitle = "Rows ordered by hierarchical clustering;  tile values show mean log2FC per condition") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, color = "grey35"),
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 8),
        axis.text.y   = element_text(size = 7),
        legend.key.height = unit(20, "pt"))

# ==================================================================
# Compose: single column, 3 panels stacked
# ==================================================================
layout <- pA / pB / pC +
  plot_layout(heights = c(1.1, 0.75, 1.3)) +
  plot_annotation(
    title = "Heat-stress disruption: sHSP suppression and HSF coordination",
    subtitle = "115-gene heat panel (9 categories);  HSF + sHSP module z-scores per sample;  rows in C clustered",
    theme = theme(plot.title    = element_text(family = "Arial",
                                               face = "bold", size = 13,
                                               hjust = 0.5,
                                               margin = margin(b = 2)),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6)))
  )

# Widen the composite beyond journal double-column (183 mm) so the 1x4
# facets in panel B fit full stage labels ("Pre-véraison", "Post-véraison")
# without clipping. At 230 mm wide, each B-facet gets ~52 mm.
save_grbv(layout, OUT_STEM,
          width_mm = 230, height_mm = 310)

cat(sprintf("[main04] panel A balloon cells: %d categories \u00d7 %d conditions\n",
            length(cat_order), 8))
cor_all <- cor(mod$HSF_score, mod$sHSP_score, use = "pairwise.complete.obs")
cat(sprintf("[main04] HSF-sHSP r (all) = %.3f\n", cor_all))
cat(sprintf("[main04] panel C heatmap: %d genes \u00d7 8 conditions\n",
            length(gene_order)))
