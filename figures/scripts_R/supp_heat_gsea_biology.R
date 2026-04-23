# Supplementary: heat stress + GSEA detail + biology cartoon — S34..S39.
#
#   S34  heat-gene landscape (DEG counts per category x condition)
#   S35  HSF <-> sHSP axis scatter (full per-sample)
#   S36  biology cartoon (cascade, same content as Main 10 standalone)
#   S37  GSEA running enrichment plots for top-4 gene sets
#   S38  GSEA vs Fisher ORA concordance
#   S39  GSEA NES heatmap (standalone, all 20 top sets)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table)
  library(grid); library(ggrepel)
})
HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))
DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- file.path(DATA_ROOT, "figures", "supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stages <- GRBV_STAGES

# -----------------------------------------------------------------
# S34 — heat-gene landscape (deg counts per category x stage x season)
# -----------------------------------------------------------------
heat_de <- fread(file.path(DATA_ROOT, "heat_stress_DEGs_all_conditions.csv")) %>%
  as.data.frame() %>%
  mutate(Stage = factor(Stage, levels = stages),
         Season = factor(Season),
         is_DEG = is_DEG == "True" | is_DEG == TRUE,
         dir = ifelse(log2FC > 0, "up", "down"))

cat_order <- heat_de %>% group_by(Category) %>%
  summarise(n = n_distinct(Gene)) %>% arrange(desc(n)) %>% pull(Category)

cat_cnt <- heat_de %>% filter(is_DEG) %>%
  count(Category, Stage, Season, dir) %>%
  group_by(Category, Stage, Season) %>%
  summarise(up = sum(n[dir == "up"]),
            down = sum(n[dir == "down"]),
            total = up + down, .groups = "drop") %>%
  complete(Category = cat_order,
           Stage = stages, Season = c("2021","2022"),
           fill = list(up = 0, down = 0, total = 0)) %>%
  pivot_longer(c(up, down), names_to = "dir", values_to = "n")

cat_cnt$Category <- factor(cat_cnt$Category, levels = rev(cat_order))

S34 <- ggplot(cat_cnt, aes(x = Stage, y = n, fill = dir)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75,
           color = "white", linewidth = 0.2) +
  facet_grid(Category ~ Season, scales = "free_y", switch = "y") +
  scale_x_discrete(labels = substr(GRBV_STAGE_LABEL, 1, 4)) +
  scale_fill_manual(values = c(up = unname(GRBV_COLORS["up"]),
                               down = unname(GRBV_COLORS["down"])),
                    labels = c(up = "Up", down = "Down"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "DEGs (padj<0.05)",
       title = "Supplementary S34. Heat-stress DEG landscape per category x stage x season") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6),
        strip.text  = element_text(size = 7.5, face = "bold"),
        strip.placement = "outside",
        legend.position = "top")
save_grbv(S34, file.path(OUT_DIR, "S34_heat_stress_gene_landscape"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 230)

# -----------------------------------------------------------------
# S35 — HSF-sHSP axis (more detail than Main 4 panel B)
# -----------------------------------------------------------------
counts <- fread(file.path(DATA_ROOT, "counts_matrix.txt"),
                sep = "\t", header = TRUE, data.table = FALSE)
rownames(counts) <- counts[,1]; counts[,1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts); counts <- counts[keep, ]
lib <- colSums(counts); lcpm <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)
meta <- read_tsv(file.path(DATA_ROOT, "sample_metadata.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
rownames(meta) <- meta$Sample; meta <- meta[colnames(lcpm), ]

shsp_genes <- intersect(unique(heat_de$Gene[grepl("sHSP", heat_de$Category)]),
                        rownames(lcpm))
hsf_genes  <- intersect(unique(heat_de$Gene[grepl("HSF", heat_de$Category)]),
                        rownames(lcpm))
z <- function(x) { m <- mean(x); s <- sd(x); if (s == 0) rep(0, length(x)) else (x - m)/s }
shsp_z <- t(apply(lcpm[shsp_genes, ], 1, z))
hsf_z  <- t(apply(lcpm[hsf_genes, ], 1, z))
mod <- data.frame(Sample = colnames(lcpm),
                  sHSP = colMeans(shsp_z),
                  HSF  = colMeans(hsf_z)) %>%
  left_join(meta, by = "Sample") %>%
  mutate(Infection_Group = factor(Infection_Group, levels = c("RB(-)","RB(+)")),
         Maturity_stage = factor(Maturity_stage, levels = stages))

S35 <- ggplot(mod, aes(x = HSF, y = sHSP, color = Infection_Group)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.6, size = 1.1) +
  geom_smooth(aes(group = Infection_Group), method = "lm",
              se = TRUE, linewidth = 0.6, alpha = 0.15) +
  facet_wrap(~ Maturity_stage, ncol = 4, labeller = as_labeller(GRBV_STAGE_LABEL)) +
  scale_color_manual(values = GRBV_INFECTION_COLOR, name = NULL) +
  labs(x = "HSF module z-score", y = "sHSP module z-score",
       title = "Supplementary S35. HSF <-> sHSP axis per maturity stage") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "top")
save_grbv(S35, file.path(OUT_DIR, "S35_heat_vs_virus_disentanglement"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 120)

# -----------------------------------------------------------------
# S36 — biology cartoon (reuse Main 10 logic as standalone)
# -----------------------------------------------------------------
box_h <- 1.45
y_centers <- seq(10.5 - box_h/2, 1.0 + box_h/2, length.out = 5)
nodes <- data.frame(
  y = y_centers,
  fill = c("#5d3a9b","#3498db","#c0392b","#e67e22","#27ae60"),
  title = c("1. GRBV infection establishes",
            "2. Canonical plant immunity dampened",
            "3. HSF <-> sHSP thermal-defence axis suppressed",
            "4. Sugar-metabolism transcription rewired",
            "5. Disease phenotype: red blotch + sugar deficit")
)
box <- data.frame(
  xmin = 0.4, xmax = 9.6,
  ymin = nodes$y - box_h/2, ymax = nodes$y + box_h/2,
  fill = nodes$fill
)
arrows <- data.frame(
  x = 5, xend = 5,
  y = head(nodes$y, -1) - box_h/2,
  yend = tail(nodes$y, -1) + box_h/2
)

S36 <- ggplot() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 11), expand = FALSE) +
  geom_rect(data = box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = box$fill, color = "black", linewidth = 0.5, alpha = 0.95) +
  geom_text(data = nodes,
            aes(x = 0.7, y = y, label = title),
            family = "Arial", color = "white", fontface = "bold",
            hjust = 0, vjust = 0.5, size = 4.0) +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.7, color = "grey25",
               arrow = arrow(length = unit(7, "pt"), type = "closed")) +
  labs(title = "Supplementary S36. GRBV molecular disruption cascade",
       caption = "Robust temperature -> titer mediators: MBF1c (LOC100249249), SPS1 (LOC100232974), uncharacterised LOC100266309") +
  theme_void(base_family = "Arial") +
  theme(plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.caption = element_text(size = 8, color = "grey35", hjust = 0.5),
        plot.margin  = margin(6, 6, 6, 6))
save_grbv(S36, file.path(OUT_DIR, "S36_GRBV_biology_diagram"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 170)

# -----------------------------------------------------------------
# GSEA base
# -----------------------------------------------------------------
gsea <- fread(file.path(DATA_ROOT, "gsea_results_all_conditions.csv")) %>%
  as.data.frame() %>%
  mutate(condition = gsub("_", " ", condition)) %>%
  separate(condition, into = c("season","stage_a","stage_b"),
           sep = " ", remove = FALSE, fill = "right") %>%
  mutate(stage = ifelse(is.na(stage_b), stage_a,
                        paste(stage_a, stage_b, sep = "_")),
         stage = factor(stage, levels = stages),
         season = factor(season, levels = c("2021","2022")),
         sig = !is.na(fdr) & fdr < 0.05)

# -----------------------------------------------------------------
# S37 — GSEA running ES for top-4 gene sets (approximation: ES bar profile)
# -----------------------------------------------------------------
top4 <- gsea %>% filter(sig) %>% group_by(gene_set) %>%
  summarise(mean_abs = mean(abs(nes), na.rm = TRUE),
            n_cond = n_distinct(condition[sig]),
            .groups = "drop") %>%
  arrange(desc(n_cond), desc(mean_abs)) %>% slice_head(n = 4)

top4_long <- gsea %>% filter(gene_set %in% top4$gene_set) %>%
  mutate(cond = paste(season,
                      substr(GRBV_STAGE_LABEL[as.character(stage)], 1, 4)),
         cond = factor(cond, levels = unique(cond[order(season, stage)])))

S37 <- ggplot(top4_long, aes(x = cond, y = es, fill = es)) +
  geom_col(color = "white", linewidth = 0.2, width = 0.75) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  facet_wrap(~ gene_set, ncol = 2, scales = "free_y") +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white", high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, guide = "none",
                       limits = c(-1, 1), oob = scales::squish) +
  labs(x = NULL, y = "Enrichment score (ES)",
       title = "Supplementary S37. GSEA ES across conditions for the four most perturbed sets") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        strip.text = element_text(size = 8, face = "bold"))
save_grbv(S37, file.path(OUT_DIR, "S37_gsea_enrichment_plots"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 130)

# -----------------------------------------------------------------
# S38 — GSEA vs Fisher ORA concordance (use NES vs p-value strength)
# -----------------------------------------------------------------
gsea_summ <- gsea %>%
  group_by(gene_set) %>%
  summarise(mean_nes = mean(nes, na.rm = TRUE),
            best_fdr = min(fdr, na.rm = TRUE),
            n_sig = sum(sig), .groups = "drop") %>%
  filter(!is.infinite(best_fdr)) %>%
  mutate(neg_log_fdr = -log10(pmax(best_fdr, 1e-30)))

S38 <- ggplot(gsea_summ, aes(x = mean_nes, y = neg_log_fdr,
                             color = n_sig)) +
  geom_point(alpha = 0.7, size = 1.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), color = "grey60",
             linetype = "dashed", linewidth = 0.3) +
  geom_text_repel(data = gsea_summ %>%
                    filter(n_sig >= 6 | neg_log_fdr > 8),
                  aes(label = gene_set), size = 2.2, family = "Arial",
                  segment.size = 0.2, max.overlaps = 20,
                  box.padding = 0.3) +
  scale_color_gradient(low = unname(GRBV_COLORS["neutral"]),
                       high = unname(GRBV_COLORS["RB_pos"]),
                       name = "sig conditions") +
  labs(x = "Mean NES across 8 conditions",
       y = expression(-log[10]*"(best FDR)"),
       title = "Supplementary S38. GSEA effect (mean NES) vs best-condition FDR") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.key.height = unit(15, "pt"))
save_grbv(S38, file.path(OUT_DIR, "S38_gsea_vs_fisher_concordance"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 130)

# -----------------------------------------------------------------
# S39 — GSEA NES heatmap (all top 30 gene sets, full detail)
# -----------------------------------------------------------------
top30 <- gsea %>% filter(sig) %>% group_by(gene_set) %>%
  summarise(max_abs = max(abs(nes), na.rm = TRUE),
            n_cond = n_distinct(condition[sig]),
            .groups = "drop") %>%
  arrange(desc(n_cond), desc(max_abs)) %>%
  slice_head(n = 30)

heat39 <- gsea %>% filter(gene_set %in% top30$gene_set) %>%
  mutate(cond = paste(season, GRBV_STAGE_LABEL[as.character(stage)]),
         cond = factor(cond, levels = unique(cond[order(season, stage)])),
         gene_set = factor(gene_set, levels = rev(top30$gene_set)))

S39 <- ggplot(heat39, aes(x = cond, y = gene_set, fill = nes)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(data = subset(heat39, sig),
            aes(label = "*"), family = "Arial", size = 2.6, color = "white") +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white", high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish, name = "NES") +
  labs(x = NULL, y = NULL,
       title = "Supplementary S39. GSEA NES heatmap (top 30 gene sets; * = FDR<0.05)") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6),
        legend.key.height = unit(22, "pt"))
save_grbv(S39, file.path(OUT_DIR, "S39_gsea_nes_heatmap"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 180)

cat("[supp_heat_gsea_biology] wrote S34, S35, S36, S37, S38, S39\n")
