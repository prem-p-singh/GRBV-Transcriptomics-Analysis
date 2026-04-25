# S81 — Carbon translocation / phloem-loading module: log2 fold-change
# heatmap across every (Season, Maturity_stage, Years-Since-Infection cohort)
# cell.  Companion to S75's summary view.
#
# Cell value = mean log2(CPM+1) of RB(+) at (season, stage, cohort)
#              minus mean log2(CPM+1) of RB(-) at matched (season, stage).

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(data.table); library(scales); library(stringr)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
RESULTS  <- file.path(ROOT, "farm_pipeline", "results", "15_carbon_translocation")
OUT_DIR  <- file.path(ROOT, "figures", "supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S81_sugar_full_heatmap")

STG_LEVELS <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STG_LABELS <- c("Pre-ver","Ver","Post-ver","Harvest")

# Order subclasses functionally: source (loading) -> transport -> sink (utilisation)
SUB_LEVELS <- c("SUT_SUC","STP_HXT","SWEET","Cell_wall_INV","Vacuolar_INV",
                "Sucrose_synth","SPS","Sucrose_phosphat","Hexokinase",
                "Trehalose","Phloem_marker")
SUB_LABELS <- c("SUT/SUC","STP/HXT","SWEET","Cell-wall INV","Vacuolar INV",
                "Sucrose synthase","SPS","Sucrose-P phosphatase","Hexokinase",
                "Trehalose","Phloem marker")

# ---- load --------------------------------------------------------------
counts <- as.matrix(fread(file.path(ROOT, "counts_matrix.txt"), data.table = FALSE)
                    %>% { rownames(.) <- .[[1]]; .[, -1, drop = FALSE] })
meta <- fread(file.path(ROOT, "sample_metadata_v2.txt"), data.table = FALSE)
meta$Sampling_Season       <- as.character(meta$Sampling_Season)
meta$Years_Since_Infection <- suppressWarnings(as.integer(meta$Years_Since_Infection))

sub_df <- fread(file.path(RESULTS, "sugar_module_genes.csv"), data.table = FALSE)

common <- intersect(colnames(counts), meta$Sample)
counts <- counts[, common]; meta <- meta[match(common, meta$Sample), ]
keep_genes <- which(rowSums(counts >= 10) >= 0.10 * ncol(counts))
counts <- counts[keep_genes, , drop = FALSE]
lib    <- colSums(counts)
logcpm <- log2(t(t(counts) / lib) * 1e6 + 1)

genes_in <- intersect(rownames(logcpm), sub_df$gene_id)
logcpm   <- logcpm[genes_in, , drop = FALSE]
sub_df   <- sub_df[match(genes_in, sub_df$gene_id), ]
cat(sprintf("[S81] %d sugar-transport genes in expression matrix\n", length(genes_in)))

# ---- aggregate ---------------------------------------------------------
agg_one <- function(gn) {
  expr <- logcpm[gn, ]
  data.frame(gene_id = gn, expr = expr,
             Sampling_Season = meta$Sampling_Season,
             Maturity_stage  = meta$Maturity_stage,
             Years_Since_Infection = meta$Years_Since_Infection,
             Infection_Group = meta$Infection_Group,
             stringsAsFactors = FALSE)
}
long <- as.data.table(do.call(rbind, lapply(genes_in, agg_one)))

rb_pos <- long[Infection_Group == "RB(+)" & !is.na(Years_Since_Infection),
               .(mean_pos = mean(expr, na.rm = TRUE), n_pos = .N),
               by = .(gene_id, Sampling_Season, Maturity_stage, Years_Since_Infection)]
rb_neg <- long[Infection_Group == "RB(-)",
               .(mean_neg = mean(expr, na.rm = TRUE), n_neg = .N),
               by = .(gene_id, Sampling_Season, Maturity_stage)]
fc <- merge(rb_pos, rb_neg,
            by = c("gene_id","Sampling_Season","Maturity_stage"),
            all.x = TRUE)
fc[, log2FC := mean_pos - mean_neg]
fc <- merge(fc, sub_df[, c("gene_id","subclass")], by = "gene_id", all.x = TRUE)

# Filter subclasses present + apply functional ordering
fc[, subclass := factor(subclass,
                        levels = SUB_LEVELS[SUB_LEVELS %in% unique(fc$subclass)],
                        labels = SUB_LABELS[match(SUB_LEVELS[SUB_LEVELS %in% unique(fc$subclass)],
                                                  SUB_LEVELS)])]

# Gene order: by subclass then by mean log2FC
gene_order <- fc[, .(mean_fc = mean(log2FC, na.rm = TRUE)), by = .(gene_id, subclass)]
gene_order <- gene_order[order(subclass, -mean_fc)]
fc[, gene_id := factor(gene_id, levels = gene_order$gene_id)]

fc[, Maturity_stage := factor(Maturity_stage, levels = STG_LEVELS, labels = STG_LABELS)]
fc[, cohort := factor(Years_Since_Infection)]
fc_keep <- fc[!is.na(log2FC)]

# ---- main heatmap (all 33 genes, faceted by Season x Stage with cohort cols) ---
p_main <- ggplot(fc_keep, aes(x = cohort, y = gene_id, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.1) +
  facet_grid(subclass ~ Sampling_Season + Maturity_stage,
             scales = "free_y", space = "free_y",
             labeller = label_value, switch = "y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-3, 3), oob = squish,
                       name = expression(log[2]*FC)) +
  labs(title = "S81 — Carbon translocation / phloem-loading module across season x cohort x stage",
       subtitle = sprintf("%d genes; cell = mean log2(CPM+1) RB(+) at (season, stage, cohort) - mean log2(CPM+1) RB(-) at (season, stage)",
                          length(unique(fc_keep$gene_id))),
       x = "Years-Since-Infection cohort", y = NULL) +
  theme_grbv(base_size = 8) +
  theme(axis.text.y  = element_text(size = 6, family = "mono"),
        axis.text.x  = element_text(size = 6),
        strip.text.x = element_text(size = 6),
        strip.text.y.left = element_text(size = 7, angle = 0, face = "bold", hjust = 1),
        panel.spacing.x = unit(0.05, "lines"),
        panel.spacing.y = unit(0.10, "lines"),
        legend.position = "right",
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.5, color = "grey30"))

ggsave(paste0(OUT_STEM, ".png"), p_main,
       width = GRBV$W_DOUBLE_MM, height = 230, units = "mm", dpi = 600,
       limitsize = FALSE)
ggsave(paste0(OUT_STEM, ".svg"), p_main,
       width = GRBV$W_DOUBLE_MM, height = 230, units = "mm", limitsize = FALSE)
cat("[S81] saved\n")
