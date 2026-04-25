# S80 — Phenylpropanoid + PR + Cell-wall genes: log2 fold-change heatmap
# resolved across every (Sampling_Season, Maturity_stage, Years-Since-Infection
# cohort) cell.  Companion to S74's summary view.
#
# Cell value = mean log2(CPM+1) of RB(+) at (season, stage, cohort)
#              minus mean log2(CPM+1) of RB(-) at matched (season, stage).
# RB(-) controls don't have a cohort, so the same RB(-) baseline is reused
# across cohorts within each (season, stage) cell.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(data.table); library(scales); library(patchwork); library(stringr)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
RESULTS  <- file.path(ROOT, "farm_pipeline", "results", "14_pathway_per_gene")
OUT_DIR  <- file.path(ROOT, "figures", "supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S80_pathway_full_heatmap")

CAT_LEVELS <- c("Phenylpropanoid", "PR_proteins", "Cell_wall")
CAT_LABELS <- c("Phenylpropanoid", "PR proteins", "Cell wall")
STG_LEVELS <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STG_LABELS <- c("Pre-ver","Ver","Post-ver","Harvest")

# ---- load ----------------------------------------------------------------
counts <- as.matrix(fread(file.path(ROOT, "counts_matrix.txt"), data.table = FALSE)
                    %>% { rownames(.) <- .[[1]]; .[, -1, drop = FALSE] })
meta <- fread(file.path(ROOT, "sample_metadata_v2.txt"), data.table = FALSE)
meta$Sampling_Season       <- as.character(meta$Sampling_Season)
meta$Years_Since_Infection <- suppressWarnings(as.integer(meta$Years_Since_Infection))

cur <- fread(file.path(RESULTS, "curated_gene_table.csv"), data.table = FALSE)
cur <- cur[!duplicated(cur$gene_id), ]   # one row per gene

# Filter to genes in expression matrix and meta-known samples
common <- intersect(colnames(counts), meta$Sample)
counts <- counts[, common]; meta <- meta[match(common, meta$Sample), ]
keep_genes <- which(rowSums(counts >= 10) >= 0.10 * ncol(counts))
counts <- counts[keep_genes, , drop = FALSE]

lib    <- colSums(counts)
logcpm <- log2(t(t(counts) / lib) * 1e6 + 1)

# Restrict to curated genes
genes_in <- intersect(rownames(logcpm), cur$gene_id)
logcpm <- logcpm[genes_in, , drop = FALSE]
cur <- cur[match(genes_in, cur$gene_id), ]
cur$category <- factor(cur$category, levels = CAT_LEVELS)
cat(sprintf("[S80] %d curated genes in expression matrix\n", length(genes_in)))

# ---- aggregate to (Season, Stage, Cohort, Infection) ---------------------
agg <- function() {
  out <- list()
  for (gn in genes_in) {
    expr <- logcpm[gn, ]
    df <- data.frame(
      gene_id = gn, expr = expr,
      Sampling_Season = meta$Sampling_Season,
      Maturity_stage  = meta$Maturity_stage,
      Years_Since_Infection = meta$Years_Since_Infection,
      Infection_Group = meta$Infection_Group,
      stringsAsFactors = FALSE
    )
    out[[gn]] <- df
  }
  rbindlist(out)
}
long <- agg()

# RB(+) means per (season, stage, cohort)
rb_pos <- long[Infection_Group == "RB(+)" & !is.na(Years_Since_Infection),
               .(mean_pos = mean(expr, na.rm = TRUE), n_pos = .N),
               by = .(gene_id, Sampling_Season, Maturity_stage, Years_Since_Infection)]

# RB(-) means per (season, stage)  — controls don't have a cohort
rb_neg <- long[Infection_Group == "RB(-)",
               .(mean_neg = mean(expr, na.rm = TRUE), n_neg = .N),
               by = .(gene_id, Sampling_Season, Maturity_stage)]

# Join: for each (gene, season, stage, cohort), compute log2 FC
fc <- merge(rb_pos, rb_neg,
            by = c("gene_id","Sampling_Season","Maturity_stage"),
            all.x = TRUE)
fc[, log2FC := mean_pos - mean_neg]
fc <- merge(fc, cur[, c("gene_id","category")], by = "gene_id", all.x = TRUE)

# Order genes within each category by mean log2FC (most-up at top)
gene_order <- fc[, .(mean_fc = mean(log2FC, na.rm = TRUE)), by = .(gene_id, category)]
gene_order <- gene_order[order(category, -mean_fc)]
fc[, gene_id := factor(gene_id, levels = gene_order$gene_id)]
fc[, Maturity_stage := factor(Maturity_stage, levels = STG_LEVELS, labels = STG_LABELS)]
fc[, cohort := factor(Years_Since_Infection)]

# Drop genes that are NA everywhere (no RB+ data)
fc_keep <- fc[!is.na(log2FC)]

# ---- plot per category ---------------------------------------------------
make_panel <- function(category_label, fc_sub, ngene_label) {
  ggplot(fc_sub, aes(x = cohort, y = gene_id, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.05) +
    facet_grid(. ~ Sampling_Season + Maturity_stage,
               labeller = label_both, switch = "x") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-3, 3), oob = squish,
                         name = expression(log[2]*FC)) +
    labs(title = sprintf("%s (n=%d genes)", category_label, ngene_label),
         x = "Years-Since-Infection cohort", y = NULL) +
    theme_grbv(base_size = 7) +
    theme(axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          strip.text  = element_text(size = 6),
          panel.spacing = unit(0.05, "lines"),
          legend.position = "right",
          plot.title = element_text(size = 9, face = "bold"))
}

panels <- list()
for (i in seq_along(CAT_LEVELS)) {
  cat_df <- fc_keep[category == CAT_LEVELS[i]]
  if (nrow(cat_df) == 0) next
  ngene <- length(unique(cat_df$gene_id))
  panels[[i]] <- make_panel(CAT_LABELS[i], cat_df, ngene)
}

fig <- (panels[[1]] / panels[[2]] / panels[[3]]) +
  plot_layout(heights = c(0.35, 0.30, 0.55)) +
  plot_annotation(
    title = "S80 — Phenylpropanoid + PR + Cell-wall log2 FC across season x cohort x stage",
    subtitle = "Cell = mean log2(CPM+1) RB(+) at (season, stage, cohort) - mean log2(CPM+1) RB(-) at (season, stage)",
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7.5, color = "grey30")))

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = GRBV$W_DOUBLE_MM, height = 320, units = "mm",
       dpi = 600, limitsize = FALSE)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = GRBV$W_DOUBLE_MM, height = 320, units = "mm",
       limitsize = FALSE)
cat("[S80] saved (3 categories,", length(unique(fc_keep$gene_id)), "genes total)\n")
