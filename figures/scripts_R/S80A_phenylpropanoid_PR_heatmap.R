# S80A — Phenylpropanoid + PR proteins log2 FC heatmap.
# Two stacked panels (Phenylpropanoid n=51, PR n=40), full season x cohort
# x stage resolution.  Hierarchical clustering (Euclidean / Ward.D2) applied
# per-panel for explicit row order.
#
# Cell value = mean log2(CPM+1) RB(+) at (season, stage, cohort) -
#              mean log2(CPM+1) RB(-) at matched (season, stage).

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(data.table); library(scales); library(patchwork); library(ggh4x)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
RESULTS  <- file.path(ROOT, "farm_pipeline", "results", "14_pathway_per_gene")
OUT_DIR  <- file.path(ROOT, "figures", "supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S80A_phenylpropanoid_PR_heatmap")

CAT_LEVELS <- c("Phenylpropanoid", "PR_proteins")
CAT_LABELS <- c("Phenylpropanoid (n_genes)", "PR proteins (n_genes)")
STG_LEVELS <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STG_LABELS <- c("Pre-vér","Vér","Post-vér","Harvest")
COHORT_LABELS <- c("2016"="'16","2019"="'19","2020"="'20","2021"="'21","2022"="'22")

# ---- load --------------------------------------------------------------
counts_df <- fread(file.path(ROOT, "counts_matrix.txt"), data.table = FALSE)
rownames(counts_df) <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])
meta <- fread(file.path(ROOT, "sample_metadata_v2.txt"), data.table = FALSE)
meta$Sampling_Season       <- as.character(meta$Sampling_Season)
meta$Years_Since_Infection <- suppressWarnings(as.integer(meta$Years_Since_Infection))

cur <- fread(file.path(RESULTS, "curated_gene_table.csv"), data.table = FALSE)
cur <- cur[!duplicated(cur$gene_id), ]
cur <- cur[cur$category %in% CAT_LEVELS, ]

common <- intersect(colnames(counts), meta$Sample)
counts <- counts[, common]; meta <- meta[match(common, meta$Sample), ]
keep_genes <- which(rowSums(counts >= 10) >= 0.10 * ncol(counts))
counts <- counts[keep_genes, , drop = FALSE]
lib    <- colSums(counts)
logcpm <- log2(t(t(counts) / lib) * 1e6 + 1)

genes_in <- intersect(rownames(logcpm), cur$gene_id)
logcpm <- logcpm[genes_in, , drop = FALSE]
cur    <- cur[match(genes_in, cur$gene_id), ]

# ---- aggregate to (Season, Stage, Cohort) -----------------------------
build_long <- function(genes, lcpm, m) {
  rbindlist(lapply(genes, function(gn) {
    data.table(gene_id = gn, expr = lcpm[gn, ],
               Sampling_Season = m$Sampling_Season,
               Maturity_stage  = m$Maturity_stage,
               Years_Since_Infection = m$Years_Since_Infection,
               Infection_Group = m$Infection_Group)
  }))
}
long <- build_long(genes_in, logcpm, meta)

rb_pos <- long[Infection_Group == "RB(+)" & !is.na(Years_Since_Infection),
               .(mean_pos = mean(expr, na.rm = TRUE)),
               by = .(gene_id, Sampling_Season, Maturity_stage, Years_Since_Infection)]
rb_neg <- long[Infection_Group == "RB(-)",
               .(mean_neg = mean(expr, na.rm = TRUE)),
               by = .(gene_id, Sampling_Season, Maturity_stage)]

fc <- merge(rb_pos, rb_neg,
            by = c("gene_id","Sampling_Season","Maturity_stage"),
            all.x = TRUE)
fc[, log2FC := mean_pos - mean_neg]
fc <- merge(fc, cur[, c("gene_id","category")], by = "gene_id", all.x = TRUE)

fc[, Maturity_stage := factor(Maturity_stage, levels = STG_LEVELS, labels = STG_LABELS)]
fc[, cohort := factor(Years_Since_Infection,
                      levels = c(2016,2019,2020,2021,2022),
                      labels = c("'16","'19","'20","'21","'22"))]

# ---- hierarchical clustering per category -----------------------------
cluster_genes <- function(cat) {
  fc_cat <- fc[category == cat]
  # Build (gene x cell) matrix, NA -> 0 for clustering only
  fc_cat[, cell := paste(Sampling_Season, Maturity_stage, cohort, sep="|")]
  mat <- dcast(fc_cat, gene_id ~ cell, value.var = "log2FC")
  rn <- mat$gene_id; mat <- as.matrix(mat[, -1])
  rownames(mat) <- rn
  mat[is.na(mat)] <- 0
  if (nrow(mat) < 2) return(rn)
  d <- dist(mat, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  rn[hc$order]
}

ord_phen <- cluster_genes("Phenylpropanoid")
ord_pr   <- cluster_genes("PR_proteins")

# ---- panel builder -----------------------------------------------------
make_panel <- function(cat_key, cat_label, gene_order, ngene) {
  fc_cat <- fc[category == cat_key & !is.na(log2FC)]
  fc_cat[, gene_id := factor(gene_id, levels = rev(gene_order))]

  ggplot(fc_cat, aes(x = cohort, y = gene_id, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.05) +
    facet_nested(. ~ Sampling_Season + Maturity_stage,
                 nest_line = element_line(color = "grey50", linewidth = 0.4),
                 strip = strip_nested(size = "variable"),
                 scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0,
                         limits = c(-3, 3), oob = scales::squish,
                         name = expression(log[2]*FC)) +
    scale_x_discrete(expand = expansion(0)) +
    labs(title = sub("n_genes", ngene, cat_label),
         x = NULL, y = NULL) +
    theme_grbv(base_size = 8) +
    theme(
      axis.text.y    = element_text(size = 6),
      axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      strip.text     = element_text(size = 8, face = "bold"),
      panel.spacing.x = unit(0.05, "lines"),
      legend.position = "right",
      plot.title     = element_text(size = 11, face = "bold")
    )
}

p_phen <- make_panel("Phenylpropanoid", CAT_LABELS[1], ord_phen, length(ord_phen))
p_pr   <- make_panel("PR_proteins",     CAT_LABELS[2], ord_pr,   length(ord_pr))

fig <- (p_phen / p_pr) +
  plot_layout(heights = c(length(ord_phen), length(ord_pr))) +
  plot_annotation(
    title = "S80A — Phenylpropanoid + PR proteins log2 FC across season x cohort x stage",
    subtitle = "Cell = mean log2(CPM+1) RB(+) at (season, stage, cohort) - mean log2(CPM+1) RB(-) at (season, stage). Diverging blue/white/red, capped at +/-3.",
    caption = paste0("Cohort label = year of infection.  '16 = 5-6 yrs since infection at 2021 sampling, 9-10 yrs at 2022 sampling.  ",
                     "'22 = 0-1 yr.  Rows hierarchically clustered (Euclidean / Ward.D2) per panel."),
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 8, color = "grey30"),
                  plot.caption = element_text(size = 7, color = "grey40", hjust = 0)))

ggsave(paste0(OUT_STEM, ".png"), fig,
       width = 10, height = 14, units = "in", dpi = 600, limitsize = FALSE)
ggsave(paste0(OUT_STEM, ".svg"), fig,
       width = 10, height = 14, units = "in", limitsize = FALSE)
cat("[S80A] saved (Phenylpropanoid n=", length(ord_phen),
    ", PR n=", length(ord_pr), ")\n", sep = "")
