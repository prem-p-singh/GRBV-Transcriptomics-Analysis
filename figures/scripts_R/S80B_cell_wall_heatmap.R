# S80B — Cell-wall genes log2 FC heatmap.
#
# 161 cell-wall genes across (season x stage x cohort).  Subcategorised
# functionally (pectin / hemicellulose / cellulose / lignin / expansin /
# XTH / PG / PME / other) via regex on the GFF description column.  Row
# order = hierarchical clustering (Euclidean / Ward.D2) within each
# subcategory; subcategory blocks separated by horizontal rules.
#
# Gene labelling strategy: only label genes whose |log2FC| >= 1 in at
# least 4 of the 8 (season, stage) outer cells (averaged across cohorts).
# All other genes shown as anonymous tile rows for the global pattern.

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
OUT_STEM <- file.path(OUT_DIR, "S80B_cell_wall_heatmap")

STG_LEVELS <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STG_LABELS <- c("Pre-vér","Vér","Post-vér","Harvest")

# Sub-category regex patterns (cell-wall functional groups)
SUBCAT_RULES <- list(
  "Pectin"        = "pectin methylesterase|pectinesterase|polygalacturonase|pectate lyase|\\bPME\\b|\\bPG\\b",
  "Hemicellulose" = "xyloglucan|XTH|beta-galactosidase|alpha-arabinofuranosidase",
  "Cellulose"     = "cellulose synthase|\\bCESA\\b|\\bCSL\\b|endoglucanase|cellulase|COBRA",
  "Lignin"        = "laccase|cinnamoyl-CoA reductase|\\bCCR\\b|cinnamyl alcohol dehydrogenase|\\bCAD\\b|caffeoyl-CoA O-methyltransferase|\\bCCoAOMT\\b|peroxidase",
  "Expansin"      = "expansin|\\bEXP\\b|\\bEXLA\\b|\\bEXLB\\b"
)
SUBCAT_LEVELS <- c("Pectin","Hemicellulose","Cellulose","Lignin","Expansin","Other")

# ---- load ----------------------------------------------------------
counts_df <- fread(file.path(ROOT, "counts_matrix.txt"), data.table = FALSE)
rownames(counts_df) <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])
meta <- fread(file.path(ROOT, "sample_metadata_v2.txt"), data.table = FALSE)
meta$Sampling_Season       <- as.character(meta$Sampling_Season)
meta$Years_Since_Infection <- suppressWarnings(as.integer(meta$Years_Since_Infection))

cur <- fread(file.path(RESULTS, "curated_gene_table.csv"), data.table = FALSE)
cur <- cur[!duplicated(cur$gene_id), ]
cur <- cur[cur$category == "Cell_wall", ]

# Assign sub-category by regex on description
cur$subcat <- "Other"
for (sc in names(SUBCAT_RULES)) {
  hits <- grepl(SUBCAT_RULES[[sc]], cur$description, ignore.case = TRUE)
  cur$subcat[cur$subcat == "Other" & hits] <- sc
}
cur$subcat <- factor(cur$subcat, levels = SUBCAT_LEVELS)
cat("[S80B] Cell-wall subcategory counts:\n"); print(table(cur$subcat))

common <- intersect(colnames(counts), meta$Sample)
counts <- counts[, common]; meta <- meta[match(common, meta$Sample), ]
keep_genes <- which(rowSums(counts >= 10) >= 0.10 * ncol(counts))
counts <- counts[keep_genes, , drop = FALSE]
lib    <- colSums(counts)
logcpm <- log2(t(t(counts) / lib) * 1e6 + 1)

genes_in <- intersect(rownames(logcpm), cur$gene_id)
logcpm   <- logcpm[genes_in, , drop = FALSE]
cur      <- cur[match(genes_in, cur$gene_id), ]
cat(sprintf("[S80B] %d cell-wall genes in expression matrix\n", length(genes_in)))

# ---- aggregate to (Season, Stage, Cohort) -----------------------
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
fc <- merge(fc, cur[, c("gene_id","subcat")], by = "gene_id", all.x = TRUE)
fc[, Maturity_stage := factor(Maturity_stage, levels = STG_LEVELS, labels = STG_LABELS)]
fc[, cohort := factor(Years_Since_Infection,
                      levels = c(2016,2019,2020,2021,2022),
                      labels = c("'16","'19","'20","'21","'22"))]

# ---- |log2FC| >= 1 in at least 4 of 8 (season, stage) outer cells ------
cell_means <- fc[!is.na(log2FC),
                 .(mean_fc = mean(log2FC, na.rm = TRUE)),
                 by = .(gene_id, Sampling_Season, Maturity_stage)]
cell_means[, big := abs(mean_fc) >= 1]
big_count <- cell_means[, .(n_big = sum(big)), by = gene_id]
labelled_genes <- big_count[n_big >= 4]$gene_id
cat(sprintf("[S80B] Genes meeting |FC|>=1 in >=4/8 outer cells: %d\n",
            length(labelled_genes)))

# ---- per-subcategory hierarchical clustering --------------------------
cluster_within <- function(sc) {
  fc_sub <- fc[subcat == sc]
  if (nrow(fc_sub) == 0) return(character(0))
  fc_sub[, cell := paste(Sampling_Season, Maturity_stage, cohort, sep="|")]
  mat <- dcast(fc_sub, gene_id ~ cell, value.var = "log2FC")
  rn <- mat$gene_id; mat <- as.matrix(mat[, -1])
  rownames(mat) <- rn
  mat[is.na(mat)] <- 0
  if (nrow(mat) < 2) return(rn)
  d <- dist(mat, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  rn[hc$order]
}

ordered_genes <- unlist(lapply(SUBCAT_LEVELS, cluster_within))
ordered_genes <- ordered_genes[ordered_genes %in% unique(fc$gene_id)]

fc[, gene_id := factor(gene_id, levels = rev(ordered_genes))]
fc_keep <- fc[!is.na(log2FC)]

# Build y-axis label: blank for unlabeled, gene_id for labelled
y_labels <- setNames(
  ifelse(ordered_genes %in% labelled_genes, ordered_genes, ""),
  ordered_genes
)

# ---- main heatmap ----------------------------------------------------
p_main <- ggplot(fc_keep, aes(x = cohort, y = gene_id, fill = log2FC)) +
  geom_tile(color = NA, linewidth = 0) +
  facet_nested(subcat ~ Sampling_Season + Maturity_stage,
               nest_line = element_line(color = "grey50", linewidth = 0.4),
               strip = strip_nested(size = "variable"),
               scales = "free", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish, name = expression(log[2]*FC)) +
  scale_y_discrete(breaks = names(y_labels), labels = y_labels) +
  scale_x_discrete(expand = expansion(0)) +
  labs(title = sprintf("S80B — Cell-wall genes log2 FC across season x cohort x stage  (n=%d genes)",
                       length(genes_in)),
       subtitle = paste0("Cell = mean log2(CPM+1) RB(+) at (season, stage, cohort) - mean log2(CPM+1) RB(-) at (season, stage). ",
                         sprintf("%d genes meet |log2FC|>=1 in >=4/8 (season, stage) cells -> labelled.",
                                 length(labelled_genes))),
       caption = paste0("Cohort label = year of infection.  '16 = 5-6 yrs since infection at 2021 sampling, 9-10 yrs at 2022 sampling.  ",
                        "'22 = 0-1 yr.  Rows clustered (Euclidean / Ward.D2) within each functional sub-category."),
       x = NULL, y = NULL) +
  theme_grbv(base_size = 8) +
  theme(
    axis.text.y    = element_text(size = 5),
    axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text.x   = element_text(size = 8, face = "bold"),
    strip.text.y   = element_text(size = 8, face = "bold", angle = 0),
    panel.spacing  = unit(0.06, "lines"),
    legend.position = "right",
    plot.title    = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    plot.caption  = element_text(size = 7, color = "grey40", hjust = 0)
  )

ggsave(paste0(OUT_STEM, ".png"), p_main,
       width = 11, height = 18, units = "in", dpi = 600, limitsize = FALSE)
ggsave(paste0(OUT_STEM, ".svg"), p_main,
       width = 11, height = 18, units = "in", limitsize = FALSE)
cat("[S80B] saved\n")
