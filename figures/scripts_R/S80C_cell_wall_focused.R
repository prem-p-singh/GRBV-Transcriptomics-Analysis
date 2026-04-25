# S80C — Cell-wall genes with consistent GRBV response.
# Companion to S80B (full landscape).  Filters the 161 cell-wall genes to
# the subset showing a substantive, consistent infection effect, then plots
# them at higher resolution with cell-text annotation and nested y-axis
# (subcategory-grouped row labels via ggh4x::scale_y_discrete(guide="axis_nested")).
#
# Selection criteria (gene retained if EITHER):
#   (i)  |mean log2FC| >= 1     in >=4 of 8 (season x stage) outer cells
#   (ii) |mean log2FC| >= 0.75  in >=6 of 8 outer cells
# Target n: 20-40.  Threshold (ii) is auto-tuned within {0.60, 0.70, 0.75,
# 0.80, 0.90} to land in the 20-40 window.
#
# Row labels: "Symbol — short_function" when an NCBI symbol is available,
# "LOC###### — best NCBI description (<=35 char)" otherwise.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(data.table); library(scales); library(patchwork); library(ggh4x)
  library(stringr)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

ROOT     <- "/Users/prempratapsingh/Data_Analysis"
RESULTS  <- file.path(ROOT, "farm_pipeline", "results", "14_pathway_per_gene")
ANNOTS   <- file.path(ROOT, "farm_pipeline", "results", "annotations", "gene_annotations.csv")
OUT_DIR  <- file.path(ROOT, "figures", "supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_STEM <- file.path(OUT_DIR, "S80C_cell_wall_focused")

STG_LEVELS <- c("Pre_veraison","Veraison","Post_veraison","Harvest")
STG_LABELS <- c("Pre-vér","Vér","Post-vér","Harvest")

# Sub-category regex on the description field.  Order matters — earlier
# rules win.
SUBCAT_RULES <- list(
  "Pectin"        = "pectin methylesterase|pectinesterase|polygalacturonase|pectate lyase|\\bPME\\b|\\bPG\\b|exopolygalacturonase",
  "Hemicellulose" = "xyloglucan|XTH|beta-galactosidase|alpha-arabinofuranosidase",
  "Cellulose"     = "cellulose synthase|\\bCESA\\b|\\bCSL\\b|endoglucanase|cellulase|COBRA",
  "Lignin"        = "laccase|cinnamoyl-CoA reductase|\\bCCR\\b|cinnamyl alcohol dehydrogenase|\\bCAD\\b|caffeoyl-CoA O-methyltransferase|\\bCCoAOMT\\b|peroxidase",
  "Expansin"      = "expansin|\\bEXP\\b|\\bEXLA\\b|\\bEXLB\\b"
)
SUBCAT_LEVELS <- c("Pectin","Hemicellulose","Cellulose","Lignin","Expansin","Other")

# ---- load expression -------------------------------------------------
counts_df <- fread(file.path(ROOT, "counts_matrix.txt"), data.table = FALSE)
rownames(counts_df) <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])
meta <- fread(file.path(ROOT, "sample_metadata_v2.txt"), data.table = FALSE)
meta$Sampling_Season       <- as.character(meta$Sampling_Season)
meta$Years_Since_Infection <- suppressWarnings(as.integer(meta$Years_Since_Infection))

cur <- fread(file.path(RESULTS, "curated_gene_table.csv"), data.table = FALSE)
cur <- cur[!duplicated(cur$gene_id), ]
cur <- cur[cur$category == "Cell_wall", ]

# Subcategorise via description regex
cur$subcat <- "Other"
for (sc in names(SUBCAT_RULES)) {
  hits <- grepl(SUBCAT_RULES[[sc]], cur$description, ignore.case = TRUE)
  cur$subcat[cur$subcat == "Other" & hits] <- sc
}

common <- intersect(colnames(counts), meta$Sample)
counts <- counts[, common]; meta <- meta[match(common, meta$Sample), ]
keep_genes <- which(rowSums(counts >= 10) >= 0.10 * ncol(counts))
counts <- counts[keep_genes, , drop = FALSE]
lib    <- colSums(counts)
logcpm <- log2(t(t(counts) / lib) * 1e6 + 1)

genes_in <- intersect(rownames(logcpm), cur$gene_id)
logcpm <- logcpm[genes_in, , drop = FALSE]
cur    <- cur[match(genes_in, cur$gene_id), ]
cat(sprintf("[S80C] %d cell-wall genes in expression matrix\n", length(genes_in)))

# ---- aggregate to (Season, Stage, Cohort) ---------------------------
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
fc <- merge(fc, cur[, c("gene_id","subcat","description")],
            by = "gene_id", all.x = TRUE)

# Outer-cell mean log2FC (averaged across cohorts) per (gene, season, stage)
cell_means <- fc[!is.na(log2FC),
                 .(mean_fc_cell = mean(log2FC, na.rm = TRUE)),
                 by = .(gene_id, Sampling_Season, Maturity_stage)]

n_strict <- cell_means[, .(n_big = sum(abs(mean_fc_cell) >= 1)), by = gene_id]
strict_set <- n_strict[n_big >= 4]$gene_id
cat(sprintf("[S80C] strict criterion |FC|>=1 in >=4/8 cells: %d genes\n",
            length(strict_set)))

# Auto-tune the broader criterion to land in the 20-40 window.  Cell-wall
# genes show generally muted GRBV effects, so we scan a (amplitude, min-cells)
# grid: amplitudes 0.30-0.90 in 0.05 steps, min-cells {6, 5, 4} in priority
# order.  Accept the first (amplitude, min-cells) pair giving 20-40 total.
amp_grid     <- seq(0.95, 0.30, by = -0.05)
mincell_grid <- c(6L, 5L, 4L)

best <- NULL
for (mc in mincell_grid) {
  for (thresh in amp_grid) {
    n_broad <- cell_means[, .(n_med = sum(abs(mean_fc_cell) >= thresh)), by = gene_id]
    broad_set <- n_broad[n_med >= mc]$gene_id
    sel <- union(strict_set, broad_set)
    if (length(sel) >= 20 && length(sel) <= 40) {
      best <- list(thresh = thresh, mincells = mc, selected = sel,
                   broad = broad_set)
      break
    }
  }
  if (!is.null(best)) break
}
if (is.null(best)) {
  # Fall back: pick the (mc=4, lowest amp) producing >=20
  for (thresh in amp_grid) {
    n_broad <- cell_means[, .(n_med = sum(abs(mean_fc_cell) >= thresh)), by = gene_id]
    broad_set <- n_broad[n_med >= 4]$gene_id
    sel <- union(strict_set, broad_set)
    if (length(sel) >= 20) {
      best <- list(thresh = thresh, mincells = 4, selected = sel, broad = broad_set)
      break
    }
  }
}
if (is.null(best)) stop("Could not select 20-40 genes; check criteria.")
selected     <- best$selected
final_thresh <- best$thresh
final_mc     <- best$mincells
cat(sprintf("[S80C] FINAL: %d genes selected (strict 13 + broad |FC|>=%.2f in >=%d/8 cells)\n",
            length(selected), final_thresh, final_mc))

# ---- annotate selected genes ---------------------------------------
ann <- fread(ANNOTS, data.table = FALSE)
# ann columns: gene_id, ncbi_gene_name, synonyms, product, biotype, ...

# Parse the description field into a "function string" — first segment of
# "{desc} | {Name} | {biotype}" pipeline-separated record.
parse_desc <- function(d) {
  if (is.na(d) || d == "") return("")
  parts <- str_trim(str_split(d, "\\|")[[1]])
  parts[1]   # the descriptive part comes first in rule 14's CSV
}
make_label <- function(gene_id, desc, ann_row) {
  fn <- parse_desc(desc)
  if (nchar(fn) > 35) fn <- paste0(substr(fn, 1, 32), "...")
  sym <- NA_character_
  if (!is.null(ann_row) && nrow(ann_row) > 0) {
    sym <- ann_row$ncbi_gene_name[1]
    if (is.na(sym) || sym == "" || sym == gene_id) sym <- NA_character_
  }
  if (!is.na(sym)) paste0(sym, " — ", fn) else paste0(gene_id, " — ", fn)
}

selected_meta <- cur[match(selected, cur$gene_id), c("gene_id","subcat","description")]
selected_meta$label <- vapply(seq_len(nrow(selected_meta)), function(i) {
  ar <- ann[ann$gene_id == selected_meta$gene_id[i], ]
  make_label(selected_meta$gene_id[i], selected_meta$description[i], ar)
}, character(1))

# Make labels unique (collisions on identical descriptions)
dup_idx <- which(duplicated(selected_meta$label))
if (length(dup_idx)) {
  selected_meta$label[dup_idx] <- paste0(
    sub(" —.*$", "", selected_meta$label[dup_idx]),
    " (", selected_meta$gene_id[dup_idx], ") —",
    sub("^[^—]+—", "", selected_meta$label[dup_idx])
  )
}

# ---- per-subcategory hierarchical clustering ----------------------
cluster_within <- function(genes_sc) {
  if (length(genes_sc) < 2) return(genes_sc)
  fc_sub <- fc[gene_id %in% genes_sc]
  fc_sub[, cell := paste(Sampling_Season, Maturity_stage,
                         Years_Since_Infection, sep="|")]
  mat <- dcast(fc_sub, gene_id ~ cell, value.var = "log2FC")
  rn <- mat$gene_id; mat <- as.matrix(mat[, -1])
  rownames(mat) <- rn; mat[is.na(mat)] <- 0
  d <- dist(mat, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  rn[hc$order]
}

ordered_genes <- character(0)
present_subcats <- intersect(SUBCAT_LEVELS, unique(selected_meta$subcat))
for (sc in present_subcats) {
  genes_sc <- selected_meta$gene_id[selected_meta$subcat == sc]
  ordered_genes <- c(ordered_genes, cluster_within(genes_sc))
}

# Build label vector (gene_id -> "Subcat > label" for axis_nested)
# axis_nested splits on a single character (default ">") so we use that separator.
selected_meta <- selected_meta[match(ordered_genes, selected_meta$gene_id), ]
selected_meta$nested_label <- paste0(selected_meta$subcat, " > ", selected_meta$label)

# ---- subset fc to selected genes ----------------------------------
fc_sel <- fc[gene_id %in% ordered_genes]
fc_sel[, gene_id := factor(gene_id, levels = rev(ordered_genes))]
fc_sel[, Maturity_stage := factor(Maturity_stage, levels = STG_LEVELS, labels = STG_LABELS)]
fc_sel[, cohort := factor(Years_Since_Infection,
                          levels = c(2016,2019,2020,2021,2022),
                          labels = c("'16","'19","'20","'21","'22"))]

# Map gene_id (factor) -> nested_label for axis_nested
y_lookup <- setNames(selected_meta$nested_label, selected_meta$gene_id)
fc_sel$y_nested <- factor(y_lookup[as.character(fc_sel$gene_id)],
                          levels = rev(y_lookup[ordered_genes]))

# ---- cell text overlay (|log2FC| >= 1.5) --------------------------
fc_sel[, cell_label := ifelse(!is.na(log2FC) & abs(log2FC) >= 1.5,
                              sprintf("%.1f", log2FC), "")]
fc_sel[, cell_text_color := ifelse(abs(log2FC) > 2, "white", "black")]

fc_keep <- fc_sel[!is.na(log2FC)]

# ---- main heatmap ------------------------------------------------
p <- ggplot(fc_keep, aes(x = cohort, y = y_nested, fill = log2FC)) +
  geom_tile(color = "grey90", linewidth = 0.1) +
  geom_text(aes(label = cell_label, color = cell_text_color),
            size = 2.4, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  facet_nested(. ~ Sampling_Season + Maturity_stage,
               nest_line = element_line(color = "grey50", linewidth = 0.4),
               strip = strip_nested(size = "variable"),
               scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish, name = expression(log[2]*FC)) +
  scale_y_discrete(guide = guide_axis_nested(delim = ">")) +
  scale_x_discrete(expand = expansion(0)) +
  labs(
    title = sprintf("S80C — Cell wall genes with consistent GRBV response (subset of S80B)  (n=%d)",
                    length(ordered_genes)),
    subtitle = sprintf("%d genes selected from 161 by |log2FC| consistency criteria.  See S80B for full landscape.",
                       length(ordered_genes)),
    caption = paste0(
      "Selection: |mean log2FC| >= 1 in >=4 of 8 (season x stage) cells  OR  ",
      sprintf("|mean log2FC| >= %.2f in >=%d of 8 cells (auto-tuned to land in the 20-40 window).  ",
              final_thresh, final_mc),
      "Cell value = mean log2(CPM+1) RB(+) at (season, stage, cohort) - ",
      "mean log2(CPM+1) RB(-) at (season, stage).  ",
      "Sub-category order top->bottom: ",
      paste(present_subcats, collapse = " -> "),
      ".  Rows clustered (Euclidean / Ward.D2) within sub-category.  ",
      "Cell text shown for |log2FC| >= 1.5 (white when |log2FC| > 2)."
    ),
    x = NULL, y = NULL
  ) +
  theme_grbv(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    ggh4x.axis.nesttext.y = element_text(size = 9, face = "bold", angle = 90, hjust = 0.5),
    ggh4x.axis.nestline.y = element_line(color = "grey40", linewidth = 0.6),
    strip.text.x = element_text(size = 8, face = "bold"),
    panel.spacing.x = unit(0.05, "lines"),
    legend.position = "right",
    plot.title    = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    plot.caption  = element_text(size = 7.5, color = "grey40", hjust = 0)
  )

ggsave(paste0(OUT_STEM, ".png"), p,
       width = 10, height = 9, units = "in", dpi = 600, limitsize = FALSE)
ggsave(paste0(OUT_STEM, ".svg"), p,
       width = 10, height = 9, units = "in", limitsize = FALSE)
cat("[S80C] saved\n")
