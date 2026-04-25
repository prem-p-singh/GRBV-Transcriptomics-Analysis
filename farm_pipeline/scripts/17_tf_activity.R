#!/usr/bin/env Rscript
# 17. TF activity inference via decoupleR + GENIE3 regulons.
#
# Expression-level TF DE (267 TFs, rule 09) doesn't equal regulatory activity.
# Here we infer per-sample TF activity scores from each TF's downstream
# regulon (its inferred target set + edge weights from GENIE3 / rule 01),
# using decoupleR's univariate linear-model method (ulm) and weighted-sum (wsum).
#
# This addresses agenda gap #5 and de-risks the HSFA2 VIGS proposal: if
# HSFA2 expression is flat but its regulon-inferred activity is dampened
# in RB(+), that strengthens the passive-failure mechanism.
#
# Inputs:
#   --counts    : counts_matrix.txt
#   --meta      : sample_metadata.txt (or v2)
#   --grn_edges : results/01_grn/grn_edges.tsv.gz (TF, target, weight)
#   --out_dir   : results/17_tf_activity/
#
# Outputs:
#   - tf_activity_ulm.csv         (samples x TFs activity score)
#   - tf_activity_de.csv          (per-TF Welch-t RB+ vs RB-, BH FDR)
#   - top_tf_activity_heatmap.png

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
})

opt_list <- list(
  make_option("--counts",    type="character"),
  make_option("--meta",      type="character"),
  make_option("--grn_edges", type="character"),
  make_option("--out_dir",   type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$out_dir, recursive=TRUE, showWarnings=FALSE)

# Try decoupleR (newer); fall back to a hand-rolled univariate linear-model
# implementation if decoupleR not available.
have_decoupleR <- suppressMessages(requireNamespace("decoupleR", quietly=TRUE))
cat("decoupleR available:", have_decoupleR, "\n")

# ---- load data --------------------------------------------------------
counts <- as.matrix(read.delim(opt$counts, row.names=1, check.names=FALSE))
meta   <- read.delim(opt$meta, stringsAsFactors=FALSE)
grn    <- read.delim(gzfile(opt$grn_edges), header=TRUE, stringsAsFactors=FALSE)
# expected columns: TF (or regulator/regulatoryGene), target/targetGene, weight
names(grn) <- tolower(names(grn))
rename_map <- c("regulator"="tf", "regulatorygene"="tf",
                "targetgene"="target", "score"="weight")
for (k in names(rename_map)) {
  if (k %in% names(grn) && !rename_map[[k]] %in% names(grn)) {
    names(grn)[names(grn)==k] <- rename_map[[k]]
  }
}
stopifnot(all(c("tf","target","weight") %in% names(grn)))

# Keep only top edges per TF (match GENIE3 typical thresholding)
grn <- grn %>% group_by(tf) %>% slice_max(order_by=weight, n=200) %>% ungroup()
cat("Regulons retained:", length(unique(grn$tf)), "TFs;",
    nrow(grn), "edges\n")

# Filter low-count genes consistent with main analysis
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]

# log2 CPM
lib  <- colSums(counts)
cpm  <- t(t(counts) / lib) * 1e6
expr <- log2(cpm + 1)

# Drop TFs with too few targets in expression matrix
grn <- grn %>% filter(target %in% rownames(expr))
n_per_tf <- grn %>% count(tf, name="n_targets")
grn <- grn %>% inner_join(n_per_tf %>% filter(n_targets >= 5), by="tf")
cat("Final regulons after expression filter:", length(unique(grn$tf)), "TFs\n")

# ---- compute activity ------------------------------------------------
if (have_decoupleR) {
  library(decoupleR)
  net <- grn %>% mutate(mor = sign(weight)) %>%
    select(source = tf, target, weight, mor)
  res_ulm <- run_ulm(mat = expr, net = net,
                     .source = "source", .target = "target", .mor = "mor",
                     minsize = 5)
  act_ulm <- res_ulm %>% select(condition, source, score) %>%
    pivot_wider(names_from = source, values_from = score)
  write.csv(act_ulm, file.path(opt$out_dir, "tf_activity_ulm.csv"),
            row.names = FALSE)
  cat("decoupleR ulm activity matrix written.\n")
} else {
  # Fallback: per-sample weighted mean of target expression (z-scored).
  cat("Falling back to weighted-mean activity (no decoupleR).\n")
  expr_z <- t(scale(t(expr)))
  expr_z[is.na(expr_z)] <- 0
  tfs <- unique(grn$tf)
  act <- matrix(NA_real_, nrow=length(tfs), ncol=ncol(expr_z),
                dimnames=list(tfs, colnames(expr_z)))
  for (tf in tfs) {
    sub <- grn[grn$tf == tf, , drop = FALSE]
    if (nrow(sub) < 5) next
    sub <- sub[sub$target %in% rownames(expr_z), , drop = FALSE]
    if (nrow(sub) < 5) next
    act[tf, ] <- colSums(expr_z[sub$target, , drop=FALSE] * sub$weight) / sum(abs(sub$weight))
  }
  act <- act[apply(act, 1, function(x) !all(is.na(x))), ]
  write.csv(t(act), file.path(opt$out_dir, "tf_activity_ulm.csv"))
  cat("Fallback activity matrix written.\n")
}

# ---- per-TF DE: RB(+) vs RB(-) on activity ---------------------------
act <- read.csv(file.path(opt$out_dir, "tf_activity_ulm.csv"),
                stringsAsFactors=FALSE, check.names=FALSE)
# First column is sample id (header is "" when written via write.csv with rownames)
samples <- act[[1]]                # positional index avoids empty-name issue
mat <- as.matrix(act[, -1, drop=FALSE])
rownames(mat) <- samples
mat <- t(mat)  # TFs x samples
cat("Activity matrix dim (TFs x samples):", dim(mat), "\n")

mm <- meta %>% filter(Sample %in% colnames(mat))
mat <- mat[, mm$Sample]
group <- factor(mm$Infection_Group, levels=c("RB(-)", "RB(+)"))
res_de <- data.frame(tf=rownames(mat), stringsAsFactors=FALSE)
res_de$mean_pos <- rowMeans(mat[, group=="RB(+)"], na.rm=TRUE)
res_de$mean_neg <- rowMeans(mat[, group=="RB(-)"], na.rm=TRUE)
res_de$diff     <- res_de$mean_pos - res_de$mean_neg
res_de$t_stat   <- NA_real_; res_de$p_value <- NA_real_
for (i in seq_len(nrow(mat))) {
  x <- mat[i, group=="RB(+)"]; y <- mat[i, group=="RB(-)"]
  if (sum(!is.na(x))<3 || sum(!is.na(y))<3) next
  tt <- try(t.test(x, y), silent=TRUE)
  if (inherits(tt,"try-error")) next
  res_de$t_stat[i] <- tt$statistic
  res_de$p_value[i] <- tt$p.value
}
res_de$padj <- p.adjust(res_de$p_value, method="BH")
res_de <- res_de %>% arrange(p_value)
write.csv(res_de, file.path(opt$out_dir, "tf_activity_de.csv"), row.names=FALSE)
cat("TF DE rows:", nrow(res_de), "; top 10:\n"); print(head(res_de, 10))

# ---- heatmap of top TFs ----------------------------------------------
suppressPackageStartupMessages({ library(ggplot2) })
# Filter to TFs with non-NaN p-values; cap at 30
res_keep <- res_de %>% filter(!is.na(p_value))
top_tf <- head(res_keep$tf, min(30, nrow(res_keep)))
mat_top <- mat[top_tf, , drop=FALSE]
mat_top <- mat_top[, !apply(mat_top, 2, function(x) all(is.na(x))), drop=FALSE]
cat("Heatmap will plot", nrow(mat_top), "TFs x", ncol(mat_top), "samples\n")
if (nrow(mat_top) >= 2 && ncol(mat_top) >= 2) {
  df <- as.data.frame(mat_top) %>%
    tibble::rownames_to_column("tf") %>%
    pivot_longer(-tf, names_to="Sample", values_to="activity") %>%
    inner_join(mm %>% select(Sample, Infection_Group, Maturity_stage, Sampling_Season),
               by="Sample") %>%
    mutate(label = paste(Sampling_Season, Maturity_stage, Infection_Group, sep="|"))
  p <- ggplot(df, aes(x=Sample, y=tf, fill=activity)) +
    geom_tile() +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B") +
    facet_grid(. ~ Infection_Group + Sampling_Season, scales="free_x", space="free_x") +
    theme_minimal(base_size=8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text = element_text(size=7)) +
    labs(title=paste0("Top ", nrow(mat_top), " TFs by activity DE (RB+ vs RB-)"), y=NULL, x=NULL)
  ggsave(file.path(opt$out_dir, "top_tf_activity_heatmap.png"), p,
         width=10, height=7, dpi=150)
} else {
  cat("Insufficient TFs/samples for heatmap; skipping plot.\n")
}
cat("Done.\n")
