#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("GRBV_RLIBS","~/grbv_extended/Rlibs"), .libPaths()))
# Fix #1b: re-run DE with a *repeated-measures* model (limma-voom + duplicate
# Correlation on vine_id). Paired with de_deseq2_voom, this is the reviewer-proof
# version of Q6 pooled and Q2 per-stage comparisons.
#
# Design:  ~ Maturity_stage + Sampling_Season + Infection_Group
# Blocking: vine_id (limma's duplicateCorrelation)
# Same thresholds: padj < 0.05 AND |log2FC| >= 1.

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma)
  library(data.table); library(optparse)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--counts", default = "data/counts_matrix.txt"),
  make_option("--meta",   default = "data/sample_metadata_v2.txt"),
  make_option("--out_dir",default = "results/de_repeated_measures"),
  make_option("--padj_thresh", type="double", default=0.05),
  make_option("--lfc_thresh",  type="double", default=1.0)
)))
dir.create(opts$out_dir, recursive=TRUE, showWarnings=FALSE)

cat("[de-rm] loading\n")
counts <- fread(opts$counts, data.table=FALSE)
rownames(counts) <- counts[,1]; counts[,1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"

meta <- fread(opts$meta, data.table=FALSE)
rownames(meta) <- meta$Sample
meta <- meta[colnames(counts), ]
meta$Infection_Group <- ifelse(meta$Infection_Group == "RB(-)", "neg", "pos")
meta$Infection_Group <- factor(meta$Infection_Group, levels=c("neg","pos"))
meta$Sampling_Season <- factor(meta$Sampling_Season)
meta$Maturity_stage  <- factor(meta$Maturity_stage,
  levels=c("Pre_veraison","Veraison","Post_veraison","Harvest"))

cat(sprintf("[de-rm] %d samples; %d unique vines\n",
            nrow(meta), length(unique(meta$vine_id))))

keep <- rowSums(counts >= 10) >= 0.10*ncol(counts)
counts <- counts[keep,]
cat(sprintf("[de-rm] %d genes x %d samples after filter\n",
            nrow(counts), ncol(counts)))

run_voom_rm <- function(counts_sub, meta_sub, design_formula, contrast_name, comp_name) {
  # voom with duplicateCorrelation blocking on vine_id
  dge <- DGEList(counts = counts_sub)
  dge <- calcNormFactors(dge, method="TMM")
  design_mat <- model.matrix(design_formula, data=meta_sub)
  colnames(design_mat) <- make.names(colnames(design_mat))
  v <- voom(dge, design_mat, plot=FALSE)
  corfit <- duplicateCorrelation(v, design_mat, block=meta_sub$vine_id)
  cat(sprintf("  [%s] consensus within-vine correlation = %.3f\n",
              comp_name, corfit$consensus.correlation))
  v2 <- voom(dge, design_mat, block=meta_sub$vine_id,
             correlation=corfit$consensus.correlation, plot=FALSE)
  corfit2 <- duplicateCorrelation(v2, design_mat, block=meta_sub$vine_id)
  fit <- lmFit(v2, design_mat, block=meta_sub$vine_id,
               correlation=corfit2$consensus.correlation)
  # Contrast: coefficient of interest
  coef_idx <- grep(contrast_name, colnames(design_mat))
  if (length(coef_idx) == 0) stop("no contrast matched ", contrast_name)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=coef_idx[1], number=Inf, sort.by="none")
  data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    pvalue = tt$P.Value,
    padj   = tt$adj.P.Val,
    method = "limma-voom-RM",
    comparison = comp_name,
    within_vine_cor = corfit2$consensus.correlation,
    stringsAsFactors=FALSE
  )
}

# Q6 pooled
cat("\n[de-rm] Q6 pooled RB+ vs RB- (repeated-measures)\n")
q6 <- run_voom_rm(counts, meta,
                  ~ Maturity_stage + Sampling_Season + Infection_Group,
                  "Infection_Grouppos", "Q6_pooled_RM")
fwrite(q6, file.path(opts$out_dir, "Q6_pooled_RM.csv"))

# Q2 per-stage
cat("\n[de-rm] Q2 per-stage (repeated-measures)\n")
for (st in levels(meta$Maturity_stage)) {
  ss <- which(meta$Maturity_stage == st)
  meta_s <- droplevels(meta[ss, ])
  if (length(unique(meta_s$Infection_Group)) < 2 ||
      min(table(meta_s$Infection_Group)) < 3) next
  # Need at least 2 vines per group for duplicateCorrelation
  if (length(unique(meta_s$vine_id)) < 4) { cat("   skip", st, "too few vines\n"); next }
  comp <- sprintf("Q2_%s_RM", st)
  cat(sprintf("   %s: n=%d samples, %d vines\n",
              comp, length(ss), length(unique(meta_s$vine_id))))
  cnt_s <- counts[, ss]
  res <- tryCatch(
    run_voom_rm(cnt_s, meta_s, ~ Sampling_Season + Infection_Group,
                "Infection_Grouppos", comp),
    error = function(e) { cat("     ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(res)) fwrite(res, file.path(opts$out_dir, paste0(comp, ".csv")))
}

cat("\n[de-rm] summary\n")
files <- list.files(opts$out_dir, pattern="\\.csv$", full.names=TRUE)
summ <- rbindlist(lapply(files, function(f) {
  d <- fread(f)
  data.table(file=basename(f), comparison=d$comparison[1],
             n_tested=nrow(d),
             within_vine_cor=round(d$within_vine_cor[1], 3),
             n_padj05=sum(!is.na(d$padj) & d$padj < opts$padj_thresh),
             n_padj05_lfc=sum(!is.na(d$padj) & d$padj < opts$padj_thresh &
                              abs(d$log2FC) >= opts$lfc_thresh))
}))
fwrite(summ, file.path(opts$out_dir, "SUMMARY.csv"))
print(summ)
cat("\n[de-rm] done\n")
