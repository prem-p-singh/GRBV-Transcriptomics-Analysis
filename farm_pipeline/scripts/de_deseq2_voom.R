#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("GRBV_RLIBS","~/grbv_extended/Rlibs"), .libPaths()))
# Fix #1: refit differential expression with DESeq2 AND limma-voom,
# in parallel to the published Welch-t results. Addresses the reviewer
# concern that Welch t on log2CPM is not the field standard for bulk
# RNA-seq.
#
# Covers all 6 core research questions:
#   Q1  infection duration (RB+ vs RB- across YSI bins)
#   Q2  maturity stage (per-stage RB+ vs RB-)
#   Q3  seasonal effects (2021 vs 2022 in RB+)
#   Q4  environmental modulation (Season x Infection)
#   Q5  baseline environment (RB- seasons)
#   Q6  pooled RB+ vs RB-
#
# Writes per-gene DE tables in a common schema:
#   gene, log2FC, logFC_shrunk (DESeq2 apeglm), pvalue, padj, method,
#   comparison
# plus a summary CSV with DEG counts per comparison per method.

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma)
  library(data.table); library(optparse)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--counts", default = "data/counts_matrix.txt"),
  make_option("--meta",   default = "data/sample_metadata.txt"),
  make_option("--out_dir",default = "results/de_redo"),
  make_option("--padj_thresh", type = "double", default = 0.05),
  make_option("--lfc_thresh",  type = "double", default = 1.0),
  make_option("--cores", type = "integer", default = 8)
)))

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

cat("[de] loading counts + metadata\n")
counts <- fread(opts$counts, data.table = FALSE)
rownames(counts) <- counts[, 1]; counts[, 1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"

meta <- fread(opts$meta, data.table = FALSE)
rownames(meta) <- meta$Sample
meta <- meta[colnames(counts), ]
# DESeq2 / limma require "clean" factor levels (letters + numbers only).
# Rename RB(-) / RB(+) -> neg / pos; remember the mapping for output labels.
meta$Infection_Group <- ifelse(meta$Infection_Group == "RB(-)", "neg", "pos")
meta$Infection_Group <- factor(meta$Infection_Group, levels = c("neg", "pos"))
meta$Sampling_Season  <- factor(meta$Sampling_Season)
meta$Maturity_stage   <- factor(meta$Maturity_stage,
                                levels = c("Pre_veraison", "Veraison",
                                           "Post_veraison", "Harvest"))

# Gene-level filter: >=10 counts in >=10% samples (same as main paper)
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]
cat(sprintf("[de] %d genes x %d samples after filter\n",
            nrow(counts), ncol(counts)))

run_DESeq2 <- function(counts_sub, meta_sub, design_formula, contrast, comp_name) {
  dds <- DESeqDataSetFromMatrix(counts_sub, colData = meta_sub, design = design_formula)
  dds <- DESeq(dds, parallel = FALSE, quiet = TRUE)
  res <- results(dds, contrast = contrast, alpha = opts$padj_thresh)
  # Shrink with apeglm if the coefficient name is unambiguous, else ashr
  res_s <- tryCatch(
    lfcShrink(dds, contrast = contrast, type = "ashr"),
    error = function(e) res
  )
  data.frame(
    gene      = rownames(res),
    log2FC    = res$log2FoldChange,
    logFC_shrunk = res_s$log2FoldChange,
    pvalue    = res$pvalue,
    padj      = res$padj,
    method    = "DESeq2",
    comparison = comp_name,
    stringsAsFactors = FALSE
  )
}

run_voom <- function(counts_sub, meta_sub, design_mat, contrast_vec, comp_name) {
  dge <- DGEList(counts = counts_sub)
  dge <- calcNormFactors(dge, method = "TMM")
  v   <- voom(dge, design_mat, plot = FALSE)
  fit <- lmFit(v, design_mat)
  fit <- contrasts.fit(fit, contrast_vec)
  fit <- eBayes(fit)
  tt  <- topTable(fit, number = Inf, sort.by = "none")
  data.frame(
    gene      = rownames(tt),
    log2FC    = tt$logFC,
    logFC_shrunk = tt$logFC,
    pvalue    = tt$P.Value,
    padj      = tt$adj.P.Val,
    method    = "limma-voom",
    comparison = comp_name,
    stringsAsFactors = FALSE
  )
}

write_results <- function(df, comp_name, method) {
  out <- file.path(opts$out_dir, sprintf("%s_%s.csv", comp_name, method))
  fwrite(df, out)
}

# -----------------------------------------------------------------------------
# Q6 : pooled RB+ vs RB- (headline comparison)
# -----------------------------------------------------------------------------
cat("\n[de] ==== Q6: pooled RB+ vs RB- ====\n")
design_q6 <- ~ Maturity_stage + Sampling_Season + Infection_Group
ds_q6 <- run_DESeq2(counts, meta, design_q6,
                    c("Infection_Group", "pos", "neg"), "Q6_pooled")
write_results(ds_q6, "Q6_pooled", "deseq2")

mm <- model.matrix(design_q6, data = meta)
colnames(mm) <- make.names(colnames(mm))
contr_q6 <- rep(0, ncol(mm)); contr_q6[grep("Infection_Grouppos", colnames(mm))] <- 1
vm_q6 <- run_voom(counts, meta, mm, contr_q6, "Q6_pooled")
write_results(vm_q6, "Q6_pooled", "voom")

# -----------------------------------------------------------------------------
# Q2 : per-stage RB+ vs RB-
# -----------------------------------------------------------------------------
cat("\n[de] ==== Q2: per-stage RB+ vs RB- ====\n")
for (st in levels(meta$Maturity_stage)) {
  ss <- which(meta$Maturity_stage == st)
  if (length(unique(meta$Infection_Group[ss])) < 2 ||
      min(table(meta$Infection_Group[ss])) < 3) next
  cnt_s <- counts[, ss]; meta_s <- droplevels(meta[ss, ])
  comp <- sprintf("Q2_%s", st)
  cat(sprintf("   %s: n=%d\n", comp, length(ss)))
  ds <- run_DESeq2(cnt_s, meta_s,
                   ~ Sampling_Season + Infection_Group,
                   c("Infection_Group","pos","neg"), comp)
  write_results(ds, comp, "deseq2")
  mm_s <- model.matrix(~ Sampling_Season + Infection_Group, data = meta_s)
  colnames(mm_s) <- make.names(colnames(mm_s))
  c_s <- rep(0, ncol(mm_s))
  c_s[grep("Infection_Grouppos", colnames(mm_s))] <- 1
  vm <- run_voom(cnt_s, meta_s, mm_s, c_s, comp)
  write_results(vm, comp, "voom")
}

# -----------------------------------------------------------------------------
# Q3 : 2021 vs 2022 in RB+
# -----------------------------------------------------------------------------
cat("\n[de] ==== Q3: seasonal effects in RB+ ====\n")
ss3 <- which(meta$Infection_Group == "pos")
meta_3 <- droplevels(meta[ss3, ]); cnt_3 <- counts[, ss3]
ds3 <- run_DESeq2(cnt_3, meta_3,
                  ~ Maturity_stage + Sampling_Season,
                  c("Sampling_Season","2022","2021"), "Q3_seasonRBpos")
write_results(ds3, "Q3_seasonRBpos", "deseq2")
mm3 <- model.matrix(~ Maturity_stage + Sampling_Season, data = meta_3)
colnames(mm3) <- make.names(colnames(mm3))
c3 <- rep(0, ncol(mm3)); c3[grep("Sampling_Season2022", colnames(mm3))] <- 1
vm3 <- run_voom(cnt_3, meta_3, mm3, c3, "Q3_seasonRBpos")
write_results(vm3, "Q3_seasonRBpos", "voom")

# -----------------------------------------------------------------------------
# Q5 : baseline seasonal effects in RB- controls
# -----------------------------------------------------------------------------
cat("\n[de] ==== Q5: seasonal effects in RB- ====\n")
ss5 <- which(meta$Infection_Group == "neg")
if (length(ss5) > 10 && length(unique(meta$Sampling_Season[ss5])) > 1) {
  meta_5 <- droplevels(meta[ss5, ]); cnt_5 <- counts[, ss5]
  ds5 <- run_DESeq2(cnt_5, meta_5,
                    ~ Maturity_stage + Sampling_Season,
                    c("Sampling_Season","2022","2021"), "Q5_seasonRBneg")
  write_results(ds5, "Q5_seasonRBneg", "deseq2")
  mm5 <- model.matrix(~ Maturity_stage + Sampling_Season, data = meta_5)
  colnames(mm5) <- make.names(colnames(mm5))
  c5 <- rep(0, ncol(mm5)); c5[grep("Sampling_Season2022", colnames(mm5))] <- 1
  vm5 <- run_voom(cnt_5, meta_5, mm5, c5, "Q5_seasonRBneg")
  write_results(vm5, "Q5_seasonRBneg", "voom")
}

# -----------------------------------------------------------------------------
# Q4 : Season x Infection interaction (pooled)
# -----------------------------------------------------------------------------
cat("\n[de] ==== Q4: Season x Infection interaction ====\n")
dds4 <- DESeqDataSetFromMatrix(counts, colData = meta,
                               design = ~ Maturity_stage + Sampling_Season * Infection_Group)
dds4 <- DESeq(dds4, test = "LRT",
              reduced = ~ Maturity_stage + Sampling_Season + Infection_Group,
              parallel = FALSE, quiet = TRUE)
res4 <- results(dds4, alpha = opts$padj_thresh)
q4 <- data.frame(
  gene     = rownames(res4),
  log2FC   = res4$log2FoldChange,
  logFC_shrunk = res4$log2FoldChange,
  pvalue   = res4$pvalue,
  padj     = res4$padj,
  method   = "DESeq2_LRT",
  comparison = "Q4_seasonXinfection",
  stringsAsFactors = FALSE
)
write_results(q4, "Q4_seasonXinfection", "deseq2")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n[de] ==== summary ====\n")
files <- list.files(opts$out_dir, pattern = "\\.csv$", full.names = TRUE)
summ <- rbindlist(lapply(files, function(f) {
  d <- fread(f)
  data.table(
    file = basename(f),
    comparison = d$comparison[1],
    method = d$method[1],
    n_tested = nrow(d),
    n_padj05 = sum(!is.na(d$padj) & d$padj < opts$padj_thresh),
    n_padj05_lfc = sum(!is.na(d$padj) & d$padj < opts$padj_thresh &
                       abs(d$log2FC) >= opts$lfc_thresh)
  )
}))
fwrite(summ, file.path(opts$out_dir, "SUMMARY.csv"))
print(summ)
cat("\n[de] done; outputs in", opts$out_dir, "\n")
