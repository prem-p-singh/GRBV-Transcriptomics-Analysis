#!/usr/bin/env Rscript
# 02. WGCNA module preservation — modules built in RB-, tested in RB+.
# Output: Zsummary per module + medianRank; PNG of Zsummary vs module size.

suppressPackageStartupMessages({
  library(WGCNA); library(data.table); library(optparse); library(ggplot2)
})
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

opts <- parse_args(OptionParser(option_list=list(
  make_option("--counts"), make_option("--meta"),
  make_option("--out_csv"), make_option("--out_plot"),
  make_option("--n_perm",    type="integer", default=100),
  make_option("--soft_power", type="integer", default=6),
  make_option("--cores",     type="integer", default=8)
)))

counts <- fread(opts$counts, data.table=FALSE); rownames(counts) <- counts[,1]; counts[,1] <- NULL
meta   <- fread(opts$meta,   data.table=FALSE)

expr <- as.matrix(counts)
# Retain genes with count >= 10 in >= 10% samples, then log2CPM
keep_g <- rowSums(expr >= 10) >= 0.10 * ncol(expr)
expr   <- expr[keep_g, ]
libsize <- colSums(expr)
lcpm <- log2(sweep(expr, 2, libsize, "/") * 1e6 + 1)

rb_neg <- meta$Sample[meta$Infection_Group == "RB(-)"]
rb_pos <- meta$Sample[meta$Infection_Group == "RB(+)"]
rb_neg <- intersect(rb_neg, colnames(lcpm))
rb_pos <- intersect(rb_pos, colnames(lcpm))
message(sprintf("RB-: %d   RB+: %d", length(rb_neg), length(rb_pos)))

# Top variable genes (use RB- variance as reference structure)
v_neg <- apply(lcpm[, rb_neg], 1, var)
keep_v <- names(sort(v_neg, decreasing=TRUE))[1:min(10000, sum(v_neg > 0))]

datRef <- t(lcpm[keep_v, rb_neg])
datTst <- t(lcpm[keep_v, rb_pos])

# Build modules in RB- reference
net <- blockwiseModules(datRef, power=opts$soft_power,
                        TOMType="signed", minModuleSize=30,
                        reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=FALSE, pamRespectsDendro=FALSE,
                        verbose=3, nThreads=opts$cores, maxBlockSize=20000)

multiExpr  <- list(ref=list(data=datRef), test=list(data=datTst))
multiColor <- list(ref=net$colors)

mp <- modulePreservation(multiExpr, multiColor,
                         referenceNetworks=1, nPermutations=opts$n_perm,
                         randomSeed=1, quickCor=0, verbose=3)

Z  <- mp$preservation$Z$ref$inColumnsAlsoPresentIn.test
obs<- mp$preservation$observed$ref$inColumnsAlsoPresentIn.test
out <- data.frame(module    = rownames(Z),
                  moduleSize= obs$moduleSize,
                  Zsummary  = Z$Zsummary.pres,
                  medianRank= obs$medianRank.pres)
out <- out[!out$module %in% c("gold","grey"), ]
out <- out[order(-out$Zsummary), ]

dir.create(dirname(opts$out_csv), recursive=TRUE, showWarnings=FALSE)
fwrite(out, opts$out_csv)

p <- ggplot(out, aes(moduleSize, Zsummary, label=module)) +
  geom_point(aes(color=module), size=3) +
  geom_text(vjust=-0.8, size=3) +
  geom_hline(yintercept=c(2,10), linetype="dashed", color=c("orange","darkgreen")) +
  scale_x_log10() + scale_color_identity() +
  labs(title="Module preservation (RB- built, RB+ tested)",
       x="Module size (log)", y="Zsummary.pres") + theme_bw()
ggsave(opts$out_plot, p, width=8, height=6, dpi=200)
message("Done.")
