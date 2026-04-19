#!/usr/bin/env Rscript
# 06. Pseudotime trajectory on RB+ samples, anchored by Years_Since_Infection.
# Uses slingshot on PCA; tests temporally variable genes with tradeSeq GAMs.

suppressPackageStartupMessages({
  library(data.table); library(optparse); library(ggplot2)
  library(SingleCellExperiment); library(slingshot); library(tradeSeq)
  library(BiocParallel); library(mgcv)
})

opts <- parse_args(OptionParser(option_list=list(
  make_option("--counts"), make_option("--meta"),
  make_option("--out_genes"), make_option("--out_plot"),
  make_option("--top_var", type="integer", default=3000),
  make_option("--n_knots", type="integer", default=6),
  make_option("--cores",   type="integer", default=4)
)))

counts <- fread(opts$counts, data.table=FALSE); rownames(counts) <- counts[,1]; counts[,1] <- NULL
meta   <- fread(opts$meta,   data.table=FALSE)

expr <- as.matrix(counts)
keep_g <- rowSums(expr >= 10) >= 0.10 * ncol(expr)
expr <- expr[keep_g, ]
libsize <- colSums(expr)
lcpm <- log2(sweep(expr, 2, libsize, "/") * 1e6 + 1)

rb_pos <- meta$Sample[meta$Infection_Group == "RB(+)" & meta$Sample %in% colnames(lcpm)]
ysi    <- meta$Years_Since_Infection[match(rb_pos, meta$Sample)]
stage  <- meta$Maturity_stage[match(rb_pos, meta$Sample)]

X <- lcpm[, rb_pos]
v <- apply(X, 1, var)
top <- names(sort(v, decreasing=TRUE))[1:opts$top_var]
Xt <- X[top, ]

pc <- prcomp(t(Xt), scale.=TRUE)
rd <- pc$x[, 1:4]

sce <- SingleCellExperiment(assays=list(counts=expr[top, rb_pos],
                                        logcounts=Xt),
                            reducedDims=list(PCA=rd))
colData(sce)$ysi   <- ysi
colData(sce)$stage <- stage

# Cluster along PCA then fit slingshot anchored by YSI correlation with pseudotime
km <- kmeans(rd, centers=5, nstart=25)
sce$cluster <- paste0("c", km$cluster)
sce <- slingshot(sce, clusterLabels="cluster", reducedDim="PCA")

pt <- slingPseudotime(sce)[,1]
# Orient so positive pseudotime correlates with YSI
if (cor(pt, ysi, use="complete.obs") < 0) pt <- max(pt, na.rm=TRUE) - pt

# tradeSeq: associationTest for genes changing along pseudotime
sce$pt <- pt
cw <- matrix(1, nrow=length(pt), ncol=1); colnames(cw) <- "lin1"
BPPARAM <- MulticoreParam(opts$cores)
set.seed(1)
sce2 <- fitGAM(counts=assay(sce,"counts"), pseudotime=pt, cellWeights=cw,
               nknots=opts$n_knots, verbose=TRUE, parallel=TRUE, BPPARAM=BPPARAM)
at  <- associationTest(sce2)
at$gene <- rownames(at)
at <- at[order(at$pvalue), ]

dir.create(dirname(opts$out_genes), recursive=TRUE, showWarnings=FALSE)
fwrite(at, opts$out_genes)

df <- data.frame(PC1=rd[,1], PC2=rd[,2], pt=pt, ysi=ysi, stage=stage)
p <- ggplot(df, aes(PC1, PC2, color=pt)) +
  geom_point(size=2) + scale_color_viridis_c() +
  labs(title="RB+ trajectory (slingshot)", color="Pseudotime") + theme_bw()
ggsave(opts$out_plot, p, width=7, height=5, dpi=200)
message("Done.")
