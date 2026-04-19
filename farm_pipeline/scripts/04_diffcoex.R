#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("GRBV_RLIBS","~/grbv_extended/Rlibs"), .libPaths()))
# 04. Differential co-expression (DGCA) — RB+ vs RB-.
# Outputs: rewired edges with pValDiff + z-scores, and gene-level rewiring counts.

suppressPackageStartupMessages({
  library(DGCA); library(data.table); library(optparse)
})

opts <- parse_args(OptionParser(option_list=list(
  make_option("--counts"), make_option("--meta"),
  make_option("--out_edges"), make_option("--out_hubs"),
  make_option("--top_var", type="integer", default=3000),
  make_option("--n_perm",  type="integer", default=100)
)))

counts <- fread(opts$counts, data.table=FALSE); rownames(counts) <- counts[,1]; counts[,1] <- NULL
meta   <- fread(opts$meta,   data.table=FALSE)

expr <- as.matrix(counts)
keep_g <- rowSums(expr >= 10) >= 0.10 * ncol(expr)
expr   <- expr[keep_g, ]
libsize <- colSums(expr)
lcpm <- log2(sweep(expr, 2, libsize, "/") * 1e6 + 1)

v <- apply(lcpm, 1, var)
keep <- names(sort(v, decreasing=TRUE))[1:opts$top_var]
M <- lcpm[keep, ]

# Design matrix (required by DGCA): columns = conditions, rows = samples
samples <- meta$Sample[meta$Sample %in% colnames(M)]
M <- M[, samples]
cond <- ifelse(meta$Infection_Group[match(samples, meta$Sample)] == "RB(+)", "infected", "control")
design <- model.matrix(~0 + cond)
colnames(design) <- levels(factor(cond))

res <- ddcorAll(inputMat=M, design=design,
                compare=c("control","infected"),
                adjust="BH", nPerm=opts$n_perm,
                classify=TRUE, corrType="pearson",
                heatmapPlot=FALSE, nPairs="all")

res <- as.data.table(res)
# keep most rewired
res <- res[order(pValDiff)]
top  <- head(res, 50000)

dir.create(dirname(opts$out_edges), recursive=TRUE, showWarnings=FALSE)
fwrite(top, opts$out_edges)

# Gene-level rewiring: count significant pairs per gene
sig <- res[pValDiff < 0.05]
hub_counts <- sort(table(c(sig$Gene1, sig$Gene2)), decreasing=TRUE)
hubs <- data.table(Gene=names(hub_counts), rewired_partners=as.integer(hub_counts))
fwrite(head(hubs, 500), opts$out_hubs)
message("Done.")
