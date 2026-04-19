#!/usr/bin/env Rscript
# 01. Gene regulatory network inference with GENIE3.
# Regulators = PlantTFDB V. vinifera TFs (matched to LOC IDs in counts matrix).

suppressPackageStartupMessages({
  library(GENIE3); library(data.table); library(optparse)
})

opts <- parse_args(OptionParser(option_list=list(
  make_option("--counts"),    make_option("--meta"),
  make_option("--tfs"),       make_option("--out_edges"),
  make_option("--out_hubs"),  make_option("--n_trees",  type="integer", default=1000),
  make_option("--top_edges",  type="integer", default=100000),
  make_option("--top_var",    type="integer", default=5000),
  make_option("--cores",      type="integer", default=16)
)))

counts <- fread(opts$counts, data.table=FALSE)
rownames(counts) <- counts[,1]; counts[,1] <- NULL
expr <- as.matrix(counts)

# log2(CPM+1), keep top-variable genes
libsize <- colSums(expr)
lcpm <- log2(sweep(expr, 2, libsize, "/") * 1e6 + 1)
v <- apply(lcpm, 1, var)
keep <- names(sort(v, decreasing=TRUE))[1:opts$top_var]
X <- lcpm[keep, , drop=FALSE]

# TF list → regulators = intersection with `keep`
tf_raw <- fread(opts$tfs, header=TRUE, fill=TRUE, data.table=FALSE)
tf_ids <- unique(unlist(tf_raw[, intersect(c("Gene_ID","Gene","Symbol","TF_ID"), colnames(tf_raw))]))
tf_ids <- tf_ids[nchar(tf_ids) > 0]
regs   <- intersect(tf_ids, rownames(X))
message(sprintf("Regulators in expression matrix: %d / %d TFs", length(regs), length(tf_ids)))
if (length(regs) < 20) stop("Too few regulators mapped — check TF list ID format.")

set.seed(42)
w <- GENIE3(X, regulators=regs, nTrees=opts$n_trees, nCores=opts$cores, verbose=TRUE)
links <- getLinkList(w)
links <- head(links[order(-links$weight), ], opts$top_edges)

# Hub TF ranking: sum of outgoing edge weights
hubs <- aggregate(weight ~ regulatoryGene, data=links, FUN=sum)
hubs <- hubs[order(-hubs$weight), ]
colnames(hubs) <- c("TF","total_weight")
hubs$n_targets <- as.integer(table(links$regulatoryGene)[hubs$TF])

dir.create(dirname(opts$out_edges), recursive=TRUE, showWarnings=FALSE)
fwrite(links, opts$out_edges, sep="\t")
fwrite(hubs,  opts$out_hubs,  sep=",")
message("Done.")
