#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("GRBV_RLIBS","~/grbv_extended/Rlibs"), .libPaths()))
# 01. Gene regulatory network inference with GENIE3.
# TFs identified by keyword search on gene products in the NCBI GFF
# (PlantTFDB VIT_* IDs don't match counts_matrix LOC* IDs).

suppressPackageStartupMessages({
  library(GENIE3); library(data.table); library(optparse)
})

opts <- parse_args(OptionParser(option_list=list(
  make_option("--counts"),    make_option("--meta"),
  make_option("--tfs"),       make_option("--gff",     default="resources/genome/genome.gff"),
  make_option("--out_edges"), make_option("--out_hubs"),
  make_option("--n_trees",    type="integer", default=1000),
  make_option("--top_edges",  type="integer", default=100000),
  make_option("--top_var",    type="integer", default=5000),
  make_option("--cores",      type="integer", default=16)
)))

counts <- fread(opts$counts, data.table=FALSE)
rownames(counts) <- counts[,1]; counts[,1] <- NULL
expr <- as.matrix(counts)

libsize <- colSums(expr)
lcpm <- log2(sweep(expr, 2, libsize, "/") * 1e6 + 1)
v <- apply(lcpm, 1, var)
keep <- names(sort(v, decreasing=TRUE))[1:opts$top_var]
X <- lcpm[keep, , drop=FALSE]

# Identify TF genes from GFF product annotations
tf_patterns <- paste(c("WRKY","\\bbHLH\\b","\\bMYB\\b","\\bMYC","\\bAP2\\b","\\bERF\\b",
                       "\\bHSF\\b","heat shock factor","\\bbZIP\\b","\\bNAC\\b","\\bGATA\\b",
                       "\\bARF\\b","\\bARR\\b","\\bTCP\\b","\\bGRF\\b","\\bSPL\\b","\\bYABBY\\b",
                       "\\bCO-like\\b","\\bDof\\b","\\bE2F\\b","\\bTRIHELIX\\b","\\bSBP\\b",
                       "\\bHB\\b","homeobox","\\bLOB\\b","\\bBBR-BPC\\b","\\bCPP\\b",
                       "\\bTGA\\b","\\bTALE\\b","MADS","\\bSRS\\b","\\bB3\\b",
                       "transcription factor","transcription activator","transcription repressor"),
                     collapse="|")
tf_ids <- character()
gff <- readLines(opts$gff)
for (line in gff) {
  if (startsWith(line, "#")) next
  p <- strsplit(line, "\t", fixed=TRUE)[[1]]
  if (length(p) < 9 || p[3] != "gene") next
  attr <- p[9]
  m_name <- regmatches(attr, regexpr("gene=[^;]+", attr))
  if (!length(m_name)) m_name <- regmatches(attr, regexpr("Name=[^;]+", attr))
  if (!length(m_name)) next
  name <- sub("^(gene|Name)=", "", m_name)
  prod <- regmatches(attr, regexpr("product=[^;]+", attr))
  desc <- regmatches(attr, regexpr("description=[^;]+", attr))
  text <- paste(c(prod, desc), collapse=" ")
  if (grepl(tf_patterns, text, ignore.case=TRUE)) tf_ids <- c(tf_ids, name)
}
tf_ids <- unique(tf_ids)
message(sprintf("TFs identified from GFF: %d", length(tf_ids)))

regs <- intersect(tf_ids, rownames(X))
message(sprintf("Regulators in expression matrix: %d / %d TFs",
                length(regs), length(tf_ids)))
if (length(regs) < 20) stop("Too few regulators mapped — check GFF parse / gene IDs.")

set.seed(42)
w <- GENIE3(X, regulators=regs, nTrees=opts$n_trees, nCores=opts$cores, verbose=TRUE)
links <- getLinkList(w)
links <- head(links[order(-links$weight), ], opts$top_edges)

hubs <- aggregate(weight ~ regulatoryGene, data=links, FUN=sum)
hubs <- hubs[order(-hubs$weight), ]
colnames(hubs) <- c("TF","total_weight")
hubs$n_targets <- as.integer(table(links$regulatoryGene)[hubs$TF])

dir.create(dirname(opts$out_edges), recursive=TRUE, showWarnings=FALSE)
fwrite(links, opts$out_edges, sep="\t")
fwrite(hubs,  opts$out_hubs,  sep=",")
message("Done.")
