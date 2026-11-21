# Shared inputs and the MPhS projection wrapper. Sourced by steps 03-04.
.libPaths(c("Rlib", .libPaths()))
suppressPackageStartupMessages({library(MPhS); library(edgeR)})

audit <- read.csv("../simplified_reanalysis_2026-11/data/sample_audit.csv", colClasses = c(Sample = "character", vine = "character"))
counts <- read.delim("../../Raw_Data/counts_matrix.txt", check.names = FALSE, colClasses = c(Gene_ID = "character"))
cm <- as.matrix(counts[, -1]); rownames(cm) <- counts$Gene_ID
stopifnot(nrow(audit) == 327, setequal(audit$Sample, colnames(cm)), !anyDuplicated(audit$Sample))
cm <- cm[, audit$Sample]
STAGES <- c("Pre_veraison", "Veraison", "Post_veraison", "Harvest")
audit$stage_num <- match(audit$stage, STAGES) - 1
audit$healthy <- audit$infection == "RB(-)"
stopifnot(!anyNA(audit$stage_num), sum(audit$healthy) == 63, length(unique(audit$vine)) == 44)

gmap <- read.csv("data/mphs_gene_map.csv", colClasses = "character")
gmap <- gmap[gmap$passes_filter == "TRUE", ]

# Expression for the usable core genes, libraries x genes, columns named by VIT ID.
# Library sizes use all genes, so CPM is not distorted by the subset.
expr <- function(type = c("cpm", "tmm_logcpm"), genes = gmap) {
  type <- match.arg(type)
  y <- if (type == "cpm") cpm(cm) else cpm(calcNormFactors(DGEList(cm)), log = TRUE, prior.count = 2)
  x <- t(y[genes$row_id, , drop = FALSE]); colnames(x) <- genes$VIT_id
  x
}

# Project the given libraries (rows of x) in one call; returns integer MPhS stage per library.
project <- function(x, libs, scaling = "scale") {
  d <- as.data.frame(x[libs, , drop = FALSE], check.names = FALSE)
  d$lib <- libs; d$grp <- "all"
  out <- suppressWarnings(MPhStimepoints(d, strata_var = "grp", stage_var = "lib", scaling_type = scaling))
  stopifnot(identical(out$PCscores$lib, libs))
  setNames(out$PCscores$timepoint, libs)
}
