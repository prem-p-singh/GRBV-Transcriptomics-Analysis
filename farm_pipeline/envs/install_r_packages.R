#!/usr/bin/env Rscript
# Install all R / Bioconductor dependencies into a project-local library.
# Called by setup_farm.sh after the LMOD R module is loaded.
lib <- Sys.getenv("GRBV_RLIBS", "~/grbv_extended/Rlibs")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)

cran <- c("optparse","data.table","ggplot2","pheatmap","igraph",
          "matrixStats","HIMA","BiocManager")
miss <- cran[!(cran %in% installed.packages(lib.loc = lib)[, "Package"])]
if (length(miss)) install.packages(miss, lib = lib)

bioc <- c("WGCNA","GENIE3","slingshot","tradeSeq","SingleCellExperiment",
          "DESeq2","edgeR","limma","impute")
miss <- bioc[!(bioc %in% installed.packages(lib.loc = lib)[, "Package"])]
if (length(miss)) BiocManager::install(miss, lib = lib, update = FALSE, ask = FALSE)

# DGCA has been archived off CRAN — install from archive tarball
if (!("DGCA" %in% installed.packages(lib.loc = lib)[, "Package"])) {
  install.packages("https://cran.r-project.org/src/contrib/Archive/DGCA/DGCA_1.0.3.tar.gz",
                   lib = lib, repos = NULL, type = "source")
}

cat("\n[R] installed into", lib, "\n")
print(installed.packages(lib.loc = lib)[, "Package"])
