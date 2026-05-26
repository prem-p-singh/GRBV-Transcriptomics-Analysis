# figures/scripts_R/test_pathway_driver.R
# Smoke tests for pathway_overlay_driver.R helpers. Run via:
#   Rscript figures/scripts_R/test_pathway_driver.R
#
# All assertions use stopifnot(); script aborts on first failure.

source("figures/scripts_R/pathway_overlay_driver.R")

cat("Testing read_q2x_matrix() ...\n")
M <- read_q2x_matrix("farm_pipeline/results/de_extras")

# Shape checks
stopifnot("matrix should have 8 columns (4 stages x 2 seasons)" = ncol(M) == 8)
stopifnot("matrix should have >10000 rows after filter"          = nrow(M) > 10000)

# Column order: 2021 first, then 2022; within each, Pre -> Vér -> Post -> Harvest
expected_cols <- c("2021_Pre_veraison", "2021_Veraison", "2021_Post_veraison", "2021_Harvest",
                   "2022_Pre_veraison", "2022_Veraison", "2022_Post_veraison", "2022_Harvest")
stopifnot("column order must be pure stage" = identical(colnames(M), expected_cols))

# All row names should be numeric Entrez IDs (LOC prefix stripped)
stopifnot("all rownames should match ^[0-9]+$" = all(grepl("^\\d+$", rownames(M))))

# No row should be all-NA (those were dropped)
all_na_rows <- apply(M, 1, function(x) all(is.na(x)))
stopifnot("no all-NA rows allowed" = sum(all_na_rows) == 0)

# Values are log2FC in a sensible range. Most shrunken estimates sit within ±10,
# but very-low-baseline genes (one group near zero counts) can produce shrunken
# log2FC up to ~22 in this dataset (verified: 2 of 114,248 cells in 2022_Pre_veraison).
# Threshold of 25 still catches corruption (un-shrunk extremes, mis-parsed scientific
# notation, integer overflow) while accepting real DESeq2 output.
finite_vals <- M[is.finite(M)]
stopifnot("log2FC magnitudes look sane (|x| < 25)" = max(abs(finite_vals)) < 25)

cat("read_q2x_matrix(): PASS\n")
cat(sprintf("  shape: %d x %d, range: [%.2f, %.2f]\n",
            nrow(M), ncol(M), min(finite_vals), max(finite_vals)))
