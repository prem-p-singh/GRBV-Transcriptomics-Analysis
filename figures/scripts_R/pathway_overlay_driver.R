# figures/scripts_R/pathway_overlay_driver.R
# Driver for the 6 canonical-KEGG pathway figures (S87-S92).
# See docs/superpowers/specs/2026-05-26-pathway-figures-design.md for full design.

suppressPackageStartupMessages({
  library(pathview)
  library(readr)
  library(dplyr)
  library(withr)
  library(XML)
})

# -----------------------------------------------------------------
# Constants
# -----------------------------------------------------------------
PATHWAY_DRIVER_VERSION <- "1.0.0"

# Column order for the 8-cell split-box: pure stage (2021 first, then 2022)
Q2X_COLS <- c(
  "2021_Pre_veraison", "2021_Veraison", "2021_Post_veraison", "2021_Harvest",
  "2022_Pre_veraison", "2022_Veraison", "2022_Post_veraison", "2022_Harvest"
)

# Theme-aligned colour palette for pathview (blue/white/red diverging)
PV_LOW  <- "#2166ac"
PV_MID  <- "#f7f7f7"
PV_HIGH <- "#b2182b"

# log2FC clip and bins for split-box rendering
PV_LIMIT <- 3
PV_BINS  <- 20

# Coverage threshold for WARNING (below this is "partial figure")
COVERAGE_WARN_THRESHOLD <- 0.70

# Directories
KGML_DIR  <- "figures/cache/kegg_xml"
TMP_DIR   <- "figures/supp/_pathview_tmp"
SUPP_DIR  <- "figures/supp"

# -----------------------------------------------------------------
# Helper: read one Q2x DESeq2 CSV and return (gene, logFC) for the matrix
# -----------------------------------------------------------------
.read_one_q2x <- function(de_dir, season, stage) {
  fname <- file.path(de_dir, sprintf("Q2x_%s_%s_deseq2.csv", season, stage))
  if (!file.exists(fname)) stop("Missing DESeq2 CSV: ", fname)
  df <- read_csv(fname, show_col_types = FALSE,
                 col_types = cols(gene     = col_character(),
                                  log2FC   = col_double(),
                                  logFC_shrunk = col_double(),
                                  .default = col_skip()))
  # Prefer shrunken; fall back to unshrunken only where shrunken is NA
  df$lfc <- ifelse(is.na(df$logFC_shrunk), df$log2FC, df$logFC_shrunk)
  df[, c("gene", "lfc")]
}

# -----------------------------------------------------------------
# read_q2x_matrix(): build the 14,604-row x 8-column log2FC matrix
# -----------------------------------------------------------------
read_q2x_matrix <- function(de_dir = "farm_pipeline/results/de_extras") {
  parts <- lapply(Q2X_COLS, function(col) {
    fields <- strsplit(col, "_", fixed = TRUE)[[1]]
    season <- fields[1]
    stage  <- paste(fields[-1], collapse = "_")
    one <- .read_one_q2x(de_dir, season, stage)
    setNames(one[, c("gene", "lfc")], c("gene", col))
  })
  # Inner-join all 8 by gene
  M_df <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "gene"), parts)
  # Filter to ^LOC\d+$ only
  M_df <- M_df[grepl("^LOC\\d+$", M_df$gene), , drop = FALSE]
  # Drop all-NA rows
  is_all_na <- apply(M_df[, Q2X_COLS], 1, function(x) all(is.na(x)))
  M_df <- M_df[!is_all_na, , drop = FALSE]
  # Strip LOC prefix; numeric Entrez rownames on the matrix
  M <- as.matrix(M_df[, Q2X_COLS])
  rownames(M) <- sub("^LOC", "", M_df$gene)
  M
}

# -----------------------------------------------------------------
# compute_mapping_coverage(): how many enzyme nodes have >=1 LOC mapped?
# Returns a list: n_ec, n_ec_mapped, percent.
# -----------------------------------------------------------------
compute_mapping_coverage <- function(pathway_id, entrez_ids,
                                     kgml_dir = KGML_DIR) {
  # Strip "vvi" prefix if present for download.kegg()
  pid <- sub("^vvi", "", pathway_id)

  # Download if not cached
  kgml_path <- file.path(kgml_dir, sprintf("vvi%s.xml", pid))
  if (!file.exists(kgml_path)) {
    invisible(download.kegg(pathway.id = pid, species = "vvi", kegg.dir = kgml_dir))
  }
  if (!file.exists(kgml_path)) {
    stop("KGML download failed for vvi", pid)
  }

  # Parse the KGML to extract gene entries (each entry has a space-separated list of vvi:NNNN IDs)
  kgml <- XML::xmlParse(kgml_path)
  entries <- XML::getNodeSet(kgml, "//entry[@type='gene']")
  ec_node_loci <- lapply(entries, function(e) {
    name_attr <- XML::xmlGetAttr(e, "name")
    if (is.null(name_attr)) return(character(0))
    ids <- strsplit(name_attr, " ", fixed = TRUE)[[1]]
    sub("^vvi:", "", ids)
  })

  n_ec <- length(ec_node_loci)
  mapped <- vapply(ec_node_loci, function(loci) any(loci %in% entrez_ids), logical(1))
  n_ec_mapped <- sum(mapped)
  percent <- if (n_ec > 0) 100 * n_ec_mapped / n_ec else 0

  list(n_ec = n_ec, n_ec_mapped = n_ec_mapped, percent = percent,
       per_ec_loci = ec_node_loci, per_ec_mapped = mapped)
}
