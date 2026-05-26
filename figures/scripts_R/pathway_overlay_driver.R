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

# -----------------------------------------------------------------
# write_paralog_auxiliary(): write per-EC paralog inventory CSV.
# For each enzyme node in the KGML, lists the LOC genes mapped and
# the mean log2FC across paralogs for each of the 8 Q2x conditions.
# -----------------------------------------------------------------
write_paralog_auxiliary <- function(pathway_id, M, out_path,
                                    kgml_dir = KGML_DIR) {
  cov <- compute_mapping_coverage(pathway_id, rownames(M), kgml_dir)
  rows <- vector("list", cov$n_ec)
  for (i in seq_len(cov$n_ec)) {
    loci <- cov$per_ec_loci[[i]]
    mapped_loci <- intersect(loci, rownames(M))
    if (length(mapped_loci) == 0) {
      mean_fc <- rep(NA_real_, length(Q2X_COLS))
    } else {
      sub_mat <- M[mapped_loci, , drop = FALSE]
      mean_fc <- colMeans(sub_mat, na.rm = TRUE)
      mean_fc[is.nan(mean_fc)] <- NA_real_
    }
    rows[[i]] <- c(
      list(pathway_id = pathway_id,
           ec_index   = i,
           n_loci     = length(mapped_loci),
           loci       = paste(paste0("LOC", mapped_loci), collapse = ";")),
      setNames(as.list(mean_fc), Q2X_COLS)
    )
  }
  out_df <- dplyr::bind_rows(rows)
  write_csv(out_df, out_path)
  nrow(out_df)
}

# -----------------------------------------------------------------
# render_kegg_pathway(): full pipeline for one KEGG pathway.
#   pathway_id : e.g. "vvi00940"
#   supp_id    : e.g. 87 (becomes S87_*.png)
#   slug       : e.g. "phenylpropanoid_kegg"
#   M          : Q2x log2FC matrix from read_q2x_matrix()
# Side effects:
#   - writes figures/supp/S<supp_id>_<slug>.png
#   - writes figures/supp/S<supp_id>_<slug>.pdf (graphviz fallback)
#   - writes figures/supp/S<supp_id>_paralog_counts.csv
#   - writes figures/supp/S<supp_id>_mapping_coverage.txt
# Returns invisibly: list(png_path, pdf_path, coverage_pct, n_ec_mapped, n_ec)
# -----------------------------------------------------------------
render_kegg_pathway <- function(pathway_id, supp_id, slug, M) {
  pid <- sub("^vvi", "", pathway_id)
  dir.create(TMP_DIR,  recursive = TRUE, showWarnings = FALSE)
  dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(KGML_DIR, recursive = TRUE, showWarnings = FALSE)

  # 1. Coverage report
  cov <- compute_mapping_coverage(pathway_id, rownames(M))
  message(sprintf("[%s] %s: %d of %d EC nodes mapped (%.1f%%)",
                  pathway_id, slug, cov$n_ec_mapped, cov$n_ec, cov$percent))
  if (cov$percent < COVERAGE_WARN_THRESHOLD * 100) {
    warning(sprintf("Coverage %.1f%% below threshold %.0f%% for %s",
                    cov$percent, COVERAGE_WARN_THRESHOLD * 100, pathway_id))
  }

  # 2. Paralog auxiliary CSV
  aux_path <- file.path(SUPP_DIR, sprintf("S%d_paralog_counts.csv", supp_id))
  n_rows <- write_paralog_auxiliary(pathway_id, M, aux_path)

  # 3. Coverage text file
  cov_path <- file.path(SUPP_DIR, sprintf("S%d_mapping_coverage.txt", supp_id))
  writeLines(sprintf("%s %s: %d of %d EC nodes mapped (%.1f%%)",
                     pathway_id, slug, cov$n_ec_mapped, cov$n_ec, cov$percent),
             cov_path)

  # 4. Native pathview render (canonical KEGG layout PNG)
  png_native <- NULL
  withr::with_dir(TMP_DIR, {
    out_suffix_native <- "grbv_native"
    pv_result <- tryCatch(
      pathview::pathview(
        gene.data    = M,
        pathway.id   = pid,
        species      = "vvi",
        gene.idtype  = "ENTREZ",
        multi.state  = TRUE,
        same.layer   = FALSE,
        node.sum     = "mean",
        low          = list(gene = PV_LOW),
        mid          = list(gene = PV_MID),
        high         = list(gene = PV_HIGH),
        limit        = list(gene = PV_LIMIT),
        bins         = list(gene = PV_BINS),
        kegg.native  = TRUE,
        out.suffix   = out_suffix_native,
        kegg.dir     = file.path("..", "..", "cache", "kegg_xml")
      ),
      error = function(e) { message("Native pathview render failed: ", e$message); NULL }
    )
    cand <- sprintf("vvi%s.%s.multi.png", pid, out_suffix_native)
    if (file.exists(cand)) png_native <<- cand
  })

  if (is.null(png_native)) {
    stop("Pathview native render produced no PNG for ", pathway_id)
  }

  # 5. Graphviz vector pass (PDF for paper figures)
  pdf_graphviz <- NULL
  withr::with_dir(TMP_DIR, {
    out_suffix_gv <- "grbv_graphviz"
    pv_result2 <- tryCatch(
      pathview::pathview(
        gene.data    = M,
        pathway.id   = pid,
        species      = "vvi",
        gene.idtype  = "ENTREZ",
        multi.state  = TRUE,
        same.layer   = FALSE,
        node.sum     = "mean",
        low          = list(gene = PV_LOW),
        mid          = list(gene = PV_MID),
        high         = list(gene = PV_HIGH),
        limit        = list(gene = PV_LIMIT),
        bins         = list(gene = PV_BINS),
        kegg.native  = FALSE,
        out.suffix   = out_suffix_gv,
        kegg.dir     = file.path("..", "..", "cache", "kegg_xml")
      ),
      error = function(e) { message("Graphviz pathview render failed: ", e$message); NULL }
    )
    cand <- sprintf("vvi%s.%s.multi.pdf", pid, out_suffix_gv)
    if (file.exists(cand)) pdf_graphviz <<- cand
  })

  # 6. Move outputs to figures/supp/ with S<n>_<slug> naming
  final_png <- file.path(SUPP_DIR, sprintf("S%d_%s.png", supp_id, slug))
  file.rename(file.path(TMP_DIR, png_native), final_png)

  final_pdf <- NULL
  if (!is.null(pdf_graphviz)) {
    final_pdf <- file.path(SUPP_DIR, sprintf("S%d_%s.pdf", supp_id, slug))
    file.rename(file.path(TMP_DIR, pdf_graphviz), final_pdf)
  }

  # 7. Try to convert PDF -> SVG if pdf2svg is on PATH
  final_svg <- NULL
  if (!is.null(final_pdf) && nchar(Sys.which("pdf2svg")) > 0) {
    final_svg <- sub("\\.pdf$", ".svg", final_pdf)
    system2("pdf2svg", c(shQuote(final_pdf), shQuote(final_svg)))
    if (!file.exists(final_svg)) final_svg <- NULL
  }

  invisible(list(
    png_path     = final_png,
    pdf_path     = final_pdf,
    svg_path     = final_svg,
    coverage_pct = cov$percent,
    n_ec_mapped  = cov$n_ec_mapped,
    n_ec         = cov$n_ec,
    aux_path     = aux_path
  ))
}
