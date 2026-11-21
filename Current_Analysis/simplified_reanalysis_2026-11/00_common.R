# Shared code for the simplified reanalysis. Scoring logic follows Era 7
# (04_developmental_test.R), with two
# changes: reference centroids are pooled over seasons, and scores may run past
# the first and last stage. Vines enter a reference through weights, so the same
# code serves cross-fitting, the nested bootstrap and vine-label permutation.
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
OUT  <- Sys.getenv("GRBV_OUT", file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11"))
PRIOR <- file.path(ROOT, "Current_Analysis/prior_era_inputs")
STAGES <- c("Pre_veraison", "Veraison", "Post_veraison", "Harvest")
SEASONS <- c("2021", "2022")
N_MARKERS <- 1000L
MIN_DYNAMIC <- 0.5    # Era 7 marker filter
MIN_SCALE <- 0.25     # Era 7 scale floor
SEED <- 20261102L
CORES <- as.integer(Sys.getenv("GRBV_CORES", "15"))

obj <- readRDS(file.path(PRIOR, "Era7_model_objects.rds"))
E <- obj$v$E
meta <- copy(obj$meta)
meta[, `:=`(Sample = as.character(Sample), season = as.character(season),
            infection = as.character(infection), vine_id = as.character(vine_id))]
meta[, stage_number := match(as.character(Maturity_stage), STAGES) - 1L]
stopifnot(identical(meta$Sample, colnames(E)), ncol(E) == 327, nrow(E) == 14604)
folds_era7 <- fread(file.path(PRIOR, "Era7_developmental_folds.csv"), select = c("vine_id", "fold"))
folds_era7 <- unique(folds_era7[, vine_id := as.character(vine_id)])
meta[, fold := folds_era7$fold[match(vine_id, folds_era7$vine_id)]]
stopifnot(!anyNA(meta$fold))

mem <- fread(file.path(PRIOR, "Era7_pathway_membership.csv"))
branch_genes <- unique(mem[startsWith(set_id, "BRANCH:")]$gene)
viral_genes <- grep("^N761_", rownames(E), value = TRUE)
background <- setdiff(rownames(E), c(branch_genes, viral_genes))

# Reference from vine weights (named by vine). pool = TRUE averages the season
# centroids at each stage; pool = FALSE keeps one centroid set per season.
make_ref <- function(w, pool = TRUE, allowed = rep(TRUE, ncol(E))) {
  w <- w[w > 0]
  train <- which(meta$vine_id %in% names(w) & allowed)
  wl <- unname(w[meta$vine_id[train]])
  cell <- paste(meta$season[train], meta$stage_number[train])
  fitted <- E[, train, drop = FALSE]
  cent <- list()
  for (cl in unique(cell)) {
    k <- which(cell == cl)
    m <- drop(E[, train[k], drop = FALSE] %*% wl[k]) / sum(wl[k])
    fitted[, k] <- m
    cent[[cl]] <- m
  }
  scale <- pmax(sqrt(drop(((E[, train, drop = FALSE] - fitted)^2) %*% wl) / sum(wl)), MIN_SCALE)
  season_cent <- function(y) {
    sapply(0:3, function(s) { key <- paste(y, s); if (key %in% names(cent)) cent[[key]] else rep(NA_real_, nrow(E)) })
  }
  ys <- sort(unique(meta$season[train]))
  by_season <- setNames(lapply(ys, function(y) { m <- season_cent(y); rownames(m) <- rownames(E); m }), ys)
  stage_mean <- sapply(0:3, function(s) {
    v <- sapply(by_season, function(m) m[, s + 1L])
    if (is.null(dim(v))) v else rowMeans(v, na.rm = TRUE)
  })
  rownames(stage_mean) <- rownames(E)
  if (anyNA(stage_mean)) stop("A stage has no reference library")
  dynamic <- pmax(stage_mean[, 1], stage_mean[, 2], stage_mean[, 3], stage_mean[, 4]) -
    pmin(stage_mean[, 1], stage_mean[, 2], stage_mean[, 3], stage_mean[, 4])
  strength <- rowMeans((stage_mean - rowMeans(stage_mean))^2) / scale^2
  names(scale) <- rownames(E)
  list(pooled = stage_mean, by_season = by_season, scale = scale, dynamic = dynamic,
       strength = strength, pool = pool)
}

pick_markers <- function(ref, ids, n = N_MARKERS) {
  ids <- intersect(ids, names(ref$strength))
  ids <- ids[ref$dynamic[ids] >= MIN_DYNAMIC & is.finite(ref$strength[ids])]
  head(ids[order(-ref$strength[ids], ids)], n)
}

# Position on the piecewise-linear healthy path. With extrapolate = TRUE, t may
# fall below 0 on the first segment and above 1 on the last; the interior clamp
# between adjacent segments stays.
project <- function(means, obs, extrapolate = TRUE) {
  errors <- pos <- matrix(0, 3, ncol(obs))
  for (s in 1:3) {
    a <- means[, s]
    delta <- means[, s + 1] - a
    centered <- obs - a
    t <- colSums(centered * delta) / sum(delta^2)
    lo <- if (extrapolate && s == 1) -Inf else 0
    hi <- if (extrapolate && s == 3) Inf else 1
    t <- pmin(hi, pmax(lo, t))
    pos[s, ] <- s - 1 + t
    errors[s, ] <- colMeans((centered - delta %o% t)^2)
  }
  pos[cbind(max.col(-t(errors), ties.method = "first"), seq_len(ncol(obs)))]
}

predict_position <- function(ref, ids, idx, extrapolate = TRUE) {
  if (length(ids) < 5) return(rep(NA_real_, length(idx)))
  scale <- ref$scale[ids]
  if (ref$pool) return(project(ref$pooled[ids, , drop = FALSE] / scale, E[ids, idx, drop = FALSE] / scale, extrapolate))
  ans <- numeric(length(idx))
  for (y in unique(meta$season[idx])) {
    loc <- which(meta$season[idx] == y)
    m <- if (y %in% names(ref$by_season)) ref$by_season[[y]] else ref$by_season[[1]]
    ans[loc] <- project(m[ids, , drop = FALSE] / scale, E[ids, idx[loc], drop = FALSE] / scale, extrapolate)
  }
  ans
}

# Expression expected at a position on the pooled path, extrapolating linearly
# past the end stages.
expected_expression <- function(ref, ids, position) {
  s <- pmin(3L, pmax(1L, floor(position) + 1L))
  t <- position - (s - 1)
  m <- ref$pooled[ids, , drop = FALSE]
  m[, s, drop = FALSE] * (rep(1, length(ids)) %o% (1 - t)) + m[, s + 1L, drop = FALSE] * (rep(1, length(ids)) %o% t)
}

# Infected minus healthy mean score per cell. labels: named by vine, "neg"/"pos".
# wt: library weights (vine multiplicity in a bootstrap draw); 0 drops a library.
CELL_KEYS <- as.vector(outer(STAGES, SEASONS, paste, sep = "__"))
CELL_OF_LIB <- meta[, paste(Maturity_stage, season, sep = "__")]
displacement <- function(score, labels = NULL, wt = rep(1, length(score))) {
  g <- if (is.null(labels)) meta$infection else unname(labels[meta$vine_id])
  wmean <- function(i) sum(score[i] * wt[i]) / sum(wt[i])
  vapply(CELL_KEYS, function(k) wmean(which(CELL_OF_LIB == k & g == "pos")) -
           wmean(which(CELL_OF_LIB == k & g == "neg")), 0)
}
