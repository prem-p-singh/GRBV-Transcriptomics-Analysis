#!/usr/bin/env Rscript
# T5. Matched-random-set null for the headline gene groups. Each real gene gets a
# random non-branch, non-viral partner from the same bin (healthy trajectory
# direction, dynamic-range decile, mean-logCPM decile). Random sets lack the
# co-regulation of real pathways, so this null is tighter than a true pathway
# null: a pass here is necessary, not sufficient.
ROOT_INIT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
.libPaths(c(file.path(ROOT_INIT, "Reference_Resources/Rlibs_phase2"), .libPaths()))
source(file.path(ROOT_INIT, "Current_Analysis/simplified_reanalysis_2026-11/00_common.R"))
suppressPackageStartupMessages({library(parallel); library(ggplot2)})
B <- 1999L
TIGHT <- Sys.getenv("GRBV_TIGHT") == "1"  # sensitivity: partners from the 25 genes nearest in dynamic range
SUFFIX <- if (TIGHT) "_tight" else ""
N_NEAREST <- 25L
MAX_SET_MARKERS <- 100L  # Era 7 branch marker limit
GROUPS <- c(Photosystem_light_harvesting = "BRANCH:58", sHSP = "BRANCH:50",
            Cytosolic_translation = "BRANCH:55", Chloroplast_translation = "BRANCH:56")

controls <- sort(unique(meta[infection == "neg"]$vine_id))
full_ref <- make_ref(setNames(rep(1, 8), controls))
genes <- rownames(E)
info <- data.table(gene = genes, dynamic = full_ref$dynamic, mean_logcpm = rowMeans(E),
                   direction = sign(full_ref$pooled[, 4] - full_ref$pooled[, 1]))
decile <- function(x, ref) findInterval(x, unique(quantile(ref, seq(.1, .9, .1))))
info[, `:=`(dyn_bin = decile(dynamic, dynamic[gene %in% background]),
            mean_bin = decile(mean_logcpm, mean_logcpm[gene %in% background]))]
info[, in_pool := gene %in% background]
pool <- info[in_pool == TRUE]

# Partner candidates by bin, widened only when the exact bin has too few genes.
candidates <- function(row, need) {
  if (TIGHT) {
    p <- pool[direction == row$direction]
    return(list(ids = p[order(abs(dynamic - row$dynamic))]$gene[seq_len(N_NEAREST)], level = "nearest_dynamic"))
  }
  a <- pool[direction == row$direction & dyn_bin == row$dyn_bin & mean_bin == row$mean_bin]$gene
  if (length(a) >= need) return(list(ids = a, level = "direction+dynamic+mean"))
  b <- pool[direction == row$direction & dyn_bin == row$dyn_bin]$gene
  if (length(b) >= need) return(list(ids = b, level = "direction+dynamic"))
  list(ids = pool[direction == row$direction]$gene, level = "direction")
}

fold_ids <- sort(unique(meta$fold))
refs <- lapply(fold_ids, function(f) {
  drop <- unique(meta[infection == "neg" & fold == f]$vine_id)
  make_ref(setNames(as.numeric(!controls %in% drop), controls))
})
fold_libs <- lapply(fold_ids, function(f) which(meta$fold == f))
pred <- fread(file.path(OUT, "data/dev_score_predictions.csv"))[match(meta$Sample, as.character(Sample))]
stopifnot(!anyNA(pred$new_score))

set_offset <- function(ids) {
  s <- rep(NA_real_, ncol(E))
  for (k in seq_along(refs)) s[fold_libs[[k]]] <- predict_position(refs[[k]], pick_markers(refs[[k]], ids, MAX_SET_MARKERS), fold_libs[[k]])
  if (anyNA(s)) return(rep(NA_real_, 8))
  displacement(s - pred$new_score)
}

adj <- fread(file.path(OUT, "data/adjusted_dream_cell_contrasts.csv"))[model == "adjusted_residual"]
adj[, cell := paste(stage, season, sep = "__")]
effect <- as.matrix(dcast(adj, gene ~ cell, value.var = "log2FC"), rownames = "gene")[genes, CELL_KEYS]

run_group <- function(name) {
  real <- intersect(mem[set_id == GROUPS[[name]]]$gene, genes)
  rows <- info[match(real, gene)]
  rows[, bin := if (TIGHT) gene else paste(direction, dyn_bin, mean_bin)]
  need <- rows[, .N, by = bin]
  cand <- lapply(split(rows, rows$bin), function(r) candidates(r[1], nrow(r)))
  level <- rows[, .(n = .N, level = cand[[bin[1]]]$level), by = bin]
  draw <- function(b) {
    set.seed(SEED + 200000L + b + 3000L * match(name, names(GROUPS)))
    ids <- unlist(lapply(names(cand), function(k) sample(cand[[k]]$ids, sum(rows$bin == k))), use.names = FALSE)
    if (anyDuplicated(ids)) ids <- unique(ids)  # widened pools can overlap; the set shrinks slightly
    c(colMeans(effect[ids, , drop = FALSE]), set_offset(ids),
      dynamic = mean(info$dynamic[match(ids, info$gene)]), mean_logcpm = mean(info$mean_logcpm[match(ids, info$gene)]), size = length(ids))
  }
  null <- do.call(rbind, mclapply(seq_len(B), draw, mc.cores = CORES, mc.preschedule = TRUE))
  stopifnot(is.numeric(null), nrow(null) == B)
  colnames(null)[1:16] <- c(paste0("adj_", CELL_KEYS), paste0("off_", CELL_KEYS))
  list(name = name, real = c(colMeans(effect[real, , drop = FALSE]), set_offset(real)), null = null, level = level,
       quality = data.table(group = name, n_genes = length(real), real_dynamic = mean(rows$dynamic),
                            null_dynamic = mean(null[, "dynamic"]), real_mean_logcpm = mean(rows$mean_logcpm),
                            null_mean_logcpm = mean(null[, "mean_logcpm"]), null_size = mean(null[, "size"]),
                            pct_genes_exact_bin = 100 * level[level %in% c("direction+dynamic+mean", "nearest_dynamic"), sum(n)] / length(real),
                            n_null_sets_scored = sum(is.finite(null[, 9]))))
}
results <- lapply(names(GROUPS), run_group)

tab <- rbindlist(lapply(results, function(r) rbindlist(lapply(1:16, function(j) {
  x <- r$null[, j]; x <- x[is.finite(x)]; v <- r$real[j]
  data.table(group = r$name, statistic = if (j <= 8) "mean_adjusted_effect" else "score_offset",
             cell = CELL_KEYS[(j - 1) %% 8 + 1], real_value = v, null_mean = mean(x), null_sd = sd(x),
             null_q025 = quantile(x, .025), null_q975 = quantile(x, .975), percentile = 100 * mean(x < v),
             p_two_sided = min(1, 2 * min(1 + sum(x <= v), 1 + sum(x >= v)) / (length(x) + 1)),
             outside_central_95 = v < quantile(x, .025) | v > quantile(x, .975), n_null = length(x))
}))))
tab[, c("stage", "season") := tstrsplit(cell, "__", fixed = TRUE)]
fwrite(tab, file.path(OUT, sprintf("data/matched_null_percentiles%s.csv", SUFFIX)))
quality <- rbindlist(lapply(results, `[[`, "quality"))
fwrite(quality, file.path(OUT, sprintf("data/matched_null_quality%s.csv", SUFFIX)))
print(quality)
print(tab[stage == "Harvest", .(group, statistic, season, real_value, null_q025, null_q975, percentile, p_two_sided, outside_central_95)])

hist_data <- rbindlist(lapply(results, function(r) rbindlist(lapply(c(4, 8, 12, 16), function(j) {
  data.table(group = r$name, statistic = if (j <= 8) "Mean adjusted effect" else "Score offset",
             season = sub(".*__", "", CELL_KEYS[(j - 1) %% 8 + 1]), value = r$null[, j])
}))))
marks <- tab[stage == "Harvest"]
marks[, statistic := fifelse(statistic == "score_offset", "Score offset", "Mean adjusted effect")]
p <- ggplot(hist_data, aes(value)) + geom_histogram(bins = 40, fill = "#9FB8B5", colour = NA) +
  geom_vline(data = marks, aes(xintercept = real_value), colour = "#7A1F3D", linewidth = .6) +
  facet_wrap(~ group + statistic + season, scales = "free", ncol = 4) +
  labs(x = "Harvest value (line: real group; histogram: matched random sets)", y = "Random sets") +
  theme_minimal(base_size = 8)
ggsave(file.path(OUT, sprintf("figures/SI_matched_null%s.png", SUFFIX)), p, width = 11, height = 7.5, dpi = 300)
cat("T5 COMPLETE\n")
