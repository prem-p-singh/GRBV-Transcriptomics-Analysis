#!/usr/bin/env Rscript
# One vine-label permutation (T7). Runs the same draw as 07b_joint_permutation.R:
# same seed, same dream model, same voom weights. Writes one file per task.
suppressPackageStartupMessages({library(BiocParallel); library(data.table); library(limma); library(variancePartition)})
b <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
stopifnot(!is.na(b), !is.na(cores))
path <- sprintf("perm_draws/draw_%03d.csv", b)
if (file.exists(path)) quit(save = "no")  # a requeued task that already finished
x <- readRDS("perm_inputs.rds")
m <- x$m; vobj <- x$vobj
make_condition <- function(g) factor(paste(g, x$stage_codes[as.character(m$Maturity_stage)], m$season, sep = "_"))
form <- ~ 0 + condition + (1 | vine_id)
BP <- MulticoreParam(cores, progressbar = FALSE)
L <- makeContrastsDream(form, m, contrasts = x$contrasts)
two_season <- sort(as.character(as.data.table(m)[, .(n = uniqueN(season)), by = vine_id][n == 2]$vine_id))
stopifnot(length(two_season) == 40)
joint_counts <- function(fit) {
  fdr <- lapply(x$joint, function(cf) {
    tt <- variancePartition::topTable(fit, coef = cf, number = Inf, sort.by = "none", adjust.method = "none")
    p.adjust(tt$P.Value, "BH")
  })
  c(sum(fdr[[1]] < .05, na.rm = TRUE), sum(fdr[[2]] < .05, na.rm = TRUE), sum(fdr[[2]] < .05 & fdr[[1]] >= .05, na.rm = TRUE))
}
t0 <- Sys.time()
if (b == 0L) {  # draw 0: true labels, the reproducibility gate
  pseudo <- as.character(unique(m$vine_id[m$infection == "neg"]))
  m2 <- m
} else {
  set.seed(x$seed + 300000L + b)
  pseudo <- sample(two_season, 8)
  m2 <- m
  m2$condition <- make_condition(ifelse(as.character(m$vine_id) %in% pseudo, "neg", "pos"))
  if (nlevels(m2$condition) != 16) stop("Draw ", b, " lost a cell")
  m2$condition <- factor(m2$condition, levels = levels(m$condition))
}
fit <- variancePartition::eBayes(dream(vobj, form, m2, L = L, ddf = "Satterthwaite", BPPARAM = BP, quiet = TRUE))
n_fail <- length(attr(fit, "errors")) + length(fit$errors)
s <- joint_counts(fit)
row <- data.table(draw = b, pseudo_controls = paste(sort(pseudo), collapse = ";"),
                  n_real_controls_as_pseudo = sum(pseudo %in% as.character(m$vine_id[m$infection == "neg"])),
                  infection_stage = s[1], infection_stage_season = s[2], only_season_specific = s[3],
                  genes_fitted = nrow(fit), failures = n_fail,
                  seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
fwrite(row, paste0(path, ".tmp")); file.rename(paste0(path, ".tmp"), path)
print(row)
