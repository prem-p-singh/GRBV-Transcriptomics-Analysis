#!/usr/bin/env Rscript
# 19. Bayesian hierarchical model with vine random effects (brms).
#
# The frequentist repeated-measures voom found |rho_hat| <= 0.04 within-vine
# correlation, but for a flagship-journal submission a Bayesian random-intercept
# model gives ROPE / posterior credible intervals on the main infection effect
# that are more defensible than NHST. We fit a small panel of high-interest
# genes (the cross-pipeline robust mediators + the consensus biomarkers) and
# report posterior summaries.
#
# Genes tested:
#   - LOC100249249 (MBF1c)
#   - LOC100232974 (SPS1)
#   - LOC100266309 (uncharacterised, robust mediator)
#   - LOC100268011 (RNA-dep RNA polymerase 1, top biomarker)
#   - LOC100854143 (F-box/LRR top biomarker)
#   - LOC100257913 (glycine-rich cell wall, top DOWN biomarker)
#
# Outputs (results/19_bayesian/):
#   - bayes_posterior_summary.csv
#   - bayes_posterior_intervals_plot.png

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2); library(tidyr)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option("--counts",   type="character"),
  make_option("--meta",     type="character"),
  make_option("--out_dir",  type="character"),
  make_option("--n_iter",   type="integer", default=2000),
  make_option("--n_chains", type="integer", default=4)
)))
dir.create(opt$out_dir, recursive=TRUE, showWarnings=FALSE)

if (!requireNamespace("brms", quietly=TRUE)) {
  stop("install.packages('brms') required (run on FARM with Rlibs)")
}
library(brms)

# ---- load data ---------------------------------------------------------
counts <- as.matrix(read.delim(opt$counts, row.names=1, check.names=FALSE))
meta   <- read.delim(opt$meta, stringsAsFactors=FALSE)
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]
lib  <- colSums(counts); cpm <- t(t(counts) / lib) * 1e6
expr <- log2(cpm + 1)

genes <- c("LOC100249249", "LOC100232974", "LOC100266309",
           "LOC100268011", "LOC100854143", "LOC100257913")
genes <- intersect(genes, rownames(expr))
cat("Genes available:", length(genes), "\n"); print(genes)

mm <- meta[match(colnames(expr), meta$Sample), ]
mm$infection <- factor(mm$Infection_Group, levels=c("RB(-)","RB(+)"))
mm$season    <- factor(mm$Sampling_Season)
mm$stage     <- factor(mm$Maturity_stage,
                       levels=c("Pre_veraison","Veraison","Post_veraison","Harvest"))
if ("vine_id" %in% names(mm)) {
  mm$vine <- factor(mm$vine_id)
} else {
  # fallback: pseudo-vine = season+ysi+biorep cell
  mm$vine <- factor(paste(mm$Sampling_Season, mm$Years_Since_Infection,
                          mm$biorep_id, sep="_"))
}
cat("Vines:", nlevels(mm$vine), "\n")

# Weakly informative priors
priors <- c(
  prior(normal(0, 1),   class="b"),
  prior(normal(0, 5),   class="Intercept"),
  prior(student_t(3, 0, 1), class="sd"),
  prior(student_t(3, 0, 1), class="sigma")
)

results <- list()
for (g in genes) {
  cat("\n=== Fitting", g, "===\n")
  mm$y <- expr[g, mm$Sample]
  fit <- try(brm(
    formula = y ~ infection + season + stage + (1 | vine),
    data    = mm,
    family  = gaussian(),
    prior   = priors,
    chains  = opt$n_chains,
    iter    = opt$n_iter,
    cores   = min(4, opt$n_chains),
    refresh = 200,
    silent  = 2
  ), silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat("  fit failed for", g, "\n"); next
  }
  post <- as_draws_df(fit)
  inf_col <- grep("^b_infection", names(post), value=TRUE)[1]
  if (is.na(inf_col)) next
  draws <- post[[inf_col]]
  rope <- mean(abs(draws) < 0.1)         # 10% of log2-CPM ~ small effect
  ci   <- quantile(draws, c(0.025, 0.5, 0.975))
  results[[g]] <- data.frame(
    gene=g, n=nrow(mm),
    posterior_mean = mean(draws),
    ci_lo = ci[[1]], ci_med = ci[[2]], ci_hi = ci[[3]],
    p_in_rope = rope,
    rhat_max = max(rhat(fit), na.rm=TRUE),
    n_eff_min = min(neff_ratio(fit), na.rm=TRUE) * sum(post$.chain != 0)
  )
}
out <- bind_rows(results)
write.csv(out, file.path(opt$out_dir, "bayes_posterior_summary.csv"), row.names=FALSE)
print(out)

# Forest plot
if (nrow(out) > 0) {
  p <- ggplot(out, aes(x=posterior_mean, y=gene)) +
    geom_point() +
    geom_errorbarh(aes(xmin=ci_lo, xmax=ci_hi), height=0.2) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=c(-0.1, 0.1), linetype="dotted", color="grey50") +
    labs(title="Bayesian posterior intervals: infection effect on log2 CPM",
         subtitle="Vertical band = ROPE [-0.1, 0.1]",
         x="Posterior mean (95% CI)", y=NULL) + theme_bw()
  ggsave(file.path(opt$out_dir, "bayes_posterior_intervals_plot.png"), p,
         width=6, height=4, dpi=150)
}
cat("Done.\n")
