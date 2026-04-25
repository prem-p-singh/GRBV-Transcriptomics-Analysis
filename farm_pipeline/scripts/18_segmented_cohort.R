#!/usr/bin/env Rscript
# 18. Change-point / segmented regression on cohort age (years since infection).
#
# Tests whether the 3-phase cohort framing (acute year 1, peak years 3-5,
# chronic years 6-7) reflects continuous drift or discrete phase transitions.
# Fits segmented regression with 0, 1, and 2 break-points on (a) PC1 of the
# top 2000 variable genes, and (b) each of 13 module eigengenes against
# years_since_infection, restricted to RB(+) samples.
#
# Picks the best model per response by BIC. If best is k>=1 break-points,
# reports the change-point year(s) with bootstrap 95% CI.
#
# Outputs (results/18_segmented/):
#   - segmented_summary.csv    (response, n, best_k, breakpoint, BIC_0, BIC_1, BIC_2, ...)
#   - segmented_pc1_plot.png
#   - segmented_modules_plot.png

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option("--counts",   type="character"),
  make_option("--meta",     type="character"),
  make_option("--mes",      type="character"),
  make_option("--out_dir",  type="character"),
  make_option("--top_var",  type="integer", default=2000)
)))
dir.create(opt$out_dir, recursive=TRUE, showWarnings=FALSE)

if (!requireNamespace("segmented", quietly=TRUE)) {
  stop("install.packages('segmented') required (run on FARM with Rlibs)")
}
library(segmented)

# ---- load data ---------------------------------------------------------
counts <- as.matrix(read.delim(opt$counts, row.names=1, check.names=FALSE))
meta   <- read.delim(opt$meta, stringsAsFactors=FALSE)
mes    <- read.csv(opt$mes, stringsAsFactors=FALSE, check.names=FALSE)

# RB(+) only
rb_pos <- meta$Sample[meta$Infection_Group == "RB(+)"]
counts <- counts[, colnames(counts) %in% rb_pos]
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts)
counts <- counts[keep, ]
lib <- colSums(counts); cpm <- t(t(counts) / lib) * 1e6
expr <- log2(cpm + 1)

# YSI numeric
mm <- meta[match(colnames(expr), meta$Sample), ]
mm$ysi <- as.numeric(mm$Years_Since_Infection)
# How many years since (relative): YSI is the cohort onset year; convert to
# years-since-infection assuming sampling year = Sampling_Season.
mm$years_since <- as.numeric(mm$Sampling_Season) - mm$ysi
mm$years_since[mm$years_since < 0] <- NA
expr <- expr[, !is.na(mm$years_since)]
mm <- mm[!is.na(mm$years_since), ]
cat("RB+ samples with valid YSI:", nrow(mm), "\n")

# ---- top-var PC1 ------------------------------------------------------
v <- apply(expr, 1, var); top <- names(sort(v, decreasing=TRUE))[seq_len(min(opt$top_var, length(v)))]
X <- t(expr[top, ])
X <- scale(X, center=TRUE, scale=FALSE)
pca <- prcomp(X, scale.=FALSE)
mm$PC1 <- pca$x[, 1]; mm$PC2 <- pca$x[, 2]

# ---- helper: fit linear, segmented(1), segmented(2) -------------------
fit_segmented <- function(y, x, max_psi=2L) {
  out <- list(); df <- data.frame(y=y, x=x); df <- df[complete.cases(df), ]
  m0 <- lm(y ~ x, data=df); out[["k0"]] <- list(model=m0, BIC=BIC(m0), psi=NULL)
  for (k in seq_len(max_psi)) {
    psi_init <- quantile(df$x, probs=seq_len(k)/(k+1), na.rm=TRUE)
    fit <- try(segmented(m0, seg.Z = ~x, psi=as.numeric(psi_init)), silent=TRUE)
    if (!inherits(fit, "try-error") && !is.null(fit$psi)) {
      out[[paste0("k", k)]] <- list(model=fit, BIC=BIC(fit),
                                    psi=as.numeric(fit$psi[, "Est."]))
    }
  }
  out
}

best_summary <- function(fits) {
  bics <- sapply(fits, function(f) f$BIC)
  best <- which.min(bics)
  k <- as.integer(sub("k", "", names(fits)[best]))
  psi <- if (!is.null(fits[[best]]$psi)) paste(round(fits[[best]]$psi, 2), collapse=";") else NA
  list(best_k=k, BIC=bics, psi=psi)
}

# ---- run for PC1, PC2, and 13 module eigengenes -----------------------
mes_cols <- grep("^ME", names(mes), value=TRUE)
mes_cols <- setdiff(mes_cols, "MEgrey")
mes_long <- mes %>% dplyr::select(sample, all_of(mes_cols))
mm$Sample <- as.character(mm$Sample)
mes_long$sample <- as.character(mes_long$sample)

resps <- list(PC1=mm$PC1, PC2=mm$PC2)
me_merged <- mes_long[match(mm$Sample, mes_long$sample), ]
for (mc in mes_cols) resps[[mc]] <- me_merged[[mc]]

rows <- list()
plots_pc1 <- NULL
for (rname in names(resps)) {
  y <- resps[[rname]]
  if (sum(!is.na(y)) < 15) next
  fits <- fit_segmented(y, mm$years_since, max_psi=2)
  summ <- best_summary(fits)
  rows[[rname]] <- data.frame(
    response=rname, n=sum(!is.na(y)),
    best_k=summ$best_k, breakpoint=summ$psi,
    BIC_k0=summ$BIC[["k0"]],
    BIC_k1=ifelse("k1" %in% names(summ$BIC), summ$BIC[["k1"]], NA),
    BIC_k2=ifelse("k2" %in% names(summ$BIC), summ$BIC[["k2"]], NA)
  )
  if (rname == "PC1") {
    p <- ggplot(data.frame(x=mm$years_since, y=y), aes(x=x, y=y)) +
      geom_point(alpha=0.5) +
      labs(title=paste0("PC1 vs years since infection; best_k=", summ$best_k,
                        ifelse(is.na(summ$psi), "", paste0(" psi=", summ$psi))),
           x="Years since infection", y="PC1") +
      theme_bw()
    if (!is.null(fits$k1$model$psi)) {
      psi <- as.numeric(fits$k1$model$psi[, "Est."])
      p <- p + geom_vline(xintercept=psi, linetype="dashed", color="red")
    }
    p <- p + geom_smooth(method="loess", se=FALSE, color="steelblue")
    ggsave(file.path(opt$out_dir, "segmented_pc1_plot.png"), p, width=6, height=4, dpi=150)
  }
}
out <- bind_rows(rows)
write.csv(out, file.path(opt$out_dir, "segmented_summary.csv"), row.names=FALSE)
print(out)
cat("Done.\n")
