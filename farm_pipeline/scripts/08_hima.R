#!/usr/bin/env Rscript
# 08. High-dimensional mediation (HIMA).
# Exposure  : per-season temperature (mean tavg).
# Mediator  : gene expression (top-variable, log-CPM).
# Outcome   : log10 viral titer per sample.
# RB+ samples only.

suppressPackageStartupMessages({
  library(data.table); library(optparse); library(HIMA)
})

opts <- parse_args(OptionParser(option_list=list(
  make_option("--counts"), make_option("--meta"), make_option("--titer"),
  make_option("--weather_2021"), make_option("--weather_2022"),
  make_option("--out"),
  make_option("--top_var", type="integer", default=3000),
  make_option("--cores",   type="integer", default=8)
)))

counts <- fread(opts$counts, data.table=FALSE); rownames(counts) <- counts[,1]; counts[,1] <- NULL
meta   <- fread(opts$meta,   data.table=FALSE)
titer  <- fread(opts$titer)

# Season temperature exposure (mean Tavg per season)
agg_tavg <- function(path) {
  w <- fread(path)
  cn <- tolower(colnames(w))
  tmax <- colnames(w)[grep("tmax|max", cn)][1]
  tmin <- colnames(w)[grep("tmin|min", cn)][1]
  tavg <- colnames(w)[grep("tavg|mean", cn)][1]
  if (is.na(tavg)) w[, tavg := (get(tmax) + get(tmin))/2] else setnames(w, tavg, "tavg")
  mean(w$tavg, na.rm=TRUE)
}
season_temp <- c(`2021`=agg_tavg(opts$weather_2021),
                 `2022`=agg_tavg(opts$weather_2022))

# Expression matrix (RB+ only)
expr <- as.matrix(counts)
keep_g <- rowSums(expr >= 10) >= 0.10 * ncol(expr)
expr <- expr[keep_g, ]
libsize <- colSums(expr)
lcpm <- log2(sweep(expr, 2, libsize, "/") * 1e6 + 1)

rb_pos <- meta$Sample[meta$Infection_Group == "RB(+)" & meta$Sample %in% colnames(lcpm)]
mm <- meta[match(rb_pos, meta$Sample), ]

# Per-sample titer (mean of biorep/stage match)
rename_stage <- c(`Pre_véraison`="Pre_veraison", `Véraison`="Veraison",
                  `Post_véraison`="Post_veraison", `Harvest`="Harvest")
titer_long <- melt(titer, id.vars=c("Season","Year","Biorep"),
                   measure.vars=c("Pre_véraison","Véraison","Post_véraison","Harvest"),
                   variable.name="stage_raw", value.name="titer")
titer_long[, Maturity_stage := rename_stage[as.character(stage_raw)]]
titer_avg <- titer_long[, .(titer=mean(titer, na.rm=TRUE)),
                        by=.(Season, Year, Maturity_stage)]
key <- data.table(Sample=mm$Sample, Season=as.integer(mm$Sampling_Season),
                  Year=as.integer(mm$Years_Since_Infection),
                  Maturity_stage=mm$Maturity_stage)
ts <- merge(key, titer_avg, by=c("Season","Year","Maturity_stage"), all.x=TRUE, sort=FALSE)
Y  <- log10(ts$titer[match(mm$Sample, ts$Sample)] + 1)

X  <- unname(season_temp[as.character(mm$Sampling_Season)])
v  <- apply(lcpm[, rb_pos], 1, var)
top <- names(sort(v, decreasing=TRUE))[1:opts$top_var]
M   <- t(lcpm[top, rb_pos])      # samples × genes

ok <- complete.cases(X, Y) & rowSums(is.na(M)) == 0
M  <- M[ok, ]; X <- X[ok]; Y <- Y[ok]
message(sprintf("HIMA: n=%d samples, m=%d mediators", length(Y), ncol(M)))

res <- hima(X=X, Y=Y, M=M, Y.family="gaussian", M.family="gaussian", verbose=TRUE)
res <- as.data.frame(res)
res$gene <- rownames(res)
res <- res[order(res$`BH.FDR`), ]

dir.create(dirname(opts$out), recursive=TRUE, showWarnings=FALSE)
fwrite(res, opts$out)
message("Done.")
