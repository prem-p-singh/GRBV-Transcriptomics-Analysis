# Step 4b: point estimates of the per-cell delay under alternative 1 (TMM log2 CPM), to see whether
# the cell signs, in particular Harvest 2022, depend on the normalization. No CIs (sensitivity only).
source("00_common.R")
CELL_KEYS <- as.vector(outer(STAGES, c("2021", "2022"), paste, sep = "__"))
cell_of <- factor(paste(audit$stage, audit$season, sep = "__"), CELL_KEYS)
d <- function(s) as.vector(tapply(s[!audit$healthy], cell_of[!audit$healthy], mean) - tapply(s[audit$healthy], cell_of[audit$healthy], mean))
prim <- read.csv("data/mphs_all_projection.csv", colClasses = c(Sample = "character"))
stopifnot(identical(prim$Sample, audit$Sample))
alt1 <- project(expr("tmm_logcpm"), audit$Sample, "scale")
out <- data.frame(cell = CELL_KEYS, primary_cpm = d(prim$mphs), alt1_tmm_logcpm = d(alt1))
out$same_sign <- sign(out$primary_cpm) == sign(out$alt1_tmm_logcpm)
write.csv(out, "data/mphs_delay_alt_normalization.csv", row.names = FALSE); print(out, digits = 3)
