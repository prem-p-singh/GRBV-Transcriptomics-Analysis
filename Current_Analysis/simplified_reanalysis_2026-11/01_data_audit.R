#!/usr/bin/env Rscript
# T1. Data hygiene and audit. Sample IDs are read as text so 2103 stays a string.
suppressPackageStartupMessages({library(data.table); library(readxl)})
ROOT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
OUT  <- Sys.getenv("GRBV_OUT", file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11"))
STAGES <- c("Pre_veraison", "Veraison", "Post_veraison", "Harvest")
RB_NEG_MAX_CPM <- 2    # brief says about 1; V63 (1.07) and V50 (1.12) are 1-2 reads, see logs/RULINGS.md
RB_POS_MIN_CPM <- 50   # brief: all RB+ libraries about >= 50 CPM

xl <- as.data.table(read_excel(file.path(ROOT, "Raw_Data/Correct.Sample_Metadata.xlsx"), col_types = "text"))
setnames(xl, c("Infection Group", "Sampling Season", "Years Since Infection", "Maturity stage",
               "Vine ID", "Bio-replicate ID"),
         c("infection", "season", "cohort", "stage", "vine", "biorep"))
tab <- fread(file.path(ROOT, "Raw_Data/counts_matrix.txt"))
genes <- tab[[1]]; tab[[1]] <- NULL
counts <- as.matrix(tab); rownames(counts) <- genes
stopifnot(nrow(counts) == 28168, ncol(counts) == 327, !anyDuplicated(xl$Sample))

no_library <- setdiff(xl$Sample, colnames(counts))
if (!all(colnames(counts) %in% xl$Sample))
  stop("Count columns missing from metadata: ",
       paste(setdiff(colnames(counts), xl$Sample), collapse = ", "))
if (length(no_library) != 9) stop("Expected 9 metadata rows without a library, got ", length(no_library))
if (!"2103" %in% colnames(counts)) stop("Sample 2103 not matched")

a <- xl[match(colnames(counts), Sample)]
a[, library_size := colSums(counts)]
viral <- grep("^N761_", rownames(counts), value = TRUE)
stopifnot(length(viral) == 6, sum(counts["N761_gp6", ]) == 0)  # gp6 has zero counts; 5 pass the filter
a[, viral_cpm := colSums(counts[viral, , drop = FALSE]) / library_size * 1e6]
a[, row_number := as.integer(sub("^R(\\d+)V\\d+$", "\\1", vine))]
a[, vine_position := as.integer(sub("^R\\d+V(\\d+)$", "\\1", vine))]
if (anyNA(a$row_number)) stop("Vine ID does not parse as R{row}V{pos}: ",
                              paste(a[is.na(row_number)]$vine, collapse = ", "))
keep <- c("Sample", "vine", "season", "stage", "infection", "cohort", "biorep",
          "row_number", "vine_position", "library_size", "viral_cpm")
fwrite(a[, ..keep], file.path(OUT, "data/sample_audit.csv"))
writeLines(no_library, file.path(OUT, "data/metadata_rows_without_library.txt"))

# Labels against viral reads.
bad <- a[(infection == "RB(-)" & viral_cpm > RB_NEG_MAX_CPM) |
         (infection == "RB(+)" & viral_cpm < RB_POS_MIN_CPM)]
cat("Viral CPM by label:\n")
print(a[, .(n = .N, min = min(viral_cpm), median = median(viral_cpm), max = max(viral_cpm)), by = infection])
if (nrow(bad)) { print(bad[, ..keep]); stop(nrow(bad), " libraries disagree with their infection label") }

# Acceptance counts.
v <- unique(a[, .(vine, infection)])
if (anyDuplicated(v$vine)) stop("A vine carries two infection labels")
n_vine <- v[, .N, by = infection]; print(n_vine)
stopifnot(nrow(v) == 44, n_vine[infection == "RB(+)"]$N == 36, n_vine[infection == "RB(-)"]$N == 8)
cells <- a[, .(libraries = .N, vines = uniqueN(vine)), by = .(infection, season, stage)][
  order(infection, season, match(stage, STAGES))]
print(cells)
neg <- cells[infection == "RB(-)"]
stopifnot(neg[!(season == "2022" & stage == "Harvest")]$libraries == 8,
          neg[season == "2022" & stage == "Harvest"]$libraries == 7,
          max(a[, .N, by = .(vine, season, stage)]$N) == 1)
fwrite(cells, file.path(OUT, "data/sample_cell_counts.csv"))
seas <- a[, .(seasons = uniqueN(season)), by = .(vine, infection)]
print(seas[, .N, by = .(infection, seasons)])

# Titer: audit only. Not used anywhere else until Q3 is answered.
ti <- fread(file.path(ROOT, "Raw_Data/Virus titer.txt"))
stopifnot(ncol(ti) == 7)
setnames(ti, c("Season", "Year", "Biorep", STAGES))  # file order: Pre, Veraison, Post, Harvest
long <- melt(ti, id.vars = c("Season", "Year", "Biorep"), variable.name = "stage",
             value.name = "titer", variable.factor = FALSE)
long[, `:=`(season = as.character(Season), cohort = as.character(Year), biorep = as.character(Biorep))]
m <- merge(a[infection == "RB(+)", .(Sample, season, cohort, biorep, stage, viral_cpm)], long,
           by = c("season", "cohort", "biorep", "stage"))
m <- m[is.finite(titer)]
if (anyDuplicated(m$Sample)) stop("Titer join is not one-to-one")
if (anyDuplicated(m[, .(season, cohort, biorep, stage)])) stop("Titer join is not one-to-one")
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
res <- rbind(m[, .(level = "overall", season = "all", stage = "all", n = .N, spearman_rho = sp(titer, viral_cpm))],
             m[, .(level = "cell", n = .N, spearman_rho = sp(titer, viral_cpm)), by = .(season, stage)][
               order(season, match(stage, STAGES))][, .(level, season, stage, n, spearman_rho)])
print(res)
fwrite(res, file.path(OUT, "data/titer_vs_viralcpm.csv"))
cat("T1 PASS\n")
