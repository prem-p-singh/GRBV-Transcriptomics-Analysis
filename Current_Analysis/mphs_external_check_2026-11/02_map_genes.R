# Step 2: map the MPhS core genes (V1 VIT_ IDs) to count-matrix gene IDs, strict one-to-one.
.libPaths(c("Rlib", .libPaths()))
suppressPackageStartupMessages(library(MPhS))
load(system.file("data", "MPhSdata.rda", package = "MPhS"))
core <- names(MPhSdata$means)

cw <- read.delim("../prior_era_inputs/Era2_id_crosswalk.tsv", colClasses = "character")
counts <- read.delim("../../Raw_Data/counts_matrix.txt", check.names = FALSE,
                     colClasses = c(Gene_ID = "character"))
stopifnot(ncol(counts) == 328, nrow(counts) == 28168)
rows <- counts$Gene_ID

# The count matrix names most genes LOC..., but some by gene symbol; resolve each crosswalk row to a count row.
cw <- cw[cw$VIT_id != "" & cw$LOC_id != "", ]
cw$row_id <- ifelse(cw$LOC_id %in% rows, cw$LOC_id, ifelse(cw$Symbol != "" & cw$Symbol %in% rows, cw$Symbol, NA))
cw <- unique(cw[, c("VIT_id", "LOC_id", "row_id")])
n_loc_per_vit <- table(cw$VIT_id); n_vit_per_loc <- table(cw$LOC_id)

map <- data.frame(VIT_id = core)
map$n_LOC <- as.integer(n_loc_per_vit[core]); map$n_LOC[is.na(map$n_LOC)] <- 0L
hit <- cw[match(core, cw$VIT_id), ]
map$LOC_id <- ifelse(map$n_LOC == 1, hit$LOC_id, NA)
map$n_VIT_for_LOC <- ifelse(is.na(map$LOC_id), NA, as.integer(n_vit_per_loc[map$LOC_id]))
map$row_id <- ifelse(map$n_LOC == 1, hit$row_id, NA)
map$status <- with(map, ifelse(n_LOC == 0, "unmapped",
                       ifelse(n_LOC > 1, "one_to_many",
                       ifelse(n_VIT_for_LOC > 1, "many_to_one",
                       ifelse(is.na(row_id), "one_to_one_not_in_counts", "one_to_one_in_counts")))))

m <- as.matrix(counts[, -1]); rownames(m) <- rows
pass <- rowSums(m >= 10) >= 32
map$passes_filter <- map$status == "one_to_one_in_counts" & pass[map$row_id] %in% TRUE

mk <- read.csv("../prior_era_inputs/Era7_healthy_clock_markers.csv", colClasses = "character")
markers <- unique(mk$gene)
map$inhouse_marker <- map$row_id %in% markers
write.csv(map, "data/mphs_gene_map.csv", row.names = FALSE)

s <- data.frame(
  item = c("core_genes_total", "unmapped", "one_to_many_dropped", "many_to_one_dropped",
           "one_to_one", "one_to_one_present_in_counts", "passing_filter_10counts_in_32libs",
           "inhouse_markers_union_across_folds", "usable_core_genes_that_are_inhouse_markers",
           "usable_core_genes_not_inhouse_markers"),
  n = c(length(core), sum(map$status == "unmapped"), sum(map$status == "one_to_many"),
        sum(map$status == "many_to_one"), sum(grepl("^one_to_one", map$status)),
        sum(map$status == "one_to_one_in_counts"), sum(map$passes_filter), length(markers),
        sum(map$passes_filter & map$inhouse_marker), sum(map$passes_filter & !map$inhouse_marker)))
write.csv(s, "data/mphs_overlap_summary.csv", row.names = FALSE)
print(s)
cat("gate: usable >= 120 ?", s$n[s$item == "passing_filter_10counts_in_32libs"] >= 120, "\n")
