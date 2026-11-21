#!/usr/bin/env Rscript
# T12. Values for the atlas paragraph, named genes, chemistry and interaction text.
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("GRBV_ROOT", if (dir.exists("Current_Analysis")) normalizePath(".") else normalizePath("../.."))
OUT <- file.path(ROOT, "Current_Analysis/simplified_reanalysis_2026-11")
PRIOR <- file.path(ROOT, "Current_Analysis/prior_era_inputs")
un <- fread(file.path(OUT, "data/unadjusted_dream_cell_contrasts.csv"))
ad <- fread(file.path(OUT, "data/adjusted_dream_cell_contrasts.csv"))[model == "adjusted_residual"]
mem <- fread(file.path(PRIOR, "Era7_pathway_membership.csv"))
defs <- fread(file.path(PRIOR, "Era7_branch_definitions.csv"))
z <- fread(file.path(OUT, "data/zenith_adjusted.csv"))
code <- c(Pre_veraison = "Pre", Veraison = "Veraison", Post_veraison = "Post", Harvest = "Harvest")
w <- merge(un[, .(gene, stage, season, log2FC_unadj = log2FC, fdr_unadj = fdr)],
           ad[, .(gene, stage, season, log2FC_adj = log2FC, fdr_adj = fdr)], by = c("gene", "stage", "season"))
atlas <- rbindlist(lapply(defs$set_id, function(id) {
  g <- unique(mem[set_id == id]$gene)
  w[gene %in% g, .(set_id = id, name = defs[set_id == id]$name, n_genes = .N, mean_log2FC_unadjusted = mean(log2FC_unadj),
                   mean_log2FC_adjusted = mean(log2FC_adj), n_fdr05_unadjusted = sum(fdr_unadj < .05), n_fdr05_adjusted = sum(fdr_adj < .05)),
    by = .(season, stage)]
}))
atlas[, stage_code := code[stage]]
atlas <- merge(atlas, z[, .(set_id, stage_code, season, zenith_adj_p = PValue, zenith_adj_fdr_branches = FDR_branches)],
               by = c("set_id", "stage_code", "season"), all.x = TRUE)
fwrite(atlas[, !"stage_code"], file.path(OUT, "data/atlas_branch_cell_summary.csv"))
genes <- c(DCL2 = "LOC100258190", RDR1 = "LOC100268011", HSFA2e = "LOC100245950", HSP18.2 = "LOC100242821")
named <- rbindlist(lapply(names(genes), function(n) w[gene == genes[[n]]][, symbol := n]))
fwrite(named, file.path(OUT, "data/named_gene_effects.csv"))
print(named[stage %in% c("Veraison", "Harvest")][order(symbol, season, stage)])
print(atlas[stage == "Harvest" & grepl("Sucrose|Starch|Sugar transport|Stilbene|hydroxyl|Expansin|Pectin hyd|ABA rec|Lignin|Antho", name),
            .(name, season, n_genes, u = round(mean_log2FC_unadjusted, 3), a = round(mean_log2FC_adjusted, 3), zp = signif(zenith_adj_p, 2))][order(name, season)])

j <- fread(file.path(PRIOR, "Era6_gene_joint_tests.csv"))
j[, viral := grepl("^N761_", gene)]
ov <- j[, .(n_infection_stage = sum(infection_stage_fdr < .05), n_three_way = sum(infection_stage_season_fdr < .05),
            n_both = sum(infection_stage_fdr < .05 & infection_stage_season_fdr < .05),
            n_three_way_only = sum(infection_stage_season_fdr < .05 & infection_stage_fdr >= .05)), by = viral]
fwrite(ov, file.path(OUT, "data/interaction_overlap.csv")); print(ov)

ch <- fread(file.path(PRIOR, "Era7_chemistry_verified_source_values.csv"))
cat("Chemistry columns:", paste(names(ch), collapse = ", "), "\n"); print(head(ch, 2))

# Pool-level chemistry: mean of three readings of one pool, descriptive only.
ch[, group := sub("^\\d+ME\\d+PASO_(RB\\([+-]\\)\\d*)_\\d+$", "\\1", sample)]
stopifnot(uniqueN(ch[, .(season, group)]) == 11, nrow(ch) == 33)
pool <- ch[, .(Brix = mean(Brix), Brix_sd = sd(Brix), pH = mean(pH), pH_sd = sd(pH), n = .N), by = .(season, group)]
h <- pool[group == "RB(-)", .(season, healthy_Brix = Brix, healthy_pH = pH)]
pool <- merge(pool, h, by = "season")
pool[, `:=`(Brix_deficit = healthy_Brix - Brix, pH_deficit = healthy_pH - pH)]
fwrite(pool, file.path(OUT, "data/chemistry_summary.csv"))
print(pool[group != "RB(-)", .(min_deficit = min(Brix_deficit), max_deficit = max(Brix_deficit), min_pH = min(pH_deficit)), by = season])
cat("T12 values COMPLETE\n")
