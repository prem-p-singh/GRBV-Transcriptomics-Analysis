# Supplementary: viral transcriptome  (S56 - S59)
#
#   S56  Viral expression atlas: 5 ORFs x (stage x season x cohort)
#   S57  Virus-titer coherence: summed viral CPM vs RT-qPCR titer
#   S58  Per-ORF host-gene correlate landscape (top 30 per ORF)
#   S59  Viral-mediator scatter grid (5 viral x 5 anchor hosts)
#
# Adds the "virus point of view" to the analysis: viral-gene dynamics
# across stages, cohorts and seasons; internal validation against
# RT-qPCR titer; per-ORF host transcriptome programmes.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(ggrepel)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/figures/scripts_R"
source(file.path(HERE, "theme_grbv.R"))

DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- file.path(DATA_ROOT, "figures", "supp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------
# Load core inputs
# -----------------------------------------------------------------
counts <- fread(file.path(DATA_ROOT, "counts_matrix.txt"),
                sep = "\t", header = TRUE, data.table = FALSE)
rownames(counts) <- counts[, 1]; counts[, 1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"

meta <- read_tsv(file.path(DATA_ROOT, "sample_metadata_v2.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
rownames(meta) <- meta$Sample
meta <- meta[colnames(counts), ]   # align order

# TMM-ish normalisation: log2(CPM + 1)
lib   <- colSums(counts)
lcpm  <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)

viral_rows <- grep("^N761_", rownames(lcpm), value = TRUE)
viral_active <- viral_rows[rowSums(counts[viral_rows, ]) > 0]
cat("[viral] active genes:", paste(viral_active, collapse = ", "), "\n")

# -----------------------------------------------------------------
# Titer join: map each sample -> its per-biorep per-stage titer
# -----------------------------------------------------------------
titer_raw <- fread(file.path(DATA_ROOT, "Virus titer.txt"),
                   sep = "\t", header = TRUE, data.table = FALSE)
names(titer_raw) <- c("Season", "Year", "Biorep",
                      "Pre_veraison", "Veraison",
                      "Post_veraison", "Harvest")

stage_to_col <- c(Pre_veraison  = "Pre_veraison",
                  Veraison      = "Veraison",
                  Post_veraison = "Post_veraison",
                  Harvest       = "Harvest")

meta_titer <- meta %>%
  mutate(Sampling_Season = as.integer(Sampling_Season),
         Years_Since_Infection = as.integer(Years_Since_Infection),
         Biorep = biorep_id) %>%
  rowwise() %>%
  mutate(titer = {
    r <- titer_raw %>%
      filter(Season == Sampling_Season,
             Year   == Years_Since_Infection,
             Biorep == biorep_id)
    if (nrow(r) == 0) NA_real_ else r[[ stage_to_col[Maturity_stage] ]][1]
  }) %>% ungroup() %>% as.data.frame()

# per-sample summed viral log-CPM (sum across the 5 active viral genes)
sum_lcpm_viral <- colSums(lcpm[viral_active, , drop = FALSE])
meta_titer$viral_sum <- sum_lcpm_viral[meta_titer$Sample]

# per-sample individual viral log-CPM
viral_lcpm_df <- as.data.frame(t(lcpm[viral_active, , drop = FALSE]))
viral_lcpm_df$Sample <- rownames(viral_lcpm_df)
meta_titer <- meta_titer %>% left_join(viral_lcpm_df, by = "Sample")

# =================================================================
# S56 — Viral expression atlas
# =================================================================
stages <- c("Pre_veraison","Veraison","Post_veraison","Harvest")

viral_long <- meta_titer %>%
  dplyr::filter(Infection_Group == "RB(+)") %>%
  dplyr::select(Sample, Sampling_Season, Maturity_stage,
                Years_Since_Infection, all_of(viral_active)) %>%
  pivot_longer(cols = all_of(viral_active),
               names_to = "viral_gene", values_to = "logCPM") %>%
  mutate(Maturity_stage = factor(Maturity_stage, levels = stages),
         Sampling_Season = as.character(Sampling_Season),
         viral_gene = factor(viral_gene, levels = rev(viral_active)))

# Panel A: stage x season heatmap, rows = viral gene, cells = mean logCPM
grid_ss <- viral_long %>%
  group_by(viral_gene, Sampling_Season, Maturity_stage) %>%
  summarise(mean_lcpm = mean(logCPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(cond = paste(Sampling_Season,
                      GRBV_STAGE_LABEL[as.character(Maturity_stage)], sep = " "),
         cond = factor(cond,
                       levels = unique(cond[order(Sampling_Season, Maturity_stage)])))

S56_A <- ggplot(grid_ss, aes(x = cond, y = viral_gene, fill = mean_lcpm)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", mean_lcpm)),
            family = "Arial", size = 2.6,
            color = ifelse(grid_ss$mean_lcpm > 4, "white", "grey15")) +
  scale_fill_gradient(low = "#fef9e7",
                      high = unname(GRBV_COLORS["med_titer"]),
                      name = expression(log[2]*"(CPM+1)")) +
  labs(x = NULL, y = NULL, tag = "A",
       title = "A. Viral-gene expression per stage \u00d7 season",
       subtitle = "mean log2(CPM+1) across RB+ samples") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        axis.text.y = element_text(face = "bold", size = 9),
        legend.key.height = unit(12, "pt"))

# Panel B: viral CPM trajectory over infection-onset cohorts (Q1b)
cohort_long <- viral_long %>%
  group_by(viral_gene, Years_Since_Infection) %>%
  summarise(mean_lcpm = mean(logCPM, na.rm = TRUE),
            n = dplyr::n(), .groups = "drop") %>%
  mutate(Years_Since_Infection = factor(Years_Since_Infection))

tf_pal <- pal_tf_family[c("HSF","NAC","bHLH","WRKY","MYB","ERF/AP2")]
viral_pal <- setNames(head(unname(tf_pal), length(viral_active)), viral_active)

S56_B <- ggplot(cohort_long,
                aes(x = Years_Since_Infection, y = mean_lcpm,
                    group = viral_gene, color = viral_gene)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_color_manual(values = viral_pal, name = NULL) +
  labs(x = "Infection-onset cohort year", y = expression(log[2]*"(CPM+1)"),
       tag = "B",
       title = "B. Viral dynamics across infection duration",
       subtitle = "mean per-sample log-CPM, RB+ only, Q1b cohorts") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.key.height = unit(8, "pt"))

S56 <- (S56_A / S56_B) +
  plot_layout(heights = c(1, 0.95)) +
  plot_annotation(
    title = "S56. Viral transcriptome atlas",
    subtitle = "5 active viral ORFs (N761_gp1-gp5) across 8 stage \u00d7 season cells and 7 infection cohorts",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 12, hjust = 0.5),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6))))

save_grbv(S56, file.path(OUT_DIR, "S56_viral_expression_atlas"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 180)

# =================================================================
# S57 — Virus-titer coherence (sanity check)
# =================================================================
coh <- meta_titer %>%
  dplyr::filter(Infection_Group == "RB(+)", !is.na(titer), !is.na(viral_sum))

r_overall <- cor(coh$viral_sum, log10(coh$titer), method = "spearman",
                 use = "pairwise.complete.obs")

S57_A <- ggplot(coh, aes(x = log10(titer), y = viral_sum,
                         color = Maturity_stage)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "grey25",
              fill = "grey85", linewidth = 0.6, formula = y ~ x) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Spearman  r = %.2f\nn = %d samples",
                           r_overall, nrow(coh)),
           hjust = -0.05, vjust = 1.3, family = "Arial",
           size = 3, fontface = "bold",
           label.size = 0, fill = alpha("white", 0.8), color = "grey10") +
  scale_color_manual(values = setNames(unname(pal_stage), GRBV_STAGES),
                     labels = GRBV_STAGE_LABEL, name = "Stage") +
  labs(x = expression(log[10]*"(RT-qPCR titer, copies)"),
       y = expression(Sigma[5]*" viral " * log[2]*"(CPM+1)"),
       tag = "A",
       title = "A. Summed viral log-CPM vs RT-qPCR titer") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.key.height = unit(8, "pt"))

# Panel B: per-viral-gene Spearman
per_gene_r <- sapply(viral_active, function(g) {
  cor(coh[[g]], log10(coh$titer), method = "spearman",
      use = "pairwise.complete.obs")
})
per_gene_df <- data.frame(gene = names(per_gene_r), r = per_gene_r) %>%
  mutate(gene = factor(gene, levels = rev(viral_active)))

S57_B <- ggplot(per_gene_df, aes(x = r, y = gene, fill = r)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", r), x = r / 2),
            family = "Arial", size = 2.9, fontface = "bold",
            color = "grey10") +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "#F7F7F7",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-1, 1),
                       guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = "Spearman r with log10(titer)", y = NULL, tag = "B",
       title = "B. Per-ORF coherence with titer") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 9))

S57 <- (S57_A | S57_B) +
  plot_layout(widths = c(1.2, 0.8)) +
  plot_annotation(
    title = "S57. Virus-titer coherence check",
    subtitle = sprintf("RNA-Seq summed-viral-CPM and RT-qPCR titer measure the same thing;  overall Spearman r = %.2f (n = %d)",
                       r_overall, nrow(coh)),
    theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                            size = 12, hjust = 0.5),
                  plot.subtitle = element_text(family = "Arial", size = 9,
                                               hjust = 0.5, color = "grey25",
                                               margin = margin(b = 6))))

save_grbv(S57, file.path(OUT_DIR, "S57_viral_titer_coherence"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 120)

# =================================================================
# S58 — Per-ORF host-gene correlate landscape
# =================================================================
rbpos_samples <- meta_titer$Sample[meta_titer$Infection_Group == "RB(+)"]
host_rows <- setdiff(rownames(lcpm), viral_rows)

# Limit to expressed host genes with variance (speed + noise control)
host_var <- apply(lcpm[host_rows, rbpos_samples, drop = FALSE], 1, var)
host_top <- names(sort(host_var, decreasing = TRUE))[1:5000]
host_mat <- lcpm[host_top, rbpos_samples, drop = FALSE]

# Load NCBI annotations so labels are interpretable
ann_path <- file.path(DATA_ROOT,
                      "farm_pipeline/results/fix11_annotations/gene_annotations.csv")
ann <- tryCatch(fread(ann_path) %>% as.data.frame(), error = function(e) NULL)
name_map <- if (!is.null(ann)) setNames(ann$product, ann$gene_id) else character(0)

N_TOP_ORF <- 15   # top 15 per viral ORF is enough room for readable labels

per_orf <- lapply(viral_active, function(vg) {
  vvec <- lcpm[vg, rbpos_samples]
  rvec <- apply(host_mat, 1, function(h) cor(h, vvec, method = "spearman",
                                             use = "pairwise.complete.obs"))
  top_idx <- order(-abs(rvec))[1:N_TOP_ORF]
  data.frame(viral_gene = vg,
             host = rownames(host_mat)[top_idx],
             r = rvec[top_idx],
             row.names = NULL)
})
per_orf <- do.call(rbind, per_orf)
per_orf$product <- ifelse(per_orf$host %in% names(name_map),
                          substr(name_map[per_orf$host], 1, 24), "")
per_orf$label <- ifelse(nchar(per_orf$product) > 0,
                        paste0(per_orf$host, "  |  ", per_orf$product),
                        per_orf$host)
per_orf$viral_gene <- factor(per_orf$viral_gene, levels = viral_active)

# Within-facet ordering: create a synthetic factor that preserves |r|-rank
per_orf <- per_orf %>%
  group_by(viral_gene) %>%
  arrange(abs(r)) %>%
  mutate(ypos = row_number(),
         label_ordered = paste(viral_gene, label, sep = "|||")) %>%
  ungroup() %>% as.data.frame()
per_orf$label_ordered <- factor(per_orf$label_ordered,
                                levels = per_orf$label_ordered)

S58 <- ggplot(per_orf, aes(x = r, y = label_ordered, fill = r)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.2f", r),
                x = ifelse(r > 0, r + 0.04, r - 0.04),
                hjust = ifelse(r > 0, 0, 1)),
            family = "Arial", size = 2.2, fontface = "bold",
            color = "grey10") +
  facet_wrap(~ viral_gene, scales = "free_y", ncol = 1) +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "#F7F7F7",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman r") +
  scale_x_continuous(limits = c(-1.15, 1.15),
                     breaks = c(-1, -0.5, 0, 0.5, 1)) +
  # show only the host label portion of the synthetic factor (strip viral prefix)
  scale_y_discrete(labels = function(x) sub(".*\\|\\|\\|", "", x)) +
  labs(x = "Spearman r  (viral-gene logCPM vs host-gene logCPM)",
       y = NULL,
       title = "S58. Per-ORF host-gene correlate landscape",
       subtitle = sprintf("Top %d host genes by |Spearman r| per viral ORF; host universe = 5,000 most-variable genes in RB+ samples",
                          N_TOP_ORF)) +
  theme_grbv() +
  theme(plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        axis.text.y   = element_text(size = 7.5),
        strip.text    = element_text(face = "bold", size = 10,
                                     color = "white"),
        strip.background = element_rect(fill = unname(GRBV_COLORS["med_titer"]),
                                        color = NA),
        panel.spacing.y = unit(4, "pt"),
        legend.key.height = unit(16, "pt"))

# height: ~10 mm per row x 15 rows x 5 facets + header ~ 220 mm, or 280 for air
save_grbv(S58, file.path(OUT_DIR, "S58_per_ORF_host_correlates"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 300)

# =================================================================
# S59 — Viral-mediator scatter grid
# =================================================================
# Choose 5 anchor host genes: MBF1c, SPS1, LOC100266309, top biomarker,
# top sHSP (pick the most-recurrent one from heat_DEGs_recurrence.csv).
shsp_rec <- fread(file.path(DATA_ROOT, "heat_DEGs_recurrence.csv")) %>%
  as.data.frame() %>%
  dplyr::filter(grepl("sHSP", Category)) %>%
  arrange(desc(n_conditions), mean_FC)   # most recurrent, then strongest down
top_shsp <- shsp_rec$Gene[1]
cat("[viral] top recurrent sHSP chosen for S59:", top_shsp, "\n")

anchors <- c(
  "LOC100249249" = "MBF1c",
  "LOC100232974" = "SPS1",
  "LOC100266309" = "uncharacterised",
  "LOC100268011" = "RdRp-1",
  setNames(paste0("sHSP (", top_shsp, ")"), top_shsp)
)
# guard against anchors not present in lcpm
anchors <- anchors[names(anchors) %in% rownames(lcpm)]

grid_long <- expand.grid(viral = viral_active, host = names(anchors),
                         stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    r  = cor(lcpm[viral, rbpos_samples], lcpm[host, rbpos_samples],
             method = "spearman"),
    p  = suppressWarnings(cor.test(lcpm[viral, rbpos_samples],
                                    lcpm[host, rbpos_samples],
                                    method = "spearman")$p.value)
  ) %>% ungroup() %>% as.data.frame() %>%
  mutate(host_label = anchors[host],
         viral = factor(viral,       levels = viral_active),
         host_label = factor(host_label, levels = unname(anchors)))

# Long data for scatter
scatter_long <- do.call(rbind, lapply(seq_len(nrow(grid_long)), function(i) {
  vg <- as.character(grid_long$viral[i])
  hg <- as.character(grid_long$host[i])
  data.frame(Sample = rbpos_samples,
             viral  = grid_long$viral[i],
             host_label = grid_long$host_label[i],
             x = lcpm[vg, rbpos_samples],
             y = lcpm[hg, rbpos_samples])
})) %>%
  mutate(viral = factor(viral, levels = viral_active),
         host_label = factor(host_label, levels = unname(anchors)))

lab_df <- grid_long %>%
  mutate(label = sprintf("r=%+.2f\np=%s",
                         r, scales::scientific(p, digits = 1)))

S59 <- ggplot(scatter_long, aes(x = x, y = y)) +
  geom_point(color = unname(GRBV_COLORS["med_titer"]),
             alpha = 0.55, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE,
              color = unname(GRBV_COLORS["up"]),
              linewidth = 0.5, formula = y ~ x) +
  geom_label(data = lab_df,
             aes(x = -Inf, y = Inf, label = label),
             hjust = -0.05, vjust = 1.15,
             family = "Arial", size = 2.2, fontface = "bold",
             label.size = 0, fill = alpha("white", 0.8), color = "grey10",
             inherit.aes = FALSE) +
  facet_grid(host_label ~ viral, scales = "free", switch = "y") +
  labs(x = expression(log[2]*"(CPM+1)  —  viral gene"),
       y = expression(log[2]*"(CPM+1)  —  host anchor"),
       title = "S59. Viral-gene \u00d7 host-anchor co-expression",
       subtitle = "Per-sample scatter, RB+ only; 5 viral ORFs \u00d7 5 host anchors (mediators + top biomarker + top recurrent sHSP); per-cell Spearman r and p") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        strip.text    = element_text(face = "bold", size = 8),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        strip.placement = "outside",
        axis.text = element_text(size = 6))

save_grbv(S59, file.path(OUT_DIR, "S59_viral_host_scatter_grid"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 200)

cat("[supp_viral_transcriptome] S56-S59 complete.\n")
cat(sprintf("  coherence overall r = %.2f, n = %d\n", r_overall, nrow(coh)))
cat(sprintf("  S58: 30 correlates shown per viral ORF, %d ORFs\n",
            length(viral_active)))
cat(sprintf("  S59: %d grid cells (%d viral x %d host anchors)\n",
            nrow(grid_long), length(viral_active), length(anchors)))
