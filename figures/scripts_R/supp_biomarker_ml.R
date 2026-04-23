# Supplementary: biomarker detail + ML reframe — S20..S22, S46..S48.
#
#   S20  biomarker score distribution (all biomarkers) + top-50 labels
#   S21  top-16 biomarker expression profiles across stage x season
#   S22  biomarker robustness: FC consistency across conditions
#   S46  multi-method AUC comparison (random + LOSO)
#   S47  LASSO stability selection frequency distribution
#   S48  consensus biomarker panel with ROC summary
#
# Inputs:
#   GRBV_biomarkers_top50.csv, counts_matrix.txt, sample_metadata.txt
#   advanced_ml_comparison_FIX04.csv, stability_selection_FIX04.csv
#   consensus_biomarker_panel.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(ggrepel)
})
HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))
DATA_ROOT <- "/Users/prempratapsingh/Data_Analysis"
OUT_DIR   <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

bm  <- fread(file.path(DATA_ROOT, "GRBV_biomarkers_top50.csv")) %>% as.data.frame()

# -----------------------------------------------------------------
# S20 — biomarker score distribution + labelled top-50
# -----------------------------------------------------------------
bm20 <- bm %>% arrange(desc(Biomarker_Score))
bm20$rank <- seq_len(nrow(bm20))
S20 <- ggplot(bm20, aes(x = rank, y = Biomarker_Score, color = Direction)) +
  geom_point(size = 1.2, alpha = 0.85) +
  geom_text_repel(data = head(bm20, 8),
                  aes(label = Gene_ID), size = 2.4, family = "Arial",
                  box.padding = 0.4, segment.size = 0.25,
                  max.overlaps = Inf) +
  scale_color_manual(values = c(UP = unname(GRBV_COLORS["up"]),
                                DOWN = unname(GRBV_COLORS["down"])),
                     name = NULL) +
  labs(x = "Biomarker rank", y = "Composite score",
       title = "Supplementary S20. Biomarker composite-score distribution (top 50)") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1))
save_grbv(S20, file.path(OUT_DIR, "S20_biomarker_overview"),
          width_mm = GRBV$W_ONEHALF_MM, height_mm = 100)

# -----------------------------------------------------------------
# S21 — top-16 biomarker expression profiles across stage x season (boxplots)
# -----------------------------------------------------------------
counts <- fread(file.path(DATA_ROOT, "counts_matrix.txt"), sep = "\t",
                header = TRUE, data.table = FALSE)
rownames(counts) <- counts[,1]; counts[,1] <- NULL
counts <- as.matrix(counts); storage.mode(counts) <- "integer"
keep <- rowSums(counts >= 10) >= 0.10 * ncol(counts); counts <- counts[keep, ]
lib <- colSums(counts); lcpm <- log2(sweep(counts, 2, lib, "/") * 1e6 + 1)

meta <- read_tsv(file.path(DATA_ROOT, "sample_metadata.txt"),
                 show_col_types = FALSE) %>% as.data.frame()
rownames(meta) <- meta$Sample; meta <- meta[colnames(lcpm), ]
meta$Maturity_stage <- factor(meta$Maturity_stage, levels = GRBV_STAGES)

top16 <- bm %>% arrange(desc(Biomarker_Score)) %>% head(16)
top16_in <- intersect(top16$Gene_ID, rownames(lcpm))
expr_long <- as.data.frame(lcpm[top16_in, ]) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "Sample", values_to = "lcpm") %>%
  left_join(meta, by = "Sample") %>%
  left_join(bm %>% select(Gene_ID, Direction), by = c("gene" = "Gene_ID"))
expr_long$gene <- factor(expr_long$gene, levels = top16_in)

S21 <- ggplot(expr_long,
              aes(x = Maturity_stage, y = lcpm, fill = Infection_Group)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3,
               linewidth = 0.25, width = 0.75) +
  facet_wrap(~ gene, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = GRBV_INFECTION_COLOR, name = NULL) +
  scale_x_discrete(labels = substr(GRBV_STAGE_LABEL, 1, 4)) +
  labs(x = NULL, y = "log2(CPM + 1)",
       title = "Supplementary S21. Top-16 biomarker expression profiles by stage x infection") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
        strip.text  = element_text(size = 7.5, face = "bold"),
        legend.position = "top")
save_grbv(S21, file.path(OUT_DIR, "S21_biomarker_profiles"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 200)

# -----------------------------------------------------------------
# S22 — biomarker FC consistency across conditions
# -----------------------------------------------------------------
# Load DESeq2 results to compute per-biomarker FC across the 8 stage x season
q2x_dir <- file.path(DATA_ROOT, "farm_pipeline/results/fix01c_de_extras")
q2x_files <- list.files(q2x_dir, pattern = "Q2x_.*_deseq2\\.csv", full.names = TRUE)
q2x <- rbindlist(lapply(q2x_files, function(f) {
  d <- fread(f); parts <- strsplit(basename(f), "_")[[1]]
  d$season <- parts[2]; d$stage <- sub("_deseq2\\.csv$", "",
                                        paste(parts[3:(length(parts)-1)], collapse = "_"))
  d
}), fill = TRUE) %>% as.data.frame()
q2x$condition <- paste(q2x$season, q2x$stage, sep = "_")

top24 <- bm %>% arrange(desc(Biomarker_Score)) %>% head(24)
bm_fc <- q2x[q2x$gene %in% top24$Gene_ID, ]
bm_fc$gene <- factor(bm_fc$gene, levels = top24$Gene_ID)
bm_fc$stage <- factor(bm_fc$stage, levels = GRBV_STAGES)
bm_fc$season <- factor(bm_fc$season, levels = c("2021","2022"))

S22 <- ggplot(bm_fc,
              aes(x = stage, y = gene, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ season, ncol = 2) +
  scale_fill_gradient2(low = unname(GRBV_COLORS["down"]),
                       mid = "white",
                       high = unname(GRBV_COLORS["up"]),
                       midpoint = 0, limits = c(-5, 5),
                       oob = scales::squish, name = "log2FC") +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL,
       title = "Supplementary S22. Top-24 biomarker log2FC consistency across stage x season") +
  theme_grbv() +
  theme(plot.title  = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 6.5),
        strip.text  = element_text(size = 9, face = "bold"),
        legend.key.height = unit(20, "pt"))
save_grbv(S22, file.path(OUT_DIR, "S22_biomarker_robustness"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 170)

# -----------------------------------------------------------------
# S46 — ML method comparison (random + LOSO)
# -----------------------------------------------------------------
auc <- fread(file.path(DATA_ROOT, "advanced_ml_comparison_FIX04.csv")) %>% as.data.frame()
auc_long <- auc %>%
  select(panel, method, random_5fold_mean_auc, loco_season_mean_auc) %>%
  pivot_longer(c(random_5fold_mean_auc, loco_season_mean_auc),
               names_to = "cv", values_to = "auc") %>%
  mutate(cv = recode(cv,
                     "random_5fold_mean_auc" = "Random 5-fold",
                     "loco_season_mean_auc"  = "Leave-one-season-out"),
         panel = recode(panel,
                        "all_vars" = "All vars (host + viral)",
                        "host_only" = "Host only"),
         method = factor(method, levels = c("LR","RF","LASSO","SVM-RBF","GB")))

S46 <- ggplot(auc_long, aes(x = method, y = auc, fill = panel)) +
  geom_col(position = position_dodge(0.78), width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("%.2f", auc)),
            position = position_dodge(0.78),
            family = "Arial", size = 2.4, vjust = -0.6) +
  facet_wrap(~ cv, ncol = 2) +
  scale_fill_manual(values = c("All vars (host + viral)" = unname(GRBV_COLORS["RB_pos"]),
                               "Host only" = unname(GRBV_COLORS["RB_neg"])),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  geom_hline(yintercept = 0.5, color = "grey60",
             linetype = "dashed", linewidth = 0.3) +
  labs(x = NULL, y = "AUC",
       title = "Supplementary S46. Multi-method AUC: random 5-fold vs LOSO season") +
  theme_grbv() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.position = "top",
        strip.text  = element_text(size = 9, face = "bold"))
save_grbv(S46, file.path(OUT_DIR, "S46_ml_method_comparison"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 110)

# -----------------------------------------------------------------
# S47 — LASSO stability selection frequency distribution
# -----------------------------------------------------------------
stab <- fread(file.path(DATA_ROOT, "stability_selection_FIX04.csv")) %>% as.data.frame()
stab <- stab %>% filter(panel == "host_only")

pA47 <- ggplot(stab, aes(x = lasso_selection_freq)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = unname(GRBV_COLORS["RB_neg"]),
                 color = "white", alpha = 0.9) +
  geom_vline(xintercept = 0.6, color = "grey40",
             linetype = "dashed", linewidth = 0.3) +
  annotate("text", x = 0.62, y = Inf, label = "stability threshold",
           hjust = 0, vjust = 1.4, family = "Arial", size = 2.6, color = "grey25") +
  labs(x = "LASSO selection frequency (100 subsamples)",
       y = sprintf("Number of host genes (n = %s)", scales::comma(nrow(stab))),
       title = "A. Selection-frequency distribution") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"))

top_stab <- stab %>% arrange(desc(lasso_selection_freq)) %>% head(15) %>%
  mutate(gene = factor(gene, levels = rev(gene)))
pB47 <- ggplot(top_stab, aes(x = lasso_selection_freq, y = gene)) +
  geom_col(fill = unname(GRBV_COLORS["RB_neg"]), alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", lasso_selection_freq)),
            hjust = -0.15, family = "Arial", size = 2.6) +
  scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  labs(x = "Selection frequency", y = NULL,
       title = "B. Top-15 stably-selected host genes") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7.5))

S47 <- pA47 | pB47 +
  plot_annotation(title = "Supplementary S47. LASSO stability selection on host-only features",
                  theme = theme(plot.title = element_text(family = "Arial",
                                                          face = "bold", size = 12,
                                                          hjust = 0.5)))
save_grbv(S47, file.path(OUT_DIR, "S47_stability_selection"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 110)

# -----------------------------------------------------------------
# S48 — consensus biomarker panel (from consensus_biomarker_panel.csv)
# -----------------------------------------------------------------
cons_f <- file.path(DATA_ROOT, "consensus_biomarker_panel.csv")
if (file.exists(cons_f)) {
  cons <- fread(cons_f) %>% as.data.frame()
  # Choose the most informative numeric columns dynamically
  num_cols <- sapply(cons, is.numeric)
  first_num <- names(cons)[num_cols][1]
  gene_col  <- grep("gene|id|symbol", names(cons), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(gene_col) && !is.na(first_num)) {
    cons2 <- cons %>% arrange(desc(.data[[first_num]])) %>% head(15)
    cons2[[gene_col]] <- factor(cons2[[gene_col]], levels = rev(cons2[[gene_col]]))
    S48 <- ggplot(cons2, aes(x = .data[[first_num]], y = .data[[gene_col]])) +
      geom_col(fill = unname(GRBV_COLORS["highlight"]), alpha = 0.9, width = 0.7) +
      geom_text(aes(label = sprintf("%.3f", .data[[first_num]])),
                hjust = -0.15, family = "Arial", size = 2.6) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = first_num, y = NULL,
           title = sprintf("Supplementary S48. Consensus biomarker panel (top 15 by %s)",
                           first_num)) +
      theme_grbv() +
      theme(plot.title = element_text(size = 11, face = "bold"),
            axis.text.y = element_text(size = 7.5))
    save_grbv(S48, file.path(OUT_DIR, "S48_consensus_panel"),
              width_mm = GRBV$W_ONEHALF_MM, height_mm = 120)
  }
}

cat("[supp_biomarker_ml] wrote S20, S21, S22, S46, S47, S48 to\n  ", OUT_DIR, "\n")
