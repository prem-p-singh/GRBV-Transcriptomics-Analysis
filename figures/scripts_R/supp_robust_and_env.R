# Supplementary batch: pooled vs stratified robustness (S26-S29) and the
# integrated temperature-titer-DEG panel (S32) refit under DESeq2.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales); library(data.table); library(ggrepel)
})

HERE <- "/Users/prempratapsingh/Data_Analysis/.claude/worktrees/gallant-sinoussi-04711f/farm_pipeline/scripts/figures"
source(file.path(HERE, "theme_grbv.R"))

RESULTS <- "/Users/prempratapsingh/Data_Analysis/farm_pipeline/results"
ANN_F   <- file.path(RESULTS, "fix11_annotations/gene_annotations.csv")
OUT_DIR <- "/Users/prempratapsingh/Data_Analysis/figures/supp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ <- 0.05; LFC <- 1.0
stages  <- GRBV_STAGES

# --------- load ---------
all_de <- rbindlist(lapply(c(
  list.files(file.path(RESULTS, "fix01_de_redo"),
             pattern = "_deseq2\\.csv$", full.names = TRUE),
  list.files(file.path(RESULTS, "fix01c_de_extras"),
             pattern = "_deseq2\\.csv$", full.names = TRUE)
), function(f) { d <- fread(f); d$file <- basename(f); d }),
fill = TRUE) %>% as.data.frame()

ann <- fread(ANN_F) %>% as.data.frame()
name_map <- setNames(ann$product, ann$gene_id)

# ============================================================
# S26 — pooled (Q6) vs stratified (per-stage Q2) DEG concordance
# ============================================================
q6  <- all_de %>% filter(comparison == "Q6_pooled")
q6_sig <- q6$gene[!is.na(q6$padj) & q6$padj < PADJ & abs(q6$log2FC) >= LFC]

strat_sig <- lapply(stages, function(st) {
  d <- all_de %>% filter(comparison == sprintf("Q2_%s", st))
  d$gene[!is.na(d$padj) & d$padj < PADJ & abs(d$log2FC) >= LFC]
})
names(strat_sig) <- stages
strat_union <- unique(unlist(strat_sig))

overlap_df <- data.frame(
  set = c("Q6 pooled only",
          "Stratified union only",
          "Shared"),
  n = c(length(setdiff(q6_sig, strat_union)),
        length(setdiff(strat_union, q6_sig)),
        length(intersect(q6_sig, strat_union)))
) %>% mutate(set = factor(set, levels = rev(set)))

p26 <- ggplot(overlap_df,
              aes(x = n, y = set, fill = set)) +
  geom_col(alpha = 0.92, width = 0.65, color = "white") +
  geom_text(aes(label = scales::comma(n)), hjust = -0.15,
            family = "Arial", size = 3, color = "grey15") +
  scale_fill_manual(values = c("Q6 pooled only"        = unname(GRBV_COLORS["RB_pos"]),
                                "Stratified union only" = unname(GRBV_COLORS["warm"]),
                                "Shared"                = unname(GRBV_COLORS["highlight"])),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::label_comma()) +
  labs(x = "DEG count", y = NULL,
       title = "S26. Pooled Q6 vs stratified Q2 union (DESeq2)",
       subtitle = sprintf("Q6 sig.: %s;  stratified (4 stages) union: %s;  shared: %s",
                          scales::comma(length(q6_sig)),
                          scales::comma(length(strat_union)),
                          scales::comma(length(intersect(q6_sig, strat_union))))) +
  theme_grbv() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)))

save_grbv(p26, file.path(OUT_DIR, "S26_pooled_vs_stratified_comparison"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 100)

# ============================================================
# S27 — connectivity/overlap heatmap between all research questions
# ============================================================
sig_of <- function(comp) {
  d <- all_de %>% filter(comparison == comp)
  d$gene[!is.na(d$padj) & d$padj < PADJ & abs(d$log2FC) >= LFC]
}
q_sets <- list(
  sig_of("Q1a_recent_le2yr_vs_neg"),
  sig_of("Q1a_mid_3to4yr_vs_neg"),
  sig_of("Q1a_old_ge5yr_vs_neg"),
  strat_sig$Pre_veraison, strat_sig$Veraison,
  strat_sig$Post_veraison, strat_sig$Harvest,
  sig_of("Q3_seasonRBpos"),
  sig_of("Q4d_seasonXinfection_dir"),
  sig_of("Q5_seasonRBneg"),
  q6_sig)
names(q_sets) <- c("Q1a recent", "Q1a mid", "Q1a old",
                   "Q2 Pre", paste0("Q2 V", "\u00e9r"),
                   "Q2 Post", "Q2 Harvest",
                   "Q3 season", "Q4 inter.", "Q5 season", "Q6 pooled")
q_names <- names(q_sets)
ov <- expand.grid(a = q_names, b = q_names, stringsAsFactors = FALSE) %>%
  mutate(n_shared = mapply(function(x, y) length(intersect(q_sets[[x]], q_sets[[y]])), a, b),
         a_n = lengths(q_sets)[a],
         b_n = lengths(q_sets)[b],
         jac = n_shared / pmax(a_n + b_n - n_shared, 1)) %>%
  mutate(a = factor(a, levels = q_names), b = factor(b, levels = rev(q_names)))

p27 <- ggplot(ov, aes(x = a, y = b, fill = jac)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", jac)),
            family = "Arial", size = 2.2,
            color = ifelse(ov$jac > 0.35, "white", "grey15")) +
  scale_fill_gradient(low = "#fff5e1", high = unname(GRBV_COLORS["RB_pos"]),
                      name = "Jaccard") +
  labs(x = NULL, y = NULL,
       title = "S27. DEG overlap (Jaccard) across comparisons",
       subtitle = "DESeq2 DEGs, padj<0.05, |log2FC|\u22651") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 7.5),
        axis.text.y   = element_text(size = 7.5),
        legend.key.height = unit(20, "pt"))

save_grbv(p27, file.path(OUT_DIR, "S27_overlap_connectivity_heatmap"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 150)

# ============================================================
# S28 — Q1+Q2 vs Q3+Q4 cross-reference scatter
# ============================================================
q12_fc <- all_de %>% filter(comparison == "Q6_pooled") %>%
  dplyr::select(gene, log2FC_q6 = log2FC, padj_q6 = padj)
q34_fc <- all_de %>% filter(comparison == "Q4d_seasonXinfection_dir") %>%
  dplyr::select(gene, log2FC_q4 = log2FC, padj_q4 = padj)
cross <- inner_join(q12_fc, q34_fc, by = "gene") %>%
  mutate(status = case_when(
    padj_q6 < PADJ & padj_q4 < PADJ ~ "both sig.",
    padj_q6 < PADJ                  ~ "Q6 only",
    padj_q4 < PADJ                  ~ "Q4 only",
    TRUE                            ~ "neither"))
r_val <- cor(cross$log2FC_q6, cross$log2FC_q4, use = "complete.obs")

p28 <- ggplot(cross, aes(x = log2FC_q6, y = log2FC_q4, color = status)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.4, size = 0.6) +
  scale_color_manual(values = c("both sig." = unname(GRBV_COLORS["RB_pos"]),
                                 "Q6 only"  = unname(GRBV_COLORS["warm"]),
                                 "Q4 only"  = unname(GRBV_COLORS["highlight"]),
                                 "neither"  = unname(GRBV_COLORS["neutral"])),
                     name = NULL) +
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
           label = sprintf("r = %.3f", r_val),
           family = "Arial", size = 3, label.size = 0,
           fill = scales::alpha("white", 0.8)) +
  labs(x = "log2FC  Q6 (RB+ vs RB\u2212 pooled)",
       y = "log2FC  Q4 interaction direction",
       title = "S28. Q6 pooled infection effect vs Q4 season-infection interaction",
       subtitle = "Each point = one gene;  colour = which effect is significant at padj<0.05") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        legend.position = "top",
        legend.text = element_text(size = 8))

save_grbv(p28, file.path(OUT_DIR, "S28_Q1Q2_vs_Q3Q4_crossreference"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 135)

# ============================================================
# S29 — integrated summary dashboard
# ============================================================
sum_df <- data.frame(
  metric = c("Samples", "Genes tested",
             "Q6 DEGs", "Q1a DEGs (union)",
             "Q2 DEGs (union)", "Q3 DEGs", "Q4 DEGs", "Q5 DEGs"),
  value  = c(327, 14604,
             length(q6_sig),
             length(unique(c(q_sets$`Q1a recent`, q_sets$`Q1a mid`, q_sets$`Q1a old`))),
             length(unique(unlist(strat_sig))),
             length(q_sets$`Q3 season`), length(q_sets$`Q4 inter.`),
             length(q_sets$`Q5 season`))) %>%
  mutate(metric = factor(metric, levels = rev(metric)))

pS29a <- ggplot(sum_df, aes(x = value, y = metric, fill = metric)) +
  geom_col(alpha = 0.9, width = 0.72, color = "white") +
  geom_text(aes(label = scales::comma(value)),
            hjust = -0.1, family = "Arial", size = 2.8, color = "grey15") +
  scale_fill_manual(values = rep(c(unname(GRBV_COLORS["RB_pos"]),
                                    unname(GRBV_COLORS["season_2021"]),
                                    unname(GRBV_COLORS["warm"])),
                                  length.out = nrow(sum_df)),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = scales::label_comma()) +
  labs(x = "Count", y = NULL, tag = "A",
       title = "Core counts") +
  theme_grbv() +
  theme(plot.title = element_text(size = 10, face = "bold"))

# Headline findings text panel
headline <- data.frame(
  text = c(
    "Stage variance \u03b7\u00b2 = 0.426 (dominates)",
    "Inverted seasonal DEG trajectory",
    "  2021 Harvest = 1,567  |  2022 Post-v\u00e9r = 772",
    "sHSP 29/36 DOWN; HSF \u2194 sHSP r = 0.625",
    "Biomarker \u2229 QTL: 2.01x over null, p < 0.0002",
    "MBF1c (LOC100249249) = cross-pipeline mediator",
    "Host-only ML LOSO AUC = 0.76 - 0.88"
  ),
  y = seq(7, 1)
)
pS29b <- ggplot(headline, aes(x = 0, y = y, label = text)) +
  geom_text(family = "Arial", size = 3.2, hjust = 0,
            color = "grey15") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.5, 7.5)) +
  labs(tag = "B", title = "Headline findings") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 10, face = "bold",
                                   family = "Arial", hjust = 0,
                                   margin = margin(b = 6)),
        plot.tag = element_text(face = "bold", size = 12),
        plot.tag.position = c(0.02, 0.98),
        plot.margin = margin(4, 4, 4, 4))

p29 <- pS29a | pS29b
p29 <- p29 + plot_annotation(
  title = "S29. Integrated findings dashboard (DESeq2 refit)",
  theme = theme(plot.title = element_text(family = "Arial", face = "bold",
                                          size = 13, hjust = 0.5,
                                          margin = margin(b = 6))))

save_grbv(p29, file.path(OUT_DIR, "S29_comprehensive_summary_dashboard"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 125)

# ============================================================
# S32 — integrated temperature-titer-DEG per stage x season
# ============================================================
titer_raw <- read.delim("/Users/prempratapsingh/Data_Analysis/Virus titer.txt",
                        sep = "\t", check.names = FALSE)
titer_long <- titer_raw %>%
  rename(Pre_veraison  = `Pre_véraison`, Veraison = `Véraison`,
         Post_veraison = `Post_véraison`) %>%
  pivot_longer(c(Pre_veraison, Veraison, Post_veraison, Harvest),
               names_to = "Maturity_stage", values_to = "titer") %>%
  filter(!is.na(titer)) %>%
  mutate(Maturity_stage = factor(Maturity_stage, levels = stages),
         Sampling_Season = factor(Season))
titer_summary <- titer_long %>%
  group_by(Sampling_Season, Maturity_stage) %>%
  summarise(mean_titer = mean(titer),
            se = sd(titer)/sqrt(n()), .groups = "drop")

# temperature at sampling (match Main 5 computation)
sampling_dates <- list(
  `2021` = c(Pre_veraison = "2021-08-02", Veraison = "2021-08-16",
             Post_veraison = "2021-09-09", Harvest = "2021-09-21"),
  `2022` = c(Pre_veraison = "2022-07-28", Veraison = "2022-08-09",
             Post_veraison = "2022-08-24", Harvest = "2022-09-13")
)
load_cimis <- function(path) {
  w <- fread(path) %>% as.data.frame()
  colnames(w) <- gsub("\\s+", " ", colnames(w))
  w$Date <- as.Date(w$Date, format = "%m/%d/%Y")
  w$Tavg <- as.numeric(w$`Avg Air Temp (C)`)
  w[, c("Date","Tavg")]
}
w21 <- load_cimis("/Users/prempratapsingh/Data_Analysis/daily_2021.csv")
w22 <- load_cimis("/Users/prempratapsingh/Data_Analysis/daily_2022.csv")

temp_at <- function(season, date) {
  w <- if (season == "2021") w21 else w22
  d <- as.Date(date)
  mean(w$Tavg[w$Date >= (d - 7) & w$Date <= d], na.rm = TRUE)
}

s32_df <- bind_rows(lapply(names(sampling_dates), function(sn) {
  do.call(rbind, lapply(names(sampling_dates[[sn]]), function(st) {
    data.frame(Sampling_Season = sn, Maturity_stage = st,
               temp7d = temp_at(sn, sampling_dates[[sn]][[st]]))
  }))
})) %>% mutate(Maturity_stage = factor(Maturity_stage, levels = stages))

# DEG counts (Q2x)
deg_cnt <- q2x_counts <- q2x_counts_builder <- NULL  # placeholder
q2x_files <- list.files(file.path(RESULTS, "fix01c_de_extras"),
                        pattern = "^Q2x_", full.names = TRUE)
q2x_all <- rbindlist(lapply(q2x_files, function(f) {
  d <- fread(f); d$file <- basename(f); d
}), fill = TRUE) %>% as.data.frame() %>%
  mutate(file2 = sub("_deseq2\\.csv$", "", file),
         parts = strsplit(file2, "_"))
q2x_all$Sampling_Season <- sapply(q2x_all$parts, function(x) x[2])
q2x_all$Maturity_stage  <- sapply(q2x_all$parts, function(x)
  paste(x[-(1:2)], collapse = "_"))
q2x_cnt <- q2x_all %>%
  filter(!is.na(padj), padj < PADJ, abs(log2FC) >= LFC) %>%
  count(Sampling_Season, Maturity_stage, name = "n_DEG") %>%
  mutate(Maturity_stage = factor(Maturity_stage, levels = stages))

s32_all <- s32_df %>%
  left_join(titer_summary, by = c("Sampling_Season","Maturity_stage")) %>%
  left_join(q2x_cnt, by = c("Sampling_Season","Maturity_stage"))

# Normalize columns for plotting
s32_long <- s32_all %>%
  mutate(Temp = temp7d,
         Titer_M = mean_titer / 1e6,
         DEG_k   = n_DEG / 1000) %>%
  pivot_longer(c(Temp, Titer_M, DEG_k), names_to = "metric", values_to = "val") %>%
  mutate(metric = recode(metric,
                         "Temp"    = "Temperature (\u00b0C, 7-d mean)",
                         "Titer_M" = "Mean titer (\u00d7 10\u2076)",
                         "DEG_k"   = "DEGs (thousands)"),
         metric = factor(metric, levels = c("Temperature (\u00b0C, 7-d mean)",
                                            "Mean titer (\u00d7 10\u2076)",
                                            "DEGs (thousands)")))

p32 <- ggplot(s32_long,
              aes(x = Maturity_stage, y = val,
                  color = Sampling_Season, group = Sampling_Season,
                  shape = Sampling_Season)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.2, fill = "white", stroke = 1.1) +
  facet_wrap(~ metric, ncol = 3, scales = "free_y") +
  scale_color_manual(values = GRBV_SEASON_COLOR, name = "Season") +
  scale_shape_manual(values = c("2021" = 21, "2022" = 24), name = "Season") +
  scale_x_discrete(labels = GRBV_STAGE_LABEL) +
  labs(x = NULL, y = NULL,
       title = "S32. Integrated temperature - viral titer - DEG count",
       subtitle = "Three parallel trajectories across stages, by season (DESeq2 DEG counts from fix01c Q2x)") +
  theme_grbv() +
  theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey25",
                                     margin = margin(b = 6)),
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 7.5),
        strip.text    = element_text(face = "bold", size = 9),
        legend.position = "top")

save_grbv(p32, file.path(OUT_DIR, "S32_integrated_temp_titer_DEG"),
          width_mm = GRBV$W_DOUBLE_MM, height_mm = 110)

cat("[supp S26-S29, S32] done.\n")
