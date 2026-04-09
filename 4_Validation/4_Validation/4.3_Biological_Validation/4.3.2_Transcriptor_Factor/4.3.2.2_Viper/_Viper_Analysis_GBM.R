# ======================================================================
# Script: VIPER_Consensus_TF_Families_And_Correlations.R
# Author: Gustavo Bueno Moreira
# Last update: 2025-09-18 (America/Sao_Paulo)
#
# Purpose
#   Build TF tables for consensus families (5/5), visualize family composition,
#   load per-study data, run VST, prepare DoRothEA regulons, run VIPER,
#   correlate NES with LT_enrichment, and summarize consistency/meta-analysis.
#
# Notes
#   - Organization/style aligned to your previous script (numbered sections).
#   - No logic changes, except restoring the missing stacked bar chart.
# ======================================================================

#### 0) Libraries & working directory #########################################

# 0.1) Core libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(biomaRt)
  library(dorothea)
  library(viper)
  library(readr)
  library(tidyr)
  library(coin)
  library(forcats)
  library(tidyr)
  library(forcats)
})

# 0.2) Working directory
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC")


#### 1) Consensus families (5/5): inputs and outputs ###########################

# 1.1) Load inputs
metadata_Jaspar <- read.csv("Transcription _Factor/Metadata_Jaspar/tf_metadata.csv", row.names = 1)
motifs          <- read.csv("Transcription _Factor/Consensus_Motif/TF_motifs_with_family.csv")

# 1.2) Select consensus motifs (5/5) and drop Study
motifs_5of5 <- motifs %>%
  dplyr::filter(n_studies_present == 5) %>%   # or: recurrence == "5/5"
  dplyr::select(jaspar_id, jaspar_base, family) %>%
  dplyr::distinct(jaspar_base, family, .keep_all = TRUE)

# 1.3) Define families of interest from the 5/5 set
families_interest <- sort(unique(motifs_5of5$family))

# 1.4) Output A — TFs from families that reached 5/5 (family-based)
tf_in_families <- metadata_Jaspar %>%
  dplyr::filter(!is.na(family) & family %in% families_interest) %>%
  dplyr::arrange(family, TF) %>%
  dplyr::select(
    family, TF, Ensembl, gene_biotype, chromosome_name, transcript_count,
    source, complex_name, matrix_id, base_id
  )

# 1.5) Includ the NR1;NR2

# Define the TFs you want to add
extra_tfs <- c("PPARG", "NR1H2", "RARA", "NR1H4", "PPARA")

# Filter them from metadata_Jaspar and assign the custom family name
extra_rows <- metadata_Jaspar %>%
  dplyr::filter(TF %in% extra_tfs) %>%
  dplyr::mutate(
    family = "Thyroid hormone receptor-related factors (NR1);RXR-related receptors (NR2)"
  ) %>%
  dplyr::select(
    family, TF, Ensembl, gene_biotype, chromosome_name, transcript_count,
    source, complex_name, matrix_id, base_id
  )

# Append (bind) these rows to the existing table
tf_in_families <- dplyr::bind_rows(tf_in_families, extra_rows)

#### 2) Quick pie chart for families of interest ###############################

# 2.1) Family-level summary (count and proportion)
family_summary_interest <- tf_in_families %>%
  dplyr::filter(!is.na(family) & family != "") %>%
  dplyr::count(family, name = "n") %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::mutate(pct = n / sum(n))

# 2.2) Pie chart
ggplot2::ggplot(
  family_summary_interest,
  ggplot2::aes(x = "", y = pct, fill = forcats::fct_reorder(family, n, .desc = TRUE))
) +
  ggplot2::geom_col(width = 1, color = NA) +
  ggplot2::coord_polar(theta = "y") +
  ggplot2::labs(
    title = "TF families",
    fill  = "Family", x = NULL, y = NULL
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "right",
    legend.title = ggplot2::element_text(face = "bold"),
    legend.text  = ggplot2::element_text(size = 9),
    plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold")
  )

# 2.3) Save
ggplot2::ggsave("Transcription _Factor/Viper_Analysis/tf_families_interest_pie.png", width = 8, height = 5, dpi = 300)


#### 3) Load study data and compute VST per study ##############################

# 3.1) Small loader for a single study (counts + metadata)
load_study_data <- function(expression_path, metadata_path) {
  expData  <- read.csv(expression_path, row.names = 1)
  metadata <- read.csv(metadata_path,   row.names = 1)
  list(expData = expData, metadata = metadata)
}

# 3.2) Paths per study
study_paths <- list(
  Anjos_2020_T1   = list(expr = "Type_1_Studies/Anjos-Afonso_etal_2021/DATA/counts_Final.csv",
                 meta = "Type_1_Studies/Anjos-Afonso_etal_2021/DATA/metadata.csv"),
  Xie_2019_T1  = list(expr = "Type_1_Studies/Xie_etal_2019/DATA/counts_NonNormalized.csv",
                 meta = "Type_1_Studies/Xie_etal_2019/DATA/metadata.csv"),
  Fares_2016_T2   = list(expr = "Type_2_Studies/Fares_2017/DATA/counts_NonNormalized.csv",
                 meta = "Type_2_Studies/Fares_2017/DATA/metadata.csv"),
  Fares_2014_T2 = list(expr = "Type_2_Studies/Sauvageau_2014/DATA/counts_Final.csv",
                   meta = "Type_2_Studies/Sauvageau_2014/DATA/metadata.csv"),
  Xie_2019_T2  = list(expr = "Type_2_Studies/Xie_2019/DATA/counts_Final.csv",
                 meta = "Type_2_Studies/Xie_2019/DATA/metadata.csv")
)

# 3.3) Load all studies into a named list
studies <- lapply(study_paths, function(paths) load_study_data(paths$expr, paths$meta))

# 3.4) VST normalization helper (returns samples × genes)
normalize_with_vst <- function(count_data, metadata) {
  metadata <- metadata[colnames(count_data), , drop = FALSE]  # align sample order
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ 1)
  vsd <- varianceStabilizingTransformation(dds)
  vst_matrix <- assay(vsd)
  as.data.frame(t(vst_matrix))  # samples in rows, genes in columns
}

# 3.5) Compute VST per study
for (study_name in names(studies)) {
  message("Normalizing study: ", study_name)
  count_data <- studies[[study_name]]$expData
  metadata   <- studies[[study_name]]$metadata
  vst_result <- normalize_with_vst(count_data, metadata)
  studies[[study_name]]$vstData <- vst_result
}


#### 4) TFs of interest present per study ######################################

# 4.1) Master TF table (no NA Ensembl)
tf_master <- tf_in_families %>%
  dplyr::select(TF, Ensembl, family) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(Ensembl) & Ensembl != "")

# 4.2) For each study, keep TFs whose Ensembl is present in the count matrix
expr_long <- lapply(names(studies), function(st) {
  ens_in_matrix <- rownames(studies[[st]]$expData)
  present <- tf_master %>%
    dplyr::semi_join(
      dplyr::tibble(Ensembl = ens_in_matrix),
      by = "Ensembl"
    ) %>%
    dplyr::mutate(study = st)
  present
}) %>% dplyr::bind_rows() %>%
  dplyr::select(study, TF, Ensembl, family) %>%
  dplyr::arrange(study, family, TF)

# 4.3) Counts by study (total TFs)
expr_counts <- expr_long %>%
  dplyr::count(study, name = "n_TFs") %>%
  dplyr::arrange(dplyr::desc(n_TFs))


#### 5) Stacked bar chart: TF family composition per study #####################

# 5.1) Counts by study × family (this object was missing before)
expr_counts_by_family <- expr_long %>%
  dplyr::count(study, family, name = "n_TFs") %>%
  dplyr::arrange(study, dplyr::desc(n_TFs))

# 5.2) Totals per study (for top labels)
expr_totals <- expr_counts_by_family %>%
  dplyr::group_by(study) %>%
  dplyr::summarise(n_TFs = sum(n_TFs), .groups = "drop")

# 5.3) Plot (stacked bars, no bar borders, total labels above bars)
p2 <- ggplot(expr_counts_by_family, aes(x = study, y = n_TFs, fill = family)) +
  geom_col() +  # stacked bars without borders
  geom_text(
    data = expr_totals,
    aes(x = study, y = n_TFs, label = n_TFs),
    vjust = -0.3, size = 3.5, inherit.aes = FALSE
  ) +
  labs(title = "Family composition (TFs of interest per study)",
       x = "Study", y = "# TFs", fill = "Family") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(p2)
ggsave("Transcription _Factor/Viper_Analysis/tf_family_composition_stacked.png", width = 8, height = 5, dpi = 300)


#### 6) Ensembl -> HGNC mapping (via VST), collapsing by highest mean ##########

# 6.1) Ensembl mart
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 6.2) Map Ensembl -> HGNC once
ens_ids_all <- unique(unlist(lapply(studies, function(s) colnames(s$vstData))))
ens_ids_all <- ens_ids_all[!is.na(ens_ids_all) & ens_ids_all != ""]

ens_map_raw <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ens_ids_all,
  mart       = ensembl
)

ens_map <- ens_map_raw
colnames(ens_map)[1:2] <- c("Ensembl", "Symbol")
ens_map <- ens_map |>
  dplyr::filter(!is.na(Symbol) & Symbol != "") |>
  dplyr::distinct()

# 6.3) Master TFs for QC
tf_master <- tf_in_families %>%
  dplyr::select(TF, Ensembl, family) %>%
  distinct() %>%
  dplyr::filter(!is.na(Ensembl) & Ensembl != "")

# 6.4) Helper: VST (samples x Ensembl) -> expr_sym (HGNC x samples)
vst_to_symbols_max <- function(vst_df, ens_map_tbl) {
  common <- intersect(colnames(vst_df), ens_map_tbl$Ensembl)
  if (length(common) == 0) return(NULL)
  M <- as.matrix(vst_df[, common, drop = FALSE])
  sym_vec <- ens_map_tbl$Symbol[match(colnames(M), ens_map_tbl$Ensembl)]
  col_means <- colMeans(M, na.rm = TRUE)
  idx_keep <- tapply(seq_along(sym_vec), sym_vec, function(ix) ix[which.max(col_means[ix])])
  M_sel <- M[, unlist(idx_keep), drop = FALSE]
  colnames(M_sel) <- names(idx_keep)
  t(M_sel)
}

# 6.5) Apply per study + QC
qc_map_list <- list()

for (st in names(studies)) {
  vst <- studies[[st]]$vstData
  n_ens_in <- ncol(vst)
  expr_sym <- vst_to_symbols_max(vst, ens_map)
  studies[[st]]$expr_sym <- expr_sym
  
  n_sym_out <- if (is.null(expr_sym)) 0 else nrow(expr_sym)
  pct_map   <- if (n_ens_in > 0) 100 * n_sym_out / n_ens_in else 0
  
  ens_in_matrix <- colnames(vst)
  tf_in_matrix  <- tf_master %>% semi_join(tibble(Ensembl = ens_in_matrix), by = "Ensembl")
  n_tf_in       <- nrow(tf_in_matrix)
  
  tf_with_symbol <- if (is.null(expr_sym)) {
    tf_master[0, , drop = FALSE]
  } else {
    tibble(Symbol = rownames(expr_sym)) %>%
      inner_join(tf_master, by = c("Symbol" = "TF"))
  }
  n_tf_with_symbol <- nrow(tf_with_symbol)
  
  qc_map_list[[st]] <- tibble(
    study = st,
    genes_in_vst         = n_ens_in,
    genes_with_symbol    = n_sym_out,
    pct_genes_mapped     = round(pct_map, 1),
    tf_interest_in_vst   = n_tf_in,
    tf_interest_with_sym = n_tf_with_symbol
  )
}

qc_mapping <- bind_rows(qc_map_list) %>%
  arrange(desc(pct_genes_mapped))

print(qc_mapping)


#### 7) DoRothEA A–C + filter by expressed targets (>= 10) ####################

# 7.1) Load DoRothEA human A–C
data("dorothea_hs", package = "dorothea")

dorothea_df <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::select(tf, target, mor, confidence) %>%
  dplyr::distinct()

tfs_in_regulon <- unique(dorothea_df$tf)

# 7.2) Coverage QC
tf_interest_symbols <- tf_in_families %>%
  dplyr::filter(!is.na(TF) & TF != "") %>%
  dplyr::pull(TF) %>%
  unique()

qc_dorothea_overall <- tibble(
  n_tf_interest_total = length(tf_interest_symbols),
  n_tf_in_dorothea    = length(intersect(tf_interest_symbols, tfs_in_regulon)),
  n_tf_not_in_dorothea= length(setdiff(tf_interest_symbols, tfs_in_regulon))
)

print(qc_dorothea_overall)
# tf_missing_overall <- setdiff(tf_interest_symbols, tfs_in_regulon)

# 7.3) Per-study regulon (expressed targets) with minimum targets
min_targets <- 10  # your relaxed criterion

qc_dorothea_by_study <- list()

for (st in names(studies)) {
  expr_sym <- studies[[st]]$expr_sym  # HGNC x samples
  if (is.null(expr_sym) || nrow(expr_sym) == 0) {
    warning("expr_sym empty for study: ", st); next
  }
  genes_expressed <- rownames(expr_sym)
  
  dorothea_study_df <- dorothea_df %>%
    dplyr::filter(target %in% genes_expressed)
  
  target_counts <- dorothea_study_df %>%
    dplyr::count(tf, name = "n_targets") %>%
    dplyr::arrange(dplyr::desc(n_targets))
  
  keep_tfs <- target_counts %>%
    dplyr::filter(n_targets >= min_targets) %>%
    dplyr::pull(tf)
  
  dorothea_study_kept <- dorothea_study_df %>%
    dplyr::filter(tf %in% keep_tfs)
  
  studies[[st]]$regulon_filtered <- dorothea::df2regulon(dorothea_study_kept)
  studies[[st]]$regulon_counts   <- target_counts
  
  interest_in_study  <- intersect(tf_interest_symbols, target_counts$tf)
  interest_passing   <- intersect(tf_interest_symbols, keep_tfs)
  
  qc_dorothea_by_study[[st]] <- tibble(
    study                 = st,
    genes_expressed       = length(genes_expressed),
    n_edges_after_filter  = nrow(dorothea_study_kept),
    n_tf_in_study_regulon = length(unique(dorothea_study_df$tf)),
    n_tf_passing_min_k    = length(keep_tfs),
    n_interest_present    = length(interest_in_study),
    n_interest_passing    = length(interest_passing)
  )
}

qc_dorothea_by_study <- dplyr::bind_rows(qc_dorothea_by_study) %>%
  dplyr::arrange(dplyr::desc(n_interest_passing))

print(qc_dorothea_by_study)


#### 8) VIPER per study: NES (TF x samples) ####################################

for (st in names(studies)) {
  message("Running VIPER on study: ", st)
  
  expr_sym <- studies[[st]]$expr_sym            # HGNC x samples
  regulon  <- studies[[st]]$regulon_filtered    # DoRothEA A–C filtered
  
  if (is.null(expr_sym) || is.null(regulon) || nrow(expr_sym) == 0 || length(regulon) == 0) {
    warning("Skipping ", st, ": expr_sym/regulon is empty."); next
  }
  
  expr_mat <- as.matrix(expr_sym)
  storage.mode(expr_mat) <- "double"
  
  nes_mat <- viper::viper(eset = expr_mat, regulon = regulon, verbose = FALSE)
  
  studies[[st]]$NES <- nes_mat
  
  message("  NES dims: ", paste(dim(nes_mat), collapse = " x "),
          " | TFs inferred: ", nrow(nes_mat))
}

# NES per study is now available at: studies[[st]]$NES


#### 9) Correlate NES with LT_enrichment per study #############################

# 9.1) Spearman with permutation (wrapped)
calculate_spearman_with_permutation <- function(x, y, n_permutations = NULL, alternative = "two.sided") {
  spearman_test <- cor.test(x, y, method = "spearman")
  r_observed <- spearman_test$estimate
  if (is.na(r_observed)) return(list(r = NA, p = NA))
  
  statistic <- function(data) cor(data[, 1], data[, 2], method = "spearman")
  
  if (is.null(n_permutations)) {
    n <- length(x)
    if (n <= 6)      n_permutations <- factorial(n)
    else if (n <= 30)  n_permutations <- 5000
    else if (n <= 100) n_permutations <- 50000
    else               n_permutations <- 100000
  }
  
  data <- data.frame(x = x, y = y)
  permutation_result <- independence_test(x ~ y, data = data, teststat = "scalar",
                                          distribution = approximate(nresample = n_permutations),
                                          alternative = alternative)
  p_value <- pvalue(permutation_result)
  if (is.na(p_value)) return(list(r = r_observed, p = NA))
  list(r = r_observed, p = p_value)
}

# 9.2) Loop per study and save per-study tables
for (st in names(studies)) {
  message("Correlating NES with LT_enrichment_fct: ", st)
  
  nes_mat <- studies[[st]]$NES   # TF x samples
  meta    <- studies[[st]]$metadata
  if (is.null(nes_mat) || is.null(meta) || nrow(nes_mat) == 0) {
    warning("Skipping ", st, ": NES or metadata is empty."); next
  }
  
  common_samples <- intersect(colnames(nes_mat), rownames(meta))
  if (length(common_samples) < 3) {
    warning("Too few samples in common for study: ", st); next
  }
  nes_mat <- nes_mat[, common_samples, drop = FALSE]
  ranking <- meta[common_samples, "LT_enrichment"]
  
  results <- lapply(rownames(nes_mat), function(tf) {
    res <- calculate_spearman_with_permutation(
      x = as.numeric(nes_mat[tf, ]),
      y = as.numeric(ranking),
      alternative = "two.sided"
    )
    tibble(
      study = st,
      TF = tf,
      rho = as.numeric(res$r),
      p   = as.numeric(res$p)
    )
  }) %>% bind_rows() %>%
    mutate(FDR = p.adjust(p, method = "BH"))
  
  studies[[st]]$NES_cor_results <- results
}

# 9.3) Peek
if (length(studies) > 0) print(head(studies[[1]]$NES_cor_results))


#### 10) Aggregate across studies, consistency & meta-analysis #################

# 10.1) Combine all per-study results
nes_cor_all <- lapply(names(studies), function(st) {
  tbl <- studies[[st]]$NES_cor_results
  if (is.null(tbl)) return(NULL)
  tbl
}) %>% bind_rows() %>%
  mutate(
    direction = case_when(
      !is.na(FDR) & FDR < 0.05 & rho > 0 ~ "activator",
      !is.na(FDR) & FDR < 0.05 & rho < 0 ~ "repressor",
      TRUE ~ "NS"
    ),
    signif = if_else(!is.na(FDR) & FDR < 0.05, TRUE, FALSE)
  ) %>%
  arrange(study, FDR, desc(abs(rho)))

# 10.2) Focus on TFs of interest
tf_interest_symbols <- tf_in_families %>%
  filter(!is.na(TF) & TF != "") %>%
  pull(TF) %>% unique()

nes_cor_interest <- nes_cor_all %>%
  filter(TF %in% tf_interest_symbols) %>%
  arrange(study, FDR, desc(abs(rho)))

# 10.3) Consistency by direction per TF
nes_cor_interest_all <- lapply(names(studies), function(st) {
  tbl <- studies[[st]]$NES_cor_results
  if (is.null(tbl)) return(NULL)
  tbl %>% filter(TF %in% tf_interest_symbols)
}) %>%
  bind_rows() %>%
  mutate(
    direction = case_when(
      !is.na(FDR) & FDR < 0.05 & rho > 0 ~ "activator",
      !is.na(FDR) & FDR < 0.05 & rho < 0 ~ "repressor",
      TRUE ~ "NS"
    ),
    signif = !is.na(FDR) & FDR < 0.05
  ) %>%
  arrange(study, FDR, desc(abs(rho)))

# 10.4) Count studies per direction (wide)
counts_dir <- nes_cor_interest_all %>%
  dplyr::filter(signif) %>%
  dplyr::mutate(dir_bin = dplyr::if_else(rho > 0, "activator", "repressor")) %>%
  dplyr::count(TF, dir_bin, name = "n_studies") %>%
  tidyr::pivot_wider(names_from = dir_bin, values_from = n_studies, values_fill = 0)

# 10.5) Fisher meta-analysis per TF
meta_fisher <- nes_cor_interest_all %>%
  group_by(TF) %>%
  summarise(
    k_studies     = sum(!is.na(p)),
    p_meta        = ifelse(k_studies >= 4,
                           pchisq(-2 * sum(log(p[!is.na(p)])),
                                  df = 2 * k_studies, lower.tail = FALSE),
                           NA_real_),
    n_pos         = sum(rho > 0, na.rm = TRUE),
    n_neg         = sum(rho < 0, na.rm = TRUE),
    n_sig_pos     = sum(signif & rho > 0, na.rm = TRUE),
    n_sig_neg     = sum(signif & rho < 0, na.rm = TRUE),
    rho_median    = median(rho, na.rm = TRUE),
    majority_dir  = case_when(
      n_pos > n_neg ~ "activator",
      n_neg > n_pos ~ "repressor",
      TRUE          ~ "tie"
    ),
    .groups = "drop"
  ) %>%
  mutate(FDR_meta = p.adjust(p_meta, method = "BH")) %>%
  arrange(FDR_meta, desc(k_studies))

library(openxlsx)

write.xlsx(
  meta_fisher,
  file = "Transcription _Factor/Viper_Analysis/meta_fisher_results.xlsx",
  row.names = FALSE
)


meta_fisher_filtered <- meta_fisher %>%
  filter(
    !is.na(FDR_meta),
    FDR_meta < 0.05,
    abs(rho_median) >= 0.5
  )
# 10.6) Candidate TFs
tf_candidates <- meta_fisher %>%
  dplyr::filter(!is.na(FDR_meta) & FDR_meta < 0.05) %>%
  dplyr::arrange(FDR_meta, dplyr::desc(k_studies), dplyr::desc(abs(rho_median))) %>%
  dplyr::pull(TF) %>%
  unique()

# 10.7) Slim table for inspection
tf_candidates_tbl <- meta_fisher %>%
  dplyr::filter(TF %in% tf_candidates) %>%
  dplyr::select(TF, FDR_meta, p_meta, k_studies, n_pos, n_neg, n_sig_pos, n_sig_neg, majority_dir, rho_median) %>%
  dplyr::arrange(FDR_meta)

# 10.8) Stringent direction-consistency rule by k_studies
# Rule: for k=5 -> at least 4 votes to one side; for k=4 -> at least 3
tf_strict_consensus <- meta_fisher %>%
  dplyr::filter(
    # Casos k = 5 → pelo menos 4 votos iguais
    (k_studies == 5 & (n_pos >= 4 | n_neg >= 4)) |
      
      # Casos k = 4 → pelo menos 3 votos iguais
      (k_studies == 4 & (n_pos >= 3 | n_neg >= 3)) |
      
      # NOVO: Casos k = 3 → todos os 3 votos iguais
      (k_studies == 3 & (n_pos == 3 | n_neg == 3))
  ) %>%
  dplyr::mutate(
    consensus_votes = pmax(n_pos, n_neg)
  ) %>%
  dplyr::arrange(
    dplyr::desc(k_studies),
    dplyr::desc(consensus_votes),
    FDR_meta,
    dplyr::desc(abs(rho_median))
  )

write.csv(tf_candidates_tbl, "Transcription _Factor/Viper_Analysis/tf_candidates.csv", row.names = FALSE)

openxlsx::write.xlsx(
  tf_strict_consensus,
  file = "Transcription _Factor/Viper_Analysis/tf_strict_consensus.xlsx",
  asTable = TRUE
)

# 10.8) Candidate table (TF-level) joined with family from metadata_Jaspar
#       - Keep only TFs already in tf_candidates_tbl
#       - Ensure one row per TF (distinct TF–family in metadata)
#       - Fill missing family with "unannotated" for clarity

family_map <- metadata_Jaspar %>%
  dplyr::select(TF, family) %>%
  dplyr::distinct()

tf_candidates_with_family <- tf_strict_consensus %>%
  dplyr::left_join(family_map, by = "TF") %>%
  dplyr::mutate(
    family = dplyr::if_else(is.na(family) | family == "", "unannotated", family)
  ) %>%
  dplyr::relocate(family, .after = TF) %>%
  dplyr::arrange(FDR_meta, dplyr::desc(k_studies), dplyr::desc(abs(rho_median)))

# Save table
write.csv(tf_candidates_with_family, "Transcription _Factor/Viper_Analysis/tf_candidates.csv", row.names = FALSE)


#### 11) Visualization: heatmap of per-study rho (NES vs LT_enrichment) ########

# 11.1) Escolher TFs: usar o strict-consensus (regra k=5→≥4 votos; k=4→≥3)
tf_cols <- tf_strict_consensus %>%
  dplyr::arrange(dplyr::desc(k_studies),
                 dplyr::desc(consensus_votes),
                 FDR_meta,
                 dplyr::desc(abs(rho_median))) %>%
  dplyr::pull(TF) %>%
  unique()

# fallback (se por algum motivo ficar vazio)
if (length(tf_cols) == 0) {
  tf_cols <- meta_fisher %>%
    dplyr::arrange(FDR_meta, dplyr::desc(k_studies)) %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::pull(TF) %>% unique()
}

# 11.2) Colapsar para um rho por TF×estudo (menor FDR)
heat_long <- nes_cor_interest_all %>%
  dplyr::filter(TF %in% tf_cols) %>%
  dplyr::select(TF, study, rho, FDR) %>%
  dplyr::group_by(TF, study) %>%
  dplyr::arrange(FDR, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# 11.3) Ordens dos eixos
#     – TFs ordenados por mediana de rho do meta (ou pelo strict-consensus se preferir)
tf_order <- tf_strict_consensus %>%
  dplyr::filter(TF %in% tf_cols) %>%
  dplyr::arrange(rho_median) %>%            # use "desc(rho_median)" se quiser ativadores no topo
  dplyr::pull(TF)

study_order <- names(studies)

# 11.4) Build a complete TF×study grid and compute label aesthetics
#       - Fill missing combinations so NAs render as grey tiles without labels
#       - Create text labels and contrast-aware text color
df_plot <- heat_long %>%
  mutate(
    study = factor(study, levels = study_order),
    TF    = factor(TF,    levels = tf_order)
  ) %>%
  tidyr::complete(TF, study) %>%                            # ensure full grid
  mutate(
    rho_lbl = ifelse(is.na(rho), "", sprintf("%.2f", rho)), # in-cell text
    txt_col = ifelse(is.na(rho), NA_character_,
                     ifelse(abs(rho) >= 0.4, "white", "black")) # simple contrast
  )

# 11.5) Plot heatmap (white borders; diverging palette; midpoint at 0)
p_heat_vals <- ggplot(df_plot, aes(x = study, y = TF, fill = pmax(pmin(rho, 1), -1))) +
  geom_tile(color = "white", size = 0.5) +   # white border around each cell
  geom_text(aes(label = rho_lbl, color = txt_col), size = 3.2, na.rm = TRUE, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "darkblue", mid = "#f7f7f7", high = "#d7191c",
    midpoint = 0, limits = c(-1, 1), na.value = "#eeeeee",
    name = "Spearman ρ"
  ) +
  labs(x = NULL, y = NULL, title = "Correlação (ρ) NES vs LT_enrichment — TF (Y) × estudo (X)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid  = element_blank(),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )
print(p_heat_vals)

# 11.6) Save figure
ggsave("Transcription _Factor/Viper_Analysis/heatmap_rho_TF_by_study_with_values.png",
       p_heat_vals, width = 8, height = 8, dpi = 300)



##### 12) Intersections of significant TFs across study groups #################

# Helper: get unique TFs with signif == TRUE in a given study (or studies)
get_signif_TFs <- function(data, study_ids) {
  data %>%
    dplyr::filter(study %in% study_ids, signif) %>%
    dplyr::pull(TF) %>%
    unique()
}

# 12.1) Grupo 1: Anjos_2020_T1 ∩ Xie_2019_T1 (Type 1)
signif_anjos_T1 <- get_signif_TFs(nes_cor_interest_all, "Anjos_2020_T1")
signif_xie_T1   <- get_signif_TFs(nes_cor_interest_all, "Xie_2019_T1")

tf_intersect_T1 <- intersect(signif_anjos_T1, signif_xie_T1)

length(signif_anjos_T1)  # só para checar quantos
length(signif_xie_T1)
length(tf_intersect_T1)  # quantos estão na intersecção T1

# 12.2) Grupo 2: Xie_2019_T2 ∩ Fares_2016_T2 ∩ Fares_2014_T2 (Type 2)
signif_xie_T2    <- get_signif_TFs(nes_cor_interest_all, "Xie_2019_T2")
signif_fares2016 <- get_signif_TFs(nes_cor_interest_all, "Fares_2016_T2")
signif_fares2014 <- get_signif_TFs(nes_cor_interest_all, "Fares_2014_T2")

tf_intersect_T2 <- Reduce(
  intersect,
  list(signif_xie_T2, signif_fares2016, signif_fares2014)
)

length(signif_xie_T2)
length(signif_fares2016)
length(signif_fares2014)
length(tf_intersect_T2)  # quantos estão na intersecção T2

# 12.3) Interseção final: TFs que são significativos em T1 e em T2
tf_intersect_all <- intersect(tf_intersect_T1, tf_intersect_T2)
length(tf_intersect_all)

# 12.3b) TFs significativos em pelo menos 4 dos 5 estudos
signif_counts <- nes_cor_interest_all %>%
  dplyr::filter(signif) %>%             # só FDR < 0.05
  dplyr::count(TF, name = "n_studies_signif")

tf_at_least_4of5 <- signif_counts %>%
  dplyr::filter(n_studies_signif >= 4) %>%
  dplyr::pull(TF)

length(tf_at_least_4of5)   # checar quantos dão 4/5 ou 5/5

# 12.4) Tabelas para inspeção / export

# (a) TFs significativos na intersecção de T1
tf_T1_tbl <- nes_cor_interest_all %>%
  dplyr::filter(TF %in% tf_intersect_T1,
                study %in% c("Anjos_2020_T1", "Xie_2019_T1")) %>%
  dplyr::arrange(TF, study)

# (b) TFs significativos na intersecção de T2
tf_T2_tbl <- nes_cor_interest_all %>%
  dplyr::filter(TF %in% tf_intersect_T2,
                study %in% c("Xie_2019_T2", "Fares_2016_T2", "Fares_2014_T2")) %>%
  dplyr::arrange(TF, study)

# (c) TFs que estão na intersecção T1 ∩ T2
tf_all_tbl <- nes_cor_interest_all %>%
  dplyr::filter(TF %in% tf_intersect_all) %>%
  dplyr::arrange(TF, study)

# 12.5) Opcional: salvar como CSV
write.csv(tf_T1_tbl,
          "Transcription _Factor/Viper_Analysis/TF_intersection_T1_Anjos_Xie.csv",
          row.names = FALSE)

write.csv(tf_T2_tbl,
          "Transcription _Factor/Viper_Analysis/TF_intersection_T2_Xie_Fares_Fares.csv",
          row.names = FALSE)

write.csv(tf_all_tbl,
          "Transcription _Factor/Viper_Analysis/TF_intersection_T1_T2_all.csv",
          row.names = FALSE)


tf_viper_counts <- lapply(names(studies), function(st) {
  nes_mat <- studies[[st]]$NES
  if (is.null(nes_mat)) return(NULL)
  
  tibble(
    study      = st,
    n_TF_VIPER = nrow(nes_mat)   # todos os TFs inferidos pelo VIPER
  )
}) %>%
  bind_rows()

tf_viper_counts

tf_interest_symbols <- tf_in_families %>%
  filter(!is.na(TF) & TF != "") %>%
  pull(TF) %>%
  unique()

tf_interest_per_study <- lapply(names(studies), function(st) {
  nes_mat <- studies[[st]]$NES
  if (is.null(nes_mat)) return(NULL)
  
  all_tfs <- rownames(nes_mat)
  tfs_interest_here <- intersect(all_tfs, tf_interest_symbols)
  
  tibble(
    study                = st,
    n_TF_VIPER           = length(all_tfs),
    n_TF_interest_tested = length(tfs_interest_here)
  )
}) %>%
  bind_rows()

tf_interest_per_study

n_interest_tested_overall <- nes_cor_interest_all %>%
  distinct(TF) %>%
  nrow()

n_interest_tested_overall

tf_viper_all <- lapply(names(studies), function(st) {
  nes_mat <- studies[[st]]$NES
  if (is.null(nes_mat)) return(NULL)
  
  tibble(
    study     = st,
    n_TF_all  = nrow(nes_mat)     # total de TFs inferidos pelo VIPER
  )
}) %>% bind_rows()

tf_viper_all

nes_cor_interest_all %>% distinct(TF) %>% nrow()

length(common_samples)
setdiff(colnames(nes_mat), rownames(meta))
setdiff(rownames(meta), colnames(nes_mat))
