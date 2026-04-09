# Title: Principal Component Analysis and Spearman Correlation with Permutation
#
# Description:
# This script performs PCA on variance-stabilized gene expression data from 
# multiple transcriptomic studies of hematopoietic stem cells (HSCs), and calculates 
# the correlation between the first principal component (PC1) and biological traits 
# using Spearman correlation with permutation testing. 
#
# Inputs:
# - Expression matrices (.csv): Raw count data from RNA-seq experiments
# - Metadata files (.csv): Sample annotations for each study
# - Gene list file (.csv): Table of genes of interest for downstream analysis
#
# Outputs:
# - PCA results (PC1 scores)
# - Correlation values (Spearman's rho and permutation-based p-values)
# - Adjusted p-values (Benjamini-Hochberg correction)
#
# Notes:
# - RNA-seq matrices are normalized via DESeq2's VST; microarray matrices are assumed
#   already normalized (RMA) and are only transposed to samples x genes format.
# - PC1 orientation is forced to increase with LT enrichment to keep biological
#   interpretation consistent across studies and lists.
# - All per-study/per-list correlation results are saved as CSVs; plot panels are
#   saved as PNGs using patchwork.


##### 0. Load required libraries #####
# Visualization, stats, normalization, and data wrangling.
library(ggplot2)    # Data visualization
library(ggpubr)     # Publication-ready plots
library(DESeq2)     # Normalization and transformation of RNA-seq data
library(coin)       # Permutation-based statistical tests
library(patchwork)  # For combining ggplot objects
library(tibble)     # Tidyverse tibble utilities
library(dplyr)      # Data manipulation

##### 1. Load the expression and metadata for each dataset #####
# Set the working directory to where inputs/outputs live.
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022")

###### 1.1 Load the functions ######
# Helper to load paired expression/metadata CSVs into a named list.
# - expression_path: counts (RNA-seq) or normalized (microarray); genes in rows, samples in columns
# - metadata_path: sample annotations; samples in rows
load_study_data <- function(expression_path, metadata_path) {
  expData <- read.csv(expression_path, row.names = 1)
  metadata <- read.csv(metadata_path, row.names = 1)
  list(expData = expData, metadata = metadata)
}

# Perform PCA on a pre-filtered, sample x gene matrix and store PC1.
# Expects `study_data$filtered_exp` to exist (set later in pipeline).
perform_pca <- function(study_data) {
  pca_result <- prcomp(study_data$filtered_exp, center = TRUE, scale. = TRUE)
  study_data$pca <- pca_result
  study_data$pc1 <- pca_result$x[, 1]
  return(study_data)
}

# Spearman correlation with permutation testing.
# - Chooses number of permutations adaptively by sample size if not provided.
# - Returns rho and a permutation-based p-value (coin::independence_test).
calculate_spearman_with_permutation <- function(x, y, n_permutations = NULL, alternative = "greater") {
  spearman_test <- cor.test(x, y, method = "spearman")
  r_observed <- spearman_test$estimate
  
  if (is.na(r_observed)) return(list(r = NA, p = NA))
  
  if (is.null(n_permutations)) {
    n <- length(x)
    # Adaptive schedule to cover small-n exactly and scale reasonably with n.
    if (n <= 6) {
      n_permutations <- factorial(n)
    } else if (n <= 30) {
      n_permutations <- 5000
    } else if (n <= 100) {
      n_permutations <- 50000
    } else {
      n_permutations <- 100000
    }
  }
  
  data <- data.frame(x = x, y = y)
  perm_test <- independence_test(x ~ y, data = data, teststat = "scalar",
                                 distribution = approximate(nresample = n_permutations),
                                 alternative = alternative)
  p_value <- pvalue(perm_test)
  
  if (is.na(p_value)) return(list(r = r_observed, p = NA))
  
  return(list(r = r_observed, p = p_value))
}

# Generic correlation wrapper for module eigengenes vs trait (kept for flexibility).
# Here, `MEs` can be any matrix-like (samples x modules), and `samp` a numeric trait vector.
analyze_module_trait_correlations <- function(MEs, samp, alternative = "two.sided") {
  results <- data.frame(
    Module = colnames(MEs),
    Correlation = numeric(ncol(MEs)),
    P_value = numeric(ncol(MEs))
  )
  
  for (i in 1:ncol(MEs)) {
    result <- calculate_spearman_with_permutation(MEs[, i], samp, alternative)
    results$Correlation[i] <- result$r
    results$P_value[i] <- result$p
  }
  
  results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")
  return(results)
}

# VST normalization for RNA-seq (counts) followed by transpose to samples x genes.
# - Keeps columns of metadata aligned to count matrix columns by sample names.
normalize_with_vst <- function(count_data, metadata) {
  metadata <- metadata[colnames(count_data), , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ 1)
  vsd <- varianceStabilizingTransformation(dds)
  vst_matrix <- assay(vsd)
  vst_df <- as.data.frame(t(vst_matrix))  # Transpose to samples x genes
  return(vst_df)
}

###### 1.2 Load data paths for multiple studies ######
# Paths for all RNA-seq studies (counts + metadata).
study_paths_RNAseq <- list(
  Cesana_2018 = list(
    expr = "RNAseq/Cesana_etal_2018/DATA/counts_final.csv",
    meta = "RNAseq/Cesana_etal_2018/DATA/metadata.csv"
  ),
  Papa_2018 = list(
    expr = "RNAseq/Papa_etal_2018/DATA/counts_final.csv",
    meta = "RNAseq/Papa_etal_2018/DATA/metadata.csv"
  ),
  Amon_2018 = list(
    expr = "RNAseq/Amon_etal_2018/DATA/counts_final.csv",
    meta = "RNAseq/Amon_etal_2018/DATA/metadata.csv"
  ),
  Calvanese_2019 = list(
    expr = "RNAseq/Calvanese_2019/DATA/counts_Final.csv",
    meta = "RNAseq/Calvanese_2019/DATA/metadata.csv"
  ),
  Expanded_H = list(
    expr  = "inHouse/2_expanded_set/counts_Final_without_100_RUVg_p06k5.csv",
    meta = "inHouse/2_expanded_set/metadata_without_100.csv"
  ),
  Fresh_H = list(
    expr  = "inHouse/1_fresh_set/counts_Final_RUVg_p07k5 (2).csv", 
    meta = "inHouse/1_fresh_set/metadata.csv"
  )
)


# Paths for all microarray studies (RMA-normalized + metadata).
study_paths_Array <- list(
  Rundberg_2015 = list(
    expr = "Microarray/Rundberg_2015/DATA/normalized_expression_RMA.csv",
    meta = "Microarray/Rundberg_2015/DATA/metadata.csv"
  ),
  Rapin_2012= list(
    expr = "Microarray/Rapin_2014/DATA/normalized_expression_RMA.csv",
    meta = "Microarray/Rapin_2014/DATA/metadata.csv"
  ),
  Prashad_2014= list(
    expr = "Microarray/Prashad_2014/DATA/normalized_expression_RMA.csv",
    meta = "Microarray/Prashad_2014/DATA/metadata.csv"
  ), 
  Nooij_2006 = list(
    expr = "Microarray/Nooij_2006/DATA/normalized_expression_RMA.csv",
    meta = "Microarray/Nooij_2006/DATA/metadata.csv"
  ),
  Barreyro_2012 = list(
    expr = "Microarray/Barreyro_2012/DATA/normalized_expression_RMA.csv",
    meta = "Microarray/Barreyro_2012/DATA/metadata.csv"
  ),
  Gentles_2010 = list(
    expr = "Microarray/Gentles_2010/DATA/normalized_expression_RMA.csv",
    meta = "Microarray/Gentles_2010/DATA/metadata.csv"
  ),
  Souyri_2019 = list(
    expr = "Microarray/Souyri_2019/DATA/normalized_expression_RMA.csv",
    meta = "Microarray/Souyri_2019/DATA/metadata.csv"
  ),
  NI_2005 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/NI_2005/DATA/normalized_expression_RMA.csv",
    meta = "Microarray/Souyri_2019/DATA/metadata.csv"
  )
)

###### 1.3 Normalize data for each study ######
# For RNA-seq: load counts + metadata and VST-normalize, yielding samples x genes.
studies_RNAseq <- lapply(study_paths_RNAseq, function(paths) {
  load_study_data(paths$expr, paths$meta)
})

for (study_name in names(studies_RNAseq)) {
  message("Normalizing study: ", study_name)
  count_data <- studies_RNAseq[[study_name]]$expData
  metadata <- studies_RNAseq[[study_name]]$metadata
  normalized_result <- normalize_with_vst(count_data, metadata)
  studies_RNAseq[[study_name]]$normalized_result <- normalized_result
}

# For microarray: load normalized matrices and only transpose to samples x genes.
studies_Array <- lapply(study_paths_Array, function(paths) {
  load_study_data(paths$expr, paths$meta)
})

for (study_name in names(studies_Array)) {
  message("Transposing data for study: ", study_name)
  studies_Array[[study_name]]$normalized_result <- t(studies_Array[[study_name]]$expData)
}

# Merge RNA-seq and microarray studies into a single named list.
studies <- c(studies_RNAseq, studies_Array)

#### 2. Filtering the modules ####
# Description: Create multiple gene lists based on different overlap criteria 
# between cultured and uncultured datasets. Also add a small marker list.

gene_table <- read.csv("Correlation/priority_lists_long.csv")

desired_lists <- c("list1_5of5", "list3_4or5of5")
gene_table$Ensembl <- as.character(gene_table$Ensembl)

gene_sets <- gene_table %>%
  dplyr::filter(List %in% desired_lists, !is.na(Ensembl), nzchar(Ensembl)) %>%
  dplyr::group_by(List) %>%
  dplyr::summarise(genes = list(unique(Ensembl)), .groups = "drop") %>%
  { setNames(.$genes, .$List) }

# Filter each study's expression matrix to each gene set.
# Result: filtered_expression_by_list[[study]][[list]] = samples x genes matrix (subset of columns).
filtered_expression_by_list <- list()

for (study_name in names(studies)) {
  expr_matrix <- studies[[study_name]]$normalized_result
  filtered_expression_by_list[[study_name]] <- list()
  
  for (list_name in names(gene_sets)) {
    genes_of_interest <- gene_sets[[list_name]]
    common_genes <- intersect(genes_of_interest, colnames(expr_matrix))
    filtered_expr <- expr_matrix[, common_genes, drop = FALSE]
    filtered_expression_by_list[[study_name]][[list_name]] <- filtered_expr
  }
}

# QC: number of genes kept per study and per list.
sapply(filtered_expression_by_list, function(x) {
  sapply(x, ncol)
})

# Build a compact matrix (lists x studies) with retained gene counts and save.
gene_counts_matrix <- sapply(filtered_expression_by_list, function(x) {
  sapply(x, ncol)
})
gene_counts_df <- as.data.frame(t(gene_counts_matrix))
write.csv(gene_counts_df, "Correlation/Gene_Counts_Per_List_Per_Study.csv")

#### 3. Extract biological trait vector ####
# Description: Extract the LT-HSC enrichment vector from metadata.
# - Prefer 'LT' if present; otherwise use 'LT_enrichment'.
# - Store a named vector per study in `studies[[study]]$samp`.

for (study_name in names(studies)) {
  metadata <- studies[[study_name]]$metadata
  
  if ("LT" %in% colnames(metadata)) {
    samp_vector <- metadata$LT
  } else if ("LT_enrichment" %in% colnames(metadata)) {
    samp_vector <- metadata$LT_enrichment
  } else {
    warning(paste("Neither 'LT' nor 'LT_enrichment' found in metadata for study:", study_name))
    next
  }
  
  # Name the trait vector by sample IDs to ensure alignment later on.
  names(samp_vector) <- rownames(metadata)
  
  # Store for downstream correlation/PC1 orientation.
  studies[[study_name]]$samp <- samp_vector
}

# Quick peek at the first 3 studies' trait vectors.
lapply(studies[1:3], function(x) head(x$samp))

#### 4. PCA calculation per gene set ####
# Description: Run PCA and extract PC1 for each gene set and study.
# - Stores the full `prcomp` object, PC1 scores, and the variance explained by all PCs.

pca_results_by_list <- list()

for (study_name in names(filtered_expression_by_list)) {
  pca_results_by_list[[study_name]] <- list()
  
  for (list_name in names(filtered_expression_by_list[[study_name]])) {
    exp_matrix <- filtered_expression_by_list[[study_name]][[list_name]]
    
    if (ncol(exp_matrix) < 2) {
      warning(paste("Not enough genes for PCA in", study_name, "-", list_name))
      next
    }
    
    pca_result <- prcomp(exp_matrix, center = TRUE, scale. = TRUE)
    pc1 <- pca_result$x[, 1]
    variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    
    pca_results_by_list[[study_name]][[list_name]] <- list(
      pca = pca_result,
      pc1 = pc1,
      variance_explained = variance_explained
    )
  }
}

###### 4.1 Force PC1 to increase with LT enrichment ######
# Rule: If the median PC1 value among the highest LT group is negative, flip PC1.
# This ensures that "higher PC1" consistently means "more LT-like".

pc1_orientation_log <- data.frame(
  Study = character(),
  List = character(),
  LT_high_level = numeric(),
  n_high = integer(),
  median_pc1_high = numeric(),
  flipped = logical(),
  stringsAsFactors = FALSE
)

for (study_name in names(pca_results_by_list)) {
  # Get trait vector (LT or LT_enrichment) previously stored
  if (!("samp" %in% names(studies[[study_name]]))) next
  trait <- studies[[study_name]]$samp
  
  # Skip if trait is missing
  if (is.null(trait) || all(is.na(trait))) next
  
  for (list_name in names(pca_results_by_list[[study_name]])) {
    # PC1 scores (names must be sample IDs)
    pc1 <- pca_results_by_list[[study_name]][[list_name]]$pc1
    if (is.null(pc1) || length(pc1) == 0) next
    
    # Align trait to PC1 sample order and drop NAs
    trait_aligned <- trait[names(pc1)]
    keep <- !is.na(trait_aligned) & !is.na(pc1)
    if (sum(keep) < 2) next
    
    pc1 <- pc1[keep]
    trait_aligned <- trait_aligned[keep]
    
    # Identify LT-high group (max level in the aligned vector)
    max_level <- max(trait_aligned, na.rm = TRUE)
    is_high <- trait_aligned == max_level
    if (sum(is_high) == 0) next
    
    # Median PC1 among LT-high samples
    med_high <- median(pc1[is_high], na.rm = TRUE)
    
    # Flip rule: if median is negative, multiply PC1 by -1
    flipped <- FALSE
    if (is.finite(med_high) && med_high < 0) {
      pca_results_by_list[[study_name]][[list_name]]$pc1 <- -1 * pca_results_by_list[[study_name]][[list_name]]$pc1
      
      # (Optional) keep PC1 orientation consistent inside the stored PCA object:
      # p <- pca_results_by_list[[study_name]][[list_name]]$pca
      # p$x[, 1] <- -p$x[, 1]
      # p$rotation[, 1] <- -p$rotation[, 1]
      # pca_results_by_list[[study_name]][[list_name]]$pca <- p
      
      flipped <- TRUE
    }
    
    # Minimal log for reproducibility
    pc1_orientation_log <- rbind(
      pc1_orientation_log,
      data.frame(
        Study = study_name,
        List = list_name,
        LT_high_level = max_level,
        n_high = sum(is_high),
        median_pc1_high = med_high,
        flipped = flipped,
        stringsAsFactors = FALSE
      )
    )
  }
}

# Inspect what happened
print(pc1_orientation_log)

###### 4.2 Flip PC1 and loadings by LT-high rule ######
# Same principle as 4.1, but now flips BOTH scores and loadings inside the `prcomp`
# object to keep internal consistency (PC1 direction, rotation signs).

pc1_orientation_log <- tibble(
  Study = character(), List = character(),
  LT_high_level = numeric(), n_high = integer(),
  median_pc1_high = numeric(), flipped = logical()
)

for (study_name in names(pca_results_by_list)) {
  if (!("samp" %in% names(studies[[study_name]]))) next
  trait <- studies[[study_name]]$samp
  if (is.null(trait) || all(is.na(trait))) next
  
  for (list_name in names(pca_results_by_list[[study_name]])) {
    res <- pca_results_by_list[[study_name]][[list_name]]
    if (is.null(res) || is.null(res$pca)) next
    
    p <- res$pca
    # Align trait to PCA sample order
    trait_aligned <- trait[rownames(p$x)]
    keep <- !is.na(trait_aligned)
    if (sum(keep) < 2) next
    
    max_level <- max(trait_aligned[keep], na.rm = TRUE)
    is_high <- trait_aligned == max_level
    if (sum(is_high, na.rm = TRUE) == 0) next
    
    med_high <- median(p$x[is_high, 1], na.rm = TRUE)
    flipped <- FALSE
    
    # Flip scores AND loadings if LT-high median PC1 is negative
    if (is.finite(med_high) && med_high < 0) {
      p$x[, 1]        <- -p$x[, 1]
      p$rotation[, 1] <- -p$rotation[, 1]
      flipped <- TRUE
    }
    
    # Store back (and keep pc1 synchronized)
    res$pca <- p
    res$pc1 <- p$x[, 1]
    pca_results_by_list[[study_name]][[list_name]] <- res
    
    pc1_orientation_log <- bind_rows(
      pc1_orientation_log,
      tibble(Study = study_name, List = list_name,
             LT_high_level = max_level, n_high = sum(is_high, na.rm = TRUE),
             median_pc1_high = med_high, flipped = flipped)
    )
  }
}

# Sanity check: inspect PC1 vectors from the first three studies and first three lists.
lapply(names(pca_results_by_list)[1:3], function(study) {
  lapply(names(pca_results_by_list[[study]])[1:3], function(list_name) {
    head(pca_results_by_list[[study]][[list_name]]$pc1)
  })
})

# Build a long table with percent variance explained by each PC, per study/list.
# Useful to inspect whether PC1 captures a substantial fraction of variance.
variance_table <- do.call(rbind, lapply(names(pca_results_by_list), function(study) {
  do.call(rbind, lapply(names(pca_results_by_list[[study]]), function(list_name) {
    variance_vec <- pca_results_by_list[[study]][[list_name]]$variance_explained
    data.frame(
      Study = study,
      List = list_name,
      PC = paste0("PC", seq_along(variance_vec)),
      Variance_Explained = round(variance_vec * 100, 2)  # in percentage
    )
  }))
}))

# Persist the PC variance table for QA/reports.
write.csv(variance_table, "Correlation/All_PC_Variances.csv", row.names = FALSE)

#### 5. Correlation analysis between PC1 and biological trait ####
# Description: For each study/list, correlate PC1 with the trait vector (LT),
# using permutation-based p-values. Results are saved per list.

correlation_by_list <- list()

for (study_name in names(pca_results_by_list)) {
  samp_vector <- studies[[study_name]]$samp
  correlation_by_list[[study_name]] <- list()
  
  for (list_name in names(pca_results_by_list[[study_name]])) {
    pc1_scores <- pca_results_by_list[[study_name]][[list_name]]$pc1
    
    # Align sample names (intersection ensures identical ordering)
    common_samples <- intersect(names(pc1_scores), names(samp_vector))
    pc1_scores <- pc1_scores[common_samples]
    trait_vector <- samp_vector[common_samples]
    
    # Run correlation with permutation (two-sided is safer given flipping)
    result <- calculate_spearman_with_permutation(pc1_scores, trait_vector, alternative = "two.sided")
    
    # Store result for later export/plotting
    correlation_by_list[[study_name]][[list_name]] <- result
  }
}

# Quick look at the first few results for sanity.
lapply(names(correlation_by_list)[1:3], function(study) {
  lapply(names(correlation_by_list[[study]])[1:3], function(list_name) {
    correlation_by_list[[study]][[list_name]]
  })
})

# Export per-list correlation summaries as CSV.
for (list_name in names(gene_sets)) {
  result_rows <- list()
  i <- 1
  
  for (study_name in names(correlation_by_list)) {
    if (!is.null(correlation_by_list[[study_name]][[list_name]])) {
      res <- correlation_by_list[[study_name]][[list_name]]
      
      result_rows[[i]] <- data.frame(
        Study = study_name,
        r = round(res$r, 4),
        p = res$p,
        stringsAsFactors = FALSE
      )
      i <- i + 1
    }
  }
  
  correlation_df <- do.call(rbind, result_rows)
  
  # Save CSV (assuming folder already exists)
  write.csv(
    correlation_df,
    file = paste0("Correlation/Correlation_", list_name, ".csv"),
    row.names = FALSE
  )
}

#### 6. Plot the Graphs ####
# Build per-study scatter-like panels (PC1 vs LT groups) for each gene list,
# add rho/p labels, and combine with patchwork into multi-row figures.

# Helper to derive a compact source label ("BM", "CB", "BM+CB", "Mixed")
# Helper to derive a compact source label across the samples in a panel.
# Returns one of: "BM", "CB", "FL", "mPB", "PL", "EB", or combinations like "CB + FL".
# Falls back to "Mixed" when no reliable match is found.
get_source_label <- function(meta, sample_ids) {
  # Candidate columns that may contain the cell source information
  source_cols <- c("Cell_Source", "CellSource", "Source", "Tissue", "Origin", "Cell_Source_simplified")
  existing <- source_cols[source_cols %in% colnames(meta)]
  if (length(existing) == 0) return("Mixed")
  
  # Restrict to the samples included in the panel (ensures alignment)
  src <- meta[sample_ids, existing[1], drop = TRUE]
  src_chr <- as.character(src)
  
  # Vector of abbreviations (one per sample)
  abbr <- rep(NA_character_, length(src_chr))
  
  # ---- Matchers (case-insensitive), covering spelling variations ----
  # Bone Marrow
  is_BM <- grepl("bone\\s*marrow|\\bBM\\b|marrow", src_chr, ignore.case = TRUE)
  # Cord Blood
  is_CB <- grepl("cord\\s*blood|\\bCB\\b", src_chr, ignore.case = TRUE)
  # Fetal Liver (covers 'fetal' and 'foetal')
  is_FL <- grepl("f[oe]tal\\s*liver|\\bFL\\b", src_chr, ignore.case = TRUE)
  # Mobilized Peripheral Blood (mobilized/mobilised; also accepts 'mPB' and 'mpb')
  is_mPB <- grepl("(mobiliz|mobilis).*(peripheral).*(blood)|\\bmpb\\b|\\bmPB\\b", src_chr, ignore.case = TRUE)
  # Placenta
  is_PL <- grepl("placenta|\\bPL\\b", src_chr, ignore.case = TRUE)
  # Embryoid Body
  is_EB <- grepl("embryoid\\s*body|\\bEB\\b", src_chr, ignore.case = TRUE)
  
  # Assignment (if multiple matchers hit the same sample, keep the first match in the order below)
  abbr[is_BM]  <- "BM"
  abbr[is_CB]  <- ifelse(is.na(abbr[is_CB]), "CB", abbr[is_CB])
  abbr[is_FL]  <- ifelse(is.na(abbr[is_FL]), "FL", abbr[is_FL])
  abbr[is_mPB] <- ifelse(is.na(abbr[is_mPB]), "mPB", abbr[is_mPB])
  abbr[is_PL]  <- ifelse(is.na(abbr[is_PL]), "PL", abbr[is_PL])
  abbr[is_EB]  <- ifelse(is.na(abbr[is_EB]), "EB", abbr[is_EB])
  
  uniq <- unique(na.omit(abbr))
  if (length(uniq) == 0) return("Mixed")
  if (length(uniq) == 1) return(uniq)
  
  # Sort according to a canonical order before collapsing (keeps labels consistent across panels)
  order_map <- c("BM" = 1, "CB" = 2, "FL" = 3, "mPB" = 4, "PL" = 5, "EB" = 6)
  uniq <- uniq[order_map[uniq] |> order(na.last = TRUE)]
  paste(uniq, collapse = " + ")
}



plot_grid <- list()

for (list_name in names(gene_sets)) {
  plot_grid[[list_name]] <- list()
  
  for (study_name in names(studies)) {
    if (is.null(pca_results_by_list[[study_name]][[list_name]])) next
    
    pc1 <- pca_results_by_list[[study_name]][[list_name]]$pc1
    correlation <- correlation_by_list[[study_name]][[list_name]]
    metadata <- studies[[study_name]]$metadata
    
    # Choose a categorical grouping column for x-axis labels (visual strata).
    if ("Enrichment" %in% colnames(metadata)) {
      group_col <- metadata$Enrichment
    } else if ("LT_enrichment_fct" %in% colnames(metadata)) {
      group_col <- metadata$LT_enrichment_fct
    } else {
      warning(paste("No valid group column found in", study_name))
      next
    }
    
    names(group_col) <- rownames(metadata)
    
    # Align sample names between PC1 and grouping column
    common_samples <- intersect(names(pc1), names(group_col))
    pc1 <- pc1[common_samples]
    group_col <- group_col[common_samples]
    
    source_label <- get_source_label(metadata, common_samples)
    
    # Assemble plotting frame with ordered factor levels (top to bottom enrichment).
    plot_df <- data.frame(
      PC1 = pc1,
      Group = factor(group_col, levels = c(
        "Upper Enrichment", "High Enrichment", 
        "Moderate Enrichment", "Low Enrichment", "Lower Enrichment"
      ))
    )
    
    # Simple diverging palette for groups (red→blue).
    group_levels <- levels(plot_df$Group)
    group_colors <- setNames(colorRampPalette(c("red", "blue"))(length(group_levels)), group_levels)
    
    # Compose correlation label (rho and p-value)
    rho <- round(correlation$r, 2)
    pval <- ifelse(is.numeric(correlation$p), round(correlation$p, 3), correlation$p)
    label_text <- paste0("\u03C1 = ", rho, "\n", "p = ", pval)
    
    # Number of samples included in this panel.
    n_samples <- nrow(plot_df)
    sample_label <- paste0("n = ", n_samples)
    
    # Panel per study.
    p <- ggplot(plot_df, aes(x = Group, y = PC1, color = Group)) +
      geom_point(size = 5, alpha = 0.6) +
      labs(
        title = paste0(study_name, " | ", source_label),
        x = "LT Enrichment", y = "PC1"
      )+
      # Top-right correlation label
      geom_label(
        aes(x = Inf, y = Inf),
        label = label_text,
        color = "black",
        fill = "white",
        label.size = 0.2,
        label.r = unit(0.1, "lines"),
        size = 3,
        fontface = "bold",
        hjust = 1.1, vjust = 1.1,
        inherit.aes = FALSE
      ) +
      # Bottom-left sample size label (kept unobtrusive)
      annotate("text",
               x = -Inf, y = -Inf,
               label = sample_label,
               hjust = -0.1, vjust = -1,
               color = "black",
               size = 3.0) +
      scale_color_manual(values = group_colors) +
      theme_minimal(base_size = 14) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 10),  # Smaller x-axis title
        axis.title.y = element_text(size = 10),  # Smaller y-axis title
        axis.text.x = element_blank()            # Hide long x labels to save space
      )
    
    # Store the plot in the grid for this list.
    plot_grid[[list_name]][[study_name]] <- p
  }
}

# Quick visual check for a specific list before bulk-saving.
wrap_plots(plot_grid$list1_5of5, nrow = 8)

###### Combine all plots into a patchwork layout ######
# Save one PNG per list (layout: 6 rows; adjust if you change the number of studies).

for (list_name in names(plot_grid)) {
  message(paste0("Saving plot for ", list_name, " ..."))
  
  final_plot <- wrap_plots(plot_grid[[list_name]], nrow = 6) &
    theme(plot.margin = margin(5, 5, 5, 5))
  
  ggsave(
    filename = paste0("Correlation/Correlation_PC1_Trait_", list_name, ".png"),
    plot = final_plot,
    width = 21,
    height = 29.7,
    units = "cm",
    dpi = 300
  )
  
  message(paste0("Plot '", list_name, "' saved successfully!\n"))
}

# Version with alphabetical tags for panel referencing in manuscripts/reviews.
for (list_name in names(plot_grid)) {
  message(paste0("Saving plot for ", list_name, " ..."))
  
  final_plot <- wrap_plots(plot_grid[[list_name]], nrow = 6) +
    plot_annotation(tag_levels = "A") & 
    theme(plot.margin = margin(5, 5, 5, 5))
  
  ggsave(
    filename = paste0("Correlation/Letter_Correlation_PC1_Trait_", list_name, ".png"),
    plot = final_plot,
    width = 25,
    height = 33,
    units = "cm",
    dpi = 300
  )
  
  message(paste0("Plot '", list_name, "' saved successfully!\n"))
}

plot_grid <- list()

for (list_name in names(gene_sets)) {
  plot_grid[[list_name]] <- list()
  n_map <- c()   # <<--- NEW: holds n per study for this list
  
  for (study_name in names(studies)) {
    if (is.null(pca_results_by_list[[study_name]][[list_name]])) next
    
    pc1 <- pca_results_by_list[[study_name]][[list_name]]$pc1
    correlation <- correlation_by_list[[study_name]][[list_name]]
    metadata <- studies[[study_name]]$metadata
    
    # Choose a categorical grouping column for x-axis labels (visual strata).
    if ("Enrichment" %in% colnames(metadata)) {
      group_col <- metadata$Enrichment
    } else if ("LT_enrichment_fct" %in% colnames(metadata)) {
      group_col <- metadata$LT_enrichment_fct
    } else {
      warning(paste("No valid group column found in", study_name))
      next
    }
    
    names(group_col) <- rownames(metadata)
    
    # Align sample names between PC1 and grouping column
    common_samples <- intersect(names(pc1), names(group_col))
    pc1 <- pc1[common_samples]
    group_col <- group_col[common_samples]
    
    source_label <- get_source_label(metadata, common_samples)
    
    # Assemble plotting frame
    plot_df <- data.frame(
      PC1 = pc1,
      Group = factor(group_col, levels = c(
        "Upper Enrichment", "High Enrichment", 
        "Moderate Enrichment", "Low Enrichment", "Lower Enrichment"
      ))
    )
    
    group_levels <- levels(plot_df$Group)
    group_colors <- setNames(colorRampPalette(c("red", "blue"))(length(group_levels)), group_levels)
    
    rho <- round(correlation$r, 2)
    pval <- ifelse(is.numeric(correlation$p), round(correlation$p, 3), correlation$p)
    label_text <- paste0("\u03C1 = ", rho, "\n", "p = ", pval)
    
    n_samples <- nrow(plot_df)
    n_map[study_name] <- n_samples    # <<--- NEW: record n for this study
    sample_label <- paste0("n = ", n_samples)
    
    p <- ggplot(plot_df, aes(x = Group, y = PC1, color = Group)) +
      geom_point(size = 5, alpha = 0.6) +
      labs(
        title = paste0(study_name, " | ", source_label),
        x = "LT Enrichment", y = "PC1"
      ) +
      geom_label(
        aes(x = Inf, y = Inf),
        label = label_text,
        color = "black",
        fill = "white",
        label.size = 0.2,
        label.r = unit(0.1, "lines"),
        size = 3,
        fontface = "bold",
        hjust = 1.1, vjust = 1.1,
        inherit.aes = FALSE
      ) +
      annotate("text",
               x = -Inf, y = -Inf,
               label = sample_label,
               hjust = -0.1, vjust = -1,
               color = "black",
               size = 3.0) +
      scale_color_manual(values = group_colors) +
      theme_minimal(base_size = 14) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank()
      )
    
    plot_grid[[list_name]][[study_name]] <- p
  }
  
  # === NEW: reorder panels by n (descending) for this list ===
  if (length(plot_grid[[list_name]]) > 0) {
    # keep only n’s for studies that actually produced a panel
    n_vec <- n_map[names(plot_grid[[list_name]])]
    ord <- order(n_vec, decreasing = TRUE, na.last = TRUE)
    plot_grid[[list_name]] <- plot_grid[[list_name]][ord]
    
    # (opcional) salvar a tabela de n por estudo já ordenada
    n_df <- data.frame(Study = names(n_vec), n = as.integer(n_vec), row.names = NULL)
    n_df <- n_df[order(n_df$n, decreasing = TRUE), ]
    write.csv(n_df, paste0("Correlation/n_per_study_", list_name, ".csv"), row.names = FALSE)
  }
}

# Quick visual check for a specific list before bulk-saving.
wrap_plots(plot_grid$list1_5of5, nrow = 8)

###### Combine all plots into a patchwork layout ######
# Save one PNG per list (layout: 6 rows; adjust if you change the number of studies).

for (list_name in names(plot_grid)) {
  message(paste0("Saving plot for ", list_name, " ..."))
  
  final_plot <- wrap_plots(plot_grid[[list_name]], nrow = 6) &
    theme(plot.margin = margin(5, 5, 5, 5))
  
  ggsave(
    filename = paste0("Correlation/Correlation_PC1_Trait_", list_name, ".png"),
    plot = final_plot,
    width = 21,
    height = 29.7,
    units = "cm",
    dpi = 300
  )
  
  message(paste0("Plot '", list_name, "' saved successfully!\n"))
}
