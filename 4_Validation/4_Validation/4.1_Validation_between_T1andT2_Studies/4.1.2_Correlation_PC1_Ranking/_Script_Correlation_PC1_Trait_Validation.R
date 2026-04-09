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
  Anjos_2020_T1 = list(
    expr = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Anjos-Afonso_etal_2021/DATA/counts_Final.csv",
    meta = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Anjos-Afonso_etal_2021/DATA/metadata.csv"
  ),
  Xie_2019_T1 = list(    # <-- Renamed from Xie
    expr = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Xie_etal_2019/DATA/counts_Final.csv",
    meta = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Xie_etal_2019/DATA/metadata.csv"
  ),
  Fares_2016_T2 = list(
    expr = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Fares_2017/DATA/counts_Final.csv",
    meta = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Fares_2017/DATA/metadata.csv"
  ),
  Fares_2014_T2 = list(
    expr = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Sauvageau_2014/DATA/counts_Final.csv",
    meta = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Sauvageau_2014/DATA/metadata.csv"
  ),
  Xie_2019_T2 = list(    # <-- Renamed from XieCult
    expr = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Xie_2019/DATA/counts_Final.csv",
    meta = "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Xie_2019/DATA/metadata.csv"
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

# Merge RNA-seq and microarray studies into a single named list.
studies <- c(studies_RNAseq)

# Manually adjusting LT_enrichment_fct for each study using the standardized labels:

studies$Anjos_2020_T1$metadata$LT_enrichment_fct <- c(
  "Upper Enrichment", "Moderate Enrichment", "Low Enrichment", "Moderate Enrichment",
  "Upper Enrichment", "Moderate Enrichment", "Low Enrichment", "Moderate Enrichment",
  "Upper Enrichment", "Moderate Enrichment", "Low Enrichment", "Moderate Enrichment"
)

studies$Xie_2019_T1$metadata$LT_enrichment_fct <- c(
  "Upper Enrichment", "Upper Enrichment", "Upper Enrichment",
  "Moderate Enrichment", "Moderate Enrichment", "Moderate Enrichment",
  "Low Enrichment", "Low Enrichment", "Low Enrichment", "Low Enrichment",
  "Low Enrichment", "Low Enrichment", "Low Enrichment", "Low Enrichment",
  "Low Enrichment"
)

studies$Fares_2014_T2$metadata$LT_enrichment_fct <- c(
  "Upper Enrichment",  "Low Enrichment",       "Moderate Enrichment", "High Enrichment",
  "Upper Enrichment",  "Low Enrichment",       "Moderate Enrichment", "High Enrichment",
  "Upper Enrichment",  "Low Enrichment",       "Moderate Enrichment", "High Enrichment",
  "Upper Enrichment",  "Low Enrichment",       "Moderate Enrichment", "High Enrichment"
)

studies$Fares_2016_T2$metadata$LT_enrichment_fct <- c(
  "Upper Enrichment", "Upper Enrichment",
  "Low Enrichment",   "Low Enrichment", "Low Enrichment",
  "Moderate Enrichment", "Moderate Enrichment", "Moderate Enrichment",
  "Upper Enrichment",
  "Low Enrichment", "Moderate Enrichment", "Moderate Enrichment"
)

studies$Xie_2019_T2$metadata$LT_enrichment_fct <- c(
  "Low Enrichment", "Low Enrichment", "Upper Enrichment", "Upper Enrichment",
  "Low Enrichment", "Low Enrichment", "Upper Enrichment", "Upper Enrichment",
  "Low Enrichment", "Low Enrichment", "Upper Enrichment", "Upper Enrichment",
  "Low Enrichment", "Low Enrichment", "Upper Enrichment", "Upper Enrichment",
  "Low Enrichment", "Low Enrichment", "Upper Enrichment", "Upper Enrichment",
  "Low Enrichment", "Low Enrichment", "Upper Enrichment", "Upper Enrichment"
)

#### 2. Filtering the modules ####
# Description: Create multiple gene lists based on different overlap criteria 
# between cultured and uncultured datasets. Also add a small marker list.

gene_table <- read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Overlap_Uncult_Cult/Gene_Prioritization//priority_lists_long.csv")

desired_lists <- c( "list3_4or5of5")
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
    
    pca_result <- prcomp(exp_matrix, center = TRUE, scale. = FALSE)
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
    file = paste0("Artigo_GBM_VPO/Suplementar/Correlation_", list_name, ".csv"),
    row.names = FALSE
  )
}

#### 6. Plot the Graphs ####

###### 6.1 Plot the Graphs ######
# Build per-study panels (PC1 vs LT enrichment groups) for each gene list,
# add rho/p labels, and combine with patchwork into multi-row figures.

plot_grid <- list()

for (list_name in names(gene_sets)) {
  plot_grid[[list_name]] <- list()
  
  for (study_name in names(studies)) {
    # Skip if PCA results are missing for this study/list
    if (is.null(pca_results_by_list[[study_name]][[list_name]])) next
    
    pc1 <- pca_results_by_list[[study_name]][[list_name]]$pc1
    correlation <- correlation_by_list[[study_name]][[list_name]]
    metadata <- studies[[study_name]]$metadata
    
    # Choose a categorical grouping column for x-axis labels (visual strata)
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
    
    # Assemble plotting frame with ordered factor levels (top to bottom enrichment)
    plot_df <- data.frame(
      PC1 = pc1,
      Group = factor(group_col, levels = c(
        "Upper Enrichment", "High Enrichment", 
        "Moderate Enrichment", "Low Enrichment", "Lower Enrichment"
      ))
    )
    
    # Simple diverging palette for groups (red → blue)
    group_levels <- levels(plot_df$Group)
    group_colors <- setNames(
      colorRampPalette(c("red", "blue"))(length(group_levels)),
      group_levels
    )
    
    # Compose correlation label (rho and p-value)
    rho <- round(correlation$r, 2)
    pval <- ifelse(is.numeric(correlation$p), round(correlation$p, 3), correlation$p)
    label_text <- paste0("\u03C1 = ", rho, "\n", "p = ", pval)
    
    # Number of samples included in this panel
    n_samples <- nrow(plot_df)
    sample_label <- paste0("n = ", n_samples)
    
    # Panel per study (no BM/CB/mPB label anymore)
    p <- ggplot(plot_df, aes(x = Group, y = PC1, color = Group)) +
      geom_point(size = 5, alpha = 0.6) +
      labs(
        title = study_name,            # <- only the study name
        x = "LT Enrichment",
        y = "PC1"
      ) +
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
      # Bottom-left sample size label
      annotate(
        "text",
        x = -Inf, y = -Inf,
        label = sample_label,
        hjust = -0.1, vjust = -1,
        color = "black",
        size = 3.0
      ) +
      scale_color_manual(values = group_colors) +
      theme_minimal(base_size = 14) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x  = element_blank()   # hide x labels to save space
      )
    
    # Store the plot in the grid for this list
    plot_grid[[list_name]][[study_name]] <- p
  }
}

# Quick visual check for a specific list before bulk-saving
wrap_plots(plot_grid$list3_4or5of5, nrow = 1)

###### Combine all plots into a patchwork layout ######
# Save one PNG per list (layout: 6 rows; adjust if you change the number of studies).

for (list_name in names(plot_grid)) {
  message(paste0("Saving plot for ", list_name, " ..."))
  
  final_plot <- wrap_plots(plot_grid[[list_name]], nrow = 1) &
    theme(plot.margin = margin(5, 5, 5, 5))
  
  ggsave(
    filename = paste0("Artigo_GBM_VPO/Suplementar/Correlation_PC1_Trait_", list_name, ".png"),
    plot = final_plot,
    width = 30,
    height = 7,
    units = "cm",
    dpi = 300
  )
  
  message(paste0("Plot '", list_name, "' saved successfully!\n"))
}


###### 6.2 Plot the Scatter Plot ######

#### 6.2.1) Handle p = 0 → cap -log10(p) at 4 
correlation_df2 <- correlation_df %>%
  mutate(
    negLogP = -log10(p),
    negLogP = ifelse(is.infinite(negLogP) | negLogP > 4, 4, negLogP)  # cap at 4
  )

#### 6.2.2) Fixed color palette per Study (keep original Study names) 
pal_studies <- c(
  "Fares_2016_T2" = "#FF96FF",
  "Fares_2014_T2" = "#80E7E9",
  "Xie_2019_T1"   = "#FFD03D",
  "Anjos_2020_T1" = "#FF944C",
  "Xie_2019_T2"   = "#92DF71"
)

#### 6.2.3) Order legend by desired sequence 
correlation_df2$Study <- factor(
  correlation_df2$Study,
  levels = c("Xie_2019_T1","Anjos_2020_T1", "Fares_2014_T2", "Fares_2016_T2", "Xie_2019_T2")
)

med_r <- median(correlation_df2$r, na.rm = TRUE)
#### 6.2.4) Plot (with gentle jitter to reduce overlap) 
p_corr <- ggplot(correlation_df2, aes(x = r, y = negLogP, color = Study)) +
  geom_jitter(
    width  = 0.00, height = 0.10,
    size   = 10, alpha = 0.9, stroke = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype   = "dashed",
    color      = "#D73027",
    linewidth  = 0.8
  ) +
  # Add median label box (no line)
  annotate("label",
           x = 0.505, y = 4.4,   # upper-left area
           label = paste0("Median ρ = ", round(med_r, 1)),
           size = 4,
           label.size = 0.3,
           fill = "white",
           color = "black",
           fontface = "plain",
           family = "Arial",
           label.r = unit(0.15, "lines"),
           hjust = 0, vjust = 1) +
  scale_color_manual(values = pal_studies, name = "Study") +
  scale_x_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4, 1)) +
  labs(
    x = expression("Spearman Correlation (PC1 vs LT-HSC Enrichment Rank)"),
    y = expression(-log[10](FDR)),
    title = ""
  ) +
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    text             = element_text(family = "Arial"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position  = "right",
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 14),
    axis.text.x      = element_text(size = 14),
    axis.text.y      = element_text(size = 14),
    axis.title.x     = element_text(size = 14, face = "bold"),
    axis.title.y     = element_text(size = 14, face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 10, alpha = 1))
  )

print(p_corr)

#### 6.2.5) Save (PNG) 
ggsave(
  "Artigo_GBM_VPO/Suplementar/Correlation_T1_T2_V3.png",
  p_corr,
  width = 15, height = 5, units = "in",
  dpi = 300
)

