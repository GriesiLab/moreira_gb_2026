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
  pca_result <- prcomp(study_data$filtered_exp, center = TRUE, scale. = FALSE)
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

# 0) Build shared random sets from the global universe (common_genes)
build_shared_random_sets <- function(gene_universe, sizes = c(5, 50, 500), n_reps = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  gene_universe <- unique(as.character(gene_universe))
  stopifnot(length(gene_universe) > 0)
  
  shared_sets <- list()
  for (sz in sizes) {
    if (sz > length(gene_universe)) {
      stop(sprintf("Requested size %d exceeds universe size %d.", sz, length(gene_universe)))
    }
    shared_sets[[as.character(sz)]] <- replicate(
      n_reps, sample(gene_universe, size = sz, replace = FALSE), simplify = FALSE
    )
  }
  return(shared_sets)  # e.g., shared_sets[["5"]][[rep_i]] is a character vector of genes
}


run_null_distribution_with_shared_sets <- function(study_name,
                                                   expr_base,            # samples x genes (normalized_filtered)
                                                   trait_vector,         # numeric vector named by samples
                                                   shared_sets,          # output of build_shared_random_sets()
                                                   n_permutations = NULL,
                                                   alternative = "two.sided",
                                                   return_genes = FALSE) {
  
  if (is.null(expr_base) || ncol(expr_base) < 2) {
    stop(sprintf("[%s] expr_base must be a samples x genes matrix with >=2 genes.", study_name))
  }
  
  # All shared sets must be contained in expr_base columns (they came from common_genes)
  pool <- colnames(expr_base)
  for (sz in names(shared_sets)) {
    ok <- vapply(shared_sets[[sz]], function(gv) all(gv %in% pool), logical(1))
    if (!all(ok)) {
      bad <- which(!ok)[1]
      stop(sprintf("[%s] Shared set of size %s (rep %d) contains genes not in expr_base.",
                   study_name, sz, bad))
    }
  }
  
  rows <- list()
  idx <- 1L
  
  for (sz in names(shared_sets)) {
    reps_list <- shared_sets[[sz]]
    for (rep_i in seq_along(reps_list)) {
      gset <- reps_list[[rep_i]]
      Xsub <- expr_base[, gset, drop = FALSE]
      
      res <- analyze_pc1_correlation(
        study_name     = study_name,
        expr_matrix    = Xsub,
        trait_vector   = trait_vector,
        n_permutations = n_permutations,  # NULL -> adaptive inside your function
        alternative    = alternative,
        return_genes   = return_genes
      )
      
      res$Subset <- paste0("Random_", sz)
      res$Rep    <- rep_i
      rows[[idx]] <- res
      idx <- idx + 1L
    }
  }
  
  out <- do.call(rbind, rows)
  cols <- c("Study_ID", "Subset", "Rep", setdiff(colnames(out), c("Study_ID", "Subset", "Rep")))
  out <- out[, cols]
  rownames(out) <- NULL
  out
}


analyze_pc1_correlation <- function(study_name,
                                    expr_matrix,         # samples x genes
                                    trait_vector,        # numeric vector named by samples
                                    n_permutations = NULL,
                                    alternative = "two.sided",
                                    return_genes = TRUE  # paste gene IDs into 'Genes' column
) {
  # ------------------ 0) Basic guards ------------------
  if (is.null(expr_matrix) || nrow(expr_matrix) < 3 || ncol(expr_matrix) < 2) {
    warning(sprintf("[%s] Not enough data: need >=3 samples and >=2 genes.", study_name))
    return(data.frame(
      Study_ID    = study_name,
      N_genes     = ifelse(is.null(expr_matrix), 0L, ncol(expr_matrix)),
      Genes       = NA_character_,
      Rho         = NA_real_,
      P_perm      = NA_real_,
      n_samples   = ifelse(is.null(expr_matrix), 0L, nrow(expr_matrix)),
      PC1_varExp  = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(trait_vector) || all(is.na(trait_vector))) {
    warning(sprintf("[%s] Trait vector is missing or all NA.", study_name))
    return(data.frame(
      Study_ID    = study_name,
      N_genes     = ncol(expr_matrix),
      Genes       = if (return_genes) paste(colnames(expr_matrix), collapse = ",") else NA_character_,
      Rho         = NA_real_,
      P_perm      = NA_real_,
      n_samples   = nrow(expr_matrix),
      PC1_varExp  = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # ------------------ 1) Align samples ----------------
  common_samples <- intersect(rownames(expr_matrix), names(trait_vector))
  if (length(common_samples) < 3) {
    warning(sprintf("[%s] Less than 3 common samples after alignment.", study_name))
    return(data.frame(
      Study_ID    = study_name,
      N_genes     = ncol(expr_matrix),
      Genes       = if (return_genes) paste(colnames(expr_matrix), collapse = ",") else NA_character_,
      Rho         = NA_real_,
      P_perm      = NA_real_,
      n_samples   = length(common_samples),
      PC1_varExp  = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  X <- expr_matrix[common_samples, , drop = FALSE]
  trait <- as.numeric(trait_vector[common_samples])
  
  keep <- is.finite(trait)
  if (sum(keep) < 3) {
    warning(sprintf("[%s] < 3 non-NA samples in trait after filtering.", study_name))
    return(data.frame(
      Study_ID    = study_name,
      N_genes     = ncol(X),
      Genes       = if (return_genes) paste(colnames(X), collapse = ",") else NA_character_,
      Rho         = NA_real_,
      P_perm      = NA_real_,
      n_samples   = sum(keep),
      PC1_varExp  = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  X <- X[keep, , drop = FALSE]
  trait <- trait[keep]
  
  # ------------------ 2) PCA --------------------------
  pca <- tryCatch(
    expr = prcomp(X, center = TRUE, scale. = FALSE),
    error = function(e) NULL
  )
  if (is.null(pca)) {
    warning(sprintf("[%s] PCA failed.", study_name))
    return(data.frame(
      Study_ID    = study_name,
      N_genes     = ncol(X),
      Genes       = if (return_genes) paste(colnames(X), collapse = ",") else NA_character_,
      Rho         = NA_real_,
      P_perm      = NA_real_,
      n_samples   = nrow(X),
      PC1_varExp  = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  pc1_scores <- pca$x[, 1]
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1_var <- if (length(var_explained) >= 1) var_explained[1] else NA_real_
  
  # ------------------ 3) Orient PC1 -------------------
  rho_sign <- suppressWarnings(cor(pc1_scores, trait, method = "spearman", use = "complete.obs"))
  if (is.finite(rho_sign) && rho_sign < 0) {
    pc1_scores <- -pc1_scores
    # Optionally keep prcomp loadings/scores consistent if you'll reuse them:
    # pca$x[, 1]        <- -pca$x[, 1]
    # pca$rotation[, 1] <- -pca$rotation[, 1]
  }
  
  # ------------------ 4) Correlation + perm p ----------
  res <- calculate_spearman_with_permutation(
    x = pc1_scores,
    y = trait,
    n_permutations = n_permutations,  # NULL -> adaptive inside your function
    alternative = alternative
  )
  
  # ------------------ 5) Output row --------------------
  out <- data.frame(
    Study_ID    = study_name,
    N_genes     = ncol(X),
    Genes       = if (return_genes) paste(colnames(X), collapse = ",") else NA_character_,
    Rho         = unname(ifelse(is.null(res$r), NA_real_, res$r)),
    P_perm      = unname(ifelse(is.null(res$p), NA_real_, res$p)),
    n_samples   = nrow(X),
    PC1_varExp  = round(100 * pc1_var, 2),  # percent
    stringsAsFactors = FALSE
  )
  return(out)
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
    expr = "4a_Panel_validation_inSilico/RNAseq/Cesana_etal_2018/DATA/counts_final.csv",
    meta = "4a_Panel_validation_inSilico/RNAseq/Cesana_etal_2018/DATA/metadata.csv"
  ),
  Papa_2018 = list(
    expr = "4a_Panel_validation_inSilico/RNAseq/Papa_etal_2018/DATA/counts_final.csv",
    meta = "4a_Panel_validation_inSilico/RNAseq/Papa_etal_2018/DATA/metadata.csv"
  ),
  Amon_2018 = list(
    expr = "4a_Panel_validation_inSilico/RNAseq/Amon_etal_2018/DATA/counts_final.csv",
    meta = "4a_Panel_validation_inSilico/RNAseq/Amon_etal_2018/DATA/metadata.csv"
  ),
  Calvanese_2019 = list(
    expr = "4a_Panel_validation_inSilico/RNAseq/Calvanese_2019/DATA/counts_Final.csv",
    meta = "4a_Panel_validation_inSilico/RNAseq/Calvanese_2019/DATA/metadata.csv"
  ), 
  Fresh_InHouse = list(
    expr = "4a_Panel_validation_inSilico/inHouse/1_fresh_set/counts_Final_RUVg_p07k5 (2).csv",
    meta = "4a_Panel_validation_inSilico/inHouse/1_fresh_set/metadata.csv"
  ), 
  Expanded_InHouse = list(
    expr = "4a_Panel_validation_inSilico/inHouse/2_expanded_set/counts_Final_without_100_RUVg_p06k5.csv",
    meta = "4a_Panel_validation_inSilico/inHouse/2_expanded_set/metadata_without_100.csv"
  )
)



# Paths for all microarray studies (RMA-normalized + metadata).
study_paths_Array <- list(
  Rundberg_2015 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/Rundberg_2015/DATA/normalized_expression_RMA.csv",
    meta = "4a_Panel_validation_inSilico/Microarray/Rundberg_2015/DATA/metadata.csv"
  ),
  Rapin_2012 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/Rapin_2014/DATA/normalized_expression_RMA.csv",
    meta = "4a_Panel_validation_inSilico/Microarray/Rapin_2014/DATA/metadata.csv"
  ),
  Prashad_2014 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/Prashad_2014/DATA/normalized_expression_RMA.csv",
    meta = "4a_Panel_validation_inSilico/Microarray/Prashad_2014/DATA/metadata.csv"
  ),
  Nooij_2006 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/Nooij_2006/DATA/normalized_expression_RMA.csv",
    meta = "4a_Panel_validation_inSilico/Microarray/Nooij_2006/DATA/metadata.csv"
  ),
  Barreyro_2012 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/Barreyro_2012/DATA/normalized_expression_RMA.csv",
    meta = "4a_Panel_validation_inSilico/Microarray/Barreyro_2012/DATA/metadata.csv"
  ),
  Gentles_2010 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/Gentles_2010/DATA/normalized_expression_RMA.csv",
    meta = "4a_Panel_validation_inSilico/Microarray/Gentles_2010/DATA/metadata.csv"
  ),
  Souyri_2019 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/Souyri_2019/DATA/normalized_expression_RMA.csv",
    meta = "4a_Panel_validation_inSilico/Microarray/Souyri_2019/DATA/metadata.csv"
  ),
  NI_2005 = list(
    expr = "4a_Panel_validation_inSilico/Microarray/NI_2005/DATA/normalized_expression_RMA.csv",
    meta = "4a_Panel_validation_inSilico/Microarray/NI_2005/DATA/metadata.csv"
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

##### Impute missing genes within a platform (RNA-seq or microarray) #####

impute_by_platform <- function(study_list, presence_threshold = 0.5) {
  # presence_threshold: proporção mínima de estudos em que o gene deve aparecer
  #                     valor entre 0 e 1 (ex.: 0.5 = pelo menos 50% dos estudos)
  
  if (presence_threshold <= 0 || presence_threshold > 1) {
    stop("presence_threshold must be > 0 and ≤ 1 (e.g., 0.5 for 50% of studies).")
  }
  
  # 1) Collect gene sets (column names of normalized_result) for each study
  gene_sets <- lapply(study_list, function(s) colnames(s$normalized_result))
  all_genes <- sort(unique(unlist(gene_sets)))
  n_studies <- length(study_list)
  
  # 2) Count in how many studies each gene appears
  gene_counts <- sapply(all_genes, function(g) {
    sum(vapply(gene_sets, function(gs) g %in% gs, logical(1)))
  })
  
  # 3) Define minimum presence (número mínimo de estudos)
  min_presence <- ceiling(presence_threshold * n_studies)
  
  # Garante pelo menos 1 estudo
  min_presence <- max(1L, min_presence)
  
  # Select genes that appear in at least 'min_presence' studies
  keep_genes <- names(gene_counts)[gene_counts >= min_presence]
  
  message(
    "Keeping ", length(keep_genes),
    " genes present in ≥ ", min_presence,
    " of ", n_studies, " studies (threshold = ",
    presence_threshold, ")."
  )
  
  # 4) Compute a global median expression per gene (across all studies where the gene exists)
  gene_medians <- setNames(numeric(length(keep_genes)), keep_genes)
  
  for (g in keep_genes) {
    values <- c()
    for (study_name in names(study_list)) {
      mat <- study_list[[study_name]]$normalized_result
      if (g %in% colnames(mat)) {
        values <- c(values, mat[, g])
      }
    }
    gene_medians[g] <- median(values, na.rm = TRUE)
  }
  
  # 5) For each study, create an imputed matrix based on normalized_result
  for (study_name in names(study_list)) {
    mat <- study_list[[study_name]]$normalized_result  # samples x genes
    
    # Genes we WANT to have (keep_genes) but are missing in this study
    missing_genes <- setdiff(keep_genes, colnames(mat))
    
    if (length(missing_genes) > 0) {
      # Build a matrix where each missing gene has a constant value = gene_medians[g]
      add_mat <- matrix(
        gene_medians[missing_genes],
        nrow = nrow(mat),
        ncol = length(missing_genes),
        byrow = TRUE
      )
      colnames(add_mat) <- missing_genes
      rownames(add_mat) <- rownames(mat)
      
      # Bind original + imputed genes
      mat_extended <- cbind(mat, add_mat)
    } else {
      mat_extended <- mat
    }
    
    # (Optional) sort columns for consistency
    mat_extended <- mat_extended[, sort(colnames(mat_extended)), drop = FALSE]
    
    # Store the imputed matrix inside the study object
    study_list[[study_name]]$imputed <- mat_extended
  }
  
  return(study_list)
}


##### Run platform-specific imputation #####

# RNA-seq platform
studies_RNAseq <- impute_by_platform(studies_RNAseq, presence_threshold = 0.5)

# Microarray platform
studies_Array  <- impute_by_platform(studies_Array,  presence_threshold = 0.5)

# Rebuild the combined list if you still want `studies` geral
studies <- c(studies_RNAseq, studies_Array)


#### 2. Extract biological trait vector ####
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

#### 3. Correlation between different gene sets — common base ####

# 3.1 Global intersection (from raw expData; pre-normalization)
gene_lists_imputed <- lapply(studies, function(s) colnames(s$imputed))
common_genes_imputed <- Reduce(intersect, gene_lists_imputed)

message("Global intersection size after imputation = ", length(common_genes_imputed))

# 3.2 For each study, keep only the common genes in the IMPUTED matrix
#     and store as 'normalized_filtered' (samples x common_genes_imputed)

for (study_name in names(studies)) {
  # Use the imputed matrix (samples x genes)
  X <- studies[[study_name]]$imputed
  
  # Safety check: all 'common_genes_imputed' must exist in the matrix columns
  missing <- setdiff(common_genes_imputed, colnames(X))
  if (length(missing) > 0) {
    stop(paste0("Missing ", length(missing), " common genes (imputed universe) in study ",
                study_name, ". Check gene ID mapping/column names or imputation step."))
  }
  
  # Reorder columns to a consistent order across all studies
  X_filtered <- X[, common_genes_imputed, drop = FALSE]
  
  # Store inside 'studies'
  studies[[study_name]]$normalized_filtered <- X_filtered
  studies[[study_name]]$n_common_genes     <- ncol(X_filtered)
  studies[[study_name]]$n_samples          <- nrow(X_filtered)
}

# 3.3 Quick sanity check: all studies must have identical #genes (intersection)
print(table(sapply(studies, function(s) s$n_common_genes)))

# 3.4 Define the unique universe for null draws (used downstream)
gene_universe <- common_genes_imputed


#### 4. Filtering the lists of biological interest (stored inside 'studies') ####

# 4.1 Load the gene table (contains List | Ensembl | Gene)
gene_table <- read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Overlap_Uncult_Cult/Gene_Prioritization/priority_lists_long.csv")

# 4.2 Select only the target list (list3_4or5of5)
selected_list_name <- "list3_4or5of5"
selected_genes <- unique(as.character(gene_table$Ensembl[gene_table$List == selected_list_name]))

# Quick sanity check
message("Number of Ensembl IDs in selected list: ", length(selected_genes))

# 4.3 For each study, filter IMPUTED data to keep only those genes
for (study_name in names(studies)) {
  # agora usamos a matriz imputada (samples x genes)
  expr_matrix <- studies[[study_name]]$imputed
  
  # Intersect the selected genes with the available columns
  genes_in_study <- intersect(selected_genes, colnames(expr_matrix))
  
  # Subset matrix
  filtered_matrix <- expr_matrix[, genes_in_study, drop = FALSE]
  
  # Store inside the study object
  studies[[study_name]]$filtered_interest <- filtered_matrix
  
  # Optional: message about how many genes were retained
  message(
    sprintf("[%s] retained %d of %d genes from the selected list",
            study_name, length(genes_in_study), length(selected_genes))
  )
}

# 4.4 Build a summary table showing how many genes were found in each study
retained_counts <- sapply(studies, function(s) ncol(s$filtered_interest))
retained_df <- data.frame(
  Study   = names(retained_counts),
  N_genes = retained_counts
)

#### 5. Calculate the null distribution ####

shared_sets <- build_shared_random_sets(
  gene_universe = gene_universe,  # aqui usamos o universo já definido pós-imputação
  sizes         = c(70,300),
  n_reps        = 1000,
  seed          = 42
)

# Run across all studies using the same gene lists per iteration
null_results_all <- list()

for (study_name in names(studies)) {
  expr_base    <- studies[[study_name]]$imputed
  trait_vector <- studies[[study_name]]$samp
  
  if (is.null(expr_base) || is.null(trait_vector)) {
    message(sprintf("⚠️ Skipping %s (missing normalized_filtered or samp)", study_name))
    next
  }
  
  message(sprintf("▶ Running null distribution (shared sets) for %s...", study_name))
  
  null_df <- tryCatch({
    run_null_distribution_with_shared_sets(
      study_name     = study_name,
      expr_base      = expr_base,
      trait_vector   = trait_vector,
      shared_sets    = shared_sets,
      n_permutations = NULL,      # adaptive inside your perm function
      alternative    = "two.sided",
      return_genes   = TRUE
    )
  }, error = function(e) {
    message(sprintf("❌ Error in %s: %s", study_name, e$message))
    NULL
  })
  
  if (!is.null(null_df)) {
    null_results_all[[study_name]] <- null_df
  }
}

# Combine and save
null_results_combined <- dplyr::bind_rows(null_results_all, .id = "Study")

#### 6. Calculate statistics ####
# Extract the set size from "Subset" (e.g., "Random_5" -> 5)
df <- null_results_combined %>%
  mutate(
    Size   = as.integer(gsub("Random_", "", Subset)),     # Numeric set size
    Size_f = factor(Size, levels = sort(unique(Size)))    # Factor for plotting order
  )

# Remove rows where correlation is missing (safety filter)
df <- df %>% filter(is.finite(Rho))

# 5.1 Individual-level statistics
stats_individual <- df %>%
  group_by(Size) %>%
  summarise(
    n_points   = n(),                             # Total number of individual correlations
    p_gt_0_9   = mean(Rho > 0.9, na.rm = TRUE),   # Probability of ρ > 0.9
    mean_rho   = mean(Rho, na.rm = TRUE),         # Mean correlation
    median_rho = median(Rho, na.rm = TRUE),       # Median correlation
    iqr_rho    = IQR(Rho, na.rm = TRUE),          # Interquartile range
    .groups = "drop"
  )

stats_individual

# 5.2 Median correlation per set
med_by_set <- df %>%
  group_by(Size, Rep) %>%
  summarise(median_rho = median(Rho, na.rm = TRUE), .groups = "drop") %>%
  mutate(Size_f = factor(Size, levels = sort(unique(Size))))

# Compute summary statistics for the median distributions
stats_median <- med_by_set %>%
  group_by(Size) %>%
  summarise(
    n_sets     = n(),                                   # Number of random sets (repetitions)
    p_gt_0_87  = mean(median_rho > 0.87, na.rm = TRUE), # Probability of median ρ > 0.87
    p_gt_0_9   = mean(median_rho > 0.9,  na.rm = TRUE), # Probability of median ρ > 0.9
    mean_rho   = mean(median_rho, na.rm = TRUE),        # Mean of medians
    median_rho = median(median_rho, na.rm = TRUE),      # Median of medians
    iqr_rho    = IQR(median_rho, na.rm = TRUE),         # IQR of medians
    .groups = "drop"
  )

stats_median

#### 7. Create the plots ####

## Paleta fixa: 70 = vermelho, 300 = azul
fill_cols <- c("70" = "red", "300" = "blue")

# 7.3 Helper: caption com probabilidades

make_caption <- function(stats_tbl, size_col = "Size", p_col = "p_gt_0_9",
                         cutoff = 0.9, sizes_order = c(70, 300), digits = 3) {
  sizes <- sizes_order[sizes_order %in% stats_tbl[[size_col]]]
  probs <- stats_tbl[[p_col]][match(sizes, stats_tbl[[size_col]])]
  parts <- paste0(sizes, " Genes: ", formatC(probs, format = "f", digits = digits))
  paste0("p (Correlation > ", cutoff, "):\n", paste(parts, collapse = " | "))
}

cap_individual <- make_caption(stats_individual, cutoff = 0.9, sizes_order = c(70, 300))
cap_median     <- make_caption(stats_median,     cutoff = 0.9, sizes_order = c(70, 300))

# 7.3 Histogram of individual ρ with caption
p_individual_legend <- ggplot(df, aes(x = Rho, fill = Size_f)) +
  geom_histogram(
    bins = 80,
    alpha = 0.40,
    position = "identity",
    color = "black",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 0.9,
    linetype = "dashed",
    color = "red"
  ) +
  scale_fill_manual(
    values = fill_cols,
    name   = "Random Set (genes)"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title   = "Null Distribution",
    x       = "Spearman Correlation (PC1 vs LT-HSC Enrichment Rank)",
    y       = "Frequency",
    caption = cap_individual
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x     = element_text(size = 14),
    axis.title.y     = element_text(size = 14),
    axis.text        = element_text(size = 13),
    legend.position  = "right",
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    plot.caption     = element_text(size = 12, hjust = 0.5, face = "bold", color = "black"),
    plot.caption.position = "plot"
  )
p_individual_legend

ggsave(
  "Null_Distribution/hist_individual_overlay.png",
  p_individual_legend,
  width = 9,
  height = 5,
  dpi = 300
)

# 7.4 Histogram of median ρ with caption
p_median_legend <- ggplot(med_by_set, aes(x = median_rho, fill = Size_f)) +
  geom_histogram(
    bins = 60,
    alpha = 0.50,
    position = "identity",
    color = "black",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 0.9,
    linetype = "dashed",
    color = "red"
  ) +
  scale_fill_manual(
    values = fill_cols,
    name   = "Random Set (genes)"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title   = "Null Distribution — Median ρ across studies",
    x       = "Spearman Correlation (PC1 vs LT-HSC Enrichment Rank)",
    y       = "Frequency",
    caption = cap_median
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x     = element_text(size = 10, color = "black"),
    axis.title.y     = element_text(size = 10, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    legend.position  = "right",
    legend.title     = element_text(size = 10, color = "black"),
    legend.text      = element_text(size = 9, color = "black"),
    plot.caption     = element_text(size = 12, hjust = 0.5, face = "bold", color = "black"),
    plot.caption.position = "plot"
  )
p_median_legend

ggsave(
  "Artigo_GBM_VPO/Suplementar/hist_median_overlay.png",
  p_median_legend,
  width = 7,
  height = 5,
  dpi = 300
)

#Shuffed

shuffle_trait <- function(trait_vector, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # embaralha APENAS os valores, mantendo os nomes das amostras
  shuffled_values <- sample(unname(trait_vector))
  shuffled_trait  <- setNames(shuffled_values, names(trait_vector))
  
  return(shuffled_trait)
}


#### 5b. Calculate the null distribution with SHUFFLED trait ####

null_results_all_shuf <- list()

for (study_name in names(studies)) {
  expr_base    <- studies[[study_name]]$imputed
  trait_vector <- studies[[study_name]]$samp
  
  if (is.null(expr_base) || is.null(trait_vector)) {
    message(sprintf("⚠️ Skipping %s (missing imputed or samp)", study_name))
    next
  }
  
  # embaralha o fenótipo UMA vez por estudo
  trait_shuf <- shuffle_trait(trait_vector)
  
  message(sprintf("▶ Running SHUFFLED null distribution (shared sets) for %s...", study_name))
  
  null_df_shuf <- tryCatch({
    run_null_distribution_with_shared_sets(
      study_name     = study_name,
      expr_base      = expr_base,
      trait_vector   = trait_shuf,   # << AQUI entra o embaralhado
      shared_sets    = shared_sets,
      n_permutations = NULL,
      alternative    = "two.sided",
      return_genes   = TRUE
    )
  }, error = function(e) {
    message(sprintf("❌ Error in shuffled %s: %s", study_name, e$message))
    NULL
  })
  
  if (!is.null(null_df_shuf)) {
    null_results_all_shuf[[study_name]] <- null_df_shuf
  }
}

# Combina tudo da versão embaralhada
null_results_combined_shuf <- dplyr::bind_rows(null_results_all_shuf, .id = "Study")

#### 6b. Calculate statistics — SHUFFLED ####

df_shuf <- null_results_combined_shuf %>%
  mutate(
    Size   = as.integer(gsub("Random_", "", Subset)),     # Numeric set size
    Size_f = factor(Size, levels = sort(unique(Size)))    # Factor for plotting order
  ) %>%
  filter(is.finite(Rho))

# 6b.1 Individual-level statistics (shuffled)
stats_individual_shuf <- df_shuf %>%
  group_by(Size) %>%
  summarise(
    n_points   = n(),                             # Total number of individual correlations
    p_gt_0_9   = mean(Rho > 0.9, na.rm = TRUE),   # Probability of ρ > 0.9
    mean_rho   = mean(Rho, na.rm = TRUE),         # Mean correlation
    median_rho = median(Rho, na.rm = TRUE),       # Median correlation
    iqr_rho    = IQR(Rho, na.rm = TRUE),          # Interquartile range
    .groups = "drop"
  )

# 6b.2 Median correlation per set (shuffled)
med_by_set_shuf <- df_shuf %>%
  group_by(Size, Rep) %>%
  summarise(median_rho = median(Rho, na.rm = TRUE), .groups = "drop") %>%
  mutate(Size_f = factor(Size, levels = sort(unique(Size))))

stats_median_shuf <- med_by_set_shuf %>%
  group_by(Size) %>%
  summarise(
    n_sets     = n(),                                   # Number of random sets (repetitions)
    p_gt_0_87  = mean(median_rho > 0.87, na.rm = TRUE), # Probability of median ρ > 0.87
    p_gt_0_9   = mean(median_rho > 0.9,  na.rm = TRUE), # Probability of median ρ > 0.9
    mean_rho   = mean(median_rho, na.rm = TRUE),        # Mean of medians
    median_rho = median(median_rho, na.rm = TRUE),      # Median of medians
    iqr_rho    = IQR(median_rho, na.rm = TRUE),         # IQR of medians
    .groups = "drop"
  )

#### 7b. Create the plots — SHUFFLED ####

## Paleta fixa: 70 = vermelho, 300 = azul (mesma dos originais)
fill_cols <- c("70" = "red", "300" = "blue")

# Reaproveita a mesma função de caption:
make_caption <- function(stats_tbl, size_col = "Size", p_col = "p_gt_0_9",
                         cutoff = 0.9, sizes_order = c(70, 300), digits = 3) {
  sizes <- sizes_order[sizes_order %in% stats_tbl[[size_col]]]
  probs <- stats_tbl[[p_col]][match(sizes, stats_tbl[[size_col]])]
  parts <- paste0(sizes, " Genes: ", formatC(probs, format = "f", digits = digits))
  paste0("p (Correlation > ", cutoff, "):\n", paste(parts, collapse = " | "))
}

cap_individual_shuf <- make_caption(stats_individual_shuf, cutoff = 0.9, sizes_order = c(70, 300))
cap_median_shuf     <- make_caption(stats_median_shuf,     cutoff = 0.9, sizes_order = c(70, 300))

# 7b.1 Histogram of individual ρ with caption — SHUFFLED
p_individual_legend_shuf <- ggplot(df_shuf, aes(x = Rho, fill = Size_f)) +
  geom_histogram(
    bins = 80,
    alpha = 0.40,
    position = "identity",
    color = "black",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 0.9,
    linetype = "dashed",
    color = "red"
  ) +
  scale_fill_manual(
    values = fill_cols,
    name   = "Random Set (genes)"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title   = "Null Distribution (shuffled trait)",
    x       = "Spearman Correlation (PC1 vs LT-HSC Enrichment Rank)",
    y       = "Frequency",
    caption = cap_individual_shuf
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x     = element_text(size = 14),
    axis.title.y     = element_text(size = 14),
    axis.text        = element_text(size = 13),
    legend.position  = "right",
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    plot.caption     = element_text(size = 12, hjust = 0.5, face = "bold", color = "black"),
    plot.caption.position = "plot"
  )

p_individual_legend_shuf

ggsave(
  "Null_Distribution/hist_individual_overlay_shuffled.png",
  p_individual_legend_shuf,
  width = 9,
  height = 5,
  dpi = 300
)

# 7b.2 Histogram of median ρ with caption — SHUFFLED
p_median_legend_shuf <- ggplot(med_by_set_shuf, aes(x = median_rho, fill = Size_f)) +
  geom_histogram(
    bins = 60,
    alpha = 0.50,
    position = "identity",
    color = "black",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 0.9,
    linetype = "dashed",
    color = "red"
  ) +
  scale_fill_manual(
    values = fill_cols,
    name   = "Random Set (genes)"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title   = "Null Distribution — Median ρ (shuffled trait)",
    x       = "Spearman Correlation (PC1 vs LT-HSC Enrichment Rank)",
    y       = "Frequency",
    caption = cap_median_shuf
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x     = element_text(size = 10, color = "black"),
    axis.title.y     = element_text(size = 10, color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    legend.position  = "right",
    legend.title     = element_text(size = 10, color = "black"),
    legend.text      = element_text(size = 9,  color = "black"),
    plot.caption     = element_text(size = 12, hjust = 0.5, face = "bold", color = "black"),
    plot.caption.position = "plot"
  )

p_median_legend_shuf

ggsave(
  "Artigo_GBM_VPO/Suplementar/hist_median_overlay_shuffled.png",
  p_median_legend_shuf,
  width = 7,
  height = 5,
  dpi = 300
)

library(patchwork)

# Painel lado a lado (2 colunas)
combined_median_panel <- p_median_legend + p_median_legend_shuf +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Median Null Distribution — Observed vs Shuffled",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# Visualizar no R
combined_median_panel

# Salvar em alta resolução
ggsave(
  "Artigo_GBM_VPO/Suplementar/hist_median_overlay_combined_observed_vs_shuffled.png",
  combined_median_panel,
  width = 12,
  height = 5,
  dpi = 300
)
