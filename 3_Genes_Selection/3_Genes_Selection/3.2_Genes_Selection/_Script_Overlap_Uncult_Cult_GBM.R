#####Overlap Analysis - Type 1 and 2 studies####

#0. Load the data 
#1. Overlap Analysis Between Coexpression Modules

#####0. Load the data #####

# Load the libraries
# --- Core WGCNA & clustering ---
library(WGCNA)          # Weighted Gene Co-expression Network Analysis
library(dynamicTreeCut) # Dynamic tree cut for module detection
library(fastcluster)    # Fast hierarchical clustering
library(FactoMineR)     # PCA/CA/HCPC (correct spelling, not 'factorMire')

# --- RNA-seq & stats ---
library(DESeq2)         # RNA-seq normalization & differential analysis
library(multtest)       # Multiple testing procedures
library(coin)           # Permutation-based statistical tests

# --- Annotation / Bioconductor base ---
library(biomaRt)        # Ensembl/BioMart annotation access
library(BiocGenerics)   # Bioconductor S4 generics (compatibility layer)

# --- Data manipulation & utilities ---
library(dplyr)          # Data wrangling (tidyverse)
library(gtools)         # Miscellaneous R utilities

# --- Parallelization ---
library(parallel)       # Base parallel computing
library(foreach)        # Parallel for-loops
library(iterators)      # Iterators for foreach
library(doParallel)     # Parallel backend for foreach

# --- Visualization ---
library(ggplot2)        # Grammar of Graphics plotting
library(cowplot)        # Plot composition helpers
library(pheatmap)       # Pretty heatmaps
library(gplots)         # heatmap.2 and other plotting utils
library(RColorBrewer)   # Color palettes
library(VennDiagram)    # Venn diagrams

#Set the directory
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC")

#Load the functions

# Function to create and filter a module frequency table based on adjusted p-values
buildFilteredModules <- function(study, studyName) {
  # Guard clauses
  if (is.null(study$correlation) || nrow(study$correlation) == 0) {
    return(data.frame(moduleColor = character(),
                      Freq = integer(),
                      CorrelationDirection = character(),
                      stringsAsFactors = FALSE))
  }
  if (is.null(study$kme) || nrow(study$kme) == 0) {
    stop("kME table is missing or empty; cannot compute Freq per module.")
  }
  
  corr <- study$correlation
  
  # Normalize module column name to 'moduleColor'
  if ("Module" %in% names(corr) && !("moduleColor" %in% names(corr))) {
    corr$moduleColor <- corr$Module
  }
  
  # Required columns check
  required_corr_cols <- c("moduleColor", "Adjusted_P_value", "Correlation")
  missing_corr <- setdiff(required_corr_cols, names(corr))
  if (length(missing_corr) > 0) {
    stop(paste("Missing columns in correlation table:", paste(missing_corr, collapse = ", ")))
  }
  
  # Significant only (Adjusted_P_value < 0.05)
  corr_sig <- subset(
    corr,
    round(Adjusted_P_value, 2) <= 0.05 & !is.na(Correlation)
  )
  
  if (nrow(corr_sig) == 0) {
    return(data.frame(moduleColor = character(),
                      Freq = integer(),
                      CorrelationDirection = character(),
                      stringsAsFactors = FALSE))
  }
  
  # Direction
  corr_sig$CorrelationDirection <- ifelse(corr_sig$Correlation > 0, "Positive", "Negative")
  
  # kME normalization to 'moduleColor'
  kme <- study$kme
  if (!("moduleColor" %in% names(kme)) && ("Module" %in% names(kme))) {
    kme$moduleColor <- kme$Module
  }
  if (!("moduleColor" %in% names(kme))) {
    stop("Could not find 'moduleColor' (or 'Module') column in kME table.")
  }
  
  # Gene count per module (Freq)
  freq_df <- as.data.frame(table(kme$moduleColor), stringsAsFactors = FALSE)
  names(freq_df) <- c("moduleColor", "Freq")
  
  # Merge significant modules with Freq and direction
  out <- merge(
    corr_sig[, c("moduleColor", "CorrelationDirection")],
    freq_df,
    by = "moduleColor",
    all.x = TRUE
  )
  
  # Append study name to moduleColor here (after 'out' exists)
  out$moduleColor <- paste0(out$moduleColor, "_", studyName)
  
  # Final order/cleanup
  out <- out[, c("moduleColor", "Freq", "CorrelationDirection")]
  out$Freq[is.na(out$Freq)] <- 0L
  out <- out[order(-out$Freq, out$moduleColor), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Function to process overlap data for any study
processOverlapData <- function(testResults, modulesStudy, studyName, modulesAll) {
  
  # Extract the enrichment results
  overlap_data <- testResults$pValues
  
  # Remove the suffix from moduleColor to get pure color names
  significantModules <- gsub(paste0("_", studyName), "", modulesStudy$moduleColor)
  
  # Filter: keep only significant input modules from this study
  overlap_data <- overlap_data[
    overlap_data$InputCategories %in% significantModules, ]
  
  # Clean UserDefinedCategories (remove suffixes like "_genes" etc.)
  overlap_data$UserDefinedCategories <- sub("(__.*)$", "", overlap_data$UserDefinedCategories)
  
  # Filter: keep only significant user-defined modules from all studies
  overlap_data <- overlap_data[
    overlap_data$UserDefinedCategories %in% modulesAll$moduleColor, ]
  
  # Restore suffix in InputCategories to identify the study
  overlap_data$InputCategories <- paste0(overlap_data$InputCategories, "_", studyName)
  
  return(overlap_data)
}

# Function to extract study name (suffix) from module name
extract_study <- function(moduleName) {
  sub(".*_(.*)$", "\\1", moduleName)
}

loadStudyData <- function(base_path, study_name) {
  # 0) Resolve study_dir: allow passing either the study folder itself
  #    or a parent folder that contains the study folder.
  if (dir.exists(file.path(base_path, study_name))) {
    study_dir <- file.path(base_path, study_name)
  } else {
    study_dir <- base_path
  }
  message(">> Using study_dir = ", normalizePath(study_dir, winslash = "/"))
  
  # 1) Build canonical file paths
  kme_path  <- file.path(study_dir, "RESULTS", "WGCNA", "3_kmeTable.csv")
  corr_path <- file.path(study_dir, "RESULTS", "WGCNA", "5_Correlation_Values.csv")
  
  # counts can vary by study; try Final first, then NonNormalized as fallback
  counts_path_final <- file.path(study_dir, "DATA", "counts_Final.csv")
  counts_path_raw   <- file.path(study_dir, "DATA", "counts_NonNormalized.csv")
  counts_path <- if (file.exists(counts_path_final)) counts_path_final else counts_path_raw
  
  metadata_path <- file.path(study_dir, "DATA", "metadata.csv")
  
  # 2) Sanity checks (fail fast with informative messages)
  missing_files <- c()
  for (p in list("KME" = kme_path, "Correlation" = corr_path, "Counts" = counts_path, "Metadata" = metadata_path)) {
    if (!file.exists(p)) missing_files <- c(missing_files, p)
  }
  if (length(missing_files) > 0) {
    stop("Missing file(s) for study '", study_name, "':\n  - ",
         paste(missing_files, collapse = "\n  - "), call. = FALSE)
  }
  
  # 3) Read expression counts (genes in rows; samples in columns)
  #    We keep row.names = 1 to treat the first column as gene IDs (Ensembl).
  expData <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  message("   Loaded counts: ", basename(counts_path), 
          " [genes x samples = ", nrow(expData), " x ", ncol(expData), "]")
  
  # 4) Read metadata robustly:
  #    Try (A) with row.names=1 and (B) without row.names, then choose the variant
  #    whose rownames best match the expression colnames.
  metaA <- try(suppressWarnings(read.csv(metadata_path, row.names = 1, check.names = FALSE)), silent = TRUE)
  metaB <- try(suppressWarnings(read.csv(metadata_path, check.names = FALSE)), silent = TRUE)
  
  choose_meta <- function(metaA, metaB, exp_cols) {
    score <- function(df) {
      if (inherits(df, "try-error")) return(-1L)
      rn <- rownames(df)
      # Score = how many rownames match exp colnames
      if (is.null(rn)) return(-1L)
      sum(rn %in% exp_cols)
    }
    sA <- score(metaA); sB <- score(metaB)
    if (sA >= sB && sA >= 0) {
      attr(metaA, "source") <- "row.names=1"
      return(metaA)
    } else if (sB >= 0) {
      # If metaB wins, try to set rownames from first column when it looks like sample IDs
      df <- metaB
      first_col <- df[[1]]
      if (length(first_col) == nrow(df) && all(first_col %in% exp_cols)) {
        rownames(df) <- as.character(first_col)
        df <- df[, -1, drop = FALSE]
      }
      attr(df, "source") <- "no row.names (adjusted)"
      return(df)
    } else {
      return(metaA) # both failed; will be caught later
    }
  }
  
  metadata <- choose_meta(metaA, metaB, colnames(expData))
  if (inherits(metadata, "try-error") || is.null(rownames(metadata))) {
    stop("Failed to parse metadata rownames for '", study_name, 
         "'. Please check the file format: ", metadata_path, call. = FALSE)
  }
  message("   Loaded metadata (", attr(metadata, "source"), 
          "). Samples in metadata: ", nrow(metadata))
  
  # 5) Align samples: ensure counts columns and metadata rows match and are ordered
  common_samples <- intersect(colnames(expData), rownames(metadata))
  if (length(common_samples) == 0) {
    stop("No overlapping sample IDs between counts and metadata for '", study_name, "'.", call. = FALSE)
  }
  # Reorder both to the same sample order
  expData   <- expData[, common_samples, drop = FALSE]
  metadata  <- metadata[common_samples, , drop = FALSE]
  message("   Matched ", length(common_samples), " samples between counts and metadata.")
  
  # 6) Read KME and correlation tables
  kme <- read.csv(kme_path, check.names = FALSE)
  correlation <- read.csv(corr_path, check.names = FALSE)
  
  # 7) Deduplicate genes by Ensembl ID (keep first occurrence)
  kme <- dplyr::distinct(kme, Ensembl, .keep_all = TRUE)
  
  # 8) Return structured list
  return(list(
    kme = kme,
    exp = expData,
    metadata = metadata,
    correlation = correlation
  ))
}

getModulesAll <- function(study_data, reference_study, include_reference = FALSE) {
  # Define which studies to include in the module comparison
  studies_to_use <- if (include_reference) {
    names(study_data)
  } else {
    setdiff(names(study_data), reference_study)
  }
  
  # Log which studies are being used
  message(paste("Modules combined from:", paste(studies_to_use, collapse = ", ")))
  
  # Bind all significant modules from the selected studies
  modules_all <- do.call(rbind, lapply(studies_to_use, function(study) {
    study_data[[study]]$modules
  }))
  
  return(modules_all)
}

runModuleOverlap <- function(study_name, kme_data, input_file, output_file) {
  geneR <- kme_data$Ensembl
  labelR <- kme_data$moduleColor
  
  result <- userListEnrichment(
    geneR, labelR,
    fnIn = input_file,
    nameOut = output_file,
    omitCategories = c("M0", "grey"),
    outputGenes = TRUE
  )
  
  return(result)
}

.get_study_full_from_module <- function(module_id) {
  sub("^[^_]+_", "", module_id)  # strip the leading "color_"
}

.abbrev_study <- function(study_full) {
  # Take part before the first underscore, then before hyphen
  s <- sub("_.*$", "", study_full)
  s <- sub("-.*$", "", s)
  s
}

# Main: process overlap results for labeledHeatmap with short study tags
process_overlap_for_heatmap <- function(
    overlap_results,
    modules_all,
    correlation_sign = "Positive",
    study_order,               # e.g., c("Anjos","Xie","Fares","Sauvageau","XieCult")
    csv_path = NULL,
    study_abbrev_map = NULL,   # optional named vector: c("Anjos-Afonso_etal_2021"="Anjos", ...)
    show_counts = FALSE        # if TRUE, append "\nN" under the short label (module frequency)
) {
  # 1) Combine all overlap results into a single data.frame
  overlap_data_combined <- do.call(rbind, overlap_results)
  
  # 2) Keep only modules with the desired correlation sign
  target_modules <- modules_all$moduleColor[modules_all$CorrelationDirection == correlation_sign]
  
  overlap_filtered <- overlap_data_combined[
    overlap_data_combined$InputCategories %in% target_modules &
      overlap_data_combined$UserDefinedCategories %in% target_modules, ]
  
  # 3) Unique module IDs (FULL IDs, e.g., "blue_Anjos-Afonso_etal_2021")
  input_modules <- unique(overlap_filtered$InputCategories)
  user_modules  <- unique(overlap_filtered$UserDefinedCategories)
  
  # 4) Order modules by study order USING STUDY ABBREVIATIONS
  #    - First, get full study names from module IDs
  in_study_full  <- vapply(input_modules, .get_study_full_from_module, character(1))
  usr_study_full <- vapply(user_modules,  .get_study_full_from_module, character(1))
  
  #    - Then get abbreviations (from map if provided, else heuristic)
  map_or_abbrev <- function(full_names) {
    if (is.null(study_abbrev_map)) {
      vapply(full_names, .abbrev_study, character(1))
    } else {
      ifelse(full_names %in% names(study_abbrev_map),
             unname(study_abbrev_map[full_names]),
             vapply(full_names, .abbrev_study, character(1)))
    }
  }
  in_study_abbrev  <- map_or_abbrev(in_study_full)
  usr_study_abbrev <- map_or_abbrev(usr_study_full)
  
  #    - Order by factor using the provided study_order (ties broken by module id)
  ordered_input <- input_modules[order(factor(in_study_abbrev,  levels = study_order),  input_modules)]
  ordered_user  <- user_modules[ order(factor(usr_study_abbrev, levels = study_order),  user_modules)]
  
  # 5) Initialize matrices
  overlapMatrix <- matrix(0, nrow = length(ordered_input), ncol = length(ordered_user),
                          dimnames = list(ordered_input, ordered_user))
  pValueMatrix  <- matrix(1, nrow = length(ordered_input), ncol = length(ordered_user),
                          dimnames = list(ordered_input, ordered_user))
  
  # 6) Fill matrices
  for (i in seq_len(nrow(overlap_filtered))) {
    mod1 <- overlap_filtered$InputCategories[i]
    mod2 <- overlap_filtered$UserDefinedCategories[i]
    overlapMatrix[mod1, mod2] <- overlap_filtered$NumOverlap[i]
    pValueMatrix[mod1, mod2]  <- overlap_filtered$CorrectedPvalues[i]
  }
  
  # 7) -log10(FDR)
  logPValueMatrix <- -log10(pValueMatrix)
  logPValueMatrix[is.infinite(logPValueMatrix) | is.na(logPValueMatrix)] <- 130
  
  # 8) Cell text: "nOverlap\n(FDR)"
  textMatrix <- paste0(overlapMatrix, "\n(", signif(pValueMatrix, 2), ")")
  dim(textMatrix) <- dim(overlapMatrix)
  
  # 9) Build color-only labels (for color bars) and short human labels (for axes)
  #    - Colors: remove "_Study" suffix -> "blue", "brown", ...
  x_colors_only <- sub("_(.*)$", "", colnames(overlapMatrix))
  y_colors_only <- sub("_(.*)$", "", rownames(overlapMatrix))
  
  #    - Full study names from the matrix dimnames
  x_study_full <- vapply(colnames(overlapMatrix), .get_study_full_from_module, character(1))
  y_study_full <- vapply(rownames(overlapMatrix), .get_study_full_from_module, character(1))
  
  #    - Abbreviations (map first, else heuristic)
  if (is.null(study_abbrev_map)) {
    x_abbrev <- vapply(x_study_full, .abbrev_study, character(1))
    y_abbrev <- vapply(y_study_full, .abbrev_study, character(1))
  } else {
    x_abbrev <- ifelse(x_study_full %in% names(study_abbrev_map),
                       unname(study_abbrev_map[x_study_full]),
                       vapply(x_study_full, .abbrev_study, character(1)))
    y_abbrev <- ifelse(y_study_full %in% names(study_abbrev_map),
                       unname(study_abbrev_map[y_study_full]),
                       vapply(y_study_full, .abbrev_study, character(1)))
  }
  
  #    - Short symbols: "color_Abbrev" (e.g., "blue_Anjos")
  x_symbols_short <- paste0(x_colors_only, "_", x_abbrev)
  y_symbols_short <- paste0(y_colors_only, "_", y_abbrev)
  
  # 10) Optional counts underneath (uses modules_all$Freq keyed by FULL module ID)
  freq_lookup <- if ("Freq" %in% names(modules_all)) {
    stats::setNames(modules_all$Freq, modules_all$moduleColor)
  } else {
    stats::setNames(rep(NA_integer_, nrow(modules_all)), modules_all$moduleColor)
  }
  x_freq <- as.integer(freq_lookup[colnames(overlapMatrix)])
  y_freq <- as.integer(freq_lookup[rownames(overlapMatrix)])
  
  if (isTRUE(show_counts)) {
    x_symbols <- ifelse(is.na(x_freq), x_symbols_short, paste0(x_symbols_short, "\n", x_freq))
    y_symbols <- ifelse(is.na(y_freq), y_symbols_short, paste0(y_symbols_short, "\n", y_freq))
  } else {
    x_symbols <- x_symbols_short
    y_symbols <- y_symbols_short
  }
  
  # 11) "ME" labels for WGCNA color bars (keep PURE colors here!)
  x_me_labels <- paste0("Me", x_colors_only)
  y_me_labels <- paste0("Me", y_colors_only)
  
  # 12) Save CSV if requested
  if (!is.null(csv_path)) {
    utils::write.csv(overlap_filtered, file = csv_path, row.names = FALSE)
  }
  
  # 13) Return bundle
  return(list(
    logP           = logPValueMatrix,
    textMatrix     = textMatrix,
    # color-only labels (for color bars or internal use)
    xLabels        = x_colors_only,
    yLabels        = y_colors_only,
    # human-facing axis text (short study tag)
    xSymbols       = x_symbols,
    ySymbols       = y_symbols,
    # also return the short symbols without counts, in case you want to switch in plotting code
    xSymbols_short = x_symbols_short,
    ySymbols_short = y_symbols_short,
    # WGCNA-style color labels
    xColors        = x_me_labels,
    yColors        = y_me_labels,
    # matrices for plotting and auditing
    overlapMatrix  = overlapMatrix,
    pValueMatrix   = pValueMatrix
  ))
}

plot_overlap_heatmap <- function(processed_data,
                                 output_path,
                                 main_title,
                                 study_labels = list(x = NULL, y = NULL)) {
  
  # Open PNG device
  png(output_path, width = 1200, height = 900, res = 150)
  par(mar = c(8, 10, 3, 3))
  par(mgp = c(5, 2, 0))
  
  # Generate labeled heatmap
  labeledHeatmap(
    Matrix = processed_data$logP,
    xLabels = processed_data$xColors,
    yLabels = processed_data$yColors,
    xSymbols = processed_data$xSymbols,
    ySymbols = processed_data$ySymbols,
    colorLabels = FALSE,
    colors = colorRampPalette(c("white", "red"))(50),
    textMatrix = processed_data$textMatrix,
    setStdMargins = FALSE,
    main = main_title,
    cex.text = 0.8,
    cex.lab.y = 0.8,
    cex.lab.x = 0.8,
    cex.legendLabel = 0.7,
    legendLabel = expression(-log[10]("FDR"))
  )
  
  # Add study labels to the axes
  if (!is.null(study_labels$x)) mtext(study_labels$x, side = 1, line = 6, cex = 1.1, font = 2)
  if (!is.null(study_labels$y)) mtext(study_labels$y, side = 2, line = 6, cex = 1.1, font = 2)
  
  dev.off()
}

build_overlap_table <- function(filtered_kme1, filtered_kme2, study_suffix1, study_suffix2, output_file = NULL) {
  # Fixed base columns
  base_cols <- c("Ensembl", "Gene", "Gene_type")
  
  # Identify relevant module labels from both studies
  module_labels1 <- unique(filtered_kme1[[paste0("moduleLabel_", study_suffix1)]])
  module_labels2 <- unique(filtered_kme2[[paste0("moduleLabel_", study_suffix2)]])
  
  # Construct expected kME and p-value column names
  kme_cols1 <- unlist(lapply(module_labels1, function(m) paste0("kME_", m, "_", study_suffix1)))
  pval_cols1 <- unlist(lapply(module_labels1, function(m) paste0("pvalBH_", m, "_", study_suffix1)))
  
  kme_cols2 <- unlist(lapply(module_labels2, function(m) paste0("kME_", m, "_", study_suffix2)))
  pval_cols2 <- unlist(lapply(module_labels2, function(m) paste0("pvalBH_", m, "_", study_suffix2)))
  
  # Also keep module color, label, and expression columns (if available)
  module_cols1 <- c(paste0("moduleColor_", study_suffix1), paste0("moduleLabel_", study_suffix1))
  module_cols2 <- c(paste0("moduleColor_", study_suffix2), paste0("moduleLabel_", study_suffix2))
  
  expression_col1 <- paste0("Expression_", study_suffix1)
  expression_col2 <- paste0("Expression_", study_suffix2)
  
  # Determine which columns to keep from each dataframe
  keep_cols1 <- c(base_cols, module_cols1, expression_col1, kme_cols1, pval_cols1)
  keep_cols2 <- c(base_cols, module_cols2, expression_col2, kme_cols2, pval_cols2)
  
  df1 <- filtered_kme1[, intersect(keep_cols1, colnames(filtered_kme1))]
  df2 <- filtered_kme2[, intersect(keep_cols2, colnames(filtered_kme2))]
  
  # Remove duplicated columns from the second dataframe
  df2 <- df2[, !colnames(df2) %in% c("Gene", "Gene_type")]
  
  # Merge by Ensembl ID
  overlap_df <- dplyr::inner_join(df1, df2, by = "Ensembl")
  
  # If output_file is provided, save the CSV
  if (!is.null(output_file)) {
    write.csv(overlap_df, output_file, row.names = FALSE)
    message("CSV file saved to: ", output_file)
  }
  
  return(overlap_df)
}

filter_kme_by_direction <- function(kme, modules, direction = c("Positive", "Negative"), study_suffix) {
  direction <- match.arg(direction)
  
  # Remove study suffix from moduleColor in modules table, if present
  modules$moduleColor <- sub("_.*$", "", modules$moduleColor)
  
  # Select modules based on direction
  selected_modules <- modules$moduleColor[modules$CorrelationDirection == direction]
  
  # Filter kME data for the selected modules
  filtered_kme <- kme[kme$moduleColor %in% selected_modules, ]
  
  # Define fixed columns to preserve
  fixed_columns <- c("Ensembl", "Gene", "Gene_type")
  
  # Identify columns to rename (all others)
  rename_columns <- setdiff(colnames(filtered_kme), fixed_columns)
  
  # Apply suffix to rename columns
  colnames(filtered_kme)[colnames(filtered_kme) %in% rename_columns] <- 
    paste0(rename_columns, "_", study_suffix)
  
  return(filtered_kme)
}



##### 1. Overlap Analysis Between Coexpression Modules####
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC")

# 1.1 Define Paths and Study Names

# Define the working directory
studyPaths <- getwd()

# List the studies to be used
study_paths <- c(
  "Anjos-Afonso_etal_2021" = "Type_1_Studies/Anjos-Afonso_etal_2021",
  "Xie_etal_2019"          = "Type_1_Studies/Xie_etal_2019",
  "Fares_2017"             = "Type_2_Studies/Fares_2017",
  "Sauvageau_2014"         = "Type_2_Studies/Sauvageau_2014",
  "Xie_2019"               = "Type_2_Studies/Xie_2019"
)


# Load all necessary study data
study_data <- lapply(names(study_paths), function(name) {
  loadStudyData(base_path = study_paths[[name]], study_name = name)
})
names(study_data) <- names(study_paths)

# Filtering correlation data to keep only rows where Trait == "LT_enrichment"
study_data <- lapply(study_data, function(study) {
  if ("correlation" %in% names(study)) {
    study$correlation <- study$correlation[study$correlation$Trait == "LT_enrichment", ]
  }
  study
})

# Creating a specific dataframe to store the frequency and correlation direction
# of modules associated with LT-HSC enrichment
study_data <- Map(function(st, nm) {
  st$modules <- buildFilteredModules(st, nm)
  st
}, study_data, names(study_data))

# 1.2 Prepare Input File for Overlap Enrichment

# Choose the study to extract module genes from
overlap_study <- c("Anjos-Afonso_etal_2021", "Xie_etal_2019" )

# Checks: ensure all studies exist
missing <- setdiff(overlap_study, names(study_data))
if (length(missing)) stop("Studies not found in study_data: ", paste(missing, collapse = ", "))

genes_for_overlap <- dplyr::bind_rows(lapply(overlap_study, function(study) {
  kme <- study_data[[study]]$kme
  # Normalize columns
  if (!"Ensembl" %in% names(kme)) {
    kme <- tibble::rownames_to_column(kme, var = "Ensembl")
  }
  # Accept 'moduleLabel' as fallback if 'moduleColor' is absent
  if (!"moduleColor" %in% names(kme) && "moduleLabel" %in% names(kme)) {
    kme <- dplyr::rename(kme, moduleColor = moduleLabel)
  }
  if (!all(c("Ensembl","moduleColor") %in% names(kme))) {
    stop("Missing Ensembl/moduleColor in ", study)
  }
  
  dplyr::transmute(
    kme,
    Ensembl,
    modules = paste0(moduleColor, "_", study),
    Study = study
  )
}))

# Export the input file to be used in enrichment analysis
write.csv(genes_for_overlap, file = "Overlap_Uncult_Cult/Overlap/0_genes_for_overlap.csv", row.names = FALSE)

# 1.3 Run Overlap Enrichment Using Each Reference Study

# Define the studies that will serve as the reference for comparison
reference_studies <- c("Fares_2017", "Sauvageau_2014", "Xie_2019")

# Set input file path
input_file <- "Overlap_Uncult_Cult/Overlap/0_genes_for_overlap.csv"

# Initialize list to store processed results from each reference
overlap_results <- list()

# Loop through reference studies and run the enrichment pipeline
for (ref_study in reference_studies) {
  
  message(paste0("\n>>> Running overlap for: ", ref_study))
  
  # 1.3.1 Run userListEnrichment from WGCNA
  testResults <- runModuleOverlap(
    study_name = ref_study,
    kme_data = study_data[[ref_study]]$kme,
    input_file = input_file,
    output_file = paste0("Overlap_Uncult_Cult/Overlap/overlap_", ref_study, ".csv")
  )
  
  # 1.3.2 Retrieve all significant modules (excluding the reference study)
  modules_all <- getModulesAll(
    study_data = study_data,
    reference_study = ref_study,
    include_reference = FALSE
  )
  
  # 1.3.3 Format the test results
  processed <- processOverlapData(
    testResults = testResults,
    modulesStudy = study_data[[ref_study]]$modules,
    studyName = ref_study,
    modulesAll = modules_all
  )
  
  # 1.3.4 Save the result in the list
  overlap_results[[ref_study]] <- processed
}

# 1.4 Format Overlap Results for Heatmap

# Define study order for axis sorting
study_order <- c("Anjos-Afonso_etal_2021", "Xie_etal_2019", "Fares_2017", "Sauvageau_2014", "Xie_2019")

# Merge all module data into a single object
modulesAll <- do.call(rbind, lapply(study_data, function(st) st$modules))

study_abbrev_map <- c(  
  "Anjos-Afonso_etal_2021" = "Anjos",
  "Xie_etal_2019"          = "Xie_Type1",
  "Fares_2017" = "Fares",
  "Sauvageau_2014"          = "Sauvageau",
  "Xie_2019" = "Xie_Type2"
)

# Remove the modules that did not show overlaps with the individual type 1 and type 2 studies
to_remove <- c("salmon_Fares_2017",
               "yellow_Xie_etal_2019",
               "yellow_Sauvageau_2014")

# Filter the dataframe to exclude the modules listed above
modulesAll <- modulesAll[!modulesAll[[1]] %in% to_remove, ]

# Prepare overlap matrices and labels for plotting
processed_positive <- process_overlap_for_heatmap(
  overlap_results = overlap_results,
  modules_all     = modulesAll,
  correlation_sign = "Positive",
  study_order      = c("Anjos", "Xie_Type1", "Fares","Sauvageau","Xie_Type2"),
  csv_path         = "Overlap_Uncult_Cult/Overlap/1_Module_Overlap_Positive.csv",
  study_abbrev_map = study_abbrev_map,
  show_counts      = TRUE
)

# Prepare overlap matrices and labels for plotting
processed_negative <- process_overlap_for_heatmap(
  overlap_results = overlap_results,
  modules_all     = modulesAll,
  correlation_sign = "Negative",
  study_order      = c("Anjos", "Xie_Type1", "Fares","Sauvageau","Xie_Type2"),
  csv_path         = "Overlap_Uncult_Cult/Overlap/1_Module_Overlap_Negative.csv",
  study_abbrev_map = study_abbrev_map,
  show_counts      = TRUE
)

# 1.5 Plot the results 

# Plot the graph 
plot_overlap_heatmap(
  processed_data = processed_positive,
  output_path = "Overlap_Uncult_Cult/Overlap/1_Heatmap_Positive_Correlation.png",
  main_title = "Module Overlap - Positive Correlation",
  study_labels = list(x = "", y = "")
)

# Plot the graph 
plot_overlap_heatmap(
  processed_data = processed_negative,
  output_path = "Overlap_Uncult_Cult/Overlap/1_Heatmap_Negative_Correlation.png",
  main_title = "Module Overlap - Negative Correlation",
  study_labels = list(x = "", y = "")
)

#### 2. Intersection between type 2 studies####

# 2.1.1. Filter positively correlated modules

#Anjos
modulesAnjos_positive <- study_data$`Anjos-Afonso_etal_2021`$modules$moduleColor[study_data$`Anjos-Afonso_etal_2021`$modules$CorrelationDirection == "Positive"]
modulesAnjos_clean <- sub("_Anjos-Afonso_etal_2021", "", modulesAnjos_positive)
modulesAnjos_clean <- setdiff(modulesAnjos_clean, "tan")
study_data$`Anjos-Afonso_etal_2021`$kme$LTmodule_Anjos <- ifelse(study_data$`Anjos-Afonso_etal_2021`$kme$moduleColor %in% c("blue", "turquoise"), 1, 0)
study_data$`Anjos-Afonso_etal_2021`$kme$LTmodule_Anjos_0.8 <- ifelse(
  study_data$`Anjos-Afonso_etal_2021`$kme$moduleColor %in% c("blue", "turquoise") &
    (study_data$`Anjos-Afonso_etal_2021`$kme$kME_M1 > 0.8 | study_data$`Anjos-Afonso_etal_2021`$kme$kME_M2 > 0.8), 1, 0)
LTAnjos = study_data$`Anjos-Afonso_etal_2021`$kme[study_data$`Anjos-Afonso_etal_2021`$kme$moduleColor %in% c("blue", "turquoise"), 1:3]

#Xie
modulesXie_Type1_positive <- study_data$Xie_etal_2019$modules$moduleColor[study_data$Xie_etal_2019$modules$CorrelationDirection == "Positive"]
modulesXie_Type1_clean <- sub("_Xie_etal_2019", "", modulesXie_Type1_positive)
study_data$Xie_etal_2019$kme$LTmodule_XieT1 <- ifelse(study_data$Xie_etal_2019$kme$moduleColor %in% c("blue", "brown", "darkturquoise", "green"), 1, 0)
study_data$Xie_etal_2019$kme$LTmodule_XieT1_0.8 <- ifelse(
  study_data$Xie_etal_2019$kme$moduleColor %in% c("blue", "brown", "darkturquoise", "green") &
    (study_data$Xie_etal_2019$kme$kME_M2 > 0.8 | study_data$Xie_etal_2019$kme$kME_M3 > 0.8 | study_data$Xie_etal_2019$kme$kME_M5 > 0.8 | study_data$Xie_etal_2019$kme$kME_M23 > 0.8), 1, 0)
LTXie_T1 = study_data$Xie_etal_2019$kme[study_data$Xie_etal_2019$kme$moduleColor %in% c("blue", "brown", "darkturquoise", "green"), 1:3]

#Fares
modulesFares_positive <- study_data$Fares_2017$modules$moduleColor[study_data$Fares_2017$modules$CorrelationDirection == "Positive"]
modulesFares_clean <- sub("_Fares_2017", "", modulesFares_positive)
modulesFares_clean <- setdiff(modulesFares_clean, c("lightcyan", "magenta", "salmon"))
study_data$Fares_2017$kme$LTmodule_Fares <- ifelse(study_data$Fares_2017$kme$moduleColor %in% c("brown"), 1, 0)
study_data$Fares_2017$kme$LTmodule_Fares_0.8 <- ifelse(
  study_data$Fares_2017$kme$moduleColor %in% "brown" & 
    study_data$Fares_2017$kme$kME_M3 > 0.8, 
  1, 0)
LTFares = study_data$Fares_2017$kme[study_data$Fares_2017$kme$moduleColor %in% c("brown"), 1:3]

#Sauvageau 
modulesSauvageau_positive <- study_data$Sauvageau_2014$modules$moduleColor[study_data$Sauvageau_2014$modules$CorrelationDirection == "Positive"]
modulesSauvageau_clean <- sub("_Sauvageau_2014", "", modulesSauvageau_positive)
study_data$Sauvageau_2014$kme$LTmodule_Sauvageau = ifelse(study_data$Sauvageau_2014$kme$moduleColor %in% c("blue", "brown", "salmon"), 1, 0)
study_data$Sauvageau_2014$kme$LTmodule_Sauvageau_0.8 <- ifelse(
  study_data$Sauvageau_2014$kme$moduleColor %in% c("blue", "brown", "salmon") &
    (study_data$Sauvageau_2014$kme$kME_M2 > 0.8 | study_data$Sauvageau_2014$kme$kME_M3 > 0.8 | study_data$Sauvageau_2014$kme$kME_M13 > 0.8), 1, 0)
LTSauv = study_data$Sauvageau_2014$kme[study_data$Sauvageau_2014$kme$moduleColor %in% c("blue", "brown", "salmon"), 1:3]

#Xie Culture 
modulesXie_Type2_positive <- study_data$Xie_2019$modules$moduleColor[study_data$Xie_2019$modules$CorrelationDirection == "Positive"]
modulesXie_Type2_clean <- sub("_Xie_2019", "", modulesXie_Type2_positive)
study_data$Xie_2019$kme$LTmodule_XieT2 = ifelse(study_data$Xie_2019$kme$moduleColor %in% c("green", "purple", "greenyellow", "cyan", "grey60"), 1, 0)
study_data$Xie_2019$kme$LTmodule_XieT2_0.8 <- ifelse(
  study_data$Xie_2019$kme$moduleColor %in% c("green", "purple", "greenyellow", "cyan", "grey60") &
    (study_data$Xie_2019$kme$kME_M5 > 0.8 | study_data$Xie_2019$kme$kME_M10 > 0.8 | study_data$Xie_2019$kme$kME_M11 > 0.8 | study_data$Xie_2019$kme$kME_M14 > 0.8 | study_data$Xie_2019$kme$kME_M17 > 0.8), 1, 0)
LTXie_T2 = study_data$Xie_2019$kme[study_data$Xie_2019$kme$moduleColor %in% c("green", "purple", "greenyellow", "cyan", "grey60"), 1:3]

#merge all genes in LT-associated modules and remove duplicates
LTgenes=rbind(LTAnjos,LTXie_T1,LTFares,LTSauv,LTXie_T2)
LTgenes= LTgenes %>% distinct(Ensembl, .keep_all = TRUE)

#create a table holding the module assignment of the gene in each study
#plus the columns with kme info (only for the LT-associated modules in that study)
#and the column that says if the gene is or not in LT-associated module in that study

#subset columns Anjos
redkmeAnj <- study_data$`Anjos-Afonso_etal_2021`$kme[,c(1,4,5, 7:10,41:42)]
colnames(redkmeAnj)[c(2:7)]=c("moduleColor_Anjos","moduleLable_Anjos",
                              "kME_M1_Anjos","pvalBH_M1_Anjos",
                              "kME_M2_Anjos","pvalBH_M2_Anjos")

#subset columns Xie
redkmeXie <- study_data$Xie_etal_2019$kme[,c(1,4,5, 9:12,15,16,51,52,61:62)]
colnames(redkmeXie)[c(2:11)]=c("moduleColor_XieT1","moduleLable_XieT1",
                               "kME_M2_XieT1","pvalBH_M2_XieT1", 
                               "kME_M3_XieT1","pvalBH_M3_XieT1", 
                               "kME_M5_XieT1","pvalBH_M5_XieT1",
                               "kME_M23_XieT1","pvalBH_M23_XieT1")

#subset columns Fares
redkmeFares <- study_data$Fares_2017$kme[,c (1,4,5,11,12,43:44)] 
colnames(redkmeFares)[c(2:5)]=c("moduleColor_Fares","moduleLable_Fares",
                                "kME_M3_Fares","pvalBH_M3_Fares")

#subset columns Sauvageau 
redkmeSauvageau <- study_data$Sauvageau_2014$kme[,c (1,4,5,9:12,31:34)]
colnames(redkmeSauvageau)[c(2:9)]=c("moduleColor_Sauvageau","moduleLable_Sauvageau",
                                    "kME_M2_Sauvageau","pvalBH_M2_Sauvageau", 
                                    "kME_M3_Sauvageau","pvalBH_M3_Sauvageau", 
                                    "kME_M13_Sauvageau","pvalBH_M13_Sauvageau")

#subset columns Xie Culture 
redkmeXieCult = study_data$Xie_2019$kme[,c (1,4,5, 15,16,25:28,33,34,39,40,53,54)]
colnames(redkmeXieCult)[c(2:13)]=c("moduleColor_XieCult","moduleLable_XieCult",
                                   "kME_M5_XieT2","pvalBH_M5_XieT2", 
                                   "kME_M10_XieT2","pvalBH_10_XieT2", 
                                   "kME_M11_XieT2","pvalBH_11_XieT2", 
                                   "kME_M14_XieT2","pvalBH_M14_XieT2",
                                   "kME_M17_XieT2","pvalBH_17_XieT2" )




#create a table with module assigment and kme info from all the studies
LTgenes = left_join(LTgenes, redkmeAnj,by="Ensembl")
LTgenes = left_join(LTgenes,redkmeXie,by="Ensembl")
LTgenes=left_join(LTgenes,redkmeFares,by="Ensembl")
LTgenes=left_join(LTgenes,redkmeSauvageau,by="Ensembl")
LTgenes=left_join(LTgenes,redkmeXieCult,by="Ensembl")

#re-organize the table, so LTmodules columns go to the end of the table
ltmodule_cols <- grep("^LTmodule", names(LTgenes), value = TRUE)
other_cols <- setdiff(names(LTgenes), ltmodule_cols)
LTgenes=LTgenes[, c(other_cols, ltmodule_cols)]

#create columns to sum up the number of times a gene was assigned to LT-associated module
#create separate columns to Fresh and Cultured
#first, we need to change NA to 0
LTgenes[, 44:53] <- lapply(
  LTgenes[, 44:53],
  function(x) ifelse(is.na(x), 0, x)
)

# Overlap between fresh datasets
LTgenes$CountFresh = LTgenes$LTmodule_Anjos + LTgenes$LTmodule_XieT1

# Overlap between expanded datasets
LTgenes$CountExpanded = LTgenes$LTmodule_Fares +
  LTgenes$LTmodule_Sauvageau + LTgenes$LTmodule_XieT2

# Overlap between fresh and expanded datasets
LTgenes$CountTotal = ifelse(LTgenes$CountFresh == 0, "NA", LTgenes$CountFresh + LTgenes$CountExpanded)

# Overlap in fresh datasets considering kME > 0.8
LTgenes$CountFresh0.8 = LTgenes$LTmodule_Anjos_0.8 + LTgenes$LTmodule_XieT1_0.8

# Overlap in expanded datasets considering kME > 0.8
LTgenes$CountExpanded0.8 = LTgenes$LTmodule_Fares_0.8 +
  LTgenes$LTmodule_Sauvageau_0.8 + LTgenes$LTmodule_XieT2_0.8

# Overlap of fresh and expanded datasets considering kME > 0.8
LTgenes$CountTotal0.8All = ifelse(LTgenes$CountFresh0.8 == 0, "NA", LTgenes$CountFresh0.8 + LTgenes$CountExpanded0.8)

# Save the LTgenes data frame (with the new 'CountExpanded' column)
write.csv(LTgenes, "Overlap_Uncult_Cult/Overlap/2_overlap_cult_kme_info_all.csv")

# Create Venn Diagram for positive modules

# Extract Ensembl IDs from filtered data
genesAnjosListPos <- LTAnjos$Ensembl
genesXieListPos_T1 <- LTXie_T1$Ensembl
genesFaresListPos <-  LTFares$Ensembl
genesSauvageauListPos <- LTSauv$Ensembl
genesXieListPos_T2 <- LTXie_T2$Ensembl

# Create Venn Diagram for positive modules
venn.plot <- venn.diagram(
  x = list(
    Anjos      = unique(genesAnjosListPos),
    Xie_Type1  = unique(genesXieListPos_T1),
    Fares      = unique(genesFaresListPos),
    Sauvageau  = unique(genesSauvageauListPos),
    Xie_Type2  = unique(genesXieListPos_T2)  ),
  filename = "Overlap_Uncult_Cult/Overlap/Diagrama_Venn_Positive.png",
  output = TRUE,
  category.names = c(
    paste0("Anjos\n(",length(unique(genesAnjosListPos)), ")"),
    paste0("Xie_Type1\n(",length(unique(genesXieListPos_T1)), ")"),
    paste0("Fares\n(", length(unique(genesFaresListPos)), ")"),
    paste0("Sauvageau\n(", length(unique(genesSauvageauListPos)), ")"),
    paste0("Xie_Type2\n(", length(unique(genesXieListPos_T2)), ")")
  ),
  col = NA,
  fill = c("#FF0000", "#FFA500", "#FFD700", "#FF6347", "#FFF749"),
  alpha = 0.2,
  lwd = 1,
  cex = 0.9,
  fontface = "plain",
  fontfamily = "sans",
  cat.cex = 1,
  cat.col = c("#B22222", "#FF8C00", "#DAA520", "#CD5C5C", "#BDB76B"),
  cat.fontface = "plain",
  cat.fontfamily = "sans",
  cat.pos = c(-10, -90, 180, 0, 90),
  cat.dist = c(0.25, 0.3, 0.2, -0.22, 0.27),
  margin = 0.2,
  sep.dist = 0.01,
  scaled = FALSE
)

# 2.2.2. Filter negatively correlated modules

#Anjos
modulesAnjos_negative <- study_data$`Anjos-Afonso_etal_2021`$modules$moduleColor[study_data$`Anjos-Afonso_etal_2021`$modules$CorrelationDirection == "Negative"]
modulesAnjos_clean_negative <- sub("_Anjos-Afonso_etal_2021", "", modulesAnjos_negative)
modulesAnjos_clean_negative <- setdiff(modulesAnjos_clean_negative, c("yellow", "green"))
study_data$`Anjos-Afonso_etal_2021`$kme$LTmodule_Anjos_Neg <- ifelse(study_data$`Anjos-Afonso_etal_2021`$kme$moduleColor %in% c("brown"), 1, 0)
LTAnjos_Neg = study_data$`Anjos-Afonso_etal_2021`$kme[study_data$`Anjos-Afonso_etal_2021`$kme$moduleColor %in% c("brown"), 1:3]

#Xie
modulesXie_Type1_negative <- study_data$Xie_etal_2019$modules$moduleColor[study_data$Xie_etal_2019$modules$CorrelationDirection == "Negative"]
modulesXie_Type1_clean_negative <- sub("_Xie_etal_2019", "", modulesXie_Type1_negative)
modulesXie_Type1_clean_negative <- setdiff(modulesXie_Type1_clean_negative, c("yellow"))
study_data$Xie_etal_2019$kme$LTmodule_XieT1_Neg <- ifelse(study_data$Xie_etal_2019$kme$moduleColor %in% c("darkgrey", "lightgreen", "tan", "turquoise"), 1, 0)
LTXie_T1_Neg = study_data$Xie_etal_2019$kme[study_data$Xie_etal_2019$kme$moduleColor %in% c("darkgrey", "lightgreen", "tan", "turquoise"), 1:3]


#Fares
modulesFares_negative <- study_data$Fares_2017$modules$moduleColor[study_data$Fares_2017$modules$CorrelationDirection == "Negative"]
modulesFares_clean_negative <- sub("_Fares_2017", "", modulesFares_negative)
modulesFares_clean_negative <- setdiff(modulesFares_clean_negative, c("greenyellow", "red"))
study_data$Fares_2017$kme$LTmodule_Fares_Neg <- ifelse(study_data$Fares_2017$kme$moduleColor %in% c("yellow"), 1, 0)
LTFares_Neg = study_data$Fares_2017$kme[study_data$Fares_2017$kme$moduleColor %in% c("yellow"), 1:3]

#Sauvageau 
modulesSauvageau_negative <- study_data$Sauvageau_2014$modules$moduleColor[study_data$Sauvageau_2014$modules$CorrelationDirection == "Negative"]
modulesSauvageau_clean_negative <- sub("_Sauvageau_2014", "", modulesSauvageau_negative)
modulesSauvageau_clean_negative <- setdiff(modulesSauvageau_clean_negative, "yellow")
study_data$Sauvageau_2014$kme$LTmodule_Sauv_Neg = ifelse(study_data$Sauvageau_2014$kme$moduleColor %in% c("green", "turquoise"), 1, 0)
LTSauv_Neg = study_data$Sauvageau_2014$kme[study_data$Sauvageau_2014$kme$moduleColor %in% c("green", "turquoise"), 1:3]

#Xie Culture 
modulesXie_Type2_negative <- study_data$Xie_2019$modules$moduleColor[study_data$Xie_2019$modules$CorrelationDirection == "Negative"]
modulesXie_Type2_clean_negative <- sub("_Xie_2019", "", modulesXie_Type2_negative)
modulesXie_Type2_clean_negative <- setdiff(modulesXie_Type2_clean_negative, c( "darkgreen", "darkturquoise", "lightgreen", "magenta"))
study_data$Xie_2019$kme$LTmodule_XieT2_Neg = ifelse(study_data$Xie_2019$kme$moduleColor %in% c("brown"), 1, 0)
LTXie_T2_Neg = study_data$Xie_2019$kme[study_data$Xie_2019$kme$moduleColor %in% c("brown"), 1:3]

#merge all genes in LT-associated modules and remove duplicates
LTgenes_Neg =rbind(LTAnjos_Neg, LTXie_T1_Neg, LTFares_Neg,LTSauv_Neg,LTXie_T2_Neg)
LTgenes_Neg= LTgenes_Neg %>% distinct(Ensembl, .keep_all = TRUE)

#create a table holding the module assignment of the gene in each study
#plus the columns with kme info (only for the LT-associated modules in that study)
#and the column that says if the gene is or not in LT-associated module in that study


redkmeAnj_neg  <- study_data$`Anjos-Afonso_etal_2021`$kme[,c (1:5,11,12,43)]
colnames(redkmeAnj_neg)[c(2:7)]=c("Gene", "Gene_type", "moduleColor_Anjos","moduleLable_Anjos",
                                       "kME_M3_Anjos","pvalBH_M3_Anjos")

redkmeXie_T1_neg  <- study_data$Xie_etal_2019$kme[,c (1:5,7,8,29,30,41,42,53,54,63)]
colnames(redkmeXie_T1_neg )[c(2:13)]=c("Gene", "Gene_type", "moduleColor_XieT1","moduleLable_XieT1",
                                  "kME_M1_XieT1","pvalBH_M1_XieT1",
                                  "kME_M12_XieT1","pvalBH_M12_XieT1",
                                  "kME_M18_XieT1","pvalBH_M18_XieT1",
                                  "kME_M24_XieT1","pvalBH_M24_XieT1")

redkmeFares_neg  <- study_data$Fares_2017$kme [,c (1:5, 13,14,45)]
colnames(redkmeFares_neg)[c(2:7)]=c("Gene", "Gene_type", "moduleColor_Fares","moduleLable_Fares",
                                     "kME_M4_Fares","pvalBH_M4_Fares")

redkmeSauvageau_neg <- study_data$Sauvageau_2014$kme [,c (1:5, 7,8,15,16,35)]
colnames(redkmeSauvageau_neg)[c(2:9)]=c("Gene", "Gene_type", "moduleColor_Sauvageau","moduleLable_Sauvageau",
                                        "kME_M1_Sauvageau","pvalBH_M1_Sauvageau",
                                        "kME_M5_Sauvageau","pvalBH_M5_Sauvageau")

redkmeXie_T2_neg  <- study_data$Xie_2019$kme[,c (1:5,11,12,55)]
colnames(redkmeXie_T2_neg)[c(2:7)]=c("Gene", "Gene_type", "moduleColor_XieCult","moduleLable_XieCult",
                                   "kME_M3_XieT2","pvalBH_M3_XieT2")

#create a table with module assigment and kme info from all the studies
LTgenes_Neg=left_join(LTgenes_Neg,redkmeAnj_neg,by="Ensembl")
LTgenes_Neg=left_join(LTgenes_Neg,redkmeXie_T1_neg,by="Ensembl")
LTgenes_Neg=left_join(LTgenes_Neg,redkmeFares_neg,by="Ensembl")
LTgenes_Neg=left_join(LTgenes_Neg,redkmeSauvageau_neg,by="Ensembl")
LTgenes_Neg=left_join(LTgenes_Neg,redkmeXie_T2_neg,by="Ensembl")

#re-organize the table, so LTmodules columns go to the end of the table
ltmodule_cols <- grep("^LTmodule", names(LTgenes_Neg), value = TRUE)
other_cols <- setdiff(names(LTgenes_Neg), ltmodule_cols)
LTgenes_Neg=LTgenes_Neg[, c(other_cols, ltmodule_cols)]
LTgenes_Neg = LTgenes_Neg[,-c(4,5,10,11,22,23,28,29,36,37)]

#create columns to sum up the number of times a gene was assigned to LT-associated module
#create separate columns to Fresh and Cultured
#first, we need to change NA to 0
LTgenes_Neg[, 32:36] <- lapply(
  LTgenes_Neg[, 32:36],
  function(x) ifelse(is.na(x), 0, x)
)

# Overlap analysis for cultivated datasets
LTgenes_Neg$CountTotal= LTgenes_Neg$LTmodule_Anjos_Neg + LTgenes_Neg$LTmodule_XieT1_Neg +
  LTgenes_Neg$LTmodule_Fares_Neg + LTgenes_Neg$LTmodule_Sauv_Neg + LTgenes_Neg$LTmodule_XieT2_Neg


# Save the LTgenes data frame (with the new 'CountExpanded' column)
write.csv(LTgenes_Neg,file="Overlap_Uncult_Cult/Overlap/2_overlap_cult_neg_Final.csv", row.names=FALSE)

# Create Venn Diagram for negative modules

# Extract Ensembl IDs from filtered data
genesAnjosListNeg <- LTAnjos_Neg$Ensembl
genesXieListNeg_T1 <- LTXie_T1_Neg$Ensembl
genesFaresListNeg <-  LTFares_Neg$Ensembl
genesSauvageauListNeg <- LTSauv_Neg$Ensembl
genesXieListNeg_T2 <- LTXie_T2_Neg$Ensembl

# Adjusted colors for 3 sets
venn_colors <- c(
  "#1B9E77",
  "#66A61E",
  "#A6D854",
  "#377EB8", 
  "#08519C"  
)

# Create Venn Diagram for positive modules
venn.plot <- venn.diagram(
  x = list(
    Anjos      = unique(genesAnjosListNeg),
    Xie_Type1  = unique(genesXieListNeg_T1),
    Fares      = unique(genesFaresListNeg),
    Sauvageau  = unique(genesSauvageauListNeg),
    Xie_Type2  = unique(genesXieListNeg_T2)  ),
  filename = "Overlap_Uncult_Cult/Overlap/Diagrama_Venn_Negative.png",
  output = TRUE,
  category.names = c(
    paste0("Anjos\n(",length(unique(genesAnjosListNeg)), ")"),
    paste0("Xie_Type1\n(",length(unique(genesXieListNeg_T1)), ")"),
    paste0("Fares\n(", length(unique(genesFaresListNeg)), ")"),
    paste0("Sauvageau\n(", length(unique(genesSauvageauListNeg)), ")"),
    paste0("Xie_Type2\n(", length(unique(genesXieListNeg_T2)), ")")
  ),
  col = NA,
  fill = venn_colors,
  alpha = 0.2,
  lwd = 1,
  cex = 0.9,
  fontface = "plain",
  fontfamily = "sans",
  cat.cex = 1,
  cat.col = c(  "#1B9E77",
                "#66A61E",
                "#A6D854",
                "#377EB8", 
                "#08519C"  ),
  cat.fontface = "plain",
  cat.fontfamily = "sans",
  cat.pos = c(-10, -90, 180, 0, 90),
  cat.dist = c(0.25, 0.3, 0.2, -0.22, 0.27),
  margin = 0.2,
  sep.dist = 0.01,
  scaled = FALSE
)


