###############################################################################
# Script:       Quality Control - Fresh in-house dataset
# Authors:      Gustavo Bueno Moreira / Victoria de Paiva Oliveira
# Last update:  2026-02-16
#
# Purpose:
#   - Perform an exploratory analysis of RNA-seq data from selected studies,
#     summarizing sample-level structure and technical variation through
#     descriptive visualizations and QC outputs.
#
# Deliverables:
#   1) QC report (boxplot, PCA, heatmap, hierarchical clustering) exported to PDF
#   2) RUVseq exploratory normalization panel exported to PDF (p- and k-grid)
#   3) PCA outputs: variance explained, PC interpretation heatmap, scree plot,
#      and top contributing genes per PC exported as tables/figures
#
# Inputs:
#   - `expData`: raw count matrix (genes × samples) loaded from text file
#   - `metadata`: sample annotations loaded from GEO
#   - DATA/GSE154588_combined_gene.exp.txt
#   - DATA/exonicLength_hg19.rda OR DATA/exonicLength_hg38.rda
#   - DATA/HK_full_gene_list_ensembl_biomart.csv
#   - Internet access for GEOquery / biomaRt (if conversion is requested)
#
# Outputs:
#   - DATA/metadata.csv
#   - DATA/raw_counts.csv
#   - RESULTS/1_TPM_genes_retained_log_plot.png
#   - RESULTS/Quality_Control/3_QC.pdf
#   - RESULTS/Quality_Control/4_Variance explained.csv
#   - RESULTS/Quality_Control/4_PC_interpretation_plot.png
#   - RESULTS/Quality_Control/4_shapiro_tests.csv
#   - RESULTS/Quality_Control/RUVg/4_QQplots/4_qqplots_all.pdf
#   - RESULTS/Quality_Control/4_Scree_plot.png
#   - RESULTS/Quality_Control/4_Key_contributing_PC1.csv
#   - RESULTS/Quality_Control/4_Key_contributing_PC2.csv
#   - RESULTS/LowVar.csv
#   - DATA/counts_RUVseq_Normalized.csv
#   - RUVseq_results_p_<p>.pdf (multiple files)
#
# Notes:
#   - This script assumes execution from the study folder set via setwd() below.
#   - `filter_genes_by_tpm()` requires `metadata$Genome_build` containing
#     "GRCh19" or "GRCh38" (case-insensitive matches handled).
#   - `convert_identifiers()` / `get_pca_top_genes()` require an Ensembl mart
#     object (`mart`). A common default is created in 0.5 if not provided.
###############################################################################

#### 0) Setup #################################################################

# 0.1) Setup: working directory ------------------------------------------------
setwd("H:/.shortcut-targets-by-id/1VEKPzB1T9K_ZyE4xteiwNQa0Tpa34pXh/_Projeto_Pi_Fapesp_2022/2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/1_fresh_set/")

# 0.2) Options & constants -----------------------------------------------------
options(stringsAsFactors = FALSE)

# 0.3) Packages ----------------------------------------------------------------
library(grid)          # low-level grid graphics (legends, drawing)
library(tidyr)         # pivot_longer and tidy data helpers
library(ggplot2)       # plotting (PCA scatter, boxplots) and exporting
library(GEOquery)      # download and parse GEO series data
library(RColorBrewer)  # palettes for heatmaps
library(RUVSeq)        # remove unwanted variation (RUVg)
library(FactoMineR)    # PCA computation
library(factoextra)    # PCA visualizations (fviz_eig, fviz_contrib)
library(pheatmap)      # heatmap visualization
library(gridExtra)     # arrangeGrob for combined heatmap + legend
library(dplyr)         # data manipulation (mutate, join, distinct)
library(DESeq2)        # VST and DESeq differential expression results
library(biomaRt)       # identifier conversion via Ensembl BioMart

# 0.4) Functions (helpers) -----------------------------------------------------

detect_separator <- function(file_path) {
  first_line <- readLines(file_path, n = 1)
  
  if (grepl(";", first_line)) {
    return(";")
  } else if (grepl(",", first_line)) {
    return(",")
  } else if (grepl("\t", first_line)) {
    return("\t")
  } else {
    stop("Separator not recognized. Check the file format.")
  }
}

check_id_type <- function(ids) {
  if (all(grepl("^ENSG\\d{11}$", ids))) {
    return("Ensembl")
  } else if (all(grepl("^\\d+$", ids))) {
    return("Entrez")
  } else {
    return("GeneName")
  }
}

convert_to_ensembl <- function(expDataRaw) {
  id_type <- check_id_type(expDataRaw[, 1])
  message("Detected identifier type: ", id_type)
  
  # If already in Ensembl format, return the data unchanged
  if (id_type == "Ensembl") {
    return(expDataRaw)
  }
  
  message("Converting identifiers to Ensembl ID...")
  
  mirrors <- c("www", "useast", "uswest", "asia")
  mart <- NULL
  
  for (mirror in mirrors) {
    message("Trying mirror: ", mirror)
    try({
      mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = mirror)
    }, silent = TRUE)
    if (!is.null(mart)) break
  }
  
  if (is.null(mart)) {
    stop("No Ensembl mirror could be reached. Please try again later.")
  }
  
  if (id_type == "Entrez") {
    attributes <- c("ensembl_gene_id", "entrezgene_id")
    filters <- "entrezgene_id"
  } else if (id_type == "GeneName") {
    attributes <- c("ensembl_gene_id", "hgnc_symbol")
    filters <- "hgnc_symbol"
  }
  
  # Query biomart to retrieve mapping
  geneR <- getBM(attributes = attributes,
                 filters = filters,
                 values = expDataRaw[, 1],
                 mart = mart)
  
  colnames(geneR) <- c("Ensembl", "GeneID")
  geneR$GeneID <- as.character(geneR$GeneID)
  
  expDataRaw <- expDataRaw %>%
    mutate(GeneID = as.character(expDataRaw[, 1])) %>%
    left_join(geneR, by = "GeneID")
  
  message("Number of genes without a match in Ensembl: ", sum(is.na(expDataRaw$Ensembl)))
  
  expDataRaw <- expDataRaw[!is.na(expDataRaw$Ensembl), ]
  
  expr_col <- which(sapply(expDataRaw, is.numeric))[1]
  if (is.na(expr_col)) stop("No numeric column found for expression sorting.")
  
  # Sort and remove duplicates
  expDataRaw <- expDataRaw[order(expDataRaw$Ensembl, -expDataRaw[[expr_col]]), ]
  
  duplicated_ensembl <- expDataRaw[duplicated(expDataRaw$Ensembl) | duplicated(expDataRaw$Ensembl, fromLast = TRUE), ]
  expDataRaw <- expDataRaw %>% distinct(Ensembl, .keep_all = TRUE)
  message("Number of duplicates removed in Ensembl: ", nrow(duplicated_ensembl))
  
  duplicated_geneid <- expDataRaw[duplicated(expDataRaw$GeneID) | duplicated(expDataRaw$GeneID, fromLast = TRUE), ]
  expDataRaw <- expDataRaw %>% distinct(GeneID, .keep_all = TRUE)
  message("Number of duplicates removed in GeneID: ", nrow(duplicated_geneid))
  
  message("Final number of genes after conversion and cleaning: ", nrow(expDataRaw))
  
  return(expDataRaw)
}

deseq_vst <- function(expData, metadata) {
  dds <- DESeqDataSetFromMatrix(countData = expData, colData = metadata, design = ~ 1)
  vsd <- varianceStabilizingTransformation(dds)
  
  transformed_data <- assay(vsd)
  transformed_data <- as.data.frame(transformed_data)
  
  return(transformed_data)
}

heatmap_plot <- function(expData, metadata, title = "Heatmap", label_column = NULL) {
  if (is.null(label_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    label_column <- readline(prompt = "Enter the column name to use for row labels: ")
  }
  
  if (!label_column %in% colnames(metadata)) {
    stop(paste("Column", label_column, "does not exist in metadata."))
  }
  
  sampleDists <- dist(t(expData))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  rownames(sampleDistMatrix) <- metadata[[label_column]]
  colnames(sampleDistMatrix) <- NULL
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    main = title
  )
}

sampleClustering_plot <- function(expData, metadata, method = "average", cex_labels = 0.7, label_column = NULL) {
  if (is.null(label_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    label_column <- readline(prompt = "Enter the column name to use for sample labels: ")
  }
  
  if (!label_column %in% colnames(metadata)) {
    stop(paste("Column", label_column, "does not exist in metadata."))
  }
  
  sampleTrees <- hclust(dist(t(expData)), method = method)
  sampleTrees$labels <- metadata[[label_column]]
  
  par(mfrow = c(2, 1))
  par(mar = c(0, 4, 2, 0))
  
  plot(sampleTrees, main = "Sample Clustering", xlab = "", sub = "", cex = cex_labels)
}

perform_pca <- function(expData, metadata, population_column = NULL, sample_column = NULL) {
  if (!requireNamespace("FactoMineR", quietly = TRUE) || !requireNamespace("factoextra", quietly = TRUE)) {
    stop("Please install 'FactoMineR' and 'factoextra' packages to use this function.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please load or install 'dplyr'.")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Please install the 'viridis' package.")
  
  if (is.null(population_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    population_column <- readline(prompt = "Enter the column name to define colors (Population) or type NA: ")
    if (toupper(population_column) == "NA") population_column <- NA
  }
  if (is.null(sample_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    sample_column <- readline(prompt = "Enter the column name to define shapes (Sample) or type NA: ")
    if (toupper(sample_column) == "NA") sample_column <- NA
  }
  if (is.na(population_column) & is.na(sample_column)) {
    stop("You must provide at least one variable for the plot!")
  }
  
  message("Assuming input data is already normalized or transformed.")
  mtx4PCA <- t(expData)
  
  if (!all(rownames(mtx4PCA) == rownames(metadata))) {
    stop("Row names of metadata must match column names of expData!")
  }
  
  pca_results <- FactoMineR::PCA(mtx4PCA, scale.unit = FALSE, ncp = 20, graph = FALSE)
  pca_coordinates <- pca_results$ind$coord
  variance_explained <- pca_results$eig
  variance_pc1 <- round(variance_explained[1, 2], 1)
  variance_pc2 <- round(variance_explained[2, 2], 1)
  
  pca_data <- as.data.frame(pca_coordinates)
  
  if (!is.na(population_column)) {
    if (!population_column %in% colnames(metadata)) stop(paste("Column", population_column, "does not exist in metadata."))
    pca_data <- dplyr::mutate(pca_data, Population = metadata[[population_column]])
    
    population_levels <- unique(metadata[[population_column]])
    pca_data$Population <- factor(pca_data$Population, levels = population_levels)
  }
  if (!is.na(sample_column)) {
    if (!sample_column %in% colnames(metadata)) stop(paste("Column", sample_column, "does not exist in metadata."))
    pca_data <- dplyr::mutate(pca_data, Sample = metadata[[sample_column]])
  }
  
  pca_plot <- ggplot(pca_data, aes(x = Dim.1, y = Dim.2)) +
    geom_point(size = 4, alpha = 0.7) +
    labs(
      x = paste0("PC1 (", variance_pc1, "%)", sep = ""),
      y = paste0("PC2 (", variance_pc2, "%)", sep = "")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.4),
      aspect.ratio = 1,
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(b = 4)),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "lines"),
      plot.margin = margin(3, 3, 3, 3),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
    )
  
  if (!is.na(population_column)) {
    num_levels <- length(unique(pca_data$Population))
    full_palette <- viridis::viridis(num_levels + 2, option = "inferno")
    trimmed_palette <- full_palette[1:num_levels]  # remove the brightest (last) color(s)
    pca_plot <- pca_plot +
      aes(color = Population) +
      scale_color_manual(values = trimmed_palette) +
      labs(color = "% LT Enrichment")
  }
  
  if (!is.na(sample_column)) {
    sample_levels <- unique(pca_data$Sample)
    custom_shapes <- c(16, 17, 15, 18)
    shape_values <- setNames(custom_shapes[seq_along(sample_levels)], sample_levels)
    
    pca_plot <- pca_plot +
      aes(shape = Sample) +
      scale_shape_manual(values = shape_values) +
      labs(shape = "Donor")
  }
  
  return(list(pca_data = pca_data, pca_plot = pca_plot, pca_results = pca_results))
}

differential_expression_analysis <- function(expData, metadata) {
  cat("Available columns in metadata:\n")
  print(colnames(metadata))
  
  design_var <- readline(prompt = "Enter the column name to use for the design (e.g., 'LT'): ")
  design_formula <- as.formula(paste("~", design_var))
  
  dds <- DESeqDataSetFromMatrix(countData = expData, colData = metadata, design = design_formula)
  dds <- DESeq(dds)
  
  res <- results(dds)
  summary(res)
  
  res_df <- as.data.frame(res)
  return(res_df)
}

boxplot_plot <- function(ncvsd,
                         metadata,
                         box_color   = "darkorange",
                         x_column    = NULL,
                         join_column = NULL,
                         y_limits    = NULL) {
  if (is.null(x_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    x_column <- readline(prompt = "Enter the column name for the X-axis: ")
  }
  
  if (is.null(join_column)) {
    join_column <- readline(
      prompt = "Enter the metadata column name that matches sample names in ncvsd: "
    )
  }
  
  if (!(x_column %in% colnames(metadata))) {
    stop(paste("Error: Column", x_column, "not found in metadata."))
  }
  if (!(join_column %in% colnames(metadata))) {
    stop(paste("Error: Column", join_column, "not found in metadata."))
  }
  
  ncvsd_long <- ncvsd %>%
    pivot_longer(
      cols       = everything(),
      names_to   = "sample",
      values_to  = "expression"
    ) %>%
    left_join(metadata, by = c("sample" = join_column))
  
  if (is.null(y_limits)) {
    y_min <- min(ncvsd_long$expression, na.rm = TRUE)
    y_max <- max(ncvsd_long$expression, na.rm = TRUE)
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }
  
  ggplot(ncvsd_long, aes(x = .data[[x_column]], y = expression, fill = "fixed")) +
    geom_boxplot(
      outlier.shape = 16,
      outlier.size  = 1.2,
      outlier.color = "black"
    ) +
    scale_fill_manual(values = c("fixed" = box_color)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    theme_minimal() +
    theme(
      panel.border    = element_rect(color = "black", fill = NA, size = 0.8),
      legend.position = "none",
      axis.title.y    = element_text(size = 12),
      axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    xlab(x_column) +
    ylab("Expression Level") +
    ggtitle("Boxplot")
}

ruvseq_analysis <- function(res, expData, metadata, HKgenes_full, method = "RUVg") {
  library(RUVSeq)
  library(DESeq2)
  library(Biobase)
  
  cat("Available columns in metadata:\n")
  print(colnames(metadata))
  
  label_column      <- readline("Enter the name of the column to use for sample labels: ")
  population_column <- readline("Enter the name of the column to use for defining colors in PCA (Population): ")
  sample_column     <- readline("Enter the name of the column to use for defining replicate groups (Sample): ")
  boxplot_x_axis    <- readline("Enter the name of the column to use on the X-axis of the boxplot: ")
  join_column       <- readline("Enter the metadata column name that matches column names in expData: ")
  
  for (col in c(label_column, population_column, sample_column, boxplot_x_axis, join_column)) {
    if (!col %in% colnames(metadata)) stop(paste("The column", col, "does not exist in the metadata."))
  }
  
  if (!all(colnames(expData) %in% metadata[[join_column]])) {
    stop("Not all expData column names are present in metadata[[join_column]].")
  }
  
  pd <- metadata
  rownames(pd) <- pd[[join_column]]
  pd <- AnnotatedDataFrame(pd)
  
  replicate_groups <- split(seq_along(pData(pd)[[sample_column]]), pData(pd)[[sample_column]])
  replicate_groups <- Filter(function(x) length(x) >= 2, replicate_groups)
  max_size <- max(lengths(replicate_groups))
  scIdx_mat <- t(sapply(replicate_groups, function(x) c(x, rep(tail(x,1), max_size - length(x)))))
  
  for (p in seq(0.1, 0.9, by = 0.1)) {
    message(paste("Threshold p =", p))
    
    seqSet <- newSeqExpressionSet(counts = as.matrix(expData), phenoData = pd)
    
    orig_names <- rownames(seqSet)
    clean_names <- sub("\\..+$", "", orig_names)
    if (!identical(clean_names, orig_names)) {
      rownames(seqSet) <- clean_names
    }
    
    raw_genes <- rownames(res)[res$pvalue >= p]
    raw_ctrl <- intersect(raw_genes, HKgenes_full$Ensembl_biomart)
    valid_ctrl <- intersect(raw_ctrl, rownames(seqSet))
    message(paste("Control genes: found", length(valid_ctrl), "of", length(raw_ctrl)))
    if (length(valid_ctrl) == 0) stop(paste("No valid control genes for p =", p))
    
    output_pdf <- paste0("RUVseq_results_", method, "_p_", p, ".pdf")
    pdf(output_pdf)
    
    for (k_int in 1:5) {
      message(paste("  Processing k =", k_int))
      
      if (method == "RUVg") {
        ruv_set <- RUVg(x = seqSet, cIdx = valid_ctrl, k = k_int)
      } else if (method == "RUVs") {
        ruv_set <- RUVs(x = seqSet, cIdx = valid_ctrl, k = k_int, scIdx = scIdx_mat)
      } else {
        stop("Invalid method. Choose 'RUVg' or 'RUVs'.")
      }
      
      N <- normCounts(ruv_set)
      N[!is.finite(N)] <- 0
      N_round <- round(N)
      N_round[N_round < 0] <- 0
      N_round[N_round > .Machine$integer.max] <- .Machine$integer.max
      storage.mode(N_round) <- "integer"
      N_round[is.na(N_round)] <- 0
      
      coldata <- pData(pd)[colnames(N_round), , drop = FALSE]
      dds <- DESeqDataSetFromMatrix(countData = N_round, colData = coldata, design = ~1)
      
      vsd <- tryCatch({
        vst(dds, blind = TRUE, fitType = "local")
      }, error = function(e) {
        message("vst falhou, usando varianceStabilizingTransformation local: ", e$message)
        varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
      })
      
      df_log <- as.data.frame(log2(N + 1))
      print(boxplot_plot(df_log, metadata, box_color = "darkorange",
                         x_column = boxplot_x_axis, join_column = join_column))
      
      pca_res <- perform_pca(assay(vsd), metadata, population_column, sample_column)
      print(pca_res$pca_plot + ggtitle(paste("PCA - p =", p, ", k =", k_int)))
      
      print(heatmap_plot(assay(vsd), metadata,
                         title = paste("Heatmap p =", p, ", k =", k_int),
                         label_column = label_column))
      
      print(sampleClustering_plot(assay(vsd), metadata,
                                  label_column = label_column,
                                  method = "average", cex_labels = 0.7))
    }
    
    dev.off()
  }
}

normalize_data <- function(res, expData, metadata, HKgenes, method = "RUVg") {
  library(RUVSeq)
  
  p <- as.numeric(readline(prompt = "Enter the p-value you want to use (e.g., 0.1): "))
  k <- as.integer(readline(prompt = "Enter the k-value you want to use (e.g., 3): "))
  
  if (p < 0 || p > 1) stop("The p-value must be between 0 and 1.")
  if (k < 1 || k > 5) stop("The k-value must be between 1 and 5.")
  
  if (is.data.frame(HKgenes)) HKgenes <- HKgenes$Ensembl_biomart
  
  nDEGs <- rownames(res[res$pvalue >= p, ])
  controlGenes <- intersect(nDEGs, HKgenes)
  if (length(controlGenes) == 0) stop("No control genes found. Check p-value and HK gene list.")
  
  if (method == "RUVg") {
    RUV <- RUVg(as.matrix(expData), cIdx = controlGenes, k = k)
  } else if (method == "RUVs") {
    # Construir scIdx como matriz de replicatas
    replicate_groups <- split(seq_along(metadata$sample), metadata$sample)
    max_size <- max(lengths(replicate_groups))
    scIdx <- t(sapply(replicate_groups, function(x) {
      c(x, rep(NA, max_size - length(x)))
    }))
    RUV <- RUVs(x = as.matrix(expData), cIdx = controlGenes, k = k, scIdx = scIdx)
  } else {
    stop("Invalid method. Choose 'RUVg' or 'RUVs'.")
  }
  
  norm_counts <- as.data.frame(RUV$normalizedCounts)
  colnames(norm_counts) <- colnames(expData)
  rownames(norm_counts) <- rownames(expData)
  return(norm_counts)
}

drawheatmap <- function(svdPV.m, title = NULL, silent = FALSE, show_rownames = TRUE, show_colnames = TRUE, show_legend = TRUE) {
  myPalette <- c("darkred", "red", "orange", "pink", "white")
  breaks.v <- c(-10000, -10, -5, -2, log10(0.05), 0)
  
  rownames(svdPV.m) <- paste("Comp.", 1:dim(svdPV.m)[1])
  
  tmp_plot <- pheatmap(
    log10(t(svdPV.m)),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    breaks = breaks.v,
    color = myPalette,
    silent = TRUE,
    legend = FALSE,
    main = paste("Interpretation of the Components -", title),
    show_rownames = show_rownames,
    show_colnames = show_colnames
  )
  
  if (show_legend) {
    leg <- grid::legendGrob(
      c(
        expression("p < 1x" ~ 10^{-10}),
        expression("p < 1x" ~ 10^{-5}),
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"
      ),
      nrow = 5,
      pch = 15,
      gp = grid::gpar(fontsize = 10, col = myPalette)
    )
  } else {
    leg <- NULL
  }
  
  final_plot <- arrangeGrob(tmp_plot$gtable, leg, ncol = 2, widths = c(5, 1))
  
  if (!silent) {
    plot.new()
    grid::grid.draw(final_plot)
  }
  
  return(final_plot)
}

check_association_CHAMP <- function(x, y) {
  if (!is.numeric(y)) {
    return(kruskal.test(x ~ as.factor(y))$p.value)
  } else {
    return(summary(lm(x ~ y))$coeff[2, 4])
  }
}

interpret_PC <- function(component_matrix, clinical_df, title = NULL, silent = FALSE, show_rownames = TRUE, show_colnames = TRUE, show_legend = TRUE) {
  tmp <- apply(
    component_matrix, 2,
    function(x) sapply(clinical_df, function(y) check_association_CHAMP(x, y))
  )
  
  drawheatmap(
    t(tmp),
    title = title,
    silent = silent,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    show_legend = show_legend
  )
}

get_pca_top_genes <- function(results_pca, pc = 1, top = 100, mart) {
  pca_contrib <- fviz_contrib(
    results_pca$pca_results,
    choice = "var",
    axes = pc,
    top = top
  )$data
  
  pca_contrib <- pca_contrib[order(-pca_contrib$contrib), ]
  
  colnames(pca_contrib)[1] <- "Ensembl"
  pca_contrib$Ensembl <- as.character(pca_contrib$Ensembl)
  
  pca_contrib <- convert_identifiers(pca_contrib, to_id = "GeneName", mart)
  return(pca_contrib)
}

draw_pc_heatmap <- function(assoc_df, pval_adj_df, test_df, output_file = "RESULTS/PC_assoc_heatmap.png") {
  mat <- as.matrix(assoc_df)
  pval <- as.matrix(pval_adj_df)
  tests <- as.matrix(test_df)
  
  labels <- matrix("", nrow = nrow(mat), ncol = ncol(mat),
                   dimnames = dimnames(mat))
  
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      val <- mat[i, j]
      p <- pval[i, j]
      test <- tests[i, j]
      
      sig <- if (!is.na(p) && p < 0.001) "***" else if (!is.na(p) && p < 0.01) "**" else if (!is.na(p) && p < 0.05) "*" else ""
      prefix <- if (test == "Kruskal") "(K)" else if (test == "ANOVA") "(A)" else if (test == "Spearman") "(S)" else ""
      
      labels[i, j] <- paste0(prefix, "\n", round(val, 2), "\n", sig)
    }
  }
  
  logp <- -log10(pval)
  logp[is.infinite(logp)] <- 10
  logp[logp > 10] <- 10
  logp[is.na(logp)] <- 0
  
  # Refined breaks and palette
  breaks.v <- c(0, 1.301029, 2, 5, 10)
  myPalette <- c("#FFF5EB", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603")
  text_colors <- matrix(ifelse(logp >= 4.5, "white", "black"), nrow = nrow(labels), ncol = ncol(labels))
  
  # Garantir diretório
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  # Criar heatmap
  tmp_plot <- pheatmap(
    logp,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = labels,
    number_color = text_colors,
    color = myPalette,
    breaks = breaks.v,
    fontsize_number = 7,
    main = "",  
    silent = TRUE,
    legend = FALSE
  )
  
  leg <- grid::legendGrob(c(expression("p < 1x" ~ 10^{-10}),
                            expression("p < 1x" ~ 10^{-5}),
                            "p < 0.01", 
                            "p < 0.05", 
                            "p > 0.05"), 
                          nrow = 5, pch = 15, gp = grid::gpar(fontsize = 10, col = rev(myPalette)))
  
  title_text <- grid::textGrob("Correlation between PCs and experimental variables\n(* p<0.05, ** p<0.01, *** p<0.001)", gp = gpar(fontsize = 14, fontface = "bold"))
  
  final_plot <- arrangeGrob(
    title_text,
    arrangeGrob(tmp_plot$gtable, leg, ncol = 2, widths = c(5, 0.8)),
    ncol = 1,
    heights = unit(c(1.0, 20), "lines")
  )
  
  png(output_file, width = 1400, height = 1000, res = 150)
  grid::grid.draw(final_plot)
  dev.off()
  
  if (interactive()) try(utils::browseURL(normalizePath(output_file)), silent = TRUE)
  
  return(final_plot)
}

generate_shapiro_outputs <- function(clinical_df, pcs, output_csv = "RESULTS/shapiro_tests.csv", qq_pdf = "RESULTS/QQplots/qqplots_all.pdf") {
  dir.create("RESULTS/QQplots", showWarnings = FALSE, recursive = TRUE)
  dir.create("RESULTS", showWarnings = FALSE, recursive = TRUE)
  
  results <- data.frame(Variable = character(), W_statistic = numeric(), p_value = numeric(), Normality = character())
  
  pdf(qq_pdf, width = 7, height = 6)
  
  test_and_plot <- function(x, name) {
    res <- shapiro.test(x)
    norm <- ifelse(res$p.value > 0.05, "Yes", "No")
    results <<- rbind(results, data.frame(Variable = name, W_statistic = res$statistic, p_value = res$p.value, Normality = norm))
    
    qqnorm(x, main = paste("Q-Q plot for", name))
    qqline(x, col = "red")
  }
  
  for (varname in colnames(clinical_df)) {
    if (is.numeric(clinical_df[[varname]])) {
      test_and_plot(clinical_df[[varname]], varname)
    }
  }
  
  for (i in seq_len(ncol(pcs))) {
    pcname <- colnames(pcs)[i]
    if (is.null(pcname) || pcname == "") pcname <- paste0("PC", i)
    test_and_plot(pcs[, i], pcname)
  }
  
  dev.off()
  
  write.csv(results, file = output_csv, row.names = FALSE)
  
  if (interactive()) try(utils::browseURL(normalizePath(qq_pdf)), silent = TRUE)
  
  return(results)
}

# 0.5) Load data ---------------------------------------------------------------
metadata = read.csv("DATA/metadata.csv", row.names = 1)

metadata$sample_name

samples_to_keep <- c(
  "F0540.0", "F0540.025", "F0540.05", "F0540.10", "F0540.25",
  "F3060.0", "F3060.025", "F3060.05", "F3060.10", "F3060.25",
  "F5652.0", "F5652.025", "F5652.05", "F5652.10", "F5652.25",
  "F8635.0", "F8635.025", "F8635.05", "F8635.10", "F8635.25"
)
metadata <- metadata[metadata$sample_name %in% samples_to_keep, ]
metadata <- metadata[match(samples_to_keep, metadata$sample_name), ]

separator <- detect_separator("DATA/merge-STAR-RSEM-gene-level-expected_count.tsv")

expDataRaw <- read.delim("DATA/merge-STAR-RSEM-gene-level-expected_count.tsv", sep = separator)
all(samples_to_keep %in% colnames(expDataRaw))
expDataRaw <- expDataRaw[, c("id1", "id2", samples_to_keep)]

expData <- convert_to_ensembl(expDataRaw)
rownames(expData) = expData$Ensembl # Set Ensembl IDs as row names
expData = expData[,-c(1,2,23,24)]

if (all(rownames(metadata) %in% colnames(expData))) {
  # If they match, inform the user
  message("The column names in expData match the row names in metadata. Proceeding with the analysis.")
} else {
  # If they don't match, ask the user to confirm
  message("The column names in expData do not match the row names in metadata.")
  message_text <- "Check if the row names of metadata match the column names of expData. Do they match? (yes/no)"
  response <- readline(prompt = message_text)
  
  # Process user's response
  if (tolower(response) == "yes") {
    # If the user confirms, rename the columns of expData using metadata row names
    colnames(expData) <- rownames(metadata)
    message("Columns of expData successfully renamed based on metadata row names.")
  } else if (tolower(response) == "no") {
    # If the user denies, ask them to manually adjust the metadata
    message("Please manually reorganize the rows of metadata to match the columns of expData.")
    message("After reorganizing, run the script again.")
  } else {
    # Handle invalid responses
    message("Invalid response. Please answer 'yes' or 'no'. Proceeding without renaming columns.")
  }
}

expData1 <- expData

write.csv(metadata, file="DATA/metadata.csv")

write.csv(expData, file="DATA/raw_counts.csv")

#### 1) Deliverable 1: Data preprocessing #####################################

message("1) Filtering and preprocessing data...")

# 1.1) Remove genes with zero counts ------------------------------------------
expData <- as.matrix(round(expData))
dds <- DESeqDataSetFromMatrix(countData = expData, colData = metadata, design = ~ 1)
expData <- dds[ rowSums(counts(dds)) > 0, ]
expData = as.data.frame(assay(expData))

# 1.2) Remove genes with TPM < 1 ----------------------------------------------
load("DATA/exonicLength_hg38.rda")  # Use 'exonicLength_hg19.rda' if working with hg19
rownames(exonicLength) <- sapply(strsplit(rownames(exonicLength), "\\."), `[`, 1)
gene_ids <- intersect(rownames(expData), rownames(exonicLength))
expData <- expData[gene_ids, ]
gene_length <- as.numeric(exonicLength[gene_ids, 1]) / 1000
names(gene_length) <- gene_ids

rpk <- sweep(expData, 1, gene_length, FUN = "/")

libsize <- colSums(rpk) / 1e6  # Scale library size to millions
tpm <- sweep(rpk, 2, libsize, FUN = "/")

n_samples <- ncol(tpm)
gene_counts <- sapply(1:n_samples, function(n) sum(rowSums(tpm > 1) >= n))

df <- data.frame(n = 1:n_samples, genes_retained = log(gene_counts))

ggplot(df, aes(x = n, y = genes_retained)) +
  geom_line(color = "red") +
  geom_point(color = "black", size = 1) +
  labs(
    title = "log(number of genes with TPM > 1 in at least n samples)",
    x = "n (minimum number of samples)",
    y = "log(number of genes retained)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

ggsave("RESULTS/1_TPM_genes_retained_log_plot.png", width = 6, height = 4, dpi = 300)

expData <- filter_genes_by_tpm(expData, metadata)
expData2a=expData

#### 2) Deliverable 2: Quality control ########################################

message("2) Running QC (VST, boxplot, PCA, heatmap, clustering)...")

# 2.1) Boxplot -----------------------------------------------------------------
boxplot_plot_obj <- boxplot_plot(ncvsd, metadata)
print(boxplot_plot_obj)

# 2.2) PCA ---------------------------------------------------------------------
pca_results <- perform_pca(ncvsd, metadata)
pca_plot_obj <- pca_results$pca_plot
print(pca_results$pca_data)
print(pca_results$pca_plot)

# 2.3) Heatmap -----------------------------------------------------------------
heatmap_plot(ncvsd, metadata)
heatmap_plot_obj <- heatmap_plot(ncvsd, metadata)

# 2.4) Sample clustering -------------------------------------------------------
sampleClustering_plot(ncvsd, metadata)

message("2) Exporting QC panel to PDF...")
pdf("RESULTS/Quality_Control/3_QC.pdf", width = 8, height = 6)

print(boxplot_plot_obj)
print(pca_plot_obj)
grid.newpage()

print(heatmap_plot_obj)
sampleClustering_plot(ncvsd, metadata)

dev.off()

repeat {
  outlier_resp <- tolower(readline("After the analysis, do you have any outlier samples? (Yes/No): "))
  
  if (outlier_resp == "yes") {
    message("Go back to the metadata and expData, remove the outliers, and run the analysis again.")
    break
  } else if (outlier_resp == "no") {
    repeat {
      batch_resp <- tolower(readline("Does your data have batch effects? (Yes/No): "))
      
      if (batch_resp == "yes") {
        message("Proceed to Run RUVseq for batch effect removal.")
        break
      } else if (batch_resp == "no") {
        num_genes <- nrow(expData)
        
        if (num_genes > 15000) {
          message("Proceed to Remove genes with low variance")
        } else {
          message("Proceed to PCA analysis")
        }
        
        break
      } else {
        message("Invalid response. Please enter 'Yes' or 'No'.")
      }
    }
    break
  } else {
    message("Invalid response. Please enter 'Yes' or 'No'.")
  }
}

#### 3) Deliverable 3: RUVseq normalization & low-variance filtering ###########

message("3) Running DESeq-based screening and RUVseq exploration...")

metadata$LT_enrichment <- as.factor(metadata$LT_enrichment)
res <- differential_expression_analysis(expData2a, metadata)
HKgenes_full <- read.csv("DATA/HK_full_gene_list_ensembl_biomart.csv")
ruvseq_analysis(res, expData2a, metadata, HKgenes_full, method = "RUVg")

# Metrics: Batch Correction p-value Threshold: 0.7;  k Value: 5
expData <- normalize_data(res,
                          expData2a,
                          metadata,
                          HKgenes = HKgenes_full$Ensembl_biomart,
                          method = "RUVg")

message("3) (Optional) Removing low-variance genes...")
ncvsd <- deseq_vst(expData, metadata)

variance <- apply(ncvsd, 1, sd) / apply(ncvsd, 1, mean)
hist(variance)

if (nrow(expData) > 15000) {
  variance_cutoff <- as.numeric(readline(prompt = "Enter the variance cutoff value (e.g., 0.035): "))
  
  keep <- which(apply(ncvsd, 1, sd) / apply(ncvsd, 1, mean) >= variance_cutoff)
  expData <- expData[keep, ]
  expData2b <- expData
  
  dev.off()
  
  genesExp <- rownames(ncvsd)
  lowVarGenes <- setdiff(genesExp, names(keep))
  
  lowVarGenes <- as.data.frame(lowVarGenes)
  lowVarGenes[, 2] <- "grey"
  lowVarGenes[, 3] <- "M0"
  lowVarGenes[, 4] <- "lowVar"
  colnames(lowVarGenes) <- c("Ensembl", "moduleColor", "moduleLabel", "Expression")
  
  write.csv(lowVarGenes, "RESULTS/LowVar.csv")
} else {
  expData2b <- expData
  message(paste("The number of genes is", nrow(expData), "(<= 15,000). Low-variance filtering was skipped."))
}

write.csv(expData, file="DATA/counts_Final.csv")

#### 4) Deliverable 4: PCA analysis & exports #################################

message("4) Running PCA analysis and exporting key outputs...")

ncvsd <- deseq_vst(expData, metadata)
results_pca <- perform_pca(ncvsd, metadata)

# 4.1) Variance explained ------------------------------------------------------
variance_explained <- as.data.frame(results_pca$pca_results$eig)
colnames(variance_explained) <- c("Eigenvalue", "Variance_Explained", "Cumulative_Variance")
write.csv(variance_explained, "RESULTS/Quality_Control/4_Variance explained.csv")

# 4.2) PC association heatmap --------------------------------------------------
clinical = metadata[match(colnames(expData),rownames(metadata)),]

clinical_subset = as.data.frame(clinical[, c(24,28:30)])

clinical_subset$provider  <- as.factor(clinical_subset$provider)
clinical_subset$sorting_batch <- as.factor(clinical_subset$sorting_batch)
clinical_subset$library_batch <- as.factor(clinical_subset$library_batch)
clinical_subset$LT_enrichment <- as.numeric(clinical_subset$LT_enrichment)
sapply(clinical_subset, class)

pcs = results_pca$pca_data[, grep("^Dim\\.", colnames(results_pca$pca_data))]

pc_assoc <- calculate_pc_clinical_association(
  pcs = pcs,
  clinical_df = clinical_subset
)

draw_pc_heatmap(
  assoc_df = pc_assoc$association,
  pval_adj_df = pc_assoc$p_value,
  test_df = pc_assoc$test_type,
  output_file = "RESULTS/Quality_Control/4_PC_interpretation_plot.png"
)

shapiro_results <- generate_shapiro_outputs(
  clinical_df = clinical_subset,
  pcs = pcs,
  output_csv = "RESULTS/Quality_Control/4_shapiro_tests.csv",
  qq_pdf = "RESULTS/Quality_Control/RUVg/4_QQplots/4_qqplots_all.pdf"
)

# 4.3) Scree plot --------------------------------------------------------------
fviz_eig(results_pca$pca_results)

scree_plot <- fviz_eig(results_pca$pca_results)
ggsave(
  "RESULTS/Quality_Control/4_Scree_plot.png",
  plot = scree_plot,
  width = 8, height = 6, dpi = 300
)

# 4.4) Top contributing genes --------------------------------------------------
contr_pc1 <- get_pca_top_genes(results_pca, pc = 1, top = 100, mart)
fviz_contrib(results_pca$pca_results, choice = "var", axes = 1, top = 50)
write.csv(contr_pc1, "RESULTS/Quality_Control/4_Key_contributing_PC1.csv")

contr_pc2 <- get_pca_top_genes(results_pca, pc = 2, top = 100, mart)
fviz_contrib(results_pca$pca_results, choice = "var", axes = 2, top = 50)
write.csv(contr_pc2, "RESULTS/Quality_Control/4_Key_contributing_PC2.csv")

#### Software environment #####################################################
sessionInfo()
