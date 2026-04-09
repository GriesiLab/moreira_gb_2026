###############################################################################
# Script:        QC Cesana_2018_T3
# Author:        Gustavo Bueno Moreira
# Last update:   2026-02-16
#
# Purpose:
#   - Perform an exploratory analysis of RNA-seq data, summarizing sample-level
#     structure and potential technical variation through descriptive QC plots.
#
# Deliverables:
#   1) QC panel (boxplot, PCA, heatmap, hierarchical clustering) exported to PDF
#   2) Processed metadata and raw counts exported as CSV files
#   3) Filtered count matrix (final) exported as CSV file
#
# Inputs:
#   - `expData`: raw count matrix loaded from local text file
#   - `metadata`: sample annotations loaded from GEO
#   - DATA/GSE109089_raw_counts_GRCh38.p13_NCBI.tsv.gz
#   - DATA/exonicLength_hg19.rda OR DATA/exonicLength_hg38.rda
#   - Internet access for GEOquery / biomaRt (if conversion is requested)
#
# Outputs:
#   - DATA/metadata.csv
#   - DATA/raw_counts.csv
#   - RESULTS/Quality_Control/3_QC_WithoutOutlier.pdf
#   - DATA/counts_final.csv
###############################################################################

#### 0) Setup #################################################################

# 0.1) Setup: working directory ------------------------------------------------
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/4a_Panel_validation_inSilico/Studies_QC/RNAseq/Cesana_etal_2018")

# 0.2) Options & constants -----------------------------------------------------
options(stringsAsFactors = FALSE)

# 0.3) Packages ----------------------------------------------------------------
library(grid)          # low-level grid graphics (legends, drawing)
library(tidyr)         # pivot_longer and tidy data helpers
library(ggplot2)       # plotting (PCA scatter, boxplots) and exporting
library(GEOquery)      # download and parse GEO series data
library(reshape2)      # reshaping helpers (legacy, wide/long conversions)
library(RColorBrewer)  # palettes for heatmaps
library(RUVSeq)        # remove unwanted variation (RUVg)
library(FactoMineR)    # PCA computation
library(factoextra)    # PCA visualizations (fviz_eig, fviz_contrib)
library(corrplot)      # correlation matrix visualization
library(pheatmap)      # heatmap visualization
library(patchwork)     # combine ggplot2 plots
library(gridExtra)     # arrangeGrob / grid-based layouts
library(MOFA2)         # multi-omics factor analysis (not used in this script yet)
library(DT)            # interactive tables (not used in this script yet)
library(psych)         # correlation/factor helpers (not used in this script yet)
library(NMF)           # matrix factorization (not used in this script yet)
library(dplyr)         # data manipulation (mutate, join, distinct)
library(DESeq2)        # VST and DESeq-based utilities
library(biomaRt)       # identifier conversion via Ensembl BioMart
library(UpSetR)        # set intersections (not used in this script yet)
library(stringr)       # string processing helpers

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
  
  if (id_type != "Ensembl") {
    message("Converting identifiers to Ensembl ID...")
    
    mart <- useMart("ensembl")
    mart <- useDataset("hsapiens_gene_ensembl", mart)
    
    if (id_type == "Entrez") {
      attributes <- c("ensembl_gene_id", "entrezgene_id")
      filters <- "entrezgene_id"
    } else if (id_type == "GeneName") {
      attributes <- c("ensembl_gene_id", "hgnc_symbol")
      filters <- "hgnc_symbol"
    }
    
    geneR <- getBM(
      attributes = attributes,
      filters = filters,
      values = expDataRaw[, 1],
      mart = mart
    )
    
    colnames(geneR) <- c("Ensembl", "GeneID")
    geneR$GeneID <- as.character(geneR$GeneID)
    
    expDataRaw <- expDataRaw %>%
      mutate(GeneID = as.character(expDataRaw[, 1])) %>%
      left_join(geneR, by = "GeneID")
    
    genes_without_ensembl <- sum(is.na(expDataRaw$Ensembl))
    message("Number of genes without a match in Ensembl: ", genes_without_ensembl)
    
    expDataRaw <- expDataRaw[!is.na(expDataRaw$Ensembl), ]
    
    expr_col <- which(sapply(expDataRaw, is.numeric))[1]
    if (is.na(expr_col)) stop("No numeric column found for expression sorting.")
    
    expDataRaw <- expDataRaw[order(expDataRaw$Ensembl, -expDataRaw[[expr_col]]), ]
    
    duplicated_ensembl <- expDataRaw[duplicated(expDataRaw$Ensembl) | duplicated(expDataRaw$Ensembl, fromLast = TRUE), ]
    expDataRaw <- expDataRaw %>% distinct(Ensembl, .keep_all = TRUE)
    message("Number of duplicates removed in Ensembl: ", nrow(duplicated_ensembl))
    
    duplicated_geneid <- expDataRaw[duplicated(expDataRaw$GeneID) | duplicated(expDataRaw$GeneID, fromLast = TRUE), ]
    expDataRaw <- expDataRaw %>% distinct(GeneID, .keep_all = TRUE)
    message("Number of duplicates removed in GeneID: ", nrow(duplicated_geneid))
    
    message("Final number of genes after conversion and cleaning: ", nrow(expDataRaw))
  }
  
  return(expDataRaw)
}

filter_genes_by_tpm <- function(expData, metadata) {
  genome <- unique(metadata$Genome_build)
  
  if (any(grepl("GRCH19", genome, ignore.case = TRUE))) {
    load("DATA/exonicLength_hg19.rda")
    cat("Reference genome identified: GRCh19\n")
  } else if (any(grepl("GRCH38", genome, ignore.case = TRUE))) {
    load("DATA/exonicLength_hg38.rda")
    cat("Reference genome identified: GRCh38\n")
  } else {
    stop("Genome not identified in metadata. Check the 'Genome_build' column.")
  }
  
  rownames(exonicLength) <- sapply(strsplit(rownames(exonicLength), "\\."), `[`, 1)
  
  gene_length <- as.numeric(exonicLength)
  names(gene_length) <- rownames(exonicLength)
  
  matched_genes <- match(rownames(expData), names(gene_length))
  gene_length <- gene_length[matched_genes] / 1000
  
  rpk <- expData
  for (j in seq_len(ncol(rpk))) {
    rpk[, j] <- expData[, j] / gene_length
  }
  
  rpk <- rpk[!is.na(rpk[, 1]), ]
  
  libsize <- colSums(rpk) / 10^6
  tpm <- rpk
  for (j in seq_len(ncol(tpm))) {
    tpm[, j] <- rpk[, j] / libsize[j]
  }
  
  n <- as.integer(readline(prompt = "Enter the minimum number of samples with TPM > 1: "))
  if (is.na(n) || n <= 0) stop("Invalid value for n. It must be a positive integer.")
  
  tpm_filtered <- tpm[rowSums(tpm > 1) >= n, ]
  expData_filtered <- expData[rownames(expData) %in% rownames(tpm_filtered), ]
  
  cat("Number of genes after filtering:", nrow(expData_filtered), "\n")
  return(expData_filtered)
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
  
  mtx4PCA <- t(expData)
  
  pca_results <- PCA(
    mtx4PCA,
    scale.unit = FALSE,
    ncp = 20,
    graph = FALSE
  )
  
  pca_coordinates <- pca_results$ind$coord
  variance_explained <- pca_results$eig
  variance_pc1 <- round(variance_explained[1, 2], 1)
  variance_pc2 <- round(variance_explained[2, 2], 1)
  
  pca_data <- as.data.frame(pca_coordinates)
  
  if (!is.na(population_column)) {
    if (!population_column %in% colnames(metadata)) stop(paste("Column", population_column, "does not exist in metadata."))
    pca_data <- pca_data %>% mutate(Population = metadata[[population_column]])
  }
  
  if (!is.na(sample_column)) {
    if (!sample_column %in% colnames(metadata)) stop(paste("Column", sample_column, "does not exist in metadata."))
    pca_data <- pca_data %>% mutate(Sample = metadata[[sample_column]])
  }
  
  pca_plot <- ggplot(pca_data, aes(x = Dim.1, y = Dim.2)) +
    geom_point(size = 4) +
    labs(
      title = "",
      x = paste("PC1 (", variance_pc1, "%)", sep = ""),
      y = paste("PC2 (", variance_pc2, "%)", sep = "")
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.key = element_blank(),
      legend.box = "vertical"
    )
  
  if (!is.na(population_column)) {
    num_colors <- length(unique(pca_data$Population))
    color_palette <- colorRampPalette(c("red", "blue"))(num_colors)
    pca_plot <- pca_plot +
      aes(color = Population) +
      scale_color_manual(values = color_palette) +
      labs(color = "Population")
  }
  
  if (!is.na(sample_column)) {
    pca_plot <- pca_plot + aes(shape = Sample) + labs(shape = "Sample")
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

boxplot_plot <- function(ncvsd, metadata, box_color = "steelblue", x_column = NULL, join_column = NULL) {
  if (is.null(x_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    x_column <- readline(prompt = "Enter the column name for the X-axis: ")
  }
  
  if (is.null(join_column)) {
    join_column <- readline(prompt = "Enter the metadata column name that exactly matches the column names in expData:  ")
  }
  
  if (!(x_column %in% colnames(metadata))) {
    stop(paste("Error: Column", x_column, "does not exist in metadata. Please try again."))
  }
  
  if (!(join_column %in% colnames(metadata))) {
    stop(paste("Error: Column", join_column, "does not exist in metadata. Please try again."))
  }
  
  ncvsd_long <- ncvsd %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "expression") %>%
    left_join(metadata, by = c("sample" = join_column))
  
  ggplot(ncvsd_long, aes(x = .data[[x_column]], y = expression, fill = "fixed_color")) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1.8, outlier.color = "black") +
    scale_fill_manual(values = c("fixed_color" = box_color)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      legend.position = "none",
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    xlab("Samples") +
    ylab("Expression Level") +
    ggtitle("Boxplot")
}

ruvseq_analysis <- function(res, expData, metadata, HKgenes_full) {
  cat("Available columns in metadata:\n")
  print(colnames(metadata))
  
  label_column <- readline(prompt = "Enter the name of the column to use for sample labels: ")
  if (!label_column %in% colnames(metadata)) stop(paste("The column", label_column, "does not exist in the metadata."))
  
  population_column <- readline(prompt = "Enter the name of the column to use for defining colors in PCA (Population): ")
  if (!population_column %in% colnames(metadata)) stop(paste("The column", population_column, "does not exist in the metadata."))
  
  sample_column <- readline(prompt = "Enter the name of the column to use for defining shapes in PCA (Sample) or press Enter to skip: ")
  if (sample_column != "" && !sample_column %in% colnames(metadata)) stop(paste("The column", sample_column, "does not exist in the metadata."))
  if (sample_column == "") sample_column <- NA
  
  boxplot_x_axis <- readline(prompt = "Enter the name of the column to use on the X-axis of the boxplot: ")
  if (!boxplot_x_axis %in% colnames(metadata)) stop(paste("The column", boxplot_x_axis, "does not exist in the metadata."))
  
  join_column <- readline(prompt = "Enter the metadata column name that exactly matches the column names in expData:  ")
  if (!join_column %in% colnames(metadata)) stop(paste("The column", join_column, "does not exist in the metadata."))
  
  for (p in seq(0.1, 0.9, by = 0.1)) {
    output_pdf <- paste0("RUVseq_results_p_", p, ".pdf")
    pdf(output_pdf)
    
    nDEGs <- rownames(res[res$pvalue >= p, ])
    controlGenes <- intersect(nDEGs, HKgenes_full$Ensembl_biomart)
    
    for (k in 1:5) {
      ruv_result <- RUVg(as.matrix(expData), controlGenes, k = k)
      N <- ruv_result$normalizedCounts
      
      vsd_data <- deseq_vst(N, metadata)
      
      N_log <- as.data.frame(log2(N + 1))
      print(boxplot_plot(N_log, metadata, box_color = "steelblue", x_column = boxplot_x_axis, join_column = join_column))
      
      pca_result <- perform_pca(vsd_data, metadata, population_column, sample_column)
      print(pca_result$pca_plot + ggtitle(paste("PCA - p =", p, ", k =", k)))
      
      heatmap_plot(vsd_data, metadata, title = paste("Heatmap - p =", p, ", k =", k), label_column = label_column)
      
      sampleClustering_plot(vsd_data, metadata, label_column = label_column, method = "average", cex_labels = 0.7)
    }
    
    dev.off()
  }
}

normalize_data <- function(res, expData, metadata, HKgenes) {
  p <- as.numeric(readline(prompt = "Enter the p-value you want to use (e.g., 0.1): "))
  k <- as.integer(readline(prompt = "Enter the k-value you want to use (e.g., 3): "))
  
  if (p < 0 || p > 1) stop("The p-value must be between 0 and 1.")
  if (k < 1 || k > 5) stop("The k-value must be between 1 and 5.")
  
  nDEGs_deseq <- rownames(res[res$pvalue >= p, ])
  nDEGs_HKgenes <- intersect(nDEGs_deseq, HKgenes)
  controlGenes <- nDEGs_HKgenes
  
  if (length(controlGenes) == 0) {
    stop("No control genes were found. Check the p-value and the housekeeping genes list.")
  }
  
  RUV <- RUVg(as.matrix(expData), controlGenes, k = k)
  N <- RUV$normalizedCounts
  
  expData_normalized <- as.data.frame(N)
  return(expData_normalized)
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

annotate_ensembl_with_genename <- function(pca_contrib, ensembl_col = "Ensembl", mart) {
  ensembl_ids <- unique(pca_contrib[[ensembl_col]])
  
  geneR <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = mart
  )
  
  colnames(geneR) <- c(ensembl_col, "Gene")
  pca_contrib <- left_join(pca_contrib, geneR, by = ensembl_col)
  
  missing_genes <- sum(is.na(pca_contrib$Gene))
  message("Number of Ensembl IDs without Gene Name annotation: ", missing_genes)
  
  return(pca_contrib)
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
  
  pca_contrib <- annotate_ensembl_with_genename(pca_contrib, mart = mart)
  return(pca_contrib)
}

#### 1) Deliverable 1: Load and harmonize data ################################

message("1) Downloading GEO metadata and reading expression matrix...")

# 1.1) Load GEO metadata -------------------------------------------------------
geo_id <- "GSE109089"

gse <- tryCatch(
  { getGEO(geo_id, GSEMatrix = TRUE) },
  error = function(e) { stop("Failed to download GEO data. Check the geo_id and internet connection.") }
)

metadata <- pData(phenoData(gse[[1]]))

# Optional: remove rows not used in downstream steps (original behavior)
metadata <- metadata[-c(19:21), ]

# 1.2) Clean and annotate metadata --------------------------------------------
metadata$LT_enrichment_fct <- metadata$source_name_ch1
names(metadata)[names(metadata) == "data_processing.5"] <- "Genome_build"

metadata <- metadata %>%
  mutate(
    Genome_build = case_when(
      Genome_build == "Genome_build: hg19 (GRCh37)" ~ "GRCH38",
      TRUE ~ as.character(Genome_build)
    ),
    LT_enrichment = case_when(
      LT_enrichment_fct %in% c("BM PROG", "CB PROG", "FL PROG") ~ 0,
      LT_enrichment_fct %in% c("BM HSC", "CB HSC", "FL HSC") ~ 1,
      TRUE ~ NA_real_
    ),
    Imunophenotype = case_when(
      LT_enrichment_fct %in% c("BM PROG", "CB PROG", "FL PROG") ~ "CD34+/CD38-",
      LT_enrichment_fct %in% c("BM HSC", "CB HSC", "FL HSC") ~ "CD34+/CD38-/CD90+/CD45RA-",
      TRUE ~ NA_character_
    ),
    Enrichment = case_when(
      LT_enrichment_fct %in% c("BM PROG", "CB PROG", "FL PROG") ~ "Lower Enrichment",
      LT_enrichment_fct %in% c("BM HSC", "CB HSC", "FL HSC") ~ "Upper Enrichment",
      TRUE ~ NA_character_
    ),
    Cell_Source = case_when(
      str_detect(str_to_lower(`cell type:ch1`), "^bone marrow .* stem cells$")          ~ "Bone Marrow",
      str_detect(str_to_lower(`cell type:ch1`), "^cord blood .* stem cells$")          ~ "Cord Blood",
      str_detect(str_to_lower(`cell type:ch1`), "^fetal liver .* stem cells$")         ~ "Fetal Liver",
      str_detect(str_to_lower(`cell type:ch1`), "^bone marrow .* progenitor cells$")   ~ "Bone Marrow",
      str_detect(str_to_lower(`cell type:ch1`), "^cord blood .* progenitor cells$")    ~ "Cord Blood",
      str_detect(str_to_lower(`cell type:ch1`), "^fetal liver .* progenitor cells$")   ~ "Fetal Liver",
      TRUE ~ NA_character_
    )
  )

# 1.3) Read expression matrix --------------------------------------------------
exp_file <- "DATA/GSE109089_raw_counts_GRCh38.p13_NCBI.tsv.gz"
separator <- detect_separator(exp_file)

expDataRaw <- read.delim(exp_file, sep = separator)

expData <- convert_to_ensembl(expDataRaw)
rownames(expData) <- expData$Ensembl
expData <- expData[, -c(19:21)]

# 1.4) Validate metadata rownames vs expression colnames -----------------------
message("1) Validating metadata rownames vs expression colnames...")
if (all(rownames(metadata) %in% colnames(expData))) {
  message("The column names in expData match the row names in metadata. Proceeding with the analysis.")
} else {
  message("The column names in expData do not match the row names in metadata.")
  response <- readline(prompt = "Check if the row names of metadata match the column names of expData. Do they match? (yes/no) ")
  
  if (tolower(response) == "yes") {
    colnames(expData) <- rownames(metadata)
    message("Columns of expData successfully renamed based on metadata row names.")
  } else if (tolower(response) == "no") {
    message("Please manually reorganize the rows of metadata to match the columns of expData.")
    message("After reorganizing, run the script again.")
  } else {
    message("Invalid response. Please answer 'yes' or 'no'. Proceeding without renaming columns.")
  }
}

# 1.5) Save processed metadata and raw counts ----------------------------------
message("1) Saving processed metadata and raw counts...")
expData1 <- expData

write.csv(metadata, file = "DATA/metadata.csv")
write.csv(expData, file = "DATA/raw_counts.csv")

#### 2) Deliverable 2: Data preprocessing #####################################

message("2) Filtering and preprocessing data...")

# 2.1) Remove genes with zero counts ------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = expData, colData = metadata, design = ~ 1)
expData <- dds[rowSums(counts(dds)) > 0, ]
expData <- as.data.frame(assay(expData))

# 2.2) Remove genes with TPM < 1 ----------------------------------------------
# - The function `filter_genes_by_tpm()` requires `metadata$Genome_build`.
expData <- filter_genes_by_tpm(expData, metadata)
expData2a <- expData

#### 3) Deliverable 3: Quality control ########################################

message("3) Running QC (VST, boxplot, PCA, heatmap, clustering)...")

ncvsd <- deseq_vst(expData, metadata)

# 3.1) Boxplot -----------------------------------------------------------------
boxplot_plot_obj <- boxplot_plot(ncvsd, metadata)
print(boxplot_plot_obj)

# 3.2) PCA ---------------------------------------------------------------------
pca_results <- perform_pca(ncvsd, metadata)
pca_plot_obj <- pca_results$pca_plot
print(pca_results$pca_data)
print(pca_plot_obj)

# 3.3) Heatmap -----------------------------------------------------------------
heatmap_plot(ncvsd, metadata)
heatmap_plot_obj <- heatmap_plot(ncvsd, metadata)

# 3.4) Sample clustering -------------------------------------------------------
sampleClustering_plot(ncvsd, metadata)

message("3) Exporting QC panel to PDF...")
pdf("RESULTS/Quality_Control/3_QC_WithoutOutlier.pdf", width = 8, height = 6)

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
        message("Your dataset is not eligible for validation analysis due to the presence of batch effects.")
        break
      } else if (batch_resp == "no") {
        message("Quality control completed. Save the expData.")
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

# 3.5) Save final count matrix -------------------------------------------------
write.csv(expData, file = "DATA/counts_final.csv")
