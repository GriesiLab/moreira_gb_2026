###############################################################################
# Script:        Quality Control (RNA-seq) - Sauvageau_2014 (GSE57299)
# Author:        Gustavo Bueno Moreira
# Last update:   2026-02-16
#
# Purpose:
#   - Perform exploratory QC for RNA-seq count data (PCA, clustering, heatmap),
#     and optionally explore batch correction with RUVseq.
#
# Deliverables:
#   1) QC report (boxplot, PCA, heatmap, clustering) exported to PDF
#   2) RUVseq exploratory normalization panel exported to PDF (p- and k-grid)
#   3) PCA outputs: variance explained, PC interpretation heatmap, scree plot,
#      and top contributing genes per PC exported as tables/figures
#
# Inputs:
#   - GEO series: GSE57299 (metadata via GEOquery)
#   - DATA/GSE57299_raw_counts_GRCh38.p13_NCBI.tsv
#   - DATA/exonicLength_hg19.rda OR DATA/exonicLength_hg38.rda
#   - DATA/HK_full_gene_list_ensembl_biomart.csv
#   - Internet access for GEOquery / biomaRt (if conversion is requested)
#
# Outputs:
#   - DATA/metadata.csv
#   - DATA/raw_counts.csv
#   - DATA/counts_NonNormalized.csv
#   - DATA/counts_Final.csv
#   - RESULTS/Quality_Control/3_QC_withoutOutliers.pdf
#   - RESULTS/Quality_Control/4_Variance_explained.csv
#   - RESULTS/Quality_Control/4_PC_interpretation_plot.png
#   - RESULTS/Quality_Control/4_Scree_plot.png
#   - RESULTS/Quality_Control/4_Key_contributing_PC1.csv
#   - RESULTS/Quality_Control/4_Key_contributing_PC2.csv
#   - RESULTS/LowVar.csv
#   - RUVseq_results_p_<p>.pdf (multiple files)
#
# Notes:
#   - This script assumes execution from the study folder set via setwd() below.
#   - `filter_genes_by_tpm()` requires `metadata$Genome_build` containing
#     "GRCh19" or "GRCh38" (case-insensitive matches handled).
###############################################################################

#### 0) Setup #################################################################

# 0.1) Working directory -------------------------------------------------------
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Fares_2014/")

# 0.2) Options -----------------------------------------------------------------
options(stringsAsFactors = FALSE)

# 0.3) Packages ----------------------------------------------------------------
library(grid)
library(tidyr)
library(ggplot2)
library(GEOquery)
library(RColorBrewer)
library(RUVSeq)
library(FactoMineR)
library(factoextra)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(DESeq2)
library(biomaRt)
library(stringr)

#### 0.4) Functions (helpers) ##################################################

detect_separator <- function(file_path) {
  first_line <- readLines(file_path, n = 1)
  if (grepl(";", first_line)) return(";")
  if (grepl(",", first_line)) return(",")
  if (grepl("\t", first_line)) return("\t")
  stop("Separator not recognized. Check the file format.")
}

check_id_type <- function(ids) {
  if (all(grepl("^ENSG\\d{11}$", ids))) return("Ensembl")
  if (all(grepl("^\\d+$", ids))) return("Entrez")
  return("GeneName")
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
    
    geneR <- getBM(attributes = attributes, filters = filters, values = expDataRaw[, 1], mart = mart)
    colnames(geneR) <- c("Ensembl", "GeneID")
    geneR$GeneID <- as.character(geneR$GeneID)
    
    expDataRaw <- expDataRaw %>%
      mutate(GeneID = as.character(expDataRaw[, 1])) %>%
      left_join(geneR, by = "GeneID")
    
    missing_ensembl <- sum(is.na(expDataRaw$Ensembl))
    message("Number of genes without a match in Ensembl: ", missing_ensembl)
    
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
  for (j in seq_len(ncol(rpk))) rpk[, j] <- expData[, j] / gene_length
  rpk <- rpk[!is.na(rpk[, 1]), ]
  
  libsize <- colSums(rpk) / 10^6
  tpm <- rpk
  for (j in seq_len(ncol(tpm))) tpm[, j] <- rpk[, j] / libsize[j]
  
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
  transformed_data <- as.data.frame(assay(vsd))
  return(transformed_data)
}

heatmap_plot <- function(expData, metadata, title = "Heatmap", label_column = NULL) {
  if (is.null(label_column)) {
    cat("Available columns in metadata:\n"); print(colnames(metadata))
    label_column <- readline(prompt = "Enter the column name to use for row labels: ")
  }
  if (!label_column %in% colnames(metadata)) stop(paste("Column", label_column, "does not exist in metadata."))
  
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
    cat("Available columns in metadata:\n"); print(colnames(metadata))
    label_column <- readline(prompt = "Enter the column name to use for sample labels: ")
  }
  if (!label_column %in% colnames(metadata)) stop(paste("Column", label_column, "does not exist in metadata."))
  
  sampleTrees <- hclust(dist(t(expData)), method = method)
  sampleTrees$labels <- metadata[[label_column]]
  
  par(mfrow = c(2, 1))
  par(mar = c(0, 4, 2, 0))
  plot(sampleTrees, main = "Sample Clustering", xlab = "", sub = "", cex = cex_labels)
}

perform_pca <- function(expData, metadata, population_column = NULL, sample_column = NULL) {
  if (is.null(population_column)) {
    cat("Available columns in metadata:\n"); print(colnames(metadata))
    population_column <- readline(prompt = "Enter the column name to define colors (Population) or type NA: ")
    if (toupper(population_column) == "NA") population_column <- NA
  }
  
  if (is.null(sample_column)) {
    cat("Available columns in metadata:\n"); print(colnames(metadata))
    sample_column <- readline(prompt = "Enter the column name to define shapes (Sample) or type NA: ")
    if (toupper(sample_column) == "NA") sample_column <- NA
  }
  
  if (is.na(population_column) & is.na(sample_column)) stop("You must provide at least one variable for the plot!")
  
  mtx4PCA <- t(expData)
  
  pca_results <- PCA(mtx4PCA, scale.unit = FALSE, ncp = 20, graph = FALSE)
  
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
      x = paste0("PC1 (", variance_pc1, "%)"),
      y = paste0("PC2 (", variance_pc2, "%)")
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
    pop_order <- c("UM171+SR1", "UM171", "SR1", "DMSO")
    pca_data$Population <- factor(pca_data$Population, levels = pop_order)
    
    num_colors <- length(unique(na.omit(pca_data$Population)))
    color_palette <- colorRampPalette(c("red", "blue"))(num_colors)
    
    pca_plot <- pca_plot +
      aes(color = Population) +
      scale_color_manual(values = setNames(color_palette, levels(pca_data$Population))) +
      labs(color = "Population")
  }
  
  if (!is.na(sample_column)) {
    pca_plot <- pca_plot + aes(shape = Sample) + labs(shape = "Sample")
  }
  
  return(list(pca_data = pca_data, pca_plot = pca_plot, pca_results = pca_results))
}

differential_expression_analysis <- function(expData, metadata) {
  cat("Available columns in metadata:\n"); print(colnames(metadata))
  design_var <- readline(prompt = "Enter the column name to use for the design: ")
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
  
  if (p < 0 || p > 1) {
    stop("The p-value must be between 0 and 1.")
  }
  
  if (k < 1 || k > 5) {
    stop("The k-value must be between 1 and 5.")
  }
  
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


#### 0.5) Load data ############################################################

message("0.5) Downloading GEO metadata and reading expression matrix...")

geo_id <- "GSE57299"

gse <- tryCatch(
  { getGEO(geo_id, GSEMatrix = TRUE) },
  error = function(e) { stop("Failed to download GEO data. Check the geo_id and internet connection.") }
)

metadata <- pData(phenoData(gse[[1]]))

# Standardize key columns (study-specific mapping)
names(metadata)[names(metadata) == "modulating agent:ch1"] <- "LT_enrichment_fct"
names(metadata)[names(metadata) == "data_processing.3"]     <- "Genome_build"
names(metadata)[names(metadata) == "growth time:ch1"]       <- "Time"

metadata <- metadata %>%
  mutate(
    LT_enrichment = case_when(
      LT_enrichment_fct == "UM171+SR1" ~ 4,
      LT_enrichment_fct == "UM171"     ~ 3,
      LT_enrichment_fct == "SR1"       ~ 2,
      LT_enrichment_fct == "DMSO"      ~ 1,
      TRUE ~ NA_real_
    ),
    LT_enrichment_fct = if_else(LT_enrichment_fct == "none", "Fresh", LT_enrichment_fct),
    title = str_remove(title, "^CB_"),
    Time = if_else(is.na(Time), "0h", as.character(Time)),
    Culture_Time = case_when(
      Time == "12hr" ~ 1,
      Time == "24hr" ~ 2,
      Time == "48hr" ~ 3,
      Time == "72hr" ~ 4,
      TRUE ~ NA_real_
    ),
    Imunophenotype = "CD34+/CD45RA-"
  )

metadata$Genome_build <- "GRCh38"

# Expression matrix ------------------------------------------------------------
exp_file <- "DATA/GSE57299_raw_counts_GRCh38.p13_NCBI.tsv"
separator <- detect_separator(exp_file)
expDataRaw <- read.delim(exp_file, sep = separator)

expData <- convert_to_ensembl(expDataRaw)
rownames(expData) <- expData$Ensembl
expData <- expData[, -c(1, 27)]

message("0.5) Validating metadata rownames vs expression colnames...")
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

# Save inputs snapshot ---------------------------------------------------------
expData1 <- expData
write.csv(metadata, file = "DATA/metadata.csv")
write.csv(expData,   file = "DATA/raw_counts.csv")

#### 1) Deliverable 1: Data preprocessing #####################################

message("1) Filtering and preprocessing data...")

dds <- DESeqDataSetFromMatrix(countData = expData, colData = metadata, design = ~ 1)
expData <- dds[rowSums(counts(dds)) > 0, ]
expData <- as.data.frame(assay(expData))

expData <- filter_genes_by_tpm(expData, metadata)
expData2a <- expData
write.csv(expData2a, "DATA/counts_NonNormalized.csv")

#### 2) Deliverable 2: Quality control ########################################

message("2) Running QC (VST, boxplot, PCA, heatmap, clustering)...")

ncvsd <- deseq_vst(expData2a, metadata)

boxplot_plot_obj <- boxplot_plot(ncvsd, metadata)
print(boxplot_plot_obj)

pca_results <- perform_pca(ncvsd, metadata)
pca_plot_obj <- pca_results$pca_plot
print(pca_plot_obj)

heatmap_plot_obj <- heatmap_plot(ncvsd, metadata)
sampleClustering_plot(ncvsd, metadata)

message("2) Exporting QC panel to PDF...")
pdf("RESULTS/Quality_Control/3_QC_withoutOutliers.pdf", width = 8, height = 6)
print(boxplot_plot_obj)
print(pca_plot_obj)
grid.newpage()
print(heatmap_plot_obj)
sampleClustering_plot(ncvsd, metadata)
dev.off()

#### 3) Deliverable 3: RUVseq normalization & low-variance filtering ###########

message("3) Running DESeq-based screening and RUVseq exploration...")

metadata$LT_enrichment <- as.factor(metadata$LT_enrichment)
res <- differential_expression_analysis(expData2a, metadata)

HKgenes_full <- read.csv("DATA/HK_full_gene_list_ensembl_biomart.csv")

ruvseq_analysis(res, expData2a, metadata, HKgenes_full)

# Metrics: Batch Correction p-value Threshold: 0.8;  k Value: 4
expData <- normalize_data(res, expData2a, metadata, HKgenes_full$Ensembl_biomart)

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
clinical <- metadata[match(colnames(expData), rownames(metadata)), ]

clinical_subset <- as.data.frame(clinical[, c(42, 43)])
interpret_PC(
  component_matrix = results_pca$pca_results$ind$coord,
  clinical_df = clinical_subset,
  title = "Expression",
  silent = FALSE
)

pc_interpretation_plot <- interpret_PC(
  component_matrix = results_pca$pca_results$ind$coord,
  clinical_df = clinical_subset,
  title = "Expression",
  silent = FALSE
)

ggsave(
  "RESULTS/Quality_Control/4_PC_interpretation_plot.png",
  plot = pc_interpretation_plot,
  width = 8, height = 6, dpi = 300
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

