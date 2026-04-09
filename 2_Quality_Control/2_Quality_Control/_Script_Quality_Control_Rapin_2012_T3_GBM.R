###############################################################################
# Script:        QC Rapin_2014
# Author:        Gustavo Bueno Moreira
# Last update:   2026-02-19
#
# Purpose:
#   - Perform quality control analysis on Affymetrix Human Genome U133 Plus 2.0 Array
#     using BrainArray (Ensembl) custom CDF, including RMA normalization of CEL files
#     and exploratory QC plots (boxplot, PCA, heatmap, hierarchical clustering).
#
# Deliverables:
#   1) QC report exported to PDF (boxplot, PCA, heatmap, clustering)
#   2) Normalized expression matrix exported to CSV (Ensembl-mapped)
#   3) Processed metadata exported to CSV
#
# Inputs:
#   - Raw CEL files: DATA/GSE42519/*.CEL
#   - GEO metadata downloaded via GEOquery (GSE42519)
#   - BrainArray custom CDF file:
#       - DATA/HGU133Plus2_Hs_ENSG.cdf
#
# Outputs:
#   - DATA/metadata.csv
#   - DATA/normalized_expression_RMA.csv
#   - RESULTS/_Quality_Control/3_QC.pdf
#
# Notes:
#   - This dataset is U133 Plus 2.0 (3'-IVT). We normalize with affy::justRMA
#     using a BrainArray CDF (Ensembl mapping).
#   - In the original script you forced subset to 16 samples. Here I keep that
#     behavior, but do it SAFELY after checking dimensions.
#   - If your metadata rownames are not the CEL-derived sample names, you must
#     realign (or rename) before plotting.
###############################################################################

#### 0) Setup #################################################################

setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/4a_Panel_validation_inSilico/Studies_QC/Microarray/Rapin_2014")
options(stringsAsFactors = FALSE)

dir.create("RESULTS/_Quality_Control", recursive = TRUE, showWarnings = FALSE)

#### 0.1) Load libraries ######################################################
library(grid)          # Grid graphics system for R, used for low-level graphics control
library(tidyr)         # Tools for tidying data, making it easier to work with in R
library(ggplot2)       # Popular package for creating customizable, high-quality plots
library(GEOquery)      # Package for accessing and analyzing Gene Expression Omnibus (GEO) data
library(reshape2)      # Tools for reshaping data between wide and long formats
library(RColorBrewer)  # Provides color schemes for maps and other graphics
library(FactoMineR)    # Package for multivariate data analysis (e.g., PCA, CA, MCA)
library(factoextra)    # Functions for visualizing the results of multivariate data analysis
library(corrplot)      # Visualization of correlation matrices (e.g., using circles or squares)
library(pheatmap)      # Package for creating heatmaps from data, often used in clustering
library(patchwork)     # Combine multiple ggplot2 plots into a single figure
library(gridExtra)     # Tools to arrange multiple grid-based visual objects (plots)
library(ggplot2)       # Popular package for creating customizable, high-quality plots
library(MOFA2)         # Unsupervised analysis of multi-omics data (Modeling Omics Factor Analysis)
library(DT)            # Interactive tables for displaying data, commonly in Shiny apps
library(psych)         # Tools for psychological research, including correlation and factor analysis
library(NMF)           # Non-negative Matrix Factorization for data analysis (e.g., genomics, image analysis)
library(dplyr)         # Data manipulation package with verbs for filtering, mutating, etc.
library(biomaRt)       # Interface to biological databases (e.g., Ensembl) for genomic queries
library(UpSetR)        # Visualize intersecting sets and their properties
library(data.table)    # High-performance data manipulation
library(affy)          # Affymetrix microarray analysis (CEL, RMA)
library(makecdfenv)    # Custom CDF environments (BrainArray)
library(stringr)       # String manipulation helpers

# If you truly need these BrainArray packages keep them;
# otherwise the .cdf + makecdfenv workflow is enough.
library(pd.hgu133plus2.hs.ensg)
library(hgu133plus2hsensgprobe)

#### 1) Functions #############################################################
# NOTE: These functions assume you already have check_id_type() implemented somewhere.
# If not, remove convert_identifiers()/get_pca_top_genes() or add your check_id_type().

convert_identifiers <- function(expDataRaw, to_id, mart) {
  id_type <- check_id_type(expDataRaw[, 1])
  message("Detected identifier type: ", id_type)
  
  if (to_id == "Ensembl") {
    if (id_type == "Ensembl") {
      message("Identifiers are already Ensembl. No conversion needed.")
      return(expDataRaw)
    }
    
    message("Converting identifiers to Ensembl ID...")
    
    if (id_type == "Entrez") {
      attributes <- c("ensembl_gene_id", "entrezgene_id")
      filters <- "entrezgene_id"
    } else if (id_type == "GeneName") {
      attributes <- c("ensembl_gene_id", "hgnc_symbol")
      filters <- "hgnc_symbol"
    } else {
      stop("Unsupported identifier type for conversion to Ensembl.")
    }
    
    geneR <- getBM(attributes = attributes,
                   filters = filters,
                   values = expDataRaw[, 1],
                   mart = mart)
    
    colnames(geneR) <- c("Ensembl", "GeneID")
    geneR$GeneID <- as.character(geneR$GeneID)
    
    expDataRaw <- expDataRaw %>%
      mutate(GeneID = as.character(expDataRaw[, 1])) %>%
      left_join(geneR, by = "GeneID")
    
    missing_ensembl <- sum(is.na(expDataRaw$Ensembl))
    message("Number of genes without a match in Ensembl: ", missing_ensembl)
    
    expDataRaw <- expDataRaw[!is.na(expDataRaw$Ensembl), ]
    expDataRaw <- expDataRaw %>% distinct(Ensembl, .keep_all = TRUE)
    
  } else if (to_id == "GeneName") {
    if (id_type == "GeneName") {
      message("Identifiers are already Gene Names. No conversion needed.")
      return(expDataRaw)
    }
    
    message("Converting Ensembl IDs to Gene Names...")
    
    attributes <- c("ensembl_gene_id", "hgnc_symbol")
    filters <- "ensembl_gene_id"
    
    geneR <- getBM(attributes = attributes,
                   filters = filters,
                   values = expDataRaw[, 1],
                   mart = mart)
    
    colnames(geneR) <- c("Ensembl", "Gene")
    expDataRaw <- left_join(expDataRaw, geneR, by = "Ensembl")
    
    missing_genes <- sum(is.na(expDataRaw$Gene))
    message("Number of Ensembl IDs without a corresponding Gene Name: ", missing_genes)
    
  } else {
    stop("Invalid conversion target. Choose 'Ensembl' or 'GeneName'.")
  }
  
  return(expDataRaw)
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
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           main = title)
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
  
  plot(sampleTrees, main = "Sample Clustering",
       xlab = "", sub = "", cex = cex_labels)
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
  
  pca_results <- PCA(mtx4PCA,
                     scale.unit = FALSE,
                     ncp = 20,
                     graph = FALSE)
  
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
    pca_plot <- pca_plot + aes(color = Population) +
      scale_color_manual(values = color_palette) +
      labs(color = "Population")
  }
  
  if (!is.na(sample_column)) {
    pca_plot <- pca_plot + aes(shape = Sample) + labs(shape = "Sample")
  }
  
  return(list(pca_data = pca_data, pca_plot = pca_plot, pca_results = pca_results))
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

#### 2) Load metadata #########################################################

message("1) Loading GEO metadata...")

geo_id <- "GSE42519"

tryCatch({
  gse <- getGEO(geo_id, GSEMatrix = TRUE)
}, error = function(e) {
  stop("Failed to download GEO data. Check the geo_id and internet connection.")
})

metadata <- pData(phenoData(gse[[1]]))

metadata <- metadata %>%
  mutate(
    LT_enrichment_fct = gsub(" [0-9]+$", "", title),
    
    LT_enrichment = case_when(
      LT_enrichment_fct %in% c("HSC") ~ 3,
      LT_enrichment_fct %in% c("MPP") ~ 2,
      LT_enrichment_fct %in% c("CMP", "GMP", "MEP") ~ 1,
      TRUE ~ NA_real_
    ),
    
    Immunophenotype = case_when(
      LT_enrichment_fct %in% c("HSC") ~ "Lin−/CD34+/CD38−/CD90+/CD45RA−",
      LT_enrichment_fct %in% c("MPP") ~ "Lin−/CD34+/CD38−/CD90−/CD45RA−",
      LT_enrichment_fct %in% c("CMP") ~ "Lin−/CD34+/CD38+/CD45RA−/CD123+",
      LT_enrichment_fct %in% c("GMP") ~ "Lin−/CD34+/CD38+/CD45RA+/CD123+",
      LT_enrichment_fct %in% c("MEP") ~ "Lin−/CD34+/CD38+/CD45RA−/CD123−",
      TRUE ~ NA_character_
    ),
    
    Enrichment = case_when(
      LT_enrichment_fct %in% c("HSC") ~ "Upper Enrichment",
      LT_enrichment_fct %in% c("MPP") ~ "Moderate Enrichment",
      LT_enrichment_fct %in% c("CMP", "GMP", "MEP") ~ "Lower Enrichment",
      TRUE ~ NA_character_
    ),
    
    Cell_Source = "Bone Marrow"
  )

#### 3) Read CEL files and normalize (BrainArray CDF) ##########################

message("2) Reading CEL files and running RMA (BrainArray Ensembl CDF)...")

data_folder <- "DATA/GSE42519"
cel_files <- list.celfiles(path = data_folder, full.names = TRUE)
length(cel_files)

# Load the custom BrainArray CDF environment (.cdf file must exist)
cdf_env <- make.cdf.env("DATA/HGU133Plus2_Hs_ENSG.cdf")

# IMPORTANT: cdfname MUST match the CDF name created by makecdfenv().
cdf_name <- "HGU133Plus2_Hs_ENSG"

normalized_data <- justRMA(
  filenames = cel_files,
  cdfname = cdf_name
)

exp_mtx <- exprs(normalized_data)

# Keep only Ensembl-mapped features (ENSG...)
ensg_rows <- grepl("^ENSG", rownames(exp_mtx))
exp_mtx <- exp_mtx[ensg_rows, , drop = FALSE]
rownames(exp_mtx) <- gsub("_at$", "", rownames(exp_mtx))

expData <- as.data.frame(exp_mtx)

#### 4) Subset to 16 samples (original behavior, but safe) ####################

# Original script: metadata <- metadata[c(1:16),]; expData <- expData[,c(1:16)]
# Here we do it only if both objects have at least 16 samples.
n_keep <- 16
if (nrow(metadata) >= n_keep && ncol(expData) >= n_keep) {
  metadata <- metadata[seq_len(n_keep), , drop = FALSE]
  expData  <- expData[, seq_len(n_keep), drop = FALSE]
} else {
  warning("Could not subset to 16 samples because metadata/expData has fewer than 16 samples.")
}

#### 5) Validate metadata vs expression #######################################

message("3) Validating metadata and expression alignment...")

if (all(rownames(metadata) %in% colnames(expData))) {
  message("The column names in expData match the row names in metadata. Proceeding with the analysis.")
} else {
  message("The column names in expData do not match the row names in metadata.")
  message("Tip: with affy::justRMA, column names often reflect CEL filenames.")
  message_text <- "Check if the row names of metadata match the column names of expData. Do they match? (yes/no): "
  response <- readline(prompt = message_text)
  
  if (tolower(response) == "yes") {
    colnames(expData) <- rownames(metadata)
    message("Columns of expData successfully renamed based on metadata row names.")
  } else if (tolower(response) == "no") {
    stop("Please reorganize metadata rows (or rename samples) to match expData columns, then rerun.")
  } else {
    message("Invalid response. Please answer 'yes' or 'no'. Proceeding without renaming columns.")
  }
}

#### 6) Save processed data ###################################################

write.csv(metadata, file = "DATA/metadata.csv")
write.csv(expData,  file = "DATA/normalized_expression_RMA.csv")

#### 7) QC plots ##############################################################

message("4) Generating QC plots...")

boxplot_plot_obj <- boxplot_plot(expData, metadata)

pca_results <- perform_pca(expData, metadata)
pca_plot_obj <- pca_results$pca_plot

heatmap_plot_obj <- heatmap_plot(expData, metadata)
sampleClustering_plot(expData, metadata)

# Export to PDF
pdf("RESULTS/_Quality_Control/3_QC.pdf", width = 8, height = 6)

print(boxplot_plot_obj)
print(pca_plot_obj)
grid.newpage()
print(heatmap_plot_obj)
sampleClustering_plot(expData, metadata)

dev.off()

#### 8) Interactive QC decision tree #########################################

repeat {
  outlier_resp <- tolower(readline("After the analysis, do you have any outlier samples? (Yes/No): "))
  
  if (outlier_resp == "yes") {
    message("Go back to the metadata and expData, remove the outliers, and run the analysis again.")
    break
  } else if (outlier_resp == "no") {
    repeat {
      batch_resp <- tolower(readline("Does your data have batch effects? (Yes/No): "))
      
      if (batch_resp == "yes") {
        message("Proceed to batch correction (e.g., ComBat) for microarray data.")
        break
      } else if (batch_resp == "no") {
        message("Quality control completed. Proceeding.")
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

write.csv(expData, file = "DATA/counts_final.csv")