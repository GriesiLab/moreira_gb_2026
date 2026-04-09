###############################################################################
# Script:        QC_Prashad_2014.R
# Author:        Gustavo Bueno Moreira
# Last update:   2026-02-19
#
# Purpose:
#   - Perform quality control analysis on Affymetrix Human Genome U133 Plus 2.0
#     Array data, including RMA normalization of CEL files and exploratory QC plots.
#
# Deliverables:
#   1) QC report (boxplot, PCA, heatmap, hierarchical clustering) exported to PDF
#   2) Normalized expression matrix (RMA) exported to CSV
#   3) Processed metadata exported to CSV
#
# Inputs:
#   - Raw CEL files: DATA/GSE54316/*.CEL (Affymetrix Human Genome U133 Plus 2.0 Array)
#   - GEO metadata downloaded via GEOquery (GSE54316)
#   - BrainArray custom CDF environment:
#       - DATA/HGU133Plus2_Hs_ENSG.cdf
#
# Outputs:
#   - DATA/metadata.csv
#   - DATA/normalized_expression_RMA.csv
#   - RESULTS/_Quality_Control/3_QC.pdf
#
###############################################################################

#### 0) Setup #################################################################

# 0.1) Setup: working directory ------------------------------------------------
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/4a_Panel_validation_inSilico/Studies_QC/Microarray/Prashad_2014")

# 0.2) Options & constants -----------------------------------------------------
options(stringsAsFactors = FALSE)

# 0.3) Packages ----------------------------------------------------------------
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
library(MOFA2)         # Unsupervised analysis of multi-omics data (Modeling Omics Factor Analysis)
library(DT)            # Interactive tables for displaying data, commonly in Shiny apps
library(psych)         # Tools for psychological research, including correlation and factor analysis
library(NMF)           # Non-negative Matrix Factorization for data analysis (e.g., genomics, image analysis)
library(dplyr)         # Data manipulation package with verbs for filtering, mutating, etc.
library(biomaRt)       # Interface to biological databases (e.g., Ensembl) for genomic queries
library(UpSetR)        # Visualize intersecting sets and their properties
library(data.table)    # High-performance data manipulation
library(affy)          # Affymetrix preprocessing (CEL reading, RMA via justRMA)
library(makecdfenv)    # Create custom CDF environments from .cdf files
library(stringr)       # String manipulation utilities
library(pd.hgu133plus2.hs.ensg) # Platform design info (BrainArray / ENSG mapping)
library(hgu133plus2hsensgprobe) # Probe annotations (BrainArray / ENSG mapping)

#### 0.4) Functions (helpers) #################################################

# NOTE: This script assumes you have a function called `check_id_type()`.
# If not, add your own identifier detector or remove conversion steps that depend on it.

convert_identifiers <- function(expDataRaw, to_id, mart) {
  # Detect identifier type
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
    
    geneR <- getBM(
      attributes = attributes,
      filters = filters,
      values = expDataRaw[, 1],
      mart = mart
    )
    
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
  
  plot(
    sampleTrees,
    main = "Sample Clustering",
    xlab = "",
    sub = "",
    cex = cex_labels
  )
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

#### 1) Deliverable 1: Load data & preprocessing ##############################

message("1) Loading GEO metadata and reading CEL files (RMA normalization)...")

##### 1.1 Read Metadata File ##################################################

geo_id <- "GSE54316"

tryCatch({
  gse <- getGEO(geo_id, GSEMatrix = TRUE)
}, error = function(e) {
  stop("Failed to download GEO data. Check the geo_id and internet connection.")
})

metadata <- pData(phenoData(gse[[1]]))

metadata <- metadata %>%
  mutate(
    LT_enrichment_fct = case_when(
      row_number() %in% 1:3   ~ "HSC",
      row_number() %in% 4:5   ~ "Progenitor CD90-",
      row_number() %in% 6:7   ~ "Committed Progenitor",
      row_number() %in% 8:10  ~ "HSC GPI-80+",
      row_number() %in% 11:13 ~ "Progenitor CD90+",
      TRUE ~ NA_character_
    ),
    LT_enrichment = case_when(
      LT_enrichment_fct %in% c("HSC GPI-80+") ~ 5,
      LT_enrichment_fct %in% c("HSC") ~ 4,
      LT_enrichment_fct %in% c("Progenitor CD90+") ~ 3,
      LT_enrichment_fct %in% c("Progenitor CD90-") ~ 2,
      LT_enrichment_fct %in% c("Committed Progenitor") ~ 1,
      TRUE ~ NA_real_
    ),
    Imunophenotype = case_when(
      LT_enrichment_fct %in% c("HSC GPI-80+") ~ "CD34+/CD38-/CD90+/GPI-80+",
      LT_enrichment_fct %in% c("HSC") ~ "CD34+/CD38-/CD90+",
      LT_enrichment_fct %in% c("Progenitor CD90+") ~ "CD34+/CD38-/CD90+/GPI-80-",
      LT_enrichment_fct %in% c("Progenitor CD90-") ~ "CD34+/CD38-/CD90-",
      LT_enrichment_fct %in% c("Committed Progenitor") ~ "CD34+/CD38+/CD90-",
      TRUE ~ NA_character_
    ),
    Enrichment = case_when(
      LT_enrichment_fct %in% c("HSC GPI-80+") ~ "Upper Enrichment",
      LT_enrichment_fct %in% c("HSC") ~ "High Enrichment",
      LT_enrichment_fct %in% c("Progenitor CD90+") ~ "Moderate Enrichment",
      LT_enrichment_fct %in% c("Progenitor CD90-") ~ "Low Enrichment",
      LT_enrichment_fct %in% c("Committed Progenitor") ~ "Lower Enrichment",
      TRUE ~ NA_character_
    ),
    Cell_Source = "Fetal Liver"
  )

##### 1.2 Read Expression Data (CEL -> RMA) ###################################

data_folder <- "DATA/GSE54316"
cel_files <- list.celfiles(path = data_folder, full.names = TRUE)
length(cel_files)

cdf_env <- make.cdf.env("DATA/HGU133Plus2_Hs_ENSG.cdf")

normalized_data <- justRMA(
  filenames = cel_files,
  cdfname = "HGU133Plus2_Hs_ENSG"
)

print(summary(normalized_data))
head(exprs(normalized_data))

raw_matrix <- exprs(normalized_data)

ensg_rows <- grepl("^ENSG", rownames(raw_matrix))
filtered_matrix <- raw_matrix[ensg_rows, ]
rownames(filtered_matrix) <- gsub("_at$", "", rownames(filtered_matrix))

expData <- as.data.frame(filtered_matrix)

##### 1.3 Validate Metadata and Expression Data ################################

if (all(rownames(metadata) %in% colnames(expData))) {
  message("The column names in expData match the row names in metadata. Proceeding with the analysis.")
} else {
  message("The column names in expData do not match the row names in metadata.")
  message_text <- "Check if the row names of metadata match the column names of expData. Do they match? (yes/no)"
  response <- readline(prompt = message_text)
  
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

##### 1.4 Save processed inputs ###############################################

expData1 <- expData
write.csv(metadata, file = "DATA/metadata.csv")
write.csv(expData, file = "DATA/normalized_expression_RMA.csv")

#### 2) Deliverable 2: Quality control ########################################

message("2) Running QC (boxplot, PCA, heatmap, clustering)...")

##### 2.1 Boxplot #############################################################
boxplot_plot_obj <- boxplot_plot(expData, metadata)
print(boxplot_plot_obj)

##### 2.2 PCA #################################################################
pca_results <- perform_pca(expData, metadata)
pca_plot_obj <- pca_results$pca_plot
print(pca_results$pca_data)
print(pca_plot_obj)

##### 2.3 Heatmap #############################################################
heatmap_plot_obj <- heatmap_plot(expData, metadata)
print(heatmap_plot_obj)

##### 2.4 Sample clustering ###################################################
sampleClustering_plot(expData, metadata)

# Export QC panel to PDF -------------------------------------------------------
pdf("RESULTS/_Quality_Control/3_QC.pdf", width = 8, height = 6)

print(boxplot_plot_obj)
print(pca_plot_obj)
grid.newpage()
print(heatmap_plot_obj)
sampleClustering_plot(expData, metadata)

dev.off()

# Interactive QC decision tree (outliers / batch effects) ----------------------
repeat {
  outlier_resp <- tolower(readline("After the analysis, do you have any outlier samples? (Yes/No): "))
  
  if (outlier_resp == "yes") {
    message("Go back to the metadata and expData, remove the outliers, and run the analysis again.")
    break
  } else if (outlier_resp == "no") {
    repeat {
      batch_resp <- tolower(readline("Does your data have batch effects? (Yes/No): "))
      
      if (batch_resp == "yes") {
        message("Proceed to Run COMBAT for batch effect removal.")
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