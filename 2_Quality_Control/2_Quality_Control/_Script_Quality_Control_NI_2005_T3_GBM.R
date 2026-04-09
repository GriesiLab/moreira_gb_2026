###############################################################################
# Script:        QC NI_2005
# Author:        Gustavo Bueno Moreira
# Last update:   2026-02-19
#
# Purpose:
#   - Perform quality control analysis on Affymetrix Human Genome U133A Array data,
#     including RMA normalization of CEL files and exploratory QC plots.
#
# Deliverables:
#   1) QC report (boxplot, PCA, heatmap, hierarchical clustering) exported to PDF
#   2) Normalized expression matrix (RMA) exported to CSV
#   3) Processed metadata exported to CSV
#
# Inputs:
#   - Raw CEL files: DATA/GSE3823/*.CEL (Affymetrix Human Genome U133A Array)
#   - GEO metadata downloaded via GEOquery (GSE3823)
#   - BrainArray custom CDF (local file):
#       - DATA/HGU133A_Hs_ENSG.cdf
#   - BrainArray custom CDF packages (installed from local tar.gz):
#       - DATA/hgu133ahsensgprobe_25.0.0.tar.gz
#       - DATA/hgu133ahsensgcdf_25.0.0.tar.gz
#       - DATA/pd.hgfocus.hs.ensg_25.0.0.tar.gz
#
# Outputs:
#   - DATA/metadata_all.csv
#   - DATA/normalized_expression_RMA_NonNormalized.csv
#   - RESULTS/_Quality_Control/3_QC_withoutOutliers.pdf
#
###############################################################################

#### 0) Setup #################################################################

# 0.1) Setup: working directory ------------------------------------------------
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/4a_Panel_validation_inSilico/Studies_QC/Microarray/NI_2005/")

# 0.2) Options & constants -----------------------------------------------------
options(stringsAsFactors = FALSE)

# 0.3) Install platform-specific annotation packages (BrainArray) --------------
# Note: These are custom packages for the Affymetrix Human Genome U133A Array
# Note: Acess http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/ensg.asp
install.packages("DATA/hgu133ahsensgprobe_25.0.0.tar.gz", repos = NULL, type = "source")
install.packages("DATA/hgu133ahsensgcdf_25.0.0.tar.gz", repos = NULL, type = "source")
install.packages("DATA/pd.hgfocus.hs.ensg_25.0.0.tar.gz", repos = NULL, type = "source")

# 0.4) Packages ----------------------------------------------------------------
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
library(data.table)    # High-performance data manipulation (fast fread/fwrite, etc.)
library(affy)          # Affymetrix microarray processing (CEL reading, RMA/justRMA)
library(makecdfenv)    # Build/load custom CDF environments (BrainArray CDF support)
library(stringr)       # String manipulation helpers (regex extraction, etc.)
library(hgu133ahsensgprobe)  # Probe annotation package (BrainArray / U133A)
library(hgu133ahsensgcdf)    # CDF package (BrainArray / U133A)
library(pd.hgfocus.hs.ensg)  # Platform design package (BrainArray)

#### 0.5) Functions (helpers) #################################################

# Function to convert identifiers to Ensembl ID
library(biomaRt)
library(dplyr)

convert_identifiers <- function(expDataRaw, to_id, mart) {
  # Detect identifier type
  id_type <- check_id_type(expDataRaw[, 1])
  message("Detected identifier type: ", id_type)
  
  # Define attributes and filters based on conversion direction
  if (to_id == "Ensembl") {
    if (id_type == "Ensembl") {
      message("Identifiers are already Ensembl. No conversion needed.")
      return(expDataRaw)
    }
    
    message("Converting identifiers to Ensembl ID...")
    
    # Define mapping attributes
    if (id_type == "Entrez") {
      attributes <- c("ensembl_gene_id", "entrezgene_id")
      filters <- "entrezgene_id"
    } else if (id_type == "GeneName") {
      attributes <- c("ensembl_gene_id", "hgnc_symbol")
      filters <- "hgnc_symbol"
    } else {
      stop("Unsupported identifier type for conversion to Ensembl.")
    }
    
    # Retrieve Ensembl IDs
    geneR <- getBM(attributes = attributes,
                   filters = filters,
                   values = expDataRaw[, 1],
                   mart = mart)
    
    colnames(geneR) <- c("Ensembl", "GeneID")
    geneR$GeneID <- as.character(geneR$GeneID)
    
    # Merge with original data
    expDataRaw <- expDataRaw %>%
      mutate(GeneID = as.character(expDataRaw[, 1])) %>%
      left_join(geneR, by = "GeneID")
    
    # Count and report missing matches
    missing_ensembl <- sum(is.na(expDataRaw$Ensembl))
    message("Number of genes without a match in Ensembl: ", missing_ensembl)
    
    # Remove genes without an Ensembl match
    expDataRaw <- expDataRaw[!is.na(expDataRaw$Ensembl), ]
    
    # Remove duplicate Ensembl IDs
    expDataRaw <- expDataRaw %>%
      distinct(Ensembl, .keep_all = TRUE)
    
  } else if (to_id == "GeneName") {
    if (id_type == "GeneName") {
      message("Identifiers are already Gene Names. No conversion needed.")
      return(expDataRaw)
    }
    
    message("Converting Ensembl IDs to Gene Names...")
    
    # Define attributes for Ensembl → Gene Name conversion
    attributes <- c("ensembl_gene_id", "hgnc_symbol")
    filters <- "ensembl_gene_id"
    
    # Retrieve gene names
    geneR <- getBM(attributes = attributes,
                   filters = filters,
                   values = expDataRaw[, 1],
                   mart = mart)
    
    colnames(geneR) <- c("Ensembl", "Gene")
    
    # Merge without removing duplicates
    expDataRaw <- left_join(expDataRaw, geneR, by = "Ensembl")
    
    # Count and report missing matches
    missing_genes <- sum(is.na(expDataRaw$Gene))
    message("Number of Ensembl IDs without a corresponding Gene Name: ", missing_genes)
    
  } else {
    stop("Invalid conversion target. Choose 'Ensembl' or 'GeneName'.")
  }
  
  return(expDataRaw)
}

# Function to generate a heatmap of sample distances
heatmap_plot <- function(expData, metadata, title = "Heatmap", label_column = NULL) {
  # If label_column is not provided, prompt the user
  if (is.null(label_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    label_column <- readline(prompt = "Enter the column name to use for row labels: ")
  }
  
  # Check if the column exists in metadata
  if (!label_column %in% colnames(metadata)) {
    stop(paste("Column", label_column, "does not exist in metadata."))
  }
  
  # Compute distance matrix
  sampleDists <- dist(t(expData))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  # Set row names based on the specified metadata column
  rownames(sampleDistMatrix) <- metadata[[label_column]]
  colnames(sampleDistMatrix) <- NULL
  
  # Define color palette
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  # Generate heatmap
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           main = title)
}

# Function for hierarchical clustering of samples
sampleClustering_plot <- function(expData, metadata, method = "average", cex_labels = 0.7, label_column = NULL) {
  # If label_column is not provided, prompt the user
  if (is.null(label_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    label_column <- readline(prompt = "Enter the column name to use for sample labels: ")
  }
  
  # Check if the column exists in metadata
  if (!label_column %in% colnames(metadata)) {
    stop(paste("Column", label_column, "does not exist in metadata."))
  }
  
  # Perform hierarchical clustering
  sampleTrees <- hclust(dist(t(expData)), method = method)
  
  # Set new sample labels based on metadata
  sampleTrees$labels <- metadata[[label_column]]
  
  # Configure plot layout (2 rows, 1 column)
  par(mfrow = c(2, 1))
  
  # Adjust plot margins
  par(mar = c(0, 4, 2, 0))
  
  # Plot clustering tree
  plot(sampleTrees, main = "Sample Clustering",
       xlab = "", sub = "", cex = cex_labels)
}

# Perform PCA analysis
perform_pca <- function(expData, metadata, population_column = NULL, sample_column = NULL) {
  # If population_column is not provided, ask the user
  if (is.null(population_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    population_column <- readline(prompt = "Enter the column name to define colors (Population) or type NA: ")
    if (toupper(population_column) == "NA") population_column <- NA
  }
  
  # If sample_column is not provided, ask the user
  if (is.null(sample_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    sample_column <- readline(prompt = "Enter the column name to define shapes (Sample) or type NA: ")
    if (toupper(sample_column) == "NA") sample_column <- NA
  }
  
  # If both are NA, stop execution
  if (is.na(population_column) & is.na(sample_column)) {
    stop("You must provide at least one variable for the plot!")
  }
  
  # Prepare the matrix for PCA: transpose data so samples are rows
  mtx4PCA <- t(expData)
  
  # Perform PCA analysis
  pca_results <- PCA(mtx4PCA,
                     scale.unit = FALSE,
                     ncp = 20,
                     graph = FALSE)
  
  # Extract coordinates of individuals (observations)
  pca_coordinates <- pca_results$ind$coord
  
  # Extract explained variance for each principal component (PC)
  variance_explained <- pca_results$eig
  variance_pc1 <- round(variance_explained[1, 2], 1)
  variance_pc2 <- round(variance_explained[2, 2], 1)
  
  # Create base DataFrame for ggplot
  pca_data <- as.data.frame(pca_coordinates)
  
  # Add population column if not NA
  if (!is.na(population_column)) {
    if (!population_column %in% colnames(metadata)) stop(paste("Column", population_column, "does not exist in metadata."))
    pca_data <- pca_data %>% mutate(Population = metadata[[population_column]])
  }
  
  # Add sample column if not NA
  if (!is.na(sample_column)) {
    if (!sample_column %in% colnames(metadata)) stop(paste("Column", sample_column, "does not exist in metadata."))
    pca_data <- pca_data %>% mutate(Sample = metadata[[sample_column]])
  }
  
  # Create PCA plot using ggplot2
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
  
  # Add colors only if Population is provided
  if (!is.na(population_column)) {
    num_colors <- length(unique(pca_data$Population))
    color_palette <- colorRampPalette(c("red", "blue"))(num_colors)
    pca_plot <- pca_plot + aes(color = Population) +
      scale_color_manual(values = color_palette) +
      labs(color = "Population")
  }
  
  # Add shapes only if Sample is provided
  if (!is.na(sample_column)) {
    pca_plot <- pca_plot + aes(shape = Sample) + labs(shape = "Sample")
  }
  
  return(list(pca_data = pca_data, pca_plot = pca_plot, pca_results = pca_results))
}

# Function to plot boxplot
boxplot_plot <- function(ncvsd, metadata, box_color = "steelblue", x_column = NULL, join_column = NULL) {
  # If x_column is not provided, display available columns in metadata and ask the user
  if (is.null(x_column)) {
    cat("Available columns in metadata:\n")
    print(colnames(metadata))
    x_column <- readline(prompt = "Enter the column name for the X-axis: ")
  }
  
  # If join_column is not provided, ask the user for it interactively
  if (is.null(join_column)) {
    join_column <- readline(prompt = "Enter the metadata column name that exactly matches the column names in expData:  ")
  }
  
  # Check if the selected columns exist in metadata
  if (!(x_column %in% colnames(metadata))) {
    stop(paste("Error: Column", x_column, "does not exist in metadata. Please try again."))
  }
  
  if (!(join_column %in% colnames(metadata))) {
    stop(paste("Error: Column", join_column, "does not exist in metadata. Please try again."))
  }
  
  # Transform the dataframe from wide to long format and perform a left join with metadata dynamically
  ncvsd_long <- ncvsd %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "expression") %>%
    left_join(metadata, by = c("sample" = join_column))
  
  # Create the boxplot
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

drawheatmap <- function(svdPV.m, title = NULL,
                        silent = FALSE, show_rownames = TRUE, show_colnames = TRUE, show_legend = TRUE) {
  myPalette         = c("darkred", "red", "orange", "pink", "white")
  breaks.v          = c(-10000, -10, -5, -2, log10(0.05), 0)
  rownames(svdPV.m) = paste("Comp.", 1:dim(svdPV.m)[1])
  tmp_plot          = pheatmap(log10(t(svdPV.m)),
                               cluster_rows  = FALSE,
                               cluster_cols  = FALSE,
                               breaks        = breaks.v,
                               color         = myPalette,
                               silent        = TRUE,
                               legend        = FALSE,
                               main          = paste("Interpretation of the Components -", title),
                               show_rownames = show_rownames,
                               show_colnames = show_colnames)
  
  if (show_legend) {
    leg = grid::legendGrob(
      c(expression("p < 1x" ~ 10^{-10}),
        expression("p < 1x" ~ 10^{-5}),
        "p < 0.01",
        "p < 0.05",
        "p > 0.05"),
      nrow = 5, pch = 15, gp = grid::gpar(fontsize = 10, col = myPalette)
    )
  } else {
    leg = NULL
  }
  
  final_plot = arrangeGrob(tmp_plot$gtable, leg, ncol = 2, widths = c(5, 1))
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

interpret_PC <- function(component_matrix, clinical_df, title = NULL,
                         silent = FALSE, show_rownames = TRUE, show_colnames = TRUE, show_legend = TRUE) {
  tmp <- apply(component_matrix, 2, function(x) {
    sapply(clinical_df, function(y) check_association_CHAMP(x, y))
  })
  drawheatmap(t(tmp), title = title,
              silent = silent, show_rownames = show_rownames, show_colnames = show_colnames, show_legend = show_legend)
}

get_pca_top_genes <- function(results_pca, pc = 1, top = 100, mart) {
  pca_contrib <- fviz_contrib(results_pca$pca_results, choice = "var",
                              axes = pc,
                              top = top)$data
  pca_contrib <- pca_contrib[order(-pca_contrib$contrib), ]
  colnames(pca_contrib)[1] <- "Ensembl"
  pca_contrib$Ensembl <- as.character(pca_contrib$Ensembl)
  pca_contrib <- convert_identifiers(pca_contrib, to_id = "GeneName", mart)
  return(pca_contrib)
}

#### 1) Deliverable 1: Load data & preprocessing ##############################

message("1) Loading GEO metadata and reading CEL files (RMA normalization)...")

##### 1.1 Read Metadata File ##################################################

geo_id <- "GSE3823"

tryCatch({
  gse <- getGEO(geo_id, GSEMatrix = TRUE)
}, error = function(e) {
  stop("Failed to download GEO data. Check the geo_id and internet connection.")
})

metadata <- pData(phenoData(gse[[1]]))

metadata <- metadata %>%
  mutate(
    LT_enrichment_fct = case_when(
      grepl("\\+/-", title) ~ paste(sub(" .*", "", title), "HSC"),
      grepl("\\+/\\+", title) ~ paste(sub(" .*", "", title), "Progenitor"),
      TRUE ~ NA_character_
    ),
    LT_enrichment = case_when(
      LT_enrichment_fct %in% c("FB HSC", "BM HSC", "CB HSC", "MPB HSC") ~ 1,
      LT_enrichment_fct %in% c("FB Progenitor", "BM Progenitor", "CB Progenitor", "MPB Progenitor") ~ 0,
      TRUE ~ NA_real_
    ),
    Immunophenotype = case_when(
      LT_enrichment_fct %in% c("FB HSC", "BM HSC", "CB HSC", "MPB HSC") ~ "CD34+/CD38-",
      LT_enrichment_fct %in% c("FB Progenitor", "BM Progenitor", "CB Progenitor", "MPB Progenitor") ~ "CD34+/CD38+",
      TRUE ~ NA_character_
    ),
    Enrichment = case_when(
      LT_enrichment_fct %in% c("FB HSC", "BM HSC", "CB HSC", "MPB HSC") ~ "Upper Enrichment",
      LT_enrichment_fct %in% c("FB Progenitor", "BM Progenitor", "CB Progenitor", "MPB Progenitor") ~ "Lower Enrichment",
      TRUE ~ NA_character_
    )
  )

##### 1.2 Read Expression Data (CEL -> RMA) ###################################

data_folder <- "DATA/GSE3823"

cel_files <- list.celfiles(path = data_folder, full.names = TRUE)
length(cel_files)

cdf_env <- make.cdf.env("DATA/HGU133A_Hs_ENSG.cdf")

normalized_data <- justRMA(
  filenames = cel_files,
  cdfname = "hgu133ahsensgcdf"
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

write.csv(metadata, file = "DATA/metadata_all.csv")
write.csv(expData, file = "DATA/normalized_expression_RMA_NonNormalized.csv")

#### 2) Deliverable 2: Quality control ########################################

message("2) Running QC (boxplot, PCA, heatmap, clustering)...")

##### 2.1 Boxplot #############################################################
boxplot_plot_obj <- boxplot_plot(expData, metadata)
print(boxplot_plot_obj)

##### 2.2 PCA #################################################################
pca_results <- perform_pca(expData, metadata)
pca_plot_obj <- pca_results$pca_plot
print(pca_results$pca_data)
print(pca_results$pca_plot)

##### 2.3 Heatmap #############################################################
heatmap_plot_obj <- heatmap_plot(expData, metadata)
print(heatmap_plot_obj)

##### 2.4 Sample clustering ###################################################
sampleClustering_plot(expData, metadata)

# Export QC panel to PDF -------------------------------------------------------
pdf("RESULTS/_Quality_Control/3_QC_withoutOutliers.pdf", width = 8, height = 6)

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
        message("Proceed to ##### 3.5 Run RUVseq ##### for batch effect removal.")
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
#### 3) Deliverable 3: Batch effect correction (ComBat) ########################

message("3) Running batch effect correction (ComBat) and re-generating QC plots...")

##### 3.1 Remove samples not suitable for ComBat ###############################
# Rationale:
#   - ComBat requires that each batch contains both biological states (here: HSC and Progenitor).
#   - We remove MPB samples because they do not contain HSC, and remove FB samples (per your rule).

metadata <- metadata[-c(1:9), ]
expData  <- expData[, -c(1:9)]

##### 3.2 Create batch-related columns ########################################
# Batch definition:
#   - CellSource: first token of LT_enrichment_fct (e.g., BM / CB / FB / MPB)
# Biological state:
#   - State: second token of LT_enrichment_fct (e.g., HSC / Progenitor)

# Extract the cell source (first word before the space)
metadata$CellSource <- sub(" .*", "", metadata$LT_enrichment_fct)

# Extract the biological state (word after the space)
metadata$State <- sub(".* ", "", metadata$LT_enrichment_fct)

##### 3.3 Convert to factors for ComBat compatibility ##########################
batch <- factor(metadata$CellSource)  # Batch = BM, CB, FB (after filtering)
state <- factor(metadata$State)       # Biological state = HSC, Progenitor

##### 3.4 Run ComBat ###########################################################
# The model matrix preserves the biological state during batch correction
mod <- model.matrix(~ state)

combat_exp <- ComBat(
  dat         = as.matrix(expData),
  batch       = batch,
  mod         = mod,
  par.prior   = TRUE,  # Use parametric empirical Bayes (recommended)
  prior.plots = TRUE   # Show diagnostic plots
)

# Replace expData with batch-corrected expression
expData <- as.data.frame(combat_exp)

##### 3.5 QC after batch correction ###########################################

# 3.5.1 Boxplot after ComBat ---------------------------------------------------
boxplot_plot_obj <- boxplot_plot(expData, metadata)
print(boxplot_plot_obj)

# 3.5.2 PCA after ComBat -------------------------------------------------------
pca_results  <- perform_pca(expData, metadata)
pca_plot_obj <- pca_results$pca_plot

# Optional: print PCA internal results
print(pca_results$pca_data)
print(pca_results$pca_plot)

# 3.5.3 Heatmap after ComBat ---------------------------------------------------
heatmap_plot_obj <- heatmap_plot(expData, metadata)
print(heatmap_plot_obj)

# 3.5.4 Sample clustering after ComBat ----------------------------------------
sampleClustering_plot(expData, metadata)

##### 3.6 Export QC panel to PDF ##############################################

pdf("RESULTS/_Quality_Control/3_QC_batch_Corrected.pdf", width = 8, height = 6)

print(boxplot_plot_obj)
print(pca_plot_obj)
grid.newpage()
print(heatmap_plot_obj)
sampleClustering_plot(expData, metadata)

dev.off()

##### 3.7 Save processed inputs ###############################################

write.csv(metadata, file = "DATA/metadata.csv")
write.csv(expData,  file = "DATA/normalized_expression_RMA.csv")
