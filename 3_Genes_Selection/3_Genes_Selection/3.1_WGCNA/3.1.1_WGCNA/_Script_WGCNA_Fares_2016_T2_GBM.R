###############################################################################
# Script:        WGCNA Fares 2014
# Author:        Gustavo Bueno Moreira
# Last update:   2026-02-16
#
# Purpose:
#   - Run WGCNA on variance-stabilized RNA-seq counts for the Fares (2014)
#     dataset (Type 2 study), generating network modules, kME tables, module
#     eigengenes, and module-trait correlation summaries and plots.
#
# Deliverables:
#   1) WGCNA network construction outputs (soft-threshold plots, dendrograms)
#   2) kME table (module membership + BH-adjusted p-values + gene annotations)
#   3) Module eigengenes table + module-trait correlation table + volcano plot
#
# Inputs:
#   - DATA/metadata.csv
#   - DATA/counts_Final.csv
#
# Outputs:
#   - DATA/counts_WGCNA_final.csv
#   - RESULTS/WGCNA/1_SoftThresholdTable.csv
#   - RESULTS/WGCNA/1_scale_independence.png
#   - RESULTS/WGCNA/1_mean_connectivity.png
#   - RESULTS/WGCNA/2_Cluster_dendrogram.png
#   - RESULTS/WGCNA/2_Cluster_dendogram_overview.png
#   - RESULTS/WGCNA/3_kmeTable.csv
#   - RESULTS/WGCNA/4_ME.csv
#   - RESULTS/WGCNA/4_moduleBarplots.pdf
#   - RESULTS/WGCNA/5_Correlation_Values.csv
#   - RESULTS/WGCNA/5_Module-Trait_relationship_padj_Spearman.png
#
# Notes:
#   - This script assumes it is executed from the study root set in setwd().
#   - The objects `colors` and `final_colors` are used for plotting and are
#     assumed to exist exactly as in the original workflow.
#   - The analysis logic, thresholds, and outputs are preserved exactly.
###############################################################################

#### 0) Setup #################################################################

#### 0.1) Setup: working directory ####
# IMPORTANT: keep setwd() exactly as in the original script
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Fares_2014")

#### 0.2) Options & constants ####
options(stringsAsFactors = FALSE)

#### 0.3) Packages ####
library(WGCNA)        # Weighted Gene Co-expression Network Analysis
library(multtest)     # Multiple testing procedures (used by WGCNA)
library(biomaRt)      # Access to Ensembl and other biological databases
library(coin)         # Permutation-based statistical tests (e.g., exact Spearman)
library(DESeq2)       # RNA-seq normalization and differential expression analysis
library(dplyr)        # Data manipulation (tidyverse)
library(tidyr)        # Data tidying (pivot, gather, etc.)
library(tibble)       # Enhanced data frames
library(ggplot2)      # Grammar of Graphics plotting
library(ggrepel)      # Improves label placement in ggplot2
library(patchwork)    # Combines multiple ggplot2 plots
library(ggnewscale)   # Allows multiple fill/colour scales in ggplot2
library(ggtext)       # Rich text formatting for ggplot2
library(grid)         # Low-level grid graphics system
library(gridExtra)    # Arrange multiple grid-based plots
library(gtable)       # Fine control over ggplot2 layouts
library(magick)       # Image editing and manipulation
library(svglite)      # SVG graphics output

#### 0.4) Functions (helpers) ####

# Calculate Spearman correlation with permutation-based p-value (coin)
calculate_spearman_with_permutation <- function(x, y, n_permutations = NULL, alternative = "greater") {
  spearman_test <- cor.test(x, y, method = "spearman")
  r_observed <- spearman_test$estimate
  
  if (is.na(r_observed)) {
    return(list(r = NA, p = NA))
  }
  
  statistic <- function(data) {
    cor(data[, 1], data[, 2], method = "spearman")
  }
  
  if (is.null(n_permutations)) {
    n <- length(x)
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
  permutation_result <- independence_test(
    x ~ y,
    data = data,
    teststat = "scalar",
    distribution = approximate(nresample = n_permutations),
    alternative = alternative
  )
  
  p_value <- pvalue(permutation_result)
  
  if (is.na(p_value)) {
    return(list(r = r_observed, p = NA))
  }
  
  return(list(r = r_observed, p = p_value))
}

# Analyze module-trait correlations for multiple traits (wide + long outputs)
analyze_module_trait_correlations_multi <- function(MEs, traits_df, alternative = "two.sided") {
  common_samples <- intersect(rownames(MEs), rownames(traits_df))
  MEs_aligned    <- MEs[common_samples, , drop = FALSE]
  traits_aligned <- traits_df[common_samples, , drop = FALSE]
  
  mod_names   <- colnames(MEs_aligned)
  trait_names <- colnames(traits_aligned)
  
  out_list <- vector("list", length(mod_names) * length(trait_names))
  k <- 1
  for (m in mod_names) {
    for (tname in trait_names) {
      res <- calculate_spearman_with_permutation(
        MEs_aligned[, m],
        traits_aligned[, tname],
        alternative = alternative
      )
      out_list[[k]] <- tibble::tibble(
        Module      = gsub("^ME", "", m),
        Trait       = tname,
        Correlation = as.numeric(res$r),
        P_value     = as.numeric(res$p)
      )
      k <- k + 1
    }
  }
  
  results_long <- dplyr::bind_rows(out_list) %>%
    dplyr::group_by(Trait) %>%
    dplyr::mutate(
      Adjusted_P_value = p.adjust(P_value, method = "BH"),
      Color            = paste0("ME", Module),
      SignifModule     = ifelse(is.finite(Adjusted_P_value) & Adjusted_P_value <= 0.05, Module, "ns")
    ) %>%
    dplyr::ungroup()
  
  results_wide <- results_long %>%
    dplyr::select(Module, Trait, Correlation, P_value, Adjusted_P_value, Color, SignifModule) %>%
    tidyr::pivot_wider(
      id_cols = Module,
      names_from = Trait,
      values_from = c(Correlation, P_value, Adjusted_P_value, Color, SignifModule),
      names_sep = "."
    )
  
  list(long = results_long, wide = results_wide)
}

#### 0.5) Load data ####
message("0.5) Loading input data...")

metadata = read.csv("DATA/metadata.csv", row.names = 1)
expData  = read.csv("DATA/counts_Final.csv", row.names = 1)

#### 1) Deliverable 1: WGCNA network construction ################################

#### 1.1) VST normalization for WGCNA input ####
message("1.1) Running VST normalization (DESeq2) and preparing WGCNA input matrix...")

dds <- DESeqDataSetFromMatrix(countData = expData, colData = metadata, design = ~ 1)
vsd = varianceStabilizingTransformation(dds)
ncvsd = assay(vsd)
ncvsd = as.data.frame(ncvsd)

# Save the WGCNA input matrix (genes x samples)
write.csv(ncvsd, file = "DATA/counts_WGCNA_final.csv")

dat <- t(ncvsd)           # WGCNA expects samples x genes
infoData <- colnames(dat) # gene identifiers (Ensembl)

#### 1.2) Soft-threshold selection ####
message("1.2) Estimating soft-threshold power...")

powers <- c(seq(from = 10, to = 30, by = 1))
sft <- pickSoftThreshold(dat, powerVector = powers, verbose = 5, blockSize = 20000, networkType = "signed")

sftTable = as.data.frame(sft)
write.csv(sftTable, "RESULTS/WGCNA/1_SoftThresholdTable.csv")

soft_threshold_df <- data.frame(
  Power = sft$fitIndices[, 1],
  SFT_R2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  MeanConnectivity = sft$fitIndices[, 5]
)

plot_a1 <- ggplot(soft_threshold_df, aes(x = Power, y = SFT_R2)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  geom_text(aes(label = Power), color = "black", size = 5) +
  scale_y_continuous(limits = c(0.5, 0.9), breaks = seq(0, 0.9, 0.1)) +
  labs(
    title = "Scale independence",
    x = "Soft Threshold (Power)",
    y = "Scale Free Topology Model Fit, signed R²"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.6),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank()
  )

plot_a2 <- ggplot(soft_threshold_df, aes(x = Power, y = MeanConnectivity)) +
  geom_text(aes(label = Power), color = "black", size = 4, vjust = -0.5) +
  labs(
    title = "Mean Connectivity",
    x = "Soft Threshold (Power)",
    y = "Mean Connectivity"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.6),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank()
  )

ggsave(
  filename = "RESULTS/WGCNA/1_scale_independence.png",
  plot = plot_a1,
  width = 7, height = 7, units = "in", dpi = 300
)

ggsave(
  filename = "RESULTS/WGCNA/1_mean_connectivity.png",
  plot = plot_a2,
  width = 7, height = 7, units = "in", dpi = 300
)

#### 1.3) Network construction (blockwiseModules) ####
message("1.3) Building coexpression network and detecting modules...")

net <- blockwiseModules(
  dat,
  power = 16,
  numericLabels = TRUE,
  networkType = "signed",
  minModuleSize = 150,
  deepsplit = 2,
  mergeCutHeight = 0.15,
  saveTOMs = FALSE,
  verbose = 5,
  minKMEtoStay = 0.5,
  nThreads = 24,
  maxBlockSize = 20000,
  checkMissingData = FALSE
)

modules = as.data.frame(table(net$colors))
colnames(modules) = c("Label", "N")
modules$Label = as.numeric(as.character(modules$Label))
modules$Color = labels2colors(modules$Label)
modules$Label = paste("M", modules$Label, sep = "")
moduleLabel = paste("M", net$colors, sep = "")
moduleColor = modules$Color[match(moduleLabel, modules$Label)]

message("1.3) Exporting module dendrogram plots...")

png("RESULTS/WGCNA/2_Cluster_dendrogram.png", width = 1200, height = 1200, res = 150)
par(oma = c(0, 0, 0, 0), mar = c(1, 5, 4, 1))

mergedColors <- labels2colors(net$colors)

plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  groupLabels = "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  cex.colorLabels = 1.1,
  cex.axis = 1.5
)

dev.off()

MEs <- moduleEigengenes(dat, colors = mergedColors)$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

png("RESULTS/WGCNA/2_Cluster_dendogram_overview.png", width = 2200, height = 2200, res = 150)
par(mar = c(5, 5, 4, 2))
plot(
  METree,
  main = "Clustering of Module Eigengenes",
  xlab = "", sub = "",
  cex.main = 2, cex.lab = 2, cex.axis = 2, cex = 2
)
abline(h = 0.15, col = "red", lty = 2)
dev.off()

#### 2) Deliverable 2: kME table and gene annotation ################################

#### 2.1) Calculating kME (module membership) and re-assigning genes ####
message("2.1) Calculating kME table and applying reassignment rules...")

nSam <- nrow(dat)
KMEs <- signedKME(dat, net$MEs, outputColumnName = "M")
kme <- data.frame(infoData[match(colnames(dat), infoData)], moduleColor, moduleLabel, KMEs)
colnames(kme)[1] = "Ensembl"

kmeInfoCols <- c(1:3)
kmedata <- kme[, -kmeInfoCols]
pvalBH <- kmedata
pvalBH[,] = NA

for (j in c(1:ncol(pvalBH))) {
  p = mt.rawp2adjp(corPvalueStudent(kmedata[, j], nSamples = nSam), proc = "BH")
  pvalBH[, j] = p$adjp[order(p$index), 2]
}

kme$newModule = "NA"
for (j in c(1:nrow(kmedata))) {
  if (j == 1) print("Working on genes 1:10000")
  if (j == 10000) print(paste("Working on genes 10000:", nrow(kmedata)))
  m = which(kmedata[j, ] == max(kmedata[j, ]))
  if ((pvalBH[j, m] < 0.05) & (kmedata[j, m] > 0.5)) kme$newModule[j] = as.character(colnames(kmedata)[m])
  else kme$newModule[j] = as.character(kme$moduleLabel[j])
}

modulesNew = as.data.frame(table(kme$newModule))
colnames(modules) = c("Label", "N")

for (j in c(1:nrow(kmedata))) {
  if (j == 1) print("Working on genes 1:10000")
  if (j == 10000) print(paste("Working on genes 10000:", nrow(kmedata)))
  m = which(colnames(kmedata[j, ]) == (kme$newModule[j]))
  if ((pvalBH[j, m] > 0.05) || (kmedata[j, m] < 0.5)) kme$newModule[j] = "NA"
}

kme$newModule[which(kme$newModule %in% "NA")] = "M0"
modulesNew2 <- as.data.frame(table(kme$newModule))
colnames(modules) = c("Label", "N")

kme$newColor <- kme$moduleColor[match(kme$newModule, kme$moduleLabel)]
kme$moduleLabel <- kme$newModule
kme$moduleColor = kme$newColor
kme <- kme[, -grep("newModule", colnames(kme))]
kme = kme[, -grep("newColor", colnames(kme))]

mod <- modules$Label[-1]
kmeTable <- kme[, kmeInfoCols]

for (j in c(1:length(mod))) {
  kmeTable = cbind(kmeTable, kmedata[, match(mod[j], colnames(kmedata))])
  colnames(kmeTable)[ncol(kmeTable)] = paste("kME", mod[j], sep = "_")
  
  kmeTable = cbind(kmeTable, pvalBH[, match(mod[j], colnames(pvalBH))])
  colnames(kmeTable)[ncol(kmeTable)] = paste("pvalBH", mod[j], sep = "_")
}

kmeTable[, (ncol(kmeTable) + 1)] = "Exp"
kmeTable = kmeTable[, c(1:3, ncol(kmeTable), 4:(ncol(kmeTable) - 1))]
colnames(kmeTable)[4] = "Expression"
colnames = colnames(kmeTable)
kmeTableComplete = kmeTable
colnames(kmeTableComplete)[which(colnames(kmeTableComplete) == "infoData.match.colnames.dat...infoData..")] <- "Ensembl"

lista <- kmeTableComplete[, 1]
geneR <- if (interactive()) {
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://grch37.ensembl.org", path = "/biomart/martservice")
  datasets <- listDatasets(mart)
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype"),
    filters = "ensembl_gene_id",
    values = lista,
    mart = mart
  )
}
colnames(geneR) <- c("Gene", "Ensembl", "Gene_type")

kmeTableFinal <- left_join(kmeTableComplete, geneR, by = "Ensembl")
kmeTableFinal <- kmeTableFinal[, c(1, (ncol(kmeTableFinal) - 1), (ncol(kmeTableFinal)), 2:(ncol(kmeTableFinal) - 2))]
write.csv(kmeTableFinal, "RESULTS/WGCNA/3_kmeTable.csv", row.names = FALSE)

#### 3) Deliverable 3: Eigengenes + module-trait correlations ################################

#### 3.1) Saving Module Eigengenes ####
message("3.1) Computing and exporting module eigengenes (ME)...")

MEs0 <- moduleEigengenes(dat, kmeTable$moduleColor)$eigengenes
MEs <- orderMEs(MEs0)
me <- MEs
write.csv(me, "RESULTS/WGCNA/4_ME.csv", row.names = TRUE)

#### 3.2) Module eigengene barplots ####
message("3.2) Exporting module eigengene barplots...")

legend_labels <- unique(metadata$LT_enrichment_fct)
legend_colors <- colors[!duplicated(metadata$LT_enrichment_fct)]

pdf("RESULTS/WGCNA/4_moduleBarplots.pdf", height = 5, width = 15)
mod = colnames(me)
for (m in mod) {
  j = match(m, colnames(me))
  col = m
  barplot(
    me[, j],
    xlab = "Samples",
    ylab = "ME",
    col = colors,
    main = col,
    names = me[, 1],
    cex.names = 0.5,
    axisnames = FALSE
  )
  legend(
    "topleft",
    legend = legend_labels,
    fill = legend_colors,
    bty = "n",
    cex = 0.7
  )
}
dev.off()

#### 3.3) Module-trait relationship (Spearman + permutation p-values) ####
message("3.3) Computing module-trait correlations and exporting tables...")

traits_df <- metadata[, c(41, 42)]
res_multi <- analyze_module_trait_correlations_multi(MEs = MEs, traits_df = traits_df, alternative = "two.sided")
moduleTrait_long <- res_multi$long

moduleTrait_long <- moduleTrait_long %>%
  dplyr::mutate(
    NegLog10FDR_raw  = -log10(Adjusted_P_value),
    NegLog10FDR_plot = ifelse(is.infinite(NegLog10FDR_raw), 4, NegLog10FDR_raw)
  )

module_ids <- kmeTableFinal %>%
  distinct(ModuleColor = moduleColor, ModuleID = moduleLabel)

gene_counts <- modulesNew2 %>%
  rename(ModuleID = Var1, GeneCount = Freq)

module_info <- module_ids %>%
  left_join(gene_counts, by = "ModuleID")

moduleTrait_long <- moduleTrait_long %>%
  left_join(module_info, by = c("Module" = "ModuleColor"))

moduleTrait_long$Study <- "Fares"

write.csv(
  moduleTrait_long,
  "RESULTS/WGCNA/5_Correlation_Values.csv",
  row.names = FALSE
)

#### 3.4) Volcano plot (module-trait: LT_enrichment) ####
message("3.4) Creating and saving volcano plot for LT_enrichment...")

df_lt <- moduleTrait_long %>%
  filter(Trait == "LT_enrichment")

plot <- ggplot(df_lt, aes(x = Correlation, y = NegLog10FDR_plot, fill = Module, color = Module)) +
  geom_jitter(width = 0.06, height = 0.06, size = 10, shape = 21, alpha = 0.6, stroke = 1) +
  scale_fill_manual(values = final_colors, guide = "none") +
  scale_color_manual(values = final_colors, guide = "none") +
  scale_x_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1),
    limits = c(-1, 1)
  ) +
  labs(
    title = "Correlation of Eigengenes with LT-HSC samples - Sauvageu",
    x = "Spearman Correlation",
    y = expression(-log[10] ~ "(FDR)")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) +
  geom_vline(xintercept = 0, col = "blue", lty = 2)

print(plot)

ggsave(
  filename = "RESULTS/WGCNA/5_Module-Trait_relationship_padj_Spearman.png",
  plot = plot,
  width = 9,
  height = 9,
  dpi = 300
)

#### Software environment #####################################################
sessionInfo()
