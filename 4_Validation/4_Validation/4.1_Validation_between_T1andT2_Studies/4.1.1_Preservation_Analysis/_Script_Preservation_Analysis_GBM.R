###############################################################################
####  Script Title: Module Preservation Analysis — Violet Gene Set ###########
####  Description: Runs module preservation across HSC studies, generates 
####               heatmaps (Zsummary, Zconnectivity, Zdensity) for a 
####               synthetic “violet” module representing core LT-HSC genes. 
###############################################################################

#### 0. Libraries & Global Options ############################################

# 0.1) Set working directory
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022")

# 0.2) Load packages
suppressPackageStartupMessages({
  library(DESeq2)   # VST normalization
  library(dplyr)    # Data manipulation
  library(ggplot2)  # Visualization
  library(WGCNA)    # Coexpression network & preservation analysis
  library(ggrepel)  # Improved labeling in ggplot
  library(openxlsx)
})


#### 1. Load Studies ###########################################################

# 1.1) Define studies (counts, metadata, kME)
studies <- list(
  Anjos = list(
    expData  = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Anjos-Afonso_etal_2021/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Anjos-Afonso_etal_2021/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Anjos-Afonso_etal_2021/RESULTS/WGCNA/3_kmeTable.csv")
  ),
  Xie = list(
    expData  = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Xie_etal_2019/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Xie_etal_2019/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Xie_etal_2019/RESULTS/WGCNA/3_kmeTable.csv")
  ),
  Fares = list(
    expData  = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Fares_2017/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Fares_2017/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Fares_2017/RESULTS/WGCNA/3_kmeTable.csv")
  ),
  Sauvageau = list(
    expData  = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Sauvageau_2014/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Sauvageau_2014/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Sauvageau_2014/RESULTS/WGCNA/3_kmeTable.csv")
  ),
  XieCult = list(
    expData  = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Xie_2019/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Xie_2019/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_2_Studies/Xie_2019/RESULTS/WGCNA/3_kmeTable.csv")
  )
)

# 1.2) Load gene list of interest (expects 'Ensembl' & 'CountTotal'; keep CountTotal ∈ {4,5})
gene_list_tbl <- read.csv("2a_b_WGCNA_public_data_fresh_cultivated_HSC/Overlap_Uncult_Cult/Overlap/2_overlap_cult_kme_info_all.csv") %>%
  dplyr::filter(CountTotal %in% c(4, 5))

# 1.2.1) Extract Ensembl IDs and log a quick check
gene_list_ids <- unique(gene_list_tbl$Ensembl)
message("Gene list loaded — total genes: ", length(gene_list_ids))

# 1.3) VST normalization + kME de-duplication
for (study_name in names(studies)) {
  message("Processing VST: ", study_name)
  
  # 1.3.1) Raw counts & metadata
  raw_counts <- studies[[study_name]]$expData
  metadata   <- studies[[study_name]]$metadata
  
  # 1.3.2) Ensure rownames(metadata) align with columns(counts)
  if (is.null(rownames(metadata))) {
    rownames(metadata) <- colnames(raw_counts)
  }
  metadata <- metadata[colnames(raw_counts), , drop = FALSE]
  
  # 1.3.3) Build DESeq2 object (normalization-only)
  dds <- DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData   = metadata,
    design    = ~ 1
  )
  
  # 1.3.4) VST and transpose (samples as rows, genes as cols)
  vst_data   <- vst(dds, blind = TRUE)
  vst_matrix <- assay(vst_data)
  vst_df <- as.data.frame(t(vst_matrix))
  vst_df$SampleID <- rownames(vst_df)
  rownames(vst_df) <- vst_df$SampleID
  vst_df$SampleID <- NULL
  
  # 1.3.5) Unique kME per Ensembl
  filtered_kme <- studies[[study_name]]$kME %>%
    dplyr::distinct(Ensembl, .keep_all = TRUE)
  
  # 1.3.6) Store back
  studies[[study_name]]$vstData <- vst_df
  studies[[study_name]]$kME     <- filtered_kme
}


#### 2. Module Preservation (Fake "violet" Module) #############################

# 2.0) Config & output dir
nPerms          <- 200                 # Number of permutations for Z-scores
fakeModuleColor <- "violet"            # Label for the synthetic module
output_dir      <- "Overlap_Uncult_Cult/Preservation_Results_geneList"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2.0.1) Container to keep in-memory results
preservation_results <- list()

# 2.1) Run module preservation across all Ref → Test pairs
for (ref_name in names(studies)) {
  message("\n[Reference]: ", ref_name)
  
  # 2.1.1) Reference data and original module colors
  datRef <- studies[[ref_name]]$vstData
  kmeRef <- studies[[ref_name]]$kME
  colorRef <- kmeRef$moduleColor
  names(colorRef) <- kmeRef$Ensembl
  
  # 2.1.2) Paint only intersection genes as "violet"
  ref_present_ids <- intersect(gene_list_ids, names(colorRef))
  colorRef_fake   <- colorRef
  colorRef_fake[ref_present_ids] <- fakeModuleColor
  
  # 2.1.3) Prepare color list for modulePreservation (1 ref network)
  multiColor <- list(Reference = colorRef_fake)
  
  # 2.1.4) Loop over all other studies as test networks
  for (test_name in setdiff(names(studies), ref_name)) {
    message("  → Testing against: ", test_name)
    
    datTest <- studies[[test_name]]$vstData
    multiExpr <- list(
      Reference = list(data = datRef),
      Test      = list(data = datTest)
    )
    
    mp <- modulePreservation(
      multiExpr, multiColor,
      referenceNetworks = 1,
      nPermutations     = nPerms,
      randomSeed        = 1,
      quickCor          = 0,
      verbose           = 3
    )
    
    # 2.1.5) Extract metrics for this Ref → Test
    ref  <- 1
    test <- 2
    modColors <- rownames(mp$preservation$observed[[ref]][[test]])
    
    PresData <- cbind(
      mp$preservation$observed[[ref]][[test]][, 2], 
      mp$preservation$Z[[ref]][[test]][, 2:4],
      mp$preservation$log.pBonf[[ref]][[test]][, 2],
      mp$preservation$observed[[ref]][[test]][, 5],
      mp$preservation$Z[[ref]][[test]][, 5],
      mp$preservation$log.pBonf[[ref]][[test]][, 5],
      mp$preservation$Z[[ref]][[test]][, 8],
      mp$preservation$log.pBonf[[ref]][[test]][, 8],
      mp$preservation$observed[[ref]][[test]][, 11],
      mp$preservation$observed[[ref]][[test]][, 10],
      mp$preservation$Z[[ref]][[test]][, 10],
      mp$preservation$log.pBonf[[ref]][[test]][, 11],
      mp$preservation$observed[[ref]][[test]][, 13],
      mp$preservation$Z[[ref]][[test]][, 13],
      mp$preservation$log.pBonf[[ref]][[test]][, 14]
    )
    
    colnames(PresData) <- c(
      "medianRank",
      "Zsummary", "Zdensity", "Zconnectivity", "p.Zsummary",
      "propVarExpl", "Z.propVarExpl", "p.Z.propVarExpl",
      "Z.meanAdj", "p.Z.meanAdj",
      "cor.kME", "cor.kIM", "Z.cor.kIM", "p.Z.cor.kIM",
      "cor.cor", "Z.cor.cor", "p.cor.cor"
    )
    rownames(PresData) <- modColors
    
    # 2.1.6) Save results in memory and on disk
    obj_name <- paste0("PresData_", ref_name, "_vs_", test_name)
    preservation_results[[obj_name]] <- PresData
    
    output_path <- file.path(output_dir, paste0(obj_name, ".csv"))
    write.csv(PresData, output_path, quote = FALSE)
  }
}

message("\nDone: preservation results saved to: ", normalizePath(output_dir))

# 2.2) Stack ALL modules from ALL pairs into one big data frame
all_rows <- list()

for (obj_name in res_names) {
  pres_data <- fetch_fun(obj_name)          # matrix/data.frame with modules as rows
  
  # Parse "<Ref>_vs_<Test>"
  parts <- strsplit(sub("^PresData_", "", obj_name), "_vs_")[[1]]
  ref_study  <- parts[1]
  test_study <- parts[2]
  
  # Convert to data.frame and keep module names
  df <- as.data.frame(pres_data, stringsAsFactors = FALSE)
  df$ModuleColor <- rownames(pres_data)
  rownames(df) <- NULL
  
  # Add identifiers
  df$ReferenceStudy <- ref_study
  df$TestStudy      <- test_study
  
  # Reorder columns: IDs first, then metrics
  id_cols     <- c("ReferenceStudy", "TestStudy", "ModuleColor")
  metric_cols <- setdiff(colnames(df), id_cols)
  df <- df[, c(id_cols, metric_cols)]
  
  all_rows[[obj_name]] <- df
}

all_modules_summary <- dplyr::bind_rows(all_rows)

# 2.3) Filter only the "violet" module
violet_only <- subset(all_modules_summary, ModuleColor == "violet")

# 2.4) Save full and filtered preservation summaries
out_dir <- "Overlap_Uncult_Cult/Preservation_Results_geneList"
write.csv(
  all_modules_summary,
  file = file.path(output_dir, "all_modules_summary.csv"),
  row.names = FALSE, quote = FALSE
)

write.csv(
  violet_only,
  file = file.path(output_dir, "violet_only.csv"),
  row.names = FALSE, quote = FALSE
)

write.xlsx(
  violet_only,
  file = file.path(output_dir, "violet_only.xlsx"),
  rowNames = FALSE
)

message("Saved: all_modules_summary.csv and violet_only.csv in ", normalizePath(output_dir))

#### 3. Plots ##################################################################

# 3.1) Map display labels and set factor orders

lab_map <- c(
  "Anjos"     = "Anjos_2020_T1",
  "Xie"       = "Xie_2019_T1",
  "Fares"     = "Fares_2016_T2",
  "Sauvageau" = "Fares_2014_T2",
  "XieCult"   = "Xie_2019_T2"
)

ord <- c("Xie_2019_T2", "Fares_2016_T2", "Fares_2014_T2", "Xie_2019_T1", "Anjos_2020_T1")

df <- violet_only %>%
  mutate(
    Ref  = recode(ReferenceStudy, !!!lab_map),
    Test = recode(TestStudy,      !!!lab_map)
  ) %>%
  mutate(
    Ref  = factor(Ref,  levels = ord),
    Test = factor(Test, levels = ord)
  )

# 3.2) Heatmap function (design-aligned)
plot_metric <- function(df, metric, title_txt = metric) {
  pal <- colorRampPalette(c("blue", "red"))(100)
  
  # Grid for the whole Ref x Test matrix + diagonal flag
  grid_base <- expand.grid(
    Ref  = levels(df$Ref),
    Test = levels(df$Test),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      Ref  = factor(Ref,  levels = levels(df$Ref)),
      Test = factor(Test, levels = levels(df$Test)),
      is_diag = as.character(Ref) == as.character(Test)
    )
  
  # Metric values (off-diagonal only)
  pdat <- df |>
    dplyr::transmute(
      Ref, Test, value = .data[[metric]],
      is_diag = as.character(Ref) == as.character(Test)
    ) |>
    dplyr::filter(!is_diag)
  
  ggplot() +
    # Base layer: draw the diagonal as light gray + white borders
    geom_tile(
      data = grid_base,
      aes(x = Test, y = Ref),
      fill  = ifelse(grid_base$is_diag, "#F2F2F2", NA),
      color = "white", linewidth = 0.8
    ) +
    # Metric layer (off-diagonal)
    geom_tile(
      data = pdat,
      aes(x = Test, y = Ref, fill = value),
      color = "white", linewidth = 0.8
    ) +
    # Numeric labels
    geom_text(
      data = pdat,
      aes(x = Test, y = Ref, label = round(value, 2)),
      color = "white", size = 6 
    ) +
    # Color scale
    scale_fill_gradientn(colors = pal, name = metric) +
    # Labels & theme
    labs(x = NULL, y = NULL, title = title_txt) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = element_text(
        size = 13, color = "black", angle = 45,
        hjust = 1, vjust = 1, margin = margin(t = 8)
      ),
      axis.text.y = element_text(size = 13, color = "black"),
      plot.title  = element_text(hjust = 0.5, face = "bold")
    )
}

# 3.3) Generate plots
p_zsum  <- plot_metric(df, "Zsummary",      "Zsummary")
p_zconn <- plot_metric(df, "Zconnectivity", "Zconnectivity")
p_zdens <- plot_metric(df, "Zdensity",      "Zdensity")

print(p_zsum); print(p_zconn); print(p_zdens)

# 3.4) Save plots
out_dir <- "Artigo_GBM_VPO/Suplementar/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(filename = file.path(out_dir, "heatmap_Zsummary_315.png"),      plot = p_zsum,  width = 9, height = 7, dpi = 300)
ggsave(filename = file.path(out_dir, "heatmap_Zconnectivity_315.png"), plot = p_zconn, width = 9, height = 7, dpi = 300)
ggsave(filename = file.path(out_dir, "heatmap_Zdensity_315.png"),      plot = p_zdens, width = 9, height = 7, dpi = 300)
