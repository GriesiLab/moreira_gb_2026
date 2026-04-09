# Module Preservation Analysis – All Modules Across All Studies
#
# Description:
# This script performs a module preservation analysis by comparing the co-expression 
# modules identified by WGCNA in each study against all other studies.
#
# Input: raw expression data (counts), sample metadata, module membership (kME)
# Output: CSV tables containing module preservation statistics (Zsummary, density, etc.)


### 0. Load required libraries ----
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC")

library(DESeq2)     # For VST normalization
library(dplyr)      # Data manipulation
library(ggplot2)    # Plotting (optional)
library(WGCNA)      # Co-expression network and preservation analysis
library(ggrepel)    # Enhanced ggplot labeling

### 1. Load expression data, metadata, and kME tables ----
# Each study is stored as a named list inside the "studies" object

studies <- list(
  Anjos = list(
    expData  = read.csv("Type_1_Studies/Anjos-Afonso_etal_2021/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("Type_1_Studies/Anjos-Afonso_etal_2021/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("Type_1_Studies/Anjos-Afonso_etal_2021/RESULTS/WGCNA/3_kmeTable.csv")
  ), 
  Xie = list(
    expData  = read.csv("Type_1_Studies/Xie_etal_2019/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("Type_1_Studies/Xie_etal_2019/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("Type_1_Studies/Xie_etal_2019/RESULTS/WGCNA/3_kmeTable.csv")
  ),
  Fares = list(
    expData  = read.csv("Type_2_Studies/Fares_2017/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("Type_2_Studies/Fares_2017/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("Type_2_Studies/Fares_2017/RESULTS/WGCNA/3_kmeTable.csv")
  ),
  Sauvageau = list(
    expData  = read.csv("Type_2_Studies/Sauvageau_2014/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("Type_2_Studies/Sauvageau_2014/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("Type_2_Studies/Sauvageau_2014/RESULTS/WGCNA/3_kmeTable.csv")
  ),
  XieCult = list(
    expData  = read.csv("Type_2_Studies/Xie_2019/DATA/counts_Final.csv", row.names = 1),
    metadata = read.csv("Type_2_Studies/Xie_2019/DATA/metadata.csv", row.names = 1),
    kME      = read.csv("Type_2_Studies/Xie_2019/RESULTS/WGCNA/3_kmeTable.csv")
  )
)

### 2. Normalize expression data using VST (Variance Stabilizing Transformation) ----
# Also filters kME for unique Ensembl gene IDs

for (study_name in names(studies)) {
  message("Processing study: ", study_name)
  
  # Extract expression matrix and metadata
  raw_counts <- studies[[study_name]]$expData
  metadata   <- studies[[study_name]]$metadata
  
  # Ensure metadata rownames match column names of counts
  if (is.null(rownames(metadata))) {
    rownames(metadata) <- colnames(raw_counts)
  }
  metadata <- metadata[colnames(raw_counts), , drop = FALSE]
  
  # Create DESeq2 object (design is ~1 because no condition is being tested)
  dds <- DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData   = metadata,
    design    = ~ 1
  )
  
  # Perform VST normalization (blind = TRUE to ignore experimental design)
  vst_data <- vst(dds, blind = TRUE)
  vst_matrix <- assay(vst_data)
  
  # Transpose and convert to data frame (samples as rows)
  vst_df <- as.data.frame(t(vst_matrix))
  vst_df$SampleID <- rownames(vst_df)
  rownames(vst_df) <- vst_df$SampleID
  vst_df$SampleID <- NULL
  
  # Filter duplicated Ensembl entries (keep only the first occurrence)
  filtered_kme <- studies[[study_name]]$kME %>%
    distinct(Ensembl, .keep_all = TRUE)
  
  # Save processed data back into the list
  studies[[study_name]]$vstData <- vst_df
  studies[[study_name]]$kME     <- filtered_kme
}

### 3. Perform module preservation analysis across all studies ----
# For each study, use its modules as reference and evaluate preservation in all other studies

# Set parameters
nPerms <- 200
output_dir <- "Preservation_Analysis/Preservation_Results_AllModules"
dir.create(output_dir, showWarnings = FALSE)

# Initialize a list to store preservation results
preservation_results <- list()

# Loop through each study as the reference network
for (ref_name in names(studies)) {
  message("\n[Reference]: ", ref_name)
  
  # Load reference expression data and module information
  datRef <- studies[[ref_name]]$vstData
  kmeRef <- studies[[ref_name]]$kME
  
  # Define the real module colors assigned by WGCNA
  colorRef <- kmeRef$moduleColor
  names(colorRef) <- kmeRef$Ensembl
  multiColor <- list(Reference = colorRef)
  
  # Compare the reference network to each other study as test network
  for (test_name in setdiff(names(studies), ref_name)) {
    message("  → Testing against: ", test_name)
    
    # Load test expression data
    datTest <- studies[[test_name]]$vstData
    
    # Define the multiExpr input (no need to manually intersect genes)
    multiExpr <- list(
      Reference = list(data = datRef),
      Test      = list(data = datTest)
    )
    
    # Run module preservation analysis using real modules
    mp <- modulePreservation(
      multiExpr,
      multiColor,
      referenceNetworks = 1,
      nPermutations = nPerms,
      randomSeed = 1,
      quickCor = 0,
      verbose = 3
    )
    
    # Extract statistics
    ref <- 1
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
    
    colnames(PresData) <- c("medianRank", "Zsummary", "Zdensity", "Zconnectivity", "p.Zsummary",
                            "propVarExpl", "Z.propVarExpl", "p.Z.propVarExpl",
                            "Z.meanAdj", "p.Z.meanAdj",
                            "cor.kME", "cor.kIM", "Z.cor.kIM", "p.Z.cor.kIM",
                            "cor.cor", "Z.cor.cor", "p.cor.cor")
    
    rownames(PresData) <- modColors
    
    # Save results to memory and disk
    obj_name <- paste0("PresData_", ref_name, "_vs_", test_name)
    preservation_results[[obj_name]] <- PresData
    assign(obj_name, PresData, envir = .GlobalEnv)
    
    output_path <- file.path(output_dir, paste0(obj_name, ".csv"))
    write.csv(PresData, output_path, quote = FALSE)
  }
}

# List all preservation result objects in the environment
pres_obj_names <- ls(pattern = "^PresData_")

# Initialize an empty list to store each comparison's data frame
pres_all_df_list <- list()

# Loop through all preservation result objects
for (obj_name in names(preservation_results)) {
  
  # Extract the reference and test names from the object name
  parts <- strsplit(obj_name, "_vs_")[[1]]
  ref_study <- gsub("PresData_", "", parts[1])
  test_study <- parts[2]
  
  # Extract the data
  df <- preservation_results[[obj_name]]
  
  # Add metadata columns
  df$ReferenceStudy <- ref_study
  df$TestStudy <- test_study
  df$ModuleColor <- rownames(df)
  
  # Move metadata columns to the front (optional)
  df <- df %>% dplyr::relocate(ReferenceStudy, TestStudy, ModuleColor)
  
  # Add to list
  pres_all_df_list[[obj_name]] <- df
}

# Combine all into one big data frame
pres_all_df <- do.call(rbind, pres_all_df_list)
rownames(pres_all_df) <- NULL  # Optional: clean rownames


# 1. Initialize an empty numeric vector
module_sizes <- vector("numeric", length = nrow(pres_all_df))

# 2. For each row in pres_all_df, retrieve the module size from the reference study
for (i in seq_len(nrow(pres_all_df))) {
  
  ref_study <- pres_all_df$ReferenceStudy[i]
  module_color <- pres_all_df$ModuleColor[i]
  
  # Count how many genes have that color in the reference study
  size <- sum(studies[[ref_study]]$kME$moduleColor == module_color)
  
  module_sizes[i] <- size
}

# 3. Add the result to the data frame
pres_all_df$moduleSize <- module_sizes


# Save
write.csv(pres_all_df, file = file.path(output_dir, "All_Preservation_Results_Combined.csv"), row.names = FALSE)

### 4. Plot Zsummary values for all module comparisons (optional) ----
# This section generates simple barplots of Zsummary statistics for each reference vs test pair
# 0. Load required libraries


# =========================
# PANELS: preservation plots (5x5 grid com diagonal cinza)
# =========================

library(ggplot2)
library(dplyr)

# 0) Ordem fixa dos estudos (linhas = TestStudy, colunas = ReferenceStudy)
study_order <- c("Anjos", "Xie", "Fares", "Sauvageau", "XieCult")

# 1) Data pronto para facet: manter todos os pares, mas remover ref==test dos pontos
pres_plot_df <- pres_all_df %>%
  mutate(
    ReferenceStudy = factor(ReferenceStudy, levels = study_order),
    TestStudy      = factor(TestStudy,      levels = study_order)
  )

# 2) Diagonal cinza (uma célula por par ref==test)
gray_df <- expand.grid(
  ReferenceStudy = study_order,
  TestStudy      = study_order
) %>%
  dplyr::filter(ReferenceStudy == TestStudy)

# 3) Paleta: usar a própria cor dos módulos (assumindo nomes/hex válidos em ModuleColor)
module_levels  <- sort(unique(pres_plot_df$ModuleColor))
module_palette <- setNames(module_levels, module_levels)

# 4) Jitter: X (moduleSize) moderado; Y (métrica) bem pequeno
my_jitter <- position_jitter(width = 15, height = 0.05)

# 5) Tema base
base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.text   = element_text(size = 10, face = "bold"),
    axis.title   = element_text(size = 11),
    axis.text    = element_text(size = 9),
    legend.position = "none"
  )

# 6) Função para um painel (uma métrica)
make_preservation_panel <- function(metric, panel_title, filename,
                                    x_lim = c(20, 2000), y_lim = c(-2, 10),
                                    hlines = NULL, outdir = output_dir) {
  
  df_points <- pres_plot_df %>% dplyr::filter(ReferenceStudy != TestStudy)
  
  p <- ggplot(df_points, aes(x = moduleSize, y = .data[[metric]], color = ModuleColor)) +
    geom_point(size = 4, alpha = 0.5, stroke = 0.8, shape = 16, position = my_jitter) +
    facet_grid(TestStudy ~ ReferenceStudy, switch = "both") +
    scale_color_manual(values = module_palette) +
    scale_x_continuous(limits = x_lim) +
    scale_y_continuous(limits = y_lim) +
    labs(title = panel_title, x = "Module size", y = panel_title) +
    base_theme +
    theme(strip.placement = "outside")  # << nomes dos estudos fora, depois dos ticks
  
  # linhas de corte (só nos painéis não-diagonais)
  if (!is.null(hlines) && length(hlines) > 0) {
    for (h in hlines) {
      p <- p + geom_hline(yintercept = h, linetype = "dashed",
                          color = ifelse(h == 10, "darkgreen", "red"),
                          linewidth = 0.7)
    }
  }
  
  # AGORA desenha a diagonal cinza por cima (cobre as hlines nos painéis diagonais)
  p <- p + geom_rect(
    data = gray_df,
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    fill = "lightgray", color = "black", linewidth = 0.5
  )
  
  ggsave(file.path(outdir, paste0(filename, ".png")), p, width = 8.5, height = 11, units = "in")
}

# =========================
# 4 PDFs (A4 retrato)
# =========================

# 1) Zsummary: X = moduleSize (20-2000), Y = Zsummary (-2 a 10), linhas y=2 (vermelha) e y=10 (verde)
make_preservation_panel(
  metric = "Zsummary",
  panel_title = "Preservation Zsummary",
  filename = "panel_Zsummary",
  x_lim = c(20, 2000),
  y_lim = c(-2, 15),
  hlines = c(2, 10),
  outdir = output_dir
)

# 2) Zdensity
make_preservation_panel(
  metric = "Zdensity",
  panel_title = "Preservation Zdensity",
  filename = "panel_Zdensity",
  x_lim = c(20, 2000),
  y_lim = c(-2, 15),
  hlines = c(2, 10),
  outdir = output_dir
)

# 3) Zconnectivity
make_preservation_panel(
  metric = "Zconnectivity",
  panel_title = "Preservation Zconnectivity",
  filename = "panel_Zconnectivity",
  x_lim = c(20, 2000),
  y_lim = c(-2, 15),
  hlines = c(2, 10),
  outdir = output_dir
)

# 4) medianRank: X = moduleSize (20-2000), Y = medianRank (0 a 1), sem linhas de corte
make_preservation_panel(
  metric = "medianRank",
  panel_title = "Preservation median rank",
  filename = "panel_medianRank",
  x_lim = c(20, 2000),
  y_lim = c(0, 1),
  hlines = NULL,
  outdir = output_dir
)

# 1) Mapa de cores permitidas (em minúsculas para padronizar)
allowed_colors_by_ref <- list(
  "Anjos"     = c("blue", "turquoise", "tan"),
  "Xie"       = c("blue", "turquoise", "green", "darkturquoise"),
  "Fares"     = c("brown"),
  "Sauvageau" = c("blue", "brown", "salmon"),
  "XieCult"   = c("cyan", "green", "greenyellow", "grey60", "purple")
)

# 2) Expandir para uma tabela de pares (ReferenceStudy x ModuleColorLower)
allowed_pairs <- lapply(names(allowed_colors_by_ref), function(ref) {
  data.frame(
    ReferenceStudy = ref,
    ModuleColorLower = tolower(allowed_colors_by_ref[[ref]]),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# 3) Filtrar pres_plot_df usando os pares permitidos
pres_plot_df_filt <- pres_plot_df %>%
  mutate(ModuleColorLower = tolower(ModuleColor)) %>%
  inner_join(allowed_pairs, by = c("ReferenceStudy", "ModuleColorLower")) %>%
  select(-ModuleColorLower)

# 4) Atualizar paleta apenas com as cores remanescentes
module_levels_filt  <- sort(unique(pres_plot_df_filt$ModuleColor))
module_palette      <- setNames(module_levels_filt, module_levels_filt)

# 5) Substituir o data frame que a função usa (mantendo a função inalterada)
pres_plot_df <- pres_plot_df_filt

# Ordem fixa dos estudos em ambos os eixos
study_order <- c("Anjos", "Xie", "Fares", "Sauvageau", "XieCult")

# Após o filtro e antes de plotar:
pres_plot_df <- pres_plot_df %>%
  dplyr::mutate(
    ReferenceStudy = factor(ReferenceStudy, levels = study_order),
    TestStudy      = factor(TestStudy,      levels = study_order)
  )

# =========================
# Gerar novamente os 4 PDFs (A4 retrato) com os filtros aplicados
# (usa a sua função make_preservation_panel já definida)
# =========================

# 1) Zsummary
make_preservation_panel(
  metric = "Zsummary",
  panel_title = "Preservation Zsummary",
  filename = "panel_Zsummary_FILTERED",
  x_lim = c(20, 2000),
  y_lim = c(-2, 25),
  hlines = c(2, 10),
  outdir = output_dir
)

# 2) Zdensity
make_preservation_panel(
  metric = "Zdensity",
  panel_title = "Preservation Zdensity",
  filename = "panel_Zdensity_FILTERED",
  x_lim = c(20, 2000),
  y_lim = c(-2, 25),
  hlines = c(2, 10),
  outdir = output_dir
)

# 3) Zconnectivity
make_preservation_panel(
  metric = "Zconnectivity",
  panel_title = "Preservation Zconnectivity",
  filename = "panel_Zconnectivity_FILTERED",
  x_lim = c(20, 2000),
  y_lim = c(-2, 25),
  hlines = c(2, 10),
  outdir = output_dir
)

# 4) medianRank (sem linhas de corte)
make_preservation_panel(
  metric = "medianRank",
  panel_title = "Preservation median rank",
  filename = "panel_medianRank_FILTERED",
  x_lim = c(20, 2000),
  y_lim = c(0, 1),
  hlines = NULL,
  outdir = output_dir
)


# Install once (if needed):
# install.packages("writexl")

library(dplyr)
install.packages("writexl")
library(writexl)

# Ensure output dir exists (you já criou acima como `output_dir`)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1) Safety: make sure these objects exist
stopifnot(exists("pres_all_df"), exists("pres_plot_df"))

# 2) Optional: if you criou a versão filtrada (pres_plot_df_filt), adicionamos como sheet
sheets_list <- list(
  pres_all_df   = pres_all_df,         # tabela completa combinada
  pres_plot_df  = pres_plot_df         # dados usados nos painéis (sem diagonal)
)

if (exists("pres_plot_df_filt")) {
  sheets_list$pres_plot_df_filtered <- pres_plot_df_filt
}

# 3) Extra: resumo por par (ReferenceStudy × TestStudy) com estatísticas úteis
by_pair_summary <- pres_plot_df %>%
  group_by(ReferenceStudy, TestStudy) %>%
  summarise(
    n_points        = dplyr::n(),
    mean_Zsummary   = mean(Zsummary, na.rm = TRUE),
    sd_Zsummary     = sd(Zsummary, na.rm = TRUE),
    mean_Zdensity   = mean(Zdensity, na.rm = TRUE),
    mean_Zconnect   = mean(Zconnectivity, na.rm = TRUE),
    mean_medianRank = mean(medianRank, na.rm = TRUE),
    .groups = "drop"
  )

sheets_list$by_pair_summary <- by_pair_summary

# 4) Write to Excel (one workbook, multiple sheets)
xlsx_path <- file.path(output_dir, "Preservation_Plot_Data.xlsx")
write_xlsx(sheets_list, path = xlsx_path)

message("Excel saved at: ", xlsx_path)
