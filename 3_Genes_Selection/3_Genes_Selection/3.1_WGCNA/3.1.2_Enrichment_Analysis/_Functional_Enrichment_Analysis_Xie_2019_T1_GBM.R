###############################################################################
# Script:        Functional_Enrichment_GO_BP_WGCNA_Modules.R
# Author:        Gustavo Bueno Moreira
# Last update:   2026-02-16
#
# Purpose:
#   - Perform GO-BP enrichment (clusterProfiler::enrichGO + simplify) for selected
#     WGCNA modules using an expressed-gene universe, export results (CSV/RDS),
#     and generate a multi-panel dotplot summary (top GO terms per module).
#
# Deliverables:
#   1) GO-BP enrichment tables (CSV) + enrichment objects (RDS) per module
#   2) Multi-panel dotplot (top terms per module) saved as a single image
#
# Inputs:
#   - DATA/counts_NonNormalized.csv (gene universe; first column must be Ensembl)
#   - RESULTS/WGCNA/3_kmeTable.csv  (kME table with Ensembl + moduleColor)
#
# Outputs:
#   - RESULTS/Enrichment_GO/GO_Analysis_<module>_GO_BP.csv
#   - RESULTS/Enrichment_GO/GO_Analysis_<module>_GO_BP.rds
#   - RESULTS/Enrichment_GO/painel_GO_BP_AllModules.jpeg
#
# Notes:
#   - This script assumes execution from the project root after setwd() below.
#   - Required input columns:
#       * universe_genes: "Ensembl"
#       * module_genes:   "Ensembl", "moduleColor"
#   - Modules without significant enrichment after simplify() are skipped.
###############################################################################

#### 0) Setup #################################################################

# 0.1) Setup: working directory
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC/Type_1_Studies/Xie_etal_2019/")

# 0.2) Options & constants
options(stringsAsFactors = FALSE)

# 0.3) Packages
library(clusterProfiler)  # GO enrichment analysis (enrichGO, simplify)
library(org.Hs.eg.db)     # human genome annotation for ENSEMBL IDs
library(ggplot2)          # plotting (ggplot-based dotplots)
library(enrichplot)       # enrichment visualization utilities (kept as in original)
library(dplyr)            # data manipulation (filter/arrange/mutate)
library(stringr)          # text wrapping for GO term descriptions
library(patchwork)        # combine plots into a single panel

# 0.4) Functions (helpers)
run_enrichment_analysis <- function(
    universe_genes,
    module_genes,
    selected_modules = NULL,
    ontology = "BP",
    p_value_cutoff = 0.05,
    min_gene_set_size = 10,
    max_gene_set_size = 500,
    output_prefix = "GO_Analysis",
    output_dir = NULL,
    simplify_cutoff = 0.5,
    simplify_measure = "Wang"
) {
  
  # Keep outputs organized and reproducible across runs
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  } else {
    output_dir <- getwd()
  }
  
  # Fail fast if expected columns are missing
  if (!"Ensembl" %in% colnames(universe_genes) |
      !"Ensembl" %in% colnames(module_genes) |
      !"moduleColor" %in% colnames(module_genes)) {
    stop("Missing 'Ensembl' or 'moduleColor' column in input.")
  }
  
  universe_gene_list <- unique(na.omit(universe_genes$Ensembl))
  
  if (is.null(selected_modules)) {
    module_gene_list <- module_genes$Ensembl
    selected_modules <- list("AllGenes")
  } else {
    selected_modules <- as.list(selected_modules)
  }
  
  for (mod in selected_modules) {
    
    if (mod != "AllGenes") {
      module_gene_list <- module_genes %>%
        filter(moduleColor == mod) %>%
        pull(Ensembl) %>%
        unique() %>%
        na.omit()
    }
    
    if (length(module_gene_list) < min_gene_set_size) {
      message(
        paste0(
          "[INFO] Skipping module '", mod, "' - Only ", length(module_gene_list),
          " genes (min required: ", min_gene_set_size, ")."
        )
      )
      next
    }
    
    message(
      paste0(
        "[INFO] Running enrichment for module '", mod, "' with ",
        length(module_gene_list), " genes..."
      )
    )
    
    enrich_go <- tryCatch({
      enrich_result <- enrichGO(
        gene          = module_gene_list,
        universe      = universe_gene_list,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENSEMBL",
        ont           = ontology,
        minGSSize     = min_gene_set_size,
        maxGSSize     = max_gene_set_size,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_value_cutoff
      )
      
      simplify(
        enrich_result,
        cutoff     = simplify_cutoff,
        by         = "p.adjust",
        select_fun = min,
        measure    = simplify_measure
      )
    }, error = function(e) {
      message(paste("[ERROR] Module:", mod, "-", e$message))
      return(NULL)
    })
    
    if (is.null(enrich_go) || length(enrich_go) == 0 ||
        is.null(enrich_go@result) || nrow(enrich_go@result) == 0) {
      message(
        paste(
          "[INFO] No significant enrichment found for module '",
          mod,
          "' – skipping export.",
          sep = ""
        )
      )
      next
    }
    
    message(
      paste0(
        "[INFO] Enrichment completed for module '", mod, "' – ",
        nrow(enrich_go@result), " GO terms found."
      )
    )
    
    csv_path <- file.path(output_dir, paste0(output_prefix, "_", mod, "_GO_", ontology, ".csv"))
    rds_path <- file.path(output_dir, paste0(output_prefix, "_", mod, "_GO_", ontology, ".rds"))
    
    write.csv(enrich_go@result, csv_path, row.names = FALSE)
    saveRDS(enrich_go, rds_path)
    
    message(paste0("[✔️] Results saved: ", basename(csv_path)))
    message(paste0("[✔️] RDS saved: ", basename(rds_path)))
  }
  
  message("Enrichment analysis finished for all modules.")
}

# 0.5) Load data
message("[INFO] Loading input datasets...")

universe_genes <- read.csv("DATA/counts_NonNormalized.csv")
colnames(universe_genes)[1] <- "Ensembl"

module_genes <- read.csv("RESULTS/WGCNA/3_kmeTable.csv")

#### 1) Deliverable 1: GO-BP enrichment per WGCNA module ######################

# 1.1) Run enrichment and export CSV/RDS per module
message("[INFO] Starting GO enrichment (BP) for selected modules...")

# Note: If you see the message "In rep(yes, length.out = len) : 'x' is NULL so the result will be NULL",
# it means that the corresponding module does not have significant enrichment. In this case,
# you should remove that module from the analysis to avoid errors.
run_enrichment_analysis(
  universe_genes     = universe_genes,
  module_genes       = module_genes,
  selected_modules   = c("green", "blue", "brown", "darkturquoise", "turquoise", "tan", "yellow"),
  ontology           = "BP",
  p_value_cutoff     = 0.05,
  min_gene_set_size  = 10,
  max_gene_set_size  = 500,
  output_prefix      = "GO_Analysis",
  output_dir         = "RESULTS/Enrichment_GO/",
  simplify_cutoff    = 0.7,
  simplify_measure   = "Wang"
)

# 1.2) Load exported results and build per-module dotplots
message("[INFO] Building dotplot panel from exported enrichment tables...")

results_files <- list.files(
  "RESULTS/Enrichment_GO/",
  pattern = "_GO_BP.csv$",
  full.names = TRUE
)

plot_list <- list()

for (file in results_files) {
  
  result_df <- read.csv(file)
  
  module_name <- gsub("GO_Analysis_|_GO_BP.csv", "", basename(file))
  result_df$module <- module_name
  
  # Convert GeneRatio if it's stored as a string (e.g., "5/150")
  if (is.character(result_df$GeneRatio)) {
    result_df$GeneRatio <- sapply(result_df$GeneRatio, function(x) eval(parse(text = x)))
  }
  
  top_terms <- result_df %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    arrange(GeneRatio) %>%
    mutate(Description_wrapped = str_wrap(Description, width = 30))
  
  p <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description_wrapped, GeneRatio))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(name = "Adjusted p-value", low = "red", high = "blue") +
    scale_size(name = "Gene Count") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    labs(
      title = module_name,
      x = "",
      y = ""
    ) +
    guides(
      size  = guide_legend(order = 1),
      color = guide_colorbar(order = 2)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y  = element_text(size = 11, color = "black"),
      axis.text.x  = element_text(size = 11, color = "black"),
      legend.title = element_text(size = 10, color = "black"),
      legend.text  = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 7, color = "black"),
      axis.title.y = element_text(size = 7, color = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
  
  plot_list[[module_name]] <- p
}

# 1.3) Combine plots and save output
final_panel <- wrap_plots(plot_list, ncol = 2)

ggsave(
  filename = "RESULTS/Enrichment_GO/painel_GO_BP_AllModules.jpeg",
  plot     = final_panel,
  width    = 15,
  height   = 13,
  dpi      = 300
)

message("[✔️] Dotplot panel saved: painel_GO_BP_AllModules.jpeg")

#### Software environment #####################################################
# R version and package versions
sessionInfo()
