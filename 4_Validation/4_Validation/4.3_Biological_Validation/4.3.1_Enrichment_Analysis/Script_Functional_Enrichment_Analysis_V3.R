##### Functional Enrichment Analysis (GO) #####
# Purpose:
#   Perform GO enrichment (Biological Process, optionally MF/CC) for a target gene set
#   against study-specific backgrounds; reduce term redundancy with simplify(),
#   export per-study CSVs, and create a multi-study dot plot.
#
# Inputs:
#   - universe_genes: derived from each study's expression matrix rownames (ENSEMBL IDs)
#   - gene_sets: named list of target gene vectors (ENSEMBL IDs), e.g., "Positive"
#
# Outputs:
#   - One CSV per study with simplified GO results
#   - A combined results table (in-memory) filtered by FDR
#   - A dot plot summarizing top terms across studies (saved as JPEG)
#
# Notes:
#   - Assumes human data (OrgDb = org.Hs.eg.db) and keyType = "ENSEMBL"
#   - Uses BH correction (FDR) and Wang semantic similarity for simplify()
#   - Make sure your input gene IDs (target + backgrounds) are ENSEMBL and in the OrgDb

##### 0) Load required libraries #####
library(clusterProfiler)  # Gene set enrichment analysis
library(org.Hs.eg.db)     # Human genome annotation (OrgDb)
library(ggplot2)          # Plotting (dot plots)
library(enrichplot)       # Enrichment visualization helpers
library(dplyr)            # Data manipulation
library(stringr)          # String utilities (wrapping long labels)
library(dplyr)
library(writexl)

##### Helper: Enrichment + simplify wrapper #####
# Runs enrichGO with defensive checks and applies simplify() to reduce redundancy.
run_enrichGO_and_simplify <- function(
    gene_list,
    background,
    list_name,
    study_name,
    ontology         = "BP",
    pvalue_cutoff    = 0.05,
    min_gene_set_size= 10,
    max_gene_set_size= 500,
    simplify_cutoff  = 0.5,
    simplify_measure = "Wang"
) {
  # Sanity checks: keep only ENSEMBL IDs present in OrgDb universe
  gene_list   <- unique(gene_list)
  background  <- unique(background)
  
  # Run enrichGO; catch and return NULL on error (e.g., empty overlap)
  enrich_result <- tryCatch({
    enrichGO(
      gene          = gene_list,
      universe      = background,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENSEMBL",
      ont           = ontology,
      minGSSize     = min_gene_set_size,
      maxGSSize     = max_gene_set_size,
      pAdjustMethod = "BH",
      pvalueCutoff  = pvalue_cutoff
    )
  }, error = function(e) {
    message("⚠️ Enrichment failed for ", list_name, " | ", study_name, ": ", e$message)
    return(NULL)
  })
  
  # If no results (NULL or 0 rows), stop here
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    message("ℹ️ No GO terms returned for ", list_name, " | ", study_name, ".")
    return(NULL)
  }
  
  # Reduce redundancy with semantic similarity
  simplified_result <- tryCatch({
    simplify(
      enrich_result,
      cutoff     = simplify_cutoff,
      by         = "p.adjust",
      select_fun = min,
      measure    = simplify_measure
    )
  }, error = function(e) {
    message("⚠️ simplify() failed for ", list_name, " | ", study_name, ": ", e$message)
    # fall back to the original enrich_result if simplify fails
    enrich_result
  })
  
  # Convert to data.frame and annotate provenance
  result_df <- as.data.frame(simplified_result)
  if (nrow(result_df) == 0) {
    return(NULL)
  }
  result_df$List   <- list_name
  result_df$Study  <- study_name
  result_df$Ont    <- ontology
  
  return(result_df)
}

##### 1) Load datasets & define backgrounds #####
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC")

# Study expression matrices (rows = ENSEMBL IDs)
expDataAnjos      <- read.csv("Type_1_Studies/Anjos-Afonso_etal_2021/DATA/counts_Final.csv", row.names = 1)
expDataXie        <- read.csv("Type_1_Studies/Xie_etal_2019/DATA/counts_NonNormalized.csv", row.names = 1)
expDataFares      <- read.csv("Type_2_Studies/Fares_2017/DATA/counts_NonNormalized.csv", row.names = 1)
expDataSauvageau  <- read.csv("Type_2_Studies/Sauvageau_2014/DATA/counts_Final.csv", row.names = 1)
expDataXieCult    <- read.csv("Type_2_Studies/Xie_2019/DATA/counts_Final.csv", row.names = 1)

# Backgrounds per study (universe = expressed genes per dataset)
backgrounds <- list(
  Anjos_2020_T1   = rownames(expDataAnjos),
  Xie_2019_T1  = rownames(expDataXie),
  Fares_2016_T2   = rownames(expDataFares),
  Fares_2014_T2 = rownames(expDataSauvageau),
  Xie_2019_T2  = rownames(expDataXieCult)
)

# Target gene set (Positive modules, filtered by CountTotal ∈ {4,5})
positive_module_genes <- read.csv("Overlap_Uncult_Cult/Overlap/2_overlap_cult_kme_info_all.csv", row.names = 1)
positive_module_genes_filter <- positive_module_genes[positive_module_genes$CountTotal %in% c(4, 5), ]
positive_gene_list <- unique(positive_module_genes_filter$Ensembl)

# Define input gene sets to analyze (can add more lists later)
gene_sets <- list(
  Positive = positive_gene_list
)

##### 2) Run enrichment across studies #####
# Output directory for per-study CSVs and figures
output_dir <- "Overlap_Uncult_Cult/Enrichment_Analysis/Enrichment_GO/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- list()

for (list_name in names(gene_sets)) {
  gene_vector <- gene_sets[[list_name]]
  
  for (study_name in names(backgrounds)) {
    background_vector <- backgrounds[[study_name]]
    
    message("▶ Running enrichGO | List: ", list_name, " | Background: ", study_name)
    
    # Enrichment + redundancy reduction
    result_df <- run_enrichGO_and_simplify(
      gene_list          = gene_vector,
      background         = background_vector,
      list_name          = list_name,
      study_name         = study_name,
      ontology           = "BP",
      pvalue_cutoff      = 0.05,
      min_gene_set_size  = 10,
      max_gene_set_size  = 500,
      simplify_cutoff    = 0.5,
      simplify_measure   = "Wang"
    )
    
    # Save and collect if non-empty
    if (!is.null(result_df) && nrow(result_df) > 0) {
      out_csv <- file.path(output_dir, paste0(list_name, "_", study_name, ".csv"))
      write.csv(result_df, out_csv, row.names = FALSE)
      all_results[[paste0(list_name, "_", study_name)]] <- result_df
    }
  }
}

# Combine all studies into one table (if anything returned)
if (length(all_results) == 0) {
  stop("No enrichment results were returned across studies. Check ID types and filters.")
}
combined_results <- dplyr::bind_rows(all_results)

# Keep only significant terms (FDR < 0.05)
combined_results <- combined_results %>% filter(p.adjust < 0.05)

# (Optional) Frequency of terms across studies
freq_table <- combined_results %>% count(Description, sort = TRUE)

# Save the frequency
write.csv(freq_table, "Overlap_Uncult_Cult/Enrichment_Analysis/Enrichment_GO/Frequency_Table.csv", row.names = FALSE)

# Select top-N terms per (List × Study) by FDR
top_results <- combined_results %>%
  group_by(List, Study) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# Convert GeneRatio from "a/b" to numeric safely (a/b as decimal)
to_numeric_ratio <- function(x) {
  # robust parsing without eval(parse())
  parts <- strsplit(x, "/", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    num <- suppressWarnings(as.numeric(parts[1]))
    den <- suppressWarnings(as.numeric(parts[2]))
    if (!is.na(num) && !is.na(den) && den > 0) return(num / den)
  }
  return(NA_real_)
}
top_results$GeneRatio <- vapply(as.character(top_results$GeneRatio), to_numeric_ratio, numeric(1))

# Focus on the "Positive" list for plotting (tweak N as needed)
top_positive <- top_results %>%
  filter(List == "Positive") %>%
  arrange(p.adjust) %>%
  slice_head(n = 96)

# Wrap long GO labels for readability
top_positive$Description_wrapped <- stringr::str_wrap(top_positive$Description, width = 30)

# (Optional) enforce a specific study order in the plot
top_positive$Study <- factor(
  top_positive$Study,
  levels = c("Anjos_2020_T1", "Xie_2019_T1", "Fares_2016_T2", "Fares_2014_T2", "Xie_2019_T2")
)

# Dot plot: each point = GO term × study (size = GeneRatio, color = FDR)
plot_positive <- ggplot(top_positive,
                        aes(x = Study, y = reorder(Description_wrapped, p.adjust))) +
  geom_point(aes(size = GeneRatio, color = p.adjust), alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(2, 8), limits = c(0, NA)) +
  labs(
    title = "GO Biological Process — Enriched Terms (Positive Module Genes)",
    x = "Study",
    y = "GO Term",
    color = "Adjusted p-value",
    size  = "Gene Ratio"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title   = element_text(size = 12, face = "bold", color = "black", hjust = 0.5)
  )

ggsave(
  filename = file.path(output_dir, "GO_Positive.jpeg"),
  plot     = plot_positive,
  width    = 8, height = 8, dpi = 300
)

# Dot plot: each point = study x GO term   (size = GeneRatio, color = FDR)
top_positive$Study <- factor(
  top_positive$Study,
  levels = c("Anjos_2020_T1", "Xie_2019_T1", "Fares_2016_T2", "Fares_2014_T2", "Xie_2019_T2")
)

term_order <- top_positive %>%
  group_by(Description_wrapped) %>%
  summarise(min_p = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_p)

top_positive$Description_wrapped <- factor(
  top_positive$Description_wrapped,
  levels = term_order$Description_wrapped
)

plot_positive <- ggplot(
  top_positive,
  aes(x = Description_wrapped, y = Study)
) +
  geom_point(aes(size = GeneRatio, color = p.adjust), alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(2, 8), limits = c(0, NA)) +
  labs(
    title = "GO Biological Process — Enriched Terms (Positive Module Genes)",
    x = "GO Term",
    y = "Study",
    color = "Adjusted p-value",
    size  = "Gene Ratio"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(
      angle = 45, vjust = 1, hjust = 1, color = "black", size = 8 
    ),
    axis.text.y  = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title   = element_text(size = 12, face = "bold", color = "black", hjust = 0.5)
  )

ggsave(
  filename = file.path(output_dir, "GO_Positive_horizontal_terms_45deg.jpeg"),
  plot     = plot_positive,
  width    = 12, height = 6, dpi = 300
)

##### Filter and plot selected GO terms #####

# --- 1) Define the GO terms of interest ---
interest_terms_pretty <- c(
  "adaptive immune response",
  "angiogenesis",
  "antigen processing and presentation of exogenous peptide antigen",
  "positive regulation of protein phosphorylation",
  "extracellular matrix organization",
  "homeostasis of number of cells",
  "Positive regulation of protein phosphorylation",
  "response to growth factor"
)

# Create a lowercase version for case-insensitive matching
interest_terms_lc <- tolower(interest_terms_pretty)

# --- 2) Filter combined_results for the selected terms ---
df_interest <- combined_results %>%
  filter(List == "Positive") %>%                            # keep only Positive list
  mutate(desc_lc = tolower(Description)) %>%                # lowercase descriptions
  filter(desc_lc %in% interest_terms_lc)                    # match with target terms

# Report missing terms, if any
missing_terms <- setdiff(interest_terms_lc, unique(df_interest$desc_lc))
if (length(missing_terms) > 0) {
  message("⚠️ Terms not found after filtering: ",
          paste(interest_terms_pretty[match(missing_terms, interest_terms_lc)], collapse = "; "))
}

# Stop execution if no results remain
if (nrow(df_interest) == 0) {
  stop("No selected terms were found in combined_results. Check spelling/ontology/filtering.")
}

# --- 3) Convert GeneRatio from "a/b" to numeric ---
to_numeric_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    num <- suppressWarnings(as.numeric(parts[1]))
    den <- suppressWarnings(as.numeric(parts[2]))
    if (!is.na(num) && !is.na(den) && den > 0) return(num / den)
  }
  return(NA_real_)
}
df_interest$GeneRatio <- vapply(as.character(df_interest$GeneRatio), to_numeric_ratio, numeric(1))

#write.csv(df_interest, "Overlap_Uncult_Cult/Enrichment_Analysis/Enrichment_GO/df_4paper.csv", row.names = FALSE)
# --- 4) Assign pretty labels and fix order ---
# Map lowercase description back to the pretty version (respect original order)
map_pretty <- setNames(interest_terms_pretty, interest_terms_lc)  
df_interest$TermPretty <- map_pretty[df_interest$desc_lc]

# Fix the order of studies
df_interest$Study <- factor(
  df_interest$Study,
  levels = c("Anjos_2020_T1", "Xie_2019_T1", "Fares_2016_T2", "Fares_2014_T2", "Xie_2019_T2")
)

# Fix the order of terms on the X-axis (exactly as listed above)
df_interest$TermPretty <- factor(
  df_interest$TermPretty,
  levels = interest_terms_pretty
)

# Wrap long labels for better readability
df_interest$TermPretty_wrapped <- str_wrap(as.character(df_interest$TermPretty), width =20)

# --- 5) Plot: GO terms on X-axis, studies on Y-axis ---
plot_selected <- ggplot(
  df_interest,
  aes(x = TermPretty_wrapped, y = Study)
) +
  geom_point(aes(size = GeneRatio, color = p.adjust), alpha = 0.7) +
  scale_color_gradient(low = "darkblue",  high = "red") +
  scale_size_continuous(range = c(2, 8), limits = c(0, NA)) +
  labs(
    title = "Selected GO Biological Processes — Positive Module Genes",
    x = "GO Term",
    y = "Study",
    color = "FDR",
    size  = "Gene Ratio"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text        = element_text(family = "Arial"),
    axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title   = element_text(size = 12, face = "bold", color = "black", hjust = 0.5)
  )

# --- 6) Save plot ---
ggsave(
  filename = file.path(output_dir, "GO_Positive_selected_terms_45deg.svg"),
  plot     = plot_selected,
  width    = 14, height = 6, dpi = 300
)

#### Create Excel table with GO results + frequency across studies ####

# 2) Total number of distinct studies in the combined table
N_STUDIES <- combined_results %>% dplyr::pull(Study) %>% unique() %>% length()

# 3) Compute frequency per GO term across studies
#    (Key choice: use ID + Description to be robust to duplicated names)
freq_df <- combined_results %>%
  distinct(ID, Description, Study) %>%
  count(ID, Description, name = "FreqStudies") %>%
  mutate(Frequency = FreqStudies)  # now it's numeric, not "n/N"

# 4) Join frequency back to each row and reorder/select columns
final_go_table <- combined_results %>%
  left_join(freq_df, by = c("ID", "Description")) %>%
 dplyr::select(
    Study,
    ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
    geneID, Count, List, Ont,
    Frequency
  ) %>%
  arrange(desc(Frequency), p.adjust, Study)

# 5) (Optional) Also create a compact summary sheet: unique terms + frequency only
summary_by_term <- freq_df %>%
  dplyr::arrange(dplyr::desc(FreqStudies), Description)

# 6) Write Excel with two sheets
out_xlsx <- "Overlap_Uncult_Cult/Enrichment_Analysis/Enrichment_GO/GO_Combined_with_Frequency.xlsx"
writexl::write_xlsx(
  list(
    GO_results = final_go_table,
    Term_frequency = summary_by_term
  ),
  path = out_xlsx
)

message("✓ Excel saved: ", out_xlsx)

