###############################################################################
# Script:        Studies Selection
# Author:        Gustavo Bueno Moreira
# Last update:   2026-02-16
#
# Purpose:
#   - Perform an exploratory analysis of the selected studies, summarizing
#     key experimental characteristics and data availability through
#     descriptive visualizations.
#
# Deliverables:
#   1) Heatmap for all included studies (Study × Cell source; sample counts)
#   2) Pie chart for Type 3 studies (Fresh vs Expanded vs Mixed)
#   3) Stacked bar chart (Study type × Sequencing type)
#
# Inputs:
#   - 1_Data/Tabela_Estudos_Final_V7.xlsx
#     Required columns:
#       - Author, Year, GEO ID, Study type, Cell source, Samples per subpopulation
#       - Inclusion status: Paper, Reason for exclusion: Paper
#       - Cell culture status, Sequencing type
#
# Outputs:
#   - Heatmap_All_Studies.png
#   - Pie.png
#   - barplot.png
###############################################################################

#### 0) Setup #################################################################

# 0.1) Setup: working directory ------------------------------------------------
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022")

# 0.2) Options & constants -----------------------------------------------------
options(stringsAsFactors = FALSE)

# 0.3) Packages ----------------------------------------------------------------
library(readxl)   # read Excel inputs
library(dplyr)    # data manipulation (filter, mutate, summarize)
library(tidyr)    # complete missing combinations for heatmap grid
library(ggplot2)  # plotting and export (ggsave)
library(viridis)  # color utilities (kept to preserve original environment)
library(openxlsx) # Excel I/O utilities (kept to preserve original environment)

# 0.4) Functions (helpers) -----------------------------------------------------
create_study_id <- function(df) {
  # Stable visible identifier used across deliverables: Author_Year_StudyType
  df$Study_ID <- paste(df$Author, df$Year, df$`Study type`, sep = "_")
  df <- df |> dplyr::relocate(Study_ID, .after = Year)
  return(df)
}

recode_cell_source <- function(x) {
  # Visualization recode: PL/EB collapsed into "Others"
  dplyr::case_when(
    x %in% c("PL", "EB") ~ "Others",
    TRUE ~ x
  )
}

standardize_culture_status <- function(x) {
  # Standardize labels used in pie chart (preserving original mapping)
  dplyr::case_when(
    x %in% c("Fresh", "fresh") ~ "Fresh",
    x %in% c("Expanded", "Cultured", "expanded", "cultured") ~ "Expanded",
    x %in% c("Mixed", "mixed") ~ "Mixed",
    TRUE ~ "Mixed"
  )
}

# 0.5) Load data ---------------------------------------------------------------
message("0.5) Loading study table from Excel...")

input_file <- "1_Data/Tabela_Estudos_Final_V7.xlsx"
studies <- read_excel(input_file)

# Create Study_ID (needed both before and after filtering; keep as in original script)
studies <- create_study_id(studies)

# Total number of studies (based on Study_ID) - kept as in original script
total <- length(unique(studies$Study_ID))

# Keep only one row per Study_ID (for exclusion reason summary)
studies_unique <- studies |> distinct(Study_ID, .keep_all = TRUE)

# Count exclusion reasons (including empty cells) and print
message("0.5) Summarizing exclusion reasons (one row per Study_ID)...")
reason_counts <- studies_unique |>
  dplyr::count(`Reason for exclusion: Paper`, name = "Count") |>
  dplyr::arrange(desc(Count))
print(reason_counts)

# Keep only Included (as in original script)
message("0.5) Filtering to included studies...")
studies <- studies |>
  dplyr::filter(`Inclusion status: Paper` == "Included")

# Recreate Study_ID after filtering (preserve original flow)
studies <- create_study_id(studies)

#### 1) Deliverable 1: Heatmap (all studies) ##################################

message("1) Building heatmap (all included studies)...")

# 1.1) Minimal base: keys, recode PL/EB, ensure numeric ------------------------
df0 <- studies |>
  dplyr::select(
    Author, Year, `GEO ID`, `Study type`,
    `Cell source`, `Samples per subpopulation`
  ) |>
  dplyr::mutate(
    StudyUID = paste(Author, Year, `GEO ID`, sep = "_"), # internal key
    X_label  = paste(Author, Year, sep = "_"),           # visible label
    `Cell source` = recode_cell_source(`Cell source`),
    `Samples per subpopulation` = as.numeric(`Samples per subpopulation`)
  )

# 1.2) Aggregate: sum per Study × Cell source ---------------------------------
agg <- df0 |>
  dplyr::group_by(StudyUID, X_label, `Study type`, `Cell source`) |>
  dplyr::summarise(
    total_samples = sum(`Samples per subpopulation`, na.rm = TRUE),
    .groups = "drop"
  )

# 1.3) Column order: Study type → Author → Year -------------------------------
order_tbl <- df0 |>
  dplyr::distinct(StudyUID, X_label, `Study type`, Author, Year) |>
  dplyr::arrange(`Study type`, Author, Year, StudyUID)

all_uids <- order_tbl$StudyUID

all_sources <- sort(unique(df0$`Cell source`))
if ("Others" %in% all_sources) {
  all_sources <- c(setdiff(all_sources, "Others"), "Others") # keep "Others" last
}

# 1.4) Complete grid (show missing combos as NA) -------------------------------
plot_df <- agg |>
  tidyr::complete(
    StudyUID = all_uids,
    `Cell source` = all_sources,
    fill = list(total_samples = NA_real_)
  ) |>
  dplyr::left_join(order_tbl, by = "StudyUID") |>
  dplyr::mutate(
    StudyUID      = factor(StudyUID, levels = order_tbl$StudyUID),
    `Cell source` = factor(`Cell source`, levels = rev(all_sources))
  )

cols <- colorRampPalette(c("blue", "red"))(100)

# 1.5) Plot heatmap ------------------------------------------------------------
p_all <- ggplot(plot_df, aes(x = StudyUID, y = `Cell source`, fill = total_samples)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(is.na(total_samples), "", total_samples)),
    size = 7, color = "white", na.rm = TRUE
  ) +
  scale_fill_gradientn(
    colours = cols,
    na.value = "#F2F2F2", name = "Samples"
  ) +
  scale_x_discrete(labels = setNames(order_tbl$X_label, order_tbl$StudyUID)) +
  labs(x = NULL, y = NULL, title = "") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15),
    panel.grid  = element_blank(),
    plot.title  = element_text(face = "bold")
  )

p_all

out_heatmap_all <- "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Studies_Selection/Paper/Heatmap_All_Studies_dissertação.png"
ggsave(
  filename = out_heatmap_all,
  plot = p_all,
  width = 15, height = 5, dpi = 300
)

#### 2) Deliverable 2: Pie chart (Fresh vs Expanded — T3 only) #################

message("2) Building pie chart (Fresh vs Expanded vs Mixed) for Type 3 only...")

# 2.1) Keep only Type 3 rows ---------------------------------------------------
t3_only <- studies |>
  dplyr::filter(`Study type` == "T3")

# 2.2) Collapse to one row per study (Mixed if >1 status within study) ---------
per_study_t3 <- t3_only |>
  dplyr::group_by(Study_ID) |>
  dplyr::summarise(
    `Cell culture status` = ifelse(
      dplyr::n_distinct(`Cell culture status`) > 1,
      "Mixed",
      dplyr::first(`Cell culture status`)
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    `Cell culture status` = standardize_culture_status(`Cell culture status`),
    `Cell culture status` = factor(`Cell culture status`, levels = c("Fresh", "Expanded", "Mixed"))
  )

# 2.3) Count + percentages -----------------------------------------------------
pie_df_t3 <- per_study_t3 |>
  dplyr::count(`Cell culture status`, name = "n") |>
  dplyr::mutate(pct = n / sum(n))

# 2.4) Plot pie chart ----------------------------------------------------------
p_pie_t3 <- ggplot(pie_df_t3, aes(x = "", y = n, fill = `Cell culture status`)) +
  geom_col(width = 1, color = "white", linewidth = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "Fresh" = "#93006B",
    "Expanded" = "#D73027",
    "Mixed" = "#D73027"
  )) +
  geom_text(
    aes(label = paste0(round(pct * 100), "%")),
    position = position_stack(vjust = 0.5),
    color = "white", fontface = "bold", size = 5
  ) +
  theme_void() +
  labs(title = "", fill = "Cell culture status") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

p_pie_t3

out_pie_t3_only <- "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Studies_Selection/Paper/Pie_CultureCondition_T3_only.png"
ggsave(
  filename = out_pie_t3_only,
  plot = p_pie_t3,
  width = 15, height = 15, units = "cm", dpi = 600
)

#### 3) Deliverable 3: Stacked bar (Study type × Sequencing type) ##############

message("3) Building stacked bar chart (Study type × Sequencing type)...")

per_study <- studies |>
  dplyr::distinct(Study_ID, .keep_all = TRUE) |>
  dplyr::select(Study_ID, `Study type`, `Sequencing type`)

bar_df <- per_study |>
  dplyr::count(`Study type`, `Sequencing type`) |>
  dplyr::mutate(`Study type` = factor(`Study type`, levels = c("T1", "T2", "T3")))

# Keep the original standalone call (even though it is not assigned) -----------
scale_x_discrete(labels = c(
  "T1" = "T1",
  "T2" = "T2",
  "T3" = "Validation datasets"
))

p_bar <- ggplot(bar_df, aes(x = `Study type`, y = n, fill = `Sequencing type`)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(
    "RNA-Seq" = "#2B00D3",
    "Microarray" = "#670097"
  )) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  scale_x_discrete(labels = c(
    "T1" = "T1",
    "T2" = "T2",
    "T3" = "Validation datasets"
  )) +
  labs(
    x = "Study type",
    y = "Number of studies",
    title = "",
    fill = "Sequencing type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9)
  )

out_barplot <- "2a_b_WGCNA_public_data_fresh_cultivated_HSC/Studies_Selection/Paper/barplot.png"
ggsave(
  filename = out_barplot,
  plot = p_bar,
  width = 15, height = 10, units = "cm", dpi = 600
)

#### Software environment #####################################################
sessionInfo()


