# ======================================================================
# Title: TF_Motif_Consistency_All_Studies.R
# Author: Gustavo Bueno Moreira
#
# Purpose:
#   Identify transcription factor (TF) motifs consistently enriched across
#   multiple studies. The pipeline loads per-study JASPAR enrichment results,
#   filters significant hits (OR > 1 & p < 0.05), visualizes overlap with a
#   5-set Venn diagram (by JASPAR ID), and exports lists of motifs/IDs that
#   are present in all studies.
#
# Notes:
#   - JASPAR source standardized to JASPAR2022 (human CORE).
#   - Family annotations and dendrogram both rely on JASPAR2022.
# ======================================================================

#### 0) Libraries, options & working directory ################################

# 0.1) Set working directory
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC")

# 0.2) Load libraries (no duplicates; JASPAR2022 only)
suppressPackageStartupMessages({
  library(dplyr)          # data manipulation
  library(purrr)          # functional programming helpers
  library(tidyr)          # data reshaping
  library(tibble)         # tidy data frames
  library(stringr)        # string operations
  library(ggplot2)        # plotting
  library(ggraph)         # dendrogram/graph viz
  library(VennDiagram)    # Venn diagrams
  library(TFBSTools)      # TF motif handling
  library(JASPAR2022)     # JASPAR CORE (human) - 2022 release
  library(universalmotif) # motif comparison & conversion
})

#### 1) Load JASPAR enrichment results (CSV) ##################################

# 1.1) Type 1 studies
Anjos  <- read.csv("Type_1_Studies/Anjos-Afonso_etal_2021/RESULTS/Transcription _Factor/jaspar2022_human_fg_bg_fisher_minScore80_BG_GCmatched_5to1.csv")
Xie_T1 <- read.csv("Type_1_Studies/Xie_etal_2019/RESULTS/Transcription _Factor/jaspar2022_human_fg_bg_fisher_minScore80_BG_GCmatched_5to1.csv")

# 1.2) Type 2 studies
Fares     <- read.csv("Type_2_Studies/Fares_2017/RESULTS/Transcription _Factor/jaspar2022_human_fg_bg_fisher_minScore80_BG_GCmatched_5to1.csv")
Sauvageau <- read.csv("Type_2_Studies/Sauvageau_2014/RESULTS/Transcription _Factor/jaspar2022_human_fg_bg_fisher_minScore80_BG_GCmatched_5to1.csv")
Xie_T2    <- read.csv("Type_2_Studies/Xie_2019/RESULTS/Transcription _Factor/jaspar2022_human_fg_bg_fisher_minScore80_BG_GCmatched_5to1.csv")


#### 2) Helper functions ######################################################

# 2.1) Pick (jaspar_id, motif_id). If motif_id is missing, fall back to tf_name.
pick_pair <- function(df) {
  if ("motif_id" %in% names(df)) {
    dplyr::transmute(df,
                     jaspar_id = as.character(.data$jaspar_id),
                     motif_id  = as.character(.data$motif_id))
  } else if ("tf_name" %in% names(df)) {
    dplyr::transmute(df,
                     jaspar_id = as.character(.data$jaspar_id),
                     motif_id  = as.character(.data$tf_name))
  } else {
    stop("Neither 'motif_id' nor 'tf_name' were found in this dataframe.")
  }
}

# 2.2) Filter significant hits: OR > 1 and p < 0.05, then keep unique pairs
filter_hits <- function(df) {
  df %>%
    dplyr::mutate(
      odds_ratio = suppressWarnings(as.numeric(.data$odds_ratio)),
      p_value    = suppressWarnings(as.numeric(.data$p_value))
    ) %>%
    dplyr::filter(!is.na(.data$odds_ratio) & !is.na(.data$p_value) &
                    .data$odds_ratio > 1 & .data$p_value < 0.05) %>%
    pick_pair() %>%
    dplyr::distinct()
}


#### 3) Bundle studies ########################################################

# 3.1) Create a named list with all study data.frames
studies <- list(
  Anjos    = Anjos,
  Xie_T1   = Xie_T1,
  Fares    = Fares,
  Sauvageau= Sauvageau,
  Xie_T2   = Xie_T2
)


#### 4) Apply filtering per study ############################################

# 4.1) Apply the significance filter to each study
filtered <- purrr::map(studies, filter_hits)


#### 5) Combine filtered data #################################################

# 5.1) Stack all filtered rows and track study of origin
combined_df <- dplyr::bind_rows(filtered, .id = "Study")

# 5.2) Keep only valid rows (defensive)
combined_df <- combined_df %>%
  dplyr::filter(!is.na(jaspar_id), !is.na(motif_id))

# 5.3) Count how many studies are present
n_studies <- length(filtered)

# 5.4) Bar plot: total number of motifs per study (thin borders, x labels at 45°)
motifs_per_study <- combined_df %>%
  dplyr::distinct(Study, jaspar_id, motif_id) %>%
  dplyr::count(Study, name = "n_motifs")

p_motifs <- ggplot(motifs_per_study, aes(x = Study, y = n_motifs)) +
  geom_col(fill = "#08306b", color = "black", linewidth = 0.25) + # dark blue bars, thin plot border only
  labs(x = "Study", y = "Total motifs (filtered)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25), # thin border around plot
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


# Save plot (keeps your folder structure)
ggsave("Transcription _Factor/Consensus_Motif/Bar_Total_Motifs_Per_Study.png", p_motifs, width = 7, height = 4.5, dpi = 300)

#### 6) Venn diagram (JASPAR IDs only) #######################################

# 6.1) Build per‑study JASPAR‑ID vectors (unique, drop NAs just in case)
venn_list <- list(
  Anjos      = unique(na.omit(filtered$Anjos$jaspar_id)),
  Xie_Type1  = unique(na.omit(filtered$Xie_T1$jaspar_id)),
  Fares      = unique(na.omit(filtered$Fares$jaspar_id)),
  Sauvageau  = unique(na.omit(filtered$Sauvageau$jaspar_id)),
  Xie_Type2  = unique(na.omit(filtered$Xie_T2$jaspar_id))
)

# 6.2) Dynamic category labels with counts
category_labels <- c(
  paste0("Anjos\n(",      length(venn_list$Anjos),      ")"),
  paste0("Xie_Type1\n(",  length(venn_list$Xie_Type1),  ")"),
  paste0("Fares\n(",      length(venn_list$Fares),      ")"),
  paste0("Sauvageau\n(",  length(venn_list$Sauvageau),  ")"),
  paste0("Xie_Type2\n(",  length(venn_list$Xie_Type2),  ")")
)

# 6.3) Output path
out_file <- "Transcription _Factor/Consensus_Motif/Diagrama_Venn_TF_JASPAR.png"

# 6.4) Draw & save
venn.plot <- VennDiagram::venn.diagram(
  x               = venn_list,
  filename        = out_file,
  output          = TRUE,
  
  # Styling (kept as provided)
  category.names  = category_labels,
  col             = NA,
  fill            = c("#FF0000", "#FFA500", "#FFD700", "#FF6347", "#FFF749"),
  alpha           = 0.2,
  lwd             = 1,
  cex             = 0.9,
  fontface        = "plain",
  fontfamily      = "sans",
  cat.cex         = 1,
  cat.col         = c("#B22222", "#FF8C00", "#DAA520", "#CD5C5C", "#BDB76B"),
  cat.fontface    = "plain",
  cat.fontfamily  = "sans",
  cat.pos         = c(-10, -90, 180, 0, 90),
  cat.dist        = c(0.25, 0.3, 0.2, -0.22, 0.27),
  margin          = 0.2,
  sep.dist        = 0.01,
  scaled          = FALSE
)


#### 7) Annotate motifs with TF families (JASPAR 2022) ########################

# 7.1) Fetch human CORE PFMs from JASPAR 2022 (no SQL connection)
opts <- list(species = 9606, collection = "CORE")
pfm_set <- TFBSTools::getMatrixSet(JASPAR2022, opts)

# 7.2) Build metadata table (motif-level)
get_tag <- function(x, keys) {
  for (k in keys) {
    v <- tryCatch(x@tags[[k]], error = function(e) NULL)
    if (!is.null(v)) return(paste(as.character(v), collapse = ";"))
  }
  NA_character_
}

meta_tbl <- purrr::map_df(pfm_set, function(pfm) {
  tibble::tibble(
    matrix_id  = TFBSTools::ID(pfm),
    name       = TFBSTools::name(pfm),
    geneSymbol = get_tag(pfm, c("geneSymbol","gene_name","symbol","TF")),
    family     = get_tag(pfm, c("family","Family")),
    species    = get_tag(pfm, c("species","tax_id","tax_group"))
  )
})

# 7.3) Version-insensitive join (ignore .version in matrix_id)
strip_version <- function(x) sub("\\.\\d+$", "", as.character(x))

meta_tbl_base <- meta_tbl %>%
  dplyr::mutate(matrix_base = strip_version(matrix_id))

combined_base <- combined_df %>%
  dplyr::mutate(jaspar_base = strip_version(jaspar_id))

# Join by base IDs; mantém repetições de combined_df
combined_with_family <- combined_base %>%
  dplyr::inner_join(meta_tbl_base %>% dplyr::select(matrix_base, family, geneSymbol),
                    by = c("jaspar_base" = "matrix_base"))

# 7.4) Recorrência por jaspar_base (quantos estudos contêm)
recurrence_tbl <- combined_base %>%
  dplyr::distinct(Study, jaspar_base) %>%
  dplyr::count(jaspar_base, name = "n_studies_present") %>%
  dplyr::mutate(recurrence = paste0(n_studies_present, "/", n_studies))

combined_with_family <- combined_with_family %>%
  dplyr::left_join(recurrence_tbl, by = "jaspar_base") %>%
  dplyr::select(Study, jaspar_id, jaspar_base, family, n_studies_present, recurrence)

# 7.5) Save & console summary
write.csv(combined_with_family,
          "Transcription _Factor/Consensus_Motif/TF_motifs_with_family.csv",
          row.names = FALSE)

message("Annotated motifs with family info (JASPAR 2022): ", nrow(combined_with_family),
        "; distinct JASPAR base IDs: ",
        combined_with_family %>% dplyr::distinct(jaspar_base) %>% nrow())

#### 8) Circular dendrogram of motifs (version-insensitive; study-level recurrence) ########################

##### 8.0) Output directory ####################################################
# 8.0.1) Define output dir and ensure it exists
out_dir <- file.path(getwd(), "Transcription _Factor/Consensus_Motif/")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

##### 8.1) Minimal master table and presence matrix ###########################
# 8.1.1) Keep one row per (Study × jaspar_base) and preserve family
motifs_clean <- combined_with_family %>%
  dplyr::select(Study, jaspar_base, family) %>%
  dplyr::filter(!is.na(Study), !is.na(jaspar_base)) %>%
  dplyr::distinct(Study, jaspar_base, .keep_all = TRUE)

# 8.1.2) Build presence matrix (motif × study, 0/1)
presence_wide <- motifs_clean %>%
  dplyr::select(Study, jaspar_base) %>%
  dplyr::mutate(present = 1L) %>%
  tidyr::pivot_wider(
    names_from  = Study,
    values_from = present,
    values_fill = 0L
  )

# 8.1.3) Collapse to a single family per jaspar_base (first non-NA)
family_tbl <- motifs_clean %>%
  dplyr::group_by(jaspar_base) %>%
  dplyr::summarise(
    family = { x <- family[!is.na(family)]; if (length(x) == 0) NA_character_ else x[1] },
    .groups = "drop"
  )

# 8.1.4) Join family and compute recurrence (n_studies_present and factor 1/5..5/5)
motif_study_mat <- presence_wide %>%
  dplyr::left_join(family_tbl, by = "jaspar_base") %>%
  dplyr::relocate(family, .after = jaspar_base)

study_cols <- setdiff(colnames(motif_study_mat), c("jaspar_base", "family"))

motif_study_mat <- motif_study_mat %>%
  dplyr::mutate(
    n_studies_present = rowSums(as.matrix(dplyr::select(., dplyr::all_of(study_cols)))),
    recurrence = factor(
      paste0(n_studies_present, "/5"),
      levels = c("1/5","2/5","3/5","4/5","5/5"),
      ordered = TRUE
    )
  ) %>%
  dplyr::arrange(dplyr::desc(n_studies_present), jaspar_base)

# 8.1.5) Master table (for downstream plotting)
motif_master <- motif_study_mat %>%
  dplyr::select(jaspar_base, family, n_studies_present, recurrence)

##### 8.2) Retrieve latest PFMs per jaspar_base (JASPAR CORE) #################
# 8.2.1) Bases of interest
bases <- unique(motif_master$jaspar_base)

# 8.2.2) Query JASPAR (human, all_versions) and convert to universalmotif
opts <- list(species = 9606, collection = "CORE", all_versions = TRUE)
pfm_list <- getMatrixSet(JASPAR2022, opts)

id_full_vec <- vapply(pfm_list, ID, character(1))               # e.g., "MA1976.1"
um_list     <- universalmotif::convert_motifs(pfm_list)          # keep order

# 8.2.3) Table with ids and versions; filter to our bases
motif_tbl <- tibble::tibble(
  motif_obj = um_list,
  id_full   = id_full_vec
) %>%
  dplyr::mutate(
    jaspar_base = stringr::str_remove(id_full, "\\.\\d+$"),
    version_num = as.integer(stringr::str_extract(id_full, "(?<=\\.)\\d+$"))
  ) %>%
  dplyr::filter(jaspar_base %in% bases)

# 8.2.4) Keep the latest version per jaspar_base
motif_latest <- motif_tbl %>%
  dplyr::group_by(jaspar_base) %>%
  dplyr::slice_max(order_by = version_num, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# 8.2.5) Warn on missing bases
missing_bases <- setdiff(bases, motif_latest$jaspar_base)
if (length(missing_bases) > 0) {
  warning("PWMs not found for ", length(missing_bases), " bases: ",
          paste(missing_bases, collapse = ", "))
}

##### 8.3) Build similarity, cluster, and dendrogram layout ###################
# 8.3.1) Align motif list to master table order and name it by jaspar_base
motifs_um <- motif_master %>%
  dplyr::filter(jaspar_base %in% motif_latest$jaspar_base) %>%
  dplyr::arrange(jaspar_base) %>%
  dplyr::left_join(dplyr::select(motif_latest, jaspar_base, motif_obj), by = "jaspar_base")

motif_list <- motifs_um$motif_obj
names(motif_list) <- motifs_um$jaspar_base

# 8.3.2) Pairwise similarity (PCC; allow reverse-complement)
sim_matrix <- universalmotif::compare_motifs(
  motif_list,
  method = "PCC",
  tryRC = TRUE,
  min.overlap = 6,
  score.strat = "a.mean",
  use.type = "PPM"
)

# 8.3.3) Normalize to [0,1], derive distance, cluster (average linkage)
sim_norm <- (sim_matrix + 1) / 2
diag(sim_norm) <- 1
dist_matrix <- 1 - sim_norm
hc <- hclust(as.dist(dist_matrix), method = "average")
leaf_order <- hc$labels[hc$order]

# 8.3.4) Map dendrogram labels → jaspar_base (prefer motif @name, fallback @altname)
get_label <- function(m) {
  lbl <- m@name
  if (is.null(lbl) || is.na(lbl) || lbl == "") {
    if (!is.null(m@altname) && !is.na(m@altname) && m@altname != "") lbl <- m@altname
  }
  lbl
}

mapping_tbl <- motif_latest %>%
  dplyr::mutate(label = vapply(motif_obj, get_label, character(1))) %>%
  dplyr::select(label, jaspar_base)

# 8.3.5) Leaf order to jaspar_base order; warn if any label is unmapped
leaf_tbl <- tibble::tibble(label = leaf_order) %>%
  dplyr::mutate(leaf_index = dplyr::row_number()) %>%
  dplyr::left_join(mapping_tbl, by = "label")

n_na <- sum(is.na(leaf_tbl$jaspar_base))
if (n_na > 0) {
  warning("Unmapped dendrogram labels: ", n_na,
          ". Inspect `leaf_tbl %>% dplyr::filter(is.na(jaspar_base))`.")
}

leaf_order_jaspar <- leaf_tbl$jaspar_base %>%
  as.character() %>%
  { .[!is.na(.)] } %>%
  unique()

# 8.3.6) Reorder master table by dendrogram order
motif_master_ordered <- motif_master %>%
  dplyr::filter(jaspar_base %in% leaf_order_jaspar) %>%
  dplyr::mutate(jaspar_base = factor(jaspar_base, levels = leaf_order_jaspar)) %>%
  dplyr::arrange(jaspar_base)

##### 8.4) Compute circular layout and label geometry #########################
# 8.4.1) Layout
lay <- ggraph::create_layout(hc, layout = "dendrogram", circular = TRUE)

# 8.4.2) Join leaf positions with recurrence and prepare label/tick geometry
leaves <- lay %>%
  dplyr::filter(leaf) %>%
  dplyr::select(x, y, label) %>%
  dplyr::left_join(mapping_tbl, by = "label") %>%
  dplyr::left_join(motif_master_ordered, by = "jaspar_base") %>%
  dplyr::filter(!is.na(jaspar_base)) %>%
  dplyr::mutate(
    label_plot = jaspar_base,
    theta      = atan2(y, x),
    r_leaf     = sqrt(x^2 + y^2)
  )

# 8.4.3) Fixed ring radius for labels + tick length/padding (consistent with your style)
tick_len <- 0.035 * max(leaves$r_leaf, na.rm = TRUE)
pad      <- 0.04  * max(leaves$r_leaf, na.rm = TRUE)
R_label  <- max(leaves$r_leaf + tick_len, na.rm = TRUE) + pad

# 8.4.4) Tick coordinates and label angles (left side flipped 180°, right side normal)
leaves <- leaves %>%
  dplyr::mutate(
    x_tick0 = r_leaf              * cos(theta),
    y_tick0 = r_leaf              * sin(theta),
    x_tick1 = (r_leaf + tick_len) * cos(theta),
    y_tick1 = (r_leaf + tick_len) * sin(theta),
    x_lab   = R_label * cos(theta),
    y_lab   = R_label * sin(theta),
    angle_deg_raw = theta * 180 / pi,
    angle_lab     = ifelse(cos(theta) < 0, angle_deg_raw + 180, angle_deg_raw),
    hjust_lab     = ifelse(cos(theta) < 0, 1, 0)
  )

##### 8.5) Plot(s): basic and with recurrence symbols #########################
# 8.5.1) Common sizing (kept from your code)
label_size_pt <- 18
line_width_mm <- 0.6
shape_size_mm <- 8

# 8.5.2) Ensure ordered factor for clean legend
leaves$recurrence <- factor(leaves$recurrence, levels = c("1/5","2/5","3/5","4/5","5/5"))

# 8.5.3) BASIC circular dendrogram (edges + ticks + outer labels)
p_circ_basic <- ggraph(lay) +
  geom_edge_elbow(linewidth = line_width_mm, colour = "grey35") +
  geom_segment(
    data = leaves,
    aes(x = x_tick0, y = y_tick0, xend = x_tick1, yend = y_tick1),
    linewidth = line_width_mm, colour = "grey35"
  ) +
  geom_text(
    data = leaves,
    aes(x = x_lab, y = y_lab, label = label_plot, angle = angle_lab, hjust = hjust_lab),
    size = label_size_pt / ggplot2::.pt,
    colour = "black"
  ) +
  theme_void() +
  theme(
    plot.margin   = margin(20, 28, 20, 28, unit = "mm"),
    panel.spacing = unit(0, "mm")
  ) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_equal(clip = "off")

# 8.5.4) WITH RECURRENCE SYMBOLS (only 5/5 filled red; others white)
p_circ <- ggraph(lay) +
  geom_edge_elbow(linewidth = line_width_mm, colour = "grey35") +
  geom_segment(
    data = leaves,
    aes(x = x_tick0, y = y_tick0, xend = x_tick1, yend = y_tick1),
    linewidth = line_width_mm, colour = "grey35"
  ) +
  geom_point(
    data = leaves,
    aes(x = x_tick1, y = y_tick1, shape = recurrence, fill = recurrence),
    size = shape_size_mm, stroke = 0.5, colour = "grey20"
  ) +
  geom_text(
    data = leaves,
    aes(x = x_lab, y = y_lab, label = label_plot, angle = angle_lab, hjust = hjust_lab),
    size = label_size_pt / ggplot2::.pt,
    colour = "black"
  ) +
  theme_void() +
  theme(
    plot.margin   = margin(20, 28, 20, 28, unit = "mm"),
    panel.spacing = unit(0, "mm")
  ) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_equal(clip = "off") +
  scale_shape_manual(
    name   = "Recurrence",
    values = c("1/5" = 22, "2/5" = 24, "3/5" = 25, "4/5" = 23, "5/5" = 21)
  ) +
  scale_fill_manual(
    name   = "Recurrence",
    values = c("1/5" = "white", "2/5" = "white", "3/5" = "white", "4/5" = "white", "5/5" = "#D32F2F")
  )

##### 8.6) Save figures (SVG + PNG; transparent) ##############################
# 8.6.1) Filenames
f_basic_svg <- file.path(out_dir, "circular_dendrogram_motifs_basic2.svg")
f_basic_png <- file.path(out_dir, "circular_dendrogram_motifs_basic2.png")
f_svg       <- file.path(out_dir, "circular_dendrogram_motifs2.svg")
f_png       <- file.path(out_dir, "circular_dendrogram_motifs2.png")

# 8.6.2) Save BASIC
ggplot2::ggsave(filename = f_basic_svg, plot = p_circ_basic, device = "svg",
                width = 22, height = 22, units = "cm", bg = "transparent")
ggplot2::ggsave(filename = f_basic_png, plot = p_circ_basic,
                width = 4500, height = 4500, units = "px", dpi = 320,
                bg = "transparent", limitsize = FALSE)

# 8.6.3) Save WITH RECURRENCE
ggplot2::ggsave(filename = f_svg, plot = p_circ, device = "svg",
                width = 22, height = 22, units = "cm", bg = "transparent")
ggplot2::ggsave(filename = f_png, plot = p_circ,
                width = 4500, height = 4000, units = "px", dpi = 320,
                bg = "transparent", limitsize = FALSE)
