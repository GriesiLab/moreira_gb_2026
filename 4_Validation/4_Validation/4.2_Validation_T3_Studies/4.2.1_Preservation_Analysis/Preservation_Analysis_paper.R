# ============================================================
# Module Preservation: 5 reference studies vs 11 test studies
# Mixed platforms (RNA-seq VST; Microarray pass-through)
# Cleaning + Resume support
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(WGCNA)
  library(glue)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# -----------------------------
# 0) Working directory (adjust)
# -----------------------------
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022")

# -----------------------------
# 1) Define 'studies'
#    (genes in rows, samples in cols in CSVs)
# -----------------------------
studies <- list(
  Xie_2019_T1 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Fresh/Xie_2019/counts_Final.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Fresh/Xie_2019/metadata.csv", row.names = 1)
  ),
  Anjos_2021_T1 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Fresh/Anjos-Afonso_2020/counts_Final.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Fresh/Anjos-Afonso_2020/metadata.csv", row.names = 1)
  ),
  Xie_2019_T2 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Expanded/Xie_exp_2019/counts_Final.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Expanded/Xie_exp_2019/metadata.csv", row.names = 1)
  ),
  Fares_2014_T2 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Expanded/Sauvageau_2014/counts_Final.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Expanded/Sauvageau_2014/metadata.csv", row.names = 1)
  ),
  Fares_2017_T2 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Expanded/Fares_2017/counts_Final.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Expanded/Fares_2017/metadata.csv", row.names = 1)
  ),
  Expanded_H = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/2_expanded_set/DATA/counts_Final_without_100_RUVg_p06k5.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/2_expanded_set/DATA/metadata_without_100.csv", row.names = 1)
  ),
  Fresh_H = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/1_fresh_set/DATA/counts_Final_RUVg_p07k5 (2).csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/1_fresh_set/DATA/metadata.csv", row.names = 1)
  ),
  Barreyro_2012_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Barreyro_2012/normalized_expression_RMA.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Barreyro_2012/metadata.csv", row.names = 1)
  ),
  Rapin_2014_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Rapin_2014/normalized_expression_RMA.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Rapin_2014/metadata.csv", row.names = 1)
  ),
  Rundberg_2015_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Rundberg_2015/normalized_expression_RMA.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Rundberg_2015/metadata.csv", row.names = 1)
  ),
  Gentles_2010_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Gentles_2010/normalized_expression_RMA.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Gentles_2010/metadata.csv", row.names = 1)
  ),
  Prashad_2014_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Prashad_2014/normalized_expression_RMA.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Prashad_2014/metadata.csv", row.names = 1)
  ),
  Amon_2018_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Fresh/Amon_2018/counts_final.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Fresh/Amon_2018/metadata.csv", row.names = 1)
  ),
  Calvanese_2019_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Fresh/Calvanese_2019/counts_Final.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Fresh/Calvanese_2019/metadata.csv", row.names = 1)
  ),
  Papa_2018_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Expanded/Papa_2018/counts_Final.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/RNA-seq/Expanded/Papa_2018/metadata.csv", row.names = 1)
  ),
  Subramaniam_2019_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Subramaniam_2019/normalized_expression_RMA.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Subramaniam_2019/metadata.csv", row.names = 1)
  ),
  Souyri_2019_T3 = list(
    expData  = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Souyri_2019/normalized_expression_RMA.csv", row.names = 1),
    metadata = read.csv("2c_WGCNA_proteomic_transcriptomic_inHouse/_Victoria_Paiva/analises/validation_insilico/studies_data/Microarray/Souyri_2019/metadata.csv", row.names = 1)
  )
)

# -----------------------------
# 2) Reference and test sets
# -----------------------------
all_names <- names(studies)
stopifnot(length(all_names) >= 17)

ref_names  <- all_names[1:5]
test_names <- all_names[6:17]

message(glue("Ref studies (n={length(ref_names)}): {paste(ref_names, collapse=', ')}"))
message(glue("Test studies (n={length(test_names)}): {paste(test_names, collapse=', ')}"))
message(glue("Total pairwise analyses: {length(ref_names)*length(test_names)}"))

# -----------------------------
# 3) Fake module genes (must have 'Ensembl')
# -----------------------------
fakeModule_path <- "4a_Panel_validation_inSilico/Validantion_Listas_Aleatorias/PPI/Random300_PPI_SelectedSets_10percent.csv"
stopifnot(file.exists(fakeModule_path))
genes_overlap <- read.csv(fakeModule_path, stringsAsFactors = FALSE)

stopifnot("Genes" %in% colnames(genes_overlap))

# Pegue os genes de UMA das listas (ex: a linha 1 ou filtrar pelo rep/set_name)
gene_string <- genes_overlap$Genes[1]  # coloque aqui a lista que vc quer rodar

# Transformar a string em vetor de Ensembl IDs
fakeModule_ensembl_ids <- unique(unlist(strsplit(gene_string, ",")))

fakeModuleColor <- "violet"

# -----------------------------
# 4) Prepare each study
# - RNA-seq (_T1/_T2): VST via DESeq2
# - Microarray (_T3): assume normalized (RMA/log2), pass-through
# Result: studies[[name]]$expr_ready is data.frame(samples x genes)
# -----------------------------
normalize_study <- function(st_name, st_obj) {
  expData  <- st_obj$expData
  metadata <- st_obj$metadata
  
  if (is.null(rownames(metadata))) rownames(metadata) <- colnames(expData)
  metadata <- metadata[colnames(expData), , drop = FALSE]
  
  is_microarray <- grepl("_T3$", st_name)
  
  if (!is_microarray) {
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(expData),
                                  colData   = metadata,
                                  design    = ~ 1)
    vst_data   <- vst(dds, blind = TRUE)
    vst_matrix <- assay(vst_data)            # genes x samples
    expr_df    <- as.data.frame(t(vst_matrix))  # samples x genes
  } else {
    expr_df <- as.data.frame(t(as.matrix(expData)))  # samples x genes
  }
  
  if (any(duplicated(colnames(expr_df)))) {
    warning(glue("[{st_name}] duplicated gene IDs; making them unique."))
    colnames(expr_df) <- make.unique(colnames(expr_df))
  }
  
  st_obj$expr_ready <- expr_df
  st_obj
}

for (nm in names(studies)) {
  message(glue("Preparing study: {nm}"))
  studies[[nm]] <- normalize_study(nm, studies[[nm]])
}

# -----------------------------
# 5) Cleaning helpers for WGCNA
# -----------------------------
clean_expr_for_wgcna <- function(expr_df,
                                 maxNAFrac_gene = 0.5,
                                 maxNAFrac_samp = 0.5,
                                 verbose = 1,
                                 tag = "") {
  stopifnot(is.data.frame(expr_df) || is.matrix(expr_df))
  if (is.matrix(expr_df)) expr_df <- as.data.frame(expr_df)
  
  rn <- rownames(expr_df)
  expr_df <- as.data.frame(data.matrix(expr_df))  # coerce to numeric
  rownames(expr_df) <- rn
  
  if (verbose) message(tag, " Initial genes: ", ncol(expr_df), " | samples: ", nrow(expr_df))
  
  na_frac_gene <- colMeans(!is.finite(as.matrix(expr_df)))
  if (any(na_frac_gene > maxNAFrac_gene)) {
    drop_g <- names(na_frac_gene)[na_frac_gene > maxNAFrac_gene]
    expr_df[, drop_g] <- NULL
    if (verbose) message(tag, " Dropped genes by NA fraction: ", length(drop_g))
  }
  
  na_frac_samp <- rowMeans(!is.finite(as.matrix(expr_df)))
  if (any(na_frac_samp > maxNAFrac_samp)) {
    drop_s <- rownames(expr_df)[na_frac_samp > maxNAFrac_samp]
    expr_df <- expr_df[setdiff(rownames(expr_df), drop_s), , drop = FALSE]
    if (verbose) message(tag, " Dropped samples by NA fraction: ", length(drop_s))
  }
  
  vars <- apply(expr_df, 2, stats::var, na.rm = TRUE)
  zero_var_genes <- names(vars)[!is.finite(vars) | vars == 0]
  if (length(zero_var_genes) > 0) {
    expr_df[, zero_var_genes] <- NULL
    if (verbose) message(tag, " Dropped zero-variance genes: ", length(zero_var_genes))
  }
  
  gsg <- goodSamplesGenes(expr_df, verbose = verbose)
  if (!gsg$allOK) {
    if (verbose) {
      message(tag, " goodSamplesGenes flagged: badGenes=", sum(!gsg$goodGenes),
              " badSamples=", sum(!gsg$goodSamples))
    }
    expr_df <- expr_df[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
  
  stopifnot(nrow(expr_df) >= 4, ncol(expr_df) >= 4)
  expr_df
}

build_color_vector <- function(expr_df, module_genes,
                               fake_color = "violet", bg_color = "grey") {
  genes <- colnames(expr_df)
  if (any(duplicated(genes))) {
    colnames(expr_df) <- make.unique(genes)
    genes <- colnames(expr_df)
  }
  col_vec <- rep(bg_color, length(genes))
  names(col_vec) <- genes
  keep <- intersect(genes, module_genes)
  if (length(keep) == 0) {
    warning("No fake module genes found in this dataset's columns.")
  } else {
    col_vec[keep] <- fake_color
  }
  col_vec
}

# -----------------------------
# 6) One ref–test run with cleaning + resume checks
# -----------------------------
run_preservation_pair <- function(ref_name, test_name,
                                  module_genes,
                                  nPerms = 200,
                                  fake_color = "violet",
                                  out_dir = "Preservation_Analysis/Overlap",
                                  err_log = file.path(out_dir, "preservation_error_log.txt")) {
  obj_name  <- paste0("PresData_", ref_name, "_vs_", test_name)
  out_file  <- file.path(out_dir, paste0(obj_name, ".csv"))
  
  # Resume: skip if already computed
  if (file.exists(out_file)) {
    message(glue("[SKIP] {obj_name} already exists."))
    return(invisible(NULL))
  }
  
  message(glue("\n[Reference]: {ref_name}  ->  [Test]: {test_name}"))
  
  datRef_raw  <- studies[[ref_name]]$expr_ready
  datTest_raw <- studies[[test_name]]$expr_ready
  
  # Try-catch full run to log errors and continue
  res <- try({
    datRef  <- clean_expr_for_wgcna(datRef_raw,  tag = glue("[{ref_name}] "))
    datTest <- clean_expr_for_wgcna(datTest_raw, tag = glue("[{test_name}] "))
    
    colorRef_fake  <- build_color_vector(datRef,  module_genes, fake_color = fake_color, bg_color = "grey")
    colorTest_fake <- build_color_vector(datTest, module_genes, fake_color = fake_color, bg_color = "grey")
    
    if (!any(colorRef_fake == fake_color)) {
      warning(glue("[{ref_name} vs {test_name}] No fake-module genes present after cleaning in REF. Skipping."))
      return(invisible(NULL))
    }
    
    multiExpr  <- list(Reference = list(data = datRef),
                       Test      = list(data = datTest))
    multiColor <- list(Reference = colorRef_fake,
                       Test      = colorTest_fake)
    
    # Sanity
    stopifnot(length(multiColor$Reference) == ncol(multiExpr$Reference$data))
    stopifnot(length(multiColor$Test)      == ncol(multiExpr$Test$data))
    stopifnot(identical(names(multiColor$Reference), colnames(multiExpr$Reference$data)))
    stopifnot(identical(names(multiColor$Test),      colnames(multiExpr$Test$data)))
    
    mp <- modulePreservation(
      multiExpr, multiColor,
      referenceNetworks = 1,
      nPermutations     = nPerms,
      randomSeed        = 1,
      quickCor          = 0,
      verbose           = 3
    )
    
    ref  <- 1
    test <- 2
    modColors <- rownames(mp$preservation$observed[[ref]][[test]])
    
    PresData <- cbind(
      mp$preservation$observed[[ref]][[test]][, 2],       # medianRank
      mp$preservation$Z[[ref]][[test]][, 2:4],            # Zsummary, Zdensity, Zconnectivity
      mp$preservation$log.pBonf[[ref]][[test]][, 2],      # p.Zsummary (log p)
      mp$preservation$observed[[ref]][[test]][, 5],       # propVarExpl
      mp$preservation$Z[[ref]][[test]][, 5],              # Z.propVarExpl
      mp$preservation$log.pBonf[[ref]][[test]][, 5],      # p.Z.propVarExpl (log p)
      mp$preservation$Z[[ref]][[test]][, 8],              # Z.meanAdj
      mp$preservation$log.pBonf[[ref]][[test]][, 8],      # p.Z.meanAdj (log p)
      mp$preservation$observed[[ref]][[test]][, 11],      # cor.kME
      mp$preservation$observed[[ref]][[test]][, 10],      # cor.kIM
      mp$preservation$Z[[ref]][[test]][, 10],             # Z.cor.kIM
      mp$preservation$log.pBonf[[ref]][[test]][, 11],     # p.Z.cor.kIM (log p)
      mp$preservation$observed[[ref]][[test]][, 13],      # cor.cor
      mp$preservation$Z[[ref]][[test]][, 13],             # Z.cor.cor
      mp$preservation$log.pBonf[[ref]][[test]][, 14]      # p.cor.cor (log p)
    )
    
    colnames(PresData) <- c("medianRank","Zsummary","Zdensity","Zconnectivity","p.Zsummary",
                            "propVarExpl","Z.propVarExpl","p.Z.propVarExpl",
                            "Z.meanAdj","p.Z.meanAdj",
                            "cor.kME","cor.kIM","Z.cor.kIM","p.Z.cor.kIM",
                            "cor.cor","Z.cor.cor","p.cor.cor")
    rownames(PresData) <- modColors
    
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(PresData, out_file, quote = FALSE)
    invisible(TRUE)
  }, silent = TRUE)
  
  if (inherits(res, "try-error")) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    msg <- paste0(Sys.time(), " | ", obj_name, " | ERROR: ",
                  conditionMessage(attr(res, "condition")), "\n")
    cat(msg, file = err_log, append = TRUE)
    message(glue("[ERROR logged] {obj_name}"))
  }
  
  invisible(NULL)
}

# -----------------------------
# 7) Run full grid with resume
# -----------------------------
nPerms  <- 200
out_dir <- "Preservation_Analysis/Overlap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------
# 7) Run full grid for ALL random lists
# -----------------------------------
nPerms <- 200
base_out_dir <- "4a_Panel_validation_inSilico/Validantion_Listas_Aleatorias/Preservation_Analysis/Overlap"
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)

# Loop over each random set (each row = one list)
for (i in seq_len(nrow(genes_overlap))) {
  
  # Identify this list
  list_id     <- genes_overlap$Set_name[i]   # ex: "Random_300_Rep415"
  list_rep    <- genes_overlap$Rep[i]        # se você quiser usar também
  gene_string <- genes_overlap$Genes[i]      # string com ENSG separados por vírgula
  
  # Turn the string into a vector of Ensembl IDs
  fakeModule_ensembl_ids <- unique(unlist(strsplit(gene_string, ",")))
  fakeModuleColor <- "violet"
  
  # Create a specific output folder for this list
  out_dir <- file.path(base_out_dir, list_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("\n===============================")
  message("Running preservation for list: ", list_id,
          " (Rep = ", list_rep, 
          ", n_genes = ", length(fakeModule_ensembl_ids), ")")
  message("===============================\n")
  
  # Your original grid: 5 refs × 11 tests
  for (r in ref_names) {
    for (t in test_names) {
      run_preservation_pair(
        ref_name    = r,
        test_name   = t,
        module_genes = fakeModule_ensembl_ids,
        nPerms       = nPerms,
        fake_color   = fakeModuleColor,
        out_dir      = out_dir   # <- agora salva separado por lista
      )
    }
  }
}


# -----------------------------
# 8) Rebuild summary from all existing per-pair results (resume-friendly)
# -----------------------------
summary_path <- file.path(out_dir, "preservation_summary_refs_vs_tests.csv")
csvs <- list.files(out_dir, pattern = "^PresData_.*_vs_.*\\.csv$", full.names = TRUE)

summary_rows <- lapply(csvs, function(fp) {
  df <- tryCatch(read.csv(fp, row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  if (!"violet" %in% rownames(df)) return(NULL)
  
  violet_data <- df["violet", , drop = FALSE]
  bn <- basename(fp)
  parts <- sub("^PresData_", "", sub("\\.csv$", "", bn))
  spl <- strsplit(parts, "_vs_")[[1]]
  ref <- spl[1]; test <- spl[2]
  
  # Module size: count violet genes on reference side (approx from column names that matched fake list)
  # We can’t reconstruct color vector here without data; keep NA if not available
  data.frame(
    ref_study     = ref,
    test_study    = test,
    Module        = "violet",
    medianRank    = as.numeric(violet_data[1, "medianRank"]),
    Zsummary      = as.numeric(violet_data[1, "Zsummary"]),
    Zdensity      = as.numeric(violet_data[1, "Zdensity"]),
    Zconnectivity = as.numeric(violet_data[1, "Zconnectivity"]),
    moduleSize    = NA_real_,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df) > 0) {
  write.csv(summary_df, summary_path, row.names = FALSE)
  message(glue("Done. Summary saved to: {summary_path}"))
} else {
  warning("Finished, but no 'violet' module rows were found. Check ID concordance or logs.")
}

message("If needed, see error log at: ", file.path(out_dir, "preservation_error_log.txt"))
