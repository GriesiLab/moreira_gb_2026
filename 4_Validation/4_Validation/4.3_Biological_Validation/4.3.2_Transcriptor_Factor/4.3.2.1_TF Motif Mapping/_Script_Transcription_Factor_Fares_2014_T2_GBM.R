############################################
# TF Motif Mapping & Enrichment — R Pipeline
# Author: Gustavo 
# Purpose: Map known TF motifs to promoters of a target gene set and test
#          motif enrichment against a well-matched background.
############################################

##### 0. Loading the libraries ------------------------------------------------------------
# Load all necessary libraries for data wrangling, genomic annotation, 
# promoter extraction, and motif analysis.

# Data manipulation and tidy utilities
library(dplyr)          # data manipulation (select, filter, mutate, etc.)
library(readr)          # fast reading/writing of tabular data (CSV, TSV)
library(stringr)        # string/text manipulation
library(tibble)         # modern data.frame structure
library(purrr)          # functional programming (map, reduce, etc.)

# Genomic data structures and annotations
library(GenomicRanges)   # representation of genomic intervals and ranges
library(GenomicFeatures) # handling transcript and gene models
library(IRanges)         # range-based operations, core for GenomicRanges
library(BSgenome.Hsapiens.UCSC.hg38) # human genome build hg38 (DNA sequences)
library(EnsDb.Hsapiens.v86) # Ensembl gene annotations (GRCh38, Ensembl v86)
library(biomaRt)         # querying Ensembl BioMart (genes, transcripts, etc.)
library(Biostrings)      # efficient manipulation of DNA/RNA/protein sequences

# Motif databases and tools
library(JASPAR2022)      # JASPAR 2022 motif collection (TF binding profiles)
library(TFBSTools)       # handling PFMs/PWMs and motif scanning
library(universalmotif)  # unified motif representation across packages
library(motifmatchr)     # fast motif scanning and matching to sequences
library(PWMEnrich)       # motif enrichment analysis framework
library(PWMEnrich.Hsapiens.background) # precomputed human motif backgrounds


##### 1. Definition of gene sets (Foreground - FG and Background - BG) ---------------------
# [Purpose]
# Define Foreground (FG) = candidate genes of interest, and 
# Background (BG) = all expressed genes in the study.
# Intersection ensures all FG genes belong to the BG universe.

# 1.1 Set working directory
setwd("~/Insync/cgustavobm2018@gmail.com/Google Drive - Shared with me/_Projeto_Pi_Fapesp_2022/2a_b_WGCNA_public_data_fresh_cultivated_HSC")

# 1.2 Load the background expression dataset (all expressed genes)
expData <- read.csv("Type_2_Studies/Sauvageau_2014/DATA/counts_NonNormalized.csv")

# 1.3 Extract Ensembl IDs for all expressed genes (background universe)
genes <- as.character(expData[[1]])

# --- QC logs ---
cat("Class of 'genes':", class(genes), "\n")        # Should be "character"
cat("Total number of background genes:", length(genes), "\n")
cat("Duplicated background genes:", sum(duplicated(genes)), "\n")

# 1.4 Load the foreground gene list (candidate genes of interest)
positive_module_genes <- read.csv("Overlap_Uncult_Cult/2_overlap_Fresh_cult_kme_info_all.csv")

# 1.5 Filter candidate genes present in at least 4 or 5 studies
positive_module_genes_filter <- positive_module_genes[positive_module_genes$CountTotal %in% c(4, 5), ]

# 1.6 Extract Ensembl IDs from the filtered foreground list
positive_gene_list <- positive_module_genes_filter$Ensembl

# --- QC logs ---
cat("Class of 'positive_gene_list':", class(positive_gene_list), "\n")
cat("Total number of candidate genes:", length(positive_gene_list), "\n")
cat("Duplicated candidate genes:", sum(duplicated(positive_gene_list)), "\n")

# 1.7 Intersect FG and BG (retain only candidates present in the background universe)
genes_intersect <- intersect(genes, positive_gene_list)

# --- QC logs ---
cat("Candidate genes present in the background:", length(genes_intersect), "\n")
cat("Percentage of candidate genes retained in the background:",
    round((length(genes_intersect) / length(positive_gene_list)) * 100, 2), "%\n")


##### 2. Identification of promoter regions (−1000/+100 bp) -------------------------------
# [Purpose]
# Retrieve TSS coordinates for each gene and construct promoter windows.
# Promoter = [−1000 bp upstream, +100 bp downstream] around TSS.
# This captures core and proximal promoter regions.

# [TSS selection priority]
# 1. MANE Select (most stable and consensus isoform)
# 2. Ensembl canonical transcript
# 3. Fallback = first available transcript (deterministic given ordering)

# [Genome build]
# Ensembl GRCh38 (BioMart) for annotation, UCSC hg38 for sequence extraction.

# 2.1 Connect to Ensembl BioMart (GRCh38 / human)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 38)

# 2.2 Fetch TSS coordinates for Foreground (FG) genes
tss_info <- getBM(
  attributes = c("ensembl_gene_id",
                 "chromosome_name",
                 "start_position", 
                 "end_position",
                 "strand",
                 "transcription_start_site",
                 "transcript_mane_select",  
                 "transcript_is_canonical"),
  filters = "ensembl_gene_id",
  values  = genes_intersect,
  mart    = ensembl
)

# --- QC logs ---
cat("Class of 'tss_info':", class(tss_info), "\n")
cat("Dimensions (rows, cols):", dim(tss_info), "\n")
cat("Returned columns:", colnames(tss_info), "\n")

# 2.3 Select one TSS per gene (priority: MANE > canonical > fallback)
mane_flag      <- !is.na(tss_info$transcript_mane_select) & nzchar(trimws(tss_info$transcript_mane_select))
canonical_flag <- tss_info$transcript_is_canonical %in% c(1, "1", TRUE)
priority       <- ifelse(mane_flag, 1L, ifelse(canonical_flag, 2L, 3L))

ord        <- order(tss_info$ensembl_gene_id, priority)
tss_ordered <- tss_info[ord, ]
tss_unique  <- tss_ordered[!duplicated(tss_ordered$ensembl_gene_id), ]

# --- QC logs ---
cat("Selected by MANE:     ", sum(!is.na(tss_unique$transcript_mane_select) & nzchar(trimws(tss_unique$transcript_mane_select))), "\n")
cat("Selected by canonical:", sum(is.na(tss_unique$transcript_mane_select) | !nzchar(trimws(tss_unique$transcript_mane_select)) &
                                    (tss_unique$transcript_is_canonical %in% c(1, '1', TRUE))), "\n")
cat("Selected by fallback: ", sum((is.na(tss_unique$transcript_mane_select) | !nzchar(trimws(tss_unique$transcript_mane_select))) &
                                    !(tss_unique$transcript_is_canonical %in% c(1, '1', TRUE))), "\n")
cat("Unique genes with TSS:", length(unique(tss_unique$ensembl_gene_id)), "\n")

# 2.3.4 Create promoter windows: −1000/+100 bp relative to TSS (strand-aware)
promoter_ranges <- promoters(
  GRanges(seqnames = tss_unique$chromosome_name,
          ranges   = IRanges(start = tss_unique$transcription_start_site,
                             end   = tss_unique$transcription_start_site),
          strand   = ifelse(tss_unique$strand == 1, "+", "-")),
  upstream  = 1000,
  downstream = 100
)

# --- QC logs ---
cat("Class of 'promoter_ranges':", class(promoter_ranges), "\n")
cat("Number of promoters created:", length(promoter_ranges), "\n")

# 2.3.5 Harmonize chromosome naming to UCSC style (chr1, chrX, chrM)
seqlevels(promoter_ranges) <- paste0("chr", seqlevels(promoter_ranges))
seqlevels(promoter_ranges) <- gsub("chrMT", "chrM", seqlevels(promoter_ranges))

# 2.3.6 Extract promoter sequences from hg38 genome
promoter_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, promoter_ranges)

# --- QC logs ---
cat("Number of sequences extracted:", length(promoter_seqs), "\n")
cat("Mean sequence length:", mean(width(promoter_seqs)), "bp\n")
cat("Class of 'promoter_seqs':", class(promoter_seqs), "\n")

# 2.3.7 Save FG promoter sequences to FASTA
names(promoter_seqs) <- tss_unique$ensembl_gene_id
writeXStringSet(promoter_seqs, "Type_2_Studies/Sauvageau_2014/RESULTS/Transcription _Factor/promoters_foreground.fa")


# ---------------------------- BACKGROUND (BG) -------------------------------------------
# Repeat the same TSS → promoter pipeline for the entire background gene universe.

tss_bg <- getBM(
  attributes = c("ensembl_gene_id",
                 "chromosome_name",
                 "start_position",
                 "end_position",
                 "strand",
                 "transcription_start_site",
                 "transcript_mane_select",
                 "transcript_is_canonical"),
  filters = "ensembl_gene_id",
  values  = genes,   # all expressed genes = background universe
  mart    = ensembl
)

# --- QC logs ---
cat("Class of 'tss_bg':", class(tss_bg), "\n")
cat("Dimensions (rows, cols):", dim(tss_bg), "\n")
cat("Unique genes in tss_bg:", length(unique(tss_bg$ensembl_gene_id)), "\n")

# Apply same priority: MANE > canonical > fallback
mane_flag_bg      <- !is.na(tss_bg$transcript_mane_select) & nzchar(trimws(tss_bg$transcript_mane_select))
canonical_flag_bg <- tss_bg$transcript_is_canonical %in% c(1, "1", TRUE)
priority_bg       <- ifelse(mane_flag_bg, 1L, ifelse(canonical_flag_bg, 2L, 3L))

ord_bg        <- order(tss_bg$ensembl_gene_id, priority_bg)
tss_bg_ordered <- tss_bg[ord_bg, ]
tss_bg_unique  <- tss_bg_ordered[!duplicated(tss_bg_ordered$ensembl_gene_id), ]

# --- QC logs ---
cat("BG selected by MANE:     ",
    sum(!is.na(tss_bg_unique$transcript_mane_select) & nzchar(trimws(tss_bg_unique$transcript_mane_select))), "\n")
cat("BG selected by canonical:",
    sum((is.na(tss_bg_unique$transcript_mane_select) | !nzchar(trimws(tss_bg_unique$transcript_mane_select))) &
          (tss_bg_unique$transcript_is_canonical %in% c(1, '1', TRUE))), "\n")
cat("BG selected by fallback: ",
    sum((is.na(tss_bg_unique$transcript_mane_select) | !nzchar(trimws(tss_bg_unique$transcript_mane_select))) &
          !(tss_bg_unique$transcript_is_canonical %in% c(1, '1', TRUE))), "\n")
cat("Unique BG genes with TSS:", length(unique(tss_bg_unique$ensembl_gene_id)), "\n")

# Create BG promoter windows (−1000/+100 bp, strand-aware)
promoter_ranges_bg <- promoters(
  GRanges(seqnames = tss_bg_unique$chromosome_name,
          ranges   = IRanges(start = tss_bg_unique$transcription_start_site,
                             end   = tss_bg_unique$transcription_start_site),
          strand   = ifelse(tss_bg_unique$strand == 1, "+", "-")),
  upstream  = 1000,
  downstream = 100
)

# --- QC logs ---
cat("Class of 'promoter_ranges_bg':", class(promoter_ranges_bg), "\n")
cat("Number of BG promoters created:", length(promoter_ranges_bg), "\n")

# Harmonize chromosome naming
seqlevels(promoter_ranges_bg) <- paste0("chr", seqlevels(promoter_ranges_bg))
seqlevels(promoter_ranges_bg) <- gsub("chrMT", "chrM", seqlevels(promoter_ranges_bg))

# Extract BG promoter sequences from hg38
promoter_seqs_bg <- getSeq(BSgenome.Hsapiens.UCSC.hg38, promoter_ranges_bg)

# --- QC logs ---
cat("BG sequences - Class:", class(promoter_seqs_bg), "\n")
cat("BG sequences - N:", length(promoter_seqs_bg), "\n")
cat("BG sequences - mean length:", mean(width(promoter_seqs_bg)), "bp\n")

# Save BG promoter sequences to FASTA
names(promoter_seqs_bg) <- tss_bg_unique$ensembl_gene_id
writeXStringSet(promoter_seqs_bg, "Type_2_Studies/Sauvageau_2014/RESULTS/Transcription _Factor/promoters_background.fa")

##### 3. Preparation and GC-matched Background (FG vs BG) ---------------------------------
# [Purpose]
# Build a GC-matched background so that FG and BG have comparable GC distributions,
# removing GC-driven motif bias. We filter sequences symmetrically, bin by FG GC%,
# and sample BG at a fixed ratio within each GC stratum.

# [Defaults for this study]
# - Symmetric 'N' filter: <= 30% N
# - GC strata: FG-based quintiles (5 bins)
# - BG sampling ratio: 5:1 (BG:FG) per bin
# - Expand GC range minimally if a bin lacks enough BG (±2.5pp, then ±5pp)
# - Verify GC equivalence via KS test (expect high p-values)

# --- Helpers -----------------------------------------------------------------------------

# Strip Ensembl version suffix (e.g., ENSG...12 -> ENSG...)
strip_version <- function(x) sub("\\..*$", "", x)

# GC proportion ignoring 'N'
gc_prop <- function(dna_set) {
  af <- alphabetFrequency(dna_set, baseOnly = FALSE)
  if (!("N" %in% colnames(af))) af <- cbind(af, N = 0)
  denom <- rowSums(af[, c("A","C","G","T"), drop = FALSE])
  gc    <- af[, "G"] + af[, "C"]
  as.numeric(gc / pmax(denom, 1L))
}

# Symmetric 'N' filter (apply equally to FG and BG)
filter_by_N <- function(dna_set, n_thresh = 0.30) {
  af <- alphabetFrequency(dna_set, baseOnly = FALSE)
  if (!("N" %in% colnames(af))) af <- cbind(af, N = 0)
  onlyN <- rowSums(af[, c("A","C","G","T"), drop = FALSE]) == 0
  propN <- af[, "N"] / pmax(rowSums(af[, c("A","C","G","T","N"), drop = FALSE]), 1L)
  keep  <- (!onlyN) & (propN <= n_thresh)
  list(dna = dna_set[keep], keep_idx = which(keep), propN = propN)
}

# Bin by FG GC% quantiles
assign_bins_by_quantiles <- function(values, probs) {
  qs <- unique(quantile(values, probs = probs, na.rm = TRUE, type = 7))
  if (length(qs) < 3) {
    # Fallback: evenly spaced breaks (5 bins)
    rng <- range(values, na.rm = TRUE)
    qs <- seq(rng[1], rng[2], length.out = 6)
  }
  bins <- cut(values, breaks = qs, include.lowest = TRUE, right = FALSE)
  list(bins = bins, breaks = qs)
}

# Expand GC interval until enough BG is available; sample k without replacement
expand_and_sample <- function(pool_gc, target_range, expand_steps_pp = c(0, 0.025, 0.05), k) {
  lo0 <- target_range[1]; hi0 <- target_range[2]
  for (e in expand_steps_pp) {
    lo <- max(0, lo0 - e); hi <- min(1, hi0 + e)
    idx <- which(pool_gc >= lo & pool_gc < hi)
    if (length(idx) >= k) {
      pick <- sample(idx, k, replace = FALSE)
      return(list(pick_idx = pick, lo = lo, hi = hi, expanded = e))
    }
  }
  # If still insufficient, take whatever is available in the original range (logged)
  idx <- which(pool_gc >= lo0 & pool_gc < hi0)
  warning(sprintf("Insufficient BG in bin [%.3f, %.3f). Using %d/%d available.", lo0, hi0, length(idx), k))
  list(pick_idx = idx, lo = lo0, hi = hi0, expanded = NA_real_)
}

# Core builder: GC-matched BG at fixed ratio per FG GC bin
build_gc_matched_bg <- function(promoters_fg,
                                promoters_bg,
                                n_thresh = 0.30,
                                probs = seq(0, 1, by = 0.20),        # quintiles (5 bins)
                                expand_steps_pp = c(0, 0.025, 0.05), # ±2.5pp then ±5pp
                                ratio = 5,                           # BG:FG per bin (default)
                                seed = 42) {
  set.seed(seed)
  
  # Names & version hygiene
  if (is.null(names(promoters_fg)) || is.null(names(promoters_bg))) {
    stop("Both promoters_fg and promoters_bg must have names (Ensembl IDs).")
  }
  names(promoters_fg) <- strip_version(names(promoters_fg))
  names(promoters_bg) <- strip_version(names(promoters_bg))
  
  # Symmetric 'N' filter
  fg_filt <- filter_by_N(promoters_fg, n_thresh = n_thresh)
  bg_filt <- filter_by_N(promoters_bg, n_thresh = n_thresh)
  promoters_fg <- fg_filt$dna
  promoters_bg <- bg_filt$dna
  
  message(sprintf("[QC] FG kept: %d | BG kept: %d (N-threshold = %.0f%%)",
                  length(promoters_fg), length(promoters_bg), 100 * n_thresh))
  
  # GC%
  fg_gc <- gc_prop(promoters_fg)
  bg_gc <- gc_prop(promoters_bg)
  
  # Remove BG entries that overlap FG (by Ensembl ID)
  bg_only_idx <- which(!(names(promoters_bg) %in% names(promoters_fg)))
  promoters_bg <- promoters_bg[bg_only_idx]
  bg_gc        <- bg_gc[bg_only_idx]
  
  # Define FG GC bins (quintiles) and apply the same breaks to BG
  bin_info   <- assign_bins_by_quantiles(fg_gc, probs = probs)
  fg_bins    <- bin_info$bins
  breaks     <- bin_info$breaks
  bin_levels <- levels(fg_bins)
  bg_bins    <- cut(bg_gc, breaks = breaks, include.lowest = TRUE, right = FALSE)
  
  # FG counts per bin (to determine BG demand)
  fg_tbl <- dplyr::count(dplyr::tibble(bin = fg_bins), .data$bin, name = "n_fg")
  fg_tbl$bin <- factor(fg_tbl$bin, levels = bin_levels)
  
  picked_bg_idx <- integer(0)
  log_expansion <- dplyr::tibble(bin = character(), lo = numeric(), hi = numeric(),
                                 expanded_pp = numeric(), need = integer(), got = integer())
  
  # Numeric ranges for logs
  bin_ranges <- lapply(seq_len(length(breaks) - 1), function(i) c(breaks[i], breaks[i+1]))
  names(bin_ranges) <- bin_levels
  
  # Sample BG within each FG bin at the chosen ratio
  for (b in bin_levels) {
    n_fg <- dplyr::filter(fg_tbl, .data$bin == b)$n_fg
    if (length(n_fg) == 0 || is.na(n_fg) || n_fg == 0) next
    
    n_needed <- n_fg * ratio
    pool_idx <- which(bg_bins == b)
    
    if (length(pool_idx) >= n_needed) {
      pick <- sample(pool_idx, n_needed, replace = FALSE)
      picked_bg_idx <- c(picked_bg_idx, pick)
      rng <- bin_ranges[[as.character(b)]]
      log_expansion <- dplyr::bind_rows(
        log_expansion,
        dplyr::tibble(bin = as.character(b), lo = rng[1], hi = rng[2],
                      expanded_pp = 0, need = n_needed, got = n_needed)
      )
    } else {
      rng <- bin_ranges[[as.character(b)]]
      res <- expand_and_sample(pool_gc = bg_gc,
                               target_range = rng,
                               expand_steps_pp = expand_steps_pp,
                               k = n_needed)
      picked_bg_idx <- c(picked_bg_idx, res$pick_idx)
      log_expansion <- dplyr::bind_rows(
        log_expansion,
        dplyr::tibble(bin = as.character(b), lo = res$lo, hi = res$hi,
                      expanded_pp = ifelse(is.na(res$expanded), NA_real_, res$expanded),
                      need = n_needed, got = length(res$pick_idx))
      )
    }
  }
  
  # Final GC-matched BG
  bg_matched <- promoters_bg[picked_bg_idx]
  
  # QC: GC distribution equivalence (KS test)
  fg_gc_final <- gc_prop(promoters_fg)
  bg_gc_final <- gc_prop(bg_matched)
  ks <- suppressWarnings(ks.test(fg_gc_final, bg_gc_final))
  
  message(sprintf("[QC] GC%% FG mean=%.3f (sd=%.3f) | BG mean=%.3f (sd=%.3f) | KS p=%.3g",
                  mean(fg_gc_final), sd(fg_gc_final), mean(bg_gc_final), sd(bg_gc_final), unname(ks$p.value)))
  
  list(
    fg            = promoters_fg,
    bg_matched    = bg_matched,
    fg_gc         = fg_gc_final,
    bg_gc         = bg_gc_final,
    ks            = ks,
    fg_bin_table  = fg_tbl,
    expansion_log = log_expansion,
    breaks        = breaks,
    ratio         = ratio,
    seed          = seed
  )
}

# --- Build GC-matched BG (5:1; FG-defined quintiles) -------------------------------------
bgm <- build_gc_matched_bg(
  promoters_fg    = promoter_seqs,
  promoters_bg    = promoter_seqs_bg,
  n_thresh        = 0.30,
  probs           = seq(0, 1, by = 0.20),  # 5 GC bins (quintiles)
  expand_steps_pp = c(0, 0.025, 0.05),     # minimal expansion if needed
  ratio           = 5,                     # BG:FG per bin
  seed            = 42
)

# QC logs: KS test and expansion summary (inspect for large expansions or shortfalls)
print(bgm$ks)              # KS test FG vs BG GC%
print(bgm$expansion_log)   # per-bin expansion; expect expanded_pp = 0 in most bins

# 3.1 Expose GC-matched FG and BG sequences for downstream motif scanning
promoters_fg <- bgm$fg
promoters_bg <- bgm$bg_matched

# --- QC logs ---
cat("Foreground - Class:", class(promoters_fg), " | # sequences:", length(promoters_fg), "\n")
cat("Background - Class:", class(promoters_bg), " | # sequences:", length(promoters_bg), "\n")
cat("BG:FG ratio realized ≈", round(length(promoters_bg) / length(promoters_fg), 3), "\n")

##### 4. Loading and conversion of JASPAR motifs (PFM → PWM) ------------------------------

# 4.1 Load human JASPAR CORE PFMs
#    - species = 9606 (Homo sapiens)
#    - collection = "CORE" (curated high-quality TF motifs)
#    - all_versions = FALSE (keep only the most recent version per motif ID)
# Note: PFMs = position frequency matrices (counts of A/C/G/T at each position)
opts_human <- list(species = 9606, collection = "CORE", all_versions = FALSE)
pfm_human  <- TFBSTools::getMatrixSet(JASPAR2022, opts_human)

# Basic sanity check: number of motifs loaded should be > 0
cat("[JASPAR] Loaded human PFMs:", length(pfm_human), "\n")
stopifnot(length(pfm_human) > 0)

# 4.2 Prepare container for PWMs (position weight matrices, log-odds form)
# Naming scheme:
#   - use stable JASPAR IDs (e.g., "MA0039.4") as list names
#   - each entry will hold a PWMatrix object after conversion
pwm_list <- vector("list", length(pfm_human))
names(pwm_list) <- sapply(seq_along(pfm_human), function(i) ID(pfm_human[[i]]))

# 4.3 Convert PFMs to PWMs (log-odds matrices)
#    - toPWM(): computes log2( freq / bg_freq ) per base
#    - default pseudocount ≈ 0.8–1.0 → prevents log(0) issues
#    - unless you have evidence of bias, keep the default
ok_pwm <- 0L
for (i in seq_along(pfm_human)) {
  m <- pfm_human[[i]]
  pwm <- try(TFBSTools::toPWM(m), silent = TRUE)
  if (inherits(pwm, "try-error")) {
    # If conversion fails, drop motif (rare but possible)
    pwm_list[[i]] <- NULL
  } else {
    pwm_list[[i]] <- pwm
    ok_pwm <- ok_pwm + 1L
  }
}

# Drop failed conversions (NULLs) to keep a clean list
pwm_list <- pwm_list[!sapply(pwm_list, is.null)]
cat("[JASPAR] Successfully converted PFMs to PWMs:", ok_pwm, "/", length(pfm_human), "\n")

# 4.4 Extra QC (optional)
# - check that all PWM IDs are unique (expected in JASPAR CORE)
# - duplicate IDs would be unexpected and should be logged
cat("[JASPAR] Unique PWM IDs:", length(unique(names(pwm_list))), "\n")
if (length(unique(names(pwm_list))) != length(pwm_list)) {
  warning("[JASPAR] Duplicate IDs detected in PWM list (unexpected).")
}

##### 5. Scanning motifs in FG and BG promoters + enrichment stats ------------------------

# 5.1 Scanning parameters
# - min_score_rel: relative score threshold for motif hits (80% of the max score by default)
#                  Higher thresholds (85–90%) increase specificity but reduce sensitivity.
# - strand_mode: "*" scans both strands of the DNA sequence.
min_score_rel <- "80%"
strand_mode   <- "*"

# 5.2 Select background set
# If a GC-matched background was built in Step 3, use that.
promoters_bg_use <- promoters_bg  
bg_label <- "BG_GCmatched_5to1"   # label for downstream audit

# 5.3 Basic counts and containers
n_motifs   <- length(pwm_list)     # total number of motifs to test
n_fg_total <- length(promoters_fg) # total number of FG promoters
n_bg_total <- length(promoters_bg_use) # total number of BG promoters

# ---------------- Sanity checks ("guard-rails") ----------------
stopifnot(exists("bgm"))   # GC-matching object must exist
stopifnot(identical(promoters_fg,     bgm$fg))          # FG consistency
stopifnot(identical(promoters_bg_use, bgm$bg_matched))  # BG consistency

# Ensure no overlap between FG and BG gene IDs
stopifnot(length(intersect(names(promoters_fg), names(promoters_bg_use))) == 0)

# Check BG ratio (e.g., 3:1, 5:1, etc.) is consistent
stopifnot(!is.null(bgm$ratio) && bgm$ratio == 5)
stopifnot(length(promoters_bg_use) == sum(bgm$fg_bin_table$n_fg) * bgm$ratio)

# KS test on GC% distributions (FG vs BG) — must be non-significant (p > 0.5)
ks_p <- as.numeric(bgm$ks$p.value)
stopifnot(ks_p > 0.5)

# Log dataset info
cat("[AUDIT] FG:", length(promoters_fg),
    "| BG:", length(promoters_bg_use),
    "| ratio:", bgm$ratio,
    "| KS(GC) p=", signif(ks_p, 4), "\n")

cat("[JASPAR] Motifs to test:", n_motifs, 
    "| FG promoters:", n_fg_total, 
    "| BG promoters:", n_bg_total, 
    "| BG choice:", bg_label, "\n")

# ---------------- Result containers ----------------
# One row per motif
res_jaspar_id   <- character(n_motifs)   # motif ID (e.g., "MA0139.1")
res_tf_name     <- character(n_motifs)   # TF symbol (if annotated)
res_n_fg_hit    <- integer(n_motifs)     # # FG promoters with ≥1 motif hit
res_n_bg_hit    <- integer(n_motifs)     # # BG promoters with ≥1 motif hit
res_p_value     <- numeric(n_motifs)     # Fisher exact test p-value
res_odds_ratio  <- numeric(n_motifs)     # Fisher odds ratio

# Optional: mapping motif ID → TF name from JASPAR
pfm_ids  <- sapply(pfm_human, ID)
pfm_name <- sapply(pfm_human, name)
id2name  <- stats::setNames(pfm_name, pfm_ids)

# ---------------- Main scan loop ----------------
# Strategy:
# - For each PWM, scan all FG and BG promoters.
# - Define "hit" = ≥1 match above threshold in a promoter.
# - Build a 2x2 contingency table and test enrichment (Fisher).
motif_counter <- 0L
for (mot_id in names(pwm_list)) {
  motif_counter <- motif_counter + 1L
  pwm <- pwm_list[[mot_id]]
  
  # Record ID and TF name (if annotated)
  res_jaspar_id[motif_counter] <- mot_id
  res_tf_name[motif_counter]   <- if (!is.na(id2name[mot_id])) id2name[mot_id] else NA_character_
  
  # Count FG hits
  n_hit_fg <- 0L
  for (s in seq_along(promoters_fg)) {
    seq_fg <- promoters_fg[[s]]
    hits_fg <- try(
      TFBSTools::searchSeq(pwm, seq_fg, min.score = min_score_rel, strand = strand_mode),
      silent = TRUE
    )
    if (!inherits(hits_fg, "try-error") && length(hits_fg) > 0) n_hit_fg <- n_hit_fg + 1L
  }
  
  # Count BG hits
  n_hit_bg <- 0L
  for (s in seq_along(promoters_bg_use)) {
    seq_bg <- promoters_bg_use[[s]]
    hits_bg <- try(
      TFBSTools::searchSeq(pwm, seq_bg, min.score = min_score_rel, strand = strand_mode),
      silent = TRUE
    )
    if (!inherits(hits_bg, "try-error") && length(hits_bg) > 0) n_hit_bg <- n_hit_bg + 1L
  }
  
  # Store counts
  res_n_fg_hit[motif_counter] <- n_hit_fg
  res_n_bg_hit[motif_counter] <- n_hit_bg
  
  # Fisher's exact test (one-sided, enrichment in FG)
  # Table:
  #            |   hit        no_hit
  #   FG       |    a           b
  #   BG       |    c           d
  a <- n_hit_fg
  b <- n_fg_total - n_hit_fg
  c <- n_hit_bg
  d <- n_bg_total - n_hit_bg
  
  if (any(c(a, b, c, d) < 0)) {
    res_p_value[motif_counter]    <- NA_real_
    res_odds_ratio[motif_counter] <- NA_real_
  } else {
    tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    ft  <- fisher.test(tab, alternative = "greater") # one-sided: FG enriched
    res_p_value[motif_counter]    <- ft$p.value
    res_odds_ratio[motif_counter] <- as.numeric(ft$estimate)
  }
  
  # Log progress every 50 motifs
  if (motif_counter %% 50 == 0) {
    cat(sprintf("[JASPAR] Processed %d/%d motifs...\n", motif_counter, n_motifs))
  }
}

# ---------------- Build results table ----------------
res_df <- data.frame(
  jaspar_id   = res_jaspar_id,
  tf_name     = res_tf_name,
  n_fg_hit    = res_n_fg_hit,
  n_fg_total  = n_fg_total,
  n_bg_hit    = res_n_bg_hit,
  n_bg_total  = n_bg_total,
  odds_ratio  = res_odds_ratio,
  p_value     = res_p_value,
  stringsAsFactors = FALSE
)

# Multiple testing correction (Benjamini–Hochberg FDR)
#res_df$FDR <- p.adjust(res_df$p_value, method = "BH")

# Calculate hit proportions (optional QC)
res_df$prop_fg <- res_df$n_fg_hit / res_df$n_fg_total
res_df$prop_bg <- res_df$n_bg_hit / res_df$n_bg_total

# Record parameters for reproducibility
res_df$scan_min_score <- min_score_rel
res_df$scan_strand    <- strand_mode
res_df$bg_choice      <- bg_label

# Order results by FDR, then raw p-value
ord <- order(res_df$p_value, res_df$p_value, na.last = TRUE)
res_df <- res_df[ord, ]

cat("[JASPAR] Motifs tested:", nrow(res_df),
    "| Significant (p-Value < 0.05):", sum(res_df$p_value < 0.05, na.rm = TRUE), "\n")

# ---------------- Save outputs ----------------
# Main results
outdir_jas <- "Type_2_Studies/Sauvageau_2014/RESULTS/Transcription _Factor"
if (!dir.exists(outdir_jas)) dir.create(outdir_jas, recursive = TRUE)

outfile <- file.path(
  outdir_jas,
  sprintf("jaspar2022_human_fg_bg_fisher_minScore%s_%s.csv",
          gsub("%","", min_score_rel), bg_label)
)
utils::write.csv(res_df, outfile, row.names = FALSE)
cat("[JASPAR] Results saved at:", outfile, "\n")

# Metadata (parameters + versions) for reproducibility
meta <- list(
  time            = as.character(Sys.time()),
  genome          = "hg38 (UCSC)",
  promoter_window = "-1000/+100 bp",
  fg_n            = n_fg_total,
  bg_n            = n_bg_total,
  bg_choice       = bg_label,
  bg_ratio        = if (!is.null(bgm$ratio)) bgm$ratio else NA,
  ks_gc_pvalue    = if (exists("ks_p")) ks_p else NA,
  min_score_rel   = min_score_rel,
  strand_mode     = strand_mode,
  jaspar_set      = "JASPAR2022 CORE (human)",
  R_version       = R.version.string,
  package_versions = sapply(c("TFBSTools","JASPAR2022","Biostrings","universalmotif"),
                            function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_)
)

meta_file <- file.path(outdir_jas, "run_metadata.txt")
writeLines(capture.output(str(meta)), con = meta_file)
cat("[JASPAR] Metadata saved at:", meta_file, "\n")
