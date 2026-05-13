# Internal helper: resolve genome shortcut to TxDb and OrgDb.
# Separates package names from object names for maintainability rather than
# hard-coding :: constructs in switch branches.
.get_txdb <- function(genome) {
  if (!is.character(genome)) {
    stop("genome must be a character string ('hg38', 'hg19', 'mm10') or 'custom'.", call. = FALSE)
  }

  config <- list(
    hg38 = list(
      txdb_pkg = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      txdb_obj = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      anno_pkg = "org.Hs.eg.db"
    ),
    hg19 = list(
      txdb_pkg = "TxDb.Hsapiens.UCSC.hg19.knownGene",
      txdb_obj = "TxDb.Hsapiens.UCSC.hg19.knownGene",
      anno_pkg = "org.Hs.eg.db"
    ),
    mm10 = list(
      txdb_pkg = "TxDb.Mmusculus.UCSC.mm10.knownGene",
      txdb_obj = "TxDb.Mmusculus.UCSC.mm10.knownGene",
      anno_pkg = "org.Mm.eg.db"
    )
  )

  cfg <- config[[genome]]
  if (is.null(cfg)) {
    stop("Unsupported genome shortcut. Use 'hg38', 'hg19', 'mm10', or 'custom'.", call. = FALSE)
  }

  if (!requireNamespace(cfg$txdb_pkg, quietly = TRUE) ||
      !requireNamespace(cfg$anno_pkg, quietly = TRUE))
    stop(sprintf("Install '%s' and '%s'.", cfg$txdb_pkg, cfg$anno_pkg))

  txdb <- get(cfg$txdb_obj, envir = asNamespace(cfg$txdb_pkg))
  list(txdb = txdb, annoDb = cfg$anno_pkg)
}

#' Stitch Proximal Epigenetic Peaks
#'
#' @description Merges closely spaced peaks into larger continuous macro-domains based on
#' a user-defined stitching distance. This is crucial for analyzing Super-Elements
#' or broad heterochromatin spreading in epiPortrait.
#'
#' @param gr A GRanges object of peaks (typically the output from get_consensus_peaks).
#' @param stitch_distance Numeric. Maximum distance (in base pairs) between peaks
#' to be stitched together. For Super-Elements, 12500 (12.5 kb) is standard.
#' For standard broad marks, you might use 3000 to 5000.
#' @return A GRanges object of stitched domains, with a metadata column
#' 'Constituent_Peaks' recording how many original peaks were merged into each domain.
#'
#' @import GenomicRanges
#' @export
stitch_epi_peaks <- function(gr, stitch_distance = 12500) {

  if (!inherits(gr, "GRanges")) stop("Input 'gr' must be a GRanges object.")

  if (length(gr) == 0) {
    warning("Input GRanges is empty. Returning as is.")
    return(gr)
  }

  if (!is.numeric(stitch_distance) || stitch_distance < 0) {
    stop("Parameter 'stitch_distance' must be a non-negative number.")
  }

  message(sprintf("Stitching peaks within %d bp of each other...", stitch_distance))

  # Core algorithm: use reduce with min.gapwidth for physical stitching
  stitched_gr <- GenomicRanges::reduce(gr, min.gapwidth = stitch_distance)

  # Count how many original peaks are merged into each macro-domain
  mcols(stitched_gr)$Constituent_Peaks <- countOverlaps(stitched_gr, gr)

  message(sprintf("Stitching complete: %d original peaks were stitched into %d continuous domains.",
                  length(gr), length(stitched_gr)))

  return(stitched_gr)
}


#' Filter Peaks by Genomic Annotations (Promoter Exclusion)
#'
#' @description Filters out peaks that overlap with specified genomic regions,
#' typically used to remove promoter-proximal peaks before Super-Element analysis.
#'
#' @param gr A GRanges object of peaks.
#' @param genome A character string ("hg38", "hg19", "mm10", etc.) or a TxDb object.
#' @param upstream Number of bp upstream of TSS to define promoter (default: 2000).
#' @param downstream Number of bp downstream of TSS to define promoter (default: 2000).
#' @return A filtered GRanges object.
#' @import GenomicRanges
#' @export
filter_promoter_peaks <- function(gr, genome = "hg38", upstream = 2000, downstream = 2000) {

  if (length(gr) == 0) return(gr)

  message("Identifying promoter regions for exclusion...")

  # 1. Resolve TxDb (supports shortcuts and custom TxDb objects)
  if (is.character(genome)) {
    txdb <- .get_txdb(genome)$txdb
  } else if (is(genome, "TxDb")) {
    txdb <- genome
  } else {
    stop("Parameter 'genome' must be a valid shortcut string or a TxDb object.")
  }

  # 2. Define promoter regions
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) stop("Please install 'GenomicFeatures'.")
  promoters_gr <- GenomicFeatures::promoters(txdb, upstream = upstream, downstream = downstream)

  # Guard: catch UCSC vs Ensembl chromosome naming mismatches before silent failure
  common_seqs <- intersect(seqlevels(gr), seqlevels(promoters_gr))
  if (length(common_seqs) == 0) {
    warning("NO overlapping chromosome names found between your peaks (e.g., '",
            seqlevels(gr)[1], "') and the reference genome (e.g., '", seqlevels(promoters_gr)[1], "').\n",
            "This usually happens when mixing UCSC ('chr1') and Ensembl ('1') formats. ",
            "Promoter filtering might silently fail to remove any peaks!")
  }

  # 3. Exclude overlapping promoter regions (invert = TRUE keeps non-overlapping)
  filtered_gr <- subsetByOverlaps(gr, promoters_gr, invert = TRUE)

  message(sprintf("Excluded %d peaks overlapping with promoters. %d peaks remaining.",
                  length(gr) - length(filtered_gr), length(filtered_gr)))

  return(filtered_gr)
}





#' Normalize Portrait Assays
#'
#' @description Corrects Intensity and Height assays in a SummarizedExperiment object
#' to account for differences in sequencing depth or distribution across samples.
#'
#' @param se A SummarizedExperiment object from build_portrait_matrix().
#' @param method Normalization method. Options are:
#'   \itemize{
#'     \item \code{"TotalSignal"}: Scales libraries to the mean total signal.
#'     \item \code{"TMM"}: Trimmed Mean of M-values. Designed for count data; use
#'           \code{force_TMM = TRUE} to apply to continuous BigWig signals anyway.
#'     \item \code{"Quantile"}: Forces identical distributions across samples (uses limma).
#'     \item \code{"Z-score"}: Standardizes each peak across samples to mean=0, sd=1 (for clustering, not DE).
#'     \item \code{"None"}: Skips normalization (use if BigWigs are already CPM/Spike-in scaled).
#'   }
#' @param force_TMM Logical. If \code{TRUE}, allows TMM normalization on continuous
#'   BigWig signals despite TMM being designed for count data. Default is \code{FALSE}.
#' @return A normalized SummarizedExperiment object.
#' @import SummarizedExperiment
#' @export
normalize_portrait <- function(se, method = "TotalSignal", force_TMM = FALSE) {

  if (method == "None") {
    message("Method set to 'None'. Skipping normalization (assuming input BigWigs are pre-normalized).")
    return(se)
  }

  valid_methods <- c("TotalSignal", "TMM", "Quantile", "Z-score")
  if (!method %in% valid_methods) {
    stop(sprintf("Invalid method. Choose from: 'None', '%s'", paste(valid_methods, collapse = "', '")))
  }

  message(sprintf("Normalizing 4D features using '%s' method...", method))

  # Extract intensity matrices (Width and Skewness are scale-invariant; skip normalisation)
  int_mat <- assay(se, "Intensity")
  h_mat <- assay(se, "Height")

  if (method == "TotalSignal") {
    sample_sums <- colSums(int_mat, na.rm = TRUE)
    target_scale <- mean(sample_sums)
    scaling_factors <- target_scale / sample_sums

    # Guard against zero-sum samples producing Inf
    bad <- !is.finite(scaling_factors)
    if (any(bad)) {
      warning(sprintf("%d sample(s) have zero total signal; scaling factors set to 1.",
                      sum(bad)))
      scaling_factors[bad] <- 1
    }

    norm_int <- sweep(int_mat, 2, scaling_factors, FUN = "*")
    norm_h <- sweep(h_mat, 2, scaling_factors, FUN = "*")

    colData(se)$ScalingFactor <- scaling_factors

  } else if (method == "TMM") {
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install 'edgeR' to use TMM.")

    if (!force_TMM) {
      stop(
        "TMM normalization is designed for count data and may produce unreliable ",
        "scaling factors on continuous BigWig signals. ",
        "Use 'TotalSignal' or 'Quantile' for continuous data, or set force_TMM = TRUE ",
        "if you understand the risks and still want to proceed."
      )
    }

    warning(
      "TMM normalization was forced on continuous BigWig data. ",
      "Scaling factors may be unreliable. ",
      "Consider TotalSignal or Quantile normalization instead."
    )

    # Use Intensity matrix as a counts proxy for TMM scaling factors
    lib_sizes <- colSums(int_mat, na.rm = TRUE)
    norm_factors <- edgeR::calcNormFactors(int_mat, method = "TMM", lib.size = lib_sizes)

    # Combine library size and TMM factors into final scaling factors
    eff_lib_sizes <- lib_sizes * norm_factors
    target_scale <- mean(eff_lib_sizes)
    scaling_factors <- target_scale / eff_lib_sizes

    norm_int <- sweep(int_mat, 2, scaling_factors, FUN = "*")
    norm_h <- sweep(h_mat, 2, scaling_factors, FUN = "*")

    colData(se)$TMM_NormFactor <- norm_factors
    colData(se)$ScalingFactor <- scaling_factors

  } else if (method == "Quantile") {
    if (!requireNamespace("limma", quietly = TRUE)) stop("Please install 'limma' to use Quantile normalization.")

    # Force distribution alignment; suitable for severe batch effects
    norm_int <- limma::normalizeBetweenArrays(int_mat, method = "quantile")
    norm_h <- limma::normalizeBetweenArrays(h_mat, method = "quantile")

  } else if (method == "Z-score") {
    message("Note: Z-score is recommended for heatmaps/clustering, but NOT for downstream differential testing (limma).")

    # Row-wise (per-peak) standardization: (x - mean) / sd
    norm_int <- t(scale(t(int_mat)))
    norm_h <- t(scale(t(h_mat)))

    # Replace NAs from zero-variance rows (all-zero peaks)
    norm_int[is.na(norm_int)] <- 0
    norm_h[is.na(norm_h)] <- 0
  }

  # Update SummarizedExperiment object
  assay(se, "Intensity") <- norm_int
  assay(se, "Height") <- norm_h

  message(sprintf("Normalization complete! Adjusted matrices written to the SummarizedExperiment object."))
  return(se)
}


#' Identify Super-Element or Broad-Domain
#'
#' @description Implements the ROSE-style ranking algorithm to identify "Super-Element" 
#' (based on Intensity), "Broad-Domain" (based on Width), or "Steep-Peak" (based on Height).
#' It automatically calculates the inflection point of the ranking curve.
#'
#' @param se A SummarizedExperiment object.
#' @param feature Character. The assay to rank: "Intensity" (to identify Super-Elements), 
#' "Width" (to identify Broad-Domains), or "Height" (to identify Steep-Peaks).
#' @param sample_id Character. Which sample to use for ranking? Default is "Mean" 
#' (averages all samples in the SE).
#'
#' @return A SummarizedExperiment object with additional rowData columns: 
#' '{feature}_Domain_Type' and '{feature}_Rank'.
#' @import SummarizedExperiment
#' @importFrom stats setNames
#' @export
call_super_domains <- function(se, feature = "Intensity", sample_id = "Mean") {
  
  if (!feature %in% assayNames(se)) stop("Feature not found in assays.")
  
  mat <- assay(se, feature)
  
  # 1. Compute ranking vector
  if (sample_id == "Mean") {
    val_vector <- rowMeans(mat, na.rm = TRUE)
  } else {
    if (!sample_id %in% colnames(mat)) stop("sample_id not found in samples.")
    val_vector <- mat[, sample_id]
  }
  
  # 2. Sort ascending
  rank_df <- data.frame(
    Peak_ID = rownames(se),
    Value = val_vector,
    Rank = rank(val_vector, ties.method = "first")
  )
  rank_df <- rank_df[order(rank_df$Value), ]
  rank_df$Cumulative_Rank <- seq_len(nrow(rank_df))
  
  # 3. Detect geometric inflection point
  x <- rank_df$Cumulative_Rank
  y <- rank_df$Value
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  dist_to_line <- (x_norm - y_norm) / sqrt(2)
  inflection_idx <- which.max(dist_to_line)
  thresh_val <- y[inflection_idx]
  
  # 4. Map feature to domain-class terminology
  target_label <- switch(feature,
                         "Intensity" = "Super-Element",
                         "Width"     = "Broad-Domain",
                         "Height"    = "Steep-Peak",
                         "Super-Domain") # Fallback
                         
  message(sprintf("Inflection point detected at Rank %d. Identifying %ss based on %s...", 
                  inflection_idx, target_label, feature))
  
  # Label domains
  rank_df$Domain_Type <- ifelse(rank_df$Value >= thresh_val, target_label, "Typical")
  
  # 5. Write results back to SE rowData
  res_vec <- stats::setNames(rank_df$Domain_Type, rank_df$Peak_ID)
  rank_vec <- stats::setNames(rank_df$Cumulative_Rank, rank_df$Peak_ID)
  
  col_name_type <- paste0(feature, "_Domain_Type")
  col_name_rank <- paste0(feature, "_Rank")
  
  rowData(se)[[col_name_type]] <- res_vec[rownames(se)]
  rowData(se)[[col_name_rank]] <- rank_vec[rownames(se)]
  
  return(se)
}


#' Export Shifted Peaks to BED Format
#'
#' @description Extracts significantly shifted peaks (or all analyzed peaks) from the
#' epiPortrait results and writes them to a standard 6-column BED file. This file
#' is ready for downstream tools like HOMER (motif discovery), MEME, or IGV.
#'
#' @param res_df A data.frame generated by \code{test_global_shape_shift}.
#' @param file Character. Path to the output BED file.
#' @param fdr_cutoff Numeric. Only export peaks with adj.P.Val below this threshold (default: 0.05).
#' Set to 1.0 to export all peaks.
#' @param top_n Integer. If specified, only exports the top N most significant peaks.
#'
#' @return Invisible NULL. Writes a file to disk.
#' @importFrom utils write.table
#' @export
export_shifted_bed <- function(res_df, file = "epiPortrait_shifted_peaks.bed",
                               fdr_cutoff = 0.05, top_n = NULL) {

  # 1. Input validation
  req_cols <- c("seqnames", "start", "end", "Peak_ID", "Shape_Shift_Score", "adj.P.Val")
  if (!all(req_cols %in% colnames(res_df))) {
    stop("Input data.frame must be the output from test_global_shape_shift (missing required columns).")
  }

  # 2. Filter by significance threshold
  bed_data <- res_df[res_df$adj.P.Val <= fdr_cutoff, ]

  if (nrow(bed_data) == 0) {
    warning("No peaks found passing the adj.P.Val cutoff. No file written.")
    return(invisible(NULL))
  }

  # 3. Slice top-N if specified (results are already sorted by significance)
  if (!is.null(top_n)) {
    bed_data <- bed_data[seq_len(min(top_n, nrow(bed_data))), ]
  }

  # 4. Build standard 6-column BED: chr, start, end, name, score, strand
  # BED is 0-based; R/Bioc is 1-based -- subtract 1 from start on export
  export_df <- data.frame(
    chrom = bed_data$seqnames,
    start = as.integer(bed_data$start) - 1,
    end   = as.integer(bed_data$end),
    name  = bed_data$Peak_ID,
    score = bed_data$Shape_Shift_Score,
    strand = "."
  )

  # 5. Write tab-delimited file without header
  utils::write.table(export_df, file = file, quote = FALSE, sep = "\t",
                     row.names = FALSE, col.names = FALSE)

  message(sprintf("Success: Exported %d peaks to %s.", nrow(export_df), file))
  return(invisible(NULL))
}
