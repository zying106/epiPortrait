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
#' @import GenomicRanges
#' @examples
#' gr <- GenomicRanges::GRanges("chr1",
#'        IRanges::IRanges(start = c(100, 10000), end = c(200, 10100)))
#' stitch_epi_peaks(gr, stitch_distance = 12500)
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

  # reduce() merges intervals whose gap is
  # STRICTLY LESS than min.gapwidth, so a gap exactly equal to stitch_distance
  # would not merge. Using min.gapwidth = stitch_distance + 1 makes the
  # documented "maximum distance to stitch" (gap <= stitch_distance) exact.
  stitched_gr <- GenomicRanges::reduce(gr, min.gapwidth = stitch_distance + 1)

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
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 5000000), width = 200))
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   filter_promoter_peaks(gr, genome = "hg38")
#' }
#' @export
filter_promoter_peaks <- function(gr, genome = "hg38", upstream = 2000, downstream = 2000) {

  if (length(gr) == 0) return(gr)

  message("Identifying promoter regions for exclusion...")

  res <- .resolve_genome_resources(genome)
  if (is.null(res$txdb)) {
    stop("Promoter filtering requires a TxDb (genome = 'hg38'/'hg19'/'mm10' or a TxDb object).")
  }
  txdb <- res$txdb

  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) stop("Please install 'GenomicFeatures'.")

  # A chr-vs-1 naming mismatch would silently remove 0 promoter overlaps
  # and leave every peak in place; enforce seqlevel compatibility like
  # annotate_epi_domains() does.
  .check_seqlevel_compatibility(gr, txdb = txdb, enforce = TRUE)

  # Validate promoter window parameters using the annotation rules.
  for (nm in c("upstream", "downstream")) {
    v <- get(nm)
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v < 0) {
      stop(sprintf("%s must be a finite non-negative number.", nm))
    }
  }

  promoters_gr <- GenomicFeatures::promoters(txdb, upstream = upstream, downstream = downstream)

  filtered_gr <- IRanges::subsetByOverlaps(gr, promoters_gr, invert = TRUE)

  message(sprintf("Excluded %d peaks overlapping with promoters. %d peaks remaining.",
                  length(gr) - length(filtered_gr), length(filtered_gr)))

  return(filtered_gr)
}


#' Normalize Portrait Assays
#'
#' @description Normalizes the Intensity assay while preserving
#'   SignalDispersion in its native bp-scale units. Corrects the Intensity
#'   matrix to account for differences in sequencing depth or distribution
#'   across samples.
#'
#' @param se A SummarizedExperiment object from build_portrait_matrix().
#' @param method Normalization method. Options are:
#'   \itemize{
#'     \item \code{"None"} (default): Skips normalization. Recommended —
#'           epiPortrait expects BigWigs that are already quantitatively
#'           comparable (same pipeline, genome, normalization; e.g. CPM/RPGC/
#'           spike-in). Post-hoc rescaling risks removing genuine global
#'           biological shifts.
#'     \item \code{"TotalSignal"}: Scales libraries to the mean total signal of
#'           the analyzed domain set. WARNING: this assumes approximately
#'           conserved total signal across samples and may remove real global
#'           gains/losses (e.g. drug-induced chromatin loss).
#'     \item \code{"TMM"}: Trimmed Mean of M-values. Designed for count data;
#'           use \code{force_TMM = TRUE} for continuous BigWig signals.
#'     \item \code{"Quantile"}: Forces identical distributions across samples
#'           (uses limma). Not recommended for between-condition comparison.
#'   }
#'   Row-wise Z-score scaling is not offered here: it is a display/clustering
#'   transform that would destroy the cross-domain magnitude ranking used for
#'   Super calling. Use \code{plot_portrait_pca()} or plotting layers for
#'   display-only scaling.
#' @param force_TMM Logical. If \code{TRUE}, allows TMM normalization on continuous
#'   BigWig signals despite TMM being designed for count data. Default is \code{FALSE}.
#' @return A normalized SummarizedExperiment object.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' normalize_portrait(example_se, method = "None")
#' @export
normalize_portrait <- function(se, method = "None", force_TMM = FALSE) {

  if (method == "None") {
    message("Method set to 'None'. Skipping normalization (assuming input BigWigs are pre-normalized).")
    return(se)
  }

  valid_methods <- c("TotalSignal", "TMM", "Quantile")
  if (!method %in% valid_methods) {
    stop(sprintf("Invalid method. Choose from: 'None', '%s'", paste(valid_methods, collapse = "', '")))
  }

  message(sprintf("Normalizing features using '%s' method...", method))

  int_mat <- assay(se, .resolve_assay(se, "Intensity"))
  disp_mat <- assay(se, .resolve_assay(se, "SignalDispersion"))

  if (method == "TotalSignal") {
    warning(
      "TotalSignal normalization scales samples to the mean total signal of the ",
      "analyzed domain set. It assumes approximately conserved total signal and ",
      "may remove genuine global biological shifts (e.g. drug-induced chromatin ",
      "loss). Consider using pre-normalized BigWigs (method = 'None') for ",
      "between-condition comparison."
    )
    sample_sums <- colSums(int_mat, na.rm = TRUE)
    target_scale <- mean(sample_sums)
    scaling_factors <- target_scale / sample_sums

    bad <- !is.finite(scaling_factors)
    if (any(bad)) {
      warning(sprintf("%d sample(s) have zero total signal; scaling factors set to 1.",
                      sum(bad)))
      scaling_factors[bad] <- 1
    }

    norm_int <- sweep(int_mat, 2, scaling_factors, FUN = "*")
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

    lib_sizes <- colSums(int_mat, na.rm = TRUE)
    norm_factors <- edgeR::calcNormFactors(int_mat, method = "TMM", lib.size = lib_sizes)

    eff_lib_sizes <- lib_sizes * norm_factors
    target_scale <- mean(eff_lib_sizes)
    scaling_factors <- target_scale / eff_lib_sizes

    norm_int <- sweep(int_mat, 2, scaling_factors, FUN = "*")
    colData(se)$TMM_NormFactor <- norm_factors
    colData(se)$ScalingFactor <- scaling_factors

  } else if (method == "Quantile") {
    if (!requireNamespace("limma", quietly = TRUE)) stop("Please install 'limma' to use Quantile normalization.")

    norm_int <- limma::normalizeBetweenArrays(int_mat, method = "quantile")
  }

  assay(se, .resolve_assay(se, "Intensity")) <- norm_int

  # SignalDispersion is not rescaled. It is a spatial measure in bp that
  # is invariant under a uniform multiplicative rescaling of the signal
  # (x_i -> c*x_i leaves the weighted genomic SD unchanged). Applying library
  # scaling factors or quantile transforms would change its unit and break its
  # biological interpretation. Row-wise Z-score scaling is likewise not applied:
  # it would destroy the cross-domain magnitude ranking used for Super calling;
  # plotting / PCA layers apply their own display scaling.
  if (!identical(disp_mat, assay(se, "SignalDispersion"))) {
    stop("SignalDispersion was unexpectedly modified during normalization.")
  }

  message("Normalization complete! Intensity adjusted; SignalDispersion ",
          "(bp-scale spatial descriptor) is left unchanged.")
  return(se)
}

# Run an expression with an OPTIONAL local random seed. When `seed` is NULL
# (default) the expression is evaluated without touching the global RNG state;
# when an integer is supplied it is wrapped in withr::with_seed() so the seed
# is applied only for the duration of `expr` and the caller's RNG stream is
# restored afterwards. This keeps the package free of global set.seed() calls
# (BiocCheck), while still offering reproducible bootstrap / sampling when the
# user explicitly asks for it.
.with_opt_seed <- function(seed, expr) {
  if (is.null(seed)) {
    force(expr)
  } else {
    withr::with_seed(seed, force(expr))
  }
}

# ---- Internal: inflection point detection ----------------------------------
