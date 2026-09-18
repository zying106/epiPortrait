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
#' 'Constituent_Peaks' recording how many original peaks were merged into each
#' domain. The call is recorded in \code{S4Vectors::metadata()} under
#' \code{stitch_provenance} (input/output counts, distance, timestamp), and is
#' propagated into \code{metadata(se)$stitch_provenance} by
#' \code{build_portrait_matrix()} when the stitched domains are used there.
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

  # Stitching determines the domain geometry used by every downstream layer,
  # so the call itself must be auditable: record what was merged and how far.
  S4Vectors::metadata(stitched_gr)$stitch_provenance <- list(
    stitch_distance_bp = stitch_distance,
    min_gapwidth_used = stitch_distance + 1,
    n_input_peaks = length(gr),
    n_output_domains = length(stitched_gr),
    rule = "GenomicRanges::reduce with gap <= stitch_distance merged",
    call = match.call(),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE))

  message(sprintf("Stitching complete: %d original peaks were stitched into %d continuous domains.",
                  length(gr), length(stitched_gr)))

  return(stitched_gr)
}


#' Filter Peaks by Genomic Annotations (Promoter Exclusion)
#'
#' @description Filters out peaks that overlap with specified genomic regions,
#' typically used to remove promoter-proximal peaks before Super-Element analysis.
#'
#' This is an explicit, opt-in universe-definition step: the rest of the
#' package never applies it automatically. For distal enhancer marks
#' (\code{get_mark_preset("H3K27ac")$exclude_promoter} is \code{TRUE}), call it
#' on the consensus peaks \emph{before} \code{\link{stitch_epi_peaks}()}, which
#' matches the order of ROSE's \code{-t} option. Note that ROSE's own default is
#' no TSS exclusion (\code{-t 0}); set \code{upstream}/\code{downstream} (e.g.
#' 2500 for a ROSE-like exclusion) only when a promoter-excluded universe is
#' intended, and use a matching external reference.
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
#' @description Optionally normalizes the Intensity assay while preserving
#'   SignalDispersion in its native bp-scale units. The recommended workflow
#'   uses quantitatively comparable BigWigs and \code{method = "None"}.
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
#'     \item \code{"Quantile"}: Forces identical distributions across samples
#'           (uses limma). This is an explicitly warned sensitivity-analysis
#'           option, not a recommended between-condition normalization.
#'   }
#'   TMM is intentionally not offered because it is a count-composition method,
#'   whereas epiPortrait operates on continuous integrated BigWig signal.
#'   Row-wise Z-score scaling is not offered here: it is a display/clustering
#'   transform that would destroy the cross-domain magnitude ranking used for
#'   Super calling. Use \code{plot_portrait_pca()} or plotting layers for
#'   display-only scaling.
#' @return A SummarizedExperiment. The decision and, when applicable, scaling
#'   factors and before/after sample totals are stored in
#'   \code{metadata(se)$normalization}. A second post-hoc normalization of an
#'   object already marked as normalized is rejected.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' normalize_portrait(example_se, method = "None")
#' @export
normalize_portrait <- function(se, method = "None") {

  valid_methods <- c("None", "TotalSignal", "Quantile")
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose from: 'None', 'TotalSignal', or 'Quantile'.")
  }

  if (method == "None") {
    if (is.null(S4Vectors::metadata(se)$normalization)) {
      S4Vectors::metadata(se)$normalization <- list(
        method = "None",
        applied = FALSE,
        input_assay = .resolve_assay(se, "Intensity"),
        assumption = paste(
          "Input BigWigs are already quantitatively comparable",
          "(same processing and CPM/RPGC/spike-in normalization).")
      )
    }
    message("Method set to 'None'. Skipping normalization; input BigWigs are assumed to be quantitatively comparable.")
    return(se)
  }

  previous <- S4Vectors::metadata(se)$normalization
  if (!is.null(previous) && isTRUE(previous$applied)) {
    stop("Intensity is already marked as normalized by normalize_portrait() ",
         "using method = '", previous$method, "'. Return to the pre-normalized ",
         "object instead of applying a second post-hoc normalization.")
  }

  message(sprintf("Normalizing features using '%s' method...", method))

  intensity_assay <- .resolve_assay(se, "Intensity")
  dispersion_assay <- .resolve_assay(se, "SignalDispersion")
  int_mat <- assay(se, intensity_assay)
  disp_mat <- assay(se, dispersion_assay)
  input_totals <- stats::setNames(colSums(int_mat, na.rm = TRUE), colnames(se))
  scaling_factors <- NULL
  assumption <- NULL

  if (method == "TotalSignal") {
    warning(
      "TotalSignal normalization scales samples to the mean total signal of the ",
      "analyzed domain set. It assumes approximately conserved total signal and ",
      "may remove genuine global biological shifts (e.g. drug-induced chromatin ",
      "loss). Consider using pre-normalized BigWigs (method = 'None') for ",
      "between-condition comparison.", call. = FALSE
    )
    sample_sums <- colSums(int_mat, na.rm = TRUE)
    target_scale <- mean(sample_sums)
    scaling_factors <- target_scale / sample_sums

    bad <- !is.finite(scaling_factors)
    if (any(bad)) {
      warning(sprintf("%d sample(s) have zero total signal; scaling factors set to 1.",
                      sum(bad)), call. = FALSE)
      scaling_factors[bad] <- 1
    }

    norm_int <- sweep(int_mat, 2, scaling_factors, FUN = "*")
    colData(se)$ScalingFactor <- scaling_factors
    assumption <- paste(
      "The total signal over the analyzed domain universe is approximately",
      "conserved across samples.")

  } else if (method == "Quantile") {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Please install 'limma' to use Quantile normalization.")
    }
    warning(
      "Quantile normalization forces identical Intensity distributions across ",
      "samples and can erase genuine global biological shifts. Use it only as ",
      "an explicitly justified sensitivity analysis, not as the default for ",
      "between-condition comparisons.", call. = FALSE
    )

    norm_int <- limma::normalizeBetweenArrays(int_mat, method = "quantile")
    assumption <- paste(
      "Differences in marginal Intensity distributions are treated as",
      "technical rather than biological.")
  }

  assay(se, intensity_assay) <- norm_int

  # SignalDispersion is not rescaled. It is a spatial measure in bp that
  # is invariant under a uniform multiplicative rescaling of the signal
  # (x_i -> c*x_i leaves the weighted genomic SD unchanged). Applying library
  # scaling factors or quantile transforms would change its unit and break its
  # biological interpretation. Row-wise Z-score scaling is likewise not applied:
  # it would destroy the cross-domain magnitude ranking used for Super calling;
  # plotting / PCA layers apply their own display scaling.
  if (!identical(disp_mat, assay(se, dispersion_assay))) {
    stop("SignalDispersion was unexpectedly modified during normalization.")
  }

  S4Vectors::metadata(se)$normalization <- list(
    method = method,
    applied = TRUE,
    input_assay = intensity_assay,
    transformed_assay = intensity_assay,
    signal_type = "continuous integrated BigWig signal",
    scaling_factors = if (is.null(scaling_factors)) NULL else
      stats::setNames(as.numeric(scaling_factors), colnames(se)),
    input_sample_totals = input_totals,
    output_sample_totals = stats::setNames(
      colSums(norm_int, na.rm = TRUE), colnames(se)),
    assumption = assumption,
    SignalDispersion_modified = FALSE
  )

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
