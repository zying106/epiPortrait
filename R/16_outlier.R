#' Detect 4D Conformational Outliers via Mahalanobis Distance
#'
#' @description Identifies genomic domains whose 4D geometric conformation
#' (Width, Intensity, Height, Skewness) deviates from the population distribution
#' using the Mahalanobis distance. Supports both classic (sample covariance) and
#' robust (Minimum Covariance Determinant) estimators. Unlike MANOVA-based shape
#' shift testing, this does not require group labels — it detects intrinsically
#' unusual conformations, making it suitable for single-sample analysis and QC.
#'
#' @param se A \code{SummarizedExperiment} from \code{build_portrait_matrix()}.
#' @param method Character. Covariance estimator: \code{"robust"} (MCD, default)
#'   or \code{"classic"} (sample mean and covariance).
#' @param ref_samples Character vector or NULL. Sample names to use as the
#'   reference population. If \code{NULL} (default), uses all samples. Supplying
#'   a subset enables query-vs-reference analysis (e.g., normal cohort as
#'   reference, tumour sample as query).
#' @param collapse_fun Function. How to collapse a domain's values across
#'   samples before computing distances. Default \code{rowMeans} (with
#'   \code{na.rm = TRUE}). Use \code{rowMedians} for skewed data, or a custom
#'   function that takes a matrix row and returns a scalar.
#'
#' @return A \code{data.frame} with columns:
#'   \item{seqnames, start, end, Peak_ID}{Genomic coordinates.}
#'   \item{Mahalanobis_D2}{Squared Mahalanobis distance.}
#'   \item{P.Value}{P-value from χ²(4).}
#'   \item{adj.P.Val}{BH-adjusted P-value.}
#'   \item{Outlier_Flag}{Logical; TRUE if adj.P.Val < 0.05.}
#'   Sorted by ascending adj.P.Val.
#'
#' @details
#' The four features are first collapsed per-domain across samples (via
#' \code{collapse_fun}) to produce an N×4 reference matrix. The covariance
#' structure is then estimated using either the classic MLE or the MCD
#' (Minimum Covariance Determinant) robust estimator. The Mahalanobis
#' distance of each domain to the centroid is computed, and a P-value is
#' derived from the χ²(4) distribution. For the robust estimator, a
#' correction factor based on the MCD consistency factor is applied.
#'
#' When \code{ref_samples} is provided, the covariance structure and
#' centroid are estimated from the reference subset, but distances are
#' computed for all domains. This enables single-sample or single-group
#' abnormality detection against a known reference distribution.
#'
#' @import SummarizedExperiment
#' @importFrom stats pchisq p.adjust cov mahalanobis
#' @export
detect_conformational_outliers <- function(se,
                                            method = c("robust", "classic"),
                                            ref_samples = NULL,
                                            collapse_fun = NULL) {

  method <- match.arg(method)

  required_assays <- c("Width", "Intensity", "Height", "Skewness")
  if (!all(required_assays %in% assayNames(se))) {
    stop("SummarizedExperiment must contain all 4 geometric assays.")
  }

  if (is.null(collapse_fun)) {
    collapse_fun <- function(x) mean(x, na.rm = TRUE)
  }

  meta <- as.data.frame(colData(se))
  all_samples <- colnames(se)

  # ---- Resolve reference samples -----------------------------------------
  if (!is.null(ref_samples)) {
    unknown <- setdiff(ref_samples, all_samples)
    if (length(unknown) > 0) {
      stop(sprintf("ref_samples not found in se: %s",
                   paste(unknown, collapse = ", ")))
    }
    ref_idx <- match(ref_samples, all_samples)
    message(sprintf("Using %d reference sample(s): %s",
                    length(ref_samples), paste(ref_samples, collapse = ", ")))
  } else {
    ref_idx <- seq_len(ncol(se))
    message(sprintf("Using all %d samples as reference.", ncol(se)))
  }

  # ---- Build collapsed 4D matrix (N domains × 4 features) ----------------
  features <- c("Width", "Intensity", "Height", "Skewness")
  n_peaks <- nrow(se)

  ref_mat <- vapply(features, function(f) {
    mat <- assay(se, f)[, ref_idx, drop = FALSE]
    apply(mat, 1, collapse_fun)
  }, numeric(n_peaks))

  colnames(ref_mat) <- features

  # Filter rows with NA or zero variance in reference
  row_vars <- apply(ref_mat, 1, stats::var, na.rm = TRUE)
  valid_rows <- !is.na(row_vars) & row_vars > 1e-12
  if (any(!valid_rows)) {
    message(sprintf("Excluding %d domain(s) with NA or zero variance in reference.",
                    sum(!valid_rows)))
  }

  ref_mat_valid <- ref_mat[valid_rows, , drop = FALSE]
  n_valid <- nrow(ref_mat_valid)

  if (n_valid < 8) {
    stop(sprintf("Only %d valid domains; need at least 8 for reliable covariance estimation.", n_valid))
  }

  # ---- Covariance estimation ----------------------------------------------
  if (method == "robust") {
    if (!requireNamespace("robustbase", quietly = TRUE)) {
      stop("Package 'robustbase' is required for robust (MCD) estimation. ",
           "Install it with: install.packages('robustbase')")
    }
    message(sprintf("Computing MCD on %d domains × 4 features...", n_valid))
    mcd_fit <- robustbase::covMcd(ref_mat_valid, alpha = 0.75)
    center  <- mcd_fit$center
    cov_mat <- mcd_fit$cov
    # Consistency correction for chi-squared approximation
    # MCD consistency factor for h/n ≈ 0.75, p = 4
    cons_factor <- mcd_fit$cnp2
    message(sprintf("MCD consistency factor: %.3f", cons_factor))
  } else {
    message("Computing classic (MLE) covariance...")
    center  <- colMeans(ref_mat_valid, na.rm = TRUE)
    cov_mat <- stats::cov(ref_mat_valid)
    cons_factor <- 1.0
  }

  # ---- Mahalanobis distance for all domains --------------------------------
  D2 <- rep(NA_real_, n_peaks)
  D2[valid_rows] <- stats::mahalanobis(ref_mat_valid, center, cov_mat)

  # Apply consistency correction for robust estimator
  D2 <- D2 / cons_factor

  # P-values from chi-squared(4)
  pvals <- rep(NA_real_, n_peaks)
  pvals[valid_rows] <- stats::pchisq(D2[valid_rows], df = 4, lower.tail = FALSE)

  adj_pvals <- stats::p.adjust(pvals, method = "BH")

  # ---- Assemble output ----------------------------------------------------
  gr_df <- as.data.frame(rowRanges(se))
  if (!"seqnames" %in% colnames(gr_df)) colnames(gr_df)[1] <- "seqnames"

  res <- data.frame(
    seqnames        = gr_df$seqnames,
    start           = gr_df$start,
    end             = gr_df$end,
    Peak_ID         = rownames(se),
    Mahalanobis_D2  = round(D2, 4),
    P.Value         = pvals,
    adj.P.Val       = adj_pvals,
    Outlier_Flag    = !is.na(adj_pvals) & adj_pvals < 0.05,
    stringsAsFactors = FALSE
  )

  res <- res[order(res$adj.P.Val, -res$Mahalanobis_D2, na.last = TRUE), ]
  rownames(res) <- NULL

  n_outliers <- sum(res$Outlier_Flag, na.rm = TRUE)
  message(sprintf("Done. %d conformational outlier(s) detected (adj.P.Val < 0.05, %d domains tested).",
                  n_outliers, n_peaks))

  return(res)
}
