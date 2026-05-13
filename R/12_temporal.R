#' Detect Temporal or Dose-Response Shape Drift
#'
#' @description Tests whether 4D geometric features drift systematically across
#' ordered conditions (timepoints, doses). Uses polynomial contrasts in MANOVA
#' to partition multivariate variance into linear and quadratic trend components,
#' detecting monotonic broadening/sharpening and non-linear patterns (e.g.,
#' peak-then-collapse) that two-group comparisons miss.
#'
#' @param se A \code{SummarizedExperiment} object from \code{build_portrait_matrix()}.
#' @param time_var Character. Column in \code{colData(se)} with ordered
#'   numeric timepoints or doses (e.g., \code{c(0, 2, 6, 24)}).
#' @param block_var Character or NULL. Optional blocking variable (e.g., patient).
#' @param poly_order Integer. Maximum polynomial order to test. Default 2
#'   (linear + quadratic). Set to 3 to also test cubic trends.
#' @param workers Number of threads for parallel processing (default: 1).
#'
#' @return A \code{data.frame} with columns:
#'   \item{Linear_Score}{Pillai trace for the linear trend.}
#'   \item{Linear_Pval}{P-value for the linear trend.}
#'   \item{Quadratic_Score}{Pillai trace for the quadratic trend (if \code{poly_order >= 2}).}
#'   \item{Quadratic_Pval}{P-value for the quadratic trend.}
#'   ... plus genomic coordinates, adj.P.Val, and Status.
#'
#' @import SummarizedExperiment
#' @import BiocParallel
#' @importFrom stats manova summary.aov p.adjust var contr.poly
#' @export
test_temporal_shape_shift <- function(se,
                                       time_var = "Timepoint",
                                       block_var = NULL,
                                       poly_order = 2,
                                       workers = 1) {

  required_assays <- c("Width", "Intensity", "Height", "Skewness")
  if (!all(required_assays %in% assayNames(se))) {
    stop("SummarizedExperiment must contain all 4 geometric assays.")
  }

  meta_data <- as.data.frame(colData(se))
  if (!time_var %in% colnames(meta_data)) {
    stop(sprintf("Time variable '%s' not found in colData(se).", time_var))
  }

  time_vals <- meta_data[[time_var]]
  if (is.character(time_vals)) time_vals <- as.numeric(time_vals)
  if (!is.numeric(time_vals)) {
    stop(sprintf("Time variable '%s' must be numeric or convertible to numeric.", time_var))
  }

  unique_times <- sort(unique(time_vals))
  if (length(unique_times) < 3) {
    stop("At least 3 distinct timepoints/doses are required for trend testing.")
  }

  # Build ordered factor with polynomial contrasts
  time_ordered <- factor(time_vals, levels = unique_times, ordered = TRUE)
  contrasts(time_ordered) <- stats::contr.poly(length(unique_times))
  # Only keep contrasts up to poly_order
  if (ncol(contrasts(time_ordered)) > poly_order) {
    contrasts(time_ordered) <- contrasts(time_ordered)[, 1:poly_order, drop = FALSE]
  }

  if (!is.null(block_var)) {
    if (!block_var %in% colnames(meta_data)) {
      stop(sprintf("Block variable '%s' not found.", block_var))
    }
    blocks <- as.factor(meta_data[[block_var]])
    model_df <- data.frame(time = time_ordered, block = blocks)
  } else {
    model_df <- data.frame(time = time_ordered)
  }

  num_peaks <- nrow(se)
  contrast_labels <- colnames(contrasts(time_ordered))

  message(sprintf("Testing temporal shape drift for %d peaks across %d timepoints (%s)...",
                  num_peaks, length(unique_times),
                  paste(unique_times, collapse = ", ")))

  # ---- Parallel backend -----------------------------------------------------
  if (workers > 1) {
    param <- if (.Platform$OS.type == "unix") BiocParallel::MulticoreParam(workers) else BiocParallel::SnowParam(workers)
  } else {
    param <- BiocParallel::SerialParam()
  }

  w_mat <- assay(se, "Width")
  i_mat <- assay(se, "Intensity")
  h_mat <- assay(se, "Height")
  s_mat <- assay(se, "Skewness")

  res_list <- BiocParallel::bplapply(seq_len(num_peaks), function(i) {

    Y <- cbind(W = w_mat[i, ], I = i_mat[i, ], H = h_mat[i, ], S = s_mat[i, ])
    col_vars <- apply(Y, 2, stats::var, na.rm = TRUE)
    valid_cols <- col_vars > 1e-8

    if (sum(valid_cols) < 2) {
      out <- c(Score = 0, Pval = 1, Eta2 = 0, Status = "Low_Variance")
      # Add NAs for all contrast columns
      for (cl in contrast_labels) {
        out <- c(out, stats::setNames(c(0, 1, 0),
                paste0(cl, c("_Score", "_Pval", "_Eta2"))))
      }
      return(out)
    }

    Y_valid <- Y[, valid_cols, drop = FALSE]
    k_valid <- sum(valid_cols)

    tryCatch({
      if (is.null(block_var)) {
        fit <- stats::manova(Y_valid ~ time, data = model_df)
      } else {
        fit <- stats::manova(Y_valid ~ time + block, data = model_df)
      }
      res <- summary(fit, test = "Pillai")

      out <- c(Status = "OK")
      for (cl in contrast_labels) {
        term_name <- paste0("time", cl)  # e.g., "time.L", "time.Q"
        if (term_name %in% rownames(res$stats)) {
          score <- res$stats[term_name, "Pillai"]
          pval  <- res$stats[term_name, "Pr(>F)"]
          df_eff <- res$stats[term_name, "Df"]
          s_val <- min(k_valid, df_eff)
          eta2 <- if (s_val > 0) score / s_val else 0
        } else {
          score <- 0; pval <- 1; eta2 <- 0
        }
        out <- c(out, stats::setNames(c(score, pval, eta2),
                paste0(cl, c("_Score", "_Pval", "_Eta2"))))
      }
      return(out)
    }, error = function(e) {
      out <- c(Status = "MANOVA_Error")
      for (cl in contrast_labels) {
        out <- c(out, stats::setNames(c(0, 1, 0),
                paste0(cl, c("_Score", "_Pval", "_Eta2"))))
      }
      return(out)
    })

  }, BPPARAM = param)

  # ---- Assemble results -----------------------------------------------------
  res_mat <- do.call(rbind, res_list)
  res_mat <- as.data.frame(res_mat, stringsAsFactors = FALSE)
  for (col in colnames(res_mat)) {
    if (col != "Status") res_mat[[col]] <- as.numeric(res_mat[[col]])
  }

  gr_df <- as.data.frame(rowRanges(se))
  if (!"seqnames" %in% colnames(gr_df)) colnames(gr_df)[1] <- "seqnames"

  res_df <- data.frame(
    seqnames = gr_df$seqnames,
    start    = gr_df$start,
    end      = gr_df$end,
    Peak_ID  = rownames(se),
    stringsAsFactors = FALSE
  )

  # Bind trend columns
  res_df <- cbind(res_df, res_mat)

  # BH correction on linear P-value
  if ("Linear_Pval" %in% colnames(res_df) || ".L_Pval" %in% colnames(res_df)) {
    pval_col <- if ("Linear_Pval" %in% colnames(res_df)) "Linear_Pval" else ".L_Pval"
    res_df$adj.P.Val <- stats::p.adjust(res_df[[pval_col]], method = "BH")
  }

  res_df <- res_df[order(res_df$adj.P.Val,
                         -res_df[[grep("_Score$", colnames(res_df))[1]]]), ]
  rownames(res_df) <- NULL

  n_sig <- sum(res_df$adj.P.Val < 0.05, na.rm = TRUE)
  message(sprintf("Temporal shape drift complete. %d peaks with significant trend (adj.P.Val < 0.05).", n_sig))

  return(res_df)
}
