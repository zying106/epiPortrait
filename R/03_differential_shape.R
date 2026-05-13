#' Perform Global Multivariate Shape Shift Test
#'
#' @description Evaluates the simultaneous shift in multiple geometric dimensions
#' (Width, Intensity, Height, Skewness) between experimental conditions using MANOVA.
#' Supports both simple two-group comparisons (\code{group_var}) and complex
#' multi-factor designs (\code{formula}).
#'
#' @param se A SummarizedExperiment object from build_portrait_matrix().
#' @param formula A formula specifying the full model design (e.g.,
#'   \code{~ genotype * treatment + batch}). If \code{NULL} (default), falls
#'   back to the simple \code{group_var} + \code{block_var} interface.
#' @param group_var Character. Column in \code{colData(se)} with the grouping
#'   variable. Only used when \code{formula = NULL}.
#' @param block_var Character or NULL. Optional blocking variable for paired
#'   designs. Only used when \code{formula = NULL}.
#' @param term Character. The model term to extract as the Shape_Shift_Score.
#'   Default \code{NULL} auto-selects the last term in the formula (typically
#'   the effect of interest).
#' @param test Character. MANOVA test statistic: \code{"Pillai"} (default),
#'   \code{"Wilks"}, \code{"Hotelling-Lawley"}, or \code{"Roy"}.
#' @param workers Number of threads for parallel processing (default: 1).
#' @return A data.frame containing \code{Shape_Shift_Score}, \code{partial_eta_sq}
#'   (effect size), multivariate P-values, and BH-adjusted P-values.
#' @import SummarizedExperiment
#' @import BiocParallel
#' @importFrom stats manova summary.aov p.adjust var as.formula terms
#' @export
test_global_shape_shift <- function(se,
                                    formula = NULL,
                                    group_var = "Condition",
                                    block_var = NULL,
                                    term = NULL,
                                    test = "Pillai",
                                    workers = 1) {

  required_assays <- c("Width", "Intensity", "Height", "Skewness")
  if (!all(required_assays %in% assayNames(se))) {
    stop("SummarizedExperiment must contain all 4 geometric assays (epiPortrait).")
  }

  meta_data <- as.data.frame(colData(se))

  # ---- Resolve model specification -------------------------------------------
  if (!is.null(formula)) {
    # Formula-based interface
    if (!inherits(formula, "formula")) {
      formula <- stats::as.formula(formula)
    }
    # Build model frame from formula — only include variables referenced in formula
    formula_vars <- all.vars(formula)
    missing_vars <- setdiff(formula_vars, colnames(meta_data))
    if (length(missing_vars) > 0) {
      stop(sprintf("Variable(s) not found in colData(se): %s",
                   paste(missing_vars, collapse = ", ")))
    }
    model_df <- meta_data[, formula_vars, drop = FALSE]
    # Ensure factor variables are treated as such
    for (v in formula_vars) {
      if (is.character(model_df[[v]])) model_df[[v]] <- as.factor(model_df[[v]])
    }
    # Rebuild formula with the clean model frame
    formula <- stats::as.formula(paste("~", paste(formula_vars, collapse = " + ")))
    # Auto-detect the term of interest (last term in formula)
    if (is.null(term)) {
      term_labels <- attr(stats::terms(formula), "term.labels")
      term <- term_labels[length(term_labels)]
      message(sprintf("Auto-selected term '%s' as the effect of interest.", term))
    }
    # Warn on sample size relative to model complexity
    n_coef <- length(attr(stats::terms(formula), "term.labels"))
    if (nrow(model_df) < n_coef + 4) {
      warning(sprintf("Only %d samples for %d model terms. MANOVA may be underpowered.",
                      nrow(model_df), n_coef))
    }
  } else {
    # Simple group_var + block_var interface (backward-compatible)
    if (!group_var %in% colnames(meta_data)) {
      stop(sprintf("Group variable '%s' not found in sample metadata.", group_var))
    }
    groups <- as.factor(meta_data[[group_var]])
    if (length(levels(groups)) < 2) {
      stop("MANOVA requires at least two distinct groups.")
    }
    if (!is.null(block_var)) {
      if (!block_var %in% colnames(meta_data)) {
        stop(sprintf("Block variable '%s' not found.", block_var))
      }
      blocks <- as.factor(meta_data[[block_var]])
      formula <- stats::as.formula("~ groups + blocks")
      model_df <- data.frame(groups = groups, blocks = blocks)
    } else {
      blocks <- NULL
      formula <- stats::as.formula("~ groups")
      model_df <- data.frame(groups = groups)
    }
    term <- "groups"

    min_reps <- min(table(groups))
    if (min_reps < 5) {
      warning(sprintf("Minimum replicates per group is %d. MANOVA may be underpowered.", min_reps))
    }
  }

  # Ensure the target term exists in the model
  if (!term %in% attr(stats::terms(formula), "term.labels") &&
      term != "groups" && term != "(Intercept)") {
    stop(sprintf("Term '%s' not found in model formula.", term))
  }

  num_peaks <- nrow(se)
  message(sprintf("Running shape shift test for %d peaks (model: %s) using %d threads...",
                  num_peaks, deparse(formula), workers))

  # ---- Parallel backend -----------------------------------------------------
  if (workers > 1) {
    param <- if (.Platform$OS.type == "unix") BiocParallel::MulticoreParam(workers) else BiocParallel::SnowParam(workers)
  } else {
    param <- BiocParallel::SerialParam()
  }

  # Extract matrices upfront to reduce per-worker I/O overhead
  w_mat <- assay(se, "Width")
  i_mat <- assay(se, "Intensity")
  h_mat <- assay(se, "Height")
  s_mat <- assay(se, "Skewness")

  # ---- Core parallel computation --------------------------------------------
  res_list <- BiocParallel::bplapply(seq_len(num_peaks), function(i) {

    Y <- cbind(
      W = w_mat[i, ],
      I = i_mat[i, ],
      H = h_mat[i, ],
      S = s_mat[i, ]
    )

    col_vars <- apply(Y, 2, stats::var, na.rm = TRUE)
    valid_cols <- col_vars > 1e-8

    if (sum(valid_cols) < 2) {
      return(c(Score = 0, Pval = 1, Eta2 = 0, Status = "Low_Variance"))
    }

    Y_valid <- Y[, valid_cols, drop = FALSE]
    k_valid <- sum(valid_cols)
    df_effect <- 1  # will be updated from fit

    tryCatch({
      fit <- stats::manova(Y_valid ~ ., data = model_df)
      res_pillai <- summary(fit, test = "Pillai")

      if (!term %in% rownames(res_pillai$stats)) {
        return(c(Score = 0, Pval = 1, Eta2 = 0, Status = "Term_Missing"))
      }

      pillai_val <- res_pillai$stats[term, "Pillai"]
      df_effect  <- res_pillai$stats[term, "Df"]
      s <- min(k_valid, df_effect)
      eta2 <- if (s > 0) pillai_val / s else 0

      if (test != "Pillai") {
        res <- summary(fit, test = test)
        score <- res$stats[term, test]
        pval  <- res$stats[term, "Pr(>F)"]
      } else {
        score <- pillai_val
        pval  <- res_pillai$stats[term, "Pr(>F)"]
      }

      return(c(Score = score, Pval = pval, Eta2 = eta2, Status = "OK"))
    }, error = function(e) {
      return(c(Score = 0, Pval = 1, Eta2 = 0, Status = "MANOVA_Error"))
    })

  }, BPPARAM = param)

  # Rapid re-assembly of parallel results
  res_mat <- do.call(rbind, res_list)
  shift_scores  <- as.numeric(res_mat[, "Score"])
  global_pvals  <- as.numeric(res_mat[, "Pval"])
  partial_eta2  <- as.numeric(res_mat[, "Eta2"])
  manova_status <- res_mat[, "Status"]

  # Extract genomic coordinates
  gr_df <- as.data.frame(rowRanges(se))
  # Tolerate different GRanges flattening conventions
  if(!"seqnames" %in% colnames(gr_df)) colnames(gr_df)[1] <- "seqnames"

  # Assemble final results table
  res_df <- data.frame(
    seqnames = gr_df$seqnames,
    start = gr_df$start,
    end = gr_df$end,
    Peak_ID = rownames(se),
    Shape_Shift_Score = shift_scores,
    partial_eta_sq = partial_eta2,
    P.Value = global_pvals,
    adj.P.Val = stats::p.adjust(global_pvals, method = "BH"),
    Status = manova_status,
    stringsAsFactors = FALSE
  )

  # Dual-sort: ascending adj.P.Val, descending partial_eta_sq
  res_df <- res_df[order(res_df$adj.P.Val, -res_df$partial_eta_sq), ]
  rownames(res_df) <- NULL

  n_errors <- sum(manova_status == "MANOVA_Error")
  n_lowvar <- sum(manova_status == "Low_Variance")
  if (n_errors > 0) message(sprintf("Note: %d peak(s) had MANOVA computation failures.", n_errors))
  if (n_lowvar > 0) message(sprintf("Note: %d peak(s) had insufficient variable dimensions (skipped).", n_lowvar))
  message("Success: Portrait shift MANOVA completed.")
  return(res_df)
}

#' Perform Differential Test on a Single Geometric Feature
#'
#' @description Performs a fast standard differential test on a single geometric dimension.
#' Supports "limma" (recommended for small samples), "t.test", or "wilcox". When a blocking
#' variable is provided, it is included as a covariate (limma) or used to compute paired
#' statistics (t.test/wilcox).
#'
#' @param se A \code{SummarizedExperiment} object from \code{build_portrait_matrix()}.
#' @param feature Character. The assay/feature to test (default: "Intensity").
#' @param group_var Character. The column in \code{colData(se)} containing the grouping variable.
#' @param target_group Character. The treatment group (numerator for FC).
#' @param ref_group Character. The reference/control group (denominator for FC).
#' @param method Character. The statistical test: "limma" (default), "t.test", or "wilcox".
#' @param block_var Character or NULL. Optional column in \code{colData(se)} specifying a
#'   blocking/pairing variable (e.g., patient ID). For limma, added as a covariate in the
#'   design matrix. For t.test and wilcox, enables paired testing.
#' @param workers Integer. Number of threads for parallel processing (default: 1).
#'
#' @return A \code{data.frame} containing logFC, P.Value, and adj.P.Val.
#'
#' @import SummarizedExperiment
#' @import BiocParallel
#' @importFrom stats t.test wilcox.test p.adjust model.matrix
#' @export
test_differential_feature <- function(se, feature = "Intensity", group_var = "Condition",
                                      target_group, ref_group, method = "limma",
                                      block_var = NULL, workers = 1) {

  if (!feature %in% assayNames(se)) {
    stop(sprintf("Feature '%s' not found in the SummarizedExperiment assays.", feature))
  }

  meta_data <- as.data.frame(colData(se))
  if (!group_var %in% colnames(meta_data)) stop(sprintf("Group variable '%s' not found.", group_var))

  # 1. Limma engine (efficient, no bplapply needed)
  if (method == "limma") {
    if (!requireNamespace("limma", quietly = TRUE)) stop("Please install 'limma'.")

    # Subset to relevant groups
    relevant_idx <- which(meta_data[[group_var]] %in% c(target_group, ref_group))
    se_sub <- se[, relevant_idx]
    mat <- assay(se_sub, feature)
    groups <- factor(colData(se_sub)[[group_var]], levels = c(ref_group, target_group))

    if (!is.null(block_var)) {
      if (!block_var %in% colnames(colData(se_sub))) {
        stop(sprintf("Block variable '%s' not found in sample metadata.", block_var))
      }
      blocks <- factor(colData(se_sub)[[block_var]])
      design <- stats::model.matrix(~ blocks + groups)
      group_coef <- ncol(design)  # 最后一个系数是组间效应
    } else {
      design <- stats::model.matrix(~ groups)
      group_coef <- 2
    }

    fit <- limma::lmFit(mat, design)
    fit <- limma::eBayes(fit)

    res_limma <- limma::topTable(fit, coef = group_coef, number = Inf, sort.by = "none")
    res_df <- data.frame(
      logFC = res_limma$logFC,
      P.Value = res_limma$P.Value,
      adj.P.Val = res_limma$adj.P.Val,
      row.names = rownames(se_sub),
      stringsAsFactors = FALSE
    )
    return(res_df)
  }

  # 2. Traditional t.test/Wilcoxon engine (with parallel support)
  idx_target <- which(meta_data[[group_var]] == target_group)
  idx_ref <- which(meta_data[[group_var]] == ref_group)
  mat <- assay(se, feature)

  # Paired design: align sample order by block_var, then enable paired testing
  if (!is.null(block_var)) {
    if (!block_var %in% colnames(meta_data)) {
      stop(sprintf("Block variable '%s' not found in sample metadata.", block_var))
    }
    blocks_target <- meta_data[[block_var]][idx_target]
    blocks_ref    <- meta_data[[block_var]][idx_ref]
    common_blocks <- intersect(blocks_target, blocks_ref)
    if (length(common_blocks) == 0) {
      stop("No shared block IDs found between target and reference groups. Cannot perform paired test.")
    }
    idx_target <- idx_target[match(common_blocks, blocks_target)]
    idx_ref    <- idx_ref[match(common_blocks, blocks_ref)]
    is_paired <- TRUE
  } else {
    is_paired <- FALSE
  }

  if (workers > 1) {
    param <- if (.Platform$OS.type == "unix") BiocParallel::MulticoreParam(workers) else BiocParallel::SnowParam(workers)
  } else {
    param <- BiocParallel::SerialParam()
  }

  res_list <- BiocParallel::bplapply(seq_len(nrow(mat)), function(i) {
    target_vals <- mat[i, idx_target]
    ref_vals <- mat[i, idx_ref]
    if (all(target_vals == 0) && all(ref_vals == 0)) return(c(logFC = 0, P.Value = 1))

    lfc <- log2((mean(target_vals, na.rm = TRUE) + 1) / (mean(ref_vals, na.rm = TRUE) + 1))
    pval <- 1
    tryCatch({
      if (method == "t.test") pval <- stats::t.test(target_vals, ref_vals, paired = is_paired)$p.value
      else pval <- stats::wilcox.test(target_vals, ref_vals, paired = is_paired)$p.value
    }, error = function(e) pval <- 1)

    return(c(logFC = lfc, P.Value = pval))
  }, BPPARAM = param)

  res_mat <- do.call(rbind, res_list)
  res_df <- data.frame(
    logFC = res_mat[, "logFC"],
    P.Value = res_mat[, "P.Value"],
    adj.P.Val = stats::p.adjust(res_mat[, "P.Value"], method = "BH"),
    row.names = rownames(se)
  )

  return(res_df)
}

