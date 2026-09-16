#' Differential Domain Analysis on Continuous Signal (limma)
#'
#' @description Performs a group-wise differential analysis of a continuous
#' domain assay (default \code{Intensity}, the integrated BigWig signal within
#' each candidate domain) using limma's empirical-Bayes framework.
#' This complements the Super-domain state calls: where \code{compare_superdomains()}
#' reports relative \emph{state} transitions between discrete phenotypes,
#' \code{analyze_differential_domains()} reports a statistical \emph{P-value} for a
#' continuous signal change between condition groups.
#'
#' @details
#' \strong{Quantitative unit.} epiPortrait extracts \emph{continuous} signal
#' magnitude (e.g. \code{Intensity}) from normalized BigWig tracks rather than
#' integer read counts. limma is applied to log-transformed continuous intensity;
#' unlike DESeq2 or edgeR, this is not a negative-binomial count model. Valid
#' inference therefore requires quantitatively comparable normalized tracks and
#' an appropriate experimental design.
#'
#' \strong{Explicit caveat (recorded in provenance).} The test is on the
#' integrated BigWig signal, \emph{not} an exact read-count model. If the user
#' requires DESeq2/DiffBind-style count statistics, the recommended path is to
#' obtain per-domain counts externally and run those tools; epiPortrait provides
#' the coordinate universe and continuous phenotype layer.
#'
#' \strong{Design and contrasts.} \code{design} defines the linear model on the
#' samples (e.g. \code{~ Condition} or \code{~ 0 + Condition + batch}); contrast
#' is set either via the convenience \code{ref_group}/\code{target_group} pair
#' (default two-group comparison on a \code{Condition} column) or via a
#' one-column \code{contrast} matrix /
#' \code{limma::makeContrasts}-style specification. One call represents one
#' contrast; a matrix with multiple contrast columns is rejected explicitly.
#'
#' \strong{Transform offset caveat.} The variance-stabilizing transform is
#' \code{log(x + 1)}; the pseudo-count 1 is expressed in the units of the input
#' BigWigs, so its relative weight depends on the upstream normalization scheme
#' (CPM vs RPGC vs spike-in). With CPM-scale tracks the offset is negligible for
#' enriched domains but dominates near-background ones; set \code{min_signal}
#' accordingly to exclude domains where the fit would be driven by the offset.
#'
#' \strong{Mean-variance trend.} Integrated intensity shows a strong
#' mean-variance dependence across domains; by default a limma mean-variance
#' trend is fitted (\code{trend = TRUE}), which prevents low-signal domains from
#' being over-called as significant. The automatic fallback (fewer than 20
#' tested domains) and \code{robust = TRUE} outlier-resistant moderation follow
#' limma's recommendations for non-count genomic data.
#'
#' @param se A SummarizedExperiment from \code{build_portrait_matrix()}.
#' @param feature Character. Continuous assay to test (default "Intensity").
#' @param group_var Character. colData column holding the groups
#'   (default "Condition"). Only used when \code{design} is NULL.
#' @param ref_group,target_group Character. Reference and target group labels for
#'   the convenience two-group contrast (used when \code{design} is NULL).
#' @param design A model formula (e.g. \code{~ Condition}) defining the linear
#'   model, or NULL to build from \code{group_var} + \code{ref_group}/
#'   \code{target_group}. If supplied, \code{contrast} must also be supplied.
#' @param contrast A contrast specification: NULL (uses ref/target), a length-2
#'   character vector of column names to subtract (coef1 - coef2), or a numeric
#'   vector or a one-column numeric matrix from \code{limma::makeContrasts}.
#'   Multiple contrasts should be run in separate calls so each result has an
#'   unambiguous domain-level status and provenance record.
#' @param result_name Character or NULL. Optional unique label for a custom
#'   contrast. When NULL, a label is derived from the ref/target groups,
#'   contrast expression, or contrast-matrix column name. Pair-specific result
#'   columns and provenance use this label.
#' @param transform Character. "log2" (only; v1.0). Transform applied
#'   to the (non-negative) signal before limma.
#' @param min_signal Numeric. Filter: domains whose mean transformed signal is
#'   below this are dropped from the fit (NA results), to avoid fitting noise.
#'   Default 0 (no filtering).
#' @param logFC_cutoff Numeric. |log2-fold-change| threshold (on log2 scale) used
#'   only for the descriptive \code{DiffStatus} label (default 1).
#' @param fdr_cutoff Numeric. FDR threshold for the \code{DiffStatus} label
#'   (default 0.05).
#' @param trend Logical. Passed to \code{limma::eBayes}: fit a mean-variance
#'   trend across domains before moderation. Integrated intensity exhibits a
#'   strong mean-variance dependence, so the default is TRUE; it is disabled
#'   automatically (with a warning) when fewer than 20 domains pass filtering.
#' @param robust Logical. Passed to \code{limma::eBayes}: robust empirical-Bayes
#'   moderation, resistant to outlier domains (default FALSE). Consider TRUE
#'   when a small subset of extreme domains is expected.
#' @return \code{se} with pair-specific rowData columns:
#'   \code{<feature>_Diff__<result>__logFC} and the corresponding
#'   \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val}, and
#'   \code{DiffStatus} columns. Complete provenance records are appended to
#'   \code{metadata(se)$differential_domain_analyses} under the same
#'   \code{<feature>_Diff__<result>} key.
#' @import SummarizedExperiment
#' @importFrom stats model.matrix
#' @examples
#' data(example_se)
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   se <- analyze_differential_domains(
#'     example_se, group_var = "Condition",
#'     ref_group = "Control", target_group = "Treatment")
#'   table(SummarizedExperiment::rowData(se)[[
#'     "Intensity_Diff__Control_vs_Treatment__DiffStatus"]])
#' }
#' @export
analyze_differential_domains <- function(se,
                                         feature = "Intensity",
                                         group_var = "Condition",
                                         ref_group = NULL,
                                         target_group = NULL,
                                         design = NULL,
                                         contrast = NULL,
                                         transform = c("log2"),
                                         min_signal = 0,
                                         logFC_cutoff = 1,
                                         fdr_cutoff = 0.05,
                                         trend = TRUE,
                                         robust = FALSE,
                                         result_name = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("analyze_differential_domains() requires the 'limma' package. ",
         "Install it before using this function.", call. = FALSE)
  }
  transform <- match.arg(transform)
  if (!is.logical(trend) || length(trend) != 1L || is.na(trend)) {
    stop("trend must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
    stop("robust must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(result_name) &&
      (!is.character(result_name) || length(result_name) != 1L ||
       is.na(result_name) || !nzchar(result_name))) {
    stop("result_name must be NULL or one non-empty character string.",
         call. = FALSE)
  }
  if (!feature %in% assayNames(se)) {
    stop("feature '", feature, "' not found in assays(se).", call. = FALSE)
  }
  if (ncol(se) < 2) {
    stop("Differential analysis requires at least 2 samples.", call. = FALSE)
  }
  meta <- as.data.frame(colData(se))
  sample_idx <- seq_len(ncol(se))

  # ---- design / contrast resolution -----------------------------------------
  design_used <- design
  contrast_used <- contrast
  if (is.null(design)) {
    if (!group_var %in% colnames(meta)) {
      stop("group_var '", group_var, "' not found in colData.", call. = FALSE)
    }
    if (is.null(ref_group) || is.null(target_group)) {
      stop("design is NULL: provide ref_group and target_group (two-group ",
           "comparison on group_var).", call. = FALSE)
    }
    idx <- meta[[group_var]] %in% c(ref_group, target_group)
    sample_idx <- which(idx)
    meta <- meta[sample_idx, , drop = FALSE]
    # Drop unused factor levels after selecting the requested comparison.
    # Otherwise a legitimate factor with levels A/B/C, subset to A/B, retains
    # C=0 and is falsely rejected by the replicate-count check below.
    m <- factor(as.character(meta[[group_var]]),
                levels = c(ref_group, target_group))
    if (!ref_group %in% m || !target_group %in% m) {
      stop("ref_group / target_group not found in group_var.", call. = FALSE)
    }
    # limma requires >=2 replicates per group to estimate residual variance.
    ig <- table(m)
    if (any(ig < 2)) {
      stop("Each group needs >=2 replicates for limma residual-variance ",
           "estimation. Found: ", paste(sprintf("%s=%d", names(ig), ig),
                                        collapse = ", "), ". ",
           "Use a single-sample descriptive analysis or provide replicates.",
           call. = FALSE)
    }
    grp <- m
    design_used <- stats::model.matrix(~ grp)
    # coefficient for target_group: "grp<target>" (raw label, e.g. grp72h)
    tgt_col <- colnames(design_used)[
      vapply(colnames(design_used), function(nm) {
        stripped <- sub("^grp", "", nm)
        !is.na(stripped) & stripped == target_group
      }, logical(1))][1]
    if (is.na(tgt_col)) {
      stop("Could not resolve the target-group coefficient in the design ",
           "matrix.", call. = FALSE)
    }
    contrast_used <- tgt_col
  }
  # Materialize a user-supplied formula into the actual model matrix so that
  # coefficient names are resolvable for contrasts and provenance.
  if (!is.matrix(design_used)) {
    design_used <- stats::model.matrix(design_used, data = meta)
  }
  if (is.null(contrast_used)) {
    stop("design was supplied without contrast. Provide either ref_group/",
         "target_group (omitting design), a one-element expression string, ",
         "two design-coefficient names (coef1 - coef2), or a numeric ",
         "contrast matrix.", call. = FALSE)
  }

  # ---- transform + filter ---------------------------------------------------
  # NA must stay NA (missing is not zero signal); negative
  # signal must NOT be silently clipped to 0 (would contradict
  # negative_policy = "allow" from build_portrait_matrix). EpigenDomain-level
  # differential analysis requires non-negative, non-missing intensity.
  # Fit only the samples selected by the convenience pair, while retaining the
  # complete input object in the return value. This permits sequential
  # pairwise analyses on a study with more than two groups.
  M <- assay(se, feature)[, sample_idx, drop = FALSE]
  if (any(M < 0, na.rm = TRUE)) {
    stop("Differential analysis requires a non-negative assay, but '",
         feature, "' contains negative value(s). Re-run build_portrait_matrix() ",
         "with a non-negative (CPM/RPGC/spike-in) track, or analyze a ",
         "non-negative assay.", call. = FALSE)
  }
  m0 <- M
  # Use one log2 scale so the result logFC, logFC_cutoff and volcano x-axis
  # remain semantically consistent.
  Mt <- log2(m0 + 1)
  # rows with NA signal are excluded from the fit (documented as untested);
  # rowMeans(na.rm=TRUE) would otherwise blur missing with zero.
  # drop rows below min_signal (based on mean transformed signal); these are
  # EXCLUDED from the limma fit (documented semantics) and reported as NA
  means <- rowMeans(Mt, na.rm = TRUE)
  fit_keep <- is.finite(means) & means >= min_signal
  Mt_fit <- Mt[fit_keep, , drop = FALSE]
  n_fit <- nrow(Mt_fit)

  # ---- limma fit on the filtered subset -------------------------------------
  # The fit runs on Mt_fit only; low-signal domains must not enter the
  # empirical-Bayes or trend estimation, and results are mapped back by the
  # original row index (fit_keep), never by position 1:n. A zero-row fit means
  # that the threshold is too high or all values are missing.
  if (n_fit == 0L) {
    stop("No domains passed min_signal filtering (min_signal = ",
         min_signal, "). Lower min_signal or check that the assay is present ",
         "and non-NA.", call. = FALSE)
  }
  fit <- limma::lmFit(Mt_fit, design_used)
  # Normalize a two-coefficient name pair (coef1 - coef2) into numeric weights.
  if (is.character(contrast_used) && length(contrast_used) == 2L) {
    missing_coef <- setdiff(contrast_used, colnames(design_used))
    if (length(missing_coef) > 0L) {
      stop("Unknown design coefficient(s) in contrast: ",
           paste(missing_coef, collapse = ", "), ". Available: ",
           paste(colnames(design_used), collapse = ", "), ".", call. = FALSE)
    }
    w <- stats::setNames(numeric(ncol(design_used)), colnames(design_used))
    w[[contrast_used[[1]]]] <- 1
    w[[contrast_used[[2]]]] <- -1
    contrast_used <- w
  }
  if (is.matrix(contrast_used) && ncol(contrast_used) != 1L) {
    stop("Exactly one contrast is supported per call; received ",
         ncol(contrast_used), ". Run analyze_differential_domains() once per ",
         "contrast and optionally set result_name.", call. = FALSE)
  }
  # Remaining forms: single expression string, or numeric weights / matrix.
  fit2 <- if (is.character(contrast_used)) {
    limma::contrasts.fit(fit, limma::makeContrasts(contrasts = contrast_used,
                                                   levels = design_used))
  } else {
    limma::contrasts.fit(fit, contrast_used)
  }
  if (ncol(fit2$coefficients) != 1L) {
    stop("Exactly one contrast is supported per call; received ",
         ncol(fit2$coefficients), ". Run analyze_differential_domains() ",
         "once per contrast and optionally set result_name.", call. = FALSE)
  }
  # Mean-variance trend needs a reasonable number of tested domains.
  min_trend_rows <- 20L
  trend_used <- trend
  if (trend_used && n_fit < min_trend_rows) {
    warning(
      sprintf("Only %d domain(s) passed min_signal filtering; disabling the ",
              n_fit),
      sprintf("mean-variance trend (requires >= %d). Use trend = FALSE to silence.",
              min_trend_rows), call. = FALSE)
    trend_used <- FALSE
  }
  fit2 <- limma::eBayes(fit2, trend = trend_used, robust = robust)
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")

  # ---- map back to original rows by index -----------------------------------
  # topTable(sort.by = "none") is row-aligned to Mt_fit (the filtered subset),
  # so tt row k corresponds to the k-th KEPT domain. We populate the output by
  # the original row indices (keep_idx), which handles non-contiguous fit_keep.
  out <- data.frame(matrix(NA_real_, nrow = nrow(se), ncol = 5))
  colnames(out) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
  keep_idx <- which(fit_keep)
  if (n_fit > 0) {
    out[keep_idx, ] <- tt[seq_len(n_fit), c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  }
  status <- rep(NA_character_, nrow(se))
  ok <- fit_keep & !is.na(out$logFC) & !is.na(out$adj.P.Val)
  # Every tested domain is NS unless it meets the Gain/Loss
  # criteria — a large-effect-but-nonsignificant domain is "tested, NS", not
  # NA (NA must stay reserved for "not tested / excluded").
  if (any(ok)) status[ok] <- "NS"
  status[ok & out$logFC > logFC_cutoff & out$adj.P.Val < fdr_cutoff] <- "Gain"
  status[ok & out$logFC < -logFC_cutoff & out$adj.P.Val < fdr_cutoff] <- "Loss"

  # Build a stable, filesystem-friendly label for this one contrast. A unique
  # suffix is added when the same label is run again so no earlier result or
  # provenance record is lost.
  inferred_name <- if (!is.null(result_name)) {
    result_name
  } else if (is.null(design)) {
    paste0(ref_group, "_vs_", target_group)
  } else if (is.character(contrast)) {
    paste(contrast, collapse = "_minus_")
  } else if (is.matrix(contrast) && !is.null(colnames(contrast)) &&
             nzchar(colnames(contrast)[1])) {
    colnames(contrast)[1]
  } else {
    "contrast"
  }
  inferred_name <- gsub("[^[:alnum:]_.-]+", "_", inferred_name)
  base_key <- paste0(feature, "_Diff__", inferred_name)
  analyses <- S4Vectors::metadata(se)$differential_domain_analyses
  analysis_key <- base_key
  suffix <- 2L
  while (!is.null(analyses) && analysis_key %in% names(analyses)) {
    analysis_key <- paste0(base_key, "_", suffix)
    suffix <- suffix + 1L
  }

  canonical_columns <- paste0(analysis_key, "__",
                              c(colnames(out), "DiffStatus"))
  for (i in seq_along(out)) rowData(se)[[canonical_columns[i]]] <- out[[i]]
  rowData(se)[[canonical_columns[length(canonical_columns)]]] <- status

  # ---- provenance -----------------------------------------------------------
  prov_diff <- list(
    result_name = analysis_key,
    feature = feature,
    transform = transform,
    design = design_used,
    contrast = contrast_used,
    ref_group = if (is.null(design)) ref_group else NULL,
    target_group = if (is.null(design)) target_group else NULL,
    group_var = if (is.null(design)) group_var else NULL,
    min_signal = min_signal,
    logFC_cutoff = logFC_cutoff,
    fdr_cutoff = fdr_cutoff,
    trend = trend_used,
    robust = robust,
    n_tested = sum(fit_keep),
    n_gain = sum(status == "Gain", na.rm = TRUE),
    n_loss = sum(status == "Loss", na.rm = TRUE),
    n_ns = sum(status == "NS", na.rm = TRUE),
    created_columns = canonical_columns,
    note = paste(
      "Differential test on continuous integrated BigWig signal via limma;",
      "not an exact read-count model. For DESeq2/DiffBind-style count",
      "statistics obtain per-domain counts externally.")
  )
  if (is.null(analyses)) analyses <- list()
  analyses[[analysis_key]] <- prov_diff
  S4Vectors::metadata(se)$differential_domain_analyses <- analyses
  se
}


#' Differential Domain Volcano Plot
#'
#' @description Volcano plot of a differential-domain analysis performed by
#' \code{\link{analyze_differential_domains}}, using the stored \code{DiffStatus}
#' labels with the package's publication palette. Domains are coloured by
#' Gain / Loss / NS (grey).
#'
#' @param se A SummarizedExperiment after \code{analyze_differential_domains()}.
#' @param feature Character. Differential feature column (default "Intensity").
#' @param logFC_col,padj_col Character. Column names for logFC and adjusted
#'   P-value. By default they are resolved from the selected canonical
#'   differential result.
#' @param label_n Integer. Number of top domains to label with
#'   \code{top_candidate_gene_symbol} / gene symbol if available (default 0).
#' @param result_name Character or NULL. Differential result to plot. May be
#'   either the full metadata key (for example
#'   \code{"Intensity_Diff__Control_vs_Treatment"}) or its result suffix
#'   (\code{"Control_vs_Treatment"}). When NULL, the result is selected
#'   automatically only if exactly one analysis exists for \code{feature}; with
#'   multiple analyses, supply \code{result_name} or explicit column names.
#' @import ggplot2
#' @return A ggplot object.
#' @examples
#' data(example_se)
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   se <- analyze_differential_domains(
#'     example_se, group_var = "Condition",
#'     ref_group = "Control", target_group = "Treatment")
#'   plot_differential_volcano(se, label_n = 5)
#' }
#' @export
plot_differential_volcano <- function(se, feature = "Intensity",
                                      logFC_col = NULL, padj_col = NULL,
                                      label_n = 0, result_name = NULL) {
  pfx <- feature
  rd <- as.data.frame(rowData(se), optional = TRUE)
  analyses <- S4Vectors::metadata(se)$differential_domain_analyses
  if (is.null(analyses)) analyses <- list()
  analysis_keys <- names(analyses)[vapply(analyses, function(x) {
    identical(x$feature, feature)
  }, logical(1))]
  analysis_key <- NULL

  if (!is.null(result_name)) {
    if (!is.character(result_name) || length(result_name) != 1L ||
        is.na(result_name) || !nzchar(result_name)) {
      stop("result_name must be NULL or one non-empty character string.",
           call. = FALSE)
    }
    clean_name <- gsub("[^[:alnum:]_.-]+", "_", result_name)
    requested_keys <- unique(c(
      result_name,
      paste0(feature, "_Diff__", result_name),
      paste0(feature, "_Diff__", clean_name)
    ))
    hits <- intersect(requested_keys, analysis_keys)
    if (length(hits) != 1L) {
      stop("Differential result '", result_name, "' was not found for feature '",
           feature, "'. Available: ",
           if (length(analysis_keys) == 0L) "none" else
             paste(analysis_keys, collapse = ", "), ".", call. = FALSE)
    }
    analysis_key <- hits
  }

  # Explicit canonical columns can identify their own result key. This keeps
  # the column-level override useful without depending on a latest-result alias.
  if (is.null(analysis_key) && !is.null(logFC_col) &&
      grepl("__logFC$", logFC_col)) {
    candidate <- sub("__logFC$", "", logFC_col)
    if (candidate %in% analysis_keys) analysis_key <- candidate
  }
  if (is.null(analysis_key) && !is.null(padj_col) &&
      grepl("__adj\\.P\\.Val$", padj_col)) {
    candidate <- sub("__adj\\.P\\.Val$", "", padj_col)
    if (candidate %in% analysis_keys) analysis_key <- candidate
  }
  if (is.null(analysis_key) && length(analysis_keys) == 1L) {
    analysis_key <- analysis_keys
  }

  needs_default <- is.null(logFC_col) || is.null(padj_col)
  if (is.null(analysis_key) && needs_default && length(analysis_keys) > 1L) {
    stop("Multiple differential results exist for feature '", feature,
         "'. Supply result_name or explicit logFC_col and padj_col. Available: ",
         paste(analysis_keys, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(analysis_key)) {
    if (is.null(logFC_col)) logFC_col <- paste0(analysis_key, "__logFC")
    if (is.null(padj_col)) padj_col <- paste0(analysis_key, "__adj.P.Val")
  } else if (length(analysis_keys) == 0L) {
    # Read compatibility for objects produced before canonical result keys were
    # introduced. New analyses never create these unqualified columns.
    if (is.null(logFC_col) && paste0(pfx, "_logFC") %in% colnames(rd)) {
      logFC_col <- paste0(pfx, "_logFC")
    }
    if (is.null(padj_col) && paste0(pfx, "_adj.P.Val") %in% colnames(rd)) {
      padj_col <- paste0(pfx, "_adj.P.Val")
    }
  }

  if (is.null(logFC_col) || is.null(padj_col)) {
    stop("No canonical differential result found for feature '", feature,
         "'. Run analyze_differential_domains() first or supply explicit ",
         "logFC_col and padj_col.", call. = FALSE)
  }
  if (!logFC_col %in% colnames(rd) || !padj_col %in% colnames(rd)) {
    stop("Run analyze_differential_domains(feature = '", feature, "') first. ",
         "Missing '", logFC_col, "' / '", padj_col, "'.", call. = FALSE)
  }
  status_col <- if (!is.null(analysis_key)) {
    paste0(analysis_key, "__DiffStatus")
  } else if (grepl("__logFC$", logFC_col)) {
    sub("__logFC$", "__DiffStatus", logFC_col)
  } else {
    paste0(pfx, "_DiffStatus")
  }
  lbl_col <- if ("top_candidate_gene_symbol" %in% colnames(rd))
    "top_candidate_gene_symbol" else if ("gene_symbol" %in% colnames(rd))
    "nearest_tss_gene_symbol" else "Domain_ID"
  if (!lbl_col %in% colnames(rd)) lbl_col <- NULL

  df <- data.frame(
    logFC = rd[[logFC_col]],
    # Untested domains (NA adj.P.Val) must not become 1e-300 and therefore
    # -log10 = 300 "extreme significance"; keep them NA (dropped by ggplot).
    negLog10P = {
      prv <- rd[[padj_col]]
      nlp <- rep(NA_real_, length(prv))
      fin <- is.finite(prv) & prv >= 0
      nlp[fin] <- -log10(pmax(prv[fin], 1e-300))
      nlp
    },
    Status = if (status_col %in% colnames(rd)) rd[[status_col]] else NA_character_,
    stringsAsFactors = FALSE)
  df$Status <- factor(df$Status, levels = c("Gain", "Loss", "NS"))
  dom_ids <- rownames(se)
  if (is.null(dom_ids)) dom_ids <- paste0("Domain_", seq_len(nrow(se)))
  df$Label <- if (is.null(lbl_col)) {
    dom_ids
  } else {
    lab <- rd[[lbl_col]]
    if (is.null(lab)) dom_ids else lab
  }
  df$Label <- as.character(df$Label)
  df$Label[is.na(df$Label)] <- ""

  pal <- c("Gain" = "#D55E00", "Loss" = "#0072B2", "NS" = "#B8B8B8")
  p <- ggplot2::ggplot(df, ggplot2::aes(logFC, negLog10P, colour = Status)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.6) +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE,
                                 name = "Differential status") +
    .epi_theme_publication() +
    ggplot2::labs(title = paste0("Volcano plot: ", feature, " differential domains"),
                  x = paste0("log2 fold-change (", feature, ")"),
                  y = "-log10(adjusted P)") +
    .epi_wrap_legend("colour", n_items = 3)

  if (label_n > 0 && nrow(df) > 0 && any(nzchar(df$Label))) {
    top_idx <- order(-df$negLog10P)
    top_idx <- top_idx[nzchar(df$Label[top_idx])][seq_len(min(label_n,
                               sum(nzchar(df$Label))))]
    if (length(top_idx) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = df[top_idx, , drop = FALSE], ggplot2::aes(label = Label),
        size = 3, max.overlaps = 30, box.padding = 0.4, color = "black",
        segment.color = "grey50")
    }
  }
  p
}
