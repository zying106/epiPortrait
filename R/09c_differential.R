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
#' integer read counts. limma is used directly on the log-transformed continuous
#' intensity because it does not require the negative-binomial count modelling of
#' DESeq2/edgeR. This is a distinct, defensible design: the rigor is that the bit
#' is based on library-comparable normalized tracks, not raw fragment counts.
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
#' \code{contrast} matrix / \code{limma::makeContrasts}-style specification.
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
#'   matrix from \code{limma::makeContrasts}.
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
#' @return \code{se} with rowData columns:
#'   \code{<feature>_logFC}, \code{<feature>_AveExpr},
#'   \code{<feature>_t}, \code{<feature>_P.Value}, \code{<feature>_adj.P.Val},
#'   \code{<feature>_DiffStatus} (Gain / Loss / NS / NA);
#'   and \code{metadata(se)$differential_domains} with the limma fit summary,
#'   design, contrast and provenance.
#' @import SummarizedExperiment
#' @importFrom stats model.matrix
#' @examples
#' data(example_se)
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   se <- analyze_differential_domains(
#'     example_se, group_var = "Condition",
#'     ref_group = "Control", target_group = "Treatment")
#'   table(SummarizedExperiment::rowData(se)$Intensity_DiffStatus)
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
                                         robust = FALSE) {
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
  if (!feature %in% assayNames(se)) {
    stop("feature '", feature, "' not found in assays(se).", call. = FALSE)
  }
  if (ncol(se) < 2) {
    stop("Differential analysis requires at least 2 samples.", call. = FALSE)
  }
  meta <- as.data.frame(colData(se))

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
    se  <- se[, idx]
    meta <- as.data.frame(colData(se))
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
  # P1 (review §2/§3): NA must stay NA (missing != zero signal); negative
  # signal must NOT be silently clipped to 0 (would contradict
  # negative_policy = "allow" from build_portrait_matrix). EpigenDomain-level
  # differential analysis requires non-negative, non-missing intensity.
  M <- assay(se, feature)
  if (any(M < 0, na.rm = TRUE)) {
    stop("Differential analysis requires a non-negative assay, but '",
         feature, "' contains negative value(s). Re-run build_portrait_matrix() ",
         "with a non-negative (CPM/RPGC/spike-in) track, or analyze a ",
         "non-negative assay.", call. = FALSE)
  }
  m0 <- M
  Mt <- log2(m0 + 1)  # v1.0: only log2, so <feature>_logFC / logFC_cutoff /
                      # volcano x-axis are semantically consistent (review §11).
  # rows with NA signal are excluded from the fit (documented as untested);
  # rowMeans(na.rm=TRUE) would otherwise blur missing with zero.
  # drop rows below min_signal (based on mean transformed signal); these are
  # EXCLUDED from the limma fit (documented semantics) and reported as NA
  means <- rowMeans(Mt, na.rm = TRUE)
  fit_keep <- is.finite(means) & means >= min_signal
  Mt_fit <- Mt[fit_keep, , drop = FALSE]
  n_fit <- nrow(Mt_fit)

  # ---- limma fit on the filtered subset -------------------------------------
  # P0 fix: the fit runs on Mt_fit only (low-signal domains must not enter the
  # empirical-Bayes / trend estimation), and results are mapped back by the
  # ORIGINAL row index (fit_keep), never by position 1:n.
  # P1-robust (review §12): a 0-row fit is a user error (threshold too high /
  # all NA), not a silent all-NA result.
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
  # Remaining forms: single expression string, or numeric weights / matrix.
  fit2 <- if (is.character(contrast_used)) {
    limma::contrasts.fit(fit, limma::makeContrasts(contrasts = contrast_used,
                                                   levels = design_used))
  } else {
    limma::contrasts.fit(fit, contrast_used)
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
  # the ORIGINAL row indices (keep_idx), which handles non-contiguous fit_keep
  # correctly (regression: previously tt was indexed 1:n causing misalignment).
  out <- data.frame(matrix(NA_real_, nrow = nrow(se), ncol = 5))
  colnames(out) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
  keep_idx <- which(fit_keep)
  if (n_fit > 0) {
    out[keep_idx, ] <- tt[seq_len(n_fit), c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  }
  status <- rep(NA_character_, nrow(se))
  ok <- fit_keep & !is.na(out$logFC) & !is.na(out$adj.P.Val)
  # P1-1 (review): every tested domain is NS unless it meets the Gain/Loss
  # criteria — a large-effect-but-nonsignificant domain is "tested, NS", not
  # NA (NA must stay reserved for "not tested / excluded").
  if (any(ok)) status[ok] <- "NS"
  status[ok & out$logFC > logFC_cutoff & out$adj.P.Val < fdr_cutoff] <- "Gain"
  status[ok & out$logFC < -logFC_cutoff & out$adj.P.Val < fdr_cutoff] <- "Loss"

  pfx <- feature
  for (cc in colnames(out)) rowData(se)[[paste0(pfx, "_", cc)]] <- out[[cc]]
  rowData(se)[[paste0(pfx, "_DiffStatus")]] <- status

  # ---- provenance -----------------------------------------------------------
  S4Vectors::metadata(se)$differential_domains <- list(
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
    note = paste(
      "Differential test on continuous integrated BigWig signal via limma;",
      "not an exact read-count model. For DESeq2/DiffBind-style count",
      "statistics obtain per-domain counts externally.")
  )
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
#'   P-value (default constructed from \code{feature}).
#' @param label_n Integer. Number of top domains to label with
#'   \code{top_candidate_gene_symbol} / gene symbol if available (default 0).
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
                                      label_n = 0) {
  pfx <- feature
  if (is.null(logFC_col)) logFC_col <- paste0(pfx, "_logFC")
  if (is.null(padj_col)) padj_col <- paste0(pfx, "_adj.P.Val")
  rd <- as.data.frame(rowData(se), optional = TRUE)
  if (!logFC_col %in% colnames(rd) || !padj_col %in% colnames(rd)) {
    stop("Run analyze_differential_domains(feature = '", feature, "') first. ",
         "Missing '", logFC_col, "' / '", padj_col, "'.", call. = FALSE)
  }
  status_col <- paste0(pfx, "_DiffStatus")
  lbl_col <- if ("top_candidate_gene_symbol" %in% colnames(rd))
    "top_candidate_gene_symbol" else if ("gene_symbol" %in% colnames(rd))
    "nearest_tss_gene_symbol" else "Domain_ID"
  if (!lbl_col %in% colnames(rd)) lbl_col <- NULL

  df <- data.frame(
    logFC = rd[[logFC_col]],
    # P1 (review §4): untested domains (NA adj.P.Val) must NOT become 1e-300 ->
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
