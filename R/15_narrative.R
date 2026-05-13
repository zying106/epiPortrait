#' Generate Biological Narrative from Shape Shift Results
#'
#' @description Synthesises classified shape-shift results into a structured,
#' human-readable narrative suitable for a manuscript results section. Includes
#' top-level summary statistics, per-class biological interpretation, gene-level
#' examples, and suggested follow-up experiments. No AI-generated fluff—every
#' sentence is backed by numbers from your data.
#'
#' @param shift_res A data.frame from \code{classify_shape_shift()} containing
#'   classified shape-shift results.
#' @param se A \code{SummarizedExperiment} object (optional). If provided and
#'   contains a \code{SYMBOL} column in \code{rowData}, top genes per class
#'   are included in the narrative.
#' @param enrichment_res Optional enrichment results from
#'   \code{enrich_shape_shifted()}.
#' @param regulatory_res Optional regulatory overlap results from
#'   \code{annotate_regulatory()}.
#' @param fdr_cutoff Numeric. FDR threshold for significance (default: 0.05).
#'
#' @return A list of character strings:
#'   \item{summary}{One-paragraph top-level summary.}
#'   \item{class_breakdown}{Per-class narrative with gene examples.}
#'   \item{top_peaks}{Table-ready summary of the top 5 peaks.}
#'   \item{enrichment_narrative}{Functional enrichment interpretation (if available).}
#'   \item{regulatory_narrative}{Regulatory context interpretation (if available).}
#'   \item{follow_up}{Recommended next experiments based on the observed pattern.}
#'   \item{full_text}{All sections concatenated, ready for copy-paste.}
#'
#' @export
summarize_findings <- function(shift_res, se = NULL,
                                enrichment_res = NULL,
                                regulatory_res = NULL,
                                fdr_cutoff = 0.05) {

  required_cols <- c("Peak_ID", "Shape_Class", "adj.P.Val", "Shape_Shift_Score",
                     "partial_eta_sq", "Intensity_FC", "Height_FC", "Skewness_Delta")
  missing_cols <- setdiff(required_cols, colnames(shift_res))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in shift_res: %s. Run test_global_shape_shift() then classify_shape_shift() first.",
         paste(missing_cols, collapse = ", ")))
  }

  sig <- shift_res[shift_res$adj.P.Val <= fdr_cutoff, , drop = FALSE]
  n_total <- nrow(shift_res)
  n_sig   <- nrow(sig)
  pct_sig <- round(n_sig / n_total * 100, 1)

  # ---- Top-level summary ----------------------------------------------------
  class_counts <- table(sig$Shape_Class)
  class_counts <- class_counts[class_counts > 0]
  dominant_class <- names(class_counts)[which.max(class_counts)]

  class_narratives <- list(
    Concentration  = "concentrate at the summit, suggesting tighter regulatory focus",
    Flattening     = "spread laterally across the domain, consistent with regulatory broadening or loss of focal recruitment",
    Polarity_Shift = "undergo asymmetric redistribution within their boundaries, indicating directional spreading or collapse",
    `Global Gain`  = "show coordinated increases across multiple geometric dimensions, consistent with activation or opening",
    `Global Loss`  = "show coordinated decreases, consistent with repression or closing",
    Complex        = "exhibit discordant multi-dimensional changes, suggesting complex structural remodelling"
  )

  summary_lines <- c(
    sprintf("Of %d genomic regions profiled, %d (%.1f%%) showed significant",
            n_total, n_sig, pct_sig),
    sprintf("4D geometric shape shift (MANOVA-Pillai, FDR < %.2f).", fdr_cutoff),
    sprintf("The dominant pattern was %s (%d regions; %.1f%% of significant peaks),",
            dominant_class, class_counts[dominant_class],
            round(class_counts[dominant_class] / n_sig * 100, 1)),
    sprintf("where peaks %s.",
            class_narratives[[dominant_class]]),
    ""
  )

  # ---- Class breakdown -------------------------------------------------------
  class_lines <- c()
  for (cl in names(class_counts)) {
    cl_peaks <- sig[sig$Shape_Class == cl, ]
    cl_peaks <- cl_peaks[order(cl_peaks$adj.P.Val), ]
    n_cl <- nrow(cl_peaks)
    median_eta2 <- round(median(cl_peaks$partial_eta_sq, na.rm = TRUE), 3)

    class_lines <- c(class_lines,
      sprintf("**%s** (n = %d, median partial eta^2 = %.3f):",
              cl, n_cl, median_eta2))

    if (!is.null(se) && "SYMBOL" %in% colnames(rowData(se))) {
      gene_map <- rowData(se)$SYMBOL
      names(gene_map) <- rownames(se)
      cl_genes <- unique(na.omit(gene_map[cl_peaks$Peak_ID]))
      cl_genes <- cl_genes[cl_genes != ""]
      if (length(cl_genes) > 0) {
        top_genes <- head(cl_genes, 5)
        class_lines <- c(class_lines,
          sprintf("  Top associated genes: %s", paste(top_genes, collapse = ", ")))
        if (length(cl_genes) > 5) {
          class_lines <- c(class_lines,
            sprintf("  (...and %d more)", length(cl_genes) - 5))
        }
      }
    }

    if (cl %in% names(class_narratives)) {
      class_lines <- c(class_lines,
        sprintf("  Interpretation: peaks %s.", class_narratives[[cl]]))
    }
    class_lines <- c(class_lines, "")
  }

  # ---- Top peaks ------------------------------------------------------------
  top5 <- head(sig[order(sig$adj.P.Val), ], 5)
  top_peaks_lines <- c(
    "| Peak ID | Shape Class | adj.P.Val | partial eta^2 | Intensity FC | Height FC |",
    "|----------|-------------|-----------|----------------|--------------|-----------|"
  )
  for (i in seq_len(nrow(top5))) {
    top_peaks_lines <- c(top_peaks_lines,
      sprintf("| %s | %s | %.2e | %.3f | %.2f | %.2f |",
              top5$Peak_ID[i], top5$Shape_Class[i],
              top5$adj.P.Val[i], top5$partial_eta_sq[i],
              top5$Intensity_FC[i], top5$Height_FC[i]))
  }

  # ---- Enrichment narrative --------------------------------------------------
  enrich_lines <- NULL
  if (!is.null(enrichment_res)) {
    enrich_lines <- "Functional enrichment was performed per shape class. "
    if (isTRUE(enrichment_res$by_class)) {
      n_enriched <- sum(vapply(enrichment_res$go, function(x) {
        !is.null(x) && nrow(x) > 0
      }, logical(1)))
      enrich_lines <- c(enrich_lines,
        sprintf("GO terms were enriched in %d of %d shape classes tested.",
                n_enriched, length(enrichment_res$go)),
        "See enrichment results for class-specific pathway interpretation.",
        "")
    }
  }

  # ---- Regulatory narrative --------------------------------------------------
  reg_lines <- NULL
  if (!is.null(regulatory_res) && !is.null(regulatory_res$enrichment)) {
    enrich_mat <- regulatory_res$enrichment
    reg_lines <- c(
      sprintf("Shape-shifted regions were overlaid with %s.",
              regulatory_res$catalog_summary),
      ""
    )
    max_enrich <- which(enrich_mat == max(enrich_mat, na.rm = TRUE), arr.ind = TRUE)
    if (nrow(max_enrich) > 0) {
      reg_lines <- c(reg_lines,
        sprintf("The strongest regulatory enrichment was observed for **%s** peaks",
                rownames(enrich_mat)[max_enrich[1, 1]]),
        sprintf("in **%s** elements (%.1f-fold over genomic background).",
                colnames(enrich_mat)[max_enrich[1, 2]],
                enrich_mat[max_enrich[1, 1], max_enrich[1, 2]]),
        "")
    }
  }

  # ---- Follow-up experiments -------------------------------------------------
  follow_up_lines <- c("Recommended follow-up experiments based on these results:", "")

  if ("Concentration" %in% names(class_counts) && class_counts["Concentration"] > 3) {
    follow_up_lines <- c(follow_up_lines,
      "- **Concentration events**: Motif analysis (HOMER) on Concentration peaks to",
      "  identify TFs whose binding sites are enriched at sharpening regulatory elements.")
  }
  if ("Flattening" %in% names(class_counts) && class_counts["Flattening"] > 3) {
    follow_up_lines <- c(follow_up_lines,
      "- **Flattening events**: Compare with Hi-C or HiChIP data to assess whether",
      "  broadening corresponds to changes in 3D chromatin looping.")
  }
  if ("Polarity_Shift" %in% names(class_counts) && class_counts["Polarity_Shift"] > 3) {
    follow_up_lines <- c(follow_up_lines,
      "- **Polarity Shift events**: Examine the direction of skewness relative to",
      "  annotated TSS orientation—unidirectional shifts may suggest enhancer",
      "  scanning or transcriptional read-through.")
  }
  if (n_sig > 20) {
    follow_up_lines <- c(follow_up_lines,
      "- **Global validation**: Select 5-10 top peaks spanning different shape classes",
      "  for ChIP-qPCR validation of the predicted geometric changes.")
  }
  follow_up_lines <- c(follow_up_lines, "")

  # ---- Assemble full text ----------------------------------------------------
  full_text <- c(
    "## Shape Shift Analysis Summary",
    "",
    summary_lines,
    "### Per-Class Breakdown",
    "",
    class_lines,
    "### Top 5 Most Significant Shape-Shifted Regions",
    "",
    top_peaks_lines,
    ""
  )
  if (!is.null(enrich_lines)) {
    full_text <- c(full_text, "### Functional Enrichment", "", enrich_lines)
  }
  if (!is.null(reg_lines)) {
    full_text <- c(full_text, "### Regulatory Context", "", reg_lines)
  }
  full_text <- c(full_text, "### Suggested Follow-Up", "", follow_up_lines)

  list(
    summary              = summary_lines,
    class_breakdown      = class_lines,
    top_peaks            = top_peaks_lines,
    enrichment_narrative = enrich_lines,
    regulatory_narrative = reg_lines,
    follow_up            = follow_up_lines,
    full_text            = full_text
  )
}
