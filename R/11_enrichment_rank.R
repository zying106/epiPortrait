# Domain-aware functional interpretation: gene ranking layer.
#
# The enrichment module reuses the existing annotation evidence
# (metadata(se)$domain_gene_links via get_domain_genes()) instead of creating
# a second domain-to-gene vocabulary. This file provides the pure helpers and
# the continuous-score -> gene-score ranking used by ORA foreground selection
# and by GSEA.

# Mark-aware default evidence relations. Broad repressive marks are restricted
# to promoter_overlap on purpose: their domains span large gene deserts, so
# gene-body / nearest-TSS links are dominated by domain width rather than by
# mark biology. Active marks use promoter overlap + proximal nearest TSS +
# BEDPE 3D contacts (far nearest-TSS links stay tier 0 and are dropped by the
# default min_evidence_tier).
.default_enrich_relations <- function(mark = NULL, relations = NULL) {
  valid <- c("nearest_tss", "promoter_overlap", "gene_body_overlap",
             "fully_contained", "bedpe_promoter_contact")
  if (!is.null(relations)) {
    if (!is.character(relations) || length(relations) == 0 ||
        anyNA(relations) || !all(relations %in% valid)) {
      stop("relations must be a character vector drawn from: ",
           paste(valid, collapse = ", "), ".", call. = FALSE)
    }
    return(unique(relations))
  }
  mark_class <- if (is.null(mark)) "generic" else get_mark_preset(mark)$mark_class
  if (identical(mark_class, "broad_repressive")) {
    "promoter_overlap"
  } else {
    c("promoter_overlap", "nearest_tss", "bedpe_promoter_contact")
  }
}

# Shared link selection for the enrichment module. Always returns long-format
# domain-gene pairs (unique_genes = FALSE) so aggregation / foreground counts
# can report how many domains support each gene.
.enrich_links <- function(se, domains = NULL, relations = NULL, mark = NULL,
                          max_per_domain = NULL,
                          nearest_tss_cutoff_bp =
                            getOption("epiPortrait.nearest_tss_cutoff_bp",
                                      10000),
                          min_evidence_tier = 1L,
                          unique_genes = FALSE) {
  if (length(min_evidence_tier) != 1L || !is.numeric(min_evidence_tier) ||
      !is.finite(min_evidence_tier) ||
      min_evidence_tier < 0 || min_evidence_tier > 4) {
    stop("min_evidence_tier must be a single number in [0, 4].", call. = FALSE)
  }
  links <- get_domain_genes(
    se,
    domains = domains,
    relations = .default_enrich_relations(mark, relations),
    unique_genes = unique_genes,
    nearest_tss_cutoff_bp = nearest_tss_cutoff_bp,
    max_per_domain = max_per_domain)
  if (nrow(links) == 0) return(links)
  links <- links[links$evidence_tier >= min_evidence_tier, , drop = FALSE]
  rownames(links) <- NULL
  links
}

# Intersect a gene vector with the annotation-database key space (the
# "testable" genes for an enrichment test) and report the loss.
.restrict_to_testable <- function(genes, keys) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes)]
  kept <- unique(genes[genes %in% keys])
  list(genes = kept,
       n_input = length(genes),
       n_kept = length(kept),
       n_dropped = length(genes) - length(kept))
}

# 2x2 table effect size for ORA terms.
#   a: in term & in foreground     b: in foreground only
#   c: in term & in background     d: in neither
# Returns log2 odds ratio with a 95% CI. Haldane-Anscombe 0.5 correction is
# applied per contrast when any cell is zero (standard practice; the
# uncorrected OR is undefined or infinite).
.ora_effect_size <- function(a, b, c, d) {
  corr <- ifelse(a == 0 | b == 0 | c == 0 | d == 0, 0.5, 0)
  a2 <- a + corr; b2 <- b + corr; c2 <- c + corr; d2 <- d + corr
  log2_or <- log2((a2 * d2) / (b2 * c2))
  # The standard error is derived on the natural-log odds-ratio scale.
  # Convert it to log2 units before combining it with log2_or.
  se_log2 <- sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2) / log(2)
  data.frame(effect = log2_or,
             effect_low = log2_or - 1.959964 * se_log2,
             effect_high = log2_or + 1.959964 * se_log2)
}

# Collapse multiple domain links per gene into one signed score.
# links: data.frame with gene_id / gene_symbol (from get_domain_genes())
# scores / weights: numeric vectors parallel to links rows
.aggregate_domain_scores <- function(links, scores, weights, aggregate) {
  sp <- split(seq_along(scores), links$gene_id)
  score <- vapply(sp, function(i) {
    s <- scores[i]; w <- weights[i]
    switch(aggregate,
           signed_weighted = sum(w * s) / sum(w),
           mean = mean(s),
           sum = sum(s),
           max_abs = s[which.max(abs(s))],
           best_domain = {
             o <- order(w, abs(s), decreasing = TRUE)
             s[o[1]]
           })
  }, numeric(1))
  weight_sum <- vapply(sp, function(i) sum(weights[i]), numeric(1))
  symbol <- vapply(sp, function(i) {
    x <- links$gene_symbol[i]
    x <- x[!is.na(x)]
    if (length(x) > 0) x[1] else NA_character_
  }, character(1))
  n_domains <- vapply(sp, length, integer(1))
  out <- data.frame(gene_id = names(sp),
                    gene_symbol = unname(symbol),
                    score = unname(score),
                    n_domains = unname(n_domains),
                    weight_sum = unname(weight_sum),
                    stringsAsFactors = FALSE)
  out <- out[order(-out$score, out$gene_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Rank Genes by a Continuous Domain-Level Score
#'
#' @description Converts a continuous per-domain quantity (differential
#'   statistic, width transition, or any numeric \code{rowData} column) into a
#'   signed gene-level score suitable for GSEA. This is the continuous
#'   counterpart of the discrete-phenotype ORA path: instead of comparing
#'   "Super" boxes, every domain contributes its actual value.
#'
#' @details Domain-gene links are read from the evidence stored by
#'   \code{annotate_epi_domains()} via \code{get_domain_genes()}. Each gene may
#'   be linked to several domains; \code{aggregate} defines how those values
#'   are collapsed:
#'   \itemize{
#'     \item \code{"signed_weighted"} (default): evidence-weighted mean of the
#'           signed scores (\code{sum(w * s) / sum(w)}). Weights default to
#'           \code{2^evidence_tier} (tier 0-4), so promoter-overlapping
#'           evidence dominates a distal nearest-TSS link.
#'     \item \code{"best_domain"}: score of the strongest-evidence link only
#'           (most conservative; avoids one gene accumulating many domains).
#'     \item \code{"mean"}: unweighted mean.
#'     \item \code{"sum"}: unweighted sum (genes in many domains score higher;
#'           kept for completeness, not recommended).
#'     \item \code{"max_abs"}: signed score with the largest absolute value
#'           (winner-takes-all; sensitive to sign instability, use as a
#'           sensitivity analysis).
#'   }
#'   Positive scores mean gain / expansion / higher signal in the direction of
#'   the supplied statistic. The sign convention is therefore inherited from
#'   \code{score_col}.
#'
#' @param se A SummarizedExperiment after \code{annotate_epi_domains()}.
#' @param score_col Character. A numeric \code{rowData(se)} column, for example
#'   \code{"Intensity_Diff__Control_vs_Treatment__t"} or
#'   \code{"Intensity_Diff__Control_vs_Treatment__logFC"} from
#'   \code{analyze_differential_domains()}, or
#'   \code{"log2WidthRatio__A_vs_B"} from \code{compute_width_transition()}.
#' @param relations Character or NULL. Evidence relation types passed to
#'   \code{get_domain_genes()}; NULL uses mark-aware defaults (see
#'   \code{enrich_epi_domains()}).
#' @param mark Character or NULL. Mark name used for mark-aware defaults
#'   (e.g. "H3K27ac", "H3K4me3", "H3K27me3").
#' @param max_per_domain Integer or NULL. Optional per-domain gene cap applied
#'   before aggregation (see \code{get_domain_genes()}).
#' @param nearest_tss_cutoff_bp Numeric. Proximal-TSS cutoff passed to
#'   \code{get_domain_genes()} (default 10000).
#' @param min_evidence_tier Integer in \code{0:4}. Links below this tier are
#'   dropped before aggregation (default 1, i.e. far nearest-TSS links alone do
#'   not qualify a gene).
#' @param aggregate Character. Aggregation rule (see Details).
#' @param weights Numeric vector of length 5 or NULL. Evidence-tier weights for
#'   tier 0-4; NULL uses \code{2^(0:4)}.
#' @return A data.frame with \code{gene_id}, \code{gene_symbol},
#'   \code{score}, \code{n_domains} (linked domains with a finite score) and
#'   \code{weight_sum}, ordered by descending score. Attributes record the
#'   parameter set and link-mapping statistics.
#' @examples
#' data(example_se)
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   se <- annotate_epi_domains(example_se, genome = "hg38")
#'   se <- analyze_differential_domains(se, feature = "Intensity",
#'                                      ref_group = "Control",
#'                                      target_group = "Treatment")
#'   ranked <- rank_epi_genes(
#'     se, score_col = "Intensity_Diff__Control_vs_Treatment__t")
#'   head(ranked)
#' }
#' @seealso \code{\link{enrich_epi_domains}}, \code{\link{gsea_epi_genes}}
#' @export
rank_epi_genes <- function(se, score_col, relations = NULL, mark = NULL,
                           max_per_domain = NULL,
                           nearest_tss_cutoff_bp =
                             getOption("epiPortrait.nearest_tss_cutoff_bp",
                                       10000),
                           min_evidence_tier = 1L,
                           aggregate = c("signed_weighted", "best_domain",
                                         "mean", "max_abs", "sum"),
                           weights = NULL) {
  aggregate <- match.arg(aggregate)
  if (!is.character(score_col) || length(score_col) != 1L ||
      is.na(score_col) || !nzchar(score_col)) {
    stop("score_col must be a single rowData column name.", call. = FALSE)
  }
  rd <- as.data.frame(rowData(se), optional = TRUE)
  if (!score_col %in% colnames(rd)) {
    stop(sprintf("score_col '%s' not found in rowData(se). ", score_col),
         "Available numeric columns: ",
         paste(colnames(rd)[vapply(rd, is.numeric, logical(1))], collapse = ", "),
         call. = FALSE)
  }
  scores_all <- rd[[score_col]]
  if (!is.numeric(scores_all)) {
    stop(sprintf("rowData column '%s' is not numeric.", score_col), call. = FALSE)
  }
  names(scores_all) <- rownames(se)

  links <- .enrich_links(se, domains = NULL, relations = relations,
                         mark = mark, max_per_domain = max_per_domain,
                         nearest_tss_cutoff_bp = nearest_tss_cutoff_bp,
                         min_evidence_tier = min_evidence_tier,
                         unique_genes = FALSE)
  if (nrow(links) == 0) {
    stop("No domain-gene links pass the evidence filter. Run ",
         "annotate_epi_domains() first, or lower min_evidence_tier.",
         call. = FALSE)
  }
  s <- scores_all[match(links$domain_id, names(scores_all))]
  keep <- !is.na(s) & is.finite(s)
  n_links_raw <- nrow(links)
  n_genes_raw <- length(unique(links$gene_id))
  n_domains_scored <- length(unique(links$domain_id[keep]))
  links <- links[keep, , drop = FALSE]
  s <- unname(s[keep])
  if (nrow(links) == 0) {
    stop(sprintf("No linked domain has a finite value in '%s'.", score_col),
         call. = FALSE)
  }
  tier <- as.integer(links$evidence_tier)
  if (is.null(weights)) {
    w <- 2^tier
  } else {
    if (!is.numeric(weights) || length(weights) != 5L || anyNA(weights) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be a numeric vector of length 5 (tiers 0-4).",
           call. = FALSE)
    }
    w <- weights[tier + 1L]
  }

  out <- .aggregate_domain_scores(links, s, w, aggregate)
  attr(out, "score_col") <- score_col
  attr(out, "aggregate") <- aggregate
  attr(out, "relations") <- .default_enrich_relations(mark, relations)
  attr(out, "mark") <- mark
  attr(out, "min_evidence_tier") <- min_evidence_tier
  attr(out, "max_per_domain") <- max_per_domain
  attr(out, "weights") <- if (is.null(weights)) 2^(0:4) else weights
  attr(out, "mapping") <- list(
    n_links = nrow(links),
    n_links_raw = n_links_raw,
    n_domains_scored = n_domains_scored,
    n_genes = nrow(out),
    n_genes_raw = n_genes_raw,
    dropped_infinite_links = n_links_raw - nrow(links))
  out
}
