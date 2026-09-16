# Domain-aware functional enrichment (GSEA layer).
#
# Ranked enrichment is the primary mode for continuous remodeling axes: the
# per-domain statistic (limma moderated t, log2 width ratio, ...) is aggregated
# to a signed gene score by rank_epi_genes() and tested with clusterProfiler's
# GSEA (fgsea backend). The background is implicit in the ranked list, which is
# itself restricted to the domain-derived testable universe.

.tidy_gsea_result <- function(df, ontology, padj_cutoff) {
  if (nrow(df) == 0) return(.empty_enrich_table())
  out <- data.frame(
    term_id = as.character(df$ID),
    term_name = as.character(df$Description),
    ontology = ontology,
    count = as.integer(df$setSize),
    set_size = as.integer(df$setSize),
    gene_ratio = NA_real_,
    bg_ratio = NA_real_,
    effect = as.numeric(df$NES),
    effect_low = NA_real_,
    effect_high = NA_real_,
    p_value = as.numeric(df$pvalue),
    p_adjust = as.numeric(df$p.adjust),
    q_value = if ("qvalues" %in% colnames(df)) as.numeric(df$qvalues) else NA_real_,
    significant = as.numeric(df$p.adjust) <= padj_cutoff,
    gene_id = as.character(df$core_enrichment),
    stringsAsFactors = FALSE)
  out <- out[order(out$p_adjust, out$p_value, -abs(out$effect)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Ranked Functional Enrichment of Domain-Derived Gene Scores
#'
#' @description GSEA of a signed gene-level ranking produced by
#'   \code{\link{rank_epi_genes}}. Unlike ORA on discrete phenotype classes,
#'   this uses each domain's actual continuous value, so the result does not
#'   depend on where a "Super" cutoff was drawn. The ranked list is restricted
#'   to the domain-derived testable universe, and the gene score sign is
#'   inherited from \code{score_col} (positive = gain / expansion).
#'
#' @param ranked A data.frame with \code{gene_id} and \code{score} columns
#'   (from \code{rank_epi_genes()}), or a named numeric vector.
#' @param universe Character vector of testable background genes (from
#'   \code{get_domain_gene_universe()}); the ranked list is intersected with it.
#' @param org_db,key_type OrgDb package name and keytype; see
#'   \code{\link{get_domain_gene_universe}}.
#' @param database Character. Only "GO" is implemented.
#' @param ontology Character. One of "BP", "MF", "CC".
#' @param min_gs_size,max_gs_size Integer. Gene-set size filter; \code{max_gs_size}
#'   is capped at the ranked-list size.
#' @param p_adjust_method Character. Multiple-testing method (default "BH").
#' @param padj_cutoff Numeric. Adjusted-p threshold for the \code{significant}
#'   flag; all tested terms are returned.
#' @param seed Integer. RNG seed for the permutation test, applied with
#'   \code{withr::with_seed()} so the caller's RNG state is not modified. The
#'   seed is recorded in the provenance.
#' @param min_ranked Integer. Minimum ranked genes; below it a structured empty
#'   result is returned (never an error).
#' @param verbose Logical. Passed to \code{clusterProfiler::gseGO}.
#' @return An \code{epi_enrichment} object with \code{$result} holding NES as
#'   the effect size and the clusterProfiler object in \code{$object}.
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("GO.db", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   universe <- head(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,
#'                                        keytype = "ENTREZID"), 2000)
#'   ranked <- data.frame(gene_id = universe,
#'                        score = seq_along(universe) - length(universe) / 2)
#'   res <- gsea_epi_genes(ranked, universe = universe,
#'                         org_db = "org.Hs.eg.db", key_type = "ENTREZID")
#'   head(res$result)
#' }
#' @seealso \code{\link{rank_epi_genes}}, \code{\link{enrich_epi_genes}}
#' @export
gsea_epi_genes <- function(ranked, universe, org_db = NULL, key_type = NULL,
                           database = "GO", ontology = "BP",
                           min_gs_size = 10, max_gs_size = 500,
                           p_adjust_method = "BH", padj_cutoff = 0.05,
                           seed = 1, min_ranked = 100, verbose = FALSE) {
  if (!identical(database, "GO")) {
    stop(sprintf("database = \"%s\" is not implemented yet; only \"GO\" is ",
                 database),
         "currently supported (KEGG / Reactome are planned).", call. = FALSE)
  }
  if (!ontology %in% c("BP", "MF", "CC")) {
    stop("ontology must be one of \"BP\", \"MF\", \"CC\".")
  }
  if (is.data.frame(ranked)) {
    if (!all(c("gene_id", "score") %in% colnames(ranked))) {
      stop("ranked data.frame must contain 'gene_id' and 'score' columns ",
           "(as returned by rank_epi_genes()).", call. = FALSE)
    }
    scores <- ranked$score
    names(scores) <- ranked$gene_id
  } else if (is.numeric(ranked) && !is.null(names(ranked))) {
    scores <- ranked
  } else {
    stop("ranked must be a data.frame with gene_id/score, or a named ",
         "numeric vector.", call. = FALSE)
  }
  universe <- unique(as.character(universe))
  universe <- universe[!is.na(universe) & nzchar(universe)]
  parameters <- list(
    database = database, ontology = ontology, org_db = org_db,
    key_type = key_type, min_gs_size = min_gs_size,
    max_gs_size = max_gs_size, p_adjust_method = p_adjust_method,
    padj_cutoff = padj_cutoff, seed = seed, min_ranked = min_ranked)

  db <- .resolve_enrich_db(se = NULL, org_db = org_db, key_type = key_type)
  .enrich_require(unique(c("clusterProfiler", db$org_db, "GO.db")))
  keys <- AnnotationDbi::keys(.get_orgdb(db$org_db), keytype = db$key_type)
  uni_rest <- .restrict_to_testable(universe, keys)
  universe <- uni_rest$genes

  scores <- scores[!is.na(scores) & is.finite(scores)]
  n_ranked_raw <- length(scores)
  scores <- scores[names(scores) %in% universe]
  # Duplicated gene names would make the ranking ambiguous; keep the strongest
  # absolute score and report how many were collapsed.
  if (anyDuplicated(names(scores))) {
    ord <- order(-abs(scores))
    scores <- scores[ord]
    scores <- scores[!duplicated(names(scores))]
  }
  mapping <- list(
    n_ranked_raw = n_ranked_raw,
    n_ranked_in_universe = length(scores),
    n_background_testable = length(universe),
    terms_tested = 0L,
    terms_significant = 0L)
  if (length(scores) < min_ranked) {
    return(.empty_epi_enrichment(
      method = "GSEA", reason = "ranked_list_too_small",
      parameters = parameters, foreground_genes = names(scores),
      background_genes = universe, mapping = mapping))
  }
  if (length(unique(scores)) < 2) {
    return(.empty_epi_enrichment(
      method = "GSEA", reason = "ranked_scores_constant",
      parameters = parameters, foreground_genes = names(scores),
      background_genes = universe, mapping = mapping))
  }
  v <- sort(scores, decreasing = TRUE)
  # withr::with_seed keeps the caller's RNG state untouched (BiocCheck
  # requirement) while making the permutation test reproducible.
  cp <- withr::with_seed(seed, clusterProfiler::gseGO(
    geneList = v, OrgDb = .get_orgdb(db$org_db), keyType = db$key_type,
    ont = ontology, minGSSize = as.integer(min_gs_size),
    maxGSSize = min(as.integer(max_gs_size), length(v)),
    pvalueCutoff = 1, pAdjustMethod = p_adjust_method,
    # NOTE: do NOT pass `by = "fgsea"`. Newer clusterProfiler releases removed
    # that argument (GSEA now always uses the fgsea/enrichit backend and extra
    # arguments are forwarded to enrichit::gsea_gson, which rejects `by`).
    # Passing it broke gseGO() on the Bioconductor build.
    verbose = verbose, seed = FALSE))
  result <- .tidy_gsea_result(as.data.frame(cp), ontology, padj_cutoff)
  mapping$terms_tested <- nrow(result)
  mapping$terms_significant <- sum(result$significant, na.rm = TRUE)

  structure(list(
    method = "GSEA",
    result = result,
    simplified = NULL,
    object = cp,
    simplified_object = NULL,
    foreground = list(
      genes = names(v),
      ranked = data.frame(gene_id = names(v), score = as.numeric(v))),
    background = list(genes = universe, n = length(universe)),
    gene_links = data.frame(),
    mapping = mapping,
    parameters = parameters,
    provenance = list(
      call = match.call(),
      timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
      versions = .enrich_versions(db$org_db))),
    class = "epi_enrichment")
}
