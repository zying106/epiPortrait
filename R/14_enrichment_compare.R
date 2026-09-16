# Cross-group comparison and visualization for domain-aware enrichment.
#
# Design: all groups share ONE universe (the linked-gene universe of the full
# object), one annotation rule and one FDR method. Cross-group comparison uses
# the explicit per-term effect size (log2 odds ratio with 95% CI), never
# -log10(FDR): groups with different foreground sizes have different power, and
# an FDR-only heatmap would reward the larger group.

#' Compare Functional Enrichment Across Domain Phenotypes or Transitions
#'
#' @description Runs the same ORA design for several domain selections
#'   (phenotypes, transitions, gain/loss states) and returns a long comparison
#'   table. Every group uses the same domain-derived testable universe, the same
#'   gene-link rules and the same multiple-testing method, so differences
#'   reflect biology rather than design. Per-term effect sizes (log2 odds ratio
#'   with 95\% CI) are reported alongside p-values.
#'
#' @param se A SummarizedExperiment after \code{annotate_epi_domains()}.
#' @param sets A named list of logical vectors (length nrow(se), no NAs), each
#'   selecting the domains of one group. Names become the group labels, e.g.
#'   \code{list(Intensity = ..., Breadth = ..., Dual = ...)} from
#'   \code{get_combined_class_results()}.
#' @param relations,mark,max_per_domain,nearest_tss_cutoff_bp,min_evidence_tier
#'   Gene-link selection, see \code{\link{get_domain_gene_universe}}.
#' @param org_db,key_type OrgDb / keytype; NULL resolves from annotation
#'   provenance.
#' @param database,ontology,min_gs_size,max_gs_size,p_adjust_method,padj_cutoff,simplify,simplify_cutoff,min_foreground
#'   Passed to \code{\link{enrich_epi_genes}}.
#' @param result_name Character or NULL. Metadata key; default
#'   \code{"comparison_ORA_GO_BP"} (made unique on collision).
#' @param BPPARAM A \code{BiocParallelParam} object for running the groups in
#'   parallel (default \code{BiocParallel::bpparam()}).
#' @return The SummarizedExperiment with an \code{epi_enrichment_comparison}
#'   object stored under \code{metadata(se)$enrichment_comparison[[result_name]]}.
#'   Retrieve it with \code{get_epi_enrichment(se, type = "comparison")}.
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
#'     requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   data(example_se)
#'   se <- annotate_epi_domains(example_se, genome = "hg38")
#'   se <- call_super_domains(se, feature = "Intensity", verbose = FALSE)
#'   cls <- SummarizedExperiment::rowData(se)$Intensity_Domain_Type
#'   se <- compare_epi_enrichment(
#'     se,
#'     sets = list(Super = cls == "Intensity_Super_Element",
#'                 Typical = cls == "Intensity_Typical"),
#'     min_foreground = 3)
#'   get_epi_enrichment(se, type = "comparison")
#' }
#' @seealso \code{\link{plot_epi_enrichment}}, \code{\link{enrich_epi_domains}}
#' @export
compare_epi_enrichment <- function(se, sets, relations = NULL, mark = NULL,
                                   max_per_domain = NULL,
                                   nearest_tss_cutoff_bp =
                                     getOption("epiPortrait.nearest_tss_cutoff_bp",
                                               10000),
                                   min_evidence_tier = 1L,
                                   org_db = NULL, key_type = NULL,
                                   database = "GO", ontology = "BP",
                                   min_gs_size = 10, max_gs_size = 500,
                                   p_adjust_method = "BH", padj_cutoff = 0.05,
                                   simplify = FALSE, simplify_cutoff = 0.7,
                                   min_foreground = 10,
                                   result_name = NULL,
                                   BPPARAM = BiocParallel::bpparam()) {
  if (!is.list(sets) || is.data.frame(sets) || length(sets) == 0) {
    stop("sets must be a non-empty named list of logical vectors.",
         call. = FALSE)
  }
  if (is.null(names(sets)) || any(!nzchar(names(sets))) ||
      anyDuplicated(names(sets))) {
    stop("sets must be named with unique, non-empty group labels.",
         call. = FALSE)
  }
  for (g in names(sets)) {
    v <- sets[[g]]
    if (!is.logical(v) || length(v) != nrow(se) || anyNA(v)) {
      stop(sprintf("sets[['%s']] must be a logical vector of length nrow(se) with no NAs.",
                   g), call. = FALSE)
    }
  }
  db <- .resolve_enrich_db(se, org_db = org_db, key_type = key_type)
  universe <- get_domain_gene_universe(
    se, relations = relations, mark = mark, domains = NULL,
    max_per_domain = max_per_domain,
    nearest_tss_cutoff_bp = nearest_tss_cutoff_bp,
    min_evidence_tier = min_evidence_tier,
    restrict_testable = TRUE, org_db = db$org_db, key_type = db$key_type)

  fg_list <- lapply(sets, function(sel) {
    links <- .enrich_links(se, domains = sel, relations = relations,
                           mark = mark, max_per_domain = max_per_domain,
                           nearest_tss_cutoff_bp = nearest_tss_cutoff_bp,
                           min_evidence_tier = min_evidence_tier,
                           unique_genes = FALSE)
    unique(links$gene_id)
  })
  results <- BiocParallel::bplapply(
    names(sets),
    function(g) {
      enrich_epi_genes(
        fg_list[[g]], universe, org_db = db$org_db, key_type = db$key_type,
        database = database, ontology = ontology,
        min_gs_size = min_gs_size, max_gs_size = max_gs_size,
        p_adjust_method = p_adjust_method, padj_cutoff = padj_cutoff,
        simplify = simplify, simplify_cutoff = simplify_cutoff,
        min_foreground = min_foreground)
    },
    BPPARAM = BPPARAM)
  names(results) <- names(sets)

  tables <- lapply(names(results), function(g) {
    r <- results[[g]]$result
    if (nrow(r) == 0) return(NULL)
    data.frame(group = g, n_fg = length(results[[g]]$foreground$genes),
               n_bg = results[[g]]$background$n, r,
               stringsAsFactors = FALSE)
  })
  tables <- tables[!vapply(tables, is.null, logical(1))]
  comp_table <- if (length(tables) == 0) {
    data.frame(group = character(), n_fg = integer(), n_bg = integer(),
               .empty_enrich_table(), stringsAsFactors = FALSE)
  } else {
    do.call(rbind, tables)
  }
  rownames(comp_table) <- NULL

  parameters <- list(
    database = database, ontology = ontology, org_db = db$org_db,
    key_type = db$key_type, min_gs_size = min_gs_size,
    max_gs_size = max_gs_size, p_adjust_method = p_adjust_method,
    padj_cutoff = padj_cutoff, simplify = simplify,
    simplify_cutoff = simplify_cutoff, min_foreground = min_foreground,
    relations = .default_enrich_relations(mark, relations), mark = mark,
    min_evidence_tier = min_evidence_tier, max_per_domain = max_per_domain,
    groups = vapply(sets, sum, integer(1)))
  comp <- structure(list(
    method = "ORA",
    table = comp_table,
    results = results,
    universe = universe,
    mapping = list(
      input_domains = nrow(se),
      universe_linked_genes = attr(universe, "n_linked"),
      universe_testable_genes = attr(universe, "n_testable")),
    parameters = parameters,
    provenance = list(
      call = match.call(),
      timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
      versions = .enrich_versions(db$org_db))),
    class = "epi_enrichment_comparison")

  base <- if (is.null(result_name)) {
    paste0("comparison_ORA_", database, "_", ontology)
  } else {
    result_name
  }
  existing <- S4Vectors::metadata(se)$enrichment_comparison
  name <- .unique_result_name(existing, base)
  comp$name <- name
  if (is.null(existing)) existing <- list()
  existing[[name]] <- comp
  S4Vectors::metadata(se)$enrichment_comparison <- existing
  se
}

#' Plot a Comparative Enrichment Heatmap
#'
#' @description Plots the stored \code{epi_enrichment_comparison} as a heatmap
#'   of signed per-term effect sizes (log2 odds ratio). Counts and significance
#'   are annotated in the tiles. For single-result dotplots, use
#'   \code{enrichplot::dotplot(as_enrich_result(res))} instead; this function
#'   deliberately implements only the comparative view that clusterProfiler
#'   does not provide.
#'
#' @param se A SummarizedExperiment with a stored comparison.
#' @param name Character or NULL. Comparison name; NULL uses the only stored
#'   comparison (error if several are stored).
#' @param terms Character or NULL. Term names to display. NULL selects the top
#'   \code{top_n} terms by maximum absolute effect, preferring terms significant
#'   in at least one group.
#' @param top_n Integer. Maximum number of terms when \code{terms = NULL}.
#' @param show_counts Logical. Annotate each tile with the term gene count
#'   (starred when \code{significant}).
#' @param fill_limits Numeric or NULL. Limits for the effect-size color scale.
#' @return A \code{ggplot} object.
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
#'     requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   data(example_se)
#'   se <- annotate_epi_domains(example_se, genome = "hg38")
#'   se <- call_super_domains(se, feature = "Intensity", verbose = FALSE)
#'   cls <- SummarizedExperiment::rowData(se)$Intensity_Domain_Type
#'   se <- compare_epi_enrichment(
#'     se,
#'     sets = list(Super = cls == "Intensity_Super_Element",
#'                 Typical = cls == "Intensity_Typical"),
#'     min_foreground = 3)
#'   plot_epi_enrichment(se)
#' }
#' @seealso \code{\link{compare_epi_enrichment}}
#' @export
plot_epi_enrichment <- function(se, name = NULL, terms = NULL, top_n = 20,
                                show_counts = TRUE, fill_limits = NULL) {
  comp <- get_epi_enrichment(se, name = name, type = "comparison")
  # get_epi_enrichment() with name = NULL always returns the named list (a
  # consistent contract). Unwrap the single stored comparison; require an
  # explicit `name` when several are stored so the wrong one is never plotted.
  if (!inherits(comp, "epi_enrichment_comparison")) {
    if (length(comp) != 1L) {
      stop("Several comparisons are stored; pass `name =` to select one. ",
           "Available: ", paste(names(comp), collapse = ", "), ".",
           call. = FALSE)
    }
    comp <- comp[[1]]
  }
  tab <- comp$table
  if (nrow(tab) == 0) {
    stop("The stored comparison has no terms to plot.", call. = FALSE)
  }
  if (is.null(terms)) {
    sig_terms <- unique(tab$term_name[tab$significant])
    pool <- if (length(sig_terms) > 0) {
      tab[tab$term_name %in% sig_terms, , drop = FALSE]
    } else {
      tab
    }
    ord <- order(-abs(pool$effect))
    terms <- unique(pool$term_name[ord])[seq_len(min(top_n, length(unique(pool$term_name))))]
  }
  tab <- tab[tab$term_name %in% terms, , drop = FALSE]
  tab$term_name <- factor(tab$term_name, levels = rev(unique(terms)))
  tab$group <- factor(tab$group, levels = names(comp$parameters$groups))
  tab$label <- if (isTRUE(show_counts)) {
    paste0(tab$count, ifelse(tab$significant, "*", ""))
  } else {
    ""
  }
  lim <- if (is.null(fill_limits)) {
    m <- max(abs(tab$effect), na.rm = TRUE)
    if (!is.finite(m) || m == 0) m <- 1
    c(-m, m)
  } else {
    fill_limits
  }
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = group, y = term_name,
                                         fill = effect)) +
    ggplot2::geom_tile(colour = "grey90") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                  high = "#B2182B", midpoint = 0,
                                  limits = lim, oob = scales::squish,
                                  name = "log2 OR") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 9))
  if (isTRUE(show_counts)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = label),
                                size = 3, colour = "grey20")
  }
  p
}

#' Convert an epi_enrichment to a clusterProfiler Result
#'
#' @description Returns the underlying clusterProfiler object stored at
#'   analysis time so users can reuse the whole enrichplot / clusterProfiler
#'   ecosystem (dotplot, cnetplot, emapplot, ...) without this package
#'   reimplementing plotting. Purely a converter: no recomputation.
#'
#' @param x An \code{epi_enrichment} object.
#' @param which Character. "raw" (default) or "simplified".
#' @return The stored \code{enrichResult} / \code{gseaResult} object.
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("GO.db", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   universe <- head(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,
#'                                        keytype = "ENTREZID"), 2000)
#'   res <- enrich_epi_genes(c("7157", "1956", "672", "675"), universe,
#'                           org_db = "org.Hs.eg.db", key_type = "ENTREZID",
#'                           simplify = FALSE, min_foreground = 3)
#'   str(as_enrich_result(res), max.level = 0)
#' }
#' @seealso \code{\link{enrich_epi_genes}}
#' @export
as_enrich_result <- function(x, which = c("raw", "simplified")) {
  which <- match.arg(which)
  if (!inherits(x, "epi_enrichment")) {
    stop("x must be an epi_enrichment object (see get_epi_enrichment()).",
         call. = FALSE)
  }
  obj <- if (which == "raw") x$object else x$simplified_object
  if (is.null(obj)) {
    stop(sprintf("No %s clusterProfiler object is stored ",
                 if (which == "raw") "raw" else "simplified"),
         if (which == "simplified") {
           "(simplify may have been disabled or failed)."
         } else {
           "(the result may be empty)."
         }, call. = FALSE)
  }
  obj
}

#' @export
print.epi_enrichment_comparison <- function(x, ...) {
  cat(sprintf("<epi_enrichment_comparison> %s | groups: %s\n",
              x$method, paste(names(x$parameters$groups), collapse = ", ")))
  cat(sprintf("  universe: %d testable genes | terms tested (union): %d\n",
              x$mapping$universe_testable_genes,
              length(unique(x$table$term_id))))
  for (g in names(x$results)) {
    r <- x$results[[g]]
    cat(sprintf("  [%s] %d foreground genes | %d terms | %d significant\n",
                g, length(r$foreground$genes), nrow(r$result),
                sum(r$result$significant, na.rm = TRUE)))
  }
  invisible(x)
}
