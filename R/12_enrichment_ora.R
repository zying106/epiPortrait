# Domain-aware functional enrichment (ORA layer).
#
# Design: reuse the existing domain-gene evidence instead of introducing a
# second mapping vocabulary. The foreground is selected by the caller with the
# existing getters (get_combined_class_results(), get_transition_results(),
# rowData columns), the background is the linked-gene universe of the same
# object (never the whole genome), and the per-term 2x2 effect size is computed
# explicitly so results can be compared across phenotypes without relying on
# -log10(FDR).

.enrich_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf(
      "Package(s) required for functional enrichment are not installed: %s. ",
      paste(missing, collapse = ", ")),
      "Install with BiocManager::install(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), ")).",
      call. = FALSE)
  }
}

.get_orgdb <- function(org_db) {
  get(org_db, envir = asNamespace(org_db))
}

# Resolve the OrgDb / keytype from the object's annotation provenance when not
# supplied. The annotation layer already records which OrgDb and keytype were
# used, so enrichment cannot silently drift away from it.
.resolve_enrich_db <- function(se = NULL, org_db = NULL, key_type = NULL) {
  prov <- NULL
  if (!is.null(se)) prov <- S4Vectors::metadata(se)$annotation_provenance
  if (is.null(org_db) && !is.null(prov)) org_db <- prov$annotation_db
  if (is.null(org_db)) {
    stop("No OrgDb available. Run annotate_epi_domains() with a built-in ",
         "genome, or supply `org_db` (e.g. \"org.Hs.eg.db\").",
         call. = FALSE)
  }
  if (!is.character(org_db) || length(org_db) != 1L || !nzchar(org_db)) {
    stop("org_db must be a single package name (e.g. \"org.Hs.eg.db\").",
         call. = FALSE)
  }
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop(sprintf("OrgDb package '%s' is required but not installed.", org_db),
         call. = FALSE)
  }
  if (is.null(key_type) && !is.null(prov)) key_type <- prov$gene_id_keytype
  if (is.null(key_type) &&
      org_db %in% c("org.Hs.eg.db", "org.Mm.eg.db")) {
    key_type <- "ENTREZID"
  }
  if (is.null(key_type)) {
    stop("key_type could not be inferred from the object; supply it ",
         "explicitly (e.g. \"ENTREZID\").", call. = FALSE)
  }
  list(org_db = org_db, key_type = key_type)
}

.enrich_versions <- function(org_db) {
  list(
    clusterProfiler = tryCatch(
      as.character(utils::packageVersion("clusterProfiler")),
      error = function(e) NA_character_),
    GO_db = tryCatch(
      as.character(utils::packageVersion("GO.db")),
      error = function(e) NA_character_),
    org_db = tryCatch(
      as.character(utils::packageVersion(org_db)),
      error = function(e) NA_character_))
}

.empty_enrich_table <- function() {
  data.frame(term_id = character(), term_name = character(),
             ontology = character(), count = integer(), set_size = integer(),
             gene_ratio = numeric(), bg_ratio = numeric(),
             effect = numeric(), effect_low = numeric(), effect_high = numeric(),
             p_value = numeric(), p_adjust = numeric(), q_value = numeric(),
             significant = logical(), gene_id = character(),
             stringsAsFactors = FALSE)
}

.empty_epi_enrichment <- function(method, reason, parameters = list(),
                                  foreground_genes = character(),
                                  background_genes = character(),
                                  mapping = list()) {
  mapping$reason <- reason
  structure(list(
    method = method,
    result = .empty_enrich_table(),
    simplified = NULL,
    object = NULL,
    simplified_object = NULL,
    foreground = list(genes = unique(as.character(foreground_genes)),
                      domains = character()),
    background = list(genes = unique(as.character(background_genes)),
                      n = length(unique(as.character(background_genes)))),
    gene_links = data.frame(),
    mapping = mapping,
    parameters = parameters,
    provenance = list(
      timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
      versions = list())),
    class = "epi_enrichment")
}

.parse_ratio_num <- function(x) {
  as.numeric(sub("/.*$", "", x))
}

.parse_ratio_den <- function(x) {
  as.numeric(sub("^[^/]*/", "", x))
}

# Convert a clusterProfiler ORA result into the tidy epiPortrait schema with an
# explicit 2x2 effect size. The background is restricted to the analysis
# universe, so b/c/d are always defined.
.tidy_ora_result <- function(df, ontology, n_fg, n_universe, padj_cutoff) {
  if (nrow(df) == 0) return(.empty_enrich_table())
  a <- as.numeric(df$Count)
  set_size <- .parse_ratio_num(df$BgRatio)
  # GO testing may exclude genes lacking annotation in the selected ontology.
  # Use the tested denominators reported by the enrichment engine so effect
  # sizes and ratios describe the same contingency table as the P-value.
  n_fg <- .parse_ratio_den(df$GeneRatio)
  n_universe <- .parse_ratio_den(df$BgRatio)
  n_bg <- n_universe - n_fg
  b <- pmax(n_fg - a, 0)
  cc <- pmax(set_size - a, 0)
  d <- pmax(n_bg - cc, 0)
  eff <- .ora_effect_size(a, b, cc, d)
  out <- data.frame(
    term_id = as.character(df$ID),
    term_name = as.character(df$Description),
    ontology = ontology,
    count = as.integer(a),
    set_size = as.integer(set_size),
    gene_ratio = .parse_ratio_num(df$GeneRatio) / n_fg,
    bg_ratio = set_size / n_universe,
    effect = eff$effect,
    effect_low = eff$effect_low,
    effect_high = eff$effect_high,
    p_value = as.numeric(df$pvalue),
    p_adjust = as.numeric(df$p.adjust),
    q_value = if ("qvalue" %in% colnames(df)) as.numeric(df$qvalue) else NA_real_,
    significant = as.numeric(df$p.adjust) <= padj_cutoff,
    gene_id = as.character(df$geneID),
    stringsAsFactors = FALSE)
  out <- out[order(out$p_adjust, out$p_value, -out$effect), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Build the Testable Gene Universe of an epiPortrait Object
#'
#' @description Returns the genes linked to the object's domains (the genes
#'   that "had the opportunity to be selected" under the same analysis design),
#'   optionally restricted to genes present in the annotation database. This is
#'   the background for \code{enrich_epi_genes()} / \code{enrich_epi_domains()};
#'   using all genes of an organism instead would ignore the domain universe and
#'   inflate enrichment.
#'
#' @param se A SummarizedExperiment after \code{annotate_epi_domains()}.
#' @param relations Character or NULL. Evidence relation types (see
#'   \code{\link{get_domain_genes}}); NULL uses mark-aware defaults.
#' @param mark Character or NULL. Mark name for mark-aware defaults
#'   (e.g. "H3K27ac", "H3K4me3", "H3K27me3").
#' @param domains Logical vector (length nrow(se)) or NULL. Restrict the
#'   universe to a subset of domains (default: all domains).
#' @param max_per_domain Integer or NULL. Optional per-domain gene cap.
#' @param nearest_tss_cutoff_bp Numeric. Proximal-TSS cutoff for
#'   \code{nearest_tss} links (default 10000).
#' @param min_evidence_tier Integer in \code{0:4}. Drop links below this tier
#'   (default 1: a far nearest-TSS link alone does not qualify a gene).
#' @param restrict_testable Logical. Intersect with the annotation database key
#'   space (default TRUE), so only testable genes enter the ORA background.
#' @param org_db Character or NULL. OrgDb package name; NULL is resolved from
#'   the object's annotation provenance.
#' @param key_type Character or NULL. Gene ID keytype; NULL is resolved from
#'   the annotation provenance (ENTREZID for the built-in genomes).
#' @return A character vector of gene IDs with class
#'   \code{epi_gene_universe} and attributes recording the link and restriction
#'   statistics.
#' @examples
#' data(example_se)
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   se <- annotate_epi_domains(example_se, genome = "hg38")
#'   u <- get_domain_gene_universe(se)
#'   length(u)
#' }
#' @seealso \code{\link{enrich_epi_domains}}
#' @export
get_domain_gene_universe <- function(se, relations = NULL, mark = NULL,
                                     domains = NULL, max_per_domain = NULL,
                                     nearest_tss_cutoff_bp =
                                       getOption("epiPortrait.nearest_tss_cutoff_bp",
                                                 10000),
                                     min_evidence_tier = 1L,
                                     restrict_testable = TRUE,
                                     org_db = NULL, key_type = NULL) {
  links <- .enrich_links(se, domains = domains, relations = relations,
                         mark = mark, max_per_domain = max_per_domain,
                         nearest_tss_cutoff_bp = nearest_tss_cutoff_bp,
                         min_evidence_tier = min_evidence_tier,
                         unique_genes = FALSE)
  genes <- unique(links$gene_id)
  n_linked <- length(genes)
  n_testable <- NA_integer_
  n_dropped <- NA_integer_
  db <- list(org_db = NA_character_, key_type = NA_character_)
  if (restrict_testable) {
    db <- .resolve_enrich_db(se, org_db = org_db, key_type = key_type)
    keys <- AnnotationDbi::keys(.get_orgdb(db$org_db), keytype = db$key_type)
    rest <- .restrict_to_testable(genes, keys)
    genes <- rest$genes
    n_testable <- rest$n_kept
    n_dropped <- rest$n_dropped
  }
  structure(genes,
            class = c("epi_gene_universe", "character"),
            n_linked = n_linked,
            n_testable = n_testable,
            n_dropped = n_dropped,
            org_db = db$org_db,
            key_type = db$key_type,
            relations = .default_enrich_relations(mark, relations),
            mark = mark,
            min_evidence_tier = min_evidence_tier,
            max_per_domain = max_per_domain,
            timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

#' Functional Enrichment of a Gene Set Against a Domain-Derived Universe
#'
#' @description Over-representation analysis (ORA) of a foreground gene set
#'   using the domain-derived testable universe. This is the low-level,
#'   object-free core: it accepts the genes returned by
#'   \code{get_domain_genes()} (or a plain gene vector) and the background from
#'   \code{get_domain_gene_universe()}. Each tested term receives an explicit
#'   log2 odds-ratio effect size with a 95\% confidence interval in addition to
#'   the usual p-value, so results remain comparable across groups of different
#'   sizes.
#'
#' @param genes Character vector of gene IDs, or a data.frame from
#'   \code{get_domain_genes()} (its \code{gene_id} column is used).
#' @param universe Character vector of background gene IDs, typically from
#'   \code{get_domain_gene_universe()}. Using all organism genes instead is
#'   intentionally not offered.
#' @param org_db,key_type OrgDb package name and keytype; see
#'   \code{get_domain_gene_universe()}.
#' @param database Character. Only "GO" is implemented; KEGG / Reactome are
#'   planned.
#' @param ontology Character. One of "BP", "MF", "CC". "ALL" is intentionally
#'   rejected to avoid redundant, hard-to-interpret output.
#' @param min_gs_size,max_gs_size Integer. Gene-set size filter. \code{max_gs_size}
#'   is automatically capped at the universe size.
#' @param p_adjust_method Character. Passed to \code{clusterProfiler::enrichGO}
#'   (default "BH").
#' @param padj_cutoff Numeric. Adjusted-p threshold used to flag
#'   \code{significant}; all tested terms are returned regardless.
#' @param simplify Logical. Run GO term redundancy reduction
#'   (\code{clusterProfiler::simplify}, "Wang" measure). Raw and simplified
#'   results are both retained.
#' @param simplify_cutoff Numeric. Semantic similarity cutoff for
#'   \code{simplify} (default 0.7).
#' @param simplify_scope Character. Which terms are passed to
#'   \code{simplify}: \code{"significant"} (default, terms with
#'   \code{p_adjust <= padj_cutoff}) or \code{"all"}. All tested terms are
#'   always returned in \code{$result}; \code{simplify} only produces the
#'   redundancy-reduced view. Simplifying the full tested set is both
#'   semantically unnecessary and slow, because the semantic similarity
#'   computation scales with the number of terms.
#' @param min_foreground Integer. Minimum number of foreground genes; below it
#'   a structured empty result is returned (never an error).
#' @return An \code{epi_enrichment} object (list) with the tidy
#'   \code{$result}, \code{$simplified}, mapping QC (\code{$mapping}), the
#'   foreground/background genes, and the original clusterProfiler object in
#'   \code{$object} for interoperability via \code{as_enrich_result()}.
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("GO.db", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   universe <- head(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,
#'                                        keytype = "ENTREZID"), 2000)
#'   fg <- intersect(c("7157", "1956", "672", "675"), universe)
#'   res <- enrich_epi_genes(fg, universe, org_db = "org.Hs.eg.db",
#'                           key_type = "ENTREZID", simplify = FALSE,
#'                           min_foreground = 3)
#'   res$mapping
#' }
#' @seealso \code{\link{enrich_epi_domains}}, \code{\link{get_domain_gene_universe}}
#' @export
enrich_epi_genes <- function(genes, universe, org_db = NULL, key_type = NULL,
                             database = "GO", ontology = "BP",
                             min_gs_size = 10, max_gs_size = 500,
                             p_adjust_method = "BH", padj_cutoff = 0.05,
                             simplify = TRUE, simplify_cutoff = 0.7,
                             simplify_scope = c("significant", "all"),
                             min_foreground = 10) {
  simplify_scope <- match.arg(simplify_scope)
  if (!identical(database, "GO")) {
    stop(sprintf("database = \"%s\" is not implemented yet; only \"GO\" is ",
                 database),
         "currently supported (KEGG / Reactome are planned).", call. = FALSE)
  }
  if (!ontology %in% c("BP", "MF", "CC")) {
    stop("ontology must be one of \"BP\", \"MF\", \"CC\". \"ALL\" is ",
         "intentionally not supported.", call. = FALSE)
  }
  if (is.data.frame(genes)) {
    if (!"gene_id" %in% colnames(genes)) {
      stop("genes data.frame must contain a 'gene_id' column ",
           "(as returned by get_domain_genes()).", call. = FALSE)
    }
    genes <- genes$gene_id
  }
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  universe <- unique(as.character(universe))
  universe <- universe[!is.na(universe) & nzchar(universe)]
  parameters <- list(
    database = database, ontology = ontology, org_db = org_db,
    key_type = key_type, min_gs_size = min_gs_size,
    max_gs_size = max_gs_size, p_adjust_method = p_adjust_method,
    padj_cutoff = padj_cutoff, simplify = simplify,
    simplify_cutoff = simplify_cutoff, simplify_scope = simplify_scope,
    min_foreground = min_foreground)

  db <- .resolve_enrich_db(se = NULL, org_db = org_db, key_type = key_type)
  .enrich_require(unique(c("clusterProfiler", db$org_db, "GO.db")))
  keys <- AnnotationDbi::keys(.get_orgdb(db$org_db), keytype = db$key_type)
  uni_rest <- .restrict_to_testable(universe, keys)
  universe <- uni_rest$genes
  fg_rest <- .restrict_to_testable(genes, universe)
  fg <- fg_rest$genes
  mapping <- list(
    n_foreground_raw = length(genes),
    n_foreground = length(fg),
    n_background_input = uni_rest$n_input,
    n_background_testable = length(universe),
    n_genes_not_testable = fg_rest$n_dropped,
    terms_tested = 0L,
    terms_significant = 0L)

  if (length(fg) < min_foreground) {
    return(.empty_epi_enrichment(
      method = "ORA", reason = "foreground_too_small",
      parameters = parameters, foreground_genes = fg,
      background_genes = universe, mapping = mapping))
  }
  if (length(universe) <= length(fg)) {
    return(.empty_epi_enrichment(
      method = "ORA", reason = "universe_not_larger_than_foreground",
      parameters = parameters, foreground_genes = fg,
      background_genes = universe, mapping = mapping))
  }

  max_gs <- min(as.integer(max_gs_size), length(universe))
  cp <- clusterProfiler::enrichGO(
    gene = fg, OrgDb = .get_orgdb(db$org_db), keyType = db$key_type,
    ont = ontology, pvalueCutoff = 1, pAdjustMethod = p_adjust_method,
    universe = universe, qvalueCutoff = 1,
    minGSSize = as.integer(min_gs_size), maxGSSize = max_gs,
    readable = FALSE)

  raw_df <- as.data.frame(cp)
  result <- .tidy_ora_result(raw_df, ontology, length(fg), length(universe),
                             padj_cutoff)
  simplified <- NULL
  simplified_object <- NULL
  simplify_reason <- NULL
  if (isTRUE(simplify) && nrow(result) > 0) {
    # Redundancy reduction is only meaningful for the significant subset; the
    # semantic similarity computation also scales with the number of terms, so
    # feeding all tested terms would be slow for no benefit.
    sig_idx <- which(result$p_adjust <= padj_cutoff)
    if (length(sig_idx) == 0) {
      simplify_reason <- "no_significant_terms"
    } else {
      cp_sig <- cp
      if (simplify_scope == "significant") {
        cp_sig@result <- cp@result[sig_idx, , drop = FALSE]
      }
      sm <- tryCatch(
        clusterProfiler::simplify(cp_sig, cutoff = simplify_cutoff,
                                  by = "p.adjust", select_fun = min,
                                  measure = "Wang"),
        error = function(e) e)
      if (inherits(sm, "error")) {
        simplify_reason <- conditionMessage(sm)
        warning("GO simplify failed: ", simplify_reason,
                ". Raw result retained; see provenance.", call. = FALSE)
      } else {
        simplified_object <- sm
        simplified <- .tidy_ora_result(as.data.frame(sm), ontology,
                                       length(fg), length(universe),
                                       padj_cutoff)
      }
    }
  }
  mapping$terms_tested <- nrow(result)
  mapping$terms_significant <- sum(result$significant, na.rm = TRUE)
  if (!is.null(simplify_reason)) mapping$simplify_note <- simplify_reason

  structure(list(
    method = "ORA",
    result = result,
    simplified = simplified,
    object = cp,
    simplified_object = simplified_object,
    foreground = list(genes = fg),
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

# Generate a non-colliding metadata key for stored results.
.unique_result_name <- function(existing, base) {
  if (is.null(existing) || !base %in% names(existing)) return(base)
  i <- 2L
  while (paste0(base, "_", i) %in% names(existing)) i <- i + 1L
  paste0(base, "_", i)
}

#' Domain-Aware Functional Enrichment on an epiPortrait Object
#'
#' @description Convenience entry point that selects foreground genes with the
#'   existing domain evidence and stores the result in the object. Two modes
#'   are supported:
#'   \itemize{
#'     \item \code{method = "ORA"}: discrete phenotype / transition selector
#'           (logical vector) -> ORA against the domain-derived universe.
#'     \item \code{method = "GSEA"}: continuous per-domain statistic
#'           (\code{score_col}) -> gene ranking -> ranked enrichment. This is
#'           the primary mode for remodeling gradients (delta intensity, delta
#'           width, limma statistics).
#'   }
#'   The background is always the linked-gene universe of the same object,
#'   never all organism genes.
#'
#' @param se A SummarizedExperiment after \code{annotate_epi_domains()}.
#' @param domains Logical vector (length nrow(se)) selecting the foreground
#'   domains for \code{method = "ORA"}, e.g.
#'   \code{get_combined_class_results(se)$Combined_Domain_Class == "Breadth-Super"}
#'   or a \code{Transition} column from \code{get_transition_results(se)}.
#' @param method Character. "ORA" or "GSEA".
#' @param score_col Character. Numeric \code{rowData} column used as the
#'   continuous ranking for \code{method = "GSEA"} (e.g.
#'   \code{"Intensity_Diff__Control_vs_Treatment__t"} from
#'   \code{analyze_differential_domains()} or
#'   \code{"log2WidthRatio__A_vs_B"} from \code{compute_width_transition()}).
#' @param relations,mark,max_per_domain,nearest_tss_cutoff_bp,min_evidence_tier
#'   Gene-link selection, see \code{\link{get_domain_gene_universe}}.
#' @param org_db,key_type OrgDb / keytype; NULL resolves from the object's
#'   annotation provenance.
#' @param database,ontology,min_gs_size,max_gs_size,p_adjust_method,padj_cutoff,simplify,simplify_cutoff,simplify_scope,min_foreground
#'   Passed to \code{\link{enrich_epi_genes}} (ORA) or
#'   \code{\link{gsea_epi_genes}} (GSEA).
#' @param seed Integer. RNG seed for GSEA.
#' @param result_name Character or NULL. Metadata key; defaults to
#'   \code{"<method>_<database>_<ontology>"} and is made unique on collision.
#' @return The SummarizedExperiment with the \code{epi_enrichment} object
#'   stored under \code{metadata(se)$enrichment[[result_name]]}. Retrieve it
#'   with \code{\link{get_epi_enrichment}}.
#' @examples
#' if (requireNamespace("clusterProfiler", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
#'     requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   data(example_se)
#'   se <- annotate_epi_domains(example_se, genome = "hg38")
#'   se <- call_super_domains(se, feature = "Intensity", verbose = FALSE)
#'   sel <- SummarizedExperiment::rowData(se)$Intensity_Domain_Type ==
#'     "Intensity_Super_Element"
#'   se <- enrich_epi_domains(se, domains = sel, simplify = FALSE,
#'                            min_foreground = 3)
#'   get_epi_enrichment(se)
#' }
#' @seealso \code{\link{enrich_epi_genes}}, \code{\link{gsea_epi_genes}},
#'   \code{\link{compare_epi_enrichment}}
#' @export
enrich_epi_domains <- function(se, domains = NULL,
                               method = c("ORA", "GSEA"),
                               score_col = NULL,
                               relations = NULL, mark = NULL,
                               max_per_domain = NULL,
                               nearest_tss_cutoff_bp =
                                 getOption("epiPortrait.nearest_tss_cutoff_bp",
                                           10000),
                               min_evidence_tier = 1L,
                               org_db = NULL, key_type = NULL,
                               database = "GO", ontology = "BP",
                               min_gs_size = 10, max_gs_size = 500,
                               p_adjust_method = "BH", padj_cutoff = 0.05,
                               simplify = TRUE, simplify_cutoff = 0.7,
                               simplify_scope = c("significant", "all"),
                               min_foreground = 10, seed = 1,
                               result_name = NULL) {
  method <- match.arg(method)
  simplify_scope <- match.arg(simplify_scope)
  db <- .resolve_enrich_db(se, org_db = org_db, key_type = key_type)
  universe <- get_domain_gene_universe(
    se, relations = relations, mark = mark, domains = NULL,
    max_per_domain = max_per_domain,
    nearest_tss_cutoff_bp = nearest_tss_cutoff_bp,
    min_evidence_tier = min_evidence_tier,
    restrict_testable = TRUE, org_db = db$org_db, key_type = db$key_type)

  if (method == "ORA") {
    if (is.null(domains)) {
      stop("domains (logical selector) is required for method = \"ORA\".",
           call. = FALSE)
    }
    if (!is.logical(domains) || length(domains) != nrow(se) ||
        anyNA(domains)) {
      stop("domains must be a logical vector of length nrow(se) with no NAs.",
           call. = FALSE)
    }
    links_fg <- .enrich_links(se, domains = domains, relations = relations,
                              mark = mark, max_per_domain = max_per_domain,
                              nearest_tss_cutoff_bp = nearest_tss_cutoff_bp,
                              min_evidence_tier = min_evidence_tier,
                              unique_genes = FALSE)
    if (nrow(links_fg) == 0) {
      res <- .empty_epi_enrichment(
        method = "ORA", reason = "no_links_in_selection",
        foreground_genes = character(), background_genes = universe,
        mapping = list(input_domains = sum(domains),
                       linked_domains = 0L,
                       unmapped_domains = sum(domains)))
    } else {
      res <- enrich_epi_genes(
        links_fg$gene_id, universe, org_db = db$org_db,
        key_type = db$key_type, database = database, ontology = ontology,
        min_gs_size = min_gs_size, max_gs_size = max_gs_size,
        p_adjust_method = p_adjust_method, padj_cutoff = padj_cutoff,
        simplify = simplify, simplify_cutoff = simplify_cutoff,
        simplify_scope = simplify_scope,
        min_foreground = min_foreground)
      res$gene_links <- links_fg
      res$foreground$domains <- rownames(se)[domains]
      res$mapping$input_domains <- sum(domains)
      res$mapping$linked_domains <- length(unique(links_fg$domain_id))
      res$mapping$unmapped_domains <-
        sum(domains) - res$mapping$linked_domains
      res$mapping$universe_linked_genes <- attr(universe, "n_linked")
    }
  } else {
    if (is.null(score_col)) {
      stop("score_col is required for method = \"GSEA\".", call. = FALSE)
    }
    ranked <- rank_epi_genes(
      se, score_col = score_col, relations = relations, mark = mark,
      max_per_domain = max_per_domain,
      nearest_tss_cutoff_bp = nearest_tss_cutoff_bp,
      min_evidence_tier = min_evidence_tier)
    res <- gsea_epi_genes(
      ranked, universe, org_db = db$org_db, key_type = db$key_type,
      database = database, ontology = ontology,
      min_gs_size = min_gs_size, max_gs_size = max_gs_size,
      p_adjust_method = p_adjust_method, padj_cutoff = padj_cutoff,
      seed = seed)
    res$mapping$input_domains <- nrow(se)
    res$mapping$linked_domains <- attr(ranked, "mapping")$n_domains_scored
    res$mapping$unmapped_domains <-
      nrow(se) - attr(ranked, "mapping")$n_domains_scored
    res$mapping$universe_linked_genes <- attr(universe, "n_linked")
  }

  base <- if (is.null(result_name)) {
    paste0(method, "_", database, "_", ontology)
  } else {
    result_name
  }
  existing <- S4Vectors::metadata(se)$enrichment
  name <- .unique_result_name(existing, base)
  res$name <- name
  if (is.null(existing)) existing <- list()
  existing[[name]] <- res
  S4Vectors::metadata(se)$enrichment <- existing
  se
}

#' Retrieve Stored Enrichment Results
#'
#' @param se A SummarizedExperiment with stored enrichment results.
#' @param name Character or NULL. Result name; NULL returns the named list of
#'   all stored objects of the requested type (consistent for both
#'   "enrichment" and "comparison").
#' @param type Character. "enrichment" for \code{epi_enrichment} objects from
#'   \code{enrich_epi_domains()}, or "comparison" for
#'   \code{epi_enrichment_comparison} objects from
#'   \code{compare_epi_enrichment()}.
#' @return A single stored object when \code{name} is given, otherwise a named
#'   list of all stored objects of the requested type (also when only one is
#'   stored, so the return type never depends on how many results happen to be
#'   present).
#' @examples
#' data(example_se)
#' # No enrichment stored yet:
#' is.null(S4Vectors::metadata(example_se)$enrichment)
#' @export
get_epi_enrichment <- function(se, name = NULL,
                               type = c("enrichment", "comparison")) {
  type <- match.arg(type)
  slot <- if (type == "enrichment") "enrichment" else "enrichment_comparison"
  enr <- S4Vectors::metadata(se)[[slot]]
  if (is.null(enr) || length(enr) == 0) {
    stop(sprintf("No %s results found. Run %s() first.", type,
                 if (type == "enrichment") "enrich_epi_domains" else
                   "compare_epi_enrichment"), call. = FALSE)
  }
  if (is.null(name)) {
    # Always a named list: callers that want "the" single result should index
    # it explicitly (or use get_epi_enrichment(se, name = <name>)).
    return(enr)
  }
  if (!name %in% names(enr)) {
    stop(sprintf("%s result '%s' not found. Available: %s", type, name,
                 paste(names(enr), collapse = ", ")), call. = FALSE)
  }
  enr[[name]]
}

#' @export
print.epi_enrichment <- function(x, ...) {
  ontology <- if (is.null(x$parameters$ontology)) "" else
    paste0(" | ontology: ", x$parameters$ontology)
  database <- if (is.null(x$parameters$database)) "" else
    paste0(" | database: ", x$parameters$database)
  cat(sprintf("<epi_enrichment> %s%s%s\n", x$method, database, ontology))
  cat(sprintf("  foreground: %d genes | background: %d genes | terms tested: %d | significant: %d\n",
              length(x$foreground$genes), x$background$n, nrow(x$result),
              sum(x$result$significant, na.rm = TRUE)))
  if (!is.null(x$mapping$reason)) {
    cat(sprintf("  no result: %s\n", x$mapping$reason))
  }
  n_not_testable <- x$mapping$n_genes_not_testable
  if (!is.null(n_not_testable) && !is.na(n_not_testable) && n_not_testable > 0) {
    cat(sprintf("  foreground genes not in annotation universe: %d\n",
                n_not_testable))
  }
  if (nrow(x$result) > 0) {
    show <- x$result[seq_len(min(5L, nrow(x$result))),
                     c("term_name", "count", "effect", "p_adjust")]
    print(show, row.names = FALSE)
  }
  invisible(x)
}
