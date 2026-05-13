#' Enrichment Analysis for Shape-Shifted Regions
#'
#' @description Runs GO and KEGG pathway enrichment on the target genes of
#' shape-shifted regions. Wraps \code{clusterProfiler} to provide a one-step
#' biological interpretation pipeline. By default, enrichment is run separately
#' for each Shape_Class category (when available), enabling class-specific
#' functional interpretation.
#'
#' @param shift_res A data.frame from \code{classify_shape_shift()}, or any
#'   data.frame containing a \code{SYMBOL} column (gene symbols) and optionally
#'   a \code{Shape_Class} column for stratified enrichment.
#' @param organism Character. Either \code{"human"} (default) or \code{"mouse"}.
#' @param ont Character. GO sub-ontology: \code{"BP"} (default), \code{"MF"},
#'   or \code{"CC"}. For all three, run the function separately.
#' @param p_cutoff Numeric. Adjusted P-value cutoff (default: 0.05).
#' @param q_cutoff Numeric. Q-value cutoff (default: 0.1).
#' @param split_by_class Logical. If \code{TRUE} and \code{Shape_Class} column
#'   exists, run enrichment per class (default: \code{TRUE}).
#'
#' @return A list with components:
#'   \item{go}{GO enrichment results (or list of results if split by class).}
#'   \item{kegg}{KEGG enrichment results (or list of results if split by class).}
#'   \item{by_class}{Logical indicating whether results are split by Shape_Class.}
#'
#' @export
enrich_shape_shifted <- function(shift_res,
                                  organism = c("human", "mouse"),
                                  ont = "BP",
                                  p_cutoff = 0.05,
                                  q_cutoff = 0.1,
                                  split_by_class = TRUE) {

  organism <- match.arg(organism)
  ont <- match.arg(ont, c("BP", "MF", "CC"))

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required. Install with BiocManager::install('clusterProfiler').")
  }

  org_config <- list(
    human = list(orgdb = "org.Hs.eg.db", kegg = "hsa"),
    mouse = list(orgdb = "org.Mm.eg.db", kegg = "mmu")
  )
  cfg <- org_config[[organism]]
  if (!requireNamespace(cfg$orgdb, quietly = TRUE)) {
    stop(sprintf("Install '%s' with BiocManager::install('%s').", cfg$orgdb, cfg$orgdb))
  }

  if (!"SYMBOL" %in% colnames(shift_res)) {
    stop("Input must contain a 'SYMBOL' column (run annotate_epi_peaks first).")
  }

  has_class <- split_by_class && "Shape_Class" %in% colnames(shift_res)

  .do_enrich <- function(gene_vec, label) {
    gene_vec <- unique(gene_vec[!is.na(gene_vec) & gene_vec != ""])
    if (length(gene_vec) < 5) {
      message(sprintf("Skipping '%s': fewer than 5 unique genes.", label))
      return(list(go = NULL, kegg = NULL))
    }

    message(sprintf("Running enrichment for '%s' (%d genes)...", label, length(gene_vec)))

    go_res <- clusterProfiler::enrichGO(
      gene          = gene_vec,
      OrgDb         = cfg$orgdb,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff
    )
    if (!is.null(go_res) && nrow(go_res) > 0) {
      go_res <- clusterProfiler::simplify(go_res)
    }

    kegg_res <- NULL
    tryCatch({
      entrez <- clusterProfiler::bitr(gene_vec, fromType = "SYMBOL",
                                      toType = "ENTREZID", OrgDb = cfg$orgdb)
      kegg_res <- clusterProfiler::enrichKEGG(
        gene          = entrez$ENTREZID,
        organism      = cfg$kegg,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_cutoff,
        qvalueCutoff  = q_cutoff
      )
    }, error = function(e) {
      message(sprintf("KEGG enrichment skipped for '%s': %s", label, e$message))
    })

    list(go = go_res, kegg = kegg_res)
  }

  if (has_class) {
    class_list <- split(shift_res$SYMBOL, shift_res$Shape_Class)
    # Skip Stable and classes with too few genes
    class_list <- class_list[names(class_list) != "Stable"]
    class_list <- class_list[lengths(class_list) >= 5]

    if (length(class_list) == 0) {
      stop("No class has enough genes for enrichment (>= 5).")
    }

    res <- lapply(names(class_list), function(cl) {
      .do_enrich(class_list[[cl]], cl)
    })
    names(res) <- names(class_list)

    # Flatten for return
    go_list   <- lapply(res, `[[`, "go")
    kegg_list <- lapply(res, `[[`, "kegg")

    message(sprintf("Enrichment complete for %d shape classes.", length(res)))
    return(list(go = go_list, kegg = kegg_list, by_class = TRUE))

  } else {
    res <- .do_enrich(shift_res$SYMBOL, "All shape-shifted")
    return(c(res, by_class = FALSE))
  }
}
