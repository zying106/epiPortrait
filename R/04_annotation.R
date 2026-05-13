#' Annotate epiPortrait Geometric Features with Target Genes
#'
#' @description Maps consensus peaks or shape-shifted regions (e.g., broad H3K27ac domains)
#' to their nearest genomic features and target genes using ChIPseeker. Supports smart
#' lazy-loading for common genomes and custom TxDb objects.
#'
#' @param gr A GRanges object or a data.frame with 'seqnames', 'start', 'end' columns
#' (such as the output from test_global_shape_shift).
#' @param genome A character string specifying the genome shortcut. Options are "hg38",
#' "hg19", "mm10", or "custom".
#' @param custom_txdb A TxDb object. Required only if genome = "custom".
#' @param custom_annoDb A character string specifying the annotation package name
#' (e.g., "org.At.tair.db"). Required only if genome = "custom".
#' @param tss_region A numeric vector of length 2 defining the TSS region (default: c(-3000, 3000)).
#' @param return_format Character. Either "dataframe" (default, merges annotations with your 4D features)
#' or "csAnno" (returns the raw ChIPseeker S4 object for plotting).
#'
#' @return A data.frame or a csAnno object depending on `return_format`.
#'
#' @examples
#' \dontrun{
#' # Example: Annotating significantly shifted domains
#' shift_res <- test_global_shape_shift(portrait_se, group_var = "Condition")
#' sig_peaks <- shift_res[shift_res$adj.P.Val < 0.05, ]
#'
#' # Get detailed annotation table
#' anno_df <- annotate_epi_peaks(sig_peaks, genome = "hg38")
#'
#' # Get csAnno object for plotting (e.g., ChIPseeker::plotAnnoPie(anno_obj))
#' anno_obj <- annotate_epi_peaks(sig_peaks, genome = "hg38", return_format = "csAnno")
#' }
#'
#' @importFrom ChIPseeker annotatePeak
#' @import GenomicRanges
#' @export
annotate_epi_peaks <- function(gr, genome = "hg38",
                               custom_txdb = NULL, custom_annoDb = NULL,
                               tss_region = c(-3000, 3000),
                               return_format = c("dataframe", "csAnno")) {

  return_format <- match.arg(return_format)

  # 1. Gracefully handle empty input
  if ((is.data.frame(gr) && nrow(gr) == 0) || (is(gr, "GRanges") && length(gr) == 0)) {
    warning("Input regions are empty (0 rows). Returning NULL.")
    return(NULL)
  }

  # 2. Auto-convert data.frame to GRanges
  if (is.data.frame(gr)) {
    req_cols <- c("seqnames", "start", "end")
    if (!all(req_cols %in% colnames(gr))) {
      stop("Input data.frame must contain 'seqnames', 'start', and 'end' columns.")
    }
    gr_obj <- GenomicRanges::GRanges(seqnames = gr$seqnames,
                                     ranges = IRanges::IRanges(start = gr$start, end = gr$end))

    # Preserve epiPortrait 4D feature columns (Shape_Shift_Score, adj.P.Val, etc.)
    meta_cols <- gr[, !colnames(gr) %in% req_cols, drop = FALSE]
    if (ncol(meta_cols) > 0) {
      GenomicRanges::mcols(gr_obj) <- meta_cols
    }
  } else if (is(gr, "GRanges")) {
    gr_obj <- gr
  } else {
    stop("Input must be a GRanges object or a coordinate data.frame.")
  }

  # 3. Smart genome lookup with lazy-loading TxDb
  message(sprintf("Initializing epiPortrait annotation module for %s...", genome))

  if (genome %in% c("hg38", "hg19", "mm10")) {
    db <- .get_txdb(genome)
    txdb <- db$txdb
    annoDb <- db$annoDb

  } else if (genome == "custom") {
    if (is.null(custom_txdb) || is.null(custom_annoDb)) {
      stop("When genome='custom', you must explicitly provide 'custom_txdb' and 'custom_annoDb'.")
    }
    txdb <- custom_txdb
    annoDb <- custom_annoDb
  } else {
    stop("Unsupported genome shortcut. Use 'hg38', 'hg19', 'mm10', or set genome='custom'.")
  }

  # 4. Run ChIPseeker annotation engine
  message("Mapping shifted domains to nearest features...")
  peakAnno <- ChIPseeker::annotatePeak(gr_obj,
                                       TxDb = txdb,
                                       annoDb = annoDb,
                                       tssRegion = tss_region,
                                       verbose = FALSE)

  # 5. Return requested format
  if (return_format == "csAnno") {
    message("Returning raw csAnno object. You can now use ChIPseeker::plotAnnoPie() or plotAnnoBar().")
    return(peakAnno)
  } else {
    message("Annotation complete! Returning integrated data.frame.")
    return(as.data.frame(peakAnno))
  }
}
