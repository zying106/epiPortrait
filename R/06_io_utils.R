#' Filter ENCODE Blacklist Regions
#'
#' @param gr A GRanges object.
#' @param genome Character. "hg38", "hg19", "mm10", "mm9".
#' @param blacklist_path Custom BED path or NULL (uses built-in).
#' @return Filtered GRanges.
#' @import GenomicRanges
#' @import rtracklayer
#' @importFrom IRanges subsetByOverlaps
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 5000000), width = 200))
#' filter_blacklist(gr, genome = "hg38")
#' @export
filter_blacklist <- function(gr, genome = "hg38", blacklist_path = NULL) {
  if (length(gr) == 0) return(gr)
  # Blacklist resolution is INDEPENDENT of the TxDb/OrgDb resolver: blacklist
  # filtering must never require a TxDb Suggests package, and must work for any
  # genome with a built-in blacklist (hg38/hg19/mm10/mm9) or a custom BED file.
  if (!is.null(blacklist_path)) {
    if (!file.exists(blacklist_path)) stop("File not found: ", blacklist_path)
    blacklist <- rtracklayer::import(blacklist_path)
  } else {
    bl_path <- system.file("extdata",
                           paste0(genome, "-blacklist.v2.bed"),
                           package = "epiPortrait")
    if (!nzchar(bl_path) || !file.exists(bl_path)) {
      stop("No built-in blacklist for genome '", genome,
           "'. Pass blacklist_path to use a custom BED file.")
    }
    blacklist <- rtracklayer::import(bl_path)
  }
  if (is.null(blacklist)) stop("Blacklist not found for ", genome)
  # P1-4: zero shared seqlevels between input and blacklist means a naming
  # mismatch (chr1 vs 1), not "no blacklist peaks" — fail loudly instead of
  # silently keeping everything.
  shared_sl <- intersect(as.character(GenomicRanges::seqnames(gr)),
                         as.character(GenomicRanges::seqnames(blacklist)))
  if (length(shared_sl) == 0) {
    stop("No shared seqlevels between the input regions and the blacklist. ",
         "Check the chromosome naming convention (e.g. 'chr1' vs '1').")
  }
  if (length(shared_sl) < length(unique(as.character(GenomicRanges::seqnames(gr))))) {
    warning("Some input seqlevels are absent from the blacklist; blacklist ",
            "filtering is incomplete for those contigs.", call. = FALSE)
  }
  filtered_gr <- subsetByOverlaps(gr, blacklist, invert = TRUE)
  message(sprintf("Removed %d blacklist peaks. %d remaining.",
                  length(gr) - length(filtered_gr), length(filtered_gr)))
  filtered_gr
}
