#' Generate Consensus Domain Universe for epiPortrait
#'
#' @description Merges multiple sample peak/domain sets into a unified genomic
#' coordinate system while applying a reproducibility filter. This is
#' CONSENSUS / domain-universe construction, not a peak caller (epiPortrait
#' does not call MACS2/SICER-type peaks).
#'
#' @details
#' **Union-support chaining semantics.** All sample peaks are first
#' \code{reduce()}-merged into connected intervals, then each reduced interval
#' counts how many sample peak sets contribute at least one overlapping peak
#' (\code{Support_Reps}). Connected overlapping peaks can chain into a wider
#' reduced interval; \code{Support_Reps} is an interval-level replicate
#' contribution count, NOT base-pair-wise support across replicates. A
#' \code{Support_Reps = 3} interval does not imply every base pair is covered by
#' all three replicates.
#'
#' @param gr_list A GRangesList or list of GRanges objects containing peaks/domains from all samples.
#' @param min_reps Minimum number of replicate peak sets contributing at least
#'   one overlap to a reduced consensus interval for it to be retained
#'   (default: 2).
#' @return A filtered GRanges object representing the consensus domain universe.
#' @import GenomicRanges
#' @examples
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 5000), width = 200))
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(150, 5200), width = 200))
#' get_consensus_peaks(list(rep1 = gr1, rep2 = gr2), min_reps = 2)
#' @export
get_consensus_peaks <- function(gr_list, min_reps = 2) {

  if (!is(gr_list, "GRangesList") && !is(gr_list, "list")) {
    stop("Input must be a GRangesList or a list of GRanges.")
  }

  grl <- if (is(gr_list, "GRangesList")) gr_list else GRangesList(gr_list)

  all_peaks <- unlist(grl)
  merged_peaks <- reduce(all_peaks)

  overlap_matrix <- vapply(grl, function(sample_gr) {
    countOverlaps(merged_peaks, sample_gr) > 0
  }, logical(length(merged_peaks)))

  # When only one sample, vapply returns a vector; coerce to matrix
  if (!is.matrix(overlap_matrix)) {
    overlap_matrix <- matrix(overlap_matrix, ncol = 1)
  }

  support_counts <- rowSums(overlap_matrix)

  # Filter high-confidence regions and store metadata
  consensus_peaks <- merged_peaks[support_counts >= min_reps]
  mcols(consensus_peaks)$Support_Reps <- support_counts[support_counts >= min_reps]

  return(consensus_peaks)
}
