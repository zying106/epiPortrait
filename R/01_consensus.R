#' Generate Consensus Peak Set for epiPortrait
#'
#' @description Merges multiple peak sets into a unified genomic coordinate system
#' while applying a reproducibility filter.
#' @param gr_list A GRangesList or list of GRanges objects containing peaks from all samples.
#' @param min_reps Minimum number of replicates a peak must appear in to be retained (default: 2).
#' @return A filtered GRanges object representing the consensus peaks.
#' @import GenomicRanges
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
