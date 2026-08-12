library(epiPortrait)
library(GenomicRanges)

test_that("get_consensus_peaks builds consensus with reproducibility filter", {
  gr1 <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = c(100, 300), end = c(200, 400)))
  gr2 <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = c(150, 350), end = c(250, 400)))
  gr3 <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = c(100, 500), end = c(200, 600)))

  gr_list <- GRangesList(gr1, gr2, gr3)

  # min_reps = 2: need support from at least 2 samples
  res <- get_consensus_peaks(gr_list, min_reps = 2)
  expect_s4_class(res, "GRanges")
  expect_true(length(res) < length(gr_list))
  expect_true("Support_Reps" %in% colnames(mcols(res)))
})

test_that("get_consensus_peaks handles edge cases", {
  expect_error(get_consensus_peaks("not a granges"), "GRangesList")
  expect_error(get_consensus_peaks(data.frame()), "GRangesList")
})