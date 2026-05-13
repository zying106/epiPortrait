context("Consensus Peaks")

test_that("get_consensus_peaks builds consensus with reproducibility filter", {
  gr1 <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = c(100, 300), end = c(200, 400)))
  gr2 <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = c(150, 350), end = c(250, 400)))
  gr3 <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = c(100, 500), end = c(200, 600)))

  gr_list <- GRangesList(gr1, gr2, gr3)

  # min_reps = 2: need support from at least 2 samples
  res <- get_consensus_peaks(gr_list, min_reps = 2)
  expect_s4_class(res, "GRanges")
  expect_true("Support_Reps" %in% colnames(mcols(res)))
  expect_true(all(mcols(res)$Support_Reps >= 2))

  # min_reps = 3: need all 3 samples
  res3 <- get_consensus_peaks(gr_list, min_reps = 3)
  expect_true(all(mcols(res3)$Support_Reps >= 3))
  expect_true(length(res3) <= length(res))
})

test_that("get_consensus_peaks handles edge cases", {
  res <- get_consensus_peaks(GRangesList(), min_reps = 1)
  expect_equal(length(res), 0)

  single <- GRanges("chr1", IRanges::IRanges(1, 10))
  res_single <- get_consensus_peaks(GRangesList(single), min_reps = 1)
  expect_equal(length(res_single), 1)

  expect_error(get_consensus_peaks(matrix(1:10)), "must be a GRangesList")
})
