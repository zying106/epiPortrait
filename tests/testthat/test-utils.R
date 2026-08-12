library(epiPortrait)
library(GenomicRanges)

test_that("stitch_epi_peaks merges proximal peaks", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(100, 250, 5000), end = c(200, 350, 5200)))
  res <- stitch_epi_peaks(gr, stitch_distance = 1000)
  expect_s4_class(res, "GRanges")
  expect_true("Constituent_Peaks" %in% colnames(mcols(res)))
  expect_true(length(res) < length(gr))
})

test_that("stitch_epi_peaks validates inputs", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(100, 250), end = c(200, 350)))
  expect_error(stitch_epi_peaks("not a granges", stitch_distance = 1000),
               "GRanges")
  expect_error(stitch_epi_peaks(gr, stitch_distance = -1),
               "non-negative")
})

test_that("stitch_epi_peaks handles single peaks", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  res <- stitch_epi_peaks(gr, stitch_distance = 1000)
  expect_equal(length(res), 1)
  expect_equal(width(res)[1], 10)
})
test_that("stitch_epi_peaks merges peaks with gap exactly equal to distance", {
  # reduce() merges only when gap < min.gapwidth; +1 fix makes gap <= distance exact.
  # peaks at 1-100 and 1101-1200 (gap = 1000) with stitch_distance = 1000 should merge.
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 1101), c(100, 1200)))
  res <- stitch_epi_peaks(gr, stitch_distance = 1000)
  expect_equal(length(res), 1)
  expect_equal(start(res), 1)
  expect_equal(end(res), 1200)
})
