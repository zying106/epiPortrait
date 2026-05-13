context("Annotation")

test_that("annotate_epi_peaks handles empty input", {
  expect_null(suppressWarnings(annotate_epi_peaks(data.frame(
    seqnames = character(0), start = integer(0), end = integer(0)
  ), genome = "hg38")))
})

test_that("annotate_epi_peaks auto-converts data.frame", {
  gr <- GRanges("chr1", IRanges::IRanges(c(100, 500), c(200, 600)))
  res <- annotate_epi_peaks(gr, genome = "hg38")
  expect_s3_class(res, "data.frame")
})

test_that("annotate_epi_peaks validates custom genome", {
  gr <- GRanges("chr1", IRanges::IRanges(100, 200))
  expect_error(annotate_epi_peaks(gr, genome = "custom"),
               "custom_txdb")
  expect_error(annotate_epi_peaks(gr, genome = "fake_genome"),
               "Unsupported")
})
