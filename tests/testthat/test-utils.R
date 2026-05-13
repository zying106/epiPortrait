context("Utility Functions")

test_that("stitch_epi_peaks merges proximal peaks", {
  gr <- GRanges("chr1", IRanges::IRanges(start = c(100, 250, 5000), end = c(200, 350, 5200)))
  res <- stitch_epi_peaks(gr, stitch_distance = 1000)
  expect_s4_class(res, "GRanges")
  expect_true("Constituent_Peaks" %in% colnames(mcols(res)))
  expect_true(length(res) < length(gr))
})

test_that("stitch_epi_peaks validates inputs", {
  gr <- GRanges("chr1", IRanges::IRanges(1, 10))
  expect_error(stitch_epi_peaks(data.frame()), "must be a GRanges")
  expect_error(stitch_epi_peaks(gr, stitch_distance = -5), "non-negative")
})

test_that("normalize_portrait handles all methods", {
  data(example_se)
  methods <- c("TotalSignal", "Quantile", "Z-score", "None")
  for (m in methods) {
    expect_s4_class(normalize_portrait(example_se, method = m), "SummarizedExperiment")
  }
  # TMM requires edgeR and force_TMM for continuous data
  if (requireNamespace("edgeR", quietly = TRUE)) {
    expect_s4_class(normalize_portrait(example_se, method = "TMM", force_TMM = TRUE),
                    "SummarizedExperiment")
  }
  # TMM without force_TMM should error
  if (requireNamespace("edgeR", quietly = TRUE)) {
    expect_error(normalize_portrait(example_se, method = "TMM", force_TMM = FALSE),
                 "continuous")
  }
  expect_error(normalize_portrait(example_se, method = "Invalid"))
})

test_that("call_super_domains identifies super elements", {
  data(example_se)
  se <- call_super_domains(example_se, feature = "Intensity")
  expect_true("Intensity_Domain_Type" %in% colnames(rowData(se)))
  expect_true("Intensity_Rank" %in% colnames(rowData(se)))
})

test_that("export_shifted_bed writes valid BED", {
  res_df <- data.frame(
    seqnames = "chr1", start = 100, end = 200,
    Peak_ID = "Peak_1", Shape_Shift_Score = 0.8, adj.P.Val = 0.01,
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".bed")
  export_shifted_bed(res_df, file = tmp, fdr_cutoff = 0.05)
  bed <- read.table(tmp, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(bed[1, 2], 99)  # 0-based start
  unlink(tmp)
})
