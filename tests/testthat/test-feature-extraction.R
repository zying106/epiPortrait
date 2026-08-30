library(epiPortrait)

test_that("build_portrait_matrix validates input", {
  expect_error(
    build_portrait_matrix(data.frame(x = 1:3), GenomicRanges::GRanges("chr1", IRanges(1, 10))),
    "bw_path"
  )
})

test_that("build_portrait_matrix produces valid SE object", {
  data(example_se)
  expect_s4_class(example_se, "SummarizedExperiment")
  expect_true(all(c("Intensity", "SignalDispersion",
                    "NativeMaxPeakWidth", "NativeOccupiedWidth",
                    "NativePeakCount") %in% assayNames(example_se)))
  expect_true("IntervalWidth" %in% colnames(rowData(example_se)))
  expect_true(all(c("SampleID", "Condition") %in% colnames(colData(example_se))))
  expect_false(is.null(S4Vectors::metadata(example_se)$native_peaks))
})
