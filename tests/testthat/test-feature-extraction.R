context("Feature Extraction")

test_that("build_portrait_matrix validates input", {
  expect_error(
    build_portrait_matrix(data.frame(x = 1:3), GRanges("chr1", IRanges(1, 10))),
    "SampleID"
  )
})

test_that("build_portrait_matrix produces valid SE object", {
  data(example_se)
  expect_s4_class(example_se, "SummarizedExperiment")
  expect_setequal(assayNames(example_se), c("Intensity", "Height", "Skewness", "Width"))
  expect_true(all(c("SampleID", "Condition") %in% colnames(colData(example_se))))
})

test_that("example_se has known shape-shift classes", {
  data(example_se)
  expect_true("True_Class" %in% colnames(mcols(rowRanges(example_se))))
  classes <- table(mcols(rowRanges(example_se))$True_Class)
  expect_true("Stable" %in% names(classes))
  expect_true(sum(classes[names(classes) != "Stable"]) > 50)
})
