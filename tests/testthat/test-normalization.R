
data(example_se)

test_that("normalize_portrait returns SummarizedExperiment for all methods", {
  for (method in c("TotalSignal", "Quantile", "None")) {
    se_norm <- normalize_portrait(example_se, method = method)
    expect_s4_class(se_norm, "SummarizedExperiment")
    expect_setequal(assayNames(se_norm), assayNames(example_se))
  }
})

test_that("normalize_portrait TotalSignal sets scaling factors", {
  se_norm <- normalize_portrait(example_se, method = "TotalSignal")
  expect_true("ScalingFactor" %in% colnames(colData(se_norm)))
  expect_true(all(is.finite(colData(se_norm)$ScalingFactor)))
})

test_that("normalize_portrait rejects invalid method", {
  expect_error(normalize_portrait(example_se, method = "INVALID"))
})

test_that("normalize_portrait TMM guards force_TMM", {
  skip_if_not_installed("edgeR")
  expect_error(
    normalize_portrait(example_se, method = "TMM", force_TMM = FALSE),
    "force_TMM"
  )
})

test_that("normalize_portrait Z-score is removed (cannot corrupt canonical Intensity)", {
  # Freeze audit 2026-08-10: Z-score is a row-wise display/clustering transform
  # that destroys the cross-domain magnitude ranking required by Super calling,
  # so it was removed from normalize_portrait().
  expect_error(normalize_portrait(example_se, method = "Z-score"))
})
