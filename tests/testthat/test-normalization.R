
data(example_se)

test_that("normalize_portrait returns SummarizedExperiment for all methods", {
  for (method in c("TotalSignal", "Quantile", "None")) {
    se_norm <- suppressWarnings(normalize_portrait(example_se, method = method))
    expect_s4_class(se_norm, "SummarizedExperiment")
    expect_setequal(assayNames(se_norm), assayNames(example_se))
    expect_identical(
      assay(se_norm, "SignalDispersion"),
      assay(example_se, "SignalDispersion")
    )
  }
})

test_that("normalize_portrait TotalSignal sets scaling factors", {
  expect_warning(
    se_norm <- normalize_portrait(example_se, method = "TotalSignal"),
    "conserved total signal"
  )
  expect_true("ScalingFactor" %in% colnames(colData(se_norm)))
  expect_true(all(is.finite(colData(se_norm)$ScalingFactor)))
  prov <- S4Vectors::metadata(se_norm)$normalization
  expect_identical(prov$method, "TotalSignal")
  expect_true(prov$applied)
  expect_equal(unname(prov$scaling_factors),
               unname(colData(se_norm)$ScalingFactor))
  expect_equal(unname(prov$output_sample_totals),
               rep(mean(prov$input_sample_totals), ncol(se_norm)))
})

test_that("normalize_portrait rejects invalid method", {
  expect_error(normalize_portrait(example_se, method = "INVALID"))
})

test_that("normalize_portrait excludes count-based TMM", {
  expect_error(normalize_portrait(example_se, method = "TMM"),
               "Invalid method")
})

test_that("normalize_portrait warns and records Quantile normalization", {
  expect_warning(
    se_norm <- normalize_portrait(example_se, method = "Quantile"),
    "erase genuine global biological shifts"
  )
  prov <- S4Vectors::metadata(se_norm)$normalization
  expect_identical(prov$method, "Quantile")
  expect_true(prov$applied)
  expect_null(prov$scaling_factors)
  expect_false(prov$SignalDispersion_modified)
})

test_that("normalize_portrait prevents repeated post-hoc normalization", {
  se_norm <- suppressWarnings(
    normalize_portrait(example_se, method = "TotalSignal")
  )
  expect_error(
    normalize_portrait(se_norm, method = "Quantile"),
    "already marked as normalized"
  )
})

test_that("normalize_portrait records the no-normalization assumption", {
  se_none <- normalize_portrait(example_se, method = "None")
  prov <- S4Vectors::metadata(se_none)$normalization
  expect_identical(prov$method, "None")
  expect_false(prov$applied)
})

test_that("normalize_portrait Z-score is removed (cannot corrupt canonical Intensity)", {
  # Freeze audit 2026-08-10: Z-score is a row-wise display/clustering transform
  # that destroys the cross-domain magnitude ranking required by Super calling,
  # so it was removed from normalize_portrait().
  expect_error(normalize_portrait(example_se, method = "Z-score"))
})
