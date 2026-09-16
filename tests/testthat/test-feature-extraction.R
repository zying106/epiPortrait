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


.custom_feature_fixture <- function() {
  skip_on_os("windows")
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")
  list(
    samples = data.frame(
      SampleID = "C1",
      Condition = "Control",
      bw_path = file.path(extdata, "C1.bw")
    ),
    domains = rtracklayer::import(file.path(extdata, "peaks.bed"))
  )
}


test_that("custom_features receive numeric coverage and create a new assay", {
  x <- .custom_feature_fixture()
  se <- build_portrait_matrix(
    x$samples,
    x$domains,
    workers = 1,
    custom_features = list(
      MeanSignal = function(values) {
        if (!is.numeric(values) || inherits(values, "NumericList")) {
          stop("expected one numeric coverage vector")
        }
        mean(values)
      }
    )
  )

  expect_true("MeanSignal" %in% SummarizedExperiment::assayNames(se))
  expect_equal(
    as.numeric(SummarizedExperiment::assay(se, "MeanSignal")),
    as.numeric(SummarizedExperiment::assay(se, "Intensity")) /
      GenomicRanges::width(SummarizedExperiment::rowRanges(se))
  )
})


test_that("custom_features cannot overwrite built-in or special features", {
  x <- .custom_feature_fixture()
  expect_error(
    build_portrait_matrix(
      x$samples,
      x$domains,
      workers = 1,
      custom_features = list(Intensity = function(values) mean(values))
    ),
    "reserved epiPortrait features: Intensity"
  )
  expect_error(
    build_portrait_matrix(
      x$samples,
      x$domains,
      workers = 1,
      custom_features = list(Breadth = function(values) mean(values))
    ),
    "reserved epiPortrait features: Breadth"
  )
})


test_that("custom_features require valid unique names and functions", {
  x <- .custom_feature_fixture()

  duplicate_features <- list(function(values) mean(values),
                             function(values) max(values))
  names(duplicate_features) <- c("ShapeScore", "ShapeScore")
  expect_error(
    build_portrait_matrix(x$samples, x$domains,
                          custom_features = duplicate_features),
    "must be unique"
  )
  expect_error(
    build_portrait_matrix(
      x$samples,
      x$domains,
      custom_features = list("bad feature" = function(values) mean(values))
    ),
    "syntactically valid R names"
  )
  expect_error(
    build_portrait_matrix(
      x$samples,
      x$domains,
      custom_features = list(ShapeScore = 1)
    ),
    "must be a function"
  )
})


test_that("invalid custom feature output fails even when fail_action is drop", {
  x <- .custom_feature_fixture()
  expect_error(
    build_portrait_matrix(
      x$samples,
      x$domains,
      workers = 1,
      custom_features = list(BadScore = function(values) c(1, 2)),
      fail_action = "drop"
    ),
    "failed during custom feature evaluation.*must return one finite numeric",
    fixed = FALSE
  )
})
