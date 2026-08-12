data(example_se)

test_that("plot_hockey_stick returns ggplot", {
  se_super <- call_super_domains(example_se, feature = "Intensity",
    verbose = FALSE)
  p <- plot_hockey_stick(se_super, feature = "Intensity")
  expect_s3_class(p, "ggplot")
})

test_that("plot_hockey_stick requires super_domains", {
  expect_error(
    plot_hockey_stick(example_se, feature = "Intensity"),
    "call_super_domains"
  )
})

test_that("plot_peak_track returns ggplot on real BigWig", {
  skip_on_os("windows")
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found (source checkout?)")
  samples <- data.frame(
    SampleID  = c("C1", "C2", "T1", "T2"),
    Condition = c("Control", "Control", "Treatment", "Treatment"),
    bw_path   = file.path(extdata, c("C1.bw", "C2.bw", "T1.bw", "T2.bw")))
  peaks <- rtracklayer::import(file.path(extdata, "peaks.bed"))
  se <- build_portrait_matrix(samples, consensus_peaks = peaks, workers = 1)
  p <- plot_peak_track(se, peak_id = rownames(se)[1],
    group_var = "Condition")
  expect_s3_class(p, "ggplot")
})

test_that("plot_portrait_correlation returns ggplot", {
  p <- plot_portrait_correlation(example_se, feature = "Intensity")
  expect_s3_class(p, "ggplot")
})

test_that("plot_portrait_pca returns ggplot", {
  p <- plot_portrait_pca(example_se, feature = "Intensity",
    group_var = "Condition")
  expect_s3_class(p, "ggplot")
})

test_that("per_sample stores per-sample cutoff and hockey supports it", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_sample", verbose = FALSE)
  prov <- get_call_provenance(se, "Intensity")
  expect_true(!is.null(prov$replicates))
  expect_true(is.finite(prov$replicates[["Control_1"]]$cutoff))
  p <- plot_hockey_stick(se, feature = "Intensity", group = "Control_1")
  expect_s3_class(p, "ggplot")
})
