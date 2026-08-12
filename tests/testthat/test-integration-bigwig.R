# Integration test running the full core pipeline on the tiny synthetic BigWig
# dataset shipped in inst/extdata (Bioconductor reviewer expectation: the real
# BigWig import path must be exercised by unit tests).

test_that("full pipeline runs on tiny real BigWigs", {
  skip_on_os("windows")  # rtracklayer BigWig is unsupported on Windows
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found (source checkout?)")

  samples <- data.frame(
    SampleID  = c("C1", "C2", "T1", "T2"),
    Condition = c("Control", "Control", "Treatment", "Treatment"),
    Replicate = c(1, 2, 1, 2),
    bw_path   = file.path(extdata, c("C1.bw", "C2.bw", "T1.bw", "T2.bw")),
    peak_path = file.path(extdata, c("C1_peaks.bed", "C2_peaks.bed",
                                     "T1_peaks.bed", "T2_peaks.bed"))
  )
  peaks <- rtracklayer::import(file.path(extdata, "peaks.bed"))

  se <- build_portrait_matrix(samples, consensus_peaks = peaks, workers = 1)
  expect_true(all(c("Intensity", "SignalDispersion",
                    "NativeMaxPeakWidth", "NativeOccupiedWidth",
                    "NativePeakCount") %in% assayNames(se)))
  expect_true("IntervalWidth" %in% colnames(rowData(se)))
  expect_true("Replicate" %in% colnames(colData(se)))
  expect_equal(nrow(se), length(peaks))
  expect_equal(ncol(se), 4)

  # Native peaks are stored for downstream Breadth-Super calling.
  expect_false(is.null(S4Vectors::metadata(se)$native_peaks))

  # SignalDispersion is finite for domains with signal (NA allowed for
  # zero-signal domains)
  d <- as.matrix(assay(se, "SignalDispersion"))
  expect_true(all(is.finite(d[!is.na(d)])))
  expect_true(all(d[!is.na(d)] >= 0))

  # Intensity is finite and non-negative (some domains may have no signal
  # overlap and hence 0 intensity)
  i <- as.matrix(assay(se, "Intensity"))
  expect_true(all(is.finite(i)))
  expect_true(all(i >= 0))
})

test_that("Breadth-Super calling works end-to-end on real BigWigs", {
  skip_on_os("windows")  # rtracklayer BigWig is unsupported on Windows
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")

  samples <- data.frame(
    SampleID  = c("C1", "C2", "T1", "T2"),
    Condition = c("Control", "Control", "Treatment", "Treatment"),
    bw_path   = file.path(extdata, c("C1.bw", "C2.bw", "T1.bw", "T2.bw")),
    peak_path = file.path(extdata, c("C1_peaks.bed", "C2_peaks.bed",
                                     "T1_peaks.bed", "T2_peaks.bed"))
  )
  peaks <- rtracklayer::import(file.path(extdata, "peaks.bed"))
  se <- build_portrait_matrix(samples, consensus_peaks = peaks, workers = 1)

  # Breadth-Super is a peak-level call (decoupled from consensus).
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition",
                           verbose = FALSE)
  expect_true("Breadth_Call__Control" %in% colnames(rowData(se)))
  expect_true("Breadth_Call__Treatment" %in% colnames(rowData(se)))
  # The broad native peak (domain 3) must be called Breadth-Super.
  expect_true("Breadth_Super_Element" %in% rowData(se)$Breadth_Call__Control)
  # Peak-level provenance table is stored.
  expect_false(is.null(S4Vectors::metadata(se)$breadth_peak_calls))
  expect_false(is.null(S4Vectors::metadata(se)$breadth_peak_mapping))
})

test_that("condition-aware consensus and combine work on real BigWigs", {
  skip_on_os("windows")  # rtracklayer BigWig is unsupported on Windows
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")

  samples <- data.frame(
    SampleID  = c("C1", "C2", "T1", "T2"),
    Condition = c("Control", "Control", "Treatment", "Treatment"),
    bw_path   = file.path(extdata, c("C1.bw", "C2.bw", "T1.bw", "T2.bw")),
    peak_path = file.path(extdata, c("C1_peaks.bed", "C2_peaks.bed",
                                     "T1_peaks.bed", "T2_peaks.bed"))
  )
  peaks <- rtracklayer::import(file.path(extdata, "peaks.bed"))
  se <- build_portrait_matrix(samples, consensus_peaks = peaks, workers = 1)

  se <- call_super_domains(se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  expect_true("Intensity_Support__Control" %in% colnames(rowData(se)))
  expect_true("Intensity_Support__Treatment" %in% colnames(rowData(se)))

  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition",
                           verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  expect_true("Combined_Class__Control" %in% colnames(rowData(se)))
  expect_true(all(rowData(se)$Combined_Class__Control %in%
                    c("Intensity-Super", "Breadth-Super", "Dual-Super",
                      "Typical", "Uncertain")))
})
