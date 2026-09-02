
# Input-dependency and replicate-aware evidence-chain tests (design
# discussion 2026-08-13). These lock down the CURRENT behaviour — no code
# change was made:
#   * candidate-domain universe is required and comes from peak/bed files
#   * Intensity needs the candidate universe but NOT per-sample peak_path
#   * Breadth needs per-sample native peaks (peak_path); single sample without
#     it errors
#   * per_group calling is replicate-aware: per-replicate independent calls
#     exposed via get_replicate_calls(), aggregated by the support rule
#   * native peak geometry uses the UNION of peaks clipped to the domain

data(example_se)

# ---- 1. candidate universe + single-sample Intensity ------------------------
test_that("single-sample Intensity-only runs with candidate universe only", {
  # rtracklayer local BigWig I/O fails on Windows ("UCSC library operation failed");
  # rtracklayer#52/#62/#128/#151. Narrow skip: this block only. epiPortrait logic
  # that does not read BigWigs is still executed on Windows.
  skip_on_os("windows")
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")

  # candidate universe read from the sample's own peak file (as documented)
  my_domains <- rtracklayer::import(file.path(extdata, "peaks.bed"))
  samples1 <- data.frame(
    SampleID = "S1", Condition = "Control",
    bw_path = file.path(extdata, "C1.bw"))
  se <- build_portrait_matrix(samples1, consensus_peaks = my_domains,
                              workers = 1)
  # no peak_path -> native geometry assays exist but are NA
  expect_true(all(c("Intensity", "SignalDispersion") %in%
                    SummarizedExperiment::assayNames(se)))
  expect_true(all(is.na(
    SummarizedExperiment::assay(se, "NativeOccupiedWidth"))))

  # Intensity call works per_sample (no replicate support with n=1)
  se <- call_super_domains(se, feature = "Intensity",
                           mode = "per_sample", verbose = FALSE)
  expect_true("Intensity_Call__S1" %in%
                colnames(SummarizedExperiment::rowData(se)))
})

test_that("single-sample Breadth errors without peak_path (current behaviour)", {
  # rtracklayer local BigWig I/O fails on Windows ("UCSC library operation failed");
  # rtracklayer#52/#62/#128/#151. Narrow skip: this block only. epiPortrait logic
  # that does not read BigWigs is still executed on Windows.
  skip_on_os("windows")
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")

  my_domains <- rtracklayer::import(file.path(extdata, "peaks.bed"))
  samples1 <- data.frame(
    SampleID = "S1", Condition = "Control",
    bw_path = file.path(extdata, "C1.bw"))
  se <- build_portrait_matrix(samples1, consensus_peaks = my_domains,
                              workers = 1)
  # Breadth requires per-sample native peaks: even with n=1, peak_path must be
  # supplied (candidate universe is NOT reused as native peaks).
  expect_error(
    call_super_domains(se, feature = "Breadth", mode = "per_group",
                       group_var = "Condition", verbose = FALSE),
    "No native peak calls found")
})

# ---- 2. multi-sample Intensity does not need peak_path ----------------------
test_that("multi-sample Intensity runs without per-sample native peaks", {
  # rtracklayer local BigWig I/O fails on Windows ("UCSC library operation failed");
  # rtracklayer#52/#62/#128/#151. Narrow skip: this block only. epiPortrait logic
  # that does not read BigWigs is still executed on Windows.
  skip_on_os("windows")
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")

  # candidate universe from consensus of per-sample peak files
  peak_list <- list(
    rep1 = rtracklayer::import(file.path(extdata, "peaks.bed")),
    rep2 = rtracklayer::import(file.path(extdata, "peaks.bed")))
  consensus <- get_consensus_peaks(peak_list, min_reps = 1)
  samples <- data.frame(
    SampleID = c("C1", "C2"), Condition = c("C", "C"),
    bw_path = file.path(extdata, c("C1.bw", "C2.bw")))
  se <- build_portrait_matrix(samples, consensus_peaks = consensus,
                              workers = 1)
  se <- call_super_domains(se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           verbose = FALSE)
  expect_true("Intensity_Call__C" %in%
                colnames(SummarizedExperiment::rowData(se)))
  # no per-sample peaks -> no native_peaks metadata
  expect_true(is.null(S4Vectors::metadata(se)$native_peaks))
})

# ---- 3. replicate-aware evidence chain (Intensity) --------------------------
test_that("get_replicate_calls exposes per-replicate independent calls", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           verbose = FALSE)
  m <- get_replicate_calls(se, feature = "Intensity", group = "Control")
  # domain x replicate matrix, one column per Control replicate
  expect_equal(nrow(m), nrow(example_se))
  expect_equal(colnames(m), c("Control_1", "Control_2", "Control_3"))
  # entries are the per-replicate independent call labels
  expect_true(all(unique(as.vector(m)) %in%
                    c("Intensity_Super_Element", "Intensity_Typical", NA)))
})

test_that("group support equals fraction of super replicates", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           verbose = FALSE)
  m <- get_replicate_calls(se, feature = "Intensity", group = "Control")
  n_super <- unname(rowSums(m == "Intensity_Super_Element", na.rm = TRUE))
  support <- n_super / ncol(m)
  expect_equal(
    unname(SummarizedExperiment::rowData(se)$Intensity_Support__Control),
    support)
  # group call consistent with majority rule (n=3 -> need >= 2 super)
  grp_call <- SummarizedExperiment::rowData(se)$Intensity_Call__Control
  expect_true(all((grp_call == "Intensity_Super_Element") == (n_super >= 2)))
})

# ---- 4. native peak geometry union semantics --------------------------------
test_that("native_peak_geometry uses union of clipped peaks", {
  domains <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(100, 1000))
  # two non-overlapping peaks inside the domain (1-based closed intervals)
  peaks_disjoint <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(200, 600), c(300, 700)))
  g1 <- epiPortrait:::.native_peak_geometry(peaks_disjoint, domains)
  # union width = 101 + 101 = 202
  expect_equal(g1$NWocc[1], 202)
  expect_equal(g1$NWmax[1], 101)
  expect_equal(g1$NWcnt[1], 2)

  # two overlapping peaks: union < sum of widths
  peaks_overlap <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(200, 250), c(400, 450)))
  g2 <- epiPortrait:::.native_peak_geometry(peaks_overlap, domains)
  expect_equal(g2$NWocc[1], 251)   # 200-450 reduced
  expect_equal(g2$NWcnt[1], 2)

  # a peak extending outside the domain is clipped to it
  peaks_outside <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(50, 500))
  g3 <- epiPortrait:::.native_peak_geometry(peaks_outside, domains)
  expect_equal(g3$NWocc[1], 401)   # clipped to 100-500
  expect_equal(g3$NWmax[1], 451)   # full width of the native peak
})
test_that("build_portrait_matrix rejects mixed declared genome assemblies", {
  ss <- data.frame(
    SampleID = c("s1", "s2"),
    Condition = c("A", "B"),
    bw_path = c(tempfile(fileext = ".bw"), tempfile(fileext = ".bw")),
    Genome = c("hg38", "hg19"))
  file.create(ss$bw_path)
  on.exit(unlink(ss$bw_path), add = TRUE)
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  expect_error(build_portrait_matrix(ss, gr), "same genome assembly")
})
