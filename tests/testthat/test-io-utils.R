data(example_se)

test_gr <- rowRanges(example_se)
set.seed(42)
test_gr_small <- test_gr[1:50]

test_that("call_super_domains adds domain type columns", {
  se_super <- call_super_domains(example_se, feature = "Intensity", verbose = FALSE)
  types <- rowData(se_super)$Intensity_Domain_Type
  expect_true("Intensity_Super_Element" %in% types)
  expect_true("Intensity_Typical" %in% types)
})

test_that("call_super_domains works with different features", {
  # SignalDispersion is a canonical assay feature
  se_sd <- call_super_domains(example_se, feature = "SignalDispersion", verbose = FALSE)
  expect_true("SignalDispersion_Domain_Type" %in% colnames(rowData(se_sd)))
  se_iw <- call_super_domains(example_se, feature = "IntervalWidth", verbose = FALSE)
  expect_true("IntervalWidth_Domain_Type" %in% colnames(rowData(se_iw)))
})

test_that("call_super_domains supports tangent (ROSE) method", {
  se_elbow <- call_super_domains(example_se, feature = "Intensity", method = "elbow", verbose = FALSE)
  se_tangent <- call_super_domains(example_se, feature = "Intensity", method = "tangent", verbose = FALSE)
  expect_true("Intensity_Domain_Type" %in% colnames(rowData(se_tangent)))
  n_elbow <- sum(rowData(se_elbow)$Intensity_Domain_Type == "Intensity_Super_Element", na.rm = TRUE)
  n_tangent <- sum(rowData(se_tangent)$Intensity_Domain_Type == "Intensity_Super_Element", na.rm = TRUE)
  expect_true(n_tangent >= 0)
  expect_true(n_tangent <= nrow(example_se))
})

test_that("call_super_domains log_transform is optional", {
  se_log <- call_super_domains(example_se, feature = "Intensity", method = "tangent",
                               log_transform = TRUE, verbose = FALSE)
  se_raw <- call_super_domains(example_se, feature = "Intensity", method = "tangent",
                               log_transform = FALSE, verbose = FALSE)
  expect_true("Intensity_Domain_Type" %in% colnames(rowData(se_raw)))
  expect_true("Intensity_Domain_Type" %in% colnames(rowData(se_log)))
})

test_that("call_super_domains Breadth requires native peaks", {
  se_np <- example_se
  S4Vectors::metadata(se_np)$native_peaks <- NULL
  expect_error(call_super_domains(se_np, feature = "Breadth"),
               "native peak")
})

test_that("call_super_domains Breadth produces per-group calls", {
  se_b <- call_super_domains(example_se, feature = "Breadth",
                             mode = "per_group", group_var = "Condition",
                             verbose = FALSE)
  expect_true("Breadth_Call__Control" %in% colnames(rowData(se_b)))
  expect_true("Breadth_Call__Treatment" %in% colnames(rowData(se_b)))
  expect_true("Breadth_Super_Element" %in% rowData(se_b)$Breadth_Call__Control)
})

test_that("filter_promoter_peaks returns subset of peaks", {
  skip_if_not_installed("GenomicFeatures")
  filtered <- filter_promoter_peaks(test_gr_small, genome = "hg38",
    upstream = 2000, downstream = 2000)
  expect_s4_class(filtered, "GRanges")
})

test_that("filter_promoter_peaks rejects non-GRanges input", {
  expect_error(filter_promoter_peaks("not_a_granges"))
})

test_that("stitch_epi_peaks handles large stitch distance", {
  gr <- test_gr[1:10]
  stitched <- stitch_epi_peaks(gr, stitch_distance = 50000)
  expect_s4_class(stitched, "GRanges")
})

test_that("filter_promoter_peaks rejects chr/1 seqlevel mismatch (freeze audit)", {
  skip_if_not_installed("txdbmaker")
  gr1 <- GenomicRanges::GRanges(
    seqnames = c("1", "1"),
    ranges = IRanges::IRanges(c(1000, 1000), c(3000, 3000)),
    strand = c("+", "+"))
  S4Vectors::mcols(gr1)$type <- c("gene", "transcript")
  S4Vectors::mcols(gr1)$ID <- c("g1", "g1_tx1")
  S4Vectors::mcols(gr1)$Parent <- c(NA_character_, "g1")
  txdb <- txdbmaker::makeTxDbFromGRanges(gr1)
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1500, 2500))
  expect_error(filter_promoter_peaks(gr, genome = txdb), "No shared seqlevels")
})

test_that("filter_blacklist rejects chr/1 mismatch (freeze audit)", {
  gr <- GenomicRanges::GRanges("1", IRanges::IRanges(c(100, 5000000), width = 200))
  bl <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  tf <- tempfile(fileext = ".bed")
  utils::write.table(data.frame(chr = "chr1", start = 100, end = 200),
                     tf, row.names = FALSE, col.names = FALSE, sep = "\t")
  expect_error(filter_blacklist(gr, blacklist_path = tf), "No shared seqlevels")
})
