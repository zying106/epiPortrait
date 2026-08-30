
# Tests for the heterochromatin / broad-repressive mark compatibility layer
# (plan 2026-08-10): mark-aware presets, display terminology, and the
# continuous width expansion/contraction descriptor.

data(example_se)

test_that("broad repressive marks default to no additional stitching", {
  p27 <- get_mark_preset("H3K27me3")
  p9  <- get_mark_preset("H3K9me3")
  expect_equal(p27$stitch_distance, 0L)
  expect_equal(p9$stitch_distance, 0L)
  expect_equal(p27$mark_class, "broad_repressive")
  expect_equal(p27$taxonomy_style, "repressive_remodeling")
})

test_that("active marks keep super-domain terminology and stitching", {
  pa <- get_mark_preset("H3K27ac")
  expect_equal(pa$mark_class, "active")
  expect_equal(pa$taxonomy_style, "super_domain")
  expect_equal(pa$stitch_distance, 12500L)
  expect_true("Intensity-Super" %in% pa$class_labels)
})

test_that("mark-aware display alias maps classes without changing internal values", {
  # internal classes in rowData stay canonical
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  internal <- unique(rowData(se)$Combined_Class__Control)
  expect_true(all(internal %in% c("Intensity-Super", "Breadth-Super", "Dual-Super",
                                  "Typical", "Uncertain")))

  # display mapping translates only the surface labels
  mapped <- epiPortrait:::.epi_class_display(
    c("Intensity-Super", "Breadth-Super", "Dual-Super", "Typical"),
    get_mark_preset("H3K27me3"))
  expect_equal(mapped, c("Intensity-Extreme", "Extended-Domain", "Dual-Extreme", "Typical"))
})

test_that("mark-aware landscape uses repressive display levels", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  p <- plot_domain_landscape(se, group = "Control", mark = "H3K27me3")
  expect_s3_class(p, "ggplot")
  expect_true("Extended-Domain" %in% levels(p$data$Class))
  expect_false("Breadth-Super" %in% levels(p$data$Class))
})

test_that("compute_width_transition reports continuous expansion/contraction", {
  se <- compute_width_transition(example_se, ref_group = "Control",
                                 target_group = "Treatment")
  # pair-specific column names (multiple comparisons do not overwrite)
  expect_true(all(c("WidthDelta_bp__Control_vs_Treatment",
                    "log2WidthRatio__Control_vs_Treatment",
                    "WidthDirection__Control_vs_Treatment") %in%
                    colnames(rowData(se))))
  dir_col <- rowData(se)$WidthDirection__Control_vs_Treatment
  expect_true(all(dir_col[!is.na(dir_col)] %in%
                    c("Expansion", "Contraction", "Stable")))
  # expansion must have positive delta and positive log2 ratio
  exp_idx <- which(dir_col == "Expansion")
  if (length(exp_idx) > 0) {
    expect_true(all(rowData(se)$WidthDelta_bp__Control_vs_Treatment[exp_idx] > 0))
    expect_true(all(rowData(se)$log2WidthRatio__Control_vs_Treatment[exp_idx] > 0))
  }
  # provenance stored
  expect_false(is.null(S4Vectors::metadata(se)$width_transitions))
})

test_that("compute_width_transition NULL threshold emits continuous only", {
  se <- compute_width_transition(example_se, ref_group = "Control",
                                 target_group = "Treatment",
                                 effect_threshold = NULL)
  expect_true("WidthDelta_bp__Control_vs_Treatment" %in% colnames(rowData(se)))
  expect_true("log2WidthRatio__Control_vs_Treatment" %in% colnames(rowData(se)))
  expect_false(any(grepl("^WidthDirection__", colnames(rowData(se)))))
})

test_that("compute_width_transition multiple comparisons do not overwrite", {
  # fake a third group
  se3 <- example_se
  colData(se3)$Condition <- c(rep("Control", 3), rep("Treatment", 3))
  se3 <- compute_width_transition(se3, ref_group = "Control", target_group = "Treatment")
  se3 <- compute_width_transition(se3, ref_group = "Control", target_group = "Control")
  expect_true("WidthDelta_bp__Control_vs_Treatment" %in% colnames(rowData(se3)))
  expect_true("WidthDelta_bp__Control_vs_Control" %in% colnames(rowData(se3)))
})

test_that("plot_replicate_support works for Breadth (rowData-based)", {
  se <- call_super_domains(example_se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  p <- plot_replicate_support(se, feature = "Breadth", group = "Control")
  expect_s3_class(p, "ggplot")
})

test_that("compute_width_transition requires native peaks", {
  se_np <- example_se
  SummarizedExperiment::assay(se_np, "NativeMaxPeakWidth") <- NULL
  expect_error(compute_width_transition(se_np, ref_group = "Control",
                                        target_group = "Treatment"),
               "NativeMaxPeakWidth")
})

test_that("unknown mark falls back to generic preset with warning", {
  expect_warning(p <- get_mark_preset("H4K20me3"), "Unknown mark")
  expect_equal(p$mark_class, "generic")
  expect_equal(p$taxonomy_style, "generic")
})

test_that("domain provenance is recorded when sample_sheet carries caller info", {
  # rtracklayer local BigWig I/O fails on Windows ("UCSC library operation failed");
  # rtracklayer#52/#62/#128/#151. Narrow skip: this block only. epiPortrait logic
  # that does not read BigWigs is still executed on Windows.
  skip_on_os("windows")
  # build with a tiny BigWig-only sheet carrying domain_source metadata
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")
  samples <- data.frame(
    SampleID = "S1", Condition = "C",
    bw_path = file.path(extdata, "C1.bw"),
    domain_source = "epic2", domain_caller = "epic2", caller_version = "0.0.1")
  se <- build_portrait_matrix(
    samples,
    consensus_peaks = rtracklayer::import(file.path(extdata, "peaks.bed")),
    workers = 1)
  dp <- S4Vectors::metadata(se)$domain_provenance
  expect_equal(dp$domain_source, "epic2")
  expect_equal(dp$domain_caller, "epic2")
})
