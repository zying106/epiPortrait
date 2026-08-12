# Tests for the object-contract accessors/export and the academic
# visualization suite (2026-08-08 design docs).

data(example_se)

# ---- object contract: validation -------------------------------------------
test_that("validate_epiportrait_object passes on a valid object", {
  expect_true(validate_epiportrait_object(example_se))
})

test_that("validate_epiportrait_object rejects bad rownames", {
  se <- example_se
  rownames(se)[1] <- rownames(se)[2]  # duplicate
  expect_error(validate_epiportrait_object(se), "unique")
})

# ---- object contract: accessors --------------------------------------------
test_that("get_domain_results returns coordinates + IntervalWidth", {
  dr <- get_domain_results(example_se, group_var = "Condition")
  expect_true(all(c("Domain_ID", "chr", "start", "end", "IntervalWidth") %in% colnames(dr)))
  expect_equal(nrow(dr), nrow(example_se))
  expect_true("Intensity_Mean__Control" %in% colnames(dr))
})

test_that("get_sample_results returns colData", {
  sr <- get_sample_results(example_se)
  expect_true("SampleID" %in% colnames(sr))
  expect_equal(nrow(sr), ncol(example_se))
})

test_that("get_call_results works wide and long", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  w <- get_call_results(se, "Intensity")
  expect_true("Intensity_Call__Control" %in% colnames(w))
  l <- get_call_results(se, "Intensity", long = TRUE)
  expect_true(all(c("Group", "Call") %in% colnames(l)))
})

test_that("get_call_provenance returns stored metadata", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  prov <- get_call_provenance(se, "Intensity")
  expect_equal(prov$method, "tangent")
  expect_false(prov$log_transform_used)
})

test_that("export_epiportrait_results writes tables + RDS", {
  se <- call_super_domains(example_se, feature = "Intensity", verbose = FALSE)
  out <- export_epiportrait_results(se, outdir = tempfile("epi_export"),
                                    group_var = "Condition")
  expect_true(file.exists(file.path(out, "domain_results.tsv")))
  expect_true(file.exists(file.path(out, "epiPortrait_object.rds")))
  expect_true(file.exists(file.path(out, "object_manifest.txt")))
  expect_true(file.exists(file.path(out, "assays", "Intensity.tsv")))
})

# ---- academic visualization -------------------------------------------------
test_that("plot_domain_landscape returns ggplot with axes", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  p <- plot_domain_landscape(se, group = "Control")
  expect_s3_class(p, "ggplot")
  expect_true(all(c("Intensity", "Breadth") %in% names(p$data)))
  expect_equal(nrow(p$data), nrow(se))
})

test_that("plot_domain_landscape per_class labels only super classes", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  p <- plot_domain_landscape(se, group = "Control", label_mode = "per_class",
                             label_n = 5)
  repel <- p$layers[vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1))]
  expect_length(repel, 1)
  labelled <- unique(as.character(repel[[1]]$data$Class))
  super <- epiPortrait:::.epi_super_class_labels()
  expect_true(all(labelled %in% super))
  expect_true(any(super %in% labelled))
  # legacy top-intensity mode still works
  p2 <- plot_domain_landscape(se, group = "Control", label_mode = "top_intensity",
                              label_n = 8)
  expect_s3_class(p2, "ggplot")
})

test_that("plot_domain_class_composition super focus drops Typical/Uncertain", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  p <- plot_domain_class_composition(se, focus = "super")
  classes <- unique(as.character(p$data$Class))
  super <- epiPortrait:::.epi_super_class_labels()
  expect_true(all(classes %in% super))
  expect_false("Typical" %in% classes)
  # full composition incl. Typical remains available
  p2 <- plot_domain_class_composition(se, focus = "all")
  expect_true("Typical" %in% unique(as.character(p2$data$Class)))
})

test_that("plot_replicate_support default min_support drops 0/n bin", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  p <- plot_replicate_support(se, feature = "Intensity", group = "Control")
  expect_false("0/3" %in% levels(p$data$Level))
  p2 <- plot_replicate_support(se, feature = "Intensity", group = "Control",
                               min_support = 0)
  expect_true("0/3" %in% levels(p2$data$Level))
})

test_that("plot_domain_feature_profile annotates class under axis", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  p <- plot_domain_feature_profile(se, peak_id = rownames(se)[1])
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "class label")
  # falls back gracefully without combined classes
  p2 <- plot_domain_feature_profile(example_se, peak_id = rownames(example_se)[1])
  expect_s3_class(p2, "ggplot")
})

test_that("plot_class_transition cell counts sum to domains", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  p <- plot_class_transition(se, ref_group = "Control", target_group = "Treatment")
  expect_s3_class(p, "ggplot")
})

test_that("plot_replicate_support returns ggplot", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  p <- plot_replicate_support(se, feature = "Intensity", group = "Control")
  expect_s3_class(p, "ggplot")
})

test_that("plot_portrait_pca supports combined features", {
  p <- plot_portrait_pca(example_se, feature = "Intensity",
                         features = c("Intensity", "SignalDispersion"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_portrait_correlation defaults to spearman", {
  p <- plot_portrait_correlation(example_se, feature = "SignalDispersion")
  expect_s3_class(p, "ggplot")
})

# ---- upgraded hockey / track -------------------------------------------------
test_that("plot_hockey_stick supports per-group mode with cutoff band", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE,
                           n_bootstrap = 10, verbose = FALSE)
  p <- plot_hockey_stick(se, feature = "Intensity", group = "Control",
                         label_genes = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_hockey_stick consensus shows Uncertain class", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  p <- plot_hockey_stick(se, feature = "Intensity", label_genes = FALSE)
  expect_s3_class(p, "ggplot")
  expect_true("Uncertain" %in% levels(p$data$Type_Label))
})

test_that("plot_peak_track show_native_occupancy overlays span", {
  skip_on_os("windows")
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")
  samples <- data.frame(SampleID = c("C1", "C2"), Condition = c("Control", "Control"),
                        bw_path = file.path(extdata, c("C1.bw", "C2.bw")),
                        peak_path = file.path(extdata, c("C1_peaks.bed", "C2_peaks.bed")))
  se <- build_portrait_matrix(samples,
                              consensus_peaks = rtracklayer::import(file.path(extdata, "peaks.bed")),
                              workers = 1)
  p <- plot_peak_track(se, peak_id = rownames(se)[1], show_native_occupancy = TRUE)
  expect_s3_class(p, "ggplot")
})


# ---- get_replicate_calls (design 2026-08-10) --------------------------------
test_that("get_replicate_calls returns per-domain x per-replicate matrix", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  m <- get_replicate_calls(se, feature = "Intensity", group = "Control")
  expect_equal(nrow(m), nrow(example_se))
  expect_true(all(grepl("^Control_", colnames(m))))
  expect_true(all(m[, 1] %in% c("Intensity_Super_Element", "Intensity_Typical", NA)))
  # long format
  l <- get_replicate_calls(se, feature = "Intensity", group = "Control", long = TRUE)
  expect_true(all(c("domain_id", "SampleID", "Group", "call") %in% colnames(l)))
  expect_equal(nrow(l), nrow(example_se) * ncol(m))
})

# ---- get_replicate_calls: freeze audit coverage (2026-08-10) ---------------
test_that("get_replicate_calls works for Breadth/per_group", {
  se <- call_super_domains(example_se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition",
                           verbose = FALSE)
  m <- get_replicate_calls(se, feature = "Breadth", group = "Control")
  expect_equal(nrow(m), nrow(example_se))
  expect_true(all(grepl("^Control_", colnames(m))))
  expect_true(all(m[, 1] %in% c("Breadth_Super_Element", "Breadth_Typical", NA)))
})

test_that("get_replicate_calls per_sample returns a matrix", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_sample", verbose = FALSE)
  m <- get_replicate_calls(se, feature = "Intensity")
  expect_true(is.matrix(m))
  expect_equal(nrow(m), nrow(example_se))
  expect_equal(ncol(m), ncol(example_se))
})

test_that("get_replicate_calls respects custom group_var in long output", {
  se <- example_se
  colData(se)$Arm <- c("A", "A", "A", "B", "B", "B")
  se <- call_super_domains(se, feature = "Intensity",
                           mode = "per_group", group_var = "Arm",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  l <- get_replicate_calls(se, feature = "Intensity", group = "A", long = TRUE)
  expect_true("Group" %in% colnames(l))
  expect_true(all(l$Group == "A"))
})

test_that("get_replicate_calls invalid group gives a clear error", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  expect_error(get_replicate_calls(se, feature = "Intensity", group = "XYZ"),
               "not found")
})

# ---- freeze audit: landscape custom group_var (2026-08-10) ------------------
test_that("plot_domain_landscape respects custom group_var", {
  se <- example_se
  colData(se)$Arm <- c("A", "A", "A", "B", "B", "B")
  se <- call_super_domains(se, feature = "Intensity",
                           mode = "per_group", group_var = "Arm",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Arm", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Arm")
  p <- plot_domain_landscape(se, group = "A", group_var = "Arm")
  expect_s3_class(p, "ggplot")
  expect_true(nrow(p$data) == nrow(se))
})

# ---- get_uncertain_cause (freeze 2026-08-11) --------------------------------
test_that("get_uncertain_cause classifies Uncertain origins", {
  se <- call_super_domains(example_se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  uc <- get_uncertain_cause(se, feature = "Breadth", group = "Control")
  expect_true(all(c("Domain_ID", "Group_Call", "N_Valid_Replicates",
                    "Min_Valid_Replicates", "Cause") %in% colnames(uc)))
  expect_equal(nrow(uc), nrow(se))
  # example_se has 6 no-call domains -> insufficient_valid_replicates
  expect_equal(sum(uc$Cause == "insufficient_valid_replicates", na.rm = TRUE), 6)
  # domains with a real group call have NA cause; uncertain ones have a cause
  called <- !is.na(uc$Group_Call)
  expect_true(all(is.na(uc$Cause[called])))
  expect_true(all(!is.na(uc$Cause[!called])))
})

test_that("plot_uncertain_cause returns a ggplot", {
  se <- call_super_domains(example_se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  p <- plot_uncertain_cause(se, feature = "Breadth", group = "Control")
  expect_s3_class(p, "ggplot")
})

test_that("get_domain_results includes uncertain-cause columns", {
  se <- call_super_domains(example_se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  dr <- get_domain_results(se, group_var = "Condition")
  expect_true("Breadth_Uncertain_Cause__Control" %in% colnames(dr))
})
