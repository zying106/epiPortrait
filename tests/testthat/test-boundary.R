# Boundary / defensive-branch tests for real logic paths that were
# uncovered (not added just to inflate coverage). Each case exercises an
# error or fallback branch that a real user could hit.

data(example_se)

# ---- filter_blacklist --------------------------------------------------------
test_that("filter_blacklist handles empty input", {
  gr <- GenomicRanges::GRanges()
  expect_s4_class(filter_blacklist(gr, genome = "hg38"), "GRanges")
})
test_that("filter_blacklist errors on missing custom path", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  expect_error(filter_blacklist(gr, genome = "hg38",
                                blacklist_path = "does_not_exist.bed"),
               "File not found")
})

test_that("filter_blacklist errors on genome without built-in blacklist", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  expect_error(filter_blacklist(gr, genome = "not_a_genome"),
               "No built-in blacklist")
})

# ---- plot_portrait_pca boundary ---------------------------------------------
test_that("plot_portrait_pca errors on a single sample", {
  se1 <- example_se[, 1]
  expect_error(plot_portrait_pca(se1, feature = "Intensity",
                                 group_var = "Condition"),
               "PCA requires at least 2 samples")
})

# ---- check_signal_compatibility ---------------------------------------------
test_that("check_signal_compatibility rejects non-GRanges regions", {
  ss <- data.frame(SampleID = "S1", Condition = "C",
                   bw_path = system.file("extdata", "C1.bw", package = "epiPortrait"))
  expect_error(check_signal_compatibility(ss, regions = "not_granges"),
               "GRanges")
})

test_that(".genome_tiles errors on unreadable BigWig", {
  expect_error(epiPortrait:::.genome_tiles("no_such_file.bw"),
               "Could not read BigWig header")
})

# ---- validate_epiportrait_object defensive branches -------------------------
test_that("validate_epiportrait_object rejects non-SE input", {
  expect_error(validate_epiportrait_object(data.frame(x = 1:3)), "SummarizedExperiment")
})

test_that("validate_epiportrait_object rejects non-unique rownames", {
  se <- example_se
  rownames(se)[1] <- rownames(se)[2]
  expect_error(validate_epiportrait_object(se), "unique")
})

# ---- get_call_provenance / get_call_results fallbacks -----------------------
test_that("get_call_provenance returns NULL when no calls stored", {
  expect_null(get_call_provenance(example_se, "Intensity"))
})

test_that("get_call_results long format returns Group/Call columns", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  l <- get_call_results(se, "Intensity", long = TRUE)
  expect_true(all(c("Domain_ID", "Group", "Call") %in% colnames(l)))
  expect_true("Control" %in% unique(l$Group))
})

# ---- export conditional branches --------------------------------------------
test_that("export_epiportrait_results handles no combined/transition columns", {
  se <- call_super_domains(example_se, feature = "Intensity", verbose = FALSE)
  out <- export_epiportrait_results(se, outdir = tempfile("epi_plain"),
                                    group_var = NULL)
  # no combined_classes.tsv or transition_results.tsv because none were computed
  expect_true(file.exists(file.path(out, "domain_results.tsv")))
  expect_false(file.exists(file.path(out, "combined_classes.tsv")))
  expect_false(file.exists(file.path(out, "transition_results.tsv")))
})
