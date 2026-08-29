library(testthat)
library(epiPortrait)
library(SummarizedExperiment)
library(GenomicRanges)
library(S4Vectors)

make_diff_se <- function(n_dom = 200, n_rep = 3, effect = 3, lfc_cut = 0.5) {
  set.seed(1)
  gr <- GRanges(rep("chr1", n_dom), IRanges(seq_len(n_dom) * 1000, width = 500))
  G <- factor(rep(c("Tum", "Nor"), each = n_rep))
  M <- matrix(rnorm(n_dom * (2 * n_rep), mean = 3, sd = 1), nrow = n_dom)
  # first 20 = true Gain, next 20 = true Loss
  M[1:20, G == "Tum"] <- M[1:20, G == "Tum"] + effect
  M[21:40, G == "Tum"] <- M[21:40, G == "Tum"] - effect
  M <- 2^M
  colnames(M) <- paste0("s", seq_len(2 * n_rep))
  SummarizedExperiment(
    assays = list(Intensity = M), rowRanges = gr,
    colData = DataFrame(SampleID = colnames(M), Condition = as.character(G)))
}

test_that("analyze_differential_domains detects true gains", {
  se <- make_diff_se()
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum",
                                      fdr_cutoff = 0.05, logFC_cutoff = 0.5)
  rd <- rowData(res)
  expect_true("Intensity_DiffStatus" %in% colnames(rd))
  expect_identical(sum(rd$Intensity_DiffStatus[1:20] == "Gain", na.rm = TRUE), 20L)
  expect_gte(sum(rd$Intensity_DiffStatus[21:40] == "Loss", na.rm = TRUE), 18L)
  expect_true(all(c("Intensity_logFC", "Intensity_adj.P.Val", "Intensity_P.Value",
                    "Intensity_t", "Intensity_AveExpr") %in% colnames(rd)))
})

test_that("analyze_differential_domains errors on <2 replicates per group", {
  se <- make_diff_se(n_rep = 1)
  expect_error(
    analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum"),
    "needs >=2 replicates")
})

test_that("analyze_differential_domains stored provenance", {
  se <- make_diff_se()
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum")
  prov <- S4Vectors::metadata(res)$differential_domains
  expect_equal(prov$feature, "Intensity")
  expect_true(prov$n_tested == 200)
  expect_true(prov$n_gain >= 20L)
  expect_true(prov$trend)      # mean-variance trend ON by default (200 domains)
  expect_false(prov$robust)
})

test_that("analyze_differential_domains eBayes options and fallbacks", {
  se <- make_diff_se()

  # trend = FALSE is honored and recorded.
  res <- analyze_differential_domains(se, ref_group = "Nor",
                                      target_group = "Tum", trend = FALSE)
  prov <- S4Vectors::metadata(res)$differential_domains
  expect_false(prov$trend)

  # robust = TRUE runs and is recorded.
  res <- analyze_differential_domains(se, ref_group = "Nor",
                                      target_group = "Tum", robust = TRUE)
  expect_true(S4Vectors::metadata(res)$differential_domains$robust)

  # Automatic trend fallback when few domains pass filtering (<20).
  set.seed(2)
  gr_small <- GenomicRanges::GRanges(rep("chr1", 8),
                                     IRanges::IRanges(seq_len(8) * 1000, width = 500))
  M_small <- 2^matrix(rnorm(8 * 4, mean = 3, sd = 1), nrow = 8)
  colnames(M_small) <- paste0("s", seq_len(4))
  se_small <- SummarizedExperiment(
    assays = list(Intensity = M_small), rowRanges = gr_small,
    colData = S4Vectors::DataFrame(
      SampleID = colnames(M_small),
      Condition = rep(c("Nor", "Tum"), each = 2)))
  expect_warning(
    res_small <- analyze_differential_domains(se_small, ref_group = "Nor",
                                              target_group = "Tum"),
    "mean-variance trend")
  expect_false(S4Vectors::metadata(res_small)$differential_domains$trend)

  # Invalid eBayes option values are rejected.
  expect_error(analyze_differential_domains(se, ref_group = "Nor",
                                            target_group = "Tum", trend = "yes"),
               "trend must be TRUE or FALSE")
  expect_error(analyze_differential_domains(se, ref_group = "Nor",
                                            target_group = "Tum", robust = NA),
               "robust must be TRUE or FALSE")
})

test_that("two-coefficient character contrast works and matches ref/target path", {
  se <- make_diff_se()

  res_ref <- analyze_differential_domains(se, ref_group = "Nor",
                                          target_group = "Tum", trend = FALSE)
  res_pair <- analyze_differential_domains(se, design = ~ 0 + Condition,
                                           contrast = c("ConditionTum",
                                                        "ConditionNor"),
                                           trend = FALSE)
  expect_equal(rowData(res_pair)$Intensity_logFC,
               rowData(res_ref)$Intensity_logFC)
  expect_equal(rowData(res_pair)$Intensity_P.Value,
               rowData(res_ref)$Intensity_P.Value, tolerance = 1e-10)

  # Provenance stores the normalized numeric weights.
  prov <- S4Vectors::metadata(res_pair)$differential_domains
  expect_identical(prov$contrast[["ConditionTum"]], 1)
  expect_identical(prov$contrast[["ConditionNor"]], -1)

  # Unknown coefficient names fail loudly.
  expect_error(
    analyze_differential_domains(se, design = ~ 0 + Condition,
                                 contrast = c("ConditionTum", "Wrong")),
    "Unknown design coefficient")

  # design without contrast is rejected instead of silently uncontrasted.
  expect_error(analyze_differential_domains(se, design = ~ Condition),
               "without contrast")
})

test_that("volcano plot returns a ggplot", {
  skip_if_not_installed("ggplot2")
  se <- make_diff_se()
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum")
  p <- plot_differential_volcano(res)
  expect_s3_class(p, "ggplot")
})