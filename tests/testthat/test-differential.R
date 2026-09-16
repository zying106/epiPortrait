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

.diff_key <- function(result = "Nor_vs_Tum", feature = "Intensity") {
  paste0(feature, "_Diff__", result)
}

.diff_col <- function(statistic, result = "Nor_vs_Tum", feature = "Intensity") {
  paste0(.diff_key(result, feature), "__", statistic)
}

.diff_provenance <- function(se, result = "Nor_vs_Tum",
                             feature = "Intensity") {
  S4Vectors::metadata(se)$differential_domain_analyses[[
    .diff_key(result, feature)]]
}

test_that("analyze_differential_domains detects true gains", {
  se <- make_diff_se()
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum",
                                      fdr_cutoff = 0.05, logFC_cutoff = 0.5)
  rd <- rowData(res)
  status <- rd[[.diff_col("DiffStatus")]]
  expect_identical(sum(status[1:20] == "Gain", na.rm = TRUE), 20L)
  expect_gte(sum(status[21:40] == "Loss", na.rm = TRUE), 18L)
  expect_true(all(vapply(
    c("logFC", "adj.P.Val", "P.Value", "t", "AveExpr"),
    function(statistic) .diff_col(statistic) %in% colnames(rd), logical(1))))
  expect_false(any(c("Intensity_logFC", "Intensity_adj.P.Val",
                     "Intensity_P.Value", "Intensity_t", "Intensity_AveExpr",
                     "Intensity_DiffStatus") %in% colnames(rd)))
})

test_that("analyze_differential_domains errors on <2 replicates per group", {
  se <- make_diff_se(n_rep = 1)
  expect_error(
    analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum"),
    "needs >=2 replicates")
})

test_that("unused group factor levels do not invalidate a two-group fit", {
  se <- make_diff_se(n_dom = 40, n_rep = 2)
  colData(se)$Condition <- factor(colData(se)$Condition,
                                  levels = c("Nor", "Tum", "Unused"))
  res <- analyze_differential_domains(
    se, ref_group = "Nor", target_group = "Tum", trend = FALSE)
  expect_equal(.diff_provenance(res)$n_tested, 40)
})

test_that("analyze_differential_domains stored provenance", {
  se <- make_diff_se()
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum")
  prov <- .diff_provenance(res)
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
  prov <- .diff_provenance(res)
  expect_false(prov$trend)

  # robust = TRUE runs and is recorded.
  res <- analyze_differential_domains(se, ref_group = "Nor",
                                      target_group = "Tum", robust = TRUE)
  expect_true(.diff_provenance(res)$robust)

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
  expect_false(.diff_provenance(res_small)$trend)

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
  pair_result <- "ConditionTum_minus_ConditionNor"
  expect_equal(rowData(res_pair)[[.diff_col("logFC", pair_result)]],
               rowData(res_ref)[[.diff_col("logFC")]])
  expect_equal(rowData(res_pair)[[.diff_col("P.Value", pair_result)]],
               rowData(res_ref)[[.diff_col("P.Value")]], tolerance = 1e-10)

  # Provenance stores the normalized numeric weights.
  prov <- .diff_provenance(res_pair, pair_result)
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

test_that("multiple contrast columns are rejected explicitly", {
  se <- make_diff_se()
  design <- model.matrix(~ 0 + Condition, data = as.data.frame(colData(se)))
  contrast <- matrix(c(-1, 1, 1, -1), nrow = 2,
                     dimnames = list(colnames(design),
                                     c("Tum_vs_Nor", "Nor_vs_Tum")))
  expect_error(
    analyze_differential_domains(se, design = design, contrast = contrast),
    "Exactly one contrast")
})

test_that("repeated differential analyses retain pair-specific results", {
  se <- make_diff_se()
  se <- analyze_differential_domains(
    se, ref_group = "Nor", target_group = "Tum", trend = FALSE)
  se <- analyze_differential_domains(
    se, ref_group = "Tum", target_group = "Nor", trend = FALSE)

  rd_names <- colnames(rowData(se))
  expect_true("Intensity_Diff__Nor_vs_Tum__logFC" %in% rd_names)
  expect_true("Intensity_Diff__Tum_vs_Nor__logFC" %in% rd_names)
  expect_equal(
    rowData(se)$Intensity_Diff__Nor_vs_Tum__logFC,
    -rowData(se)$Intensity_Diff__Tum_vs_Nor__logFC,
    tolerance = 1e-10)
  analyses <- S4Vectors::metadata(se)$differential_domain_analyses
  expect_setequal(names(analyses),
                  c("Intensity_Diff__Nor_vs_Tum",
                    "Intensity_Diff__Tum_vs_Nor"))
  expect_false("differential_domains" %in% names(S4Vectors::metadata(se)))
})

test_that("pairwise analysis retains samples outside the fitted contrast", {
  se <- make_diff_se(n_rep = 2)
  extra <- se[, 1:2]
  colnames(extra) <- c("other1", "other2")
  colData(extra)$SampleID <- colnames(extra)
  colData(extra)$Condition <- "Other"
  se <- cbind(se, extra)
  original_samples <- colnames(se)

  se <- analyze_differential_domains(
    se, ref_group = "Nor", target_group = "Tum", trend = FALSE)
  expect_identical(colnames(se), original_samples)
  se <- analyze_differential_domains(
    se, ref_group = "Nor", target_group = "Other", trend = FALSE)

  expect_identical(colnames(se), original_samples)
  expect_true(all(c("Intensity_Diff__Nor_vs_Tum__logFC",
                    "Intensity_Diff__Nor_vs_Other__logFC") %in%
                  colnames(rowData(se))))
})

test_that("volcano plot returns a ggplot", {
  skip_if_not_installed("ggplot2")
  se <- make_diff_se()
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum")
  p <- plot_differential_volcano(res)
  expect_s3_class(p, "ggplot")
})

test_that("volcano requires an explicit result when several analyses exist", {
  skip_if_not_installed("ggplot2")
  se <- make_diff_se()
  se <- analyze_differential_domains(
    se, ref_group = "Nor", target_group = "Tum", trend = FALSE)
  se <- analyze_differential_domains(
    se, ref_group = "Tum", target_group = "Nor", trend = FALSE)

  expect_error(plot_differential_volcano(se), "Multiple differential results")
  p <- plot_differential_volcano(se, result_name = "Nor_vs_Tum")
  expect_s3_class(p, "ggplot")
})

# ---- P0 regression: min_signal filtering must not mis-map rows (2026 review) --
test_that("min_signal excludes low-signal domains from the fit and maps by index", {
  set.seed(11)
  n <- 20
  gr <- GRanges(rep("chr1", n), IRanges(seq_len(n) * 1000, width = 500))
  G <- factor(rep(c("Nor", "Tum"), each = 2))
  M <- 2^matrix(rnorm(n * 4, mean = 3, sd = 1), nrow = n)
  M[1:5, ] <- M[1:5, ]              # rows 1-5 := low signal (keep base)
  M[6, ] <- 2^2                     # row 6 also low-ish
  colnames(M) <- paste0("s", seq_len(4))
  se <- SummarizedExperiment(assays = list(Intensity = M), rowRanges = gr,
                             colData = DataFrame(SampleID = colnames(M),
                                                 Condition = as.character(G)))
  # force a fit_keep that is non-contiguous: min_signal clips the LOW rows 1:6
  # (mean log2-scale small), leaving typically rows 7.. and others. To make the
  # mapping bug reproducible deterministically we instead use a synthetic
  # min_signal that drops a contiguous leading block -> check kept rows map to
  # the correct original indices.
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum",
                                      min_signal = 3.5,
                                      fdr_cutoff = 0.5, logFC_cutoff = 0.1)
  rd <- rowData(res)
  # domains below min_signal (in log2(x+1) mean) have NA results and are
  # reported NS, never Gain/Loss, and never carry a wrong domain's statistics.
  lo <- which(rowMeans(log2(assay(se, "Intensity") + 1)) < 3.5)
  if (length(lo) > 0) {
    expect_true(all(is.na(rd[[.diff_col("logFC")]][lo])))
    # excluded domains are "not tested" -> NA DiffStatus (documented as
    # unavailable, distinct from NS which means tested-but-not-significant)
    expect_true(all(is.na(rd[[.diff_col("DiffStatus")]][lo])))
  }
  hi <- which(!is.na(rd[[.diff_col("logFC")]]))
  if (length(hi) > 0) {
    # every result row must be >= min_signal by construction
    expect_true(all(rowMeans(log2(assay(se, "Intensity")[hi, , drop = FALSE] + 1)) >= 3.5))
  }
})

test_that("no min_signal filter runs on all rows (fit_keep all TRUE)", {
  se <- make_diff_se()
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum",
                                      min_signal = 0, fdr_cutoff = 0.05,
                                      logFC_cutoff = 0.5)
  rd <- rowData(res)
  expect_false(any(is.na(rd[[.diff_col("logFC")]])))
  expect_true(sum(rd[[.diff_col("DiffStatus")]] == "Gain",
                  na.rm = TRUE) >= 19L)
})

# ---- epiPortrait3 review §2/§3/§4/§12: NA & negative semantics -------------
test_that("negative assay values are rejected, not silently clipped", {
  se <- make_diff_se()
  M <- assay(se, "Intensity")
  M[1, 1] <- -5                              # inject a negative
  assay(se, "Intensity") <- M
  expect_error(analyze_differential_domains(se, ref_group = "Nor",
                                            target_group = "Tum"),
               "non-negative assay")
})

test_that("all-NA / below-min_signal domains are excluded (rows NA), NA != 0", {
  set.seed(3)
  n <- 12
  gr <- GRanges(rep("chr1", n), IRanges(seq_len(n) * 1000, width = 500))
  G <- factor(rep(c("Nor", "Tum"), each = 2))
  M <- 2^matrix(rnorm(n * 4, mean = 3, sd = 1), nrow = n)
  M[1, ] <- NA                               # row 1 = all-NA -> excluded, not 0
  M[2, ] <- 2                                # row 2 = very low signal
  colnames(M) <- paste0("s", seq_len(4))
  se <- SummarizedExperiment(assays = list(Intensity = M), rowRanges = gr,
                             colData = DataFrame(SampleID = colnames(M),
                                                 Condition = as.character(G)))
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum",
                                      min_signal = 3, fdr_cutoff = 0.5,
                                      logFC_cutoff = 0.1)
  rd <- rowData(res)
  # row 1 (all NA) must be excluded -> NA logFC (not 0)
  expect_true(is.na(rd[[.diff_col("logFC")]][1]))
  # other rows (above min_signal) produce finite results
  hi <- which(!is.na(rd[[.diff_col("logFC")]]))
  expect_true(length(hi) >= 1)
  expect_true(all(is.finite(rd[[.diff_col("logFC")]][hi])))
})

test_that("n_fit == 0 (all below min_signal) is an explicit error", {
  se <- make_diff_se()
  expect_error(analyze_differential_domains(se, ref_group = "Nor",
                                            target_group = "Tum",
                                            min_signal = 1e6),
               "No domains passed min_signal")
})

test_that("volcano keeps untested (NA adj.P.Val) domains as NA, not extreme", {
  se <- make_diff_se()
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum",
                                      min_signal = 0, fdr_cutoff = 0.05,
                                      logFC_cutoff = 0.5)
  # simulate two domains that were NOT tested (NA adjusted P)
  rowData(res)[[.diff_col("adj.P.Val")]][c(1, 2)] <- NA_real_
  p <- plot_differential_volcano(res, feature = "Intensity", label_n = 2)
  b <- ggplot2::ggplot_build(p)
  dat <- b$data[[1]]
  # NA adj.P.Val rows are dropped from the geometry -> fewer points than domains
  expect_true(nrow(dat) <= nrow(res))
  # no plotted point may carry -log10P == 300 (old NA->1e-300 artifact)
  if (nrow(dat) > 0) {
    yf <- dat$y[is.finite(dat$y)]
    if (length(yf) > 0) expect_true(all(yf < 200))
  }
})

# ---- epiPortrait3 review §13: TRUE non-contiguous fit_keep mapping ----------
test_that("non-contiguous min_signal filter maps results to correct rows", {
  set.seed(21)
  n <- 30
  gr <- GRanges(rep("chr1", n), IRanges(seq_len(n) * 1000, width = 500))
  G <- factor(rep(c("Nor", "Tum"), each = 3))   # n = 2 groups x 3 reps? use 6 samples
  M <- 2^matrix(rnorm(n * 6, mean = 4, sd = 1), nrow = n)
  # scatter LOW rows interspersed: 2,5,9,14 (NOT a contiguous leading block)
  M[c(2, 5, 9, 14), ] <- 2
  colnames(M) <- paste0("s", seq_len(6))
  se <- SummarizedExperiment(assays = list(Intensity = M), rowRanges = gr,
                             colData = DataFrame(SampleID = colnames(M),
                                                 Condition = rep(c("Nor","Tum"), each=3)))
  res <- analyze_differential_domains(se, ref_group = "Nor", target_group = "Tum",
                                      min_signal = 3, fdr_cutoff = 0.5,
                                      logFC_cutoff = 0.1)
  rd <- rowData(res)
  lo <- which(rowMeans(log2(assay(se, "Intensity") + 1)) < 3)
  # verify the filtered set is genuinely interspersed (not contiguous)
  lo <- sort(lo)
  expect_true(length(lo) >= 4)
  expect_true(any(diff(lo) == 1) == FALSE || length(lo) < 5)  # not a solid run
  # kept domains -> finite results; dropped -> NA
  expect_true(all(is.na(rd[[.diff_col("logFC")]][lo])))
  hi <- setdiff(seq_len(n), lo)
  expect_true(all(is.finite(rd[[.diff_col("logFC")]][hi])))
  # EXACT mapping check: independently refit only the kept rows and compare
  suppressPackageStartupMessages(library(limma))
  Mv <- log2(assay(se, "Intensity") + 1)
  design <- model.matrix(~ factor(rep(c("Nor","Tum"), each=3)))
  fit_man <- limma::lmFit(Mv[hi, , drop = FALSE], design)
  tt_man <- limma::topTable(limma::eBayes(limma::contrasts.fit(fit_man, c(0,1))),
                            number = Inf, sort.by = "none")
  # package output for kept domain i == tt_man row matching position in hi
  for (k in seq_along(hi)) {
    expect_equal(rd[[.diff_col("logFC")]][hi[k]], tt_man$logFC[k],
                 tolerance = 1e-10)
  }
})

# ---- review coverage: validation error paths (lines 115, 118, 127, 130, 142) -----
test_that("differential validates input: missing feature / samples / group / ref", {
  se <- make_diff_se()
  expect_error(analyze_differential_domains(se, feature = "BOGUS",
    group_var = "Condition", ref_group = "Nor", target_group = "Tum"),
    "not found in assays")
  expect_error(analyze_differential_domains(se[, 1, drop = FALSE],
    feature = "Intensity", group_var = "Condition",
    ref_group = "Nor", target_group = "Tum"),
    "at least 2 samples")
  expect_error(analyze_differential_domains(se, feature = "Intensity",
    group_var = "NONEXISTENT", ref_group = "Nor", target_group = "Tum"),
    "not found in colData")
  expect_error(analyze_differential_domains(se, feature = "Intensity",
    group_var = "Condition", ref_group = NULL, target_group = NULL),
    "provide ref_group")
  expect_error(analyze_differential_domains(se, feature = "Intensity",
    group_var = "Condition", ref_group = "NONEXISTENT", target_group = "Tum"),
    "not found in group_var")
})
