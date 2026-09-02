# Tests for cross-mark integration (integrate_cross_mark) and per-sample
# pairwise transitions (transition_matrix_per_sample), 2026-08-26.

library(epiPortrait)
library(GenomicRanges)
library(SummarizedExperiment)
library(S4Vectors)

# two objects on chr1; object A: 4 domains (0-40kb), object B: 3 overlapping
# domains (slightly shifted universe, as happens with different stitching)
make_cross_se <- function(type) {
  suppressWarnings({
    if (type == "A") {
      dom <- GRanges("chr1", IRanges::IRanges(c(1, 101, 201, 301),
                                              c(300, 300, 300, 600)))
      se <- SummarizedExperiment(
        assays = list(Intensity = matrix(1, length(dom), 2)),
        rowRanges = dom)
      rownames(se) <- sprintf("A_%02d", seq_len(length(dom)))
      colnames(se) <- c("t0", "t1")
      colData(se)$Condition <- c("t0", "t1")
      rowData(se)$Breadth_Call__t0 <- c("Breadth_Super_Element",
                                        "Breadth_Typical",
                                        "Breadth_Super_Element",
                                        "Breadth_Typical")
      rowData(se)$Breadth_Call__t1 <- c("Breadth_Typical",
                                        "Breadth_Super_Element",
                                        "Breadth_Typical", NA)
    } else {
      # B universe shifted: only A4 (301-600) overlaps two B domains
      dom <- GRanges("chr1", IRanges::IRanges(c(1, 302, 601),
                                              c(400, 500, 900)))
      se <- SummarizedExperiment(
        assays = list(Intensity = matrix(1, length(dom), 2)),
        rowRanges = dom)
      rownames(se) <- sprintf("B_%02d", seq_len(length(dom)))
      colnames(se) <- c("t0", "t1")
      colData(se)$Condition <- c("t0", "t1")
      rowData(se)$Intensity_Call__t0 <- c("Intensity_Super_Element",
                                          "Intensity_Typical",
                                          "Intensity_Super_Element")
      rowData(se)$Intensity_Call__t1 <- c(NA,
                                          "Intensity_Typical",
                                          "Intensity_Typical")
    }
  })
  se
}

test_that("transition_matrix_per_sample compares per-sample call columns", {
  se <- make_cross_se("A")
  se <- transition_matrix_per_sample(se, feature = "Breadth",
                                     timepoints = c("t0", "t1"), verbose = FALSE)
  col <- "Breadth_SampleTransition__t0_vs_t1"
  expect_true(col %in% colnames(rowData(se)))
  tr <- rowData(se)[[col]]
  # domain 1: Super -> Typical
  expect_equal(tr[1], "Breadth_Super_Element_to_Breadth_Typical")
  # domain 2: Typical -> Super
  expect_equal(tr[2], "Breadth_Typical_to_Breadth_Super_Element")
  # domain 3: Super -> Typical
  expect_equal(tr[3], "Breadth_Super_Element_to_Breadth_Typical")
  # domain 4: Typical -> NA (Uncertain propagation)
  expect_equal(tr[4], "Uncertain")
  expect_false(is.null(metadata(se)$transitions[["t0_vs_t1"]]$counts))
  # auto timepoint discovery from column suffixes
  se2 <- make_cross_se("A")
  se2 <- transition_matrix_per_sample(se2, feature = "Breadth",
                                      verbose = FALSE)
  expect_true("Breadth_SampleTransition__t0_vs_t1" %in%
                colnames(rowData(se2)))
})

test_that("transition_matrix_per_sample propagates explicit Uncertain", {
  se <- make_cross_se("A")
  rowData(se)$Breadth_Call__t0[1] <- "Uncertain"
  se <- transition_matrix_per_sample(
    se, feature = "Breadth", timepoints = c("t0", "t1"), verbose = FALSE)
  expect_equal(rowData(se)$Breadth_SampleTransition__t0_vs_t1[1], "Uncertain")
})

test_that("transition_matrix_per_sample gives all ordered pairs when ref=NULL", {
  se <- make_cross_se("A")
  rowData(se)$Breadth_Call__t2 <- rowData(se)$Breadth_Call__t0
  se <- transition_matrix_per_sample(se, feature = "Breadth",
                                     timepoints = c("t0","t1","t2"),
                                     verbose = FALSE)
  cols <- grep("SampleTransition", colnames(rowData(se)), value = TRUE)
  expect_equal(length(cols), 3)  # t0-t1, t0-t2, t1-t2
  expect_true(all(c("t0_vs_t1", "t0_vs_t2", "t1_vs_t2") %in%
                    names(metadata(se)$transitions)))
  expect_error(transition_matrix_per_sample(se, feature = "Breadth",
                                            ref = "nope",
                                            timepoints = c("t0","t1")), "`ref`")
})

test_that("integrate_cross_mark overlays second-mark calls per time point", {
  a <- make_cross_se("A")
  b <- make_cross_se("B")
  res <- integrate_cross_mark(a, b,
                              call_fmt_a = "Breadth_Call__%s",
                              call_fmt_b = "Intensity_Call__%s",
                              mark_a = "Breadth", mark_b = "Intensity",
                              timepoints = c("t0", "t1"),
                              aggregate_ov = "any_super")
  rd <- rowData(res)
  # A domain 1 (1-300) overlaps only B_01 (1-400, t0=Super) -> Super
  expect_equal(rd$Intensity_Call__t0[1], "Intensity_Super_Element")
  expect_equal(rd$Intensity_NOverlaps__t0[1], 1)
  # A domain 2 (101-300) overlaps B_01 too (1-400)
  expect_equal(rd$Intensity_Call__t0[2], "Intensity_Super_Element")
  # A domain 3 (201-300) overlaps B_01 (Super at t0)
  expect_equal(rd$Intensity_Call__t0[3], "Intensity_Super_Element")
  # A domain 4 (301-300 -> 301-600) overlaps B_01 & B_02; t0 B_01=Super,B_02=Typical
  expect_equal(rd$Intensity_Call__t0[4], "Intensity_Super_Element")
  expect_true(rd$Intensity_NOverlaps__t0[4] >= 2)
  # t1: B_01 is NA, B_02 Typical -> A domain 2 (101-300) overlaps B_01 only -> Uncertain
  expect_equal(rd$Intensity_Call__t1[2], "Uncertain")
  # A domain 4 at t1 overlaps B_01(NA) + B_02(Typical) -> any_super = Typical
  expect_equal(rd$Intensity_Call__t1[4], "Intensity_Typical")
  # long table + provenance exist
  expect_false(is.null(metadata(res)$cross_mark_integration))
  expect_equal(metadata(res)$cross_mark_provenance$n_overlap_pairs, 5)
})

test_that("integrate_cross_mark majority rule and validation", {
  a <- make_cross_se("A")
  b <- make_cross_se("B")
  res <- integrate_cross_mark(a, b,
                              call_fmt_a = "Breadth_Call__%s",
                              call_fmt_b = "Intensity_Call__%s",
                              mark_b = "Intensity",
                              timepoints = c("t0", "t1"),
                              aggregate_ov = "majority")
  # A domain 4 at t0 overlaps Super(1) + Typical(1) -> majority: not >0.5 -> Typical
  expect_equal(rowData(res)$Intensity_Call__t0[4], "Intensity_Typical")
  # errors
  expect_error(integrate_cross_mark(a, b, mark_b = "Breadth",
                                    timepoints = c("t0", "t1")),
               "same call prefix")
  expect_error(integrate_cross_mark(a, b, timepoints = "nope"), "Missing call")
})

test_that("integrate_cross_mark: mark_b label can differ from source call prefix", {
  a <- make_cross_se("A")
  b <- make_cross_se("B")
  # se_b carries "Intensity_Call__t0/t1" but we label output "H3K4me3"
  res <- integrate_cross_mark(a, b,
                              call_fmt_a = "Breadth_Call__%s",
                              call_fmt_b = "Intensity_Call__%s",
                              mark_a = "Breadth", mark_b = "H3K4me3",
                              timepoints = c("t0", "t1"),
                              aggregate_ov = "any_super")
  rd <- rowData(res)
  # A domain 1 overlaps B_01 (t0 = Intensity_Super_Element) -> output H3K4me3 label
  expect_equal(rd$H3K4me3_Call__t0[1], "H3K4me3_Super_Element")
  # t1: B_01 is NA, so overlapping -> Uncertain
  expect_equal(rd$H3K4me3_Call__t1[2], "Uncertain")
  # output column exists under mark_b name
  expect_true("H3K4me3_NOverlaps__t0" %in% colnames(rd))
})

test_that("integrate_cross_mark is strand-agnostic and rejects assembly mismatch", {
  a <- make_cross_se("A")[1, ]
  b <- make_cross_se("B")[1, ]
  GenomicRanges::strand(rowRanges(a)) <- "+"
  GenomicRanges::strand(rowRanges(b)) <- "-"
  metadata(a)$signal_contract <- list(genome = "hg38")
  metadata(b)$signal_contract <- list(genome = "hg38")

  res <- integrate_cross_mark(a, b, timepoints = c("t0", "t1"))
  expect_equal(rowData(res)$Intensity_NOverlaps__t0, 1)

  metadata(b)$signal_contract$genome <- "hg19"
  expect_error(integrate_cross_mark(a, b, timepoints = c("t0", "t1")),
               "different genome assemblies")
})
