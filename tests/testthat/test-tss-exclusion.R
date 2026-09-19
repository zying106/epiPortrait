test_that("filter_promoter_peaks overlap vs contained semantics", {
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(850, 980, 5000), width = 100))   # 850, 980, 5000
  tss <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 1000),
                                strand = "+")

  ov <- filter_promoter_peaks(gr, tss = tss, upstream = 100, downstream = 100,
                              mode = "overlap")
  ct <- filter_promoter_peaks(gr, tss = tss, upstream = 100, downstream = 100,
                              mode = "contained")

  # overlap removes 850 (touches window) and 980 (inside); keeps 5000
  expect_length(ov, 1L)
  expect_equal(GenomicRanges::start(ov), 5000)
  # contained keeps 850 (crosses the window edge), removes 980, keeps 5000
  expect_length(ct, 2L)
  expect_setequal(GenomicRanges::start(ct), c(850, 5000))
})

test_that("filter_promoter_peaks records provenance metadata", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, width = 1))
  tss <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 1000),
                                strand = "+")
  out <- filter_promoter_peaks(gr, tss = tss, upstream = 500, downstream = 500,
                               mode = "contained")
  prov <- S4Vectors::metadata(out)$promoter_exclusion
  expect_type(prov, "list")
  expect_equal(prov$mode, "contained")
  expect_equal(prov$n_input, 1L)
  expect_equal(prov$n_removed, 1L)
  expect_equal(prov$upstream, 500)
})

test_that("filter_promoter_peaks validates arguments", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, width = 1))
  tss <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 1000),
                                strand = "+")
  expect_error(filter_promoter_peaks("nope"), "GRanges")
  expect_error(filter_promoter_peaks(gr, tss = "nope"), "tss must be")
  expect_error(filter_promoter_peaks(gr, tss = tss, upstream = -1),
               "non-negative")
  expect_error(filter_promoter_peaks(gr, tss = tss, mode = "bogus"))
})

test_that("tss_from_rose parses a refGene table strand-aware", {
  f <- tempfile(fileext = ".ucsc")
  writeLines(c("#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tname2",
               "0\tNM_000001\tchr1\t+\t1000\t2000\tGENE1",
               "0\tNM_000002\tchr1\t-\t5000\t6000\tGENE2"), f)
  tss <- tss_from_rose("hg38", f)
  expect_s4_class(tss, "GRanges")
  expect_equal(length(tss), 2L)
  # + strand TSS = txStart + 1 ; - strand TSS = txEnd
  expect_equal(GenomicRanges::start(tss)[1], 1001L)
  expect_equal(GenomicRanges::start(tss)[2], 6000L)
  expect_equal(as.character(S4Vectors::mcols(tss)$gene), c("GENE1", "GENE2"))
  expect_equal(S4Vectors::metadata(tss)$tss_source$source,
               "ROSE refGene table")
})

test_that("tss_from_rose rejects bad input", {
  expect_error(tss_from_rose("hg38", tempfile()), "existing ROSE refGene")
  f <- tempfile(fileext = ".ucsc")
  writeLines(c("#bin\tname\tchrom", "0\tNM_1\tchr1"), f)
  expect_error(tss_from_rose("hg38", f), "refGene table with columns")
})

test_that("revert_multi_tss reverts regions spanning > max_tss genes", {
  orig <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(1000, 1100, 1900, 2800, 2900, 4100), width = 50))
  stitched <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(1000, 4000), c(3000, 4200)))
  tss <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(c(1500, 2000, 2500, 4500), width = 1))
  tss$gene <- c("A", "B", "C", "D")   # A,B,C inside first region

  out <- revert_multi_tss(stitched, tss, orig, max_tss = 2)
  prov <- S4Vectors::metadata(out)$multi_tss_revert
  expect_equal(prov$n_reverted, 1L)
  # the first region is replaced by its 5 constituent peaks; second region kept
  expect_equal(length(out), 6L)
  expect_false(any(GenomicRanges::width(out) == 2000L))
})

test_that("revert_multi_tss keeps a region spanning at most max_tss genes", {
  orig <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, width = 50))
  stitched <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 3000))
  tss <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1500, 2500),
                                                         width = 1))
  tss$gene <- c("A", "B")
  out <- revert_multi_tss(stitched, tss, orig, max_tss = 2)
  expect_equal(length(out), 1L)
  expect_equal(S4Vectors::metadata(out)$multi_tss_revert$n_reverted, 0L)
})

test_that("revert_multi_tss validates arguments", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  tss <- GenomicRanges::GRanges("chr1", IRanges::IRanges(5, 5))
  expect_error(revert_multi_tss("x", tss, gr), "GRanges")
  expect_error(revert_multi_tss(gr, tss, gr, max_tss = -1), "non-negative")
})
