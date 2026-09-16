test_that("continuous quantile calls retain interpolated cutoffs", {
  x <- setNames(1:4, paste0("d", 1:4))
  for (logged in c(FALSE, TRUE)) {
    for (policy in c("strict", "inclusive")) {
      ans <- .call_super_domains_on_vector(
        x, "Intensity", 0.5, logged, FALSE, tie_policy = policy)
      cutoff <- if (logged) {
        10^unname(stats::quantile(log10(x + 1), 0.5)) - 1
      } else 2.5
      expect_equal(unname(ans$cutoff_value), cutoff)
      expect_identical(ans$Domain_Type,
                       paste0("Intensity_", c("Typical", "Typical",
                                              "Super_Element", "Super_Element")))
    }
  }
  tied <- setNames(c(1, 2, 2, 3), names(x))
  strict <- .call_super_domains_on_vector(
    tied, "Intensity", 0.5, FALSE, FALSE, tie_policy = "strict")
  inclusive <- .call_super_domains_on_vector(
    tied, "Intensity", 0.5, FALSE, FALSE, tie_policy = "inclusive")
  expect_equal(sum(strict$Domain_Type == "Intensity_Super_Element"), 1)
  expect_equal(sum(inclusive$Domain_Type == "Intensity_Super_Element"), 3)
})

.evidence_fixture <- function() {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Intensity = matrix(1, 3, 2)),
    colData = S4Vectors::DataFrame(Condition = c("A", "B")))
  rownames(se) <- paste0("d", 1:3)
  colnames(se) <- c("s1", "s2")
  evidence <- matrix(c("Broad", "Typical", "PeakAbsent",
                       "NoCall", "PeakAbsent", "Broad"), 3, 2,
                     dimnames = dimnames(se))
  reason <- matrix(paste0("reason", 1:6), 3, 2, dimnames = dimnames(se))
  S4Vectors::metadata(se)$breadth_domain_evidence <-
    list(evidence = evidence, reason = reason, group_var = "Condition")
  se
}

test_that("Breadth access follows IDs after row and sample reordering", {
  se <- .evidence_fixture()
  reordered <- se[c(3, 1, 2), c(2, 1)]
  expected <- matrix(c("Broad", "NoCall", "PeakAbsent",
                       "PeakAbsent", "Broad", "Typical"), 3, 2,
                     dimnames = dimnames(reordered))
  expect_identical(get_breadth_evidence(reordered), expected)
  tab <- get_breadth_evidence(reordered, group = "B", long = TRUE)
  expect_identical(tab$Domain_ID, c("d3", "d1", "d2"))
  expect_identical(tab$SampleID, rep("s2", 3))
  expect_identical(tab$Group, rep("B", 3))
  expect_identical(tab$Evidence, c("Broad", "NoCall", "PeakAbsent"))
  expect_identical(tab$Reason, c("reason6", "reason4", "reason5"))
  small <- se[c(3, 1), 2, drop = FALSE]
  expect_identical(get_breadth_evidence(small, type = "reason")[, 1],
                   c(d3 = "reason6", d1 = "reason4"))
  expect_equal(nrow(get_breadth_evidence(se[FALSE, ], long = TRUE)), 0)
})

test_that("Breadth access rejects ambiguous or inconsistent identifiers", {
  se <- .evidence_fixture()
  rownames(se)[1] <- "unknown"
  expect_error(get_breadth_evidence(se), "not aligned")
  se <- .evidence_fixture()
  colnames(S4Vectors::metadata(se)$breadth_domain_evidence$reason) <- c("s2", "s1")
  expect_error(get_breadth_evidence(se), "not aligned")
  se <- .evidence_fixture()
  rownames(se) <- c("d1", "d1", "d3")
  expect_error(get_breadth_evidence(se), "not aligned")
})

test_that("called objects remain ID-aligned after domain and sample reordering", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           verbose = FALSE)
  se <- call_super_domains(se, feature = "Breadth", mode = "per_group",
                           group_var = "Condition", verbose = FALSE)
  reordered <- se[c(3, 1, 2), rev(seq_len(ncol(se)))]
  expect_true(validate_epiportrait_object(reordered))
  calls <- get_replicate_calls(reordered, feature = "Intensity")
  evidence <- get_breadth_evidence(reordered)
  expect_identical(rownames(calls), rownames(reordered))
  expect_identical(colnames(calls), colnames(reordered))
  expect_identical(dimnames(evidence), dimnames(reordered))
})

test_that("ORA ratios and odds ratios use the tested annotation universe", {
  d <- data.frame(ID = "GO:test", Description = "example", Count = 2,
                  GeneRatio = "2/5", BgRatio = "10/100", pvalue = 0.1,
                  p.adjust = 0.2, qvalue = 0.2, geneID = "a/b")
  # Supplied genes can include IDs without annotation in this GO ontology.
  ans <- .tidy_ora_result(d, "BP", n_fg = 8, n_universe = 150,
                          padj_cutoff = 0.05)
  expect_equal(ans$gene_ratio, 0.4)
  expect_equal(ans$bg_ratio, 0.1)
  expect_equal(ans$effect, log2((2 * 87) / (3 * 8)))
  expect_equal(ans$p_value, d$pvalue)
  expect_equal(ans$p_adjust, d$p.adjust)
})
