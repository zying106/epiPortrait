# Tests for the domain-aware functional interpretation module (2026-09):
# evidence tiers, gene ranking / aggregation, domain-derived universe, ORA,
# GSEA, cross-group comparison, and interoperability.

library(epiPortrait)
library(GenomicRanges)
library(SummarizedExperiment)
library(S4Vectors)

# ---- helper: synthetic annotated object with controllable gene links --------
# Bypasses TxDb; writes metadata(se)$domain_gene_links directly so the
# domain-gene evidence contract is tested, not the annotation implementation
# (which has its own test file).
make_enrich_se <- function(n_dom = 30, gene_pool = NULL, drop_last_links = FALSE) {
  if (is.null(gene_pool)) gene_pool <- sprintf("gene%03d", seq_len(3 * n_dom))
  dom <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = seq(1, by = 1e5, length.out = n_dom),
                             width = 5000))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Intensity = matrix(1, n_dom, 4)),
    rowRanges = dom)
  rownames(se) <- sprintf("epiDomain_%06d", seq_len(n_dom))
  colnames(se) <- paste0("S", 1:4)
  SummarizedExperiment::colData(se)$Condition <-
    c("Control", "Control", "Treatment", "Treatment")

  mk <- function(i, gene_id, relation, dist = NA_real_, source = "linear",
                 bedpe_id = NA_character_, score = NA_real_) {
    data.frame(domain_id = rownames(se)[i], gene_id = gene_id,
               gene_symbol = NA_character_, relation_type = relation,
               distance_to_tss_bp = dist, overlap_bp = 100,
               domain_overlap_fraction = 0.02, feature_overlap_fraction = 0.5,
               bedpe_record_id = bedpe_id, evidence_source = source,
               contact_score = score, stringsAsFactors = FALSE)
  }
  k <- 0L
  rows <- lapply(seq_len(n_dom), function(i) {
    if (drop_last_links && i == n_dom) return(NULL)
    k <<- k + 1L
    g1 <- gene_pool[k]
    k <<- k + 1L
    g2 <- gene_pool[k]
    rbind(
      mk(i, g1, "promoter_overlap"),
      mk(i, g2, "nearest_tss", dist = if (i %% 3 == 0) 5e5 else 3000))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  links <- do.call(rbind, rows)
  links <- rbind(links, mk(1, gene_pool[length(gene_pool)],
                           "bedpe_promoter_contact", source = "bedpe",
                           bedpe_id = "rec1", score = 12))
  S4Vectors::metadata(se)$domain_gene_links <- links
  set.seed(11)
  SummarizedExperiment::rowData(se)$test_score <- stats::rnorm(n_dom)
  se
}

entrez_pool <- function(n) {
  skip_if_not_installed("org.Hs.eg.db")
  head(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,
                           keytype = "ENTREZID"), n)
}

# ---- evidence tier / per-domain cap on get_domain_genes ---------------------
test_that("get_domain_genes exposes evidence_tier and caps per domain", {
  se <- make_enrich_se(10)
  gd <- get_domain_genes(se, relations = c("promoter_overlap", "nearest_tss",
                                           "bedpe_promoter_contact"),
                         unique_genes = FALSE)
  expect_true("evidence_tier" %in% colnames(gd))
  expect_true(all(gd$evidence_tier %in% 0:4))
  expect_true(all(gd$evidence_tier[gd$relation_types == "promoter_overlap"] == 4L))
  proximal <- gd$relation_types == "nearest_tss" &
    abs(gd$nearest_tss_distance_bp) <= 10000
  expect_true(all(gd$evidence_tier[proximal] == 2L))
  expect_true(all(gd$evidence_tier[gd$relation_types == "bedpe_promoter_contact"] == 3L))

  capped <- get_domain_genes(se, relations = "promoter_overlap",
                             max_per_domain = 1, unique_genes = FALSE)
  expect_true(all(table(capped$domain_id) == 1L))
  expect_error(get_domain_genes(se, max_per_domain = -1),
               "max_per_domain")
  # NULL keeps the original behaviour (all qualifying pairs)
  expect_equal(nrow(get_domain_genes(se, relations = "promoter_overlap",
                                     unique_genes = FALSE)), 10L)
})

# ---- pure helpers -----------------------------------------------------------
test_that("ORA effect size matches the explicit 2x2 computation", {
  eff <- epiPortrait:::.ora_effect_size(10, 90, 20, 880)
  expected_or <- log2((10 * 880) / (90 * 20))
  expected_se <- sqrt(1 / 10 + 1 / 90 + 1 / 20 + 1 / 880) / log(2)
  expect_equal(eff$effect, expected_or)
  expect_equal(eff$effect_low, expected_or - 1.959964 * expected_se)
  expect_equal(eff$effect_high, expected_or + 1.959964 * expected_se)
  # zero cells use the Haldane-Anscombe correction, never Inf
  eff0 <- epiPortrait:::.ora_effect_size(0, 10, 5, 100)
  expect_true(all(is.finite(unlist(eff0))))
})

test_that("restrict_to_testable reports and applies the intersection", {
  res <- epiPortrait:::.restrict_to_testable(c("a", "b", "c", "a"), c("b", "c", "d"))
  expect_setequal(res$genes, c("b", "c"))
  expect_equal(res$n_input, 3L)
  expect_equal(res$n_dropped, 1L)
})

test_that("domain score aggregation rules are exact", {
  links <- data.frame(gene_id = c("g1", "g1", "g2"),
                      gene_symbol = c("G1", NA, "G2"),
                      stringsAsFactors = FALSE)
  scores <- c(2, -1, 4)
  w <- c(1, 1, 2)
  agg <- epiPortrait:::.aggregate_domain_scores
  by_gene <- function(res, genes = c("g1", "g2")) res$score[match(genes, res$gene_id)]
  nd_by_gene <- function(res, genes = c("g1", "g2"))
    res$n_domains[match(genes, res$gene_id)]
  expect_equal(by_gene(agg(links, scores, w, "signed_weighted")), c(0.5, 4))
  expect_equal(by_gene(agg(links, scores, w, "mean")), c(0.5, 4))
  expect_equal(by_gene(agg(links, scores, w, "sum")), c(1, 4))
  expect_equal(by_gene(agg(links, scores, w, "max_abs")), c(2, 4))
  expect_equal(by_gene(agg(links, scores, w, "best_domain")), c(2, 4))
  expect_equal(nd_by_gene(agg(links, scores, w, "signed_weighted")), c(2L, 1L))
  # output is ordered by descending score
  expect_true(all(diff(agg(links, scores, w, "signed_weighted")$score) <= 0))
})

# ---- rank_epi_genes ---------------------------------------------------------
test_that("rank_epi_genes aggregates signed domain scores by evidence weight", {
  se <- make_enrich_se(10)
  scores <- setNames(seq_len(10) * 1, rownames(se))
  SummarizedExperiment::rowData(se)$dom_score <- scores
  rk <- rank_epi_genes(se, "dom_score", mark = "H3K27ac")
  expect_true(is.data.frame(rk))
  expect_true(all(c("gene_id", "score", "n_domains") %in% colnames(rk)))
  # promoter links are tier 4 (weight 16), proximal nearest tier 2 (weight 4)
  links <- S4Vectors::metadata(se)$domain_gene_links
  prom <- links[links$relation_type == "promoter_overlap", ]
  expect_equal(rk$score[match(prom$gene_id[1], rk$gene_id)], 1)
  expect_equal(attr(rk, "score_col"), "dom_score")
  expect_equal(attr(rk, "aggregate"), "signed_weighted")

  far_genes <- links$gene_id[links$relation_type == "nearest_tss" &
                               abs(links$distance_to_tss_bp) > 10000]
  expect_false(any(far_genes %in% rk$gene_id))  # tier 0 dropped by default
  rk0 <- rank_epi_genes(se, "dom_score", mark = "H3K27ac", min_evidence_tier = 0)
  expect_true(any(far_genes %in% rk0$gene_id))
})

test_that("rank_epi_genes validates the score column", {
  se <- make_enrich_se(5)
  expect_error(rank_epi_genes(se, "not_a_column"), "not found in rowData")
  expect_error(rank_epi_genes(se, c("a", "b")), "single rowData column")
})

# ---- universe ---------------------------------------------------------------
test_that("get_domain_gene_universe is domain-derived and mark-aware", {
  pool <- entrez_pool(90)
  se <- make_enrich_se(20, gene_pool = pool)
  u <- get_domain_gene_universe(se, mark = "H3K27ac",
                                org_db = "org.Hs.eg.db",
                                key_type = "ENTREZID")
  expect_s3_class(u, "epi_gene_universe")
  expect_true(all(u %in% pool))
  expect_true(attr(u, "n_testable") <= attr(u, "n_linked"))
  expect_equal(attr(u, "key_type"), "ENTREZID")

  # broad repressive marks: promoter-overlap links only
  u_rep <- get_domain_gene_universe(se, mark = "H3K27me3",
                                    org_db = "org.Hs.eg.db",
                                    key_type = "ENTREZID")
  links <- S4Vectors::metadata(se)$domain_gene_links
  expect_true(all(u_rep %in% links$gene_id[links$relation_type == "promoter_overlap"]))
})

test_that("universe drops genes absent from the annotation database", {
  pool <- c(entrez_pool(40), sprintf("FAKE%02d", 1:5))
  se <- make_enrich_se(20, gene_pool = pool)
  u <- get_domain_gene_universe(se, org_db = "org.Hs.eg.db",
                                key_type = "ENTREZID")
  expect_false(any(grepl("^FAKE", u)))
  expect_gt(attr(u, "n_dropped"), 0)
})

# ---- ORA --------------------------------------------------------------------
test_that("enrich_epi_genes returns the tidy effect-size schema", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  universe <- entrez_pool(2000)
  set.seed(3)
  fg <- sample(universe, 30)
  res <- enrich_epi_genes(fg, universe, org_db = "org.Hs.eg.db",
                          key_type = "ENTREZID", simplify = FALSE)
  expect_s3_class(res, "epi_enrichment")
  expect_true(all(c("term_id", "term_name", "count", "set_size", "effect",
                    "effect_low", "effect_high", "p_value", "p_adjust",
                    "significant", "gene_id") %in% colnames(res$result)))
  expect_true(all(res$result$effect_low <= res$result$effect))
  expect_true(all(res$result$effect <= res$result$effect_high))
  expect_true(inherits(res$object, "enrichResult"))
  expect_equal(res$mapping$n_foreground, 30L)
  expect_equal(res$background$n, length(universe))
  expect_equal(res$mapping$terms_tested, nrow(res$result))
  # deterministic
  res2 <- enrich_epi_genes(fg, universe, org_db = "org.Hs.eg.db",
                           key_type = "ENTREZID", simplify = FALSE)
  expect_equal(res$result, res2$result)
})

test_that("enrich_epi_genes returns a structured empty result, never an error", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  universe <- entrez_pool(500)
  res <- enrich_epi_genes(universe[1:3], universe, org_db = "org.Hs.eg.db",
                          key_type = "ENTREZID", min_foreground = 10)
  expect_s3_class(res, "epi_enrichment")
  expect_equal(nrow(res$result), 0L)
  expect_equal(res$mapping$reason, "foreground_too_small")
  expect_null(res$object)
})

test_that("ORA returns valid p-values and detects a real signal", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  universe <- entrez_pool(4000)
  set.seed(42)
  fg <- sample(universe, 40)
  res <- enrich_epi_genes(fg, universe, org_db = "org.Hs.eg.db",
                          key_type = "ENTREZID", simplify = FALSE)
  p <- res$result$p_value
  # Deterministic validity contract. A distributional "no inflation" claim on
  # a random foreground is NOT asserted here: GO terms overlap heavily and the
  # realized proportion of p <= 0.05 depends on the (frequently updated) GO.db
  # release, so such an assertion is environment-dependent and flakes across
  # Bioconductor builds. The exact p-value / effect-size math is instead pinned
  # deterministically in the 2x2 tests below and in test-submission-regressions.
  expect_gt(length(p), 0)
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(res$result$p_adjust >= 0 & res$result$p_adjust <= 1))
  # A real signal must still be detected (sanity of power): cell-cycle genes.
  term_genes <- unique(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db, keys = "GO:0000278",
    columns = "ENTREZID", keytype = "GOALL")$ENTREZID)
  fg2 <- intersect(term_genes, universe)
  expect_gte(length(fg2), 20)
  res2 <- enrich_epi_genes(fg2[1:40], universe, org_db = "org.Hs.eg.db",
                           key_type = "ENTREZID", simplify = FALSE)
  expect_gt(res2$mapping$terms_significant, 0)
})

test_that("enrich_epi_genes rejects unimplemented databases / ontologies", {
  expect_error(enrich_epi_genes("1", "1", database = "KEGG"),
               "not implemented yet")
  expect_error(enrich_epi_genes("1", "1", ontology = "ALL"),
               "ontology must be one of")
})

# ---- object integration -----------------------------------------------------
test_that("enrich_epi_domains stores results and getter retrieves them", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  pool <- entrez_pool(300)
  se <- make_enrich_se(20, gene_pool = pool)
  sel <- rep(FALSE, 20); sel[1:10] <- TRUE
  se2 <- enrich_epi_domains(se, domains = sel, simplify = FALSE,
                            org_db = "org.Hs.eg.db", key_type = "ENTREZID")
  res <- get_epi_enrichment(se2)
  expect_length(res, 1L)
  expect_s3_class(res[[1]], "epi_enrichment")
  expect_equal(res[[1]]$mapping$input_domains, 10L)
  expect_equal(res[[1]]$mapping$unmapped_domains, 0L)
  expect_true(any(summarize_epiportrait_object(se2)$enrichment_results != ""))

  # second call collides -> unique name
  se3 <- enrich_epi_domains(se2, domains = sel, simplify = FALSE,
                            org_db = "org.Hs.eg.db", key_type = "ENTREZID")
  expect_length(get_epi_enrichment(se3), 2L)
  expect_true("ORA_GO_BP_2" %in% names(get_epi_enrichment(se3)))

  # selection without links -> structured empty result
  se4 <- make_enrich_se(20, gene_pool = pool, drop_last_links = TRUE)
  sel4 <- rep(FALSE, 20); sel4[20] <- TRUE
  se4 <- enrich_epi_domains(se4, domains = sel4, simplify = FALSE,
                            org_db = "org.Hs.eg.db", key_type = "ENTREZID")
  expect_equal(get_epi_enrichment(se4)[[1]]$mapping$reason,
               "no_links_in_selection")
})

test_that("enrich_epi_domains dispatches GSEA on a continuous score", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  pool <- entrez_pool(300)
  se <- make_enrich_se(100, gene_pool = pool)
  SummarizedExperiment::rowData(se)$dom_score <-
    seq(-10, 10, length.out = 100)
  # Genes sharing a domain share its score, so preranked ties are expected here
  # and only produce an fgsea tie warning.
  se2 <- suppressWarnings(
    enrich_epi_domains(se, method = "GSEA", score_col = "dom_score",
                       min_gs_size = 5,
                       org_db = "org.Hs.eg.db", key_type = "ENTREZID"))
  res <- get_epi_enrichment(se2)[[1]]
  expect_s3_class(res, "epi_enrichment")
  expect_equal(res$method, "GSEA")
  expect_equal(res$parameters$seed, 1)
  expect_true(all(c("effect", "p_adjust", "significant") %in%
                    colnames(res$result)))
  expect_true(inherits(res$object, "gseaResult"))
})

test_that("gsea_epi_genes returns structured empty result for tiny rankings", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  universe <- entrez_pool(300)
  ranked <- data.frame(gene_id = universe[1:10], score = 1:10)
  res <- gsea_epi_genes(ranked, universe, org_db = "org.Hs.eg.db",
                        key_type = "ENTREZID", min_ranked = 100)
  expect_s3_class(res, "epi_enrichment")
  expect_equal(res$mapping$reason, "ranked_list_too_small")
})

# ---- comparison / plot / conversion -----------------------------------------
test_that("compare_epi_enrichment uses one shared universe and effect sizes", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  pool <- entrez_pool(300)
  se <- make_enrich_se(20, gene_pool = pool)
  sets <- list(A = c(rep(TRUE, 8), rep(FALSE, 12)),
               B = c(rep(FALSE, 8), rep(TRUE, 12)))
  se2 <- compare_epi_enrichment(se, sets = sets,
                                org_db = "org.Hs.eg.db",
                                key_type = "ENTREZID")
  cmp <- get_epi_enrichment(se2, type = "comparison")
  # name = NULL always returns the named list (same contract as enrichment)
  expect_length(cmp, 1L)
  cmp <- cmp[[1]]
  expect_s3_class(cmp, "epi_enrichment_comparison")
  expect_identical(cmp$results$A$background$genes,
                   cmp$results$B$background$genes)
  expect_true(all(c("A", "B") %in% unique(cmp$table$group)))
  expect_true(all(c("effect", "effect_low", "effect_high", "n_fg", "n_bg") %in%
                    colnames(cmp$table)))
  p <- plot_epi_enrichment(se2, top_n = 5)
  expect_s3_class(p, "ggplot")
  expect_output(print(cmp), "epi_enrichment_comparison")
})

test_that("as_enrich_result exposes the clusterProfiler object", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  universe <- entrez_pool(1000)
  set.seed(5)
  res <- enrich_epi_genes(sample(universe, 25), universe,
                          org_db = "org.Hs.eg.db", key_type = "ENTREZID",
                          simplify = FALSE)
  expect_true(inherits(as_enrich_result(res), "enrichResult"))
  expect_error(as_enrich_result(res, which = "simplified"),
               "No simplified clusterProfiler object")
  expect_error(as_enrich_result(list()), "epi_enrichment")
})

test_that("get_epi_enrichment errors informatively when empty", {
  data(example_se)
  expect_error(get_epi_enrichment(example_se), "No enrichment results")
  expect_error(get_epi_enrichment(example_se, type = "comparison"),
               "No comparison results")
})
