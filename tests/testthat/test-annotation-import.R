test_that("standardized external annotations replace and rebuild all views", {
  data(example_se)
  links <- data.frame(
    domain_id = rownames(example_se)[1:4],
    gene_id = c("g1", "g2", "g3", "g4"),
    gene_symbol = c("A", "B", "C", "D"),
    relation_type = c("promoter_overlap", "nearest_tss",
                      "gene_body_overlap", "fully_contained"),
    distance_to_tss_bp = c(NA, -4500, NA, NA),
    external_relation = c("Promoter", "Downstream", "Intron", "Exon"),
    stringsAsFactors = FALSE)

  out <- import_domain_annotations(example_se, links, source = "ChIPseeker")
  raw <- S4Vectors::metadata(out)$domain_gene_links
  dedup <- S4Vectors::metadata(out)$domain_gene_links_dedup
  summary <- S4Vectors::metadata(out)$annotation_summary

  expect_equal(nrow(raw), 4L)
  expect_equal(nrow(dedup), 4L)
  expect_equal(nrow(summary), nrow(example_se))
  expect_equal(raw$annotation_source, rep("ChIPseeker", 4L))
  expect_equal(raw$evidence_source, rep("external", 4L))
  expect_equal(raw$external_relation,
               c("Promoter", "Downstream", "Intron", "Exon"))
  expect_equal(rowData(out)$primary_genomic_context[1:4],
               c("Promoter-associated", "Intergenic",
                 "Gene-body-associated", "Gene-body-associated"))
  expect_equal(rowData(out)$nearest_tss_distance_bp[2], -4500)
  expect_equal(dedup$best_tier, c(4L, 2L, 1L, 1L))
  expect_equal(S4Vectors::metadata(out)$annotation_provenance$annotation_mode,
               "external_replace")
  expect_equal(
    S4Vectors::metadata(out)$annotation_import_provenance[[1]]$source,
    "ChIPseeker")
})

test_that("append preserves existing links and records both sources", {
  data(example_se)
  first <- data.frame(
    domain_id = rownames(example_se)[1], gene_id = "g1",
    gene_symbol = "A", relation_type = "promoter_overlap")
  second <- data.frame(
    domain_id = rownames(example_se)[1], gene_id = "g2",
    gene_symbol = "B", relation_type = "gene_body_overlap")

  out <- import_domain_annotations(example_se, first, source = "tool_a")
  out <- import_domain_annotations(out, second, source = "tool_b",
                                   mode = "append")
  raw <- S4Vectors::metadata(out)$domain_gene_links

  expect_equal(nrow(raw), 2L)
  expect_setequal(raw$annotation_source, c("tool_a", "tool_b"))
  expect_equal(S4Vectors::metadata(out)$annotation_summary$n_linked_gene[1], 2L)
  expect_equal(length(S4Vectors::metadata(out)$annotation_import_provenance), 2L)
})

test_that("BEDPE-like external links retain record-aware score summaries", {
  data(example_se)
  links <- data.frame(
    domain_id = rep(rownames(example_se)[1], 2),
    gene_id = rep("g1", 2),
    gene_symbol = rep("A", 2),
    relation_type = rep("bedpe_promoter_contact", 2),
    bedpe_record_id = c("loop_1", "loop_2"),
    contact_score = c(2, 3))

  out <- import_domain_annotations(example_se, links, source = "ABC")
  dedup <- S4Vectors::metadata(out)$domain_gene_links_dedup

  expect_equal(dedup$bedpe_support_count, 2L)
  expect_equal(dedup$bedpe_contact_score, 5)
  expect_equal(dedup$bedpe_contact_score_max, 3)
  expect_equal(rowData(out)$bedpe_contact_score[1], 5)
  expect_equal(rowData(out)$n_bedpe_contact_gene[1], 1L)
})

test_that("bedpe rows require a stable bedpe_record_id", {
  data(example_se)
  no_id <- data.frame(
    domain_id = rep(rownames(example_se)[1], 2),
    gene_id = rep("g1", 2),
    relation_type = rep("bedpe_promoter_contact", 2),
    contact_score = c(1, 2))
  expect_error(
    import_domain_annotations(example_se, no_id, source = "ABC"),
    "bedpe_record_id is required")

  # A stable user-supplied id keeps record-aware support counting intact.
  with_id <- transform(no_id, bedpe_record_id = c("r1", "r2"))
  out <- import_domain_annotations(example_se, with_id, source = "ABC")
  dedup <- S4Vectors::metadata(out)$domain_gene_links_dedup
  expect_equal(dedup$bedpe_support_count, 2L)
  expect_equal(dedup$bedpe_contact_score, 3)
})

test_that("import rejects non-canonical or unmatched annotations", {
  data(example_se)
  bad_relation <- data.frame(
    domain_id = rownames(example_se)[1], gene_id = "g1",
    relation_type = "Intron")
  expect_error(
    import_domain_annotations(example_se, bad_relation, source = "tool"),
    "Unsupported relation_type")

  bad_domain <- transform(bad_relation,
                          domain_id = "not_in_object",
                          relation_type = "gene_body_overlap")
  expect_error(
    import_domain_annotations(example_se, bad_domain, source = "tool"),
    "not present in se")
  expect_error(
    import_domain_annotations(example_se, bad_relation[, -2], source = "tool"),
    "missing required")
  expect_error(
    import_domain_annotations(example_se, bad_relation, source = ""),
    "source")
  expect_error(
    import_domain_annotations(example_se, transform(bad_relation,
      relation_type = "gene_body_overlap"), source = "tool", mode = "append"),
    "requires existing")
})

test_that("file input works and imported links invalidate stored enrichment", {
  data(example_se)
  links <- data.frame(
    domain_id = rownames(example_se)[1], gene_id = "g1",
    relation_type = "nearest_tss", distance_to_tss_bp = 100)
  path <- tempfile(fileext = ".tsv")
  utils::write.table(links, path, sep = "\t", quote = FALSE, row.names = FALSE)
  S4Vectors::metadata(example_se)$enrichment <- list(old = TRUE)
  S4Vectors::metadata(example_se)$enrichment_comparison <- list(old = TRUE)

  expect_warning(
    out <- import_domain_annotations(example_se, path, source = "GREAT"),
    "changed the candidate-gene universe")
  expect_null(S4Vectors::metadata(out)$enrichment)
  expect_null(S4Vectors::metadata(out)$enrichment_comparison)
  expect_equal(S4Vectors::metadata(out)$domain_gene_links$annotation_source,
               "GREAT")
})

test_that("external annotation provenance is exported", {
  data(example_se)
  links <- data.frame(
    domain_id = rownames(example_se)[1], gene_id = "g1",
    relation_type = "promoter_overlap")
  out <- import_domain_annotations(example_se, links, source = "external_tool")
  path <- export_epiportrait_results(
    out, outdir = tempfile("epi_external_export"), save_object = FALSE)
  provenance <- readLines(file.path(path, "annotation", "provenance.txt"))
  expect_true(any(grepl("annotation_import_provenance", provenance,
                        fixed = TRUE)))
})
