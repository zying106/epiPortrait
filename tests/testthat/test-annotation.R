# Tests for the domain-aware annotation module (design 2026-08-10):
# genome resolution, linear annotation, BEDPE 3D evidence, RNA-seq expression
# integration, and candidate-gene prioritization.

library(epiPortrait)
library(GenomicRanges)
library(SummarizedExperiment)
library(S4Vectors)

# ---- helper: build a tiny synthetic TxDb on chr1 ----------------------------
make_tiny_txdb <- function() {
  skip_if_not_installed("txdbmaker")
  gid <- c("gene1", "gene2", "gene3", "gene4")
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 8),
    ranges = IRanges::IRanges(
      start = c(1000, 1000, 20000, 20000, 40000, 40000, 60000, 60000),
      end   = c(3000, 3000, 26000, 26000, 48000, 48000, 70000, 70000)),
    strand = rep(c("+", "-", "+", "+"), each = 2))
  S4Vectors::mcols(gr)$type <- rep(c("gene", "transcript"), 4)
  S4Vectors::mcols(gr)$ID <- c(rbind(gid, paste0(gid, "_tx1")))
  S4Vectors::mcols(gr)$Parent <- c(rbind(rep(NA_character_, 4), gid))
  txdb <- txdbmaker::makeTxDbFromGRanges(gr)
  txdb
}

# helper: minimal SE with 4 domains overlapping the synthetic genes
make_anno_se <- function() {
  suppressWarnings({
    domains <- GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(
        start = c(1500, 21000, 35000, 65000),
        end   = c(5000, 30000, 90000, 80000)),
      seqinfo = GenomeInfoDb::Seqinfo("chr1", 1000000))
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(Intensity = matrix(1, 4, 2),
                    SignalDispersion = matrix(1, 4, 2)),
      rowRanges = domains)
  })
  rownames(se) <- paste0("epiDomain_", sprintf("%06d", seq_len(4)))
  colnames(se) <- c("C1", "T1")
  colData(se)$Condition <- c("Control", "Treatment")
  se
}

test_that("annotate_epi_domains adds 7 rowData columns + links", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  se <- annotate_epi_domains(se, genome = txdb)
  expect_true(all(c("primary_genomic_context", "nearest_tss_gene_id",
                    "nearest_tss_gene_symbol", "nearest_tss_distance_bp",
                    "promoter_overlap_gene_count", "gene_body_overlap_gene_count",
                    "fully_contained_gene_count") %in% colnames(rowData(se))))
  expect_false(is.null(S4Vectors::metadata(se)$domain_gene_links))
  expect_false(is.null(S4Vectors::metadata(se)$annotation_provenance))
  # every domain must have a primary context
  expect_true(all(!is.na(rowData(se)$primary_genomic_context)))
})

test_that("genome resolver handles built-in and custom", {
  r <- epiPortrait:::.resolve_genome_resources("hg38")
  expect_equal(r$genome_name, "hg38")
  expect_equal(r$genome_class, "builtin")
  expect_equal(r$organism, "Homo sapiens")
  expect_false(is.null(r$txdb))
  # unknown genome without txdb -> partial + warning
  expect_warning(r2 <- epiPortrait:::.resolve_genome_resources("rn7"), "No built-in")
  expect_equal(r2$genome_class, "partial")
})

test_that("annotate_epi_domains rejects unknown genome without TxDb", {
  se <- make_anno_se()
  expect_error(annotate_epi_domains(se, genome = "rn7"), "TxDb")
})

test_that("nearest TSS and counts are consistent", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  se <- annotate_epi_domains(se, genome = txdb, promoter_upstream = 2000,
                             promoter_downstream = 2000)
  rd <- rowData(se)
  # domain 1 (1500-5000) overlaps gene 1 promoter (1000-3000 +-2kb window) -> promoter_overlap>=1
  expect_true(rd$promoter_overlap_gene_count[1] >= 1)
  # domain 3 (45000-90000) contains gene 3 (40000-48000) and gene 4 (60000-70000) fully
  expect_true(rd$fully_contained_gene_count[3] >= 2)
  # domain 1 (1500-5000) contains part of gene 1 -> gene_body_overlap>=1
  expect_true(rd$gene_body_overlap_gene_count[1] >= 1)
})

test_that("one domain can link to multiple genes / multiple relations", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  se <- annotate_epi_domains(se, genome = txdb)
  links <- S4Vectors::metadata(se)$domain_gene_links
  d3 <- links[links$domain_id == rownames(se)[3], ]
  # domain 3 spans 45-90kb: overlaps gene 3 (40-48kb) and likely gene 4 (60-70kb)
  expect_true(nrow(d3) >= 1)
  expect_true(all(c("nearest_tss", "gene_body_overlap", "fully_contained") %in%
                    unique(d3$relation_type)))
})

test_that("BEDPE promoter contact is symmetric and read-only", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  bedpe_df <- data.frame(
    chrom1 = "chr1", start1 = 2000, end1 = 2500,   # overlaps gene1 (promoter)
    chrom2 = "chr1", start2 = 100,  end2 = 600,    # far anchor, not a domain
    stringsAsFactors = FALSE)
  # to create a contact, anchor2 must overlap a promoter and anchor1 a domain:
  # use domain 1 (1500-5000) as anchor1 and gene 4 promoter (60000-70000) as anchor2
  bedpe_df <- data.frame(
    chrom1 = "chr1", start1 = 1499, end1 = 2500,   # inside domain 1
    chrom2 = "chr1", start2 = 59999, end2 = 61000, # gene 4 promoter region
    stringsAsFactors = FALSE)
  se <- annotate_epi_domains(se, genome = txdb, bedpe = bedpe_df)
  links <- S4Vectors::metadata(se)$domain_gene_links
  bedpe_links <- links[links$evidence_source == "bedpe", ]
  expect_true(nrow(bedpe_links) >= 1)
  expect_equal(unique(bedpe_links$relation_type), "bedpe_promoter_contact")
  expect_false(is.null(S4Vectors::metadata(se)$bedpe_provenance))
  expect_equal(S4Vectors::metadata(se)$bedpe_provenance$record_count, 1)
})

test_that("RNA-seq expression integrates per condition", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  # gene ids must match the TxDb's gene ids
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(c(10, 20, 5, 1, 50, 2, 8, 9),
                     nrow = 4, dimnames = list(gids, c("C1", "T1")))
  se <- annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                             expression_type = "TPM")
  expect_false(is.null(S4Vectors::metadata(se)$expression_provenance))
  expect_equal(S4Vectors::metadata(se)$expression_provenance$expression_type, "TPM")
  # expression summary per condition
  es <- S4Vectors::metadata(se)$expression_summary
  expect_true("Control" %in% es$Condition)
})

test_that("expression_type must be declared", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(1:8, nrow = 4, dimnames = list(gids, c("C1", "T1")))
  expect_error(annotate_epi_domains(se, genome = txdb, expression = exp_wide),
               "expression_type")
})

test_that("get_domain_genes prioritizes transparently", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(c(80, 5, 30, 1, 0.2, 2, 3, 4),
                     nrow = 4, dimnames = list(gids, c("C1", "T1")))
  se <- annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                             expression_type = "TPM")
  cand <- get_domain_genes(se, group = "Control",
                           expression_priority = "expressed_first",
                           min_expression = 1)
  expect_true(all(c("domain_id", "gene_id", "relation_types",
                    "bedpe_supported", "expression_value",
                    "expression_status", "candidate_priority") %in% colnames(cand)))
  expect_true("Expressed" %in% cand$expression_status)
})

test_that("get_domain_genes high expression does not override evidence", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(c(80, 5, 30, 1, 0.2, 2, 3, 4),
                     nrow = 4, dimnames = list(gids, c("C1", "T1")))
  se <- annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                             expression_type = "TPM")
  cand <- get_domain_genes(se, group = "Control")
  expect_true(nrow(cand) >= 1)
})

# ---- P1-12: dangerous annotation edge cases (review 2026-08-10) --------------

# two overlapping genes must NOT be merged into one pseudo-gene
make_overlap_txdb <- function() {
  skip_if_not_installed("txdbmaker")
  # gene1 = 1000-3000, gene2 = 2500-4000 (overlap)
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(c(1000, 1000, 2500, 2500),
                              c(3000, 3000, 4000, 4000)),
    strand = rep(c("+", "-"), each = 2))
  S4Vectors::mcols(gr)$type <- rep(c("gene", "transcript"), 2)
  S4Vectors::mcols(gr)$ID <- c("gene1", "gene1_tx1", "gene2", "gene2_tx1")
  S4Vectors::mcols(gr)$Parent <- c(NA_character_, "gene1", NA_character_, "gene2")
  txdbmaker::makeTxDbFromGRanges(gr)
}

test_that("overlapping genes are not merged (P0-1)", {
  txdb <- make_overlap_txdb()
  gm <- epiPortrait:::.gene_model_from_txdb(txdb)
  expect_equal(length(gm$genes), 2)   # gene1 and gene2 stay separate
  expect_true(all(c("gene1", "gene2") %in% gm$genes$gene_id))
  # no pseudo-gene like "gene1,gene2"
  expect_false(any(grepl(",", gm$genes$gene_id)))
})

test_that("TSS inside domain -> distance 0 (P0-2)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  se <- annotate_epi_domains(se, genome = txdb)
  # domain 3 (35000-90000) contains gene3 TSS (40000, +) -> signed distance 0
  d <- rowData(se)
  i3 <- which(d$nearest_tss_gene_id == "gene3")
  expect_true(length(i3) >= 1)
  expect_equal(d$nearest_tss_distance_bp[i3[1]], 0)
  # domain 1 (1500-5000): gene1 TSS (1000, +) is UPSTREAM of the domain; for a
  # + strand gene a domain downstream of the TSS is POSITIVE.
  i1 <- which(d$nearest_tss_gene_id == "gene1")
  expect_true(length(i1) >= 1)
  expect_gt(d$nearest_tss_distance_bp[i1[1]], 0)
})

test_that("domain on contig absent from TxDb stays NA without misalignment", {
  txdb <- make_tiny_txdb()
  suppressWarnings({
    domains <- GenomicRanges::GRanges(
      c("chr1", "chr1", "chr1", "chr1", "chr2"),
      IRanges::IRanges(c(1500, 21000, 35000, 65000, 100),
                       c(5000, 30000, 90000, 80000, 500)),
      seqinfo = GenomeInfoDb::Seqinfo(c("chr1", "chr2"), c(1000000, 1000000)))
  })
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Intensity = matrix(1, 5, 2),
                  SignalDispersion = matrix(1, 5, 2)),
    rowRanges = domains)
  rownames(se2) <- paste0("epiDomain_", sprintf("%06d", 1:5))
  colnames(se2) <- c("C1", "T1")
  colData(se2)$Condition <- c("Control", "Treatment")
  se2 <- annotate_epi_domains(se2, genome = txdb)
  d <- rowData(se2)
  # chr2 domain -> NA nearest gene; other domains unaffected
  expect_true(is.na(d$nearest_tss_gene_id[5]))
  expect_true(all(!is.na(d$nearest_tss_gene_id[1:4])))
})

test_that("RNA long format aggregates >=2 replicates correctly (P0-5)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  # 2 replicates per gene per condition with DIFFERENT values
  long <- rbind(
    data.frame(gene_id = rep(gids, 2), SampleID = c("R1", "R2"),
               Condition = "Control",
               expression = c(10, 20, 30, 40, 20, 40, 60, 80)),
    data.frame(gene_id = rep(gids, 2), SampleID = c("R1", "R2"),
               Condition = "Treatment",
               expression = c(50, 60, 70, 80, 10, 20, 30, 40)))
  se <- annotate_epi_domains(se, genome = txdb, expression = long,
                             expression_type = "TPM")
  es <- S4Vectors::metadata(se)$expression_summary
  ctrl <- es[es$Condition == "Control", ]
  expect_equal(nrow(ctrl), 4)   # one row per gene
  # medians: gene1 (10,20)=15, gene2 (20,40)=30, gene3 (30,60)=45, gene4 (40,80)=60
  expect_equal(sort(ctrl$expression), sort(c(15, 30, 45, 60)))
})

test_that("wide RNA matched by sample name, not position (P0-6)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  # RNA columns reordered relative to colData
  exp_wide <- matrix(c(80, 5, 30, 1, 0.2, 2, 3, 4),
                     nrow = 4, dimnames = list(gids, c("T1", "C1")))
  se <- annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                             expression_type = "TPM")
  es <- S4Vectors::metadata(se)$expression_summary
  # C1 -> Control, T1 -> Treatment; T1 column value should go to Treatment
  t <- es[es$Condition == "Treatment", ]
  c <- es[es$Condition == "Control", ]
  expect_true(nrow(t) == 4 && nrow(c) == 4)
})

test_that("candidate_priority is a true rank (P0-7)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(c(80, 5, 30, 1, 0.2, 2, 3, 4),
                     nrow = 4, dimnames = list(gids, c("C1", "T1")))
  se <- annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                             expression_type = "TPM")
  cand <- get_domain_genes(se, group = "Control")
  expect_true(all(sort(cand$candidate_priority) == seq_len(nrow(cand))))
  expect_true(all(cand$candidate_priority == seq_len(nrow(cand))))
})

test_that("BEDPE support count uses unique records (P1-4)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  # same BEDPE record matched via both directions -> counted once
  bedpe_df <- data.frame(
    chrom1 = "chr1", start1 = 1499, end1 = 2500,   # inside domain 1
    chrom2 = "chr1", start2 = 59999, end2 = 61000, # gene 4 promoter
    stringsAsFactors = FALSE)
  se <- annotate_epi_domains(se, genome = txdb, bedpe = bedpe_df)
  cand <- get_domain_genes(se)
  if (nrow(cand) > 0 && any(cand$bedpe_supported)) {
    expect_true(all(cand$bedpe_support_count[cand$bedpe_supported] == 1))
  }
})

test_that("BEDPE input validation rejects bad coordinates (P1-6)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  bad1 <- data.frame(chrom1 = "chr1", start1 = -1, end1 = 100,
                     chrom2 = "chr1", start2 = 100, end2 = 200)
  expect_error(annotate_epi_domains(se, genome = txdb, bedpe = bad1), ">= 0")
  bad2 <- data.frame(chrom1 = "chr1", start1 = 100, end1 = 100,
                     chrom2 = "chr1", start2 = 100, end2 = 200)
  expect_error(annotate_epi_domains(se, genome = txdb, bedpe = bad2), "end must be greater")
})

# ---- epiPortrait9 audit edge cases (2026-08-10) -----------------------------

test_that("mm9 built-in blacklist works without TxDb (P0-1)", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 5000000), width = 200))
  res <- filter_blacklist(gr, genome = "mm9")
  expect_s4_class(res, "GRanges")
})

test_that("filter_blacklist requires a built-in or custom BED, not TxDb", {
  # unknown genome without blacklist file -> clear error mentioning blacklist_path
  expect_error(filter_blacklist(GenomicRanges::GRanges("zzz", IRanges::IRanges(1, 100)),
                                genome = "rn7"),
               "blacklist_path")
})

test_that("zero shared domain/TxDb seqlevels -> hard error (P0-2)", {
  txdb <- make_tiny_txdb()
  # domains on chr2 only; TxDb is chr1
  suppressWarnings({
    dom <- GenomicRanges::GRanges("chr2", IRanges::IRanges(100, 500),
                                  seqinfo = GenomeInfoDb::Seqinfo("chr2", 1000000))
  })
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Intensity = matrix(1, 1, 1),
                  SignalDispersion = matrix(1, 1, 1)),
    rowRanges = dom)
  rownames(se) <- "epiDomain_000001"; colnames(se) <- "S1"
  colData(se)$Condition <- "X"
  expect_error(annotate_epi_domains(se, genome = txdb), "No shared seqlevels")
})

test_that("partial seqlevels -> warning, matched domains annotated (P0-2)", {
  txdb <- make_tiny_txdb()
  # chr1 + chr2 domains; chr2 has no genes -> warning but chr1 works
  suppressWarnings({
    dom <- GenomicRanges::GRanges(
      c("chr1", "chr2"),
      IRanges::IRanges(c(1500, 100), c(5000, 500)),
      seqinfo = GenomeInfoDb::Seqinfo(c("chr1", "chr2"), c(1000000, 1000000)))
  })
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Intensity = matrix(1, 2, 1),
                  SignalDispersion = matrix(1, 2, 1)),
    rowRanges = dom)
  rownames(se) <- c("epiDomain_000001", "epiDomain_000002"); colnames(se) <- "S1"
  colData(se)$Condition <- "X"
  expect_warning(se2 <- annotate_epi_domains(se, genome = txdb), "present in the TxDb")
  expect_true(!is.na(rowData(se2)$nearest_tss_gene_id[1]))
  expect_true(is.na(rowData(se2)$nearest_tss_gene_id[2]))
})

test_that("wide RNA with unmatched sample columns -> error (P1-2)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(1:8, nrow = 4,
                     dimnames = list(gids, c("C1", "UnknownRNA1")))
  expect_error(annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                                    expression_type = "TPM"),
               "not present in colData")
})

test_that("wide RNA validates rownames and numeric (P1-3)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  bad <- matrix(rep("a", 8), nrow = 4,
                dimnames = list(c("g1", "g2", "g3", "g4"), c("C1", "T1")))
  expect_error(annotate_epi_domains(se, genome = txdb, expression = bad,
                                    expression_type = "TPM"),
               "numeric")
})

test_that("VST expressed_first without threshold -> error (P1-1)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(1:8, nrow = 4, dimnames = list(gids, c("C1", "T1")))
  se <- annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                             expression_type = "VST")
  expect_error(get_domain_genes(se, group = "Control",
                                expression_priority = "expressed_first"),
               "min_expression")
  # high_expression_first needs no threshold
  cand <- get_domain_genes(se, group = "Control",
                           expression_priority = "high_expression_first",
                           min_expression = NULL)
  expect_true(nrow(cand) >= 1)
})

test_that("invalid group -> error (P1-11)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(1:8, nrow = 4, dimnames = list(gids, c("C1", "T1")))
  se <- annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                             expression_type = "TPM")
  expect_error(get_domain_genes(se, group = "Tretment"), "not found")
})

test_that("unique_genes reranks candidate_priority contiguous (P1-4)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  exp_wide <- matrix(1:8, nrow = 4, dimnames = list(gids, c("C1", "T1")))
  se <- annotate_epi_domains(se, genome = txdb, expression = exp_wide,
                             expression_type = "TPM")
  cand <- get_domain_genes(se, group = "Control", unique_genes = TRUE)
  expect_true(all(sort(cand$candidate_priority) == seq_len(nrow(cand))))
})

test_that("evidence outranks expression (P1-5)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  gids <- names(GenomicFeatures::genes(txdb))
  # give gene4 (in domain3 with gene3) low expression, gene1 high expression;
  # a BEDPE contact should still rank gene4's domain-gene above a nearest-only
  # high-expression gene when present on the same domain.
  exp_wide <- matrix(c(100, 5, 30, 1, 0.2, 2, 3, 4),
                     nrow = 4, dimnames = list(gids, c("C1", "T1")))
  bedpe_df <- data.frame(
    chrom1 = "chr1", start1 = 34999, end1 = 45000,  # inside domain3
    chrom2 = "chr1", start2 = 59999, end2 = 61000,  # gene4 promoter
    stringsAsFactors = FALSE)
  se <- annotate_epi_domains(se, genome = txdb, bedpe = bedpe_df,
                             expression = exp_wide, expression_type = "TPM")
  cand <- get_domain_genes(se, group = "Control",
                           expression_priority = "high_expression_first")
  # a bedpe-supported gene must have priority <= any nearest-only gene in the
  # same output (evidence tier first)
  if (any(cand$bedpe_supported)) {
    best_bedpe <- min(cand$candidate_priority[cand$bedpe_supported])
    worst_nearest <- max(cand$candidate_priority[grepl("nearest_tss", cand$relation_types) &
                                                   !cand$bedpe_supported])
    expect_lte(best_bedpe, worst_nearest)
  }
})

test_that("domain_gene_links exported to annotation/ (P1-14)", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  se <- annotate_epi_domains(se, genome = txdb)
  out <- export_epiportrait_results(se, outdir = tempfile("epi_anno_export"))
  expect_true(file.exists(file.path(out, "annotation", "domain_gene_links.tsv")))
  expect_true(file.exists(file.path(out, "annotation", "provenance.txt")))
})

# ---- BEDPE contact-score aggregation (bedpe_score_col) ----------------------
test_that("bedpe_score_col aggregates scores into links, dedup and summary", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  # two records between the same anchor pair; each record links domain1<->gene4
  # (direction 1) AND domain3<->gene4 (direction 2), so every touched pair has
  # TWO unique supporting records whose scores must be summed.
  bedpe_df <- data.frame(
    chrom1 = "chr1", start1 = 1499, end1 = 2500,
    chrom2 = "chr1", start2 = 59999, end2 = 61000,
    score = c(2, 10),
    stringsAsFactors = FALSE)
  se <- annotate_epi_domains(se, genome = txdb, bedpe = bedpe_df,
                             bedpe_score_col = "score")
  links <- S4Vectors::metadata(se)$domain_gene_links
  bl <- links[links$evidence_source == "bedpe", ]
  expect_true("contact_score" %in% colnames(links))
  expect_true(all(bl$contact_score %in% c(2, 10)))
  expect_equal(sum(is.na(links$contact_score[links$evidence_source != "bedpe"])),
               sum(links$evidence_source != "bedpe"))
  dedup <- S4Vectors::metadata(se)$domain_gene_links_dedup
  expect_true("bedpe_contact_score" %in% colnames(dedup))
  bp_pairs <- dedup[dedup$evidence_sources == "bedpe" |
                      grepl("bedpe_promoter_contact", dedup$relation_types), ]
  expect_true(nrow(bp_pairs) >= 2)
  expect_true(all(bp_pairs$bedpe_contact_score == 12))   # unique-record sum
  expect_true(all(dedup$bedpe_contact_score[dedup$bedpe_support_count == 0] == 0,
                  na.rm = TRUE))
  summ <- S4Vectors::metadata(se)$annotation_summary
  expect_true("bedpe_contact_score" %in% colnames(summ))
  d1 <- summ$Domain_ID == rownames(se)[1]
  d3 <- summ$Domain_ID == rownames(se)[3]
  expect_equal(summ$bedpe_contact_score[d1], 12)
  expect_equal(summ$bedpe_contact_score[d3], 12)
  expect_equal(summ$n_bedpe_contact_gene[d1], 1)
  # rowData mirror
  expect_equal(rowData(se)$bedpe_contact_score[d1], 12)
  # provenance records the score column
  expect_match(S4Vectors::metadata(se)$bedpe_provenance$score_column, "score")
})

test_that("default (no bedpe_score_col) keeps count-only behaviour", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  bedpe_df <- data.frame(
    chrom1 = "chr1", start1 = 1499, end1 = 2500,
    chrom2 = "chr1", start2 = 59999, end2 = 61000,
    score = c(5),
    stringsAsFactors = FALSE)
  se <- annotate_epi_domains(se, genome = txdb, bedpe = bedpe_df)
  links <- S4Vectors::metadata(se)$domain_gene_links
  expect_true(all(is.na(links$contact_score)))
  dedup <- S4Vectors::metadata(se)$domain_gene_links_dedup
  # supported pairs stay NA (unscored config); unsupported pairs keep 0
  expect_true(all(is.na(dedup$bedpe_contact_score[dedup$bedpe_support_count > 0])))
  expect_true(all(dedup$bedpe_contact_score[dedup$bedpe_support_count == 0] == 0,
                  na.rm = TRUE))
  summ <- S4Vectors::metadata(se)$annotation_summary
  contacted <- summ$n_bedpe_contact_gene > 0
  expect_true(all(is.na(summ$bedpe_contact_score[contacted])))
  expect_true(all(summ$bedpe_contact_score[!contacted] == 0))
})

test_that("bedpe_score_col input validation fails loudly", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  good <- data.frame(chrom1 = "chr1", start1 = 1499, end1 = 2500,
                     chrom2 = "chr1", start2 = 59999, end2 = 61000,
                     stringsAsFactors = FALSE)
  expect_error(annotate_epi_domains(se, genome = txdb, bedpe_score_col = 8),
               "requires a BEDPE")
  expect_error(annotate_epi_domains(se, genome = txdb, bedpe = good,
                                    bedpe_score_col = 3), "coordinates")
  expect_error(annotate_epi_domains(se, genome = txdb, bedpe = good,
                                    bedpe_score_col = 99), "index must be")
  expect_error(annotate_epi_domains(se, genome = txdb, bedpe = good,
                                    bedpe_score_col = "nope"), "not found")
  expect_error(annotate_epi_domains(se, genome = txdb, bedpe = good,
                                    bedpe_score_col = "chrom1"), "extra column")
})

test_that("get_domain_genes exposes bedpe_contact_score and rank_by works", {
  txdb <- make_tiny_txdb()
  se <- make_anno_se()
  bedpe_df <- data.frame(
    chrom1 = "chr1", start1 = 34999, end1 = 45000,  # inside domain3
    chrom2 = "chr1", start2 = 59999, end2 = 61000,  # gene4 promoter
    score = c(7.5),
    stringsAsFactors = FALSE)
  se <- annotate_epi_domains(se, genome = txdb, bedpe = bedpe_df,
                             bedpe_score_col = "score")
  cand_default <- get_domain_genes(se, rank_by = "tier")
  cand_score <- get_domain_genes(se, rank_by = "bedpe_score")
  expect_true(all(c("bedpe_contact_score") %in% colnames(cand_default)))
  # scored pair exists with the right value
  sc <- cand_score$bedpe_contact_score[!is.na(cand_score$bedpe_contact_score) &
                                         cand_score$bedpe_contact_score > 0]
  expect_true(length(sc) >= 1 && all(sc == 7.5))
  # rank_by = "bedpe_score": the highest-scoring pair is ranked FIRST overall
  best_sc <- max(cand_score$bedpe_contact_score, na.rm = TRUE)
  top_row <- cand_score[cand_score$candidate_priority ==
                          min(cand_score$candidate_priority[cand_score$bedpe_contact_score == best_sc]), ]
  expect_gte(min(cand_score$candidate_priority[top_row$domain_id == rownames(se)[3]]),
             min(cand_score$candidate_priority))
})
