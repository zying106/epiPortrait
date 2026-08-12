# Domain-aware annotation and candidate-gene prioritization.
#
# Design (2026-08-10): epiPortrait's analysis object is the shared FINAL
# domain (H3K4me3 broad domains, H3K27ac stitched SE-like domains, H3K27me3 /
# H3K9me3 repressive domains), which can be tens to hundreds of kb. Traditional
# peak-level annotation (ChIPseeker nearest-gene) is information-poor for these.
# This module provides genome-flexible, domain-aware annotation with optional
# BEDPE 3D-contact and RNA-seq expression evidence, plus transparent
# candidate-gene prioritization. It does NOT call peaks, loops, or fit
# differential expression models.

# Resolve genome resources (TxDb / OrgDb / blacklist) from a genome shortcut or
# user-supplied objects (design §8). annotate_epi_domains() and
# filter_promoter_peaks() share the TxDb resolution via this resolver;
# filter_blacklist() intentionally resolves blacklists independently of the
# TxDb/OrgDb.
.resolve_genome_resources <- function(genome = "hg38", txdb = NULL,
                                      anno_db = NULL, blacklist = NULL) {
  if (is(genome, "TxDb")) {
    # infer seqstyle from the TxDb seqlevels (may be UCSC 'chr1' or '1')
    sl <- GenomeInfoDb::seqlevels(genome)
    seqstyle <- if (length(sl) > 0 && all(grepl("^chr", sl))) "UCSC" else
      if (length(sl) > 0) "custom" else "unknown"
    return(list(genome_name = "custom", genome_class = "custom",
                organism = NA_character_, txdb = genome,
                anno_db = anno_db, blacklist = blacklist,
                seqstyle = seqstyle,
                resource_source = "user-supplied TxDb"))
  }
  if (!is.character(genome)) {
    stop("genome must be a character shortcut ('hg38','hg19','mm10','mm9') or a TxDb.", call. = FALSE)
  }
  cfg <- switch(genome,
    hg38 = list(txdb_pkg = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                txdb_obj = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                anno_pkg = "org.Hs.eg.db",
                organism = "Homo sapiens"),
    hg19 = list(txdb_pkg = "TxDb.Hsapiens.UCSC.hg19.knownGene",
                txdb_obj = "TxDb.Hsapiens.UCSC.hg19.knownGene",
                anno_pkg = "org.Hs.eg.db",
                organism = "Homo sapiens"),
    mm10 = list(txdb_pkg = "TxDb.Mmusculus.UCSC.mm10.knownGene",
                txdb_obj = "TxDb.Mmusculus.UCSC.mm10.knownGene",
                anno_pkg = "org.Mm.eg.db",
                organism = "Mus musculus"),
    NULL)
  if (is.null(cfg)) {
    # mm9 / other: partial support (blacklist only) unless TxDb supplied.
    if (!is.null(txdb)) {
      return(list(genome_name = genome, genome_class = "custom",
                  organism = NA_character_, txdb = txdb, anno_db = anno_db,
                  blacklist = blacklist, seqstyle = "custom",
                  resource_source = "custom"))
    }
    warning(sprintf(
      "No built-in TxDb/OrgDb shortcut for genome '%s'. Provide txdb / anno_db for annotation; core quantitative analysis is genome-agnostic.",
      genome), call. = FALSE)
    return(list(genome_name = genome, genome_class = "partial",
                organism = NA_character_, txdb = NULL, anno_db = anno_db,
                blacklist = blacklist, seqstyle = "unknown",
                resource_source = "none"))
  }
  # Full built-in support: resolve TxDb and OrgDb lazily.
  tx <- txdb
  if (is.null(tx)) {
    if (!requireNamespace(cfg$txdb_pkg, quietly = TRUE))
      stop(sprintf("Install '%s' for built-in genome '%s'.", cfg$txdb_pkg, genome))
    tx <- get(cfg$txdb_obj, envir = asNamespace(cfg$txdb_pkg))
  }
  adb <- anno_db
  if (is.null(adb)) adb <- cfg$anno_pkg
  bl <- blacklist
  if (is.null(bl)) {
    bl_path <- system.file("extdata", paste0(genome, "-blacklist.v2.bed"),
                           package = "epiPortrait")
    if (file.exists(bl_path)) bl <- rtracklayer::import(bl_path)
  }
  list(genome_name = genome, genome_class = "builtin",
       organism = cfg$organism, txdb = tx, anno_db = adb, blacklist = bl,
       seqstyle = "UCSC",
       resource_source = sprintf("built-in (%s)", cfg$txdb_pkg))
}


# Check seqlevel compatibility across BigWig / domains / TxDb / BEDPE without
# silently converting 1<->chr1 or MT<->chrM (design §10). When enforce = TRUE
# (used by annotate_epi_domains), zero shared domain/TxDb seqlevels is a hard
# error and partial overlap is a warning, so a genome-build / naming mismatch
# is never silently reinterpreted as "intergenic domains" (P0-2).
.check_seqlevel_compatibility <- function(domains, txdb = NULL, bedpe = NULL,
                                          bw_seqlevels = NULL, enforce = FALSE) {
  out <- list(shared = NULL, signal_only = NULL, domain_only = NULL,
              txdb_only = NULL, bedpe_only = NULL,
              domain_without_txdb = NULL)
  dom_sl <- as.character(GenomicRanges::seqnames(domains))
  dom_set <- unique(dom_sl)
  out$shared <- dom_set

  if (!is.null(bw_seqlevels)) {
    bw_set <- unique(as.character(bw_seqlevels))
    out$shared <- intersect(out$shared, bw_set)
    out$signal_only <- setdiff(bw_set, dom_set)
    out$domain_only <- setdiff(dom_set, bw_set)
  }
  if (!is.null(txdb)) {
    tx_sl <- as.character(GenomeInfoDb::seqlevels(txdb))
    out$txdb_only <- setdiff(tx_sl, dom_set)
    out$domain_without_txdb <- setdiff(dom_set, tx_sl)
    out$shared <- intersect(out$shared, tx_sl)
    if (enforce) {
      if (length(out$shared) == 0) {
        stop("No shared seqlevels between the domains and the TxDb. ",
             "Check the genome build and chromosome naming convention ",
             "(e.g. 'chr1' vs '1').",
             if (length(out$domain_without_txdb) > 0)
               paste0(" Domain-only seqlevels: ",
                      paste(head(out$domain_without_txdb, 10), collapse = ", ")) else "")
      }
      if (length(out$domain_without_txdb) > 0) {
        warning(sprintf(
          "Only %d/%d domain seqlevels are present in the TxDb; %d domain(s) will have NA annotation (domain-only: %s).",
          length(out$shared), length(dom_set), length(out$domain_without_txdb),
          paste(head(out$domain_without_txdb, 10), collapse = ", ")),
          call. = FALSE)
      }
    }
  }
  if (!is.null(bedpe) && is(bedpe, "GRanges")) {
    bp_sl <- unique(as.character(GenomicRanges::seqnames(bedpe)))
    out$bedpe_only <- setdiff(bp_sl, dom_set)
  }
  out
}


# Build a gene-level GRanges model from a TxDb: one range per gene with a
# strand-aware TSS, gene body, and a per-gene promoter window (design §13-20).
# Transcripts are collapsed to gene level (multi-transcript -> gene).
.gene_model_from_txdb <- function(txdb, promoter_upstream = 3000,
                                  promoter_downstream = 3000) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE))
    stop("Please install 'GenomicFeatures'.")
  # one TxDb gene = one gene-level feature. NEVER reduce() across genes: two
  # overlapping genes must remain two separate features (P0-1).
  gr <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
  if (length(gr) == 0) {
    stop("TxDb contains no genes.")
  }
  gene_ids <- names(gr)
  if (is.null(gene_ids) || any(gene_ids == "")) {
    gene_ids <- as.character(S4Vectors::mcols(gr)$gene_id)
  }
  if (is.null(gene_ids) || any(is.na(gene_ids))) {
    gene_ids <- as.character(seq_along(gr))
  }
  names(gr) <- gene_ids
  S4Vectors::mcols(gr)$gene_id <- gene_ids
  # strand-aware per-gene TSS: start for +, end for -; midpoint if unstranded
  gene_strand <- as.character(GenomicRanges::strand(gr))
  gene_tss <- ifelse(gene_strand == "+", GenomicRanges::start(gr),
              ifelse(gene_strand == "-", GenomicRanges::end(gr),
                     round((GenomicRanges::start(gr) + GenomicRanges::end(gr)) / 2)))
  gene_tss <- as.integer(gene_tss)
  S4Vectors::mcols(gr)$tss <- gene_tss
  S4Vectors::mcols(gr)$gene_strand <- gene_strand
  # promoter window: TSS +/- up/downstream, strand-aware
  prom_start <- ifelse(gene_strand == "-", gene_tss - promoter_downstream,
                       gene_tss - promoter_upstream)
  prom_end <- ifelse(gene_strand == "-", gene_tss + promoter_upstream,
                     gene_tss + promoter_downstream)
  promoters_gr <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr),
    ranges = IRanges::IRanges(pmin(prom_start, prom_end),
                              pmax(prom_start, prom_end)),
    gene_id = gene_ids, tss = gene_tss)
  list(genes = gr, promoters = promoters_gr, tss = gene_tss,
       gene_id_type = "TxDb gene_id")
}


#' Domain-Aware Annotation of epiPortrait Domains
#'
#' @description Annotates the final shared domains (not individual peaks) with
#' genome-aware linear evidence plus optional BEDPE 3D-contact and RNA-seq
#' expression evidence. Adds 7 compact rowData columns and stores a long-format
#' domain-gene link table in \code{metadata(se)$domain_gene_links}.
#'
#' @details
#' **Representative gene-level TSS.** epiPortrait uses one strand-aware
#' 5' boundary per gene derived from \code{GenomicFeatures::genes(txdb)} (start
#' for +, end for -) as the representative TSS for compact domain annotation. It
#' does not enumerate every alternative-transcript TSS of a gene; the compact
#' \code{nearest_tss_distance_bp} is therefore a gene-level proximity measure.
#'
#' **Promoter models.** \code{filter_promoter_peaks()} uses
#' \code{GenomicFeatures::promoters(txdb)} (transcript-level promoter set), while
#' \code{annotate_epi_domains()} uses a gene-level representative-TSS promoter
#' window. The two are not the same promoter universe; in an enhancer workflow,
#' promoter exclusion and \code{promoter_overlap_gene_count} should not be
#' assumed to use identical definitions.
#'
#' @param se A SummarizedExperiment from \code{build_portrait_matrix()}.
#' @param genome Character. Built-in shortcut ("hg38","hg19","mm10","mm9") or a
#'   TxDb object. For non-built-in species supply \code{txdb}.
#' @param txdb A TxDb object or NULL (resolved from genome shortcut).
#' @param anno_db Character. OrgDb package name, or NULL (resolved).
#' @param promoter_upstream Numeric. Promoter window upstream of the TSS
#'   (default 3000 bp).
#' @param promoter_downstream Numeric. Promoter window downstream of the TSS
#'   (default 3000 bp).
#' @param bedpe Character path to a BEDPE file, a data.frame, or NULL. Only
#'   \code{bedpe_promoter_contact} evidence is added (3D promoter-contact);
#'   records are read-only (no calling / filtering / merging).
#' @param expression A gene-level expression object: a wide matrix (genes x
#'   samples), a long data.frame (gene_id, SampleID, Condition, expression), or
#'   NULL.
#' @param expression_type Character. "TPM", "CPM", "normalized_counts", "VST",
#'   "rlog". Must be declared; never auto-guessed.
#' @param aggregate_fun Character. Per-condition aggregation over RNA replicates
#'   (default "median").
#' @param gene_id_keytype Character or NULL. OrgDb keytype for the gene IDs
#'   (default "ENTREZID" for built-in UCSC genomes; required for custom TxDbs
#'   using other ID systems such as ENSEMBL / TAIR / FlyBase).
#' @return \code{se} with rowData columns: primary_genomic_context,
#'   nearest_tss_gene_id, nearest_tss_gene_symbol, nearest_tss_distance_bp,
#'   promoter_overlap_gene_count, gene_body_overlap_gene_count,
#'   fully_contained_gene_count; metadata(se)$domain_gene_links,
#'   metadata(se)$annotation_provenance, and (if provided)
#'   metadata(se)$bedpe_provenance / metadata(se)$expression_provenance.
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   se <- annotate_epi_domains(example_se, genome = "hg38")
#'   head(SummarizedExperiment::rowData(se)$primary_genomic_context)
#' }
#' @export
annotate_epi_domains <- function(se, genome = "hg38", txdb = NULL, anno_db = NULL,
                                 promoter_upstream = 3000, promoter_downstream = 3000,
                                 bedpe = NULL, expression = NULL,
                                 expression_type = NULL,
                                 aggregate_fun = "median",
                                 gene_id_keytype = NULL) {
  # P1-10: validate promoter window parameters (negative values would be
  # "corrected" by pmin/pmax but are biologically meaningless).
  for (nm in c("promoter_upstream", "promoter_downstream")) {
    v <- get(nm)
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v < 0) {
      stop(sprintf("%s must be a finite non-negative number.", nm))
    }
  }
  # ---- genome resources ----------------------------------------------------
  res <- .resolve_genome_resources(genome, txdb = txdb, anno_db = anno_db)
  if (is.null(res$txdb)) {
    stop("Annotation requires a TxDb. Provide genome = 'hg38'/'hg19'/'mm10' or a TxDb object.")
  }
  gm <- .gene_model_from_txdb(res$txdb,
                              promoter_upstream = promoter_upstream,
                              promoter_downstream = promoter_downstream)
  genes_gr <- gm$genes
  promoters_gr <- gm$promoters
  gene_ids <- gm$gene_id_type

  domains <- rowRanges(se)
  domain_ids <- rownames(se)
  n_dom <- length(domains)

  # ---- seqlevel compatibility (P0-2: enforced, never silently ignored) ------
  sl_check <- .check_seqlevel_compatibility(domains, txdb = res$txdb,
                                            enforce = TRUE)

  # ---- linear evidence ------------------------------------------------------
  # nearest TSS: absolute shortest interval-to-point distance, backfilled by
  # queryHits() so domains without a hit (e.g. on contigs absent from the TxDb)
  # remain NA instead of misaligning (P0-2, P0-3).
  tss_gr <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(genes_gr),
    ranges = IRanges::IRanges(mcols(genes_gr)$tss, width = 1))
  d2t <- GenomicRanges::distanceToNearest(domains, tss_gr)
  qh_t <- S4Vectors::queryHits(d2t)
  sh_t <- S4Vectors::subjectHits(d2t)
  nearest_gene_id <- rep(NA_character_, n_dom)
  nearest_dist <- rep(NA_real_, n_dom)
  nearest_gene_id[qh_t] <- mcols(genes_gr)$gene_id[sh_t]
  nearest_dist[qh_t] <- mcols(d2t)$distance
  # Strand-aware SIGNED shortest distance: domain interval to TSS. When the
  # TSS lies inside the domain, distance is 0 (NOT the midpoint offset, P0-2).
  # Sign: + strand -> positive downstream of TSS, negative upstream; - strand
  # is mirrored.
  signed_dist <- rep(NA_real_, n_dom)
  for (i in qh_t) {
    g <- genes_gr[mcols(genes_gr)$gene_id == nearest_gene_id[i]]
    if (length(g) == 0) next
    tss_i <- mcols(g)$tss[1]
    s <- mcols(g)$gene_strand[1]
    # shortest absolute distance already from distanceToNearest; recompute sign
    d <- nearest_dist[i]
    if (s == "+") {
      signed_dist[i] <- if (end(domains[i]) < tss_i) -d else d
    } else if (s == "-") {
      signed_dist[i] <- if (start(domains[i]) > tss_i) -d else d
    } else {
      signed_dist[i] <- d
    }
  }

  # promoter overlap / gene body overlap / fully contained
  prom_ov <- GenomicRanges::findOverlaps(domains, promoters_gr)
  body_ov <- GenomicRanges::findOverlaps(domains, genes_gr, ignore.strand = TRUE)
  prom_count <- tabulate(S4Vectors::queryHits(prom_ov), nbins = n_dom)
  body_count <- tabulate(S4Vectors::queryHits(body_ov), nbins = n_dom)
  # fully contained: gene body entirely inside domain
  contained <- GenomicRanges::findOverlaps(genes_gr, domains, ignore.strand = TRUE,
                                           type = "within")
  contained_count <- tabulate(S4Vectors::subjectHits(contained), nbins = n_dom)

  # primary genomic context (P0-4): simplified to 3 values consistent with the
  # count columns. Detailed exon/intron/downstream relations live in
  # domain_gene_links; they are not re-encoded here.
  ctx <- rep("Intergenic", n_dom)
  prom_dom <- unique(S4Vectors::queryHits(prom_ov))
  body_dom <- unique(S4Vectors::queryHits(body_ov))
  ctx[prom_dom] <- "Promoter-associated"
  ctx[setdiff(body_dom, prom_dom)] <- "Gene-body-associated"

  # ---- rowData columns ------------------------------------------------------
  # P1-2: built-in UCSC knownGene uses ENTREZID; custom TxDbs may use other ID
  # systems, so require an explicit gene_id_keytype for custom genomes.
  if (is.null(gene_id_keytype) && res$genome_class == "builtin") {
    gene_id_keytype <- "ENTREZID"
  }
  rowData(se)$primary_genomic_context <- ctx
  rowData(se)$nearest_tss_gene_id <- nearest_gene_id
  rowData(se)$nearest_tss_distance_bp <- signed_dist
  rowData(se)$promoter_overlap_gene_count <- prom_count
  rowData(se)$gene_body_overlap_gene_count <- body_count
  rowData(se)$fully_contained_gene_count <- contained_count

  # ---- long-format domain-gene links (design §22-34) ------------------------
  links <- data.frame()
  # nearest_tss rows
  if (n_dom > 0 && any(!is.na(nearest_gene_id))) {
    links <- rbind(links, data.frame(
      domain_id = domain_ids[!is.na(nearest_gene_id)],
      gene_id = nearest_gene_id[!is.na(nearest_gene_id)],
      gene_symbol = NA_character_,   # backfilled below (P1-1)
      relation_type = "nearest_tss",
      distance_to_tss_bp = signed_dist[!is.na(nearest_gene_id)],
      overlap_bp = NA_real_,
      domain_overlap_fraction = NA_real_,
      feature_overlap_fraction = NA_real_,
      bedpe_record_id = NA_character_,
      evidence_source = "linear",
      stringsAsFactors = FALSE))
  }
  # promoter_overlap rows
  if (length(prom_ov) > 0) {
    qh <- S4Vectors::queryHits(prom_ov); sh <- S4Vectors::subjectHits(prom_ov)
    ov <- GenomicRanges::pintersect(domains[qh], promoters_gr[sh])
    links <- rbind(links, data.frame(
      domain_id = domain_ids[qh],
      gene_id = mcols(promoters_gr)$gene_id[sh],
      gene_symbol = NA_character_,
      relation_type = "promoter_overlap",
      distance_to_tss_bp = NA_real_,
      overlap_bp = GenomicRanges::width(ov),
      domain_overlap_fraction = GenomicRanges::width(ov) /
        GenomicRanges::width(domains[qh]),
      feature_overlap_fraction = GenomicRanges::width(ov) /
        GenomicRanges::width(promoters_gr[sh]),
      bedpe_record_id = NA_character_,
      evidence_source = "linear",
      stringsAsFactors = FALSE))
  }
  # gene_body_overlap rows
  if (length(body_ov) > 0) {
    qh <- S4Vectors::queryHits(body_ov); sh <- S4Vectors::subjectHits(body_ov)
    ov <- GenomicRanges::pintersect(domains[qh], genes_gr[sh], ignore.strand = TRUE)
    links <- rbind(links, data.frame(
      domain_id = domain_ids[qh],
      gene_id = mcols(genes_gr)$gene_id[sh],
      gene_symbol = NA_character_,
      relation_type = "gene_body_overlap",
      distance_to_tss_bp = NA_real_,
      overlap_bp = GenomicRanges::width(ov),
      domain_overlap_fraction = GenomicRanges::width(ov) /
        GenomicRanges::width(domains[qh]),
      feature_overlap_fraction = GenomicRanges::width(ov) /
        GenomicRanges::width(genes_gr[sh]),
      bedpe_record_id = NA_character_,
      evidence_source = "linear",
      stringsAsFactors = FALSE))
  }
  # fully_contained rows
  if (length(contained) > 0) {
    qh <- S4Vectors::queryHits(contained); sh <- S4Vectors::subjectHits(contained)
    links <- rbind(links, data.frame(
      domain_id = domain_ids[sh],
      gene_id = mcols(genes_gr)$gene_id[qh],
      gene_symbol = NA_character_,
      relation_type = "fully_contained",
      distance_to_tss_bp = NA_real_,
      overlap_bp = GenomicRanges::width(genes_gr[qh]),
      domain_overlap_fraction = GenomicRanges::width(genes_gr[qh]) /
        GenomicRanges::width(domains[sh]),
      feature_overlap_fraction = 1,
      bedpe_record_id = NA_character_,
      evidence_source = "linear",
      stringsAsFactors = FALSE))
  }

  # ---- BEDPE 3D evidence (design §39-46) ------------------------------------
  bedpe_prov <- NULL
  if (!is.null(bedpe)) {
    bp_res <- .load_bedpe(bedpe)
    bedpe_prov <- bp_res$provenance
    # P1-5: seqlevel compatibility must actually be enforced for BEDPE. A
    # chr1-vs-1 style mismatch would silently yield zero contacts.
    dom_sl <- unique(as.character(GenomicRanges::seqnames(domains)))
    bp_sl <- unique(c(as.character(GenomicRanges::seqnames(bp_res$gr)),
                      as.character(GenomicRanges::seqnames(bp_res$anchor2))))
    shared <- intersect(dom_sl, bp_sl)
    if (length(shared) == 0) {
      stop("No shared seqlevels between the domains and the BEDPE anchors. ",
           "Check that both use the same chromosome naming convention ",
           "(e.g. chr1 vs 1).")
    }
    if (length(shared) < length(dom_sl)) {
      warning(sprintf(
        "BEDPE seqlevels match %d/%d domain seqlevels. Contacts on ",
        length(shared), length(dom_sl),
        "unmatched contigs will be absent."), call. = FALSE)
    }
    bedpe_links <- .bedpe_promoter_contacts(bp_res$gr, domains, domain_ids,
                                            promoters_gr)
    if (nrow(bedpe_links) > 0) links <- rbind(links, bedpe_links)
  }

  # P1-1: resolve gene symbols ONCE over all linked gene ids (nearest_tss,
  # promoter, gene body, fully contained, BEDPE), then backfill into the link
  # table and the rowData nearest_tss_gene_symbol column.
  symbol_map <- NULL
  if (!is.null(gene_id_keytype) && !is.null(res$anno_db) &&
      requireNamespace(res$anno_db, quietly = TRUE)) {
    all_gids <- unique(c(nearest_gene_id[!is.na(nearest_gene_id)], links$gene_id))
    all_gids <- unique(all_gids[!is.na(all_gids)])
    if (length(all_gids) > 0) {
      # P1-8: a user-provided keytype that fails must NOT be silently ignored.
      res_map <- tryCatch({
        sym <- AnnotationDbi::mapIds(get(res$anno_db, asNamespace(res$anno_db)),
                                     keys = all_gids, column = "SYMBOL",
                                     keytype = gene_id_keytype, multiVals = "first")
        stats::setNames(as.character(sym), all_gids)
      }, error = function(e) e)
      if (inherits(res_map, "error")) {
        warning(sprintf(
          "gene symbol mapping failed with keytype '%s' on %s: %s. ",
          gene_id_keytype, res$anno_db, conditionMessage(res_map),
          "gene_symbol will be NA."), call. = FALSE)
        symbol_map <- NULL
      } else {
        symbol_map <- res_map
      }
    }
  }
  symbol_of <- function(gid) {
    # P1-7: gene_symbol must be a real symbol or NA, never a gene_id masquerading
    # as a symbol.
    if (is.null(symbol_map)) return(rep(NA_character_, length(gid)))
    out <- symbol_map[as.character(gid)]
    out[is.na(out)] <- NA_character_
    unname(out)
  }
  if (nrow(links) > 0) {
    links$gene_symbol <- symbol_of(links$gene_id)
  }
  rowData(se)$nearest_tss_gene_symbol <- symbol_of(nearest_gene_id)

  # ---- RNA-seq expression evidence (design §47-65) --------------------------
  expr_prov <- NULL
  expr_data <- NULL
  if (!is.null(expression)) {
    if (is.null(expression_type)) {
      stop("expression_type must be declared (e.g. 'TPM'). Never auto-guessed.")
    }
    expr <- .build_expression_matrix(expression, colData(se),
                                     aggregate_fun = aggregate_fun)
    # P1-10: explicit gene-ID match QC. TxDb gene_ids and RNA gene_ids are
    # often on different ID systems (Entrez vs Ensembl); report the match
    # fraction and fail on zero matches instead of silently returning NA.
    link_gids <- unique(links$gene_id[!is.na(links$gene_id)])
    expr_gids <- unique(expr$summary$gene_id[!is.na(expr$summary$gene_id)])
    n_matched <- length(intersect(link_gids, expr_gids))
    match_frac <- if (length(link_gids) > 0) n_matched / length(link_gids) else 0
    if (length(link_gids) > 0 && n_matched == 0) {
      stop(sprintf(
        "No RNA gene IDs match the annotation gene IDs (0/%d). Check that the ",
        "RNA gene_id system matches the TxDb gene_id system (e.g. Entrez vs ",
        "Ensembl), or map IDs before passing expression.",
        length(link_gids)))
    }
    if (length(link_gids) > 0 && match_frac < 0.5) {
      warning(sprintf(
        "Only %.0f%% of annotation gene IDs matched the RNA expression IDs. ",
        100 * match_frac,
        "Check that the gene_id systems are consistent."), call. = FALSE)
    }
    expr_prov <- list(expression_type = expression_type,
                      gene_id_type = expr$gene_id_type,
                      samples = expr$rna_samples,   # P1-9: RNA sample IDs
                      aggregation_method = aggregate_fun,
                      min_expression = NULL,
                      n_annotation_gene_ids = length(link_gids),
                      n_expression_gene_ids = length(expr_gids),
                      n_matched_gene_ids = n_matched,
                      gene_id_match_fraction = match_frac)
    expr_data <- expr$summary  # gene_id x Condition median matrix
  }

  S4Vectors::metadata(se)$domain_gene_links <- links
  S4Vectors::metadata(se)$annotation_provenance <- list(
    genome = res$genome_name,
    organism = res$organism,
    resource_mode = res$genome_class,
    txdb = res$resource_source,
    annotation_db = res$anno_db,
    gene_id_type = gm$gene_id_type,
    gene_id_keytype = gene_id_keytype,   # P1-8
    promoter_upstream_bp = promoter_upstream,
    promoter_downstream_bp = promoter_downstream,
    gene_model_source = "TxDb genes()",
    seqlevel_style = res$seqstyle,
    shared_seqlevels = sl_check$shared,
    domain_without_txdb = sl_check$domain_without_txdb,
    txdb_only_seqlevels = sl_check$txdb_only)
  if (!is.null(bedpe_prov)) S4Vectors::metadata(se)$bedpe_provenance <- bedpe_prov
  if (!is.null(expr_prov)) {
    S4Vectors::metadata(se)$expression_provenance <- expr_prov
    S4Vectors::metadata(se)$expression_summary <- expr_data
  }
  se
}


# Load a BEDPE file / data.frame into GRanges of the two anchors (design §40).
# Read-only: no filtering, merging, or calling. Original records are not
# modified or re-scored; bedpe_record_id links derived evidence back to the
# source file.
.load_bedpe <- function(bedpe) {
  if (is.character(bedpe)) {
    if (!file.exists(bedpe)) stop("BEDPE file not found: ", bedpe)
    tab <- utils::read.table(bedpe, header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE, fill = TRUE,
                             comment.char = "")
    src <- basename(bedpe)
  } else if (is.data.frame(bedpe)) {
    tab <- bedpe
    src <- "data.frame"
  } else {
    stop("bedpe must be a file path or data.frame.")
  }
  req <- seq_len(6)
  if (ncol(tab) < 6) stop("BEDPE requires at least 6 columns (chrom1..end2).")
  extra <- if (ncol(tab) > 6) tab[, 7:ncol(tab), drop = FALSE] else NULL
  n <- nrow(tab)
  # P1-6: strict input validation (chrom non-empty, start/end numeric, start>=0,
  # end>start). BEDPE is 0-based half-open.
  for (cc in c(1, 4)) {
    if (any(is.na(tab[[cc]]) | !nzchar(as.character(tab[[cc]])))) {
      stop("BEDPE chrom columns must be non-empty.")
    }
  }
  for (cc in c(2, 3, 5, 6)) {
    if (!is.numeric(tab[[cc]]) && !is.integer(tab[[cc]])) {
      stop(sprintf("BEDPE column %d (start/end) must be numeric.", cc))
    }
  }
  if (any(tab[[2]] < 0 | tab[[5]] < 0)) {
    stop("BEDPE start coordinates must be >= 0 (0-based).")
  }
  if (any(tab[[3]] <= tab[[2]] | tab[[6]] <= tab[[5]])) {
    stop("BEDPE end must be greater than start for both anchors.")
  }
  anchor1 <- GenomicRanges::GRanges(
    tab[, 1], IRanges::IRanges(tab[, 2] + 1, tab[, 3]))  # 0-based -> 1-based
  anchor2 <- GenomicRanges::GRanges(
    tab[, 4], IRanges::IRanges(tab[, 5] + 1, tab[, 6]))
  bedpe_id <- sprintf("BEDPE_%06d", seq_len(n))
  S4Vectors::mcols(anchor1)$bedpe_record_id <- bedpe_id
  S4Vectors::mcols(anchor2)$bedpe_record_id <- bedpe_id
  if (!is.null(extra)) {
    S4Vectors::mcols(anchor1) <- cbind(S4Vectors::mcols(anchor1), extra)
    S4Vectors::mcols(anchor2) <- cbind(S4Vectors::mcols(anchor2), extra)
  }
  # carry anchor2 into anchor1's metadata for symmetric contact linking
  S4Vectors::mcols(anchor1)$.anchor2 <- anchor2
  # P1-9: infer seqlevel style from the anchors instead of hard-coding "UCSC".
  sl1 <- as.character(GenomicRanges::seqnames(anchor1))
  seqstyle <- if (length(sl1) > 0 && all(grepl("^chr", sl1))) "UCSC" else "custom/unknown"
  list(gr = anchor1, anchor2 = anchor2, bedpe_id = bedpe_id,
       provenance = list(file_name = src, record_count = n,
                         column_names = colnames(tab),
                         coordinate_convention = "BEDPE 0-based half-open -> 1-based closed",
                         seqlevel_style = seqstyle,
                         note = paste(
                           "Original BEDPE records are not modified or re-scored;",
                           "bedpe_record_id links derived evidence back to the",
                           "source file.")))
}


# bedpe_promoter_contact evidence (design §41-45): symmetric / bidirectional.
# A domain-gene 3D contact is recorded when the domain overlaps one anchor AND
# the gene promoter overlaps the OTHER anchor (either orientation). Only
# promoter contacts are used as 3D evidence in v1.0.
.bedpe_promoter_contacts <- function(anchor1, domains, domain_ids, promoters_gr) {
  a1 <- anchor1
  a2 <- S4Vectors::mcols(a1)$.anchor2
  rec1 <- S4Vectors::mcols(a1)$bedpe_record_id
  rec2 <- S4Vectors::mcols(a2)$bedpe_record_id

  # P1-8: range-vectorized. For each direction, one findOverlaps for
  # domain<->anchor and one for promoter<->anchor, then join by
  # bedpe_record_id. No per-hit inner findOverlaps (scales to 10^4-10^6
  # interactions).
  .link_side <- function(d_hits, other_anchor, other_rec) {
    if (length(d_hits) == 0) return(data.frame())
    qh <- S4Vectors::queryHits(d_hits)
    sh <- S4Vectors::subjectHits(d_hits)
    # records whose OTHER anchor overlaps a promoter
    p_hits <- GenomicRanges::findOverlaps(promoters_gr, other_anchor)
    if (length(p_hits) == 0) return(data.frame())
    p_rec <- other_rec[S4Vectors::subjectHits(p_hits)]
    p_gene <- S4Vectors::mcols(promoters_gr)$gene_id[S4Vectors::queryHits(p_hits)]
    # join domain-side records to promoter-side records by bedpe_record_id
    d_rec <- S4Vectors::mcols(a1)$bedpe_record_id[sh]
    d_dom <- domain_ids[qh]
    d_df <- data.frame(domain_id = d_dom, bedpe_record_id = d_rec,
                       stringsAsFactors = FALSE)
    p_df <- data.frame(gene_id = p_gene, bedpe_record_id = p_rec,
                       stringsAsFactors = FALSE)
    merged <- merge(d_df, p_df, by = "bedpe_record_id",
                    all = FALSE, stringsAsFactors = FALSE)
    if (nrow(merged) == 0) return(data.frame())
    data.frame(
      domain_id = merged$domain_id,
      gene_id = merged$gene_id,
      gene_symbol = NA_character_,
      relation_type = "bedpe_promoter_contact",
      distance_to_tss_bp = NA_real_,
      overlap_bp = NA_real_,
      domain_overlap_fraction = NA_real_,
      feature_overlap_fraction = NA_real_,
      bedpe_record_id = merged$bedpe_record_id,
      evidence_source = "bedpe",
      stringsAsFactors = FALSE)
  }

  # direction 1: domain overlaps anchor1, promoter overlaps anchor2
  d1 <- GenomicRanges::findOverlaps(domains, a1)
  out <- .link_side(d1, a2, rec2)
  # direction 2: domain overlaps anchor2, promoter overlaps anchor1 (symmetric)
  d2 <- GenomicRanges::findOverlaps(domains, a2)
  out2 <- .link_side(d2, a1, rec1)
  rbind(out, out2)
}


# Build a per-condition per-gene expression summary from wide matrix or long
# table, aggregating RNA replicates by median (design §49-52).
.build_expression_matrix <- function(expression, col_data, aggregate_fun = "median") {
  agg <- match.fun(aggregate_fun)
  if (is.data.frame(expression)) {
    req <- c("gene_id", "SampleID", "Condition", "expression")
    if (!all(req %in% colnames(expression))) {
      stop("Long expression table requires columns: gene_id, SampleID, Condition, expression.")
    }
    # P0-5: aggregate with stats::aggregate so gene x replicate rows collapse
    # correctly to one row per gene per condition (tapply length mismatch fixed).
    summ <- stats::aggregate(expression ~ gene_id + Condition,
                             data = expression, FUN = agg)
    colnames(summ) <- c("gene_id", "Condition", "expression")
    rownames(summ) <- NULL
    return(list(summary = summ, gene_id_type = "user gene_id",
                rna_samples = unique(as.character(expression$SampleID))))
  }
  # wide matrix / data.frame: rows = genes, cols = RNA samples. Columns MUST be
  # matched by name against the epiPortrait sample IDs; positional matching is
  # unsafe for an independent RNA experiment (P0-6). If RNA sample IDs differ
  # from ChIP colnames, require the long-format input.
  m <- as.matrix(expression)
  # P1-3: validate gene rownames and numeric matrix BEFORE any aggregation.
  if (is.null(rownames(m)) || any(rownames(m) == "") || any(is.na(rownames(m)))) {
    stop("Wide expression matrix must have non-empty, unique gene rownames.")
  }
  if (anyDuplicated(rownames(m))) {
    stop("Wide expression matrix gene rownames must be unique.")
  }
  if (!is.numeric(m)) {
    stop("Wide expression matrix must be numeric.")
  }
  meta <- as.data.frame(col_data)
  if (is.null(colnames(m)) || any(is.na(colnames(m)))) {
    stop("Wide expression matrix must have column names matching rownames(colData(se)) (or use long format for independent RNA sample IDs).")
  }
  # P1-2: ANY unmatched RNA sample column is an error (never silently drop
  # columns). A subset of the ChIP samples is allowed, but unknown RNA samples
  # must be flagged so the user switches to long format for independent data.
  unmatched <- setdiff(colnames(m), rownames(meta))
  if (length(unmatched) > 0) {
    stop("Wide expression matrix contains RNA sample IDs not present in ",
         "colData(se): ", paste(head(unmatched, 5), collapse = ", "),
         ". Use the long-format input (with an explicit Condition column) ",
         "for independent RNA-seq samples.")
  }
  matched <- colnames(m)
  if (length(matched) == 0) {
    stop("For independent RNA-seq sample IDs, use the long-format input with an explicit Condition column (wide matrix columns must match rownames(colData(se))).")
  }
  cond_col <- intersect(c("Condition"), colnames(meta))
  conds <- if (length(cond_col) > 0) unique(meta[[cond_col[1]]]) else "ALL"
  summ <- do.call(rbind, lapply(conds, function(cn) {
    idx <- if (length(cond_col) > 0) {
      which(meta[[cond_col[1]]] == cn & rownames(meta) %in% matched)
    } else seq_len(ncol(m))
    idx <- idx[rownames(meta)[idx] %in% matched]
    if (length(idx) == 0) return(NULL)
    # re-map RNA columns by sample name
    rna_cols <- rownames(meta)[idx]
    rna_idx <- match(rna_cols, colnames(m))
    rna_idx <- rna_idx[!is.na(rna_idx)]
    data.frame(gene_id = rownames(m),
               Condition = cn,
               expression = apply(m[, rna_idx, drop = FALSE], 1, agg),
               stringsAsFactors = FALSE)
  }))
  rownames(summ) <- NULL
  list(summary = summ, gene_id_type = "user gene_id",
       rna_samples = matched)
}


#' Extract and Prioritize Candidate Genes for Domains
#'
#' @description Extracts candidate (domain-associated) genes from the stored
#'   domain-gene evidence (metadata(se)$domain_gene_links), optionally applying
#'   RNA-seq expression prioritization. Transparent ordinal prioritization:
#'   evidence relationships come first, expression is a secondary tie-breaker;
#'   no black-box composite score is computed.
#'
#' @param se A SummarizedExperiment after \code{annotate_epi_domains()}.
#' @param domains Logical vector (length nrow(se)) selecting which domains to
#'   consider, or NULL for all.
#' @param relations Character. Evidence relation types to include
#'   (default all: nearest_tss, promoter_overlap, gene_body_overlap,
#'   fully_contained, bedpe_promoter_contact).
#' @param group Character or NULL. Condition group whose expression is used for
#'   prioritization (NULL = no expression).
#' @param expression_priority Character. "none" (default), "expressed_first",
#'   or "high_expression_first".
#' @param min_expression Numeric or NULL. Descriptive threshold for
#'   \code{expression_status} = "Expressed". For TPM/CPM with
#'   \code{expression_priority = "expressed_first"}, a value like 1 is a
#'   convenience default. For VST / rlog / normalized_counts, expressed_first
#'   requires an explicit user threshold (no universal meaning of 1). For
#'   \code{"high_expression_first"} no threshold is needed.
#' @param unique_genes Logical. Collapse to unique genes (default TRUE),
#'   retaining the highest-priority domain association per gene.
#' @return A data.frame: domain_id, gene_id, gene_symbol, relation_types,
#'   bedpe_supported, bedpe_support_count, expression_value,
#'   expression_status, expression_rank, candidate_priority.
#' @examples
#' data(example_se)
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   se <- annotate_epi_domains(example_se, genome = "hg38")
#'   head(get_domain_genes(se))
#' }
#' @export
get_domain_genes <- function(se, domains = NULL,
                             relations = c("nearest_tss", "promoter_overlap",
                                           "gene_body_overlap", "fully_contained",
                                           "bedpe_promoter_contact"),
                             group = NULL,
                             expression_priority = c("none", "expressed_first",
                                                     "high_expression_first"),
                             min_expression = NULL, unique_genes = TRUE) {
  expression_priority <- match.arg(expression_priority)
  # ---- input validation (P1-11) --------------------------------------------
  if (!is.null(domains)) {
    if (length(domains) != nrow(se) || any(is.na(domains))) {
      stop("domains must be a logical vector of length nrow(se) with no NAs.")
    }
  }
  expr_prov <- S4Vectors::metadata(se)$expression_provenance
  has_expr <- !is.null(expr_prov)
  if (expression_priority != "none" && is.null(group)) {
    stop("expression_priority != 'none' requires a `group` (condition group).")
  }
  if (expression_priority != "none" && !has_expr) {
    stop("expression_priority != 'none' but the object has no expression data. ",
         "Run annotate_epi_domains() with `expression` first.")
  }
  if (!is.null(group) && has_expr) {
    es_all <- S4Vectors::metadata(se)$expression_summary
    if (!is.null(es_all) && !group %in% es_all$Condition) {
      stop(sprintf("group '%s' not found in the expression conditions (%s).",
                   group, paste(unique(es_all$Condition), collapse = ", ")))
    }
  }
  # P1-1: resolve min_expression semantics from the declared expression_type.
  etype <- if (has_expr) expr_prov$expression_type else NULL
  threshold <- min_expression
  if (expression_priority == "expressed_first") {
    if (is.null(threshold)) {
      if (!is.null(etype) && tolower(etype) %in% c("tpm", "cpm")) {
        threshold <- 1   # descriptive convenience default (documented)
      } else {
        stop("expression_priority = 'expressed_first' requires an explicit ",
             "min_expression threshold for expression_type '", etype, "' ",
             "(no universal meaning of 1 for VST/rlog/normalized_counts).")
      }
    }
  }
  links <- S4Vectors::metadata(se)$domain_gene_links
  if (is.null(links) || nrow(links) == 0) {
    stop("No domain-gene links found. Run annotate_epi_domains() first.")
  }
  if (!is.null(domains)) {
    keep_dom <- rownames(se)[as.logical(domains)]
    links <- links[links$domain_id %in% keep_dom, , drop = FALSE]
  }
  links <- links[links$relation_type %in% relations, , drop = FALSE]
  if (nrow(links) == 0) {
    return(data.frame())
  }

  # collapse per domain-gene
  key <- paste(links$domain_id, links$gene_id, sep = "|")
  dd <- unique(key)
  out <- do.call(rbind, lapply(dd, function(k) {
    sub <- links[key == k, ]
    bedpe_sub <- sub[sub$evidence_source == "bedpe", ]
    data.frame(
      domain_id = sub$domain_id[1],
      gene_id = sub$gene_id[1],
      gene_symbol = if (any(!is.na(sub$gene_symbol))) sub$gene_symbol[!is.na(sub$gene_symbol)][1] else NA_character_,
      relation_types = paste(sort(unique(sub$relation_type)), collapse = ";"),
      bedpe_supported = any(sub$evidence_source == "bedpe"),
      # P1-4: count UNIQUE BEDPE records (a single record may match via both
      # directions -> two rows; must not be double counted).
      bedpe_support_count = if (nrow(bedpe_sub) > 0)
        length(unique(stats::na.omit(bedpe_sub$bedpe_record_id))) else 0L,
      stringsAsFactors = FALSE)
  }))

  # expression integration
  out$expression_value <- NA_real_
  out$expression_status <- NA_character_
  out$expression_rank <- NA_integer_
  if (!is.null(group)) {
    expr_summ <- S4Vectors::metadata(se)$expression_summary
    if (!is.null(expr_summ)) {
      es <- expr_summ[expr_summ$Condition == group, , drop = FALSE]
      if (nrow(es) > 0) {
        m <- match(out$gene_id, es$gene_id)
        out$expression_value <- es$expression[m]
        # P1-1: expression_status only when a threshold is in effect (i.e. for
        # expressed_first); high_expression_first needs no threshold and leaves
        # status NA while still ranking by expression_value.
        if (!is.null(threshold)) {
          out$expression_status <- ifelse(is.na(out$expression_value), NA_character_,
                                          ifelse(out$expression_value >= threshold,
                                                 "Expressed", "Not_detected"))
        }
      }
    }
  }

  # ordinal prioritization (design §57-61): evidence hierarchy first, then
  # expression as secondary. Use pmax() so the strongest evidence is NEVER
  # overwritten by a weaker relation on the same domain-gene pair (P0-8).
  order_score <- rep(1L, nrow(out))
  order_score <- pmax(order_score,
                      ifelse(grepl("gene_body_overlap|fully_contained",
                                   out$relation_types), 1L, 1L))
  order_score <- pmax(order_score,
                      ifelse(grepl("nearest_tss", out$relation_types), 2L, 1L))
  order_score <- pmax(order_score,
                      ifelse(grepl("promoter_overlap", out$relation_types), 2L, 1L))
  order_score <- pmax(order_score,
                      ifelse(grepl("bedpe_promoter_contact", out$relation_types),
                             3L, 1L))
  # secondary: expression
  expr_tie <- rep(0, nrow(out))
  if (expression_priority == "expressed_first") {
    expr_tie[out$expression_status == "Expressed"] <- 1
  } else if (expression_priority == "high_expression_first") {
    expr_tie <- ifelse(is.na(out$expression_value), 0, out$expression_value)
  }
  # P0-7: order() returns the permutation; assign ranks correctly.
  ord <- order(order_score, expr_tie, decreasing = TRUE, na.last = TRUE)
  priority <- integer(nrow(out))
  priority[ord] <- seq_along(ord)
  out$candidate_priority <- priority
  # expression_rank: descending expression
  if (any(!is.na(out$expression_value))) {
    expr_ord <- order(-out$expression_value, na.last = TRUE)
    out$expression_rank <- match(seq_len(nrow(out)), expr_ord)
  }
  # P1-3: when collapsing to unique genes, keep the BEST-SUPPORTED domain-gene
  # association (highest priority), not the first row.
  if (unique_genes) {
    out <- out[order(out$candidate_priority, out$domain_id), ]
    out <- out[!duplicated(out$gene_id), , drop = FALSE]
  }
  # P1-4: candidate_priority is always a contiguous ordinal rank of the FINAL
  # output table (a gap can appear when unique_genes drops intermediate rows).
  out <- out[order(out$candidate_priority, out$domain_id), ]
  out$candidate_priority <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}
