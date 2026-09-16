# Domain-aware annotation and candidate-gene prioritization.
#
# epiPortrait annotates the shared final
# domain (H3K4me3 broad domains, H3K27ac stitched SE-like domains, H3K27me3 /
# H3K9me3 repressive domains), which can span tens to hundreds of kb. A single
# nearest-gene label does not capture all relationships for such intervals.
# This module provides genome-flexible, domain-aware annotation with optional
# BEDPE 3D-contact and RNA-seq expression evidence, plus explicit
# candidate-gene prioritization. It does NOT call peaks, loops, or fit
# differential expression models.

# Resolve genome resources (TxDb / OrgDb / blacklist) from a genome shortcut or
# user-supplied objects. annotate_epi_domains() and
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
    stop("genome must be a character shortcut ('hg38','hg19','mm10') or a TxDb.", call. = FALSE)
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
# silently converting 1<->chr1 or MT<->chrM. When enforce = TRUE
# (used by annotate_epi_domains), zero shared domain/TxDb seqlevels is a hard
# error and partial overlap is a warning, so a genome-build / naming mismatch
# is never silently reinterpreted as "intergenic domains".
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
# strand-aware TSS, gene body, and a per-gene promoter window.
# Transcripts are collapsed to gene level (multi-transcript -> gene).
.gene_model_from_txdb <- function(txdb, promoter_upstream = 3000,
                                  promoter_downstream = 3000) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE))
    stop("Please install 'GenomicFeatures'.")
  # one TxDb gene = one gene-level feature. NEVER reduce() across genes: two
  # overlapping genes must remain separate features.
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


# Typed empty raw-link table so an object annotated on domains with NO gene
# evidence still carries the documented schema instead of a 0x0 data.frame.
.empty_links_table <- function() {
  data.frame(domain_id = character(0), gene_id = character(0),
             gene_symbol = character(0), relation_type = character(0),
             distance_to_tss_bp = numeric(0), overlap_bp = numeric(0),
             domain_overlap_fraction = numeric(0),
             feature_overlap_fraction = numeric(0),
             bedpe_record_id = character(0), evidence_source = character(0),
             contact_score = numeric(0), annotation_source = character(0),
             stringsAsFactors = FALSE)
}


# Single implementation of the domain-gene link -> export-view semantics.
#
# Both the native annotator (annotate_epi_domains) and the external importer
# (import_domain_annotations) funnel their raw long-table `links` through this
# function. It is the ONE place that defines the evidence hierarchy, the
# per-pair dedup collapse and the per-domain master summary, so the two entry
# points cannot drift. `links` must use the documented link schema; the only
# requirement is that relation_type / evidence_source / domain_id / gene_id are
# present (plus the standard columns produced by .empty_links_table()).
.rebuild_annotation_views <- function(se, links, nearest_tss_cutoff_bp) {
  dom_ids <- rownames(se)
  n_dom <- nrow(se)
  split_domain <- split(seq_len(nrow(links)), links$domain_id)

  count_relation <- function(types) {
    out <- stats::setNames(integer(n_dom), dom_ids)
    keep <- links$relation_type %in% types
    if (any(keep)) {
      values <- vapply(split(links$gene_id[keep], links$domain_id[keep]),
                       function(z) length(unique(z)), integer(1))
      out[names(values)] <- values
    }
    unname(out)
  }
  promoter_count <- count_relation("promoter_overlap")
  body_count <- count_relation("gene_body_overlap")
  contained_count <- count_relation("fully_contained")

  nearest_id <- nearest_symbol <- rep(NA_character_, n_dom)
  nearest_distance <- rep(NA_real_, n_dom)
  names(nearest_id) <- names(nearest_symbol) <- names(nearest_distance) <- dom_ids
  nearest_rows <- which(links$relation_type == "nearest_tss")
  if (length(nearest_rows) > 0L) {
    nearest_split <- split(nearest_rows, links$domain_id[nearest_rows])
    for (domain in names(nearest_split)) {
      i <- nearest_split[[domain]]
      finite <- which(is.finite(links$distance_to_tss_bp[i]))
      chosen <- if (length(finite) > 0L) {
        i[finite[which.min(abs(links$distance_to_tss_bp[i][finite]))]]
      } else {
        i[1L]
      }
      nearest_id[domain] <- links$gene_id[chosen]
      nearest_symbol[domain] <- links$gene_symbol[chosen]
      nearest_distance[domain] <- links$distance_to_tss_bp[chosen]
    }
  }

  context <- rep("Intergenic", n_dom)
  context[body_count > 0L | contained_count > 0L] <- "Gene-body-associated"
  context[promoter_count > 0L] <- "Promoter-associated"
  rowData(se)$primary_genomic_context <- context
  rowData(se)$nearest_tss_gene_id <- unname(nearest_id[dom_ids])
  rowData(se)$nearest_tss_gene_symbol <- unname(nearest_symbol[dom_ids])
  rowData(se)$nearest_tss_distance_bp <- unname(nearest_distance[dom_ids])
  rowData(se)$promoter_overlap_gene_count <- promoter_count
  rowData(se)$gene_body_overlap_gene_count <- body_count
  rowData(se)$fully_contained_gene_count <- contained_count

  row_tier <- integer(nrow(links))
  row_tier[links$relation_type == "promoter_overlap"] <- 4L
  row_tier[links$relation_type == "bedpe_promoter_contact"] <- 3L
  proximal <- links$relation_type == "nearest_tss" &
    is.finite(links$distance_to_tss_bp) &
    abs(links$distance_to_tss_bp) <= nearest_tss_cutoff_bp
  row_tier[proximal] <- 2L
  body <- links$relation_type %in% c("gene_body_overlap", "fully_contained")
  row_tier[body] <- pmax(row_tier[body], 1L)

  key <- paste(links$domain_id, links$gene_id, sep = "|")
  pairs <- split(seq_len(nrow(links)), key)
  pair_value <- function(i, column) {
    z <- links[[column]][i]
    z <- z[!is.na(z) & nzchar(as.character(z))]
    if (length(z) == 0L) NA_character_ else as.character(z[1L])
  }
  pair_distance <- function(i) {
    z <- links$distance_to_tss_bp[i]
    z <- z[is.finite(z)]
    if (length(z) == 0L) NA_real_ else z[which.min(abs(z))]
  }
  bedpe_summary <- function(i) {
    b <- i[links$relation_type[i] == "bedpe_promoter_contact"]
    if (length(b) == 0L) return(c(0, 0, 0, 0, 0))
    keep <- !duplicated(links$bedpe_record_id[b])
    b <- b[keep]
    score <- links$contact_score[b]
    if (all(is.na(score))) {
      c(length(b), NA_real_, NA_real_, NA_real_, length(b))
    } else {
      c(length(b), sum(score, na.rm = TRUE), max(score, na.rm = TRUE),
        mean(score, na.rm = TRUE), length(b))
    }
  }
  bedpe_values <- t(vapply(pairs, bedpe_summary, numeric(5)))
  best_index <- vapply(pairs, function(i) i[which.max(row_tier[i])], integer(1))
  dedup <- data.frame(
    domain_id = vapply(pairs, function(i) links$domain_id[i[1L]], character(1)),
    gene_id = vapply(pairs, function(i) links$gene_id[i[1L]], character(1)),
    gene_symbol = vapply(pairs, pair_value, character(1), column = "gene_symbol"),
    relation_types = vapply(pairs, function(i) {
      paste(sort(unique(links$relation_type[i])), collapse = ";")
    }, character(1)),
    evidence_sources = vapply(pairs, function(i) {
      paste(sort(unique(links$evidence_source[i])), collapse = ";")
    }, character(1)),
    annotation_sources = vapply(pairs, function(i) {
      z <- unique(stats::na.omit(links$annotation_source[i]))
      paste(sort(z), collapse = ";")
    }, character(1)),
    best_relation = links$relation_type[best_index],
    best_evidence_source = links$evidence_source[best_index],
    best_tier = as.integer(vapply(pairs, function(i) max(row_tier[i]), integer(1))),
    nearest_tss_distance_bp = vapply(pairs, pair_distance, numeric(1)),
    bedpe_support_count = as.integer(bedpe_values[, 1L]),
    bedpe_contact_score = as.numeric(bedpe_values[, 2L]),
    bedpe_contact_score_max = as.numeric(bedpe_values[, 3L]),
    bedpe_contact_score_mean = as.numeric(bedpe_values[, 4L]),
    bedpe_contact_score_n = as.integer(bedpe_values[, 5L]),
    stringsAsFactors = FALSE)
  rownames(dedup) <- NULL

  n_linked <- stats::setNames(integer(n_dom), dom_ids)
  linked_values <- vapply(split(dedup$gene_id, dedup$domain_id),
                          function(z) length(unique(z)), integer(1))
  n_linked[names(linked_values)] <- linked_values
  n_bedpe_gene <- stats::setNames(integer(n_dom), dom_ids)
  bedpe_gene_values <- vapply(
    split(dedup$bedpe_support_count > 0L, dedup$domain_id), sum, integer(1))
  n_bedpe_gene[names(bedpe_gene_values)] <- bedpe_gene_values

  dom_bedpe <- lapply(split_domain, function(i) {
    b <- i[links$relation_type[i] == "bedpe_promoter_contact"]
    if (length(b) == 0L) return(c(0, 0, 0, 0))
    b <- b[!duplicated(links$bedpe_record_id[b])]
    score <- links$contact_score[b]
    if (all(is.na(score))) return(c(NA_real_, NA_real_, NA_real_, length(b)))
    c(sum(score, na.rm = TRUE), max(score, na.rm = TRUE),
      mean(score, na.rm = TRUE), length(b))
  })
  dom_score <- dom_max <- dom_mean <- stats::setNames(numeric(n_dom), dom_ids)
  dom_n <- stats::setNames(integer(n_dom), dom_ids)
  for (domain in names(dom_bedpe)) {
    z <- dom_bedpe[[domain]]
    dom_score[domain] <- z[1L]
    dom_max[domain] <- z[2L]
    dom_mean[domain] <- z[3L]
    dom_n[domain] <- as.integer(z[4L])
  }

  ordered <- dedup[order(match(dedup$domain_id, dom_ids), -dedup$best_tier,
                         is.na(dedup$gene_symbol), dedup$gene_symbol,
                         dedup$gene_id), , drop = FALSE]
  top <- ordered[!duplicated(ordered$domain_id), , drop = FALSE]
  top_symbol <- stats::setNames(top$gene_symbol, top$domain_id)
  top_id <- stats::setNames(top$gene_id, top$domain_id)

  summary <- data.frame(
    Domain_ID = dom_ids,
    primary_genomic_context = context,
    nearest_tss_gene_symbol = unname(nearest_symbol[dom_ids]),
    nearest_tss_distance_bp = unname(nearest_distance[dom_ids]),
    n_promoter_overlap_gene = promoter_count,
    n_gene_body_overlap_gene = body_count,
    n_fully_contained_gene = contained_count,
    n_linked_gene = unname(n_linked[dom_ids]),
    n_bedpe_contact_gene = unname(n_bedpe_gene[dom_ids]),
    bedpe_contact_score = unname(dom_score[dom_ids]),
    bedpe_contact_score_max = unname(dom_max[dom_ids]),
    bedpe_contact_score_mean = unname(dom_mean[dom_ids]),
    bedpe_contact_score_n = unname(dom_n[dom_ids]),
    top_candidate_gene_symbol = unname(top_symbol[dom_ids]),
    top_candidate_gene_id = unname(top_id[dom_ids]),
    stringsAsFactors = FALSE)

  rowData(se)$n_bedpe_contact_gene <- summary$n_bedpe_contact_gene
  rowData(se)$bedpe_contact_score <- summary$bedpe_contact_score
  rowData(se)$bedpe_contact_score_max <- summary$bedpe_contact_score_max
  rowData(se)$bedpe_contact_score_mean <- summary$bedpe_contact_score_mean
  rowData(se)$bedpe_contact_score_n <- summary$bedpe_contact_score_n
  rowData(se)$top_candidate_gene_symbol <- summary$top_candidate_gene_symbol
  S4Vectors::metadata(se)$domain_gene_links <- links
  S4Vectors::metadata(se)$domain_gene_links_dedup <- dedup
  S4Vectors::metadata(se)$annotation_summary <- summary
  se
}


# Assemble the linear (TxDb) domain-gene evidence rows: nearest-TSS,
# promoter-overlap, gene-body-overlap and fully-contained relationships.
# Returns the raw long-table links (through `evidence_source`); gene_symbol is
# filled later by .resolve_link_gene_symbols() and the per-domain context /
# counts / summary are derived by .rebuild_annotation_views().
.linear_domain_gene_links <- function(domains, domain_ids, genes_gr,
                                      promoters_gr, nearest_gene_id,
                                      signed_dist, n_dom) {
  prom_ov <- GenomicRanges::findOverlaps(domains, promoters_gr)
  body_ov <- GenomicRanges::findOverlaps(domains, genes_gr, ignore.strand = TRUE)
  contained <- GenomicRanges::findOverlaps(genes_gr, domains,
                                           ignore.strand = TRUE, type = "within")

  links <- data.frame()
  # nearest_tss rows
  if (n_dom > 0 && any(!is.na(nearest_gene_id))) {
    links <- rbind(links, data.frame(
      domain_id = domain_ids[!is.na(nearest_gene_id)],
      gene_id = nearest_gene_id[!is.na(nearest_gene_id)],
      gene_symbol = NA_character_,
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
  links
}


# Load optional BEDPE 3D evidence, enforce seqlevel compatibility, and derive
# the promoter-contact link rows. Returns NULL components when no BEDPE is
# supplied. A chr1-vs-1 style seqlevel mismatch is a hard error because it would
# otherwise silently yield zero contacts.
.bedpe_domain_gene_links <- function(domains, domain_ids, promoters_gr, bedpe,
                                     bedpe_score_col, min_anchor_overlap_bp) {
  if (is.null(bedpe)) {
    return(list(bp_res = NULL, bedpe_prov = NULL, bedpe_links = NULL))
  }
  bp_res <- .load_bedpe(bedpe, score_col = bedpe_score_col)
  bedpe_prov <- bp_res$provenance
  bedpe_prov$min_anchor_overlap_bp <-
    if (is.null(min_anchor_overlap_bp)) 1 else min_anchor_overlap_bp
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
      "BEDPE seqlevels match %d/%d domain seqlevels. Contacts on unmatched contigs will be absent.",
      length(shared), length(dom_sl)), call. = FALSE)
  }
  bedpe_links <- .bedpe_promoter_contacts(bp_res$gr, domains, domain_ids,
                                          promoters_gr,
                                          min_overlap_bp =
                                            if (is.null(min_anchor_overlap_bp)) 1 else
                                              min_anchor_overlap_bp)
  list(bp_res = bp_res, bedpe_prov = bedpe_prov, bedpe_links = bedpe_links)
}


# Resolve gene symbols once over all linked gene IDs using the TxDb's associated
# OrgDb. gene_symbol is a real symbol or NA, never a gene_id masquerading as a
# symbol. A user-supplied keytype that fails mapping is a warning, not a silent
# NA.
.resolve_link_gene_symbols <- function(links, nearest_gene_id, gene_id_keytype,
                                       anno_db) {
  symbol_map <- NULL
  if (!is.null(gene_id_keytype) && !is.null(anno_db) &&
      requireNamespace(anno_db, quietly = TRUE)) {
    all_gids <- unique(c(nearest_gene_id[!is.na(nearest_gene_id)], links$gene_id))
    all_gids <- unique(all_gids[!is.na(all_gids)])
    if (length(all_gids) > 0) {
      res_map <- tryCatch({
        sym <- AnnotationDbi::mapIds(get(anno_db, asNamespace(anno_db)),
                                     keys = all_gids, column = "SYMBOL",
                                     keytype = gene_id_keytype, multiVals = "first")
        stats::setNames(as.character(sym), all_gids)
      }, error = function(e) e)
      if (inherits(res_map, "error")) {
        warning(sprintf(
          "gene symbol mapping failed with keytype '%s' on %s: %s. gene_symbol will be NA.",
          gene_id_keytype, anno_db, conditionMessage(res_map)), call. = FALSE)
        symbol_map <- NULL
      } else {
        symbol_map <- res_map
      }
    }
  }
  if (nrow(links) > 0) {
    links$gene_symbol <- if (is.null(symbol_map)) {
      rep(NA_character_, nrow(links))
    } else {
      out <- symbol_map[as.character(links$gene_id)]
      out[is.na(out)] <- NA_character_
      unname(out)
    }
  }
  links
}


# Join optional per-record contact scores onto the raw links by
# bedpe_record_id. Linear rows keep contact_score = NA. A zero-row links table
# is returned untouched so the caller can substitute the typed empty schema.
.attach_bedpe_contact_scores <- function(links, bp_res) {
  if (nrow(links) == 0) return(links)
  links$contact_score <- NA_real_
  if (!is.null(bp_res) && !is.null(bp_res$score)) {
    is_bp_row <- !is.na(links$bedpe_record_id) & links$evidence_source == "bedpe"
    score_lookup <- stats::setNames(bp_res$score, bp_res$bedpe_id)
    hit <- is_bp_row & links$bedpe_record_id %in% names(score_lookup)
    links$contact_score[hit] <-
      unname(score_lookup[links$bedpe_record_id[hit]])
  }
  links
}


# Build the per-condition expression summary and its provenance, with explicit
# gene-ID matching diagnostics against the annotation gene IDs (TxDb and RNA
# IDs are often on different ID systems).
.build_expression_evidence <- function(links, expression, expression_type,
                                       aggregate_fun, col_data) {
  if (is.null(expression)) {
    return(list(expr_prov = NULL, expr_data = NULL))
  }
  if (is.null(expression_type)) {
    stop("expression_type must be declared (e.g. 'TPM'). Never auto-guessed.")
  }
  expr <- .build_expression_matrix(expression, col_data,
                                   aggregate_fun = aggregate_fun)
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
      "Only %.0f%% of annotation gene IDs matched the RNA expression IDs. Check that the gene_id systems are consistent.",
      100 * match_frac), call. = FALSE)
  }
  expr_prov <- list(expression_type = expression_type,
                    gene_id_type = expr$gene_id_type,
                    samples = expr$rna_samples,
                    aggregation_method = aggregate_fun,
                    min_expression = NULL,
                    n_annotation_gene_ids = length(link_gids),
                    n_expression_gene_ids = length(expr_gids),
                    n_matched_gene_ids = n_matched,
                    gene_id_match_fraction = match_frac)
  list(expr_prov = expr_prov, expr_data = expr$summary)
}


#' Domain-Aware Annotation of epiPortrait Domains
#'
#' @description Annotates the final shared domains (not individual peaks) with
#' genome-aware linear evidence plus optional BEDPE 3D-contact and RNA-seq
#' expression evidence. Adds 13 compact rowData columns (7 linear/symbol
#' columns plus 6 BEDPE-contact and candidate columns) and stores a long-format
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
#' @param genome Character. Built-in shortcut ("hg38","hg19","mm10") or a
#'   TxDb object. For other genomes (e.g. mm9, rn7) supply \code{txdb}
#'   explicitly; core quantitative analysis is genome-agnostic.
#' @param txdb A TxDb object or NULL (resolved from genome shortcut).
#' @param anno_db Character. OrgDb package name, or NULL (resolved).
#' @param promoter_upstream Numeric. Promoter window upstream of the TSS
#'   (default 3000 bp).
#' @param promoter_downstream Numeric. Promoter window downstream of the TSS
#'   (default 3000 bp).
#' @param nearest_tss_cutoff_bp Numeric. Distance (bp) below which a
#'   \code{nearest_tss} relation counts as strong proximal evidence (tier 2)
#'   when \code{best_relation} / \code{best_tier} are computed for the dedup
#'   and summary tables. Must match the value used by
#'   \code{get_domain_genes()} (which reads the same
#'   \code{epiPortrait.nearest_tss_cutoff_bp} option by default) so the two
#'   layers cannot drift apart. Default 10000.
#' @param bedpe Character path to a BEDPE file, a data.frame, or NULL. Only
#'   \code{bedpe_promoter_contact} evidence is added (3D promoter-contact);
#'   records are read-only (no calling / filtering / merging).
#' @param bedpe_score_col Numeric or character or NULL. Optional column holding
#'   a per-record contact STRENGTH (e.g. the FitHiChIP TMM-normalized contact
#'   frequency in column 8). Integer column index (>= 7; columns 1-6 are the
#'   coordinates) or a column name (file inputs expose V1..Vn names; a
#'   data.frame input can carry arbitrary names). When given, quantitative
#'   outputs are produced IN ADDITION to the count-based evidence:
#'   \code{metadata(se)$domain_gene_links$contact_score} (per evidence row), and
#'   four per-record aggregates (\code{sum} / \code{max} / \code{mean} /
#'   \code{n}) over the UNIQUE supporting records, exposed as
#'   \code{bedpe_contact_score(_max/_mean/_n)} in both the dedup pair table and
#'   the annotation summary (and mirrored on rowData). Rationale: FitHiChIP-style
#'   TMM BEDPEs often share an identical loop SET across conditions and vary only
#'   in score, so count-based integration alone cannot capture contact dynamics;
#'   \code{sum} reflects total contact burden, \code{max} a single dominant
#'   interaction, \code{mean} per-record intensity, and \code{n} record count.
#'   NULL (default) keeps the previous count-only behaviour unchanged.
#' @param min_anchor_overlap_bp Numeric or NULL. Minimum base pairs of overlap
#'   required between a domain anchor and a BEDPE anchor for the 3D contact to
#'   be counted (default NULL = any overlap, i.e. >= 1 bp). A value of
#'   \code{1} requires at least 1 bp, \code{2} at least 2 bp, \code{100} at
#'   least 100 bp of anchor overlap; typical meaningful thresholds are tens to
#'   hundreds of base pairs depending on anchor resolution. Does NOT reject
#'   1-bp overlaps (use >= 2 for that). Applies to BOTH anchor sides
#'   (domain-anchor and promoter-anchor) of the domain-gene contact.
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
#'   nearest_tss_gene_id, nearest_tss_gene_symbol,
#'   nearest_tss_distance_bp (strand-aware SIGNED distance from the domain to
#'   the representative TSS; negative means the domain is upstream of a + gene
#'   or downstream of a - gene; 0 when the TSS is inside the domain),
#'   promoter_overlap_gene_count, gene_body_overlap_gene_count,
#'   fully_contained_gene_count, n_bedpe_contact_gene, bedpe_contact_score,
#'   bedpe_contact_score_max, bedpe_contact_score_mean, bedpe_contact_score_n.
#'   BEDPE score-summary columns are always present:
#'   \code{bedpe_contact_score}/\code{_max}/\code{_mean} are NA for supported
#'   contacts without a usable score (including runs without
#'   \code{bedpe_score_col}) and 0 for domains with no BEDPE evidence;
#'   \code{bedpe_contact_score_n} is the number of UNIQUE supporting BEDPE
#'   records (0 when there is no BEDPE evidence), so it mirrors
#'   \code{bedpe_support_count} and stays interpretable without a score column.
#'   Also adds top_candidate_gene_symbol; the annotation result is exposed as a
#'   three-level structure in metadata:
#'   \itemize{
#'     \item \code{metadata(se)$annotation_summary} — ONE ROW PER DOMAIN master
#'           table (nearest gene, per-evidence overlap/gene counts, number of
#'           BEDPE-contacted genes, and a representative top candidate gene).
#'     \item \code{metadata(se)$domain_gene_links_dedup} — ONE ROW PER
#'           domain-gene pair, with summarised \code{relation_types},
#'           \code{evidence_sources} and a \code{best_relation} label.
#'     \item \code{metadata(se)$domain_gene_links} — raw per-relationship
#'           detail (for auditing).
#'   }
#'   \code{metadata(se)$annotation_provenance}, and (if provided)
#'   \code{metadata(se)$bedpe_provenance} / \code{metadata(se)$expression_provenance}.
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
                                  nearest_tss_cutoff_bp =
                                    getOption("epiPortrait.nearest_tss_cutoff_bp",
                                              10000),
                                  bedpe = NULL, bedpe_score_col = NULL,
                                  min_anchor_overlap_bp = NULL,
                                  expression = NULL,
                                  expression_type = NULL,
                                  aggregate_fun = "median",
                                  gene_id_keytype = NULL) {
  # Validate promoter window parameters; negative values would be
  # "corrected" by pmin/pmax but are biologically meaningless).
  for (nm in c("promoter_upstream", "promoter_downstream")) {
    v <- get(nm)
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v < 0) {
      stop(sprintf("%s must be a finite non-negative number.", nm))
    }
  }
  # The proximal-TSS cutoff feeds best_tier / best_relation in the stored
  # tables; get_domain_genes() must use the same value or the two layers drift.
  if (length(nearest_tss_cutoff_bp) != 1L ||
      !is.numeric(nearest_tss_cutoff_bp) ||
      !is.finite(nearest_tss_cutoff_bp) || nearest_tss_cutoff_bp < 0) {
    stop("nearest_tss_cutoff_bp must be a finite non-negative number.")
  }
  # BEDPE contact-score column: only meaningful together with a BEDPE input.
  # Type/index validation happens inside .load_bedpe() where the table shape
  # is known; here we only catch the obvious misuse early.
  if (!is.null(bedpe_score_col) && is.null(bedpe)) {
    stop("bedpe_score_col requires a BEDPE input (`bedpe`).")
  }
  if (!is.null(bedpe_score_col) &&
      !is.numeric(bedpe_score_col) && !is.character(bedpe_score_col)) {
    stop("bedpe_score_col must be an integer column index or a column name.")
  }
  if (!is.null(min_anchor_overlap_bp)) {
    if (length(min_anchor_overlap_bp) != 1L || !is.numeric(min_anchor_overlap_bp) ||
        !is.finite(min_anchor_overlap_bp) || min_anchor_overlap_bp < 1) {
      stop("min_anchor_overlap_bp must be a finite number >= 1 (bp of anchor ",
           "overlap), or NULL for the default (any overlap).", call. = FALSE)
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

  domains <- rowRanges(se)
  domain_ids <- rownames(se)
  n_dom <- length(domains)

  # ---- seqlevel compatibility -----------------------------------------------
  sl_check <- .check_seqlevel_compatibility(domains, txdb = res$txdb,
                                            enforce = TRUE)

  # ---- linear evidence ------------------------------------------------------
  # nearest TSS: absolute shortest interval-to-point distance, backfilled by
  # queryHits() so domains without a hit (e.g. on contigs absent from the TxDb)
  # remain NA instead of misaligning.
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
  # TSS lies inside the domain, distance is 0, not the midpoint offset.
  # Sign: + strand -> positive downstream of TSS, negative upstream; - strand
  # is mirrored.
  # Vectorized to avoid an O(n_domains x n_genes) per-hit subset operation.
  # match() maps each hit to its gene's TSS and strand in one pass.
  signed_dist <- rep(NA_real_, n_dom)
  if (length(qh_t) > 0) {
    g_idx <- match(nearest_gene_id[qh_t], mcols(genes_gr)$gene_id)
    hit_ok <- !is.na(g_idx)
    if (any(hit_ok)) {
      d_hit <- nearest_dist[qh_t[hit_ok]]
      tss_hit <- mcols(genes_gr)$tss[g_idx[hit_ok]]
      s_hit <- mcols(genes_gr)$gene_strand[g_idx[hit_ok]]
      sd <- d_hit  # default: unstranded (or TSS inside domain -> distance 0)
      plus <- s_hit == "+"
      sd[plus] <- ifelse(GenomicRanges::end(domains)[qh_t[hit_ok][plus]] < tss_hit[plus],
                         -d_hit[plus], d_hit[plus])
      minus <- s_hit == "-"
      sd[minus] <- ifelse(GenomicRanges::start(domains)[qh_t[hit_ok][minus]] > tss_hit[minus],
                          -d_hit[minus], d_hit[minus])
      signed_dist[qh_t[hit_ok]] <- sd
    }
  }

  # Built-in UCSC knownGene uses ENTREZID; custom TxDbs may use other ID
  # systems, so require an explicit gene_id_keytype for custom genomes.
  if (is.null(gene_id_keytype) && res$genome_class == "builtin") {
    gene_id_keytype <- "ENTREZID"
  }

  # ---- long-format domain-gene links ----------------------------------------
  # Linear evidence rows (nearest-TSS, promoter overlap, gene body, fully
  # contained) are assembled by an internal helper. Primary genomic context,
  # per-domain counts and nearest-TSS columns are NOT computed here: they are
  # derived from the links by .rebuild_annotation_views() (the single shared
  # builder), so the native and external-import paths cannot drift.
  links <- .linear_domain_gene_links(domains, domain_ids, genes_gr, promoters_gr,
                                     nearest_gene_id, signed_dist, n_dom)

  # ---- BEDPE 3D evidence -----------------------------------------------------
  bedpe_res <- .bedpe_domain_gene_links(domains, domain_ids, promoters_gr, bedpe,
                                        bedpe_score_col, min_anchor_overlap_bp)
  bp_res <- bedpe_res$bp_res
  bedpe_prov <- bedpe_res$bedpe_prov
  if (!is.null(bedpe_res$bedpe_links) && nrow(bedpe_res$bedpe_links) > 0) {
    links <- rbind(links, bedpe_res$bedpe_links)
  }

  # Resolve gene symbols once over all linked gene IDs (nearest_tss, promoter,
  # gene body, fully contained, BEDPE).
  links <- .resolve_link_gene_symbols(links, nearest_gene_id, gene_id_keytype,
                                      res$anno_db)

  # ---- BEDPE record-level contact scores (optional, bedpe_score_col) --------
  # Scores live on the raw evidence rows and are joined by bedpe_record_id
  # (identical for both anchors of a record). Linear rows keep contact_score NA.
  links <- .attach_bedpe_contact_scores(links, bp_res)

  # ---- RNA-seq expression evidence ------------------------------------------
  expr_res <- .build_expression_evidence(links, expression, expression_type,
                                         aggregate_fun,
                                         SummarizedExperiment::colData(se))
  expr_prov <- expr_res$expr_prov
  expr_data <- expr_res$expr_data

  # A zero-link annotation (e.g. every domain on a contig that carries no
  # genes and no BEDPE contacts) is legitimate boundary data: store the
  # documented typed schema instead of a 0x0 data.frame so downstream column
  # access and exports stay type-stable.
  if (nrow(links) == 0) {
    links <- .empty_links_table()
  } else {
    # Provenance label for every raw evidence row. Linear TxDb evidence and
    # BEDPE contacts are both produced natively; import_domain_annotations()
    # records its external source in the SAME column.
    links$annotation_source <- ifelse(links$evidence_source == "bedpe",
                                      "epiPortrait:BEDPE", "epiPortrait:native")
  }

  # ---- user-facing export tables (single shared builder) --------------------
  # The annotation result is exposed at three levels:
  #   * metadata(se)$annotation_summary      : ONE ROW PER DOMAIN (master table)
  #   * metadata(se)$domain_gene_links_dedup : ONE ROW PER domain-gene PAIR
  #   * metadata(se)$domain_gene_links       : raw per-relationship detail (audit)
  # .rebuild_annotation_views() is the SINGLE implementation of the evidence
  # tier, dedup collapse and master-summary semantics. The external importer
  # calls the exact same function, so the two entry points cannot drift.
  se <- .rebuild_annotation_views(
    se, links, nearest_tss_cutoff_bp = nearest_tss_cutoff_bp)

  S4Vectors::metadata(se)$annotation_provenance <- list(
    genome = res$genome_name,
    organism = res$organism,
    resource_mode = res$genome_class,
    txdb = res$resource_source,
    annotation_db = res$anno_db,
    gene_id_type = gm$gene_id_type,
    gene_id_keytype = gene_id_keytype,
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


# Load a BEDPE file or data.frame into GRanges for the two anchors.
# Read-only: no filtering, merging, or calling. Original records are not
# modified or re-scored; bedpe_record_id links derived evidence back to the
# source file.
.load_bedpe <- function(bedpe, score_col = NULL) {
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
  # ---- optional contact-score column (bedpe_score_col) ----------------------
  # Resolved against the RAW table so both numeric indices and names work
  # (file inputs expose V1..Vn names; data.frame inputs keep their own).
  score_vals <- NULL
  score_label <- NULL
  if (!is.null(score_col)) {
    if (is.numeric(score_col)) {
      if (length(score_col) != 1L || is.na(score_col) ||
          floor(score_col) != score_col || score_col < 7 || score_col > ncol(tab)) {
        stop("bedpe_score_col index must be a single integer in [7, ",
             ncol(tab), "] (columns 1-6 are the coordinates).")
      }
      sidx <- as.integer(score_col)
    } else {
      cn <- colnames(tab)
      if (is.null(cn) || length(score_col) != 1L || !score_col %in% cn) {
        stop(sprintf("bedpe_score_col '%s' not found among BEDPE column names.",
                     paste(score_col, collapse = ",")), call. = FALSE)
      }
      sidx <- match(score_col, cn)
      if (sidx < 7) {
        stop("bedpe_score_col must point to an extra column (>= 7); ",
             "columns 1-6 are the coordinates.")
      }
    }
    score_label <- sprintf("column %d%s", sidx,
                           if (!is.null(colnames(tab)))
                             paste0(" (", colnames(tab)[sidx], ")") else "")
    score_vals <- as.numeric(tab[[sidx]])
    # as.numeric("Inf") gives Inf rather than NA; convert
    # ANY non-finite value (NA, NaN, Inf, -Inf) to NA so downstream sum/max/
    # mean/n never propagate Inf.
    bad <- !is.finite(score_vals)
    if (any(bad)) {
      score_vals[bad] <- NA_real_
      warning(sprintf(
        "%d/%d non-finite values in the BEDPE contact-score column were converted to NA.",
        sum(bad), length(score_vals)), call. = FALSE)
    }
  }
  # Strict input validation: chromosome names are non-empty, coordinates are
  # numeric and non-negative, and each end is greater than its start. BEDPE is
  # interpreted as 0-based half-open.
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
  # Infer seqlevel style from the anchors instead of hard-coding "UCSC".
  sl1 <- as.character(GenomicRanges::seqnames(anchor1))
  seqstyle <- if (length(sl1) > 0 && all(grepl("^chr", sl1))) "UCSC" else "custom/unknown"
  list(gr = anchor1, anchor2 = anchor2, bedpe_id = bedpe_id,
       score = score_vals,
       provenance = list(file_name = src, record_count = n,
                         column_names = colnames(tab),
                         coordinate_convention = "BEDPE 0-based half-open -> 1-based closed",
                         seqlevel_style = seqstyle,
                         score_column = score_label,
                         note = paste(
                           "Original BEDPE records are not modified or re-scored;",
                           "bedpe_record_id links derived evidence back to the",
                           "source file.")))
}


# bedpe_promoter_contact evidence is symmetric and bidirectional.
# A domain-gene 3D contact is recorded when the domain overlaps one anchor AND
# the gene promoter overlaps the OTHER anchor (either orientation). Only
# promoter contacts are used as 3D evidence in v1.0.
.bedpe_promoter_contacts <- function(anchor1, domains, domain_ids, promoters_gr,
                                     min_overlap_bp = 1) {
  a1 <- anchor1
  a2 <- S4Vectors::mcols(a1)$.anchor2
  rec1 <- S4Vectors::mcols(a1)$bedpe_record_id
  rec2 <- S4Vectors::mcols(a2)$bedpe_record_id
  if (length(min_overlap_bp) != 1L || !is.numeric(min_overlap_bp) ||
      !is.finite(min_overlap_bp) || min_overlap_bp < 1) {
    stop("min_overlap_bp must be a finite number >= 1 (bp of anchor overlap).")
  }

  # Range-vectorized: for each direction, one findOverlaps for
  # domain<->anchor and one for promoter<->anchor, then join by
  # bedpe_record_id. No per-hit inner findOverlaps (scales to 10^4-10^6
  # interactions).
  .link_side <- function(d_hits, dom_anchor, other_anchor, other_rec) {
    if (length(d_hits) == 0) return(data.frame())
    qh <- S4Vectors::queryHits(d_hits)
    sh <- S4Vectors::subjectHits(d_hits)
    # enforce a MINIMUM overlap width on the DOMAIN-anchor side (dom_anchor,
    # which is the anchor this direction overlaps: a1 for dir1, a2 for dir2)
    # and on the PROMOTER-anchor side, so edge-touching / 1-bp anchors do not
    # create spurious contacts.
    ov_w <- GenomicRanges::width(GenomicRanges::pintersect(
      domains[qh], dom_anchor[sh]))
    keep_d <- ov_w >= min_overlap_bp
    qh <- qh[keep_d]; sh <- sh[keep_d]
    if (length(qh) == 0) return(data.frame())
    # records whose OTHER anchor overlaps a promoter (with min overlap)
    p_hits <- GenomicRanges::findOverlaps(promoters_gr, other_anchor)
    if (length(p_hits) == 0) return(data.frame())
    p_qh <- S4Vectors::queryHits(p_hits); p_sh <- S4Vectors::subjectHits(p_hits)
    p_ov <- GenomicRanges::width(GenomicRanges::pintersect(
      promoters_gr[p_qh], other_anchor[p_sh]))
    p_keep <- p_ov >= min_overlap_bp
    p_hits <- p_hits[p_keep]
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
  out <- .link_side(d1, a1, a2, rec2)
  # direction 2: domain overlaps anchor2, promoter overlaps anchor1 (symmetric)
  d2 <- GenomicRanges::findOverlaps(domains, a2)
  out2 <- .link_side(d2, a2, a1, rec1)
  rbind(out, out2)
}


# Build a per-condition per-gene expression summary from wide matrix or long
# table, aggregating RNA replicates by median.
.build_expression_matrix <- function(expression, col_data, aggregate_fun = "median") {
  agg <- match.fun(aggregate_fun)
  if (is.data.frame(expression)) {
    req <- c("gene_id", "SampleID", "Condition", "expression")
    if (!all(req %in% colnames(expression))) {
      stop("Long expression table requires columns: gene_id, SampleID, Condition, expression.")
    }
    # Collapse gene-by-replicate rows to one row per gene and condition.
    summ <- stats::aggregate(expression ~ gene_id + Condition,
                             data = expression, FUN = agg)
    colnames(summ) <- c("gene_id", "Condition", "expression")
    rownames(summ) <- NULL
    return(list(summary = summ, gene_id_type = "user gene_id",
                rna_samples = unique(as.character(expression$SampleID))))
  }
  # wide matrix / data.frame: rows = genes, cols = RNA samples. Columns MUST be
  # matched by name against the epiPortrait sample IDs; positional matching is
  # unsafe for an independent RNA experiment. If RNA sample IDs differ
  # from ChIP colnames, require the long-format input.
  m <- as.matrix(expression)
  # Validate gene row names and the numeric matrix before aggregation.
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
  # Unmatched RNA sample columns are rejected rather than silently dropped.
  # A subset of the ChIP samples is allowed, but unknown RNA samples
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
#'   RNA-seq expression prioritization. Ordinal prioritization is explicit:
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
#' @param nearest_tss_cutoff_bp Numeric. Distance (bp) below which a
#'   \code{nearest_tss} relation is treated as a strong PROXIMAL candidate and
#'   ranked above \code{gene_body_overlap}/distal links; a currently-nearest
#'   but FAR TSS (beyond this cutoff) is downgraded to the lowest evidence tier.
#'   This prevents a "nearest but remote" gene (e.g. 500 kb away) from ranking
#'   equal to \code{promoter_overlap} evidence. Defaults to the option
#'   \code{epiPortrait.nearest_tss_cutoff_bp}, or 10000 if unset. Pass the SAME
#'   value to \code{annotate_epi_domains()} (which stores \code{best_tier} in
#'   the dedup table) so the stored and recomputed tiers cannot disagree.
#' @param unique_genes Logical. Collapse to unique genes (default TRUE),
#'   retaining the highest-priority domain association per gene.
#' @param max_per_domain Integer or NULL. When set, keep at most this many
#'   genes per domain, selected by \code{candidate_priority}. This is the
#'   transparent alternative to permutation-based background control for
#'   functional enrichment: wide domains would otherwise link more genes and
#'   inflate their opportunity to be selected. NULL (default) keeps all
#'   domain-gene pairs, preserving the original behaviour.
#' @param rank_by Character. Primary ordering key for \code{candidate_priority}.
#'   \itemize{
#'     \item \code{"tier"} (default): evidence hierarchy first (see Details),
#'           expression as tie-breaker — the original behaviour.
#'     \item \code{"bedpe_score"}: descending 3D contact strength
#'           (\code{bedpe_contact_score}) first, then the evidence tier and
#'           expression as tie-breakers. Requires that
#'           \code{annotate_epi_domains()} was run with
#'           \code{bedpe_score_col}; pairs without a usable score fall back to
#'           the tier ordering. This implements contact-strength-first target
#'           prioritization (dominant-loop style).
#'   }
#' @return A data.frame: domain_id, gene_id, gene_symbol, relation_types,
#'   evidence_tier (0-4 ordinal evidence strength), bedpe_supported,
#'   bedpe_support_count, bedpe_contact_score,
#'   nearest_tss_distance_bp (strand-aware SIGNED distance, matching
#'   \code{annotation_summary}; NA when the pair has no nearest-TSS evidence),
#'   expression_value, expression_status, expression_rank (contiguous rank over
#'   the returned table), candidate_priority.
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
                             min_expression = NULL, unique_genes = TRUE,
                             nearest_tss_cutoff_bp =
                               getOption("epiPortrait.nearest_tss_cutoff_bp",
                                         10000),
                             max_per_domain = NULL,
                             rank_by = c("tier", "bedpe_score")) {
  expression_priority <- match.arg(expression_priority)
  rank_by <- match.arg(rank_by)
  if (!is.null(max_per_domain)) {
    if (length(max_per_domain) != 1L || !is.numeric(max_per_domain) ||
        !is.finite(max_per_domain) || max_per_domain < 1) {
      stop("max_per_domain must be a single finite number >= 1, or NULL.")
    }
    max_per_domain <- as.integer(max_per_domain)
  }
  # ---- input validation -----------------------------------------------------
  if (!is.null(domains)) {
    if (!is.logical(domains) || length(domains) != nrow(se) ||
        any(is.na(domains))) {
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
  # Resolve min_expression semantics from the declared expression_type.
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
  # Vectorised with split (single pass) instead of a per-pair rbind loop,
  # so it scales to 10^5-10^6 links.
  sp <- split(seq_len(nrow(links)), key)
  rel_collapse <- function(field) unlist(lapply(sp, function(i) {
    paste(sort(unique(links[[field]][i])), collapse = ";")
  }))
  relation_types <- rel_collapse("relation_type")
  dom_part  <- sub("\\|.*$", "", names(sp))
  gene_part <- sub("^[^\\|]*\\|", "", names(sp))
  symbol_of <- function(i) {
    s <- links$gene_symbol[i]; s[!is.na(s)][1]
  }
  bedpe_of <- function(i) {
    b <- links$evidence_source[i] == "bedpe"
    if (!any(b)) return(0L)
    length(unique(stats::na.omit(links$bedpe_record_id[i][b])))
  }
  # pair-level contact score: sum over UNIQUE supporting records (0 when the
  # pair has no BEDPE evidence; NA when evidence exists but scores are absent,
  # e.g. objects annotated without bedpe_score_col).
  has_score_col <- "contact_score" %in% colnames(links)
  bedpe_score_of <- function(i) {
    b <- which(links$evidence_source[i] == "bedpe" &
                 !is.na(links$bedpe_record_id[i]))
    if (!has_score_col || length(b) == 0) return(0)
    s <- links$contact_score[i][b]
    s <- s[!duplicated(links$bedpe_record_id[i][b])]
    if (all(is.na(s))) return(NA_real_)
    sum(s, na.rm = TRUE)
  }
  dist_value <- function(i) {
    d <- stats::na.omit(links$distance_to_tss_bp[i])
    if (length(d) == 0) return(NA_real_)
    # A pair has at most one nearest_tss row; keep the raw SIGNED distance
    # (same convention as rowData/annotation_summary).
    d <- d[which.min(abs(d))]
    if (!is.finite(d)) NA_real_ else d
  }
  out <- data.frame(
    domain_id = dom_part,
    gene_id   = gene_part,
    gene_symbol = vapply(sp, symbol_of, character(1)),
    relation_types = relation_types,
    bedpe_supported = vapply(sp, function(i) any(links$evidence_source[i] == "bedpe"), logical(1)),
    # Count unique BEDPE records; a single record may match via both
    # directions -> two rows; must not be double counted).
    bedpe_support_count = vapply(sp, bedpe_of, integer(1)),
    bedpe_contact_score = vapply(sp, bedpe_score_of, numeric(1)),
    # Use the closest absolute TSS distance across all linear links of
    # the pair (raw SIGNED value; see dist_value). Only nearest_tss rows carry
    # distance_to_tss_bp; for pairs without any nearest_tss row this stays NA.
    nearest_tss_distance_bp = vapply(sp, dist_value, numeric(1)),
    stringsAsFactors = FALSE)
  # NA symbol: fall back to any symbol lookup across all rows
  na_sym <- which(is.na(out$gene_symbol))
  if (length(na_sym) > 0) {
    sym_map <- stats::setNames(links$gene_symbol, links$gene_id)
    out$gene_symbol[na_sym] <- sym_map[out$gene_id[na_sym]]
  }
  rownames(out) <- NULL

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
        # expression_status is populated only when a threshold is in effect.
        # high_expression_first needs no threshold and leaves status NA while
        # still ranking by expression_value.
        if (!is.null(threshold)) {
          out$expression_status <- ifelse(is.na(out$expression_value), NA_character_,
                                          ifelse(out$expression_value >= threshold,
                                                 "Expressed", "Not_detected"))
        }
      }
    }
  }

  # Ordinal prioritization uses the evidence hierarchy first,
  # then expression as secondary. Use pmax() so the strongest evidence is NEVER
  # overwritten by a weaker relation on the same domain-gene pair.
  # A distal nearest TSS must not rank equally to a
  # promoter-overlapping one. Promoter overlap (TSS inside/at the domain) is now
  # the strongest linear evidence; nearest_tss is graded by proximity against
  # nearest_tss_cutoff_bp (proximal ranks above gene-body links, distal drops to
  # the lowest tier).
  has_promoter  <- grepl("promoter_overlap", out$relation_types)
  has_bedpe     <- grepl("bedpe_promoter_contact", out$relation_types)
  has_nearest   <- grepl("nearest_tss", out$relation_types)
  has_body      <- grepl("gene_body_overlap|fully_contained", out$relation_types)
  nearest_dist  <- abs(out$nearest_tss_distance_bp)
  near_nearest  <- has_nearest & is.finite(nearest_dist) &
                   nearest_dist <= nearest_tss_cutoff_bp
  # Note: "far nearest-TSS" (has_nearest & !near_nearest) needs no explicit
  # variable - it is intentionally absorbed by tier 0 (order_score stays 0), so
  # a distal nearest gene ranks below gene-body links rather than equal to them.

  order_score <- rep(0L, nrow(out))
  # Tier 4: promoter overlap (strongest linear evidence)
  order_score <- pmax(order_score, ifelse(has_promoter, 4L, 0L))
  # Tier 3: BEDPE promoter contact (3D evidence)
  order_score <- pmax(order_score, ifelse(has_bedpe, 3L, 0L))
  # Tier 2: proximal nearest-TSS (strong linear proximity)
  order_score <- pmax(order_score, ifelse(near_nearest, 2L, 0L))
  # Tier 1: gene-body overlap / fully contained
  order_score <- pmax(order_score, ifelse(has_body, 1L, 0L))
  # Tier 0: far nearest-TSS (nearest but remote) and anything else
  # secondary: expression
  expr_tie <- rep(0, nrow(out))
  if (expression_priority == "expressed_first") {
    expr_tie[out$expression_status == "Expressed"] <- 1
  } else if (expression_priority == "high_expression_first") {
    expr_tie <- ifelse(is.na(out$expression_value), 0, out$expression_value)
  }
  # order() returns a permutation; convert it to per-row ranks.
  # rank_by = "bedpe_score": contact STRENGTH is the primary key (dominant-
  # loop style prioritization). Pairs without a usable score (NA, i.e. BEDPE
  # evidence present but no score column / all-NA scores) are treated like
  # pairs with no contact (0) and fall through to the evidence-tier ordering,
  # exactly as documented. Using -Inf here would rank unscored contacts BELOW
  # contacts with no evidence at all, which is not what "fall back to tier"
  # means.
  if (rank_by == "bedpe_score") {
    sc_key <- out$bedpe_contact_score
    sc_key[is.na(sc_key)] <- 0
    ord <- order(sc_key, order_score, expr_tie, decreasing = TRUE,
                 na.last = TRUE)
  } else {
    ord <- order(order_score, expr_tie, decreasing = TRUE, na.last = TRUE)
  }
  priority <- integer(nrow(out))
  priority[ord] <- seq_along(ord)
  out$candidate_priority <- priority
  # Ordinal evidence tier (0 = far nearest-TSS / other, 4 = promoter overlap)
  # exposed so downstream interpretation (enrichment, ranking) does not have
  # to re-parse relation_types.
  out$evidence_tier <- as.integer(order_score)
  # Optional per-domain cap: wide domains naturally link more genes; keeping
  # only the top-k pairs per domain bounds that opportunity without changing
  # the evidence hierarchy itself.
  if (!is.null(max_per_domain)) {
    out <- out[order(out$domain_id, out$candidate_priority), , drop = FALSE]
    dom_rank <- stats::ave(seq_len(nrow(out)), out$domain_id,
                           FUN = seq_along)
    out <- out[dom_rank <= max_per_domain, , drop = FALSE]
  }
  # When collapsing to unique genes, keep the best-supported domain-gene
  # association (highest priority), not the first row.
  if (unique_genes) {
    out <- out[order(out$candidate_priority, out$domain_id), ]
    out <- out[!duplicated(out$gene_id), , drop = FALSE]
  }
  # expression_rank: contiguous descending rank over the FINAL returned table.
  # Recomputed after max_per_domain / unique_genes filtering so the ranks have
  # no gaps and stay comparable.
  if (any(!is.na(out$expression_value))) {
    expr_ord <- order(-out$expression_value, na.last = TRUE)
    out$expression_rank <- match(seq_len(nrow(out)), expr_ord)
  }
  # candidate_priority is a contiguous ordinal rank of the final
  # output table (a gap can appear when unique_genes drops intermediate rows).
  out <- out[order(out$candidate_priority, out$domain_id), ]
  out$candidate_priority <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}
