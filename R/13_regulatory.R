#' Annotate Shape-Shifted Regions with Regulatory Element Catalogs
#'
#' @description Overlaps shape-shifted regions with regulatory element
#' annotations (ENCODE cCREs, Roadmap chromatin states, or custom catalogs)
#' and computes per-class fold enrichment over the genomic background.
#' This answers the question: "Are Concentration events enriched in enhancers?"
#'
#' @param shift_res A data.frame from \code{classify_shape_shift()} containing
#'   genomic coordinates and Shape_Class labels.
#' @param catalog A \code{GRanges} object with regulatory annotations, or a
#'   character path to a BED file, or \code{"encode"} to auto-fetch the
#'   ENCODE SCREEN cCREs via AnnotationHub.
#' @param genome Character. Reference genome (e.g., \code{"hg38"}). Only
#'   required when \code{catalog = "encode"}.
#' @param background A \code{GRanges} of background regions for enrichment
#'   calculation. If \code{NULL} (default), uses all peaks in \code{shift_res}.
#' @param fdr_cutoff Numeric. Only test peaks with adj.P.Val below this
#'   threshold for enrichment (default: 0.05).
#'
#' @return A list with components:
#'   \item{overlap_table}{Per-class per-category overlap counts.}
#'   \item{enrichment}{Per-class per-category fold enrichment (observed/expected).}
#'   \item{catalog_summary}{Summary of the regulatory catalog used.}
#'
#' @import GenomicRanges
#' @export
annotate_regulatory <- function(shift_res,
                                 catalog,
                                 genome = "hg38",
                                 background = NULL,
                                 fdr_cutoff = 0.05) {

  # ---- 1. Resolve regulatory catalog ----------------------------------------
  if (is.character(catalog) && length(catalog) == 1 && catalog == "encode") {
    catalog_gr <- .fetch_encode_ccres(genome)
    catalog_name <- sprintf("ENCODE SCREEN cCREs (%s)", genome)
  } else if (is.character(catalog) && length(catalog) == 1 && file.exists(catalog)) {
    catalog_gr <- rtracklayer::import(catalog)
    catalog_name <- basename(catalog)
  } else if (is(catalog, "GRanges")) {
    catalog_gr <- catalog
    catalog_name <- "user-provided catalog"
  } else {
    stop("catalog must be a GRanges, a BED file path, or 'encode'.")
  }

  # Ensure catalog has a 'type' or 'name' column; fall back to generic labels
  type_col <- if ("type" %in% colnames(mcols(catalog_gr))) "type"
              else if ("name" %in% colnames(mcols(catalog_gr))) "name"
              else NULL
  if (is.null(type_col)) {
    mcols(catalog_gr)$type <- "regulatory_element"
    type_col <- "type"
  }
  reg_types <- unique(mcols(catalog_gr)[[type_col]])
  message(sprintf("Regulatory catalog: %s (%d elements, %d types)",
                  catalog_name, length(catalog_gr), length(reg_types)))

  # ---- 2. Build query GRanges ------------------------------------------------
  sig <- shift_res[shift_res$adj.P.Val <= fdr_cutoff, , drop = FALSE]
  if (nrow(sig) == 0) {
    warning("No significant peaks at the current FDR cutoff.")
    return(list(overlap_table = NULL, enrichment = NULL,
                catalog_summary = catalog_name))
  }

  query_gr <- GRanges(
    seqnames = sig$seqnames,
    ranges   = IRanges(start = sig$start, end = sig$end)
  )
  mcols(query_gr)$Shape_Class <- sig$Shape_Class

  if (is.null(background)) {
    bg_gr <- GRanges(
      seqnames = shift_res$seqnames,
      ranges   = IRanges(start = shift_res$start, end = shift_res$end)
    )
  } else {
    bg_gr <- background
  }

  # ---- 3. Compute per-class overlap and enrichment --------------------------
  classes <- setdiff(unique(sig$Shape_Class), "Stable")
  if (length(classes) == 0) classes <- unique(sig$Shape_Class)

  overlap_list <- lapply(classes, function(cl) {
    cl_gr <- query_gr[mcols(query_gr)$Shape_Class == cl]
    hits <- findOverlaps(cl_gr, catalog_gr)
    if (length(hits) == 0) return(NULL)
    hit_types <- mcols(catalog_gr)[subjectHits(hits), type_col]
    table(factor(hit_types, levels = reg_types))
  })
  names(overlap_list) <- classes
  overlap_list <- overlap_list[!vapply(overlap_list, is.null, logical(1))]

  if (length(overlap_list) == 0) {
    warning("No overlaps found between significant peaks and regulatory catalog.")
    return(list(overlap_table = NULL, enrichment = NULL,
                catalog_summary = catalog_name))
  }

  overlap_mat <- do.call(rbind, overlap_list)

  # Compute background overlap rates for expected counts
  bg_hits <- findOverlaps(bg_gr, catalog_gr)
  bg_type_counts <- integer(length(reg_types))
  names(bg_type_counts) <- reg_types
  if (length(bg_hits) > 0) {
    bg_hit_types <- table(mcols(catalog_gr)[subjectHits(bg_hits), type_col])
    for (nm in names(bg_hit_types)) {
      bg_type_counts[nm] <- bg_hit_types[nm]
    }
  }
  total_bg <- length(bg_gr)

  enrichment_mat <- overlap_mat
  for (i in seq_len(nrow(overlap_mat))) {
    cl <- rownames(overlap_mat)[i]
    cl_size <- sum(sig$Shape_Class == cl)
    for (j in seq_len(ncol(overlap_mat))) {
      cat_name <- colnames(overlap_mat)[j]
      bg_rate <- bg_type_counts[cat_name] / total_bg
      expected <- cl_size * max(bg_rate, 1 / total_bg)
      enrichment_mat[i, j] <- overlap_mat[i, j] / max(expected, 0.01)
    }
  }

  message(sprintf("Regulatory overlap analysis complete for %d classes.", nrow(overlap_mat)))

  list(
    overlap_table = as.data.frame(overlap_mat),
    enrichment    = as.data.frame(enrichment_mat),
    catalog_summary = catalog_name
  )
}


# ---- Internal: fetch ENCODE SCREEN cCREs via AnnotationHub ------------------
.fetch_encode_ccres <- function(genome = "hg38") {
  if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
    stop("Package 'AnnotationHub' is required to fetch ENCODE data. ",
         "Install with BiocManager::install('AnnotationHub'), or provide a ",
         "local BED file as the 'catalog' argument.")
  }

  # Map genome to ENCODE SCREEN AH record
  genome_map <- list(
    hg38 = c("AH95512", "ENCODE cCREs hg38 v3"),
    hg19 = c("AH95511", "ENCODE cCREs hg19 v3"),
    mm10 = c("AH95513", "ENCODE cCREs mm10 v3")
  )
  if (!genome %in% names(genome_map)) {
    stop(sprintf("ENCODE cCREs not available for '%s'. Provide a local BED file.", genome))
  }

  ah_id <- genome_map[[genome]][1]
  ah <- AnnotationHub::AnnotationHub()
  catalog_gr <- ah[[ah_id]]

  # Standardise type column: ENCODE uses 'label' or 'type'
  if ("label" %in% colnames(mcols(catalog_gr))) {
    mcols(catalog_gr)$type <- mcols(catalog_gr)$label
  }
  message(sprintf("Fetched %s (%d elements) from AnnotationHub.",
                  genome_map[[genome]][2], length(catalog_gr)))

  return(catalog_gr)
}
