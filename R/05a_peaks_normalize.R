#' Stitch Proximal Epigenetic Peaks
#'
#' @description Merges closely spaced peaks into larger continuous macro-domains based on
#' a user-defined stitching distance. This is crucial for analyzing Super-Elements
#' or broad heterochromatin spreading in epiPortrait.
#'
#' @param gr A GRanges object of peaks (typically the output from get_consensus_peaks).
#' @param stitch_distance Numeric. Maximum distance (in base pairs) between peaks
#' to be stitched together. For Super-Elements, 12500 (12.5 kb) is standard.
#' For standard broad marks, you might use 3000 to 5000.
#' @return A GRanges object of stitched domains, with a metadata column
#' 'Constituent_Peaks' recording how many original peaks were merged into each
#' domain. The call is recorded in \code{S4Vectors::metadata()} under
#' \code{stitch_provenance} (input/output counts, distance, timestamp), and is
#' propagated into \code{metadata(se)$stitch_provenance} by
#' \code{build_portrait_matrix()} when the stitched domains are used there.
#' @import GenomicRanges
#' @examples
#' gr <- GenomicRanges::GRanges("chr1",
#'        IRanges::IRanges(start = c(100, 10000), end = c(200, 10100)))
#' stitch_epi_peaks(gr, stitch_distance = 12500)
#' @export
stitch_epi_peaks <- function(gr, stitch_distance = 12500) {

  if (!inherits(gr, "GRanges")) stop("Input 'gr' must be a GRanges object.")

  if (length(gr) == 0) {
    warning("Input GRanges is empty. Returning as is.")
    return(gr)
  }

  if (!is.numeric(stitch_distance) || stitch_distance < 0) {
    stop("Parameter 'stitch_distance' must be a non-negative number.")
  }

  message(sprintf("Stitching peaks within %d bp of each other...", stitch_distance))

  # reduce() merges intervals whose gap is
  # STRICTLY LESS than min.gapwidth, so a gap exactly equal to stitch_distance
  # would not merge. Using min.gapwidth = stitch_distance + 1 makes the
  # documented "maximum distance to stitch" (gap <= stitch_distance) exact.
  stitched_gr <- GenomicRanges::reduce(gr, min.gapwidth = stitch_distance + 1)

  mcols(stitched_gr)$Constituent_Peaks <- countOverlaps(stitched_gr, gr)

  # Stitching determines the domain geometry used by every downstream layer,
  # so the call itself must be auditable: record what was merged and how far.
  S4Vectors::metadata(stitched_gr)$stitch_provenance <- list(
    stitch_distance_bp = stitch_distance,
    min_gapwidth_used = stitch_distance + 1,
    n_input_peaks = length(gr),
    n_output_domains = length(stitched_gr),
    rule = "GenomicRanges::reduce with gap <= stitch_distance merged",
    call = match.call(),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE))

  message(sprintf("Stitching complete: %d original peaks were stitched into %d continuous domains.",
                  length(gr), length(stitched_gr)))

  return(stitched_gr)
}


#' Filter Peaks by Genomic Annotations (Promoter / TSS Exclusion)
#'
#' @description Filters out peaks that overlap or are contained within specified
#' promoter-proximal regions. This is an explicit, opt-in universe-definition
#' step: the rest of the package never applies it automatically. For distal
#' enhancer marks (\code{get_mark_preset("H3K27ac")$exclude_promoter} is
#' \code{TRUE}), call it on the consensus peaks \emph{before}
#' \code{\link{stitch_epi_peaks}()}, which matches the order of ROSE's
#' \code{-t} option. Note that ROSE's own default is no TSS exclusion
#' (\code{-t 0}); set \code{upstream}/\code{downstream} (e.g. 2500 for a
#' ROSE-like exclusion) only when a promoter-excluded universe is intended, and
#' use a matching external reference.
#'
#' \strong{ROSE equivalence.} ROSE removes constituent peaks that are
#' \emph{contained within} a TSS +/- \code{tssWindow} zone and then, after
#' stitching, reverts any stitched region spanning more than two gene TSS
#' (see \code{\link{revert_multi_tss}}). To approximate ROSE: use a RefSeq
#' (\code{refGene}) TSS set (\code{\link{tss_from_refgene}} or
#' \code{\link{tss_from_rose}}), \code{mode = "contained"} and
#' \code{upstream = downstream = 2500}. The annotation source matters: the
#' built-in \code{"hg38"}/\code{"hg19"}/\code{"mm10"} shortcuts resolve to
#' UCSC \emph{knownGene}, not RefSeq, and give a different TSS set.
#'
#' @param gr A GRanges object of peaks (or consensus domains).
#' @param genome A character string ("hg38", "hg19", "mm10", etc.) or a TxDb
#'   object. Ignored when \code{tss} is supplied.
#' @param upstream Number of bp upstream of TSS to define the window
#'   (default: 2000).
#' @param downstream Number of bp downstream of TSS to define the window
#'   (default: 2000).
#' @param mode \code{"overlap"} (default, legacy behaviour) removes any peak
#'   that intersects a promoter/TSS window. \code{"contained"} removes only
#'   peaks fully contained within a window (ROSE semantics), so large peaks
#'   that merely cross a TSS are kept.
#' @param tss Optional GRanges of TSS positions (1 bp, strand-aware) or a TxDb.
#'   When supplied it overrides \code{genome} and lets the caller reuse an
#'   annotation source (for example RefSeq via \code{\link{tss_from_refgene}})
#'   across many samples.
#' @return A filtered GRanges object. Provenance (mode, window, source, counts)
#'   is stored in \code{S4Vectors::metadata(x)$promoter_exclusion}.
#' @import GenomicRanges
#' @examples
#' gr <- GenomicRanges::GRanges("chr1",
#'   IRanges::IRanges(c(850, 980, 5000), width = 100))
#' tss <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 1000),
#'                               strand = "+")
#' # "overlap" removes the peaks at 850 and 980 (both intersect TSS +/- 100)
#' length(filter_promoter_peaks(gr, tss = tss, upstream = 100,
#'                              downstream = 100, mode = "overlap"))
#' # "contained" keeps the 850 peak (it crosses the window edge)
#' length(filter_promoter_peaks(gr, tss = tss, upstream = 100,
#'                              downstream = 100, mode = "contained"))
#' @export
filter_promoter_peaks <- function(gr, genome = "hg38", upstream = 2000,
                                  downstream = 2000,
                                  mode = c("overlap", "contained"),
                                  tss = NULL) {
  mode <- match.arg(mode)
  if (length(gr) == 0) return(gr)
  if (!methods::is(gr, "GRanges")) {
    stop("gr must be a GRanges object.", call. = FALSE)
  }
  for (nm in c("upstream", "downstream")) {
    v <- get(nm)
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v < 0) {
      stop(sprintf("%s must be a finite non-negative number.", nm),
           call. = FALSE)
    }
  }

  src <- "genome shortcut"
  if (!is.null(tss)) {
    if (methods::is(tss, "GRanges")) {
      windows <- GenomicRanges::promoters(tss, upstream = upstream,
                                          downstream = downstream)
      src <- "user-supplied TSS GRanges"
    } else if (methods::is(tss, "TxDb")) {
      windows <- GenomicFeatures::promoters(tss, upstream = upstream,
                                            downstream = downstream)
      src <- "user-supplied TxDb"
    } else {
      stop("tss must be NULL, a GRanges of TSS positions, or a TxDb.",
           call. = FALSE)
    }
  } else {
    res <- .resolve_genome_resources(genome)
    if (is.null(res$txdb)) {
      stop("Promoter filtering requires a TxDb (genome = 'hg38'/'hg19'/'mm10' ",
           "or a TxDb object).", call. = FALSE)
    }
    if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
      stop("Please install 'GenomicFeatures'.")
    }
    # A chr-vs-1 naming mismatch would silently remove 0 promoter overlaps and
    # leave every peak in place; enforce seqlevel compatibility.
    .check_seqlevel_compatibility(gr, txdb = res$txdb, enforce = TRUE)
    windows <- GenomicFeatures::promoters(res$txdb, upstream = upstream,
                                          downstream = downstream)
    src <- as.character(genome)[1]
  }

  shared <- intersect(as.character(GenomeInfoDb::seqlevels(gr)),
                      as.character(GenomeInfoDb::seqlevels(windows)))
  if (length(shared) == 0L) {
    stop("No shared seqlevels between gr and the promoter/TSS windows. ",
         "Check the genome build and chromosome naming (e.g. chr1 vs 1).",
         call. = FALSE)
  }

  # "overlap" = any intersection; "contained" = the peak lies fully inside a
  # promoter/TSS window (ROSE semantics). findOverlaps(type=) is vectorised.
  type <- if (identical(mode, "contained")) "within" else "any"
  hits <- GenomicRanges::findOverlaps(gr, windows, type = type,
                                      ignore.strand = TRUE)
  removed <- sort(unique(S4Vectors::queryHits(hits)))
  filtered_gr <- gr[setdiff(seq_len(length(gr)), removed)]

  S4Vectors::metadata(filtered_gr)$promoter_exclusion <- list(
    mode = mode, upstream = upstream, downstream = downstream,
    source = src, n_input = length(gr), n_removed = length(removed),
    order = "applied to peaks BEFORE stitch_epi_peaks() (ROSE -t order)",
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE))

  message(sprintf(
    "Excluded %d peaks (%s; %s, window -%d/+%d bp). %d peaks remaining.",
    length(removed), mode, src, upstream, downstream, length(filtered_gr)))
  filtered_gr
}


#' Build a RefSeq (refGene) TSS set
#'
#' @description Returns a 1 bp, strand-aware transcription-start-site (TSS)
#' GRanges from the UCSC \code{refGene} (RefSeq) track. This is the annotation
#' family used by ROSE's bundled \code{<genome>_refseq.ucsc} tables, and is the
#' appropriate source when a promoter/TSS exclusion must be comparable to ROSE.
#' It is \strong{not} the UCSC \code{knownGene} set used by the built-in
#' \code{genome} shortcuts of \code{\link{filter_promoter_peaks}}, so the two
#' give different TSS sets.
#'
#' @param genome Genome assembly passed to
#'   \code{txdbmaker::makeTxDbFromUCSC} (default "hg38").
#' @param txdb Optional TxDb to use instead of downloading (must be RefSeq
#'   based for ROSE comparability). When supplied, \code{genome} is only used
#'   for provenance.
#' @return A GRanges of 1 bp TSS positions with a \code{transcript} column.
#' @examples
#' \donttest{
#' if (requireNamespace("txdbmaker", quietly = TRUE) &&
#'     interactive()) {
#'   tss <- tss_from_refgene("hg38")
#'   length(tss)
#' }
#' }
#' @export
tss_from_refgene <- function(genome = "hg38", txdb = NULL) {
  if (is.null(txdb)) {
    if (!requireNamespace("txdbmaker", quietly = TRUE)) {
      stop("tss_from_refgene() needs the 'txdbmaker' package to build a ",
           "refGene TxDb; alternatively supply txdb = <TxDb>.", call. = FALSE)
    }
    txdb <- tryCatch(
      txdbmaker::makeTxDbFromUCSC(genome = genome, tablename = "refGene"),
      error = function(e) stop(
        sprintf("Could not build a refGene TxDb for '%s': %s. ", genome,
                conditionMessage(e)),
        "Pass txdb = <TxDb>, or use tss_from_rose() with a local ROSE ",
        "annotation file instead.", call. = FALSE))
  } else if (!methods::is(txdb, "TxDb")) {
    stop("txdb must be a TxDb object or NULL.", call. = FALSE)
  }
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop("Please install 'GenomicFeatures'.")
  }
  tx <- GenomicFeatures::transcripts(txdb)
  tss <- GenomicRanges::promoters(tx, upstream = 0, downstream = 1)
  S4Vectors::metadata(tss)$tss_source <- list(
    source = "UCSC refGene (RefSeq)", genome = genome,
    definition = "transcript 5' end, 1 bp, strand-aware", n = length(tss))
  tss
}


#' Parse a ROSE RefSeq annotation file into a TSS set
#'
#' @description Reads ROSE's bundled UCSC refGene table
#' (\code{<genome>_refseq.ucsc}, columns \code{name, chrom, strand, txStart,
#' txEnd, ...}) and returns a 1 bp, strand-aware TSS GRanges using ROSE's
#' definition (positive strand: \code{txStart}; negative strand: \code{txEnd}).
#' Reusing the exact file that a ROSE run used guarantees an identical TSS set
#' (see \href{https://github.com/younglab/ROSE}{younglab/ROSE}).
#'
#' @param genome Genome label, used for provenance only.
#' @param file Path to a ROSE refGene annotation file (e.g.
#'   \code{"hg38_refseq.ucsc"}).
#' @return A GRanges of 1 bp TSS positions with \code{transcript} (and
#'   \code{gene} when \code{name2} is present) columns.
#' @examples
#' f <- tempfile(fileext = ".ucsc")
#' writeLines(c("#bin\tname\tchrom\tstrand\ttxStart\ttxEnd",
#'              "0\tNM_000001\tchr1\t+\t1000\t2000",
#'              "0\tNM_000002\tchr1\t-\t5000\t6000"), f)
#' tss_from_rose("hg38", f)
#' @export
tss_from_rose <- function(genome, file) {
  if (missing(genome) || !is.character(genome) || length(genome) != 1L) {
    stop("genome must be a single character string.", call. = FALSE)
  }
  if (missing(file) || !is.character(file) || length(file) != 1L ||
      !file.exists(file)) {
    stop("file must be a path to an existing ROSE refGene annotation ",
         "(e.g. 'hg38_refseq.ucsc').", call. = FALSE)
  }
  raw <- utils::read.delim(file, sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE, comment.char = "",
                           check.names = FALSE, quote = "")
  if (nrow(raw) > 0L && grepl("^#", as.character(raw[1, 1]))) {
    colnames(raw) <- sub("^#", "", as.character(raw[1, ]))
    raw <- raw[-1, , drop = FALSE]
  }
  req <- c("name", "chrom", "strand", "txStart", "txEnd")
  if (!all(req %in% colnames(raw))) {
    stop("The ROSE annotation must be a UCSC refGene table with columns ",
         "'name','chrom','strand','txStart','txEnd'.", call. = FALSE)
  }
  pos <- ifelse(raw$strand == "-", as.numeric(raw$txEnd),
                as.numeric(raw$txStart) + 1L)
  tss <- GenomicRanges::GRanges(
    seqnames = raw$chrom,
    ranges = IRanges::IRanges(pos, pos),
    strand = raw$strand,
    transcript = raw$name)
  if ("name2" %in% colnames(raw)) tss$gene <- raw$name2
  S4Vectors::metadata(tss)$tss_source <- list(
    source = "ROSE refGene table", genome = genome, file = file,
    definition = "txStart (+ strand) / txEnd (- strand)", n = length(tss))
  tss
}


#' Revert stitched regions that span more than two TSS (ROSE semantics)
#'
#' @description Implements the ROSE post-stitching safeguard: a stitched region
#' whose span contains the TSS (+/- \code{tss_span}) of more than \code{max_tss}
#' distinct genes is replaced by its original constituent peaks (ROSE records
#' these as \code{MULTIPLE_TSS}). This prevents gene-dense loci from being
#' merged into a single artefactual macro-domain. It is only relevant when TSS
#' exclusion is enabled; ROSE applies it only for \code{-t != 0}.
#'
#' @param gr Stitched GRanges (e.g. from \code{\link{stitch_epi_peaks}}).
#' @param tss TSS GRanges (e.g. from \code{\link{tss_from_refgene}} /
#'   \code{\link{tss_from_rose}}). A \code{gene} or \code{transcript} column is
#'   used to count distinct genes.
#' @param original The constituent GRanges from which \code{gr} was stitched;
#'   flagged regions are replaced by their overlapping constituents.
#' @param max_tss Maximum number of distinct TSS allowed per stitched region
#'   (default 2, i.e. revert when > 2 = 3 or more).
#' @param tss_span Half-window (bp) used only for the TSS-counting overlap
#'   (default 50, matching ROSE).
#' @return A GRanges with flagged regions reverted to constituents. Provenance
#'   is stored in \code{S4Vectors::metadata(x)$multi_tss_revert}.
#' @examples
#' orig <- GenomicRanges::GRanges("chr1",
#'   IRanges::IRanges(c(1000, 1100, 1900, 2800, 2900, 4100), width = 50))
#' stitched <- GenomicRanges::GRanges("chr1",
#'   IRanges::IRanges(c(1000, 4000), c(3000, 4200)))
#' tss <- GenomicRanges::GRanges("chr1",
#'   IRanges::IRanges(c(1500, 2000, 2500, 4500), width = 1))
#' tss$gene <- c("A", "B", "C", "D")
#' # the first stitched region spans 3 TSS (A, B, C) -> reverted
#' length(revert_multi_tss(stitched, tss, orig, max_tss = 2))
#' @export
revert_multi_tss <- function(gr, tss, original, max_tss = 2L, tss_span = 50L) {
  if (!methods::is(gr, "GRanges") || !methods::is(tss, "GRanges") ||
      !methods::is(original, "GRanges")) {
    stop("gr, tss and original must be GRanges objects.", call. = FALSE)
  }
  if (length(max_tss) != 1L || !is.numeric(max_tss) || !is.finite(max_tss) ||
      max_tss < 0 || max_tss != floor(max_tss)) {
    stop("max_tss must be a single non-negative integer.", call. = FALSE)
  }
  if (length(tss_span) != 1L || !is.numeric(tss_span) ||
      !is.finite(tss_span) || tss_span < 0) {
    stop("tss_span must be a single non-negative number.", call. = FALSE)
  }
  if (length(gr) == 0L) return(gr)

  tss_w <- GenomicRanges::promoters(tss, upstream = tss_span,
                                    downstream = tss_span)
  gene <- if (!is.null(tss$gene)) as.character(tss$gene)
          else if (!is.null(tss$transcript)) as.character(tss$transcript)
          else as.character(seq_len(length(tss)))
  hits <- GenomicRanges::findOverlaps(gr, tss_w, ignore.strand = TRUE)
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  n_gene <- integer(length(gr))
  if (length(qh) > 0L) {
    tab <- tapply(gene[sh], qh, function(x) length(unique(x)))
    n_gene[as.integer(names(tab))] <- as.integer(tab)
  }
  flag <- which(n_gene > max_tss)

  if (length(flag) == 0L) {
    S4Vectors::metadata(gr)$multi_tss_revert <- list(
      max_tss = max_tss, tss_span = tss_span, n_reverted = 0L)
    return(gr)
  }

  keep <- gr[-flag]
  names(keep) <- NULL
  reverted <- lapply(flag, function(i) {
    ov <- GenomicRanges::findOverlaps(original, gr[i], ignore.strand = TRUE)
    cons <- original[unique(S4Vectors::queryHits(ov))]
    names(cons) <- NULL
    cons
  })
  out <- do.call(c, c(list(keep), reverted))
  out <- GenomicRanges::sort(out)
  S4Vectors::metadata(out)$multi_tss_revert <- list(
    max_tss = max_tss, tss_span = tss_span, n_reverted = length(flag),
    definition = paste("stitched regions spanning > max_tss distinct TSS",
                       "were reverted to their constituent peaks"))
  message(sprintf(
    "Reverted %d stitched region(s) spanning > %d TSS; %d regions remain.",
    length(flag), max_tss, length(out)))
  out
}



#' Normalize Portrait Assays
#'
#' @description Optionally normalizes the Intensity assay while preserving
#'   SignalDispersion in its native bp-scale units. The recommended workflow
#'   uses quantitatively comparable BigWigs and \code{method = "None"}.
#'
#' @param se A SummarizedExperiment object from build_portrait_matrix().
#' @param method Normalization method. Options are:
#'   \itemize{
#'     \item \code{"None"} (default): Skips normalization. Recommended —
#'           epiPortrait expects BigWigs that are already quantitatively
#'           comparable (same pipeline, genome, normalization; e.g. CPM/RPGC/
#'           spike-in). Post-hoc rescaling risks removing genuine global
#'           biological shifts.
#'     \item \code{"TotalSignal"}: Scales libraries to the mean total signal of
#'           the analyzed domain set. WARNING: this assumes approximately
#'           conserved total signal across samples and may remove real global
#'           gains/losses (e.g. drug-induced chromatin loss).
#'     \item \code{"Quantile"}: Forces identical distributions across samples
#'           (uses limma). This is an explicitly warned sensitivity-analysis
#'           option, not a recommended between-condition normalization.
#'   }
#'   TMM is intentionally not offered because it is a count-composition method,
#'   whereas epiPortrait operates on continuous integrated BigWig signal.
#'   Row-wise Z-score scaling is not offered here: it is a display/clustering
#'   transform that would destroy the cross-domain magnitude ranking used for
#'   Super calling. Use \code{plot_portrait_pca()} or plotting layers for
#'   display-only scaling.
#' @return A SummarizedExperiment. The decision and, when applicable, scaling
#'   factors and before/after sample totals are stored in
#'   \code{metadata(se)$normalization}. A second post-hoc normalization of an
#'   object already marked as normalized is rejected.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' normalize_portrait(example_se, method = "None")
#' @export
normalize_portrait <- function(se, method = "None") {

  valid_methods <- c("None", "TotalSignal", "Quantile")
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose from: 'None', 'TotalSignal', or 'Quantile'.")
  }

  if (method == "None") {
    if (is.null(S4Vectors::metadata(se)$normalization)) {
      S4Vectors::metadata(se)$normalization <- list(
        method = "None",
        applied = FALSE,
        input_assay = .resolve_assay(se, "Intensity"),
        assumption = paste(
          "Input BigWigs are already quantitatively comparable",
          "(same processing and CPM/RPGC/spike-in normalization).")
      )
    }
    message("Method set to 'None'. Skipping normalization; input BigWigs are assumed to be quantitatively comparable.")
    return(se)
  }

  previous <- S4Vectors::metadata(se)$normalization
  if (!is.null(previous) && isTRUE(previous$applied)) {
    stop("Intensity is already marked as normalized by normalize_portrait() ",
         "using method = '", previous$method, "'. Return to the pre-normalized ",
         "object instead of applying a second post-hoc normalization.")
  }

  message(sprintf("Normalizing features using '%s' method...", method))

  intensity_assay <- .resolve_assay(se, "Intensity")
  dispersion_assay <- .resolve_assay(se, "SignalDispersion")
  int_mat <- assay(se, intensity_assay)
  disp_mat <- assay(se, dispersion_assay)
  input_totals <- stats::setNames(colSums(int_mat, na.rm = TRUE), colnames(se))
  scaling_factors <- NULL
  assumption <- NULL

  if (method == "TotalSignal") {
    warning(
      "TotalSignal normalization scales samples to the mean total signal of the ",
      "analyzed domain set. It assumes approximately conserved total signal and ",
      "may remove genuine global biological shifts (e.g. drug-induced chromatin ",
      "loss). Consider using pre-normalized BigWigs (method = 'None') for ",
      "between-condition comparison.", call. = FALSE
    )
    sample_sums <- colSums(int_mat, na.rm = TRUE)
    target_scale <- mean(sample_sums)
    scaling_factors <- target_scale / sample_sums

    bad <- !is.finite(scaling_factors)
    if (any(bad)) {
      warning(sprintf("%d sample(s) have zero total signal; scaling factors set to 1.",
                      sum(bad)), call. = FALSE)
      scaling_factors[bad] <- 1
    }

    norm_int <- sweep(int_mat, 2, scaling_factors, FUN = "*")
    colData(se)$ScalingFactor <- scaling_factors
    assumption <- paste(
      "The total signal over the analyzed domain universe is approximately",
      "conserved across samples.")

  } else if (method == "Quantile") {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Please install 'limma' to use Quantile normalization.")
    }
    warning(
      "Quantile normalization forces identical Intensity distributions across ",
      "samples and can erase genuine global biological shifts. Use it only as ",
      "an explicitly justified sensitivity analysis, not as the default for ",
      "between-condition comparisons.", call. = FALSE
    )

    norm_int <- limma::normalizeBetweenArrays(int_mat, method = "quantile")
    assumption <- paste(
      "Differences in marginal Intensity distributions are treated as",
      "technical rather than biological.")
  }

  assay(se, intensity_assay) <- norm_int

  # SignalDispersion is not rescaled. It is a spatial measure in bp that
  # is invariant under a uniform multiplicative rescaling of the signal
  # (x_i -> c*x_i leaves the weighted genomic SD unchanged). Applying library
  # scaling factors or quantile transforms would change its unit and break its
  # biological interpretation. Row-wise Z-score scaling is likewise not applied:
  # it would destroy the cross-domain magnitude ranking used for Super calling;
  # plotting / PCA layers apply their own display scaling.
  if (!identical(disp_mat, assay(se, dispersion_assay))) {
    stop("SignalDispersion was unexpectedly modified during normalization.")
  }

  S4Vectors::metadata(se)$normalization <- list(
    method = method,
    applied = TRUE,
    input_assay = intensity_assay,
    transformed_assay = intensity_assay,
    signal_type = "continuous integrated BigWig signal",
    scaling_factors = if (is.null(scaling_factors)) NULL else
      stats::setNames(as.numeric(scaling_factors), colnames(se)),
    input_sample_totals = input_totals,
    output_sample_totals = stats::setNames(
      colSums(norm_int, na.rm = TRUE), colnames(se)),
    assumption = assumption,
    SignalDispersion_modified = FALSE
  )

  message("Normalization complete! Intensity adjusted; SignalDispersion ",
          "(bp-scale spatial descriptor) is left unchanged.")
  return(se)
}

# Run an expression with an OPTIONAL local random seed. When `seed` is NULL
# (default) the expression is evaluated without touching the global RNG state;
# when an integer is supplied it is wrapped in withr::with_seed() so the seed
# is applied only for the duration of `expr` and the caller's RNG stream is
# restored afterwards. This keeps the package free of global set.seed() calls
# (BiocCheck), while still offering reproducible bootstrap / sampling when the
# user explicitly asks for it.
.with_opt_seed <- function(seed, expr) {
  if (is.null(seed)) {
    force(expr)
  } else {
    withr::with_seed(seed, force(expr))
  }
}

# ---- Internal: inflection point detection ----------------------------------
