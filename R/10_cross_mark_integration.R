# Cross-mark integration and per-sample transition helpers.
#
# Two design gaps closed (GSE251898 two-axis benchmark, 2026-08-26):
#   * A single epiPortrait object carries ONE mark's axes (e.g. H3K4me3
#     Breadth or H3K27ac Intensity). The full two-axis picture (intensity
#     on the enhancer mark, width on the promoter mark) requires aligning
#     two objects whose universes usually DIFFER (different stitching /
#     promoter filtering). integrate_cross_mark() overlays the second
#     mark's per-sample calls onto the first mark's domains.
#   * Time-course workflows called in per_sample mode (1 sample per
#     condition) have no transition helper: compare_superdomains() /
#     compare_superdomain_classes() require per-GROUP columns.
#     transition_matrix_per_sample() fills this gap.


#' Integrate a Second Mark's Super-Domain Call onto Another Object
#'
#' @description Overlays the per-sample super-domain calls of object B
#'   (\code{se_b}) onto the shared domains of object A (\code{se_a}). The two
#'   objects may have DIFFERENT domain universes (e.g. a promoter-filtered +
#'   12.5 kb-stitched H3K27ac universe vs a 5 kb-stitched H3K4me3 universe);
#'   each \code{se_a} domain is classified from the \code{se_b} calls of the
#'   \code{se_b} domains that overlap it. This enables the two-axis design:
#'   e.g. whether promoter-proximal ARBS that an enhancer-mark pipeline
#'   deliberately excludes are captured by the promoter-mark WIDTH axis
#'   ("boundary rescue"), or the orthogonal Intensity x Width combinations.
#'
#' @param se_a A \code{SummarizedExperiment} with per-sample calls
#'   (the LEFT object whose domains are annotated).
#' @param se_b A \code{SummarizedExperiment} with per-sample calls
#'   (the RIGHT object whose evidence is overlaid).
#' @param call_fmt_a,call_fmt_b Character. \code{sprintf} formats producing the
#'   per-sample call column names, e.g. \code{"Breadth_Call__\%s"} /
#'   \code{"Intensity_Call__\%s"}. At least the union of used time points must
#'   exist in each object.
#' @param mark_a,mark_b Character or NULL. Short labels (e.g. "Width",
#'   "Intensity") prepended to the new rowData column names. When NULL they
#'   are derived from the call column prefix (text before \code{"_Call__"}).
#' @param timepoints Character or NULL. Time points to integrate. When NULL,
#'   the intersection of the column-derived time points of both objects is
#'   used (falling back to \code{colData$TimePoint}).
#' @param aggregate_ov Character. How to collapse MULTIPLE overlapping
#'   \code{se_b} domains onto one \code{se_a} domain:
#'   \itemize{
#'     \item \code{"any_super"} (default): Super if ANY overlapping \code{se_b}
#'           domain is super; otherwise Typical if any overlapping domain is
#'           called (non-NA); Uncertain if all overlapping calls are NA; NA
#'           (outside-universe) when there is no overlap at all.
#'     \item \code{"majority"}: strict majority (> half) of the overlapping
#'           non-NA \code{se_b} calls are super.
#'   }
#' @param min_overlap_bp Numeric. Minimum base pairs of overlap required for a
#'   \code{se_b} domain to count (default 0 = any overlap within 1 bp).
#' @return \code{se_a} with per-time-point rowData columns:
#'   \itemize{
#'     \item \code{<mark_b>_Call__<tp>}: aggregated \code{se_b} call
#'     \item \code{<mark_b>_NOverlaps__<tp>}: number of overlapping \code{se_b}
#'           domains (0 = not covered by the \code{se_b} universe)
#'   }
#'   plus a long-format integration table stored in
#'   \code{metadata(se_a)$cross_mark_integration} and a provenance list in
#'   \code{metadata(se_a)$cross_mark_provenance}.
#' @import GenomicRanges
#' @importFrom utils combn
#' @examples
#' data(example_se)
#' ## both objects must carry per_sample calls; see examples of
#' ## transition_matrix_per_sample() for building a two-object setup.
#' @export
integrate_cross_mark <- function(se_a, se_b,
                                 call_fmt_a = "Breadth_Call__%s",
                                 call_fmt_b = "Intensity_Call__%s",
                                 mark_a = NULL, mark_b = NULL,
                                 timepoints = NULL,
                                 aggregate_ov = c("any_super", "majority"),
                                 min_overlap_bp = 0) {
  aggregate_ov <- match.arg(aggregate_ov)
  if (!inherits(se_a, "SummarizedExperiment") ||
      !inherits(se_b, "SummarizedExperiment")) {
    stop("se_a and se_b must be SummarizedExperiment objects.")
  }
  if (length(min_overlap_bp) != 1L || !is.numeric(min_overlap_bp) ||
      !is.finite(min_overlap_bp) || min_overlap_bp < 0) {
    stop("min_overlap_bp must be a finite non-negative number.")
  }

  .derive_tps <- function(se, fmt) {
    pre <- sub("_Call__%s$", "", fmt)
    cols <- colnames(rowData(se))
    grp <- cols[grepl(paste0("^", pre, "_Call__"), cols)]
    if (length(grp) > 0) {
      unique(sub(paste0("^", pre, "_Call__"), "", grp))
    } else if ("TimePoint" %in% colnames(colData(se))) {
      unique(as.character(colData(se)$TimePoint))
    } else character(0)
  }
  tps_a <- .derive_tps(se_a, call_fmt_a)
  tps_b <- .derive_tps(se_b, call_fmt_b)
  tps <- if (is.null(timepoints)) intersect(tps_a, tps_b) else timepoints
  if (length(tps) == 0) {
    stop(sprintf("No shared time points between the two objects. ",
                 "Object A has: %s; object B has: %s.",
                 paste(tps_a, collapse = ","), paste(tps_b, collapse = ",")))
  }
  missing_a <- tps[!sprintf(call_fmt_a, tps) %in% colnames(rowData(se_a))]
  missing_b <- tps[!sprintf(call_fmt_b, tps) %in% colnames(rowData(se_b))]
  if (length(missing_a) > 0 || length(missing_b) > 0) {
    stop(sprintf("Missing call columns: A: %s; B: %s.",
                 paste(sprintf(call_fmt_a, missing_a), collapse = ","),
                 paste(sprintf(call_fmt_b, missing_b), collapse = ",")))
  }

  if (is.null(mark_a)) mark_a <- sub("_Call__%s$", "", call_fmt_a)
  if (is.null(mark_b)) mark_b <- sub("_Call__%s$", "", call_fmt_b)
  if (identical(mark_a, mark_b)) {
    stop("se_a and se_b use the same call prefix (", mark_a,
         "). Pass explicit mark_b / a different call prefix so the overlaid ",
         "columns do not overwrite se_a's own calls.")
  }

  dom_a <- rowRanges(se_a)
  dom_b <- rowRanges(se_b)

  # Chromosome names alone cannot establish coordinate compatibility: hg19
  # and hg38 both commonly use chr1..chrY. Prefer the assembly recorded by
  # build_portrait_matrix(), falling back to a unique colData Genome value.
  .object_genome <- function(se) {
    g <- S4Vectors::metadata(se)$signal_contract$genome
    if (is.null(g) && "Genome" %in% colnames(colData(se))) {
      g <- unique(stats::na.omit(as.character(colData(se)$Genome)))
    }
    g <- unique(as.character(g))
    g[nzchar(g) & !is.na(g)]
  }
  genome_a <- .object_genome(se_a)
  genome_b <- .object_genome(se_b)
  if (length(genome_a) == 1L && length(genome_b) == 1L &&
      !identical(genome_a, genome_b)) {
    stop("se_a and se_b use different genome assemblies ('", genome_a,
         "' vs '", genome_b, "'). Lift coordinates to one assembly before ",
         "cross-mark integration.", call. = FALSE)
  }
  sl_a <- unique(as.character(GenomicRanges::seqnames(dom_a)))
  sl_b <- unique(as.character(GenomicRanges::seqnames(dom_b)))
  shared <- intersect(sl_a, sl_b)
  if (length(shared) == 0) {
    stop("se_a and se_b share no seqlevels (chromosome naming mismatch?).")
  }
  if (length(shared) < length(sl_a)) {
    warning(sprintf("Integrating over %d/%d se_a seqlevels present in se_b.",
                    length(shared), length(sl_a)), call. = FALSE)
  }

  # Source call-value prefix comes from call_fmt_b (the REAL values in se_b's
  # rowData), NOT from mark_b. mark_b is only the OUTPUT label. Deriving from
  # mark_b would silently produce "H3K4me3_Super_Element" when the data carry
  # "Breadth_Super_Element" (review 2026: P1 cross-mark prefix/label mixing).
  source_prefix_b <- sub("_Call__%s$", "", call_fmt_b)
  super_b <- paste0(source_prefix_b, "_Super_Element")
  typ_b   <- paste0(source_prefix_b, "_Typical")
  out_super_b <- paste0(mark_b, "_Super_Element")
  out_typ_b   <- paste0(mark_b, "_Typical")

  # Histone-mark / accessibility domains are genomic intervals rather than
  # strand-specific transcript features. Explicitly ignore strand so GRanges
  # imported from different sources behave consistently with BED intervals.
  hits <- GenomicRanges::findOverlaps(dom_a, dom_b, ignore.strand = TRUE)
  if (min_overlap_bp > 0) {
    ovw <- GenomicRanges::width(
      GenomicRanges::pintersect(dom_a[S4Vectors::queryHits(hits)],
                                dom_b[S4Vectors::subjectHits(hits)],
                                ignore.strand = TRUE))
    hits <- hits[ovw >= min_overlap_bp]
  }

  long_rows <- list()
  for (tp in tps) {
    cl_a <- rowData(se_a)[[sprintf(call_fmt_a, tp)]]
    cl_b <- rowData(se_b)[[sprintf(call_fmt_b, tp)]]
    is_sup_b <- !is.na(cl_b) & cl_b == super_b
    n_ov <- integer(length(dom_a))
    agg <- rep(NA_character_, length(dom_a))  # NA = no se_b coverage

    if (length(hits) > 0) {
      qh <- S4Vectors::queryHits(hits)
      sh <- S4Vectors::subjectHits(hits)
      n_ov <- tabulate(qh, nbins = length(dom_a))
      spl <- split(seq_along(sh), qh)
      agg <- vapply(seq_len(length(dom_a)), function(i) {
        idx <- spl[[as.character(i)]]
        if (is.null(idx)) return(NA_character_)
        cb <- cl_b[sh[idx]]
        n_sup <- sum(cb == super_b, na.rm = TRUE)
        n_typ <- sum(cb == typ_b, na.rm = TRUE)
        if (n_sup + n_typ == 0) return("Uncertain")
        if (aggregate_ov == "any_super") {
          if (n_sup > 0) out_super_b else {
            # typical evidence exists; label with mark_b nomenclature
            out_typ_b
          }
        } else {
          if (n_sup / (n_sup + n_typ) > 0.5) out_super_b else out_typ_b
        }
      }, character(1))
    }

    rowData(se_a)[[sprintf("%s_Call__%s", mark_b, tp)]] <- unname(agg)
    rowData(se_a)[[sprintf("%s_NOverlaps__%s", mark_b, tp)]] <- n_ov
    long_rows[[tp]] <- data.frame(
      Domain_ID = rownames(se_a),
      seqnames = as.character(GenomicRanges::seqnames(dom_a)),
      start = GenomicRanges::start(dom_a),
      end = GenomicRanges::end(dom_a),
      TimePoint = tp,
      Call_A = cl_a,
      Call_B = unname(agg),
      NOverlaps_B = n_ov,
      stringsAsFactors = FALSE)
  }

  long <- do.call(rbind, long_rows)
  S4Vectors::metadata(se_a)$cross_mark_integration <- long
  S4Vectors::metadata(se_a)$cross_mark_provenance <- list(
    se_a_mark = mark_a, se_b_mark = mark_b,
    call_fmt_a = call_fmt_a, call_fmt_b = call_fmt_b,
    timepoints = tps,
    aggregate_ov = aggregate_ov,
    min_overlap_bp = min_overlap_bp,
    n_se_b_domains = length(dom_b),
    n_overlap_pairs = length(hits),
    note = "NA call = no overlapping se_b domain (outside se_b universe); 'Uncertain' = overlap exists but all overlapping se_b calls are NA.")
  se_a
}


#' Pairwise Per-Sample Transition Table Across Time Points
#'
#' @description Computes super-domain / class transitions between per-sample
#'   call columns (mode \code{"per_sample"}, one sample per time point), which
#'   \code{compare_superdomains()} / \code{compare_superdomain_classes()} do
#'   not support (they require per-GROUP columns). Adds one transition column
#'   per (ref, target) time-point pair and a labelled counts table to
#'   \code{metadata(se)$transitions}.
#'
#' @param se A \code{SummarizedExperiment} carrying per-sample call columns.
#' @param feature Character or NULL. Feature prefix for the call columns
#'   (e.g. \code{"Breadth"}, \code{"Intensity"}, or \code{"Combined_Class"}).
#'   Ignored when \code{call_fmt} is given.
#' @param call_fmt Character or NULL. \code{sprintf} format of the per-sample
#'   call columns, e.g. \code{"Breadth_Call__\%s"}. Defaults to
#'   \code{sprintf("\%s_Call__\%s", feature, "\%s")}.
#' @param timepoints Character or NULL. Ordered time points. When NULL,
#'   derived from the call column suffixes, preserving their discovery order
#'   (use \code{colData$TimePoint} order when present).
#' @param ref Character or NULL. Reference time point for transitions; when
#'   given, all other time points are targets. When NULL, all ordered pairs
#'   (a, b) with a before b are computed.
#' @param verbose Logical. Print the transition counts (default TRUE).
#' @return \code{se} with rowData columns
#'   \code{<prefix>_SampleTransition__<ref>_vs_<target>} (labels:
#'   \code{Persistent_<value>}, \code{<from>_to_<to>}, \code{Uncertain},
#'   matching \code{compare_superdomain_classes()} semantics) and
#'   \code{metadata(se)$transitions[[key]]} entries with the counts tables.
#' @importFrom stats setNames
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_sample", verbose = FALSE)
#' se <- transition_matrix_per_sample(se, feature = "Intensity",
#'                                    verbose = FALSE)
#' SummarizedExperiment::rowData(se)$Intensity_SampleTransition__C1_vs_T1
#' @export
transition_matrix_per_sample <- function(se, feature = NULL, call_fmt = NULL,
                                         timepoints = NULL, ref = NULL,
                                         verbose = TRUE) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment.")
  }
  if (is.null(call_fmt)) {
    if (is.null(feature)) {
      stop("Provide either `feature` or `call_fmt`.")
    }
    call_fmt <- sprintf("%s_Call__%%s", feature)
  }
  prefix <- sub("_Call__%s$", "", call_fmt)
  rd <- rowData(se)

  if (is.null(timepoints)) {
    cols <- colnames(rd)
    grp <- cols[grepl(paste0("^", prefix, "_Call__"), cols)]
    if (length(grp) == 0) {
      stop(sprintf("No columns matching '%s' found.",
                   sprintf("^%s_Call__", prefix)))
    }
    tps <- sub(paste0("^", prefix, "_Call__"), "", grp)
    if ("TimePoint" %in% colnames(colData(se))) {
      ord <- as.character(colData(se)$TimePoint)
      tps <- unique(c(intersect(ord, tps), setdiff(tps, ord)))
    }
    timepoints <- unique(tps)
  }
  timepoints <- unique(as.character(timepoints))
  missing <- timepoints[!sprintf(call_fmt, timepoints) %in% colnames(rd)]
  if (length(missing) > 0) {
    stop(sprintf("Missing per-sample call column(s): %s.",
                 paste(sprintf(call_fmt, missing), collapse = ", ")))
  }

  refs_targets <- if (!is.null(ref)) {
    if (!ref %in% timepoints) {
      stop("`ref` must be one of the time points.")
    }
    lapply(setdiff(timepoints, ref), function(t) c(ref, t))
  } else {
    cm <- combn(timepoints, 2)
    lapply(seq_len(ncol(cm)), function(i) cm[, i])
  }

  if (is.null(S4Vectors::metadata(se)$transitions)) {
    S4Vectors::metadata(se)$transitions <- list()
  }
  for (pair in refs_targets) {
    r_tp <- pair[1]; t_tp <- pair[2]
    r <- rd[[sprintf(call_fmt, r_tp)]]
    t <- rd[[sprintf(call_fmt, t_tp)]]
    if (length(r) != nrow(se) || length(t) != nrow(se)) {
      stop("Call columns must be length nrow(se).")
    }
    r_uncertain <- is.na(r) | r == "Uncertain"
    t_uncertain <- is.na(t) | t == "Uncertain"
    tr <- ifelse(r_uncertain | t_uncertain, "Uncertain",
                 ifelse(r == t, paste0("Persistent_", r),
                        paste0(r, "_to_", t)))
    col <- sprintf("%s_SampleTransition__%s_vs_%s", prefix, r_tp, t_tp)
    rowData(se)[[col]] <- tr
    key <- paste0(r_tp, "_vs_", t_tp)
    S4Vectors::metadata(se)$transitions[[key]] <- list(
      type = "per_sample_transition",
      ref = r_tp, target = t_tp,
      call_column_fmt = call_fmt,
      created_columns = col,
      counts = table(tr))
    if (verbose) {
      message(sprintf("[%s] %s -> %s", prefix, r_tp, t_tp))
      message(paste(names(table(tr)), table(tr), sep = "=", collapse = ", "))
    }
  }
  se
}
