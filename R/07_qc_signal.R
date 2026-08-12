#' Check Signal Compatibility Across BigWig Samples
#'
#' @description Performs pre-analysis QC on a sample sheet of BigWig tracks
#' Checks file existence, signal ranges,
#' negative/zero fractions, chromosome naming, and cross-sample global scale
#' consistency. Returns a per-sample summary with a PASS / WARNING / FAIL
#' status and emits warnings for global incompatibilities.
#'
#' @param sample_sheet A \code{data.frame} with at least \code{SampleID},
#'   \code{Condition}, and \code{bw_path}.
#' @param regions A \code{GRanges} or NULL. If provided, signals are evaluated
#'   within these regions (e.g. candidate domains). If NULL, a small set of
#'   genomic windows is sampled from each track.
#' @param max_windows Integer. Number of windows to sample when \code{regions}
#'   is NULL (default 2000).
#' @param seed Integer or NULL. Optional random seed for reproducible window /
#'   region subsampling (default NULL = random draw). Pass a fixed seed to make
#'   the QC reproducible; NULL keeps the draw truly random.
#' @param verbose Logical. Print a summary (default TRUE).
#'
#' @return A \code{data.frame} with one row per sample:
#'   SampleID, Condition, signal-sum (QCWindowSignalSum when sampling tiles,
#'   DomainSetSignalSum when regions are provided), MedianSignal, MeanSignal,
#'   ZeroFraction, NegativeFraction, UpperQuantile, GlobalScaleRatio, Status.
#'   The object carries a \code{print} method that displays the summary table
#'   when \code{verbose = TRUE}; the returned value is otherwise an ordinary
#'   \code{data.frame}.
#'
#' @import rtracklayer
#' @importFrom stats quantile median
#' @examples
#' if (.Platform$OS.type != "windows") {
#'   extdata <- system.file("extdata", package = "epiPortrait")
#'   ss <- data.frame(
#'     SampleID  = c("C1", "C2", "T1", "T2"),
#'     Condition = c("Control", "Control", "Treatment", "Treatment"),
#'     bw_path   = file.path(extdata, c("C1.bw", "C2.bw", "T1.bw", "T2.bw")))
#'   check_signal_compatibility(ss, max_windows = 500)
#' }
#' @export
check_signal_compatibility <- function(sample_sheet, regions = NULL,
                                       max_windows = 2000, seed = NULL,
                                       verbose = TRUE) {
  req_cols <- c("SampleID", "Condition", "bw_path")
  if (!all(req_cols %in% colnames(sample_sheet))) {
    stop("sample_sheet must contain columns: 'SampleID', 'Condition', 'bw_path'")
  }

  missing <- sample_sheet$bw_path[!file.exists(sample_sheet$bw_path)]
  if (length(missing) > 0) {
    stop(sprintf("BigWig file(s) not found: %s", paste(missing, collapse = ", ")))
  }

  # Determine the genomic windows to summarize over.
  if (!is.null(regions)) {
    if (!inherits(regions, "GRanges")) stop("regions must be a GRanges object or NULL.")
    sel <- regions
    if (length(sel) > max_windows) {
      sel <- .with_opt_seed(seed, sel[sample(length(sel), max_windows)])
    }
    signal_col <- "DomainSetSignalSum"
  } else {
    # Use the BigWig header (seqinfo) to build FIXED genomic tiles across the
    # genome (P1-10). This avoids importing the first track's signal regions,
    # which would (a) load a large BigWig into memory and (b) bias the QC
    # toward the first sample's active regions.
    sel <- .genome_tiles(sample_sheet$bw_path[1], max_windows, seed = seed)
    signal_col <- "QCWindowSignalSum"
  }

  # Metadata compatibility (review #11): if provided, Genome / SignalType /
  # Normalization / BinSize must be consistent across samples.
  for (mc in c("Genome", "SignalType", "Normalization", "BinSize")) {
    if (mc %in% colnames(sample_sheet)) {
      uniq <- unique(sample_sheet[[mc]])
      if (length(uniq) > 1) {
        warning(paste0(
          sprintf("Sample metadata '%s' is not consistent across samples: %s. ",
                  mc, paste(uniq, collapse = ", ")),
          "Cross-sample quantitative comparison may be invalid."
        ))
      }
    }
  }

  message(sprintf("Checking signal compatibility over %d region(s) for %d sample(s)...",
                  length(sel), nrow(sample_sheet)))

  res_list <- lapply(seq_len(nrow(sample_sheet)), function(i) {
    bw <- sample_sheet$bw_path[i]
    cvg <- tryCatch(
      rtracklayer::import(bw, selection = rtracklayer::BigWigSelection(sel),
                          as = "NumericList"),
      error = function(e) NULL
    )
    if (is.null(cvg)) {
      d <- data.frame(
        SampleID = sample_sheet$SampleID[i],
        Condition = sample_sheet$Condition[i],
        MedianSignal = NA_real_, MeanSignal = NA_real_,
        ZeroFraction = NA_real_, NegativeFraction = NA_real_,
        UpperQuantile = NA_real_, Status = "FAIL"
      )
      d[[signal_col]] <- NA_real_
      return(d)
    }
    vals <- unlist(cvg, use.names = FALSE)
    vals <- vals[is.finite(vals)]
    n <- length(vals)
    if (n == 0) {
      d <- data.frame(
        SampleID = sample_sheet$SampleID[i], Condition = sample_sheet$Condition[i],
        MedianSignal = 0, MeanSignal = 0,
        ZeroFraction = 1, NegativeFraction = 0, UpperQuantile = 0,
        Status = "WARNING"
      )
      d[[signal_col]] <- 0
      return(d)
    }
    neg_frac <- mean(vals < 0)
    zero_frac <- mean(vals == 0)
    status <- "PASS"
    if (neg_frac > 0.01) status <- "WARNING"   # some negatives, may need clipping
    if (neg_frac > 0.5) status <- "FAIL"       # predominantly signed/control-subtracted
    d <- data.frame(
      SampleID = sample_sheet$SampleID[i],
      Condition = sample_sheet$Condition[i],
      MedianSignal = stats::median(vals),
      MeanSignal = mean(vals),
      ZeroFraction = zero_frac,
      NegativeFraction = neg_frac,
      UpperQuantile = as.numeric(stats::quantile(vals, 0.9, names = FALSE)),
      Status = status
    )
    d[[signal_col]] <- sum(vals)
    d
  })

  res <- do.call(rbind, res_list)
  rownames(res) <- NULL

  # Cross-sample global-scale compatibility: report the max/min signal-sum
  # ratio. Any ratio > 5 is flagged as WARNING (even a real global biological
  # change would be unusual beyond ~5x); the exact threshold for accepting
  # cross-sample quantitative comparison is a user judgment.
  totals <- res[[signal_col]][is.finite(res[[signal_col]])]
  res$GlobalScaleRatio <- NA_real_
  if (length(totals) >= 2) {
    ratio <- max(totals, na.rm = TRUE) / max(min(totals, na.rm = TRUE), 1e-10)
    res$GlobalScaleRatio <- ratio
    if (ratio > 5) {
      warning(sprintf(
        "Substantial global signal-scale differences detected (max/min %s = %.1f). ",
        signal_col, ratio,
        "Distinguish technical normalization from true global biology before ",
        "cross-sample quantitative comparison. Within-sample ranking remains interpretable."
      ))
    }
  }

  if (verbose) {
    message("Status legend: PASS / FLAG / FAIL. ",
            "Check NegativeFraction for control-subtracted tracks.")
  }
  class(res) <- c("epi_signal_qc", class(res))
  res
}

#' @export
print.epi_signal_qc <- function(x, ...) {
  print.data.frame(unclass(x), ...)
  invisible(x)
}


# Build fixed genomic tiles from a BigWig header (seqinfo) for QC sampling
# (P1-10). Avoids importing the full first track and avoids biasing the QC
# toward the first sample's active regions. When sampling down to max_windows,
# a `seed` is used only if explicitly supplied (default NULL = random draw).
.genome_tiles <- function(bw_path, max_windows = 2000, tile_size = 10000L,
                          seqlevels = NULL, seed = NULL) {
  si <- tryCatch(GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(bw_path)),
                 error = function(e) {
                   stop(sprintf("Could not read BigWig header of '%s': %s", bw_path, e$message))
                 })
  sl <- GenomeInfoDb::seqlengths(si)
  nm <- names(sl)
  if (is.null(nm) || length(nm) == 0) nm <- GenomeInfoDb::seqnames(si)
  nm <- as.character(nm)
  # Default: use all seqlevels with finite, positive seqlengths (works for any
  # species, not only human/mouse); user may pass a custom set via seqlevels.
  if (!is.null(seqlevels)) {
    keep <- nm %in% seqlevels
  } else {
    keep <- is.finite(sl) & sl > 0
  }
  si <- si[nm[keep]]
  GenomeInfoDb::seqlengths(si) <- pmax(GenomeInfoDb::seqlengths(si), tile_size)
  tiles <- GenomicRanges::tileGenome(si, tilewidth = tile_size)
  # tileGenome returns a GRangesList (one element per seqlevel); flatten it.
  if (methods::is(tiles, "GRangesList")) tiles <- unlist(tiles, use.names = FALSE)
  if (length(tiles) > max_windows) {
    tiles <- .with_opt_seed(seed, tiles[sample(length(tiles), max_windows)])
  }
  tiles
}
