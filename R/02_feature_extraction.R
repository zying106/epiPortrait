#' Extract and Build Multi-dimensional Portrait Matrix
#'
#' @description Extracts domain-level quantitative features from BigWig files
#' across consensus or stitched peaks: \code{Intensity} (integrated signal
#' magnitude) and \code{SignalDispersion} (signal-weighted genomic dispersion).
#' When sample-specific native peak files (\code{peak_path}) are provided,
#' additional sample-specific breadth geometry is added as
#' \code{NativeMaxPeakWidth}, \code{NativeOccupiedWidth} and
#' \code{NativePeakCount}. This is the core engine of epiPortrait.
#'
#' @param sample_sheet A \code{data.frame} containing at least 'SampleID',
#'   'Condition', and 'bw_path'. Optional \code{peak_path}: path to the
#'   sample-specific native peak call file (BED/narrowPeak/broadPeak readable by
#'   \code{rtracklayer::import}); when absent, native breadth geometry is not
#'   computed and Breadth-Super calling is unavailable for that sample.
#'   Recommended additional metadata columns, all passed through to
#'   colData(se): 'Replicate', 'Mark', 'Genome', 'Normalization', 'SignalType',
#'   'BinSize', 'SpikeIn'. These enable provenance tracking and
#'   condition-aware replicate handling.
#' @param consensus_peaks A \code{GRanges} object of unified/stitched peaks
#'   (the shared candidate domain universe).
#' @param workers Integer. Number of threads for parallel processing (default: 1).
#' @param on_disk Logical. If TRUE, stores assays as HDF5 files on disk to save RAM,
#' which is highly recommended for mammalian-scale BigWig processing.
#' @param h5_dir Character. Directory to store HDF5 files if \code{on_disk = TRUE}.
#' @param cache_dir Character or NULL. If provided, caches per-sample RAW
#' coverage views as RDS files in this directory. Subsequent calls with the same
#' BigWig + consensus peaks combination skip BigWig I/O entirely (default: NULL).
#' The cache stores raw (un-sanitized) coverage so that negative_policy and
#' custom features can change between runs without cache pollution.
#' @param custom_features A named list of functions. Each function must accept a
#' single \code{NumericList} (the coverage view for one peak) and return a single
#' numeric value. Results are added as additional assays alongside the core
#' dimensions (e.g., \code{list(Entropy = function(x) ...)}).
#' @param negative_policy Character. How to handle negative signal values in the
#'   BigWig tracks:
#'   \itemize{
#'     \item \code{"error"} (default): fail if any negative signal is found,
#'           forcing the user to provide non-negative (e.g. clipped,
#'           CPM/RPGC/spike-in) tracks.
#'     \item \code{"clip_zero"}: set negative values to 0 and continue, recording
#'           the negative fraction in the returned object.
#'     \item \code{"allow"}: keep negatives in the Intensity assay (use only
#'           for advanced within-sample Intensity analyses with an explicitly
#'           compatible downstream transform). Note that \code{SignalDispersion}
#'           always uses non-negative weights (negative values are clipped to 0
#'           internally, since signed signal is not a valid probability-like
#'           weight), and the default \code{log10(pmax(x,0)+1)} transform /
#'           tangent inflection of Super calling also clip negatives to 0.
#'           Signed tracks are not part of the recommended canonical workflow.
#'   }
#' @param fail_action Character. Behavior when a BigWig cannot be imported
#'   \code{"stop"} (default) fails the whole run; \code{"drop"}
#'   excludes the failed sample(s) with a clear warning about the changed
#'   sample composition.
#' @return A \code{SummarizedExperiment} object. Core assays (domain x sample):
#'   \itemize{
#'     \item \code{Intensity}: integrated BigWig signal within the shared domain.
#'     \item \code{SignalDispersion}: signal-weighted genomic SD within the
#'           domain (secondary within-domain architecture descriptor).
#'     \item \code{NativeMaxPeakWidth}: max width of native peaks overlapping the
#'           domain (NA when \code{peak_path} missing).
#'     \item \code{NativeOccupiedWidth}: sum of reduced native peak widths
#'           inside the domain (NA when \code{peak_path} missing).
#'     \item \code{NativePeakCount}: number of native peaks overlapping the
#'           domain (NA when \code{peak_path} missing).
#'   }
#'   The static shared interval length is stored in \code{rowData(se)$IntervalWidth}.
#'   Imported native peaks (if provided) are stored in
#'   \code{metadata(se)$native_peaks} for downstream peak-level Breadth-Super
#'   calling.
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import SummarizedExperiment
#' @import BiocParallel
#' @importFrom stats quantile
#' @importFrom S4Vectors DataFrame
#' @examples
#' if (.Platform$OS.type != "windows") {
#'   extdata <- system.file("extdata", package = "epiPortrait")
#'   samples <- data.frame(
#'     SampleID  = c("C1", "C2", "T1", "T2"),
#'     Condition = c("Control", "Control", "Treatment", "Treatment"),
#'     bw_path   = file.path(extdata, c("C1.bw", "C2.bw", "T1.bw", "T2.bw")),
#'     peak_path = file.path(extdata, c("C1_peaks.bed", "C2_peaks.bed",
#'                                      "T1_peaks.bed", "T2_peaks.bed")))
#'   peaks <- rtracklayer::import(file.path(extdata, "peaks.bed"))
#'   build_portrait_matrix(samples, consensus_peaks = peaks, workers = 1)
#' }
#' @export
build_portrait_matrix <- function(sample_sheet, consensus_peaks,
                                   workers = 1,
                                   on_disk = FALSE,
                                   h5_dir = "epiPortrait_h5",
                                   cache_dir = NULL,
                                   custom_features = NULL,
                                   negative_policy = c("error", "clip_zero", "allow"),
                                   fail_action = c("stop", "drop")) {
  negative_policy <- match.arg(negative_policy)
  fail_action <- match.arg(fail_action)

  # 1. Input validation
  req_cols <- c("SampleID", "Condition", "bw_path")
  if (!all(req_cols %in% colnames(sample_sheet))) {
    stop("sample_sheet must contain columns: 'SampleID', 'Condition', and 'bw_path'")
  }
  # Freeze review 2026-08-11: many downstream steps index by sample name
  # (metadata(se)$native_peaks[[s]], colnames(se), replicate_call_matrix), so a
  # missing/empty/duplicated SampleID must fail loudly at the entry point
  # rather than producing ambiguous list indexing / output columns.
  if (anyNA(sample_sheet$SampleID) ||
      any(!nzchar(as.character(sample_sheet$SampleID))) ||
      anyDuplicated(sample_sheet$SampleID)) {
    stop("SampleID must be non-missing, non-empty, and unique.")
  }
  if (anyNA(sample_sheet$Condition) ||
      any(!nzchar(as.character(sample_sheet$Condition)))) {
    stop("Condition must be non-missing and non-empty.")
  }

  if (!all(file.exists(sample_sheet$bw_path))) {
    missing_files <- sample_sheet$bw_path[!file.exists(sample_sheet$bw_path)]
    stop(paste("The following BigWig files do not exist:", paste(missing_files, collapse = ", ")))
  }

  # Genome assembly is an experimental condition, not decorative metadata.
  # When supplied, require one explicit assembly across every sample and carry
  # it into the object contract for downstream cross-object integration.
  genome_id <- NULL
  if ("Genome" %in% colnames(sample_sheet)) {
    gv <- as.character(sample_sheet$Genome)
    if (anyNA(gv) || any(!nzchar(gv))) {
      stop("Genome must be non-missing and non-empty for every sample when the ",
           "column is supplied.")
    }
    gu <- unique(gv)
    if (length(gu) != 1L) {
      stop("All samples must use the same genome assembly. Found: ",
           paste(gu, collapse = ", "), ".")
    }
    genome_id <- gu
  }

  if (!inherits(consensus_peaks, "GRanges")) {
    stop("consensus_peaks must be a GRanges object.")
  }
  if (length(consensus_peaks) == 0L) {
    stop("consensus_peaks is empty; provide at least one candidate domain.")
  }

  # Header-level chromosome compatibility. A region absent from a BigWig is
  # represented by zero coverage by the sparse import path, so without this
  # check a chr/build mismatch can silently become biological low signal.
  domain_seqlevels <- unique(as.character(GenomicRanges::seqnames(consensus_peaks)))
  for (i in seq_len(nrow(sample_sheet))) {
    bw_sl <- tryCatch({
      si <- GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(sample_sheet$bw_path[i]))
      as.character(GenomeInfoDb::seqlevels(si))
    }, error = function(e) character(0))
    # If the upstream library cannot read the header, defer to the ordinary
    # import error path below (also keeps mocked unit tests independent of I/O).
    if (length(bw_sl) == 0L) next
    shared <- intersect(domain_seqlevels, bw_sl)
    if (length(shared) == 0L) {
      stop("No shared seqlevels between consensus_peaks and BigWig '",
           sample_sheet$bw_path[i], "'. Check genome build and chromosome ",
           "naming (for example chr1 vs 1).")
    }
    missing_sl <- setdiff(domain_seqlevels, bw_sl)
    if (length(missing_sl) > 0L) {
      warning(sprintf(
        "BigWig '%s' lacks %d/%d candidate-domain seqlevels (%s); domains on those seqlevels would otherwise be read as zero signal.",
        basename(sample_sheet$bw_path[i]), length(missing_sl),
        length(domain_seqlevels), paste(utils::head(missing_sl, 5),
                                       collapse = ", ")), call. = FALSE)
    }
  }

  # Native peak files are OPTIONAL. When absent for a sample, native breadth
  # geometry is NA and that sample cannot contribute to Breadth-Super calling.
  # has_peaks is ALWAYS a per-sample logical vector (length = nrow(sample_sheet))
  # so indexing has_peaks[i] is safe in both modes (P0-3).
  if ("peak_path" %in% colnames(sample_sheet)) {
    has_peaks <- !is.na(sample_sheet$peak_path) &
      nzchar(sample_sheet$peak_path)
  } else {
    has_peaks <- rep(FALSE, nrow(sample_sheet))
  }
  if (any(has_peaks)) {
    missing_peaks <- sample_sheet$peak_path[has_peaks][
      !file.exists(sample_sheet$peak_path[has_peaks])]
    if (length(missing_peaks) > 0) {
      stop(paste("The following native peak files do not exist:",
                 paste(missing_peaks, collapse = ", ")))
    }
    message(sprintf("Native peak files provided for %d/%d samples; Breadth-Super calling will be available.",
                    sum(has_peaks), nrow(sample_sheet)))
  }

  # 2. Configure parallel backend
  if (workers > 1) {
    if (.Platform$OS.type == "unix") {
      param <- BiocParallel::MulticoreParam(workers = workers)
    } else {
      param <- BiocParallel::SnowParam(workers = workers)
    }
  } else {
    param <- BiocParallel::SerialParam()
  }

  # Validate custom features
  use_custom <- !is.null(custom_features) && length(custom_features) > 0
  if (use_custom) {
    if (!is.list(custom_features) || is.null(names(custom_features))) {
      stop("custom_features must be a named list of functions.")
    }
    for (nm in names(custom_features)) {
      if (!is.function(custom_features[[nm]])) {
        stop(sprintf("custom_features[['%s']] is not a function.", nm))
      }
    }
    message(sprintf("Including %d custom feature(s): %s",
                    length(custom_features),
                    paste(names(custom_features), collapse = ", ")))
  }

  message(sprintf("Drawing domain portraits for %d peaks across %d samples...",
                  length(consensus_peaks), nrow(sample_sheet)))

  # 3. Pre-compute loop-invariant constants (Width is constant across samples)
  peak_widths <- GenomicRanges::width(consensus_peaks)

  # 3.5 Prepare caching infrastructure
  use_cache <- !is.null(cache_dir) && nzchar(cache_dir)
  if (use_cache) {
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    # Compute a stable fingerprint of the consensus peak set for cache keys
    region_hash <- digest::digest(list(
      ranges   = as.character(ranges(consensus_peaks)),
      seqnames = as.character(seqnames(consensus_peaks))
    ), algo = "xxhash64")
    message(sprintf("Cache enabled: %s (region hash: %s)", cache_dir, substr(region_hash, 1, 8)))
  }

  # Cache key per sample: full path + size + mtime + region hash.
  # Using basename() alone causes collisions across directories and stale reads
  # when a BigWig file is updated in place (P1-3).
  .cache_key <- function(bw_path) {
    finfo <- file.info(bw_path)
    digest::digest(list(
      path  = normalizePath(bw_path, mustWork = FALSE),
      size  = finfo$size,
      mtime = as.numeric(finfo$mtime),
      region_hash = region_hash
    ), algo = "xxhash64")
  }

  # 4. Parallel feature extraction core
  num_peaks <- length(consensus_peaks)
  # P1-4: process regions in batches so the base-level NumericList for all
  # domains never has to reside in memory at once (the real memory bottleneck).
  batch_size <- 5000L

  feature_list <- BiocParallel::bplapply(seq_len(nrow(sample_sheet)), function(i) {

    bw_path <- sample_sheet$bw_path[i]
    peak_path <- if (has_peaks[i]) sample_sheet$peak_path[i] else NA_character_

    tryCatch({

      # ---- load native peaks (optional) for breadth geometry ----------------
      # Native geometry is computed from the sample's own peak calls, NOT from
      # the shared consensus universe (Breadth-Super design revision).
      native_peaks <- NULL
      if (!is.na(peak_path)) {
        native_peaks <- rtracklayer::import(peak_path)
        if (length(native_peaks) == 0) native_peaks <- NULL
        if (!is.null(native_peaks)) {
          native_sl <- unique(as.character(GenomicRanges::seqnames(native_peaks)))
          if (length(intersect(domain_seqlevels, native_sl)) == 0L) {
            stop("No shared seqlevels between consensus_peaks and native peak ",
                 "file '", peak_path, "'. Check genome build and chromosome ",
                 "naming (for example chr1 vs 1).")
          }
        }
      }
      native_geom <- if (!is.null(native_peaks)) {
        .native_peak_geometry(native_peaks, consensus_peaks)
      } else {
        list(NWmax = rep(NA_real_, num_peaks),
             NWocc = rep(NA_real_, num_peaks),
             NWcnt = rep(NA_real_, num_peaks))
      }

      # ---- cache / import coverage in batches (P1-4) ---------------------
      # P0-A: the cache must store RAW imported coverage (before any
      # negative_policy / baseline handling). Sanitization happens AFTER
      # loading, using the CURRENT run's parameters — otherwise a previous
      # clip_zero run would hide negatives from a later "error" run, and a
      # later "allow" run would consume clipped (contaminated) signal.
      if (use_cache) {
        cache_file <- file.path(cache_dir, paste0(.cache_key(bw_path), ".rds"))
        if (file.exists(cache_file)) {
          bw_views <- tryCatch(readRDS(cache_file), error = function(e) NULL)
          if (!is.null(bw_views)) {
            message(sprintf("[cache hit] %s", basename(bw_path)))
          } else {
            bw_views <- NULL
          }
        } else {
          bw_views <- NULL
        }
      } else {
        bw_views <- NULL
      }

      if (!is.null(bw_views)) {
        # ---- extract from RAW cached NumericList ----------------------------
        # Sanitize with the CURRENT negative_policy (P0-A, P0-1).
        n_neg <- sum(vapply(bw_views, function(x) sum(x < 0, na.rm = TRUE), integer(1)))
        n_positions <- sum(vapply(bw_views, function(x) sum(!is.na(x)), integer(1)))
        bw_views <- lapply(bw_views, function(x) {
          if (is.null(x)) return(x)
          .sanitize_signal(as.numeric(x), negative_policy = negative_policy)
        })
        bw_views <- IRanges::NumericList(bw_views)

        intensities <- sum(bw_views, na.rm = TRUE)
        dispersions <- vapply(bw_views, .signal_dispersion, numeric(1))
      } else {
        # ---- batch import to bound memory (P1-4) ---------------------------
        # P2-fix: accumulate into DOUBLE vectors. batch_ints was integer(); R
        # silently promotes on assignment so values were never truncated, but
        # the type was misleading for continuous (e.g. CPM/RPGC) signal.
        batch_ints <- numeric(num_peaks)
        batch_disp <- rep(NA_real_, num_peaks)
        batch_neg <- integer(num_peaks)
        batch_pos <- integer(num_peaks)
        raw_views <- if (use_cache) vector("list", num_peaks) else NULL
        n_batches <- ceiling(num_peaks / batch_size)
        for (b in seq_len(n_batches)) {
          idx <- ((b - 1) * batch_size + 1):min(b * batch_size, num_peaks)
          cvg_b <- .import_bw_views(bw_path, consensus_peaks[idx])
          # Count negatives BEFORE sanitization so the fraction reflects the
          # raw input (P0-1); error policy throws inside .sanitize_signal.
          batch_neg[idx] <- vapply(cvg_b, function(x) sum(x < 0, na.rm = TRUE), integer(1))
          batch_pos[idx] <- vapply(cvg_b, function(x) sum(!is.na(x)), integer(1))
          # Store RAW views for the cache (P0-A): sanitization below must not
          # be persisted.
          if (use_cache) {
            for (j in seq_along(idx)) raw_views[[idx[j]]] <- cvg_b[[j]]
          }
          # Sanitize BEFORE any feature computation (P0-1): the same vector is
          # used for Intensity, SignalDispersion and custom features, so
          # clip_zero semantics are consistent across all of them.
          cvg_b <- lapply(cvg_b, function(x) {
            if (is.null(x)) return(x)
            .sanitize_signal(as.numeric(x), negative_policy = negative_policy)
          })
          cvg_b <- IRanges::NumericList(cvg_b)

          batch_ints[idx] <- sum(cvg_b, na.rm = TRUE)
          batch_disp[idx] <- vapply(cvg_b, .signal_dispersion, numeric(1))
        }
        intensities <- batch_ints
        dispersions <- batch_disp
        n_neg <- sum(batch_neg)
        n_positions <- sum(batch_pos)

        # Persist RAW per-sample coverage for future runs (P0-A, P1-3 key).
        # NEVER persist sanitized views.
        if (use_cache) {
          cache_file <- file.path(cache_dir, paste0(.cache_key(bw_path), ".rds"))
          tryCatch(saveRDS(IRanges::NumericList(raw_views), cache_file),
                   error = function(e) warning("Failed to write cache: ", e$message))
        }
      }

      # Correct negative-fraction denominator: negative positions over total
      # signal positions (P0-1). After clip_zero no negatives remain.
      neg_frac <- if (n_positions > 0) n_neg / n_positions else 0
      if (n_neg > 0 && negative_policy == "error") {
        stop("BigWig '", bw_path, "' contains ", n_neg,
             " negative value(s) (", format(100 * neg_frac, digits = 3),
             "% of positions). ",
             "epiPortrait expects normalized non-negative signal tracks. ",
             "Clip or provide CPM/RPGC/spike-in normalized tracks, or set ",
             "negative_policy = 'clip_zero'.")
      } else if (n_neg > 0 && negative_policy == "clip_zero") {
        # clip_zero already applied inside .sanitize_signal(); here the
        # negatives have been set to 0, so report the pre-clip fraction.
        message(sprintf("Clipped %d negative value(s) (%.3f%% of positions) to 0 in '%s'.",
                        n_neg, 100 * neg_frac, basename(bw_path)))
      }

      # Custom feature extraction: needs the per-region views. In batch mode
      # we re-import the views per region one at a time (rare feature; keeps
      # the common path memory-bounded). Same negative policy applied.
      custom_result <- list()
      if (use_custom) {
        for (fn_name in names(custom_features)) {
          fn <- custom_features[[fn_name]]
          custom_result[[fn_name]] <- vapply(seq_len(num_peaks), function(k) {
            cvg_k <- .import_bw_views(bw_path, consensus_peaks[k])[[1]]
            x_clean <- .sanitize_signal(as.numeric(cvg_k), negative_policy = negative_policy)
            x_clean <- x_clean[!is.na(x_clean)]
            if (length(x_clean) >= 3) fn(x_clean) else NA_real_
          }, numeric(1))
        }
      }

      c(list(I = intensities, D = dispersions,
             NWmax = native_geom$NWmax,
             NWocc = native_geom$NWocc,
             NWcnt = native_geom$NWcnt,
             NegFraction = neg_frac),
        custom_result,
        list(NativePeaks = native_peaks,
             Error = NULL))

    }, error = function(e) {
      err_list <- list(I = rep(NA_real_, num_peaks),
                       D = rep(NA_real_, num_peaks),
                       NWmax = rep(NA_real_, num_peaks),
                       NWocc = rep(NA_real_, num_peaks),
                       NWcnt = rep(NA_real_, num_peaks),
                       NegFraction = NA_real_)
      if (use_custom) {
        for (fn_name in names(custom_features)) {
          err_list[[fn_name]] <- rep(NA_real_, num_peaks)
        }
      }
      c(err_list, list(NativePeaks = NULL,
                       Error = paste("Sample", sample_sheet$SampleID[i], "-", e$message)))
    })

  }, BPPARAM = param)

  # Check for per-sample failures in parallel computation.
  # P1-10: default is to STOP rather than silently drop samples, because a
  # failed replicate silently changes the sample composition (e.g. Treatment
  # n=2 -> n=1) and invalidates downstream group analyses.
  error_flags <- !vapply(feature_list, function(x) is.null(x$Error), logical(1))
  if (any(error_flags)) {
    error_msgs <- vapply(feature_list[error_flags], function(x) x$Error, character(1))
    if (fail_action == "stop") {
      stop(sprintf(
        "%d sample(s) failed during BigWig import:\n  %s\n",
        sum(error_flags), paste(error_msgs, collapse = "\n  "),
        "Fix the failing input(s) or set fail_action = 'drop' to exclude them explicitly."
      ))
    }
    warning(sprintf("%d sample(s) failed during BigWig import and will be excluded:\n  %s",
                    sum(error_flags), paste(error_msgs, collapse = "\n  ")))
    feature_list <- feature_list[!error_flags]
    sample_sheet <- sample_sheet[!error_flags, , drop = FALSE]
    has_peaks <- has_peaks[!error_flags]   # P1-8: keep availability aligned
    if (nrow(sample_sheet) == 0) {
      stop("All samples failed. Please check your BigWig files.")
    }
    # Emit the changed sample composition clearly
    comp <- table(sample_sheet$Condition)
    warning(sprintf("Sample composition after drop: %s",
                    paste(sprintf("%s n=%d", names(comp), comp), collapse = ", ")))
  }

  # 5. Assemble matrices and route storage (RAM vs HDF5)
  num_peaks <- length(consensus_peaks)
  num_samples <- nrow(sample_sheet)
  sample_names <- sample_sheet$SampleID

  # Assemble matrices with do.call(cbind, ...) instead of explicit for-loops.
  # NOTE (P1-5): the static genomic interval length (IntervalWidth) is moved to
  # rowData — it is identical across samples and is not a per-sample signal
  # feature. NOTE (v1.0 design): canonical assays are Intensity,
  # SignalDispersion and (when peak files are provided) NativeMaxPeakWidth /
  # NativeOccupiedWidth / NativePeakCount. RobustHeight and EffectiveWidth have
  # been removed from the core design.
  mat_list <- list(
    Intensity = do.call(cbind, lapply(feature_list, function(x) x$I)),
    SignalDispersion = do.call(cbind, lapply(feature_list, function(x) x$D)),
    NativeMaxPeakWidth = do.call(cbind, lapply(feature_list, function(x) x$NWmax)),
    NativeOccupiedWidth = do.call(cbind, lapply(feature_list, function(x) x$NWocc)),
    NativePeakCount = do.call(cbind, lapply(feature_list, function(x) x$NWcnt))
  )

  # Append custom feature matrices
  if (use_custom) {
    for (fn_name in names(custom_features)) {
      mat_list[[fn_name]] <- do.call(cbind, lapply(feature_list, function(x) {
        if (fn_name %in% names(x)) x[[fn_name]] else rep(NA_real_, num_peaks)
      }))
    }
  }

  # Add column names to all matrices
  mat_list <- lapply(mat_list, function(m) {
    colnames(m) <- sample_names
    return(m)
  })

  # Route to HDF5 on-disk storage if requested
  if (on_disk) {
    if (!requireNamespace("HDF5Array", quietly = TRUE)) {
      stop("Please install 'HDF5Array' to use the on_disk feature.")
    }
    if (!dir.exists(h5_dir)) dir.create(h5_dir, recursive = TRUE)

    for (feat_name in names(mat_list)) {
      h5_path <- file.path(h5_dir, paste0(feat_name, ".h5"))
      mat_list[[feat_name]] <- HDF5Array::writeHDF5Array(
        mat_list[[feat_name]],
        filepath = h5_path,
        name = feat_name,
        verbose = FALSE
      )
    }
  }

  # 6. Build the final SummarizedExperiment object
  # Ensure unique Domain_IDs (P1-9): duplicate or NULL names would silently
  # mis-map during match() in downstream ranking.
  se_ranges <- consensus_peaks
  if (is.null(names(se_ranges)) || any(names(se_ranges) == "") ||
      anyDuplicated(names(se_ranges))) {
    message("Assigning unique domain IDs (epiDomain_xxxxxx).")
    names(se_ranges) <- sprintf("epiDomain_%06d", seq_len(num_peaks))
  }
  stopifnot(anyDuplicated(names(se_ranges)) == 0)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = mat_list,
    rowRanges = se_ranges,
    colData = S4Vectors::DataFrame(sample_sheet, row.names = sample_names)
  )

  # Static genomic interval length lives in rowData (P1-5), not as an assay.
  rowData(se)$IntervalWidth <- GenomicRanges::width(se_ranges)

  # Record the negative-signal fraction per sample as a QC provenance column
  # (P0-1): the fraction of raw signal positions that were negative before the
  # negative policy was applied.
  se$NegativeFraction <- vapply(feature_list, function(x) x$NegFraction, numeric(1))

  # ---- Metadata contract (object-contract design) ---------------------------
  S4Vectors::metadata(se)$signal_contract <- list(
    input_type = "BigWig",
    required_signal = "normalized non-negative continuous signal",
    negative_policy = negative_policy,
    native_peaks = if (any(has_peaks)) "provided (per-sample peak files)" else "not provided",
    fail_action = fail_action,
    genome = genome_id
  )
  # Domain-source / caller provenance (heterochromatin plan 2026-08-10):
  # broad repressive marks (H3K27me3/H3K9me3) rely on upstream broad-domain
  # callers (SICER / epic2 / RECOGNICER / MACS2_broad / precalled), whose
  # window/gap/boundary strategy determines domain width. Optional sample-sheet
  # columns (domain_source, domain_caller, caller_version, caller_mode) are
  # passed through to colData and mirrored here for auditability.
  S4Vectors::metadata(se)$domain_provenance <- list(
    domain_source = if ("domain_source" %in% colnames(sample_sheet)) {
      unique(as.character(sample_sheet$domain_source))
    } else "unspecified",
    domain_caller = if ("domain_caller" %in% colnames(sample_sheet)) {
      unique(as.character(sample_sheet$domain_caller))
    } else "unspecified",
    caller_version = if ("caller_version" %in% colnames(sample_sheet)) {
      unique(as.character(sample_sheet$caller_version))
    } else "unspecified",
    caller_mode = if ("caller_mode" %in% colnames(sample_sheet)) {
      unique(as.character(sample_sheet$caller_mode))
    } else "unspecified",
    note = paste(
      "peak_path may contain narrow peaks or broad enriched domains",
      "depending on the histone mark and upstream caller")
  )
  S4Vectors::metadata(se)$feature_definitions <- list(
    Intensity = list(
      definition = "Integrated BigWig signal within candidate domain",
      units = "signal x bp"
    ),
    SignalDispersion = list(
      definition = "Signal-weighted genomic SD within candidate domain (within-domain architecture descriptor)",
      units = "bp"
    ),
    NativeMaxPeakWidth = list(
      definition = "Max width of native peaks overlapping the domain (any-overlap geometry, not Breadth unique mapping)",
      units = "bp",
      available = any(has_peaks)
    ),
    NativeOccupiedWidth = list(
      definition = "Sum of reduced native peak widths inside the domain",
      units = "bp",
      available = any(has_peaks)
    ),
    NativePeakCount = list(
      definition = "Number of native peaks overlapping the domain",
      units = "count",
      available = any(has_peaks)
    ),
    IntervalWidth = list(
      definition = "Candidate genomic interval width",
      units = "bp",
      location = "rowData"
    )
  )
  # Imported native peaks, kept for peak-level Breadth-Super calling. Stored as
  # a named list (one GRanges per sample with peak files; NULL otherwise).
  if (any(has_peaks)) {
    np <- lapply(feature_list, function(x) x$NativePeaks)
    names(np) <- sample_names
    S4Vectors::metadata(se)$native_peaks <- np
    S4Vectors::metadata(se)$native_peak_available <- has_peaks[seq_len(nrow(sample_sheet))]
  }
  S4Vectors::metadata(se)$provenance <- list(
    epiPortrait_version = if (requireNamespace("epiPortrait", quietly = TRUE))
      as.character(utils::packageVersion("epiPortrait")) else "source",
    R_version = R.version.string,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  message(sprintf("Success: Portrait Matrix built with %d assays (On-disk = %s).",
                  length(mat_list), on_disk))
  return(se)
}


# Import per-region BigWig coverage as a NumericList of full-length per-base
# vectors, using the GRanges (sparse) import path rather than `as="NumericList"`.
#
# The GRanges path avoids the Windows-specific `NumericList` failure observed in
# some `rtracklayer` builds, but local BigWig I/O may still fail on Windows
# through the upstream UCSC backend ("UCSC library operation failed";
# rtracklayer#52/#62/#128/#151); that case is re-raised with a clear, actionable
# message (see the tryCatch below). The expansion reproduces the exact numeric
# semantics of `as="NumericList"` on platforms where BigWig reads succeed.
.import_bw_views <- function(bw_path, regions) {
  sparse <- tryCatch(
    rtracklayer::import.bw(bw_path, which = regions),
    error = function(e) {
      msg <- conditionMessage(e)
      # rtracklayer raises 'UCSC library operation failed' (and sometimes
      # 'attempt to apply non-function') when reading local BigWig on the
      # Windows binary (rtracklayer#52/#62/#128/#151). Give users a clear,
      # actionable message instead of an opaque backend error, so a
      # Windows user can ingest BigWigs on Linux/macOS or feed peak files /
      # a pre-built matrix.
      if (.Platform$OS.type == "windows" &&
          grepl("UCSC|BigWig|operation failed|non-function", msg, ignore.case = TRUE)) {
        stop(
          "BigWig region import failed on Windows through the upstream ",
          "rtracklayer/UCSC backend ('", msg, "'). This is a platform-specific ",
          "upstream limitation (rtracklayer#52/#62/#128/#151), not an ",
          "epiPortrait analysis error. Options: ingest BigWigs on ",
          "Linux/macOS, or provide per-sample native peak (BED/narrowPeak) ",
          "files / a pre-built portrait matrix so BigWig I/O is not needed.",
          call. = FALSE)
      }
      # non-Windows or unrelated error: re-raise as-is
      stop(conditionMessage(e), call. = FALSE)
    })
  region_w <- GenomicRanges::width(regions)
  region_start <- GenomicRanges::start(regions)
  views <- lapply(seq_len(length(regions)), function(i) {
    rep(0, region_w[i])
  })
  if (length(sparse) > 0) {
    ov <- GenomicRanges::findOverlaps(sparse, regions)
    qh <- S4Vectors::queryHits(ov)
    sh <- S4Vectors::subjectHits(ov)
    # a sparse bin spans width(sparse)[qh] bp and carries the bin score at
    # every covered base (matching the `as="NumericList"` expansion); clip to
    # the region boundary so a bin straddling two regions is split correctly.
    bin_start <- GenomicRanges::start(sparse)[qh]
    bin_end   <- GenomicRanges::end(sparse)[qh]
    score     <- S4Vectors::mcols(sparse)$score[qh]
    for (j in seq_along(sh)) {
      lo <- max(bin_start[j], region_start[sh[j]]) - region_start[sh[j]] + 1L
      hi <- min(bin_end[j], region_start[sh[j]] + region_w[sh[j]] - 1L) -
        region_start[sh[j]] + 1L
      if (lo <= hi) views[[sh[j]]][lo:hi] <- score[j]
    }
  }
  IRanges::NumericList(views)
}


# Sanitize a per-region signal vector according to the negative policy (P0-1).
#
# This MUST be applied once, before any feature (Intensity, SignalDispersion,
# custom) is computed, so that clip_zero semantics are consistent across all
# features. NA handling: missing positions are treated as 0 (preserving the
# genomic coordinate axis for span metrics).
.sanitize_signal <- function(x, negative_policy = c("error", "clip_zero", "allow")) {
  negative_policy <- match.arg(negative_policy)
  # keep positions: NA -> 0 so genomic span metrics are not compressed (P1-25)
  x[is.na(x)] <- 0
  if (negative_policy == "clip_zero") {
    x[x < 0] <- 0
  } else if (negative_policy == "error") {
    n_neg <- sum(x < 0, na.rm = TRUE)
    if (n_neg > 0) {
      stop("Signal contains ", n_neg, " negative value(s) (",
           format(100 * n_neg / max(sum(!is.na(x)), 1), digits = 3),
           "% of positions). ",
           "epiPortrait expects normalized non-negative signal tracks. ",
           "Clip or provide CPM/RPGC/spike-in normalized tracks, or set ",
           "negative_policy = 'clip_zero'.")
    }
  }
  x
}


# Signal-weighted genomic dispersion (SignalDispersion).
#
# For a domain with per-position signal x_i at positions p_i (relative bp
# offset within the domain), the dispersion is the signal-weighted SD:
#
#     w_i = x_i / sum(x_i)
#     mu  = sum(w_i * p_i)
#     SignalDispersion = sqrt(sum(w_i * (p_i - mu)^2))
#
# It answers "how far does the signal mass spread in genomic space", WITHOUT a
# peak-caller boundary. Two domains with the same occupied width but different
# internal architecture (compact vs multi-modal) receive different values. This
# is a SECONDARY within-domain architecture descriptor, not the Breadth-Super
# caller (Breadth-Super uses native PeakWidth; design doc 2026-08-09).
.signal_dispersion <- function(x) {
  xc <- as.numeric(x)
  xc[is.na(xc)] <- 0
  xc[xc < 0] <- 0
  total <- sum(xc)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  pos <- seq_along(xc)
  w <- xc / total
  mu <- sum(w * pos)
  sqrt(max(sum(w * (pos - mu)^2), 0))
}


# Native peak geometry within each shared domain (v1.0 breadth design).
#
# Computes per-domain sample-specific breadth quantities from the sample's own
# native peak calls, restricted to the shared domain intervals:
#
#   NativeMaxPeakWidth : max(width(overlapping native peaks))
#   NativeOccupiedWidth: sum(width(reduce(native peaks clipped to the domain)))
#   NativePeakCount    : number of native peaks overlapping the domain
#
# Any-overlap is used here because these are GEOMETRY DESCRIPTORS (for display
# and interpretation), not Breadth-Super evidence. Breadth-Super evidence uses
# the unique-assignment mapping in call_super_domains(feature = "Breadth").
.native_peak_geometry <- function(native_peaks, domains) {
  hits <- GenomicRanges::findOverlaps(domains, native_peaks)
  if (length(hits) == 0) {
    return(list(NWmax = rep(NA_real_, length(domains)),
                NWocc = rep(NA_real_, length(domains)),
                NWcnt = rep(0L, length(domains))))
  }
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  widths <- GenomicRanges::width(native_peaks)[sh]
  n_peaks <- tabulate(qh, nbins = length(domains))

  # NativeMaxPeakWidth: max width of peaks overlapping each domain. Fully
  # vectorized via split + vapply (no per-domain R loop).
  max_w <- rep(NA_real_, length(domains))
  mx <- vapply(split(widths, qh), max, numeric(1))
  max_w[as.integer(names(mx))] <- mx

  # NativeOccupiedWidth: union of native peaks clipped to the domain. Clip all
  # peaks to their overlapping domain at once (pintersect), then reduce ONLY
  # the multi-peak domains. Single-peak domains need no reduce (occupied width
  # = clipped peak width). This avoids the O(n_domains) per-domain reduce loop
  # that made mammalian-scale build_portrait_matrix() unusable.
  ov <- GenomicRanges::pintersect(native_peaks[sh], domains[qh],
                                  ignore.strand = TRUE)
  occ_w <- rep(NA_real_, length(domains))
  single_peaks <- which(n_peaks == 1L)
  if (length(single_peaks) > 0) {
    idx <- match(single_peaks, qh)
    occ_w[single_peaks] <- GenomicRanges::width(ov)[idx]
  }
  multi <- which(n_peaks > 1L)
  if (length(multi) > 0) {
    sub_ov <- ov[qh %in% multi]
    sub_qh <- qh[qh %in% multi]
    occ_list <- split(sub_ov, factor(sub_qh, levels = multi))
    occ_w[multi] <- vapply(occ_list, function(g) {
      sum(GenomicRanges::width(GenomicRanges::reduce(g)))
    }, numeric(1))
  }

  list(NWmax = max_w, NWocc = occ_w, NWcnt = n_peaks)
}
