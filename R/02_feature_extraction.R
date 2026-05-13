#' Extract and Build Multi-dimensional Portrait Matrix
#'
#' @description Extracts Width, Intensity, Height, and robust Skewness from BigWig files
#' across consensus or stitched peaks. This function is the core engine of epiPortrait,
#' generating a 4D feature representation (Portrait) for epigenetic domains.
#'
#' @param sample_sheet A \code{data.frame} containing 'SampleID', 'Condition', and 'bw_path'.
#' @param consensus_peaks A \code{GRanges} object of unified/stitched peaks.
#' @param workers Integer. Number of threads for parallel processing (default: 1).
#' @param bg_quantile Numeric (0-1). Signals below this quantile within each peak are
#' excluded from skewness calculation to filter out background noise (Valley Trap).
#' @param on_disk Logical. If TRUE, stores assays as HDF5 files on disk to save RAM,
#' which is highly recommended for mammalian-scale BigWig processing.
#' @param h5_dir Character. Directory to store HDF5 files if \code{on_disk = TRUE}.
#' @param cache_dir Character or NULL. If provided, caches per-sample coverage views
#' as RDS files in this directory. Subsequent calls with the same BigWig + consensus
#' peaks + \code{bg_quantile} combination skip BigWig I/O entirely (default: NULL).
#' @param custom_features A named list of functions. Each function must accept a
#' single \code{NumericList} (the coverage view for one peak) and return a single
#' numeric value. Results are added as additional assays alongside the core 4D
#' dimensions (e.g., \code{list(Entropy = function(x) ...)}).
#'
#' @return A \code{SummarizedExperiment} object containing 4 assays (The Portrait matrix).
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import SummarizedExperiment
#' @import BiocParallel
#' @importFrom stats quantile
#' @importFrom moments skewness
#' @importFrom S4Vectors DataFrame
#' @export
build_portrait_matrix <- function(sample_sheet, consensus_peaks,
                                  workers = 1,
                                  bg_quantile = 0.1,
                                  on_disk = FALSE,
                                  h5_dir = "epiPortrait_h5",
                                  cache_dir = NULL,
                                  custom_features = NULL) {

  # 1. Input validation
  req_cols <- c("SampleID", "Condition", "bw_path")
  if (!all(req_cols %in% colnames(sample_sheet))) {
    stop("sample_sheet must contain columns: 'SampleID', 'Condition', and 'bw_path'")
  }

  if (!all(file.exists(sample_sheet$bw_path))) {
    missing_files <- sample_sheet$bw_path[!file.exists(sample_sheet$bw_path)]
    stop(paste("The following BigWig files do not exist:", paste(missing_files, collapse = ", ")))
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

  message(sprintf("Drawing 4D epigenetic portraits for %d peaks across %d samples...",
                  length(consensus_peaks), nrow(sample_sheet)))

  # 3. Pre-compute loop-invariant constants (Width is constant across samples)
  peak_widths <- GenomicRanges::width(consensus_peaks)

  # 3.5 Prepare caching infrastructure
  use_cache <- !is.null(cache_dir) && nzchar(cache_dir)
  if (use_cache) {
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    # Compute a stable fingerprint of the consensus peak set for cache keys
    hash_input <- serialize(list(
      ranges   = as.character(ranges(consensus_peaks)),
      seqnames = as.character(seqnames(consensus_peaks)),
      bg_q     = bg_quantile
    ), connection = NULL)
    region_hash <- as.character(sum(as.integer(hash_input)))
    message(sprintf("Cache enabled: %s (region hash: %s)", cache_dir, substr(region_hash, 1, 8)))
  }

  # 4. Parallel feature extraction core
  num_peaks <- length(consensus_peaks)

  feature_list <- BiocParallel::bplapply(seq_len(nrow(sample_sheet)), function(i) {

    bw_path <- sample_sheet$bw_path[i]

    tryCatch({

      bw_views <- NULL

      # Cache lookup: load cached coverage views if available
      if (use_cache) {
        cache_file <- file.path(cache_dir,
          paste0(basename(bw_path), "_", region_hash, ".rds"))
        if (file.exists(cache_file)) {
          bw_views <- tryCatch(readRDS(cache_file), error = function(e) NULL)
          if (!is.null(bw_views)) {
            message(sprintf("[cache hit] %s", basename(bw_path)))
          }
        }
      }

      # Cache miss: import from BigWig
      if (is.null(bw_views)) {
        bw_views <- rtracklayer::import(bw_path,
                                        selection = rtracklayer::BigWigSelection(consensus_peaks),
                                        as = "NumericList")
        if (use_cache) {
          cache_file <- file.path(cache_dir,
            paste0(basename(bw_path), "_", region_hash, ".rds"))
          tryCatch(saveRDS(bw_views, cache_file), error = function(e) {
            warning("Failed to write cache: ", e$message)
          })
        }
      }

      # Feature extraction
      intensities <- sum(bw_views, na.rm = TRUE)
      heights <- max(bw_views, na.rm = TRUE)
      heights[is.infinite(heights)] <- 0

      # Robust skewness with background-noise gate (Valley Trap)
      skewnesses <- vapply(bw_views, function(x) {
        x_clean <- x[!is.na(x)]
        if (length(x_clean) >= 3) {
          # names = FALSE speeds up quantile by skipping name generation
          bg_thresh <- stats::quantile(x_clean, probs = bg_quantile, names = FALSE)
          robust_signal <- x_clean[x_clean > bg_thresh]

          if (length(robust_signal) >= 3) {
            return(moments::skewness(robust_signal, na.rm = TRUE))
          }
        }
        return(0)
      }, numeric(1))

      # Custom feature extraction
      custom_result <- list()
      if (use_custom) {
        for (fn_name in names(custom_features)) {
          fn <- custom_features[[fn_name]]
          custom_result[[fn_name]] <- vapply(bw_views, function(x) {
            x_clean <- x[!is.na(x)]
            if (length(x_clean) >= 3) fn(x_clean) else NA_real_
          }, numeric(1))
        }
      }

      c(list(I = intensities, H = heights, S = skewnesses),
        custom_result,
        list(Error = NULL))

    }, error = function(e) {
      err_list <- list(I = rep(NA_real_, num_peaks),
                       H = rep(NA_real_, num_peaks),
                       S = rep(NA_real_, num_peaks))
      if (use_custom) {
        for (fn_name in names(custom_features)) {
          err_list[[fn_name]] <- rep(NA_real_, num_peaks)
        }
      }
      c(err_list, list(Error = paste("Sample", sample_sheet$SampleID[i], "-", e$message)))
    })

  }, BPPARAM = param)

  # Check for per-sample failures in parallel computation
  error_flags <- !vapply(feature_list, function(x) is.null(x$Error), logical(1))
  if (any(error_flags)) {
    error_msgs <- vapply(feature_list[error_flags], function(x) x$Error, character(1))
    warning(sprintf("%d sample(s) failed during BigWig import and will be excluded:\n  %s",
                    sum(error_flags), paste(error_msgs, collapse = "\n  ")))
    feature_list <- feature_list[!error_flags]
    sample_sheet <- sample_sheet[!error_flags, , drop = FALSE]
    if (nrow(sample_sheet) == 0) {
      stop("All samples failed. Please check your BigWig files.")
    }
  }

  # 5. Assemble matrices and route storage (RAM vs HDF5)
  num_peaks <- length(consensus_peaks)
  num_samples <- nrow(sample_sheet)
  sample_names <- sample_sheet$SampleID

  # Assemble matrices with do.call(cbind, ...) instead of explicit for-loops
  mat_list <- list(
    Intensity = do.call(cbind, lapply(feature_list, function(x) x$I)),
    Height    = do.call(cbind, lapply(feature_list, function(x) x$H)),
    Skewness  = do.call(cbind, lapply(feature_list, function(x) x$S)),
    Width     = matrix(rep(peak_widths, num_samples), ncol = num_samples)
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
  se_ranges <- consensus_peaks
  if (is.null(names(se_ranges))) {
    names(se_ranges) <- paste0("Peak_", seq_len(num_peaks))
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = mat_list,
    rowRanges = se_ranges,
    colData = S4Vectors::DataFrame(sample_sheet, row.names = sample_names)
  )

  message(sprintf("Success: 4D Portrait Matrix built (On-disk = %s).", on_disk))
  return(se)
}
