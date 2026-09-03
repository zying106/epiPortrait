#' Validate an epiPortrait SummarizedExperiment Object
#'
#' @description Verifies the internal consistency invariants of an epiPortrait
#' object: unique row names, alignment between
#' assays and rowRanges/colData, IntervalWidth matching rowRanges widths,
#' presence of required assays, and consistent per-group class/support columns.
#'
#' @param se A SummarizedExperiment from \code{build_portrait_matrix()}.
#' @return Invisibly \code{TRUE} if all checks pass; otherwise stops with the
#'   first failing invariant.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' validate_epiportrait_object(example_se)
#' @export
validate_epiportrait_object <- function(se) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object.")
  }
  # rownames unique
  if (is.null(rownames(se)) || anyDuplicated(rownames(se))) {
    stop("rownames(se) must be present and unique (Domain_ID).")
  }
  # required canonical assays: Intensity and SignalDispersion are mandatory;
  # native geometry assays are optional (only present when peak_path provided).
  required <- c("Intensity", "SignalDispersion")
  missing <- setdiff(required, assayNames(se))
  if (length(missing) > 0) {
    stop(sprintf("Missing required epiPortrait assay(s): %s. ",
                 paste(missing, collapse = ", ")),
         "Run build_portrait_matrix().")
  }
  # assay row count aligns with rowRanges
  if (nrow(se) != length(rowRanges(se))) {
    stop("nrow(se) must equal length(rowRanges(se)).")
  }
  # IntervalWidth consistency (if present)
  if ("IntervalWidth" %in% colnames(rowData(se))) {
    iw <- as.numeric(rowData(se)$IntervalWidth)
    w <- GenomicRanges::width(rowRanges(se))
    if (any(!is.na(iw) & iw != w)) {
      stop("rowData(se)$IntervalWidth must equal width(rowRanges(se)).")
    }
  }
  # colData SampleID alignment
  if (ncol(se) != nrow(colData(se))) {
    stop("ncol(se) must equal nrow(colData(se)).")
  }
  if ("SampleID" %in% colnames(colData(se))) {
    if (!identical(colnames(se), as.character(colData(se)$SampleID))) {
      stop("colnames(se) must equal colData(se)$SampleID.")
    }
  }
  # If a per-group replicate_call_matrix is stored, verify that its rows and
  # row names match the domains and its columns are a subset of the samples.
  calls <- S4Vectors::metadata(se)$superdomain_calls
  if (!is.null(calls)) {
    for (feat in names(calls)) {
      grps <- calls[[feat]]$groups
      if (is.null(grps)) next
      for (g in names(grps)) {
        m <- grps[[g]]$replicate_call_matrix
        if (is.null(m)) next
        if (nrow(m) != nrow(se)) {
          stop(sprintf("replicate_call_matrix for %s/%s has %d rows, expected %d.",
                       feat, g, nrow(m), nrow(se)))
        }
        if (!identical(rownames(m), rownames(se))) {
          stop(sprintf("replicate_call_matrix for %s/%s rownames must match domain IDs.",
                       feat, g))
        }
        if (!all(colnames(m) %in% colnames(se))) {
          stop(sprintf("replicate_call_matrix for %s/%s has sample columns absent from the object.",
                       feat, g))
        }
      }
    }
  }
  invisible(TRUE)
}


#' Extract Per-Domain Results Table
#'
#' @description Returns a flat data.frame of domain-level results for export or
#' inspection. Includes genomic coordinates, IntervalWidth, and — when present
#' — per-group calls/support/rank, combined classes, transitions and group mean
#' quantitative values.
#'
#' @param se A SummarizedExperiment with super-domain calls.
#' @param group_var Character. Column in colData for group-mean feature columns.
#'   If NULL, no group-mean columns are added.
#' @param include_assay_means Logical. Add per-group mean of each dynamic assay
#'   (default TRUE when \code{group_var} is provided).
#' @return A \code{data.frame}.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' head(get_domain_results(example_se, group_var = "Condition"))
#' @export
get_domain_results <- function(se, group_var = "Condition",
                               include_assay_means = TRUE) {
  rd <- as.data.frame(rowData(se), optional = TRUE)
  coords <- as.data.frame(rowRanges(se))
  coords <- coords[, intersect(c("seqnames", "start", "end"), colnames(coords)),
                   drop = FALSE]
  if ("seqnames" %in% colnames(coords)) colnames(coords)[1] <- "chr"
  out <- data.frame(Domain_ID = rownames(se), coords, rd,
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!"IntervalWidth" %in% colnames(out)) {
    out$IntervalWidth <- GenomicRanges::width(rowRanges(se))
  }

  if (!is.null(group_var) && include_assay_means) {
    if (!group_var %in% colnames(colData(se))) {
      stop(sprintf("group_var '%s' not found in colData.", group_var))
    }
    meta <- as.data.frame(colData(se))
    for (grp in unique(meta[[group_var]])) {
      idx <- which(meta[[group_var]] == grp)
      for (feat in assayNames(se)) {
        out[[sprintf("%s_Mean__%s", feat, grp)]] <-
          rowMeans(assay(se, feat)[, idx, drop = FALSE], na.rm = TRUE)
      }
    }
  }

  # For every feature with per-group calls, append
  # <feature>_Uncertain_Cause__<group> so the
  # reason behind each NA/Uncertain call is directly auditable in the flat
  # table (see get_uncertain_cause()). Purely read-only; skipped on error.
  if (!is.null(group_var) && group_var %in% colnames(colData(se))) {
    meta <- as.data.frame(colData(se))
    call_cols <- grep("_Call__", colnames(rd), value = TRUE)
    feats <- unique(sub("_Call__.*$", "", call_cols))
    for (f in feats) {
      for (g in unique(meta[[group_var]])) {
        cc <- paste0(f, "_Call__", g)
        if (!cc %in% colnames(rd)) next
        out[[sprintf("%s_Uncertain_Cause__%s", f, g)]] <- tryCatch(
          get_uncertain_cause(se, feature = f, group = g,
                              group_var = group_var)$Cause,
          error = function(e) rep(NA_character_, nrow(se)))
      }
    }
  }
  out
}


#' Extract Sample-Level Results / QC Table
#'
#' @param se A SummarizedExperiment.
#' @return A \code{data.frame} with one row per sample (colData + QC).
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' get_sample_results(example_se)
#' @export
get_sample_results <- function(se) {
  meta <- as.data.frame(colData(se), optional = TRUE)
  meta
}


#' Extract Super-Domain Call Provenance
#'
#' @description Returns the provenance for a feature's super-domain calls
#' (cutoff, stability interval, quality, per-group replicate calls), stored in
#' \code{metadata(se)$superdomain_calls}.
#'
#' @param se A SummarizedExperiment after \code{call_super_domains}.
#' @param feature Character. Feature name (default "Intensity").
#' @return A list (the stored provenance), or NULL if unavailable.
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity", verbose = FALSE)
#' get_call_provenance(se, "Intensity")$method
#' @export
get_call_provenance <- function(se, feature = "Intensity") {
  calls <- S4Vectors::metadata(se)$superdomain_calls
  if (is.null(calls)) return(NULL)
  calls[[.resolve_assay(se, feature)]]
}


#' Extract Call Results as a Flat Table
#'
#' @description Returns per-domain call results for a feature in long or wide
#' form: for each condition group, the call, replicate support and rank.
#'
#' @param se A SummarizedExperiment after \code{call_super_domains}.
#' @param feature Character. Feature name (default "Intensity").
#' @param long Logical. If TRUE, returns long format (one row per domain x
#'   group); if FALSE, wide format with one column per group.
#' @return A \code{data.frame}.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' head(get_call_results(se, "Intensity"))
#' head(get_call_results(se, "Intensity", long = TRUE))
#' @export
get_call_results <- function(se, feature = "Intensity", long = FALSE) {
  feat <- .resolve_assay(se, feature)
  rd <- as.data.frame(rowData(se), optional = TRUE)
  call_cols <- grep(sprintf("^%s_Call__", feat), colnames(rd), value = TRUE)
  if (length(call_cols) == 0) {
    stop(sprintf("No per-group call columns found for '%s'. Run call_super_domains(mode='per_group').",
                 feat))
  }
  groups <- sub(sprintf("^%s_Call__", feat), "", call_cols)
  if (!long) {
    out <- data.frame(Domain_ID = rownames(se), check.names = FALSE)
    for (g in groups) {
      out[[sprintf("%s_Call__%s", feat, g)]] <- rd[[sprintf("%s_Call__%s", feat, g)]]
      if (sprintf("%s_Support__%s", feat, g) %in% colnames(rd)) {
        out[[sprintf("%s_Support__%s", feat, g)]] <- rd[[sprintf("%s_Support__%s", feat, g)]]
      }
      if (sprintf("%s_Rank__%s", feat, g) %in% colnames(rd)) {
        out[[sprintf("%s_Rank__%s", feat, g)]] <- rd[[sprintf("%s_Rank__%s", feat, g)]]
      }
    }
    return(out)
  }
  # long format
  parts <- lapply(groups, function(g) {
    data.frame(
      Domain_ID = rownames(se),
      Group = g,
      Call = rd[[sprintf("%s_Call__%s", feat, g)]],
      Support = if (sprintf("%s_Support__%s", feat, g) %in% colnames(rd))
        rd[[sprintf("%s_Support__%s", feat, g)]] else NA_real_,
      Rank = if (sprintf("%s_Rank__%s", feat, g) %in% colnames(rd))
        rd[[sprintf("%s_Rank__%s", feat, g)]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parts)
}


#' Extract Combined Class Results
#'
#' @param se A SummarizedExperiment after \code{combine_superdomain_calls}.
#' @return A \code{data.frame} of combined classes (and Combined_Domain_Class if
#'   present).
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity", verbose = FALSE)
#' se <- call_super_domains(se, feature = "Breadth", verbose = FALSE)
#' se <- combine_superdomain_calls(se)
#' head(get_combined_class_results(se))
#' @export
get_combined_class_results <- function(se) {
  rd <- as.data.frame(rowData(se), optional = TRUE)
  cc <- grep("^(Combined_Domain_Class|Combined_Class__)", colnames(rd), value = TRUE)
  if (length(cc) == 0) {
    stop("No combined class columns found. Run combine_superdomain_calls().")
  }
  out <- data.frame(Domain_ID = rownames(se), check.names = FALSE)
  for (c in cc) out[[c]] <- rd[[c]]
  out
}


#' Extract Transition Results
#'
#' @param se A SummarizedExperiment after \code{compare_superdomains} /
#'   \code{compare_superdomain_classes}.
#' @return A \code{data.frame} of transition columns.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' se <- compare_superdomains(se, group_var = "Condition",
#'                            ref_group = "Control", target_group = "Treatment")
#' head(get_transition_results(se))
#' @export
get_transition_results <- function(se) {
  rd <- as.data.frame(rowData(se), optional = TRUE)
  tc <- grep("Transition", colnames(rd), value = TRUE)
  if (length(tc) == 0) {
    stop("No transition columns found. Run compare_superdomains() / compare_superdomain_classes().")
  }
  out <- data.frame(Domain_ID = rownames(se), check.names = FALSE)
  for (c in tc) out[[c]] <- rd[[c]]
  out
}


#' Extract Per-Domain x Per-Replicate Call Matrix
#'
#' @description Exposes the per-domain x per-replicate super-domain evidence that
#' underlies the group-level replicate-support call. For
#' \code{mode = "per_group"} the matrix is read from stored provenance; for
#' \code{mode = "per_sample"} it is reconstructed from the
#' \code{<feature>_Call__<sample>} rowData columns. Group calls are
#' replicate-support aggregates; this matrix is the auditable evidence trail.
#'
#' @param se A SummarizedExperiment after \code{call_super_domains()}.
#' @param feature Character. Feature (default "Intensity").
#' @param group Character or NULL. For per_group calls, restrict to the samples
#'   of this condition group. NULL returns all samples.
#' @param long Logical. If TRUE, returns long format (domain_id, SampleID,
#'   Group, call); if FALSE (default), a domain x sample matrix.
#' @return A matrix (long = FALSE) or data.frame (long = TRUE).
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' m <- get_replicate_calls(se, feature = "Intensity", group = "Control")
#' head(m)
#' head(get_replicate_calls(se, feature = "Intensity", group = "Control",
#'                          long = TRUE))
#' @export
get_replicate_calls <- function(se, feature = "Intensity", group = NULL,
                                long = FALSE) {
  feat <- .resolve_assay(se, feature)
  meta <- as.data.frame(colData(se))

  # Resolve the group_var actually used by the call (from stored provenance,
  # defaulting to "Condition" for objects produced without per_group calls).
  .group_var_used <- function() {
    calls <- S4Vectors::metadata(se)$superdomain_calls
    if (!is.null(calls) && !is.null(calls[[feat]]$group_var)) {
      return(calls[[feat]]$group_var)
    }
    if ("Condition" %in% colnames(meta)) "Condition" else NULL
  }

  # 1. try the stored per-group replicate call matrices
  calls <- S4Vectors::metadata(se)$superdomain_calls
  mat <- NULL
  if (!is.null(calls) && !is.null(calls[[feat]]$groups)) {
    groups <- calls[[feat]]$groups
    mats <- lapply(names(groups), function(g) {
      m <- groups[[g]]$replicate_call_matrix
      if (is.null(m)) return(NULL)
      list(m = m, group = g)
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) > 0) {
      # bind columns across groups (samples are distinct within a group)
      mats_ordered <- mats[order(vapply(mats, function(x) x$group, character(1)))]
      cols <- unlist(lapply(mats_ordered, function(x) colnames(x$m)))
      combined <- do.call(cbind, lapply(mats_ordered, function(x) x$m))
      colnames(combined) <- cols
      if (!is.null(group)) {
        g_samples <- rownames(meta)[meta[[calls[[feat]]$group_var]] == group]
        g_samples <- intersect(g_samples, colnames(combined))
        combined <- combined[, g_samples, drop = FALSE]
      } else {
        # preserve the original sample order of the object
        # rather than the alphabetical group order used during cbind.
        combined <- combined[, intersect(colnames(se), colnames(combined)),
                             drop = FALSE]
      }
      mat <- combined
    }
  }
  # 2. fallback: per_sample rowData columns
  if (is.null(mat)) {
    rd <- colnames(rowData(se))
    s_cols <- grep(sprintf("^%s_Call__", feat), rd, value = TRUE)
    if (length(s_cols) > 0) {
      samples <- sub(sprintf("^%s_Call__", feat), "", s_cols)
      if (!is.null(group)) {
        gv <- .group_var_used()
        if (is.null(gv) || !gv %in% colnames(meta)) {
          stop("Cannot resolve the group variable used by the call to filter samples.")
        }
        keep <- samples %in% rownames(meta)[meta[[gv]] == group]
        s_cols <- s_cols[keep]; samples <- samples[keep]
      }
      if (length(s_cols) > 0) {
        # The per_sample fallback also returns a matrix, consistent with
        # the per_group branch and the documented return type).
        mat <- as.matrix(as.data.frame(rowData(se)[, s_cols, drop = FALSE]))
        colnames(mat) <- samples
        rownames(mat) <- rownames(se)
      }
    }
  }
  if (is.null(mat) || ncol(mat) == 0) {
    # Distinguish "no calls at all" from "group absent".
    if (!is.null(group) && !group %in% unique(as.character(meta[[.group_var_used()]]))) {
      stop(sprintf("group '%s' not found in colData.", group))
    }
    stop("No replicate call matrix found. Run call_super_domains() with ",
         "mode = 'per_group' or mode = 'per_sample' first.")
  }
  if (!long) return(mat)

  # long format: use the ACTUAL group_var used by the call (not a hard-coded
  # "Condition"), so custom group_var (e.g. "Subtype") is
  # respected. Column is named "Group".
  gv <- .group_var_used()
  group_vals <- if (!is.null(gv) && gv %in% colnames(meta)) {
    stats::setNames(meta[[gv]], rownames(meta))
  } else NULL
  long_df <- do.call(rbind, lapply(colnames(mat), function(s) {
    data.frame(domain_id = rownames(mat),
               SampleID = s,
               Group = if (!is.null(group_vals)) unname(group_vals[s]) else NA_character_,
               call = mat[, s],
               stringsAsFactors = FALSE)
  }))
  rownames(long_df) <- NULL
  long_df
}


#' Explain Why a Domain Call Is Uncertain
#'
#' @description For a per-group feature call, classifies the cause of each
#'   \code{NA} / \code{Uncertain} group call so the ambiguity is auditable:
#'   \itemize{
#'     \item \code{no_native_peak_in_condition} (\code{Breadth}): the domain has
#'           no native peak evidence in any replicate of the group — e.g. it is
#'           condition-specific (present only in the other condition's universe).
#'     \item \code{no_valid_signal_in_any_replicate} (signal features): no
#'           replicate produced a usable ranked distribution.
#'     \item \code{insufficient_valid_replicates}: some, but fewer than
#'           \code{min_valid_replicates}, replicates carried evidence (e.g. 1/3
#'           with a majority rule).
#'     \item \code{inflection_no_call_all_replicates}: every replicate's
#'           inflection was unreliable (\code{no_call}).
#'   }
#'   Domains whose group call is \emph{not} \code{NA} get \code{cause = NA}.
#'   This reads stored provenance only and never re-computes calls.
#'
#' @param se A SummarizedExperiment after \code{call_super_domains(mode =
#'   "per_group")}.
#' @param feature Character. Feature (default "Breadth").
#' @param group Character. Condition group.
#' @param group_var Character or NULL. Column used for grouping; resolved from
#'   stored provenance when NULL (default "Condition").
#' @return A data.frame with columns \code{Domain_ID}, \code{Group_Call},
#'   \code{N_Valid_Replicates}, \code{Min_Valid_Replicates}, \code{Cause}.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Breadth",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' get_uncertain_cause(se, feature = "Breadth", group = "Control")
#' @export
get_uncertain_cause <- function(se, feature = "Breadth", group,
                                group_var = NULL) {
  prov <- get_call_provenance(se, feature)
  if (is.null(prov)) {
    stop("No call provenance found for feature '", feature,
         "'. Run call_super_domains(mode='per_group') first.")
  }
  meta <- as.data.frame(colData(se))
  if (is.null(group_var)) {
    group_var <- if (!is.null(prov$group_var)) prov$group_var else "Condition"
  }
  if (!group_var %in% colnames(meta)) {
    stop("group_var '", group_var, "' not found in colData.")
  }
  if (!group %in% unique(meta[[group_var]])) {
    stop("group '", group, "' not found in colData.")
  }

  mat <- get_replicate_calls(se, feature = feature, group = group)
  n_reps <- ncol(mat)
  n_valid <- rowSums(!is.na(mat))

  call_col <- paste0(.resolve_assay(se, feature), "_Call__", group)
  rd <- as.data.frame(rowData(se))
  if (!call_col %in% colnames(rd)) {
    stop("Call column '", call_col,
         "' not found. Run call_super_domains(mode='per_group') first.")
  }
  grp_call <- rd[[call_col]]

  # min_valid / support rule from stored provenance (never recomputed).
  sr <- if (!is.null(prov$support_rule)) prov$support_rule else "majority"
  mvr <- prov$min_valid_replicates
  if (is.null(mvr)) {
    frac <- if (!is.null(prov$min_replicate_support)) prov$min_replicate_support else 0.5
    mvr <- switch(sr,
      majority = floor(n_reps / 2) + 1L,
      all      = n_reps,
      fraction = ceiling(frac * n_reps))
  }

  # Per-replicate call status: Breadth stores prov$replicates; per-group signal
  # features store prov$groups[[group]]$replicate_calls.
  rep_status <- NULL
  if (!is.null(prov$replicates)) {
    rep_status <- vapply(colnames(mat), function(s) {
      r <- prov$replicates[[s]]
      if (is.null(r)) NA_character_ else if (!is.null(r$call_status)) r$call_status else "called"
    }, character(1))
  } else if (!is.null(prov$groups[[group]]$replicate_calls)) {
    rep_status <- vapply(prov$groups[[group]]$replicate_calls, function(r) {
      if (!is.null(r$call_status)) r$call_status else "called"
    }, character(1))
  }
  all_no_call <- !is.null(rep_status) &&
    length(rep_status) > 0 && all(rep_status == "no_call", na.rm = TRUE)

  is_breadth <- identical(feature, "Breadth") ||
    grepl("Breadth", as.character(feature))

  cause <- rep(NA_character_, length(grp_call))
  unc <- is.na(grp_call)
  cause[unc & all_no_call] <- "inflection_no_call_all_replicates"
  cause[unc & !all_no_call & n_valid == 0] <-
    if (is_breadth) "no_native_peak_in_condition"
    else "no_valid_signal_in_any_replicate"
  cause[unc & !all_no_call & n_valid > 0] <- "insufficient_valid_replicates"

  data.frame(Domain_ID = rownames(se), Group_Call = grp_call,
             N_Valid_Replicates = n_valid, Min_Valid_Replicates = mvr,
             Cause = cause, stringsAsFactors = FALSE)
}


#' Summarize an epiPortrait Object
#'
#' @param se A SummarizedExperiment.
#' @return A list with a short text summary of the object contents.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' summarize_epiportrait_object(example_se)
#' @export
summarize_epiportrait_object <- function(se) {
  rd <- as.data.frame(rowData(se), optional = TRUE)
  classes <- NULL
  if ("Combined_Domain_Class" %in% colnames(rd)) {
    classes <- table(rd$Combined_Domain_Class)
  } else {
    cc <- grep("^Combined_Class__", colnames(rd), value = TRUE)
    if (length(cc) > 0) classes <- table(rd[[cc[1]]])
  }
  trans <- grep("Transition", colnames(rd), value = TRUE)
  list(
    domains = nrow(se),
    samples = ncol(se),
    assays = assayNames(se),
    combined_class_counts = classes,
    transition_columns = trans,
    superdomain_provenance_features =
      names(if (is.null(S4Vectors::metadata(se)$superdomain_calls)) list()
            else S4Vectors::metadata(se)$superdomain_calls)
  )
}


#' Export Complete epiPortrait Results to Disk
#'
#' @description Writes a self-contained results directory: flat TSV tables,
#' per-assay matrices, per-call tables, an object manifest, and the complete
#' SummarizedExperiment as RDS.
#'
#' @param se A SummarizedExperiment.
#' @param outdir Character. Output directory (created if missing).
#' @param save_object Logical. Also save \code{se} as RDS (default TRUE).
#' @param group_var Character. Passed to \code{get_domain_results} for group-mean
#'   columns.
#' @return Invisibly the path to the output directory.
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity", verbose = FALSE)
#' out <- export_epiportrait_results(se, outdir = tempfile("epi_export"))
#' @export
export_epiportrait_results <- function(se, outdir = "epiPortrait_results",
                                       save_object = TRUE,
                                       group_var = "Condition") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  dir.create(file.path(outdir, "assays"), showWarnings = FALSE)
  dir.create(file.path(outdir, "calls"), showWarnings = FALSE)

  # Domain results
  utils::write.table(get_domain_results(se, group_var = group_var),
                     file.path(outdir, "domain_results.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  # Sample results
  utils::write.table(get_sample_results(se),
                     file.path(outdir, "sample_results.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  # Combined classes (if present)
  if (any(grepl("^Combined", colnames(rowData(se))))) {
    utils::write.table(get_combined_class_results(se),
                       file.path(outdir, "combined_classes.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  # Transitions (if present)
  if (any(grepl("Transition", colnames(rowData(se))))) {
    utils::write.table(get_transition_results(se),
                       file.path(outdir, "transition_results.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  # Domain-gene annotation evidence: long-format links and provenance.
  if (!is.null(S4Vectors::metadata(se)$domain_gene_links)) {
    dir.create(file.path(outdir, "annotation"), showWarnings = FALSE)
    utils::write.table(S4Vectors::metadata(se)$domain_gene_links,
                       file.path(outdir, "annotation", "domain_gene_links.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    prov_lines <- c()
    for (nm in c("annotation_provenance", "bedpe_provenance",
                 "expression_provenance")) {
      p <- S4Vectors::metadata(se)[[nm]]
      if (!is.null(p)) {
        prov_lines <- c(prov_lines, paste0("==", nm, "=="),
                        utils::capture.output(dput(p)))
      }
    }
    if (length(prov_lines) > 0) {
      writeLines(prov_lines, file.path(outdir, "annotation", "provenance.txt"))
    }
  }
  # Assay matrices
  for (a in assayNames(se)) {
    utils::write.table(as.matrix(assay(se, a)),
                       file.path(outdir, "assays", paste0(a, ".tsv")),
                       sep = "\t", quote = FALSE, row.names = TRUE)
  }
  # Per-feature call tables (per_group)
  rd <- colnames(rowData(se))
  feats <- assayNames(se)
  for (f in feats) {
    if (any(grepl(sprintf("^%s_Call__", f), rd))) {
      utils::write.table(get_call_results(se, feature = f),
                         file.path(outdir, "calls", paste0(f, "_calls.tsv")),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  # Manifest
  manifest <- c(
    sprintf("epiPortrait object export"),
    sprintf("timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("package version: %s",
            if (requireNamespace("epiPortrait", quietly = TRUE))
              as.character(utils::packageVersion("epiPortrait")) else "source"),
    sprintf("domains: %d", nrow(se)),
    sprintf("samples: %d", ncol(se)),
    sprintf("assays: %s", paste(assayNames(se), collapse = ", ")),
    sprintf("output files:")
  )
  files <- list.files(outdir, recursive = TRUE)
  manifest <- c(manifest, paste0("  ", files))
  writeLines(manifest, file.path(outdir, "object_manifest.txt"))

  if (save_object) {
    saveRDS(se, file.path(outdir, "epiPortrait_object.rds"))
  }
  invisible(outdir)
}
