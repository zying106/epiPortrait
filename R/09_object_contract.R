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
  # Stored matrices retain the full call universe when a SummarizedExperiment
  # is subsequently reordered or subset by rows. Validate by stable IDs rather
  # than physical position; sample removal still requires recalling because it
  # changes group support summaries.
  calls <- S4Vectors::metadata(se)$superdomain_calls
  if (!is.null(calls)) {
    for (feat in names(calls)) {
      grps <- calls[[feat]]$groups
      if (is.null(grps)) next
      for (g in names(grps)) {
        m <- grps[[g]]$replicate_call_matrix
        if (is.null(m)) next
        if (is.null(rownames(m)) || anyNA(rownames(m)) ||
            anyDuplicated(rownames(m)) ||
            !all(rownames(se) %in% rownames(m))) {
          stop(sprintf("replicate_call_matrix for %s/%s must contain unique current domain IDs.",
                       feat, g))
        }
        if (!all(colnames(m) %in% colnames(se))) {
          stop(sprintf("replicate_call_matrix for %s/%s has sample columns absent from the object.",
                       feat, g))
        }
      }
    }
  }
  # Orthogonal evidence is also ID-addressed so row/sample reordering is safe.
  breadth_evidence <- S4Vectors::metadata(se)$breadth_domain_evidence
  if (!is.null(breadth_evidence)) {
    for (nm in c("evidence", "reason")) {
      m <- breadth_evidence[[nm]]
      if (!is.matrix(m)) {
        stop("breadth_domain_evidence$", nm,
             " must be a domain x sample matrix.")
      }
      if (is.null(rownames(m)) || is.null(colnames(m)) ||
          anyNA(rownames(m)) || anyNA(colnames(m)) ||
          anyDuplicated(rownames(m)) || anyDuplicated(colnames(m)) ||
          !all(rownames(se) %in% rownames(m)) ||
          !all(colnames(se) %in% colnames(m))) {
        stop("breadth_domain_evidence$", nm,
             " must be a domain x sample matrix containing unique ",
             "current domain and sample IDs.", call. = FALSE)
      }
    }
    if (!identical(dim(breadth_evidence$evidence),
                   dim(breadth_evidence$reason)) ||
        !identical(dimnames(breadth_evidence$evidence),
                   dimnames(breadth_evidence$reason))) {
      stop("Breadth evidence and reason matrices must have identical IDs.")
    }
    allowed <- c("Broad", "Typical", "PeakAbsent", "NoCall")
    if (!all(unique(as.character(breadth_evidence$evidence)) %in% allowed)) {
      stop("breadth_domain_evidence contains an unsupported evidence state.")
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
      if (is.null(rownames(combined)) || anyDuplicated(rownames(combined)) ||
          !all(rownames(se) %in% rownames(combined))) {
        stop("Stored replicate calls are not aligned to current domain IDs.")
      }
      combined <- combined[rownames(se), , drop = FALSE]
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


#' Extract Orthogonal Breadth Presence Evidence
#'
#' @description Returns the per-domain x per-replicate Breadth evidence stored
#'   by \code{call_super_domains(feature = "Breadth")}. Evidence is orthogonal
#'   to the canonical Breadth-Super call and has four states:
#'   \itemize{
#'     \item \code{Broad}: a uniquely assigned native peak exceeds the
#'           replicate's width cutoff;
#'     \item \code{Typical}: at least one uniquely assigned peak is present,
#'           but none is broad;
#'     \item \code{PeakAbsent}: the replicate-level width call is valid but no
#'           eligible native peak overlaps the shared domain;
#'     \item \code{NoCall}: the replicate is not callable, or overlapping
#'           peaks cannot be assigned uniquely at the requested overlap
#'           threshold.
#'   }
#'   \code{PeakAbsent} is an operational statement about the supplied peak
#'   calls, not proof that the biological chromatin domain disappeared.
#'   Evidence is matched by domain and sample IDs after subsetting or
#'   reordering. Stored group summaries are not recalculated by this accessor;
#'   rerun calling when changing the replicate composition of a group.
#'
#' @param se A SummarizedExperiment after Breadth calling.
#' @param group Character or NULL. Optionally restrict samples to one condition
#'   group using the \code{group_var} recorded in call provenance.
#' @param type Character. For matrix output, return \code{"evidence"} or its
#'   machine-readable \code{"reason"}. Ignored when \code{long = TRUE}.
#' @param long Logical. Return a long data.frame containing both Evidence and
#'   Reason instead of a matrix.
#' @return A domain x sample matrix, or a long data.frame with columns
#'   \code{Domain_ID}, \code{SampleID}, \code{Group}, \code{Evidence}, and
#'   \code{Reason}.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Breadth",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' head(get_breadth_evidence(se, group = "Control", long = TRUE))
#' @export
get_breadth_evidence <- function(se, group = NULL,
                                 type = c("evidence", "reason"),
                                 long = FALSE) {
  type <- match.arg(type)
  stored <- S4Vectors::metadata(se)$breadth_domain_evidence
  if (is.null(stored) || is.null(stored$evidence) || is.null(stored$reason)) {
    stop("No Breadth evidence found. Run call_super_domains(feature = ",
         "'Breadth') first.")
  }
  evidence <- stored$evidence
  reason <- stored$reason
  valid_ids <- function(ids) {
    !is.null(ids) && !anyNA(ids) && all(nzchar(ids)) &&
      !anyDuplicated(ids)
  }
  if (!is.matrix(evidence) || !is.matrix(reason) ||
      !identical(dim(evidence), dim(reason)) ||
      !identical(dimnames(evidence), dimnames(reason)) ||
      !valid_ids(rownames(evidence)) || !valid_ids(colnames(evidence)) ||
      !valid_ids(rownames(se)) || !valid_ids(colnames(se)) ||
      !all(rownames(se) %in% rownames(evidence)) ||
      !all(colnames(se) %in% colnames(evidence))) {
    stop("Stored Breadth evidence is not aligned to the object.")
  }
  evidence <- evidence[rownames(se), colnames(se), drop = FALSE]
  reason <- reason[rownames(se), colnames(se), drop = FALSE]

  meta <- as.data.frame(colData(se))
  group_var <- stored$group_var
  if (is.null(group_var)) {
    prov <- get_call_provenance(se, "Breadth")
    if (!is.null(prov$group_var)) group_var <- prov$group_var
  }
  keep <- seq_len(ncol(se))
  if (!is.null(group)) {
    if (length(group) != 1L || is.na(group)) {
      stop("group must be NULL or one non-missing value.")
    }
    if (is.null(group_var) || !group_var %in% colnames(meta)) {
      stop("No group variable was recorded for the Breadth evidence.")
    }
    keep <- which(as.character(meta[[group_var]]) == as.character(group))
    if (length(keep) == 0L) {
      stop("group '", group, "' not found in colData.")
    }
  }
  evidence <- evidence[, keep, drop = FALSE]
  reason <- reason[, keep, drop = FALSE]
  if (!long) return(if (type == "evidence") evidence else reason)

  group_values <- if (!is.null(group_var) && group_var %in% colnames(meta)) {
    as.character(meta[[group_var]])
  } else {
    rep(NA_character_, ncol(se))
  }
  if (nrow(se) == 0L || length(keep) == 0L) {
    return(data.frame(Domain_ID = character(), SampleID = character(),
                      Group = character(), Evidence = character(),
                      Reason = character()))
  }
  out <- do.call(rbind, lapply(seq_along(keep), function(j) {
    i <- keep[j]
    data.frame(
      Domain_ID = rownames(se),
      SampleID = colnames(se)[i],
      Group = group_values[i],
      Evidence = evidence[, j],
      Reason = reason[, j],
      stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}


#' Explain Why a Domain Call Is Uncertain
#'
#' @description For a per-group feature call, classifies the cause of each
#'   \code{NA} / \code{Uncertain} group call so the ambiguity is auditable:
#'   \itemize{
#'     \item \code{peak_absent_by_support_rule} (\code{Breadth}): enough
#'           callable replicates have no eligible native peak overlapping the
#'           domain. This is operational peak-call absence, not proof of
#'           biological disappearance.
#'     \item \code{mixed_peak_presence} (\code{Breadth}): callable replicates
#'           disagree between peak presence and peak absence.
#'     \item \code{overlap_without_unique_assignment} (\code{Breadth}): peaks
#'           overlap the domain, but none can be assigned uniquely at the
#'           selected overlap threshold.
#'     \item \code{insufficient_assessable_replicates} (\code{Breadth}): too
#'           few replicates distinguish peak presence from operational absence.
#'     \item \code{no_valid_signal_in_any_replicate} (signal features): no
#'           replicate produced a usable ranked distribution.
#'     \item \code{insufficient_valid_replicates}: some, but fewer than
#'           \code{min_valid_replicates}, replicates carried evidence (e.g. 1/3
#'           with a majority rule).
#'     \item \code{inflection_no_call_all_replicates}: every replicate's
#'           inflection was unreliable (\code{no_call}).
#'     \item \code{sharp_peak_regime} (\code{Breadth}): every replicate was
#'           withheld by the \code{min_broad_width_bp} sharp-peak guard (all
#'           eligible native peaks narrower than the floor), so no Broad
#'           evidence exists by design; use \code{feature = "Intensity"}.
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
#'   \code{N_Valid_Replicates}, \code{N_Assessable_Replicates},
#'   \code{Min_Valid_Replicates}, and \code{Cause}. For Breadth,
#'   \code{N_Valid_Replicates} counts Broad/Typical evidence, whereas
#'   \code{N_Assessable_Replicates} also counts operational
#'   \code{PeakAbsent} evidence.
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

  # Sharp-peak guard provenance (Breadth only): a replicate withheld by
  # min_broad_width_bp carries sharp_peak_regime = TRUE.
  rep_sharp <- NULL
  if (!is.null(prov$replicates)) {
    rep_sharp <- vapply(colnames(mat), function(s) {
      r <- prov$replicates[[s]]
      !is.null(r) && isTRUE(r$sharp_peak_regime)
    }, logical(1))
  }
  sharp_all <- !is.null(rep_sharp) && length(rep_sharp) > 0 &&
    all(rep_sharp)

  is_breadth <- identical(feature, "Breadth") ||
    grepl("Breadth", as.character(feature))

  cause <- rep(NA_character_, length(grp_call))
  n_assessable <- n_valid
  unc <- is.na(grp_call)
  cause[unc & all_no_call] <- if (sharp_all) {
    "sharp_peak_regime"
  } else {
    "inflection_no_call_all_replicates"
  }
  remaining <- unc & is.na(cause)
  if (is_breadth) {
    presence_col <- paste0("Breadth_PresenceStatus__", group)
    presence_status <- if (presence_col %in% colnames(rd)) {
      rd[[presence_col]]
    } else {
      rep(NA_character_, nrow(se))
    }
    cause[remaining & presence_status == "PeakAbsent"] <-
      "peak_absent_by_support_rule"
    cause[remaining & presence_status == "Mixed"] <- "mixed_peak_presence"

    evidence_mat <- get_breadth_evidence(
      se, group = group, type = "evidence", long = FALSE)
    reason_mat <- get_breadth_evidence(
      se, group = group, type = "reason", long = FALSE)
    n_assessable <- rowSums(evidence_mat != "NoCall")
    overlap_without_unique <- apply(reason_mat, 1, function(x) {
      any(x == "overlap_without_unique_assignment")
    })
    remaining <- unc & is.na(cause)
    cause[remaining & overlap_without_unique] <-
      "overlap_without_unique_assignment"
    remaining <- unc & is.na(cause)
    cause[remaining & n_valid > 0] <- "insufficient_valid_replicates"
    cause[remaining & n_assessable < mvr] <-
      "insufficient_assessable_replicates"
  } else {
    cause[remaining & n_valid == 0] <-
      "no_valid_signal_in_any_replicate"
    cause[remaining & n_valid > 0] <- "insufficient_valid_replicates"
  }

  data.frame(
    Domain_ID = rownames(se),
    Group_Call = grp_call,
    N_Valid_Replicates = n_valid,
    N_Assessable_Replicates = n_assessable,
    Min_Valid_Replicates = mvr,
    Cause = cause,
    stringsAsFactors = FALSE)
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
  enr <- S4Vectors::metadata(se)$enrichment
  enr_cmp <- S4Vectors::metadata(se)$enrichment_comparison
  list(
    domains = nrow(se),
    samples = ncol(se),
    assays = assayNames(se),
    combined_class_counts = classes,
    transition_columns = trans,
    superdomain_provenance_features =
      names(if (is.null(S4Vectors::metadata(se)$superdomain_calls)) list()
            else S4Vectors::metadata(se)$superdomain_calls),
    enrichment_results = names(if (is.null(enr)) list() else enr),
    enrichment_comparisons = names(if (is.null(enr_cmp)) list() else enr_cmp)
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
  # Orthogonal per-domain x per-replicate Breadth presence evidence.
  if (!is.null(S4Vectors::metadata(se)$breadth_domain_evidence)) {
    utils::write.table(
      get_breadth_evidence(se, long = TRUE),
      file.path(outdir, "calls", "breadth_domain_evidence.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  # Domain-gene annotation evidence: all three documented levels (per-domain
  # summary, per-domain-gene pair, raw relationship detail) plus provenance.
  if (!is.null(S4Vectors::metadata(se)$domain_gene_links)) {
    dir.create(file.path(outdir, "annotation"), showWarnings = FALSE)
    ann_summary <- S4Vectors::metadata(se)$annotation_summary
    if (!is.null(ann_summary)) {
      utils::write.table(ann_summary,
                         file.path(outdir, "annotation",
                                   "annotation_summary.tsv"),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
    ann_dedup <- S4Vectors::metadata(se)$domain_gene_links_dedup
    if (!is.null(ann_dedup)) {
      utils::write.table(ann_dedup,
                         file.path(outdir, "annotation",
                                   "domain_gene_links_dedup.tsv"),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
    utils::write.table(S4Vectors::metadata(se)$domain_gene_links,
                       file.path(outdir, "annotation", "domain_gene_links.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    prov_lines <- c()
    for (nm in c("annotation_provenance", "annotation_import_provenance",
                 "bedpe_provenance", "expression_provenance")) {
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
