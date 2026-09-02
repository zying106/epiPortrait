#' Combine Intensity and Breadth Super-Domain Calls
#'
#' @description Combines the Intensity-Super and Breadth-Super calls into the
#' final 2D architecture taxonomy (v1.0 design):
#' \code{Typical} / \code{Intensity-Super} / \code{Breadth-Super} /
#' \code{Dual-Super} / \code{Uncertain}. \code{Dual-Super} means the domain is
#' simultaneously high in integrated signal and unusually broad.
#'
#' @param se SummarizedExperiment after call_super_domains on Intensity and
#'   Breadth.
#' @param intensity_feature Character (default "Intensity").
#' @param width_feature Character. Breadth axis feature (default "Breadth").
#'   This is the peak-level native breadth call, NOT a static interval length.
#' @param call_col_suffix "Domain_Type" (consensus) or "Call__Control" (per_group).
#' @param group_var Character or NULL. If provided, the calls are combined per
#'   condition group (requires mode = "per_group"), producing
#'   \code{Combined_Class__<group>} columns used by
#'   \code{compare_superdomain_classes()}. If NULL (default), the single
#'   \code{call_col_suffix} columns are combined into \code{Combined_Domain_Class}.
#' @return se with Combined_Domain_Class column (or Combined_Class__<group>).
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          method = "tangent", log_transform = FALSE)
#' se <- call_super_domains(se, feature = "Breadth")
#' se <- combine_superdomain_calls(se)
#' table(SummarizedExperiment::rowData(se)$Combined_Domain_Class)
#' @export
combine_superdomain_calls <- function(se, intensity_feature = "Intensity",
                                       width_feature = "Breadth",
                                       call_col_suffix = "Domain_Type",
                                       group_var = NULL) {

  .combine_pair <- function(i_type, w_type, out_col) {
    i_na <- is.na(i_type)
    w_na <- is.na(w_type)
    i_s <- grepl("_Super", i_type) & !is.na(i_type)
    w_s <- grepl("_Super", w_type) & !is.na(w_type)
    cl <- rep("Typical", length(i_type))
    cl[i_s & !w_s] <- "Intensity-Super"
    cl[!i_s & w_s] <- "Breadth-Super"
    cl[i_s & w_s]  <- "Dual-Super"
    cl[i_na | w_na] <- "Uncertain"
    rowData(se)[[out_col]] <- cl
    se
  }

  # per-group mode: combine the Intensity and Breadth calls within each group
  # (P1-6). Produces Combined_Class__<group> columns so a combined-class
  # transition can be computed between conditions.
  if (!is.null(group_var)) {
    meta <- as.data.frame(colData(se))
    if (!group_var %in% colnames(meta)) stop("group_var not found in colData.")
    groups <- unique(meta[[group_var]])
    for (g in groups) {
      ic <- paste0(intensity_feature, "_Call__", g)
      wc <- paste0(width_feature, "_Call__", g)
      alt_i <- .resolve_assay(se, intensity_feature)
      alt_w <- .resolve_assay(se, width_feature)
      if (!ic %in% colnames(rowData(se))) ic <- paste0(alt_i, "_Call__", g)
      if (!wc %in% colnames(rowData(se))) wc <- paste0(alt_w, "_Call__", g)
      if (!ic %in% colnames(rowData(se)) || !wc %in% colnames(rowData(se))) {
        stop(sprintf("Run call_super_domains(mode='per_group') for %s and %s first.",
                     intensity_feature, width_feature))
      }
      se <- .combine_pair(rowData(se)[[ic]], rowData(se)[[wc]],
                          paste0("Combined_Class__", g))
    }
    return(se)
  }

  # single suffix (consensus / per_sample): current behaviour
  ic <- paste0(intensity_feature, "_", call_col_suffix)
  wc <- paste0(width_feature, "_", call_col_suffix)

  # Resolve feature aliases to the names actually used by call_super_domains().
  if (!ic %in% colnames(rowData(se))) {
    alt <- .resolve_assay(se, intensity_feature)
    ic <- paste0(alt, "_", call_col_suffix)
  }
  if (!wc %in% colnames(rowData(se))) {
    alt <- .resolve_assay(se, width_feature)
    wc <- paste0(alt, "_", call_col_suffix)
  }

  if (!ic %in% colnames(rowData(se))) stop("Run call_super_domains(", intensity_feature, ") first.")
  if (!wc %in% colnames(rowData(se))) stop("Run call_super_domains(", width_feature, ") first.")

  .combine_pair(rowData(se)[[ic]], rowData(se)[[wc]], "Combined_Domain_Class")
}


#' Compare Super-Domain Status Between Conditions
#'
#' @param se after call_super_domains(mode="per_group").
#' @param group_var Character.
#' @param ref_group Reference condition.
#' @param target_group Treatment condition.
#' @param feature Feature name (default "Intensity").
#' @param cutoff_scope Character. How the super/non-super threshold is defined
#'   The default \code{"relative"} calls each group independently, so a
#'   transition such as \code{Typical_to_Super} reflects a RELATIVE prominence
#'   change and must NOT be interpreted as an absolute signal gain.
#'   \itemize{
#'     \item \code{"relative"} (default): per-group cutoff; output is a
#'           relative class switch (rename the classes accordingly, do not call
#'           them Gain/Loss).
#'     \item \code{"reference"}: use the reference group's cutoff for both
#'           groups; closer to absolute gain/loss when BigWigs are
#'           quantitatively comparable.
#'     \item \code{"pooled"}: a common cutoff estimated on the pooled
#'           (ref + target) signal, applied to both groups.
#'   }
#'   For \code{"reference"} / \code{"pooled"} the feature must be a per-sample
#'   assay, and the two groups must be quantitatively comparable. These modes
#'   inherit the FULL replicate-support configuration of the primary call from
#'   stored provenance (support_rule, min_replicate_support, tie_policy,
#'   min_valid_replicates) so transition semantics match the group calls.
#' @return se with *_Transition column.
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          method = "tangent", log_transform = FALSE)
#' se <- compare_superdomains(se, group_var = "Condition",
#'                            ref_group = "Control", target_group = "Treatment")
#' table(SummarizedExperiment::rowData(se)$Intensity_Transition)
#' @export
compare_superdomains <- function(se, group_var = "Condition",
                                  ref_group, target_group,
                                  feature = "Intensity",
                                  cutoff_scope = c("relative", "reference", "pooled")) {
  cutoff_scope <- match.arg(cutoff_scope)
  meta <- as.data.frame(colData(se))
  if (!group_var %in% colnames(meta)) stop("group_var not found in colData.")
  idx_ref <- which(meta[[group_var]] == ref_group)
  idx_tgt <- which(meta[[group_var]] == target_group)
  if (length(idx_ref) == 0 || length(idx_tgt) == 0)
    stop("ref_group / target_group not found in colData.")

  rc <- paste0(feature, "_Call__", ref_group)
  tc <- paste0(feature, "_Call__", target_group)

  if (cutoff_scope == "relative") {
    # ---- per-group independent cutoffs (default) ---------------------------
    if (!rc %in% colnames(rowData(se))) stop("Run per_group mode first.")
    r_type <- rowData(se)[[rc]]
    t_type <- rowData(se)[[tc]]
  } else {
    # ---- common cutoff applied to raw signal values ------------------------
    if (!feature %in% assayNames(se)) {
      stop(sprintf("cutoff_scope '%s' requires a per-sample assay feature, but '%s' is not an assay.",
                   cutoff_scope, feature))
    }
    # P0-C: inherit the algorithmic settings of the primary call from the
    # stored provenance (method, log_transform_used, min_quality), instead of
    # silently falling back to the defaults (elbow + auto-log) used by
    # .call_super_domains_on_vector().
    prov <- S4Vectors::metadata(se)$superdomain_calls[[feature]]
    c_method <- if (!is.null(prov$method)) prov$method else "tangent"
    c_log <- if (!is.null(prov$log_transform_used)) prov$log_transform_used else FALSE
    c_minq <- if (!is.null(prov$min_quality)) prov$min_quality else 0.1
    c_quantile <- if (!is.null(prov$quantile_cutoff)) prov$quantile_cutoff else NULL

    mat <- assay(se, feature)
    r_raw <- rowMeans(mat[, idx_ref, drop = FALSE], na.rm = TRUE)
    t_raw <- rowMeans(mat[, idx_tgt, drop = FALSE], na.rm = TRUE)

    if (cutoff_scope == "reference") {
      cutoff <- .call_super_domains_on_vector(
        setNames(r_raw, rownames(se)), feature, c_quantile, c_log, FALSE,
        method = c_method, min_quality = c_minq)$cutoff_value
    } else { # pooled
      pooled <- rowMeans(cbind(r_raw, t_raw), na.rm = TRUE)
      cutoff <- .call_super_domains_on_vector(
        setNames(pooled, rownames(se)), feature, c_quantile, c_log, FALSE,
        method = c_method, min_quality = c_minq)$cutoff_value
    }
    if (!is.finite(cutoff)) stop("Could not estimate a common cutoff.")

    # P0-C: keep the transition replicate-aware. Each replicate is classified
    # against the SAME common cutoff, then the per-group support rule is
    # applied — so reference/pooled transitions use the same replicate-support
    # semantics as relative ones, rather than a single group-mean crossing.
    # P1-7: inherit the stored support_rule; fall back to "majority" ONLY when
    # provenance is absent. The previous condition was inverted and forced
    # "majority" whenever provenance carried a rule.
    support_rule <- if (is.null(prov) || is.null(prov$support_rule)) {
      "majority"
    } else {
      prov$support_rule
    }
    # Inherit the FULL replicate-support configuration from the primary call
    # (freeze review 2026-08-11): min_replicate_support (was hard-coded 0.5),
    # tie_policy (was always strict), and min_valid_replicates (NULL must stay
    # NULL so .replicate_support_call() re-resolves it from the per-group
    # replicate count — falling back to 1L changed semantics, especially for
    # unbalanced designs).
    support_fraction <- if (!is.null(prov) && !is.null(prov$min_replicate_support)) {
      prov$min_replicate_support
    } else {
      0.5
    }
    c_tie <- if (!is.null(prov) && !is.null(prov$tie_policy)) prov$tie_policy else "strict"
    m_valid <- if (!is.null(prov) && !is.null(prov$min_valid_replicates)) {
      prov$min_valid_replicates
    } else {
      NULL
    }
    r_type <- .replicate_support_call(
      .classify_vs_cutoff(mat[, idx_ref, drop = FALSE], cutoff, feature,
                          tie_policy = c_tie),
      feature, support_rule, support_fraction, m_valid)$group_type
    t_type <- .replicate_support_call(
      .classify_vs_cutoff(mat[, idx_tgt, drop = FALSE], cutoff, feature,
                          tie_policy = c_tie),
      feature, support_rule, support_fraction, m_valid)$group_type
  }

  r_na <- is.na(r_type)
  t_na <- is.na(t_type)
  r_s <- grepl("_Super", r_type) & !is.na(r_type)
  t_s <- grepl("_Super", t_type) & !is.na(t_type)

  tr <- rep(NA_character_, nrow(se))
  tr[!r_s & !t_s] <- "Persistent_Typical"
  tr[r_s & t_s] <- "Persistent_Super"
  tr[!r_s & t_s] <- "Typical_to_Super"
  tr[r_s & !t_s] <- "Super_to_Typical"
  tr[r_na | t_na] <- "Uncertain"  # propagate no-call (P0-8)

  # Relative transitions must not be labelled as absolute gain/loss (P0-4).
  col_name <- paste0(feature, "_Transition")
  if (cutoff_scope == "relative") {
    col_name <- paste0(feature, "_Relative_Transition")
    tr[tr == "Typical_to_Super"] <- "Relative_Prominence_Up"
    tr[tr == "Super_to_Typical"] <- "Relative_Prominence_Down"
  } else {
    col_name <- paste0(feature, "_Transition__", cutoff_scope)
    tr[tr == "Typical_to_Super"] <- "Gain"
    tr[tr == "Super_to_Typical"] <- "Loss"
  }
  rowData(se)[[col_name]] <- tr

  # Transition provenance (object-contract design). For reference/pooled modes
  # the FULL inherited replicate-support configuration is recorded so the
  # transition semantics are auditable (freeze review 2026-08-11).
  key <- paste0(ref_group, "_vs_", target_group)
  if (is.null(S4Vectors::metadata(se)$transitions)) {
    S4Vectors::metadata(se)$transitions <- list()
  }
  prov_t <- list(
    ref_group = ref_group,
    target_group = target_group,
    feature = feature,
    cutoff_scope = cutoff_scope,
    created_columns = col_name
  )
  if (cutoff_scope != "relative") {
    prov_t$support_rule <- support_rule
    prov_t$min_replicate_support <- support_fraction
    prov_t$tie_policy <- c_tie
    prov_t$min_valid_replicates <- m_valid
  }
  S4Vectors::metadata(se)$transitions[[key]] <- prov_t
  se
}


#' Compare Combined Super-Domain Classes Between Conditions
#'
#' @description Computes transitions of the combined taxonomy
#' (Intensity-Super / Breadth-Super / Dual-Super / Typical / Uncertain) between
#' two conditions, using the \code{Combined_Class__<group>} columns produced by
#' \code{combine_superdomain_calls(se, group_var = ...)}. This closes the core
#' loop: the combined architecture taxonomy is tracked across conditions
#' not just a single feature.
#'
#' @param se SummarizedExperiment after combine_superdomain_calls(group_var).
#' @param ref_group Reference condition.
#' @param target_group Treatment condition.
#' @return se with \code{Combined_Relative_Class_Transition} column.
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          method = "tangent", log_transform = FALSE)
#' se <- call_super_domains(se, feature = "Breadth",
#'                          mode = "per_group", group_var = "Condition")
#' se <- combine_superdomain_calls(se, group_var = "Condition")
#' se <- compare_superdomain_classes(se, ref_group = "Control",
#'                                   target_group = "Treatment")
#' table(SummarizedExperiment::rowData(se)$Combined_Relative_Class_Transition)
#' @export
compare_superdomain_classes <- function(se, ref_group, target_group) {
  rc <- paste0("Combined_Class__", ref_group)
  tc <- paste0("Combined_Class__", target_group)
  if (!rc %in% colnames(rowData(se)) || !tc %in% colnames(rowData(se))) {
    stop("Run combine_superdomain_calls(se, group_var = ...) first.")
  }
  r <- rowData(se)[[rc]]
  t <- rowData(se)[[tc]]
  # Combined calls encode abstention as the literal "Uncertain" whereas
  # feature-level no-calls use NA. Both representations are uncertainty and
  # must never become apparent biological transitions such as
  # Uncertain_to_Typical.
  r_uncertain <- is.na(r) | r == "Uncertain"
  t_uncertain <- is.na(t) | t == "Uncertain"
  tr <- ifelse(r_uncertain | t_uncertain, "Uncertain",
               ifelse(r == t, paste0("Persistent_", r),
                      paste0(r, "_to_", t)))
  # Relative architecture-state transition: combined classes come from
  # per-group independent cutoffs (review #10), so never label as absolute
  # gain/loss.
  rowData(se)[["Combined_Relative_Class_Transition"]] <- tr

  key <- paste0(ref_group, "_vs_", target_group)
  if (is.null(S4Vectors::metadata(se)$transitions)) {
    S4Vectors::metadata(se)$transitions <- list()
  }
  S4Vectors::metadata(se)$transitions[[key]] <- list(
    ref_group = ref_group,
    target_group = target_group,
    type = "combined_relative_class",
    cutoff_scope = "relative",
    created_columns = "Combined_Relative_Class_Transition"
  )
  se
}


#' Continuous Domain-Width Expansion / Contraction Descriptor
#'
#' @description Computes a continuous effect descriptor for domain breadth
#' between two condition groups.
#' For each shared domain it summarizes the mapped native width per group
#' (median NativeMaxPeakWidth across replicates) and reports:
#'   \itemize{
#'     \item \code{WidthDelta_bp}: median(target) - median(ref)
#'     \item \code{log2WidthRatio}: log2(median(target)/median(ref))
#'     \item \code{WidthDirection}: "Expansion" / "Contraction" / "Stable"
#'   }
#' This is a CONTINUOUS effect descriptor, not a significance test: with a
#' handful of replicates it reports the observed magnitude of width remodeling.
#' Statistical testing is left to downstream / external frameworks.
#'
#' @param se A SummarizedExperiment after build_portrait_matrix() with native
#'   peaks (requires the NativeMaxPeakWidth assay).
#' @param group_var Character. colData column for grouping (default "Condition").
#' @param ref_group,target_group Character. Reference and target condition
#'   groups.
#' @param effect_threshold Numeric or NULL. Minimum \code{|log2WidthRatio|}
#'   (fold-change on log2 scale) above which a domain is called "Expansion" or
#'   "Contraction"; below this it is "Stable" (default 0.585 ~ 1.5-fold). This
#'   is a DESCRIPTIVE effect-size label, not a statistical significance call.
#'   Set \code{NULL} to emit only the continuous columns (no hard
#'   classification).
#' @param min_replicates Integer. Minimum number of replicates per group with a
#'   finite width for the domain to be scored (default 1; NA otherwise).
#' @return A SummarizedExperiment with pair-specific rowData columns
#'   \code{WidthDelta_bp__<ref>_vs_<target>} and
#'   \code{log2WidthRatio__<ref>_vs_<target>} (and
#'   \code{WidthDirection__<ref>_vs_<target>} when a threshold is set), so
#'   multiple comparisons can be stored without overwriting each other.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' se <- compute_width_transition(example_se, ref_group = "Control",
#'                                target_group = "Treatment")
#' head(SummarizedExperiment::rowData(se)$WidthDelta_bp__Control_vs_Treatment)
#' @export
compute_width_transition <- function(se, group_var = "Condition",
                                     ref_group, target_group,
                                     effect_threshold = log2(1.5),
                                     min_replicates = 1L) {
  if (!"NativeMaxPeakWidth" %in% assayNames(se)) {
    stop("NativeMaxPeakWidth assay not found. Re-run build_portrait_matrix() ",
         "with per-sample peak files (peak_path).")
  }
  meta <- as.data.frame(colData(se))
  if (!group_var %in% colnames(meta)) stop("group_var not found in colData.")
  idx_ref <- which(meta[[group_var]] == ref_group)
  idx_tgt <- which(meta[[group_var]] == target_group)
  if (length(idx_ref) == 0 || length(idx_tgt) == 0) {
    stop("ref_group / target_group not found in colData.")
  }

  mat <- assay(se, "NativeMaxPeakWidth")
  ref_w <- mat[, idx_ref, drop = FALSE]
  tgt_w <- mat[, idx_tgt, drop = FALSE]
  med_ref <- apply(ref_w, 1, stats::median, na.rm = TRUE)
  med_tgt <- apply(tgt_w, 1, stats::median, na.rm = TRUE)
  n_ref <- rowSums(is.finite(ref_w))
  n_tgt <- rowSums(is.finite(tgt_w))
  valid <- n_ref >= min_replicates & n_tgt >= min_replicates

  delta <- rep(NA_real_, nrow(se))
  ratio <- rep(NA_real_, nrow(se))
  direction <- rep(NA_character_, nrow(se))
  ok <- valid & is.finite(med_ref) & is.finite(med_tgt) & med_ref > 0 & med_tgt > 0
  delta[ok] <- med_tgt[ok] - med_ref[ok]
  ratio[ok] <- log2(med_tgt[ok] / med_ref[ok])
  # WidthDirection is a DESCRIPTIVE thresholded label (1.5-fold by default),
  # not a statistical significance call. When effect_threshold = NULL only the
  # continuous effect columns are produced (no hard classification).
  if (!is.null(effect_threshold)) {
    direction[ok] <- ifelse(ratio[ok] > effect_threshold, "Expansion",
                            ifelse(ratio[ok] < -effect_threshold, "Contraction",
                                   "Stable"))
  }

  # Pair-specific column names so repeated comparisons (Control vs TreatmentA,
  # Control vs TreatmentB, ...) do not overwrite each other's rowData.
  key <- paste0(ref_group, "_vs_", target_group)
  col_delta <- paste0("WidthDelta_bp__", key)
  col_ratio <- paste0("log2WidthRatio__", key)
  col_dir   <- paste0("WidthDirection__", key)
  rowData(se)[[col_delta]] <- delta
  rowData(se)[[col_ratio]] <- ratio
  if (!is.null(effect_threshold)) {
    rowData(se)[[col_dir]] <- direction
  }

  if (is.null(S4Vectors::metadata(se)$width_transitions)) {
    S4Vectors::metadata(se)$width_transitions <- list()
  }
  S4Vectors::metadata(se)$width_transitions[[key]] <- list(
    ref_group = ref_group,
    target_group = target_group,
    group_var = group_var,
    width_assay = "NativeMaxPeakWidth",
    summary = "median native width per group (continuous descriptor)",
    effect_threshold_log2 = effect_threshold,
    note = "WidthDirection is a descriptive 1.5-fold threshold label, not a statistical significance call",
    created_columns = if (is.null(effect_threshold))
      c(col_delta, col_ratio) else c(col_delta, col_ratio, col_dir)
  )
  se
}


#' Get Mark-Specific Preset
#'
#' @param mark "H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3", etc.
#' @return List of preset parameters, including:
#'   \itemize{
#'     \item \code{label}: human-readable domain label
#'     \item \code{stitch_distance}: default stitching gap in bp. 0 / NULL is a
#'           SENTINEL meaning "do not automatically re-stitch at the
#'           preset-workflow level" — treat it as "skip \code{stitch_epi_peaks()}"
#'           (e.g. broad repressive marks such as H3K27me3 / H3K9me3, whose
#'           domains are already spatially clustered by the upstream broad-domain
#'           caller, and unknown/generic marks). Note that
#'           \code{stitch_epi_peaks(x, 0)} itself still merges
#'           overlapping/adjacent ranges; the sentinel must be handled by the
#'           caller (see the vignette preset workflow).
#'     \item \code{primary_feature}: ranking feature
#'     \item \code{mark_class}: "active" (H3K27ac, H3K4me3, H3K4me1) or
#'           "broad_repressive" (H3K27me3, H3K9me3)
#'     \item \code{taxonomy_style}: "super_domain" (active marks) or
#'           "repressive_remodeling" (broad repressive marks). This controls
#'           DISPLAY terminology only; the internal canonical classes
#'           (Intensity-Super / Breadth-Super / Dual-Super) are unchanged.
#'     \item \code{class_labels}: display aliases for the canonical taxonomy.
#'   }
#' @examples
#' get_mark_preset("H3K27ac")
#' get_mark_preset("H3K27me3")
#' @export
