#' Call Super-Domains with Multi-Mode Support
#'
#' @description Identifies super-domains by ranking domains on a feature and
#' separating extreme ("super") from typical domains via a hockey-stick
#' inflection, with replicate-aware support aggregation.
#' Supports three modes:
#' \itemize{
#'   \item \code{"global_consensus"} (default; alias \code{"consensus"}):
#'     calls per-sample, computes replicate support per condition group, and
#'     requires the support rule to be met in EVERY group. Identifies domains
#'     reproducibly super across all conditions (persistent/global consensus).
#'   \item \code{"per_group"}: EACH replicate is called independently using its
#'     own ranked feature distribution; replicate-level calls are then
#'     aggregated within each condition using the selected support rule.
#'     Group-mean ranks are stored for visualization only and do not determine
#'     the group call. Stores one call column per group. Recommended for the
#'     main analysis of multi-condition studies. The per-domain x per-replicate
#'     call matrix is available via \code{get_replicate_calls()}.
#'   \item \code{"per_sample"}: calls each sample independently, stores all
#'     call columns.
#' }
#'
#' @param se A SummarizedExperiment from \code{build_portrait_matrix()}.
#' @param feature Character. Ranking feature.
#'   \itemize{
#'     \item Per-sample assays: \code{"Intensity"} (integrated signal
#'           magnitude), or \code{"SignalDispersion"} (within-domain signal
#'           architecture).
#'     \item \code{"Breadth"}: peak-level Breadth-Super calling. Requires the
#'           object to carry sample-specific native peak calls (see
#'           \code{peak_path} in \code{build_portrait_matrix}; stored in
#'           \code{metadata(se)$native_peaks}). Each replicate's genome-wide
#'           eligible native PeakWidth distribution is ranked and cut by an
#'           elbow/inflection; broad peaks are then mapped to the shared domains
#'           by unique assignment and aggregated across replicates. This
#'           follows the v1.0 design (Breadth-Super calling is decoupled from
#'           consensus construction).
#'   }
#'   The static genomic interval length is available as \code{"IntervalWidth"}
#'   (alias \code{"Width"}; read from rowData, constant across samples — use
#'   only for static broadness ranking, never for per-sample calls).
#' @param mode Character. "global_consensus" (default), "per_group", or
#'   "per_sample". The alias "consensus" is also accepted.
#' @param group_var Character. Column in colData for grouping
#'   (used when mode = "per_group" or "global_consensus").
#' @param method Character. Inflection-detection heuristic:
#'   \itemize{
#'     \item \code{"elbow"} (default): point with maximum perpendicular distance
#'           to the first-last ranked line (classic knee detection).
#'     \item \code{"tangent"}: a ROSE-inspired tangent-optimization inflection
#'           (Whyte et al., 2013, Cell 153:307-319). Slides a line whose slope is
#'           \code{(max - min)/n} and minimizes the number of points below it.
#'           Note that this is a geometric variant of ROSE's
#'           \code{calculate_cutoff()} step — not a bit-exact reproduction — and
#'           implements only an inflection step, not the full ROSE pipeline.
#'   }
#'   Ignored when \code{quantile_cutoff} is provided.
#' @param log_transform Logical or NULL. Whether to \code{log10(x + 1)}-transform
#'   the feature before ranking. \code{NULL} (default) auto-selects per feature
#'   (Intensity and SignalDispersion are log-transformed, matching epiPortrait
#'   conventions). Set \code{FALSE} to rank on the raw feature scale (needed for
#'   a ROSE-style tangent benchmark on untransformed signals).
#'
#' @details
#' **Optional H3K27ac / ROSE benchmark.** If the goal is a head-to-head
#' comparison with the original ROSE method (Whyte et al., 2013, Cell
#' 153:307-319), use \code{method = "tangent"} with
#' \code{log_transform = FALSE} (raw signal, ROSE-style scale). Because the
#' tangent implementation is a geometric variant rather than a bit-exact port
#' of ROSE's \code{calculate_cutoff()}, validate any ROSE-reproduction numbers
#' against the reference implementation before drawing quantitative
#' conclusions. This is an optional benchmark setting, not the package's primary
#' workflow; the core value of epiPortrait is the Intensity x native-Breadth
#' decomposition with replicate support. The log-transformed default is retained
#' as the epiPortrait convention because it stabilises the skewed intensity
#' distribution for the elbow heuristic and downstream analysis.
#'
#' @param quantile_cutoff Numeric or NULL. If provided (e.g. 0.95), bypasses
#'   inflection detection and uses an explicit top fraction instead. For
#'   \code{"Intensity"} the cutoff is the \code{quantile_cutoff} quantile of
#'   the ranked domain values. For \code{"Breadth"} it is the
#'   \code{quantile_cutoff} quantile of each replicate's genome-wide
#'   eligible native PeakWidth distribution — i.e. the top
#'   \code{100*(1-quantile_cutoff)}\% widest native peaks are Broad (the same
#'   "top fraction" convention as Intensity: \code{0.95} = top 5\%). This is an
#'   explicit user opt-in (a common broad-domain practice) with no quality
#'   gate; it never auto-falls-back to, or from, the data-driven inflection
#'   path, so the no-call semantics of the default path are unchanged. Must be
#'   in (0, 1).
#' @param min_replicate_support Numeric. Minimum support fraction used when
#'   \code{support_rule = "fraction"} (default 0.5).
#' @param support_rule Character. Replicate-support rule for per-group /
#'   consensus calls:
#'   \itemize{
#'     \item \code{"majority"} (default): strictly more than half of the
#'           replicates must call Super (n=2 needs 2, n=3 needs 2, n=4 needs 3).
#'     \item \code{"all"}: every replicate must call Super.
#'     \item \code{"fraction"}: at least \code{min_replicate_support} fraction.
#'   }
#'   no-call replicates never count as Super and never inflate support.
#' @param min_valid_replicates Integer or NULL. Minimum number of valid
#'   (non-NA) replicates required before a domain can be called Super in a
#'   group. When \code{NULL} (default), it is set to the number required by the
#'   selected support rule (majority: floor(n/2)+1; all: n; fraction:
#'   ceiling(fraction*n)), so a "Super + no_call" pair in an n=2 group is
#'   Uncertain rather than Typical. Below the requirement the group call is
#'   \code{NA} / \code{Uncertain}. Note the deliberate asymmetry with
#'   \code{support_rule = "fraction"}: because the user explicitly chose how
#'   much support suffices, a low threshold (e.g. 0.5 with n=2 resolves the
#'   requirement to 1) can call a "Super + no_call" pair Super — the no-call
#'   replicate never counts as Super, it simply does not lower the bar the user
#'   set. Use the stricter majority rule if that behaviour is unwanted.
#' @param min_quality Numeric. Minimum inflection quality for a reliable call
#'   (see \code{find_hockey_inflection}; default 0.1).
#' @param tie_policy Character. How to classify domains whose value equals the
#'   cutoff. \code{"strict"} (default) requires value \code{>} cutoff,
#'   matching the original ROSE super-enhancer selection
#'   (\code{rankBy_vector > cutoff}); this avoids inflating the super set when
#'   discrete features tie exactly at the cutoff.
#'   \code{"inclusive"} uses \code{>=} cutoff (legacy behaviour).
#' @param n_bootstrap Integer or NULL. Bootstrap resamples for the cutoff
#'   stability interval and success rate (default NULL; set e.g. 200 to
#'   enable). The interval reflects resampling stability over the fixed
#'   ranked set, not a classical iid confidence interval. Its information
#'   content shrinks as the ranked set grows (interval width ~ 1/sqrt(n)); at
#'   genome-wide scale the primary cutoff-stability check is per-replicate
#'   cutoff consistency via \code{get_call_provenance()$replicates}. For
#'   \code{"Intensity"} the ranked set is the candidate domains; for
#'   \code{"Breadth"} it is each replicate's eligible native PeakWidth
#'   distribution (a per-replicate interval is stored in
#'   \code{prov$replicates[[s]]}), and it applies only to the data-driven
#'   inflection path (ignored with a warning when \code{quantile_cutoff} is
#'   set).
#' @param seed Integer or NULL. Random seed for the bootstrap (default NULL).
#'   A NULL seed leaves the resampling truly random (recommended for general
#'   use); pass a fixed integer to make the bootstrap reproducible.
#' @param min_peak_overlap_fraction Numeric. Only used when
#'   \code{feature = "Breadth"}. Minimum fraction of a native peak's genomic
#'   extent that must fall inside a shared domain for the peak's Broad evidence
#'   to be assigned to that domain (default 0.5). A native peak is assigned to
#'   at most ONE shared domain (unique assignment by maximum overlap); peaks
#'   not reaching the threshold are \code{Unmapped}, and tied maxima are
#'   \code{Ambiguous}. Prevents edge-overlap double counting of Broad evidence.
#' @param valid_chroms Character or NULL. Only used when
#'   \code{feature = "Breadth"}. Optional vector of allowed chromosomes for the
#'   genome-wide eligible native peak set (e.g. \code{c("chr1", ..., "chr22")}
#'   or the standard chromosomes of the assembly). NULL (default) keeps all
#'   finite-length seqlevels.
#' @param verbose Logical. Print messages (default TRUE).
#'
#' @return A SummarizedExperiment with rowData columns for each call.
#'   Call provenance (feature, method, requested/used log_transform, cutoff,
#'   cutoff_stability_interval, bootstrap_success_rate, quality_score,
#'   call_status, n_super) is stored in
#'   \code{metadata(se)$superdomain_calls}.
#' @import SummarizedExperiment
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          method = "tangent", log_transform = FALSE)
#' table(SummarizedExperiment::rowData(se)$Intensity_Domain_Type)
#' @export
call_super_domains <- function(se, feature = "Intensity",
                                mode = c("global_consensus", "consensus",
                                         "per_group", "per_sample"),
                                group_var = "Condition",
                                method = "elbow",
                                log_transform = NULL,
                                quantile_cutoff = NULL,
                                min_replicate_support = 0.5,
                                support_rule = c("majority", "all", "fraction"),
                                min_valid_replicates = NULL,
                                min_quality = 0.1,
                                tie_policy = c("strict", "inclusive"),
                                n_bootstrap = NULL,
                                seed = NULL,
                                min_peak_overlap_fraction = 0.5,
                                valid_chroms = NULL,
                                verbose = TRUE) {
  mode <- match.arg(mode)
  # "consensus" is retained as a backward-compatible alias for
  # "global_consensus" (P1-9): the mode identifies domains reproducibly super
  # in EVERY condition group, i.e. a global/persistent consensus.
  if (mode == "consensus") mode <- "global_consensus"
  method <- match.arg(method, c("elbow", "tangent"))
  support_rule <- match.arg(support_rule)
  tie_policy <- match.arg(tie_policy)
  # Freeze review 2026-08-11: min_replicate_support must be a finite fraction
  # in (0, 1] when support_rule = "fraction". 0 would set required = 0 (every
  # domain Super); > 1 or NA silently breaks ceiling(fraction * n).
  if (support_rule == "fraction" &&
      (length(min_replicate_support) != 1L || !is.numeric(min_replicate_support) ||
       !is.finite(min_replicate_support) ||
       min_replicate_support <= 0 || min_replicate_support > 1)) {
    stop("min_replicate_support must be a finite number in (0, 1] ",
         "when support_rule = 'fraction'.")
  }
  # quantile_cutoff (both Intensity and Breadth routes): an explicit top
  # fraction must be a finite number in (0, 1). 0/1 would collapse the cut to
  # the min/max of the ranked values; NA or out-of-range silently breaks
  # stats::quantile.
  if (!is.null(quantile_cutoff) &&
      (length(quantile_cutoff) != 1L || !is.numeric(quantile_cutoff) ||
       !is.finite(quantile_cutoff) ||
       quantile_cutoff <= 0 || quantile_cutoff >= 1)) {
    stop("quantile_cutoff must be a finite number in (0, 1) ",
         "(e.g. 0.95 = top 5%).")
  }

  # ---- Breadth route: peak-level native PeakWidth calling (v1.0 design) -----
  # Breadth-Super is NOT ranked on a domain x sample assay. Broad/typical is
  # decided per replicate on the genome-wide eligible native PeakWidth
  # distribution (elbow/inflection), THEN mapped to the shared domains by
  # unique assignment and aggregated across replicates. Consensus construction
  # never redefines the native broad status.
  if (feature == "Breadth") {
    # P1-5: min_peak_overlap_fraction must be a finite fraction in [0, 1].
    if (length(min_peak_overlap_fraction) != 1L ||
        !is.numeric(min_peak_overlap_fraction) ||
        !is.finite(min_peak_overlap_fraction) ||
        min_peak_overlap_fraction < 0 || min_peak_overlap_fraction > 1) {
      stop("min_peak_overlap_fraction must be a finite number in [0, 1].")
    }
    # quantile_cutoff for Breadth is an explicit user opt-in: the top fraction
    # is taken from each replicate's genome-wide eligible native PeakWidth
    # distribution (same population as the inflection path), giving a bp
    # cutoff. It NEVER auto-falls-back to (or from) the data-driven
    # inflection: when quantile_cutoff is NULL the elbow/tangent path runs with
    # its no-call semantics unchanged.
    if (!is.null(quantile_cutoff) && !is.null(n_bootstrap) && n_bootstrap > 0) {
      warning("n_bootstrap applies only to the data-driven inflection path; ",
              "it is ignored when quantile_cutoff is set.")
    }
    return(.call_breadth_super_domains(
      se, mode = mode, group_var = group_var, method = method,
      quantile_cutoff = quantile_cutoff, min_quality = min_quality,
      tie_policy = tie_policy, n_bootstrap = n_bootstrap, seed = seed,
      support_rule = support_rule,
      min_replicate_support = min_replicate_support,
      min_valid_replicates = min_valid_replicates,
      min_peak_overlap_fraction = min_peak_overlap_fraction,
      valid_chroms = valid_chroms, verbose = verbose))
  }

  # Feature may be a per-sample assay OR a static rowData column (e.g.
  # IntervalWidth, P1-5). Resolve aliases (Width -> IntervalWidth) against
  # either source.
  if (!feature %in% assayNames(se) && !feature %in% colnames(rowData(se))) {
    alt <- .resolve_assay(se, feature)
    if (alt != feature) feature <- alt
  }
  if (!feature %in% assayNames(se) && !feature %in% colnames(rowData(se))) {
    stop("Feature not found in assays or rowData: ", feature)
  }
  meta <- as.data.frame(colData(se))

  # Feature matrix: per-sample assays, OR a static rowData column such as
  # IntervalWidth (constant across samples, P1-5).
  if (feature %in% assayNames(se)) {
    mat <- assay(se, feature)
  } else if (feature %in% colnames(rowData(se))) {
    v <- as.numeric(rowData(se)[[feature]])
    mat <- matrix(rep(v, ncol(se)), ncol = ncol(se))
    colnames(mat) <- colnames(se)
  } else {
    stop(sprintf("Feature '%s' not found in assays or rowData.", feature))
  }

  .run_call <- function(val_vector) {
    .call_super_domains_on_vector(val_vector, feature, quantile_cutoff,
                                   log_transform = log_transform, verbose = verbose,
                                   n_bootstrap = n_bootstrap, seed = seed,
                                   method = method, min_quality = min_quality,
                                   tie_policy = tie_policy)
  }

  if (mode == "global_consensus") {
    # ---- Global (persistent) consensus across condition groups (P1-9) -------
    # Stratify replicate support by group so that a domain called Super in
    # 2/2 Control replicates but 0/2 Treatment replicates is NOT silently
    # called "Global Super" via a pooled 2/4 support fraction.
    #
    # P0-B: the global class must be decided by each group's GROUP_TYPE (does
    # the group pass the selected support rule?), NOT by support >= 1. Using
    # support >= 1 silently downgrades majority (2/3) groups to "all" and
    # makes support_rule meaningless.
    has_group <- group_var %in% colnames(meta)
    groups <- if (has_group) unique(meta[[group_var]]) else "ALL"

    per_sample_types <- lapply(colnames(mat), function(s) {
      .run_call(setNames(mat[, s], rownames(se)))$Domain_Type
    })
    names(per_sample_types) <- colnames(mat)
    type_mat <- do.call(cbind, per_sample_types)

    support_vec <- rep(NA_real_, nrow(se))
    min_support_vec <- rep(NA_real_, nrow(se))
    n_valid_vec <- rep(0L, nrow(se))
    n_super_reps <- rowSums(type_mat == paste0(feature, "_Super_Element"), na.rm = TRUE)
    group_type_mat <- matrix(NA_character_, nrow = nrow(se), ncol = length(groups))
    colnames(group_type_mat) <- groups

    for (g in groups) {
      idx <- if (has_group) which(meta[[group_var]] == g) else seq_len(ncol(mat))
      if (length(idx) == 0) next
      gm <- type_mat[, idx, drop = FALSE]
      g_res <- .replicate_support_call(gm, feature, support_rule,
                                       min_replicate_support, min_valid_replicates)
      group_type_mat[, g] <- g_res$group_type
      support_vec <- ifelse(is.na(support_vec), g_res$support,
                            pmin(support_vec, g_res$support, na.rm = TRUE))
      n_valid_vec <- n_valid_vec + g_res$n_valid
      rowData(se)[[paste0(feature, "_Support__", g)]] <- g_res$support
    }

    # Global class from per-group group_type (P0-B):
    #   any group NA (insufficient valid / no-call)  -> Uncertain
    #   all groups Super                             -> Super
    #   otherwise                                    -> Typical
    super_label <- paste0(feature, "_Super_Element")
    any_uncertain <- apply(group_type_mat, 1, function(r) any(is.na(r)))
    all_super <- apply(group_type_mat, 1, function(r) all(!is.na(r) & r == super_label))
    consensus_type <- rep(paste0(feature, "_Typical"), nrow(se))
    consensus_type[all_super] <- super_label
    consensus_type[any_uncertain] <- NA_character_

    rowData(se)[[paste0(feature, "_Domain_Type")]] <- consensus_type
    rowData(se)[[paste0(feature, "_Replicate_Support")]] <- support_vec
    rowData(se)[[paste0(feature, "_N_Super_Replicates")]] <- n_super_reps

    # Display/rank columns come from the group-MEAN signal (visualization only).
    # They are NOT the classification engine; the global consensus decision is
    # the per-sample -> replicate-support -> group-type chain above (review
    # §6). Keep them clearly separated from the decision provenance.
    mean_call <- .run_call(setNames(rowMeans(mat, na.rm = TRUE), rownames(se)))
    rowData(se)[[paste0(feature, "_Rank")]] <- mean_call$Rank
    rowData(se)[[paste0(feature, "_Value")]] <- mean_call$Value_Used

    n_super <- sum(consensus_type == super_label, na.rm = TRUE)
    if (verbose) message(sprintf("Global consensus: %d super (rule '%s' per group), %d domains.",
                                  n_super, support_rule, nrow(se)))
    # P1 (review §6): provenance must reflect the DECISION chain. Global
    # consensus has no single cutoff: each replicate is cut on its own ranked
    # distribution, support is aggregated per group. Store the per-sample
    # cutoffs/status/quality PROFILES for audit, the per-group support/type
    # matrices, and the final n_super from the real classification. The
    # mean-call (display) cutoff is kept under an explicit "display_rank"
    # key so it can never be mistaken for a classification cutoff.
    per_sample_prov <- lapply(colnames(mat), function(s) {
      r <- .run_call(setNames(mat[, s], rownames(se)))
      list(replicate = s, cutoff = r$cutoff_value,
           cutoff_stability_interval = r$cutoff_stability_interval,
           bootstrap_success_rate = r$bootstrap_success_rate,
           quality_score = r$quality_score,
           call_status = r$call_status,
           inflection_method = r$inflection_method)
    })
    names(per_sample_prov) <- colnames(mat)
    provenance <- list(
      mode = mode,
      decision_chain = "per-sample cut -> replicate support -> per-group type -> global consensus (all groups Super)",
      per_sample_evidence = per_sample_prov,
      group_support = if (has_group) {
        lapply(groups, function(g) rowData(se)[[paste0(feature, "_Support__", g)]])
      } else NULL,
      group_type = group_type_mat,
      support_rule = support_rule,
      min_replicate_support = min_replicate_support,
      min_valid_replicates = min_valid_replicates,
      n_super = n_super,
      n_total = nrow(se),
      cutoff_value = NA_real_,
      quality_score = NA_real_,
      call_status = "called-by-support",
      effective_log_transform = mean_call$effective_log_transform,
      display_rank = list(
        note = "group-MEAN-signal rank/display cutoff; NOT the classification cutoff",
        cutoff = mean_call$cutoff_value,
        quality_score = mean_call$quality_score,
        call_status = mean_call$call_status,
        n_super = mean_call$n_super),
      note = paste(
        "Global consensus: no single cutoff ranks the consensus; each",
        "replicate is cut on its own ranked distribution and group types are",
        "aggregated by the support rule. cutoff/quality above are per-sample",
        "profiles; the mean-signal cutoff under display_rank is for display only.")
    )

  } else if (mode == "per_group") {
    # ---- Replicate-aware per-group calling (P0-7) ---------------------------
    # Call each replicate, compute per-group replicate support via the
    # support rule, and report the group call. The group-mean rank is
    # retained for display/visualization.
    if (!group_var %in% colnames(meta)) stop("group_var not found in colData.")
    groups <- unique(meta[[group_var]])
    n_super <- NA_integer_
    last_prov <- NULL
    group_prov <- list()   # per-group replicate provenance (P1-7)
    for (g in groups) {
      idx <- which(meta[[group_var]] == g)
      if (length(idx) == 0) next
      g_calls <- lapply(idx, function(i) {
        .run_call(setNames(mat[, i], rownames(se)))
      })
      g_type <- do.call(cbind, lapply(g_calls, function(x) x$Domain_Type))
      g_res <- .replicate_support_call(g_type, feature, support_rule,
                                       min_replicate_support, min_valid_replicates)

      g_mean <- rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
      mean_call <- .run_call(setNames(g_mean, rownames(se)))

      rowData(se)[[paste0(feature, "_Call__", g)]] <- g_res$group_type
      rowData(se)[[paste0(feature, "_Support__", g)]] <- g_res$support
      rowData(se)[[paste0(feature, "_Rank__", g)]] <- mean_call$Rank
      rowData(se)[[paste0(feature, "_Value__", g)]] <- mean_call$Value_Used
      n_g_super <- sum(g_res$group_type == paste0(feature, "_Super_Element"), na.rm = TRUE)
      if (verbose) message("  ", g, ": ", n_g_super, " super (rule '", support_rule, "')")
      last_prov <- mean_call

      # Per-replicate provenance (P1-7): cutoff/quality/status for each
      # replicate in this group.
      rep_prov <- lapply(idx, function(i) {
        r <- g_calls[[match(i, idx)]]
        list(replicate = colnames(mat)[i],
             cutoff = r$cutoff_value,
             cutoff_stability_interval = r$cutoff_stability_interval,
             bootstrap_success_rate = r$bootstrap_success_rate,
             quality_score = r$quality_score,
             call_status = r$call_status,
             inflection_method = r$inflection_method)
      })
      group_prov[[g]] <- list(
        replicate_calls = rep_prov,
        # full domain x replicate call matrix for auditability (design
        # 2026-08-10: group calls are replicate-support aggregates; this matrix
        # exposes the per-domain per-replicate evidence).
        replicate_call_matrix = {
          m <- g_type
          colnames(m) <- colnames(mat)[idx]
          rownames(m) <- rownames(se)
          m
        },
        support_rule = support_rule,
        min_valid_replicates = min_valid_replicates,
        n_super = sum(g_res$group_type == paste0(feature, "_Super_Element"), na.rm = TRUE)
      )
    }
    # There is no single global cutoff in per_group mode: each replicate is cut
    # on its own ranked distribution, and per-replicate cutoffs / quality live
    # in groups[[g]]$replicate_calls. Do not expose the last group's group-mean
    # cutoff as if it were THE cutoff (P1-?): null out the ambiguous fields.
    provenance <- list(
      cutoff_value = NA_real_,
      cutoff_stability_interval = NULL,
      bootstrap_success_rate = NULL,
      quality_score = NA_real_,
      call_status = NA_character_,
      inflection_method = NA_character_,
      n_total = nrow(se),
      n_super = NA_integer_,
      effective_log_transform =
        if (!is.null(last_prov$effective_log_transform))
          last_prov$effective_log_transform else log_transform,
      note = "per_group: no single global cutoff; per-replicate cutoffs and quality are stored in metadata(se)$superdomain_calls[[feature]]$groups[[group]]$replicate_calls"
    )

  } else {
    # per_sample: call each sample independently; store the per-sample
    # cutoff/quality in provenance (symmetric with Breadth) so that
    # plot_hockey_stick(se, feature, group = <sample>) can draw the per-sample
    # cutoff line (freeze review 2026-08-11).
    rep_prov <- list()
    for (s in colnames(mat)) {
      call_res <- .run_call(setNames(mat[, s], rownames(se)))
      rowData(se)[[paste0(feature, "_Call__", s)]] <- call_res$Domain_Type
      rowData(se)[[paste0(feature, "_Rank__", s)]] <- call_res$Rank
      rep_prov[[s]] <- list(
        replicate = s,
        cutoff = call_res$cutoff_value,
        cutoff_stability_interval = call_res$cutoff_stability_interval,
        bootstrap_success_rate = call_res$bootstrap_success_rate,
        quality_score = call_res$quality_score,
        call_status = call_res$call_status,
        inflection_method = call_res$inflection_method,
        n_super = call_res$n_super)
    }
    first_call <- .run_call(setNames(mat[, 1], rownames(se)))
    n_super <- sum(first_call$Domain_Type == paste0(feature, "_Super_Element"), na.rm = TRUE)
    provenance <- first_call
  }

  # ---- Store call provenance (P0-5, P1-7) -----------------------------------
  if (is.null(S4Vectors::metadata(se)$superdomain_calls))
    S4Vectors::metadata(se)$superdomain_calls <- list()
  S4Vectors::metadata(se)$superdomain_calls[[feature]] <- list(
    feature = feature,
    mode = mode,
    group_var = if (mode %in% c("per_group", "global_consensus")) group_var else NULL,
    method = method,
    log_transform_requested = log_transform,
    log_transform_used = if (!is.null(provenance$effective_log_transform))
      provenance$effective_log_transform else log_transform,
    quantile_cutoff = quantile_cutoff,
    cutoff = provenance$cutoff_value,
    cutoff_stability_interval = provenance$cutoff_stability_interval,
    bootstrap_success_rate = provenance$bootstrap_success_rate,
    quality_score = provenance$quality_score,
    call_status = provenance$call_status,
    inflection_method = provenance$inflection_method,
    n_total = provenance$n_total,
    n_super = if (!is.null(provenance$n_super)) provenance$n_super else n_super,
    min_replicate_support = min_replicate_support,
    support_rule = support_rule,
    min_valid_replicates = min_valid_replicates,
    min_quality = min_quality,
    tie_policy = tie_policy,
    seed = seed,
    note = provenance$note,
    groups = if (exists("group_prov", inherits = FALSE)) group_prov else NULL,
    replicates = if (exists("rep_prov", inherits = FALSE)) rep_prov else NULL,
    # reviewer §6: global_consensus provenance must expose the DECISION chain
    # (per-sample evidence, per-group support/types, display cutoffs kept
    # separate). Only populated in global_consensus mode.
    decision_chain = provenance$decision_chain,
    per_sample_evidence = provenance$per_sample_evidence,
    group_support = provenance$group_support,
    group_type = provenance$group_type,
    display_rank = provenance$display_rank
  )

  se
}


# ---- Breadth-Super: peak-level native PeakWidth calling (v1.0 design) -------
#
# Implements the design decision (2026-08-09): Breadth-Super must NOT be ranked
# on a domain x sample assay. Instead:
#   1. each replicate's GENOME-WIDE eligible native peaks are filtered by fixed
#      QC (valid width, allowed chromosomes);
#   2. the replicate-specific native PeakWidth distribution is ranked and cut by
#      an elbow/inflection -> PeakBroadCall (peak-level label, independent of
#      any consensus construction);
#   3. native peaks are mapped to the shared domains by UNIQUE assignment (max
#      overlap + min_peak_overlap_fraction), so a broad peak cannot double-count
#      into neighbouring domains;
#   4. per-replicate domain broad evidence is aggregated across replicates
#      (support rule) to produce the group-level / global Breadth-Super.
.call_breadth_super_domains <- function(se, mode, group_var, method,
                                        quantile_cutoff, min_quality,
                                        tie_policy, n_bootstrap, seed,
                                        support_rule, min_replicate_support,
                                        min_valid_replicates,
                                        min_peak_overlap_fraction,
                                        valid_chroms, verbose) {
  # ---- 0. native peaks must be available -------------------------------
  np <- S4Vectors::metadata(se)$native_peaks
  if (is.null(np)) {
    stop("No native peak calls found. Breadth-Super requires peak_path in ",
         "build_portrait_matrix(). Re-run build_portrait_matrix with per-sample ",
         "peak files, or call a signal feature (e.g. 'Intensity') instead.")
  }
  meta <- as.data.frame(colData(se))
  feature <- "Breadth"
  super_label <- paste0(feature, "_Super_Element")
  typ_label <- paste0(feature, "_Typical")
  domains <- rowRanges(se)

  # ---- 1. replicate-specific peak-level broad calling ---------------------
  # Each replicate's cutoff is estimated on its OWN genome-wide eligible native
  # PeakWidth distribution, decoupled from the shared-domain universe. Two
  # mutually exclusive ways to set the cutoff (never auto-falling back between
  # them):
  #   * quantile_cutoff given (user opt-in): the (1 - quantile_cutoff)
  #     quantile of the eligible PeakWidth distribution -> top-fraction rule.
  #   * quantile_cutoff NULL (default): data-driven elbow/tangent inflection.
  #
  # A replicate whose inflection is NOT reliable (no_call) provides NO
  # broad/typical evidence at all: it must not be coerced to "all Typical"
  # (P0-1). The quantile path has no quality gate (the user owns the QC
  # choice), so it always yields a call (0 Broad for a constant-width
  # distribution). rep_call_valid lets the mapping step skip no-call
  # replicates -> Uncertain.
  peak_tables <- list()
  rep_broad <- vector("list", ncol(se))   # per-replicate logical vector over
                                          # that replicate's eligible peaks
  names(rep_broad) <- colnames(se)
  rep_call_valid <- stats::setNames(rep(FALSE, ncol(se)), colnames(se))
  prov_replicates <- list()
  for (s in colnames(se)) {
    peaks <- np[[s]]
    if (is.null(peaks) || length(peaks) == 0) {
      # Sample has no peak calls: cannot contribute broad evidence.
      rep_broad[[s]] <- rep(NA, 0L)
      prov_replicates[[s]] <- list(
        replicate = s, n_eligible = 0L, call_status = "no_call",
        reason = "no native peaks")
      next
    }
    eligible <- .eligible_native_peaks(peaks, valid_chroms)
    w <- GenomicRanges::width(eligible)
    if (length(w) == 0L) {
      # No eligible native peaks after (optional) valid_chroms filtering:
      # neither the quantile nor the elbow path can estimate a cutoff, so the
      # replicate provides no broad evidence -> Uncertain (P0-1). This must be
      # handled BEFORE quantile()/find_hockey_inflection() to avoid a runtime
      # error on an empty width vector.
      rep_broad[[s]] <- rep(NA, 0L)
      rep_call_valid[s] <- FALSE
      prov_replicates[[s]] <- list(
        replicate = s, n_eligible = 0L, call_status = "no_call",
        reason = "no eligible native peaks after filtering")
      next
    }
    if (!is.null(quantile_cutoff)) {
      # Explicit top-fraction rule, symmetric with Intensity: the cutoff is the
      # quantile_cutoff quantile of this replicate's genome-wide eligible
      # native PeakWidth distribution (bp), so the top
      # 100*(1-quantile_cutoff)% widest peaks are Broad. The user owns the QC
      # choice; no-call semantics of the data-driven path are untouched.
      cutoff_value <- as.numeric(stats::quantile(
        w, probs = quantile_cutoff, names = FALSE))
      method_label <- sprintf("quantile_%.2f", quantile_cutoff)
      quality_score <- NA_real_
    } else {
      inflect <- find_hockey_inflection(
        w, method = method, min_quality = min_quality, verbose = FALSE)
      if (inflect$call_status != "called") {
        # no reliable elbow -> no evidence for this replicate (P0-1).
        rep_broad[[s]] <- rep(NA, length(eligible))
        rep_call_valid[s] <- FALSE
        prov_replicates[[s]] <- list(
          replicate = s, n_eligible = length(eligible),
          cutoff = NA_real_, quality_score = inflect$quality_score,
          call_status = inflect$call_status, reason = inflect$reason)
        if (verbose) {
          message(sprintf("  [Breadth] %s: no reliable inflection -> no evidence (Uncertain)",
                          s))
        }
        next
      }
      cutoff_value <- inflect$cutoff_value
      method_label <- inflect$method
      quality_score <- inflect$quality_score
    }
    rep_call_valid[s] <- TRUE
    if (tie_policy == "strict") {
      is_broad <- w > cutoff_value
    } else {
      is_broad <- w >= cutoff_value
    }
    rep_broad[[s]] <- is_broad

    # Per-replicate within-set cutoff stability (symmetric with Intensity):
    # resample the replicate's eligible native PeakWidth distribution with
    # replacement and re-estimate the inflection cutoff, giving a per-replicate
    # stability interval + success rate. This is a resampling-stability
    # diagnostic (not a classical CI) and, like Intensity, its information
    # content shrinks as the ranked set grows; it only applies to the
    # data-driven inflection path (quantile_cutoff is deterministic). The
    # dominant uncertainty for Breadth is the upstream native-peak SET
    # composition, which this within-set bootstrap does not capture — see
    # validation/02 (peak-set perturbation).
    cutoff_stability_interval <- NULL
    bootstrap_success_rate <- NULL
    if (!is.null(n_bootstrap) && n_bootstrap > 0 && is.null(quantile_cutoff)) {
      boot_cutoffs <- .with_opt_seed(seed, {
        cutoffs <- numeric(n_bootstrap)
        for (b in seq_len(n_bootstrap)) {
          wr <- w[sample(length(w), replace = TRUE)]
          br <- tryCatch(
            find_hockey_inflection(wr, method = method,
                                   min_quality = min_quality, verbose = FALSE),
            error = function(e) list(call_status = "no_call",
                                     cutoff_value = NA_real_))
          # Only a "called" resample is a bootstrap success: no_call resamples
          # (quality below min_quality) can still carry a finite cutoff_value,
          # which must NOT be counted as a stable cutoff (P1, review 14 §11).
          cutoffs[b] <- if (identical(br$call_status, "called")) {
            br$cutoff_value
          } else {
            NA_real_
          }
        }
        cutoffs
      })
      finite_idx <- is.finite(boot_cutoffs)
      cutoff_stability_interval <- if (any(finite_idx)) {
        stats::quantile(boot_cutoffs[finite_idx], c(0.025, 0.975), na.rm = TRUE)
      } else NULL
      bootstrap_success_rate <- mean(finite_idx)
    }

    prov_replicates[[s]] <- list(
      replicate = s, n_eligible = length(eligible),
      cutoff = cutoff_value,
      cutoff_stability_interval = cutoff_stability_interval,
      bootstrap_success_rate = bootstrap_success_rate,
      quality_score = quality_score,
      call_status = "called",
      n_broad = sum(is_broad),
      inflection_method = method_label)
    if (verbose) {
      message(sprintf("  [Breadth] %s: %d broad / %d eligible (cutoff %s %.0f bp)",
                      s, sum(is_broad), length(eligible),
                      if (is.null(quantile_cutoff))
                        ">= " else sprintf("= top %s%%, >= ",
                                            format(100 * (1 - quantile_cutoff), digits = 2)),
                      cutoff_value))
    }
    # peak-level table (design: PeakWidth / PeakBroadCall / PeakRank)
    tab <- data.frame(
      NativePeakID = paste(s, seq_along(eligible), sep = "."),
      SampleID = s,
      seqnames = as.character(GenomicRanges::seqnames(eligible)),
      start = GenomicRanges::start(eligible),
      end = GenomicRanges::end(eligible),
      PeakWidth = as.numeric(w),
      PeakRank = rank(-w, ties.method = "min"),
      PeakBroadCall = ifelse(is_broad, "Broad", "Typical"),
      PeakWidthCutoff = cutoff_value,
      stringsAsFactors = FALSE)
    peak_tables[[s]] <- tab
  }

  # ---- 2. unique peak -> domain mapping ------------------------------------
  # A native peak is assigned to at most ONE shared domain: the overlapping
  # domain with the maximum overlap bp, provided the peak-overlap fraction
  # >= min_peak_overlap_fraction. Unmapped / Ambiguous peaks NEVER provide
  # broad OR typical evidence (P0-2): only successfully uniquely-mapped peaks
  # contribute, so an edge-overlap Unmapped peak cannot silently turn a domain
  # into Typical.
  domain_evidence <- matrix(NA_character_, nrow = nrow(se), ncol = ncol(se),
                            dimnames = list(rownames(se), colnames(se)))
  mapping_prov <- list()
  for (s in colnames(se)) {
    peaks <- np[[s]]
    if (!isTRUE(rep_call_valid[s])) {
      # no native peaks, or no reliable elbow -> no evidence -> Uncertain
      domain_evidence[, s] <- NA_character_
      next
    }
    eligible <- .eligible_native_peaks(peaks, valid_chroms)
    w <- GenomicRanges::width(eligible)
    is_broad <- rep_broad[[s]]

    # All overlapping domain x peak pairs (any overlap for geometry).
    hits <- GenomicRanges::findOverlaps(domains, eligible)
    if (length(hits) == 0) {
      domain_evidence[, s] <- NA_character_
      next
    }
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    ov <- GenomicRanges::pintersect(
      domains[qh], eligible[sh], ignore.strand = TRUE)
    ov_bp <- GenomicRanges::width(ov)
    peak_frac <- ov_bp / pmax(w[sh], 1)
    domain_frac <- ov_bp / pmax(GenomicRanges::width(domains)[qh], 1)

    # For each peak, the unique best domain (max overlap bp); ties -> Ambiguous.
    status <- rep("Unmapped", length(eligible))
    assign_domain <- rep(NA_integer_, length(eligible))
    best_ov_bp <- rep(NA_real_, length(eligible))
    best_peak_frac <- rep(NA_real_, length(eligible))
    best_domain_frac <- rep(NA_real_, length(eligible))
    for (pk in unique(sh)) {
      hits_pk <- which(sh == pk)
      best <- which.max(ov_bp[hits_pk])
      if (length(hits_pk) > 1) {
        # ties for the max overlap -> ambiguous (never double count)
        ties <- which(ov_bp[hits_pk] == ov_bp[hits_pk][best])
        if (length(ties) > 1) {
          status[pk] <- "Ambiguous"
          next
        }
      }
      d_best <- qh[hits_pk][best]
      if (peak_frac[hits_pk][best] >= min_peak_overlap_fraction) {
        status[pk] <- "Unique"
        assign_domain[pk] <- d_best
        best_ov_bp[pk] <- ov_bp[hits_pk][best]
        best_peak_frac[pk] <- peak_frac[hits_pk][best]
        best_domain_frac[pk] <- domain_frac[hits_pk][best]
      }
    }

    # Domain-level evidence for this replicate. ONLY uniquely-mapped peaks
    # contribute: broad peaks set Super, and typical evidence comes only from
    # uniquely-mapped non-broad peaks (P0-2). Unmapped / Ambiguous peaks never
    # create evidence.
    ev <- rep(NA_character_, nrow(se))
    mapped_idx <- which(!is.na(assign_domain))
    for (pk in mapped_idx[is_broad[mapped_idx]]) {
      ev[assign_domain[pk]] <- super_label
    }
    for (pk in mapped_idx[!is_broad[mapped_idx]]) {
      d <- assign_domain[pk]
      if (is.na(ev[d])) ev[d] <- typ_label   # broad evidence wins over typical
    }
    domain_evidence[, s] <- ev

    mapping_prov[[s]] <- data.frame(
      NativePeakID = paste(s, seq_along(eligible), sep = "."),
      SampleID = s,
      SharedDomainIndex = assign_domain,
      SharedDomainID = ifelse(is.na(assign_domain), NA_character_,
                              rownames(se)[assign_domain]),
      PeakWidth = as.numeric(w),
      PeakBroadCall = ifelse(is_broad, "Broad", "Typical"),
      OverlapBp = best_ov_bp,
      PeakOverlapFraction = best_peak_frac,
      DomainOverlapFraction = best_domain_frac,
      MappingStatus = status,
      stringsAsFactors = FALSE)
  }

  # ---- 3. replicate-aware aggregation ---------------------------------------
  if (mode == "global_consensus") {
    has_group <- group_var %in% colnames(meta)
    groups <- if (has_group) unique(meta[[group_var]]) else "ALL"
    group_type_mat <- matrix(NA_character_, nrow = nrow(se), ncol = length(groups),
                             dimnames = list(rownames(se), groups))
    support_vec <- rep(NA_real_, nrow(se))
    n_valid_vec <- rep(0L, nrow(se))
    for (g in groups) {
      idx <- if (has_group) which(meta[[group_var]] == g) else seq_len(ncol(se))
      g_res <- .replicate_support_call(domain_evidence[, idx, drop = FALSE],
                                       feature, support_rule = support_rule,
                                       min_replicate_support = min_replicate_support,
                                       min_valid_replicates = min_valid_replicates)
      group_type_mat[, g] <- g_res$group_type
      support_vec <- ifelse(is.na(support_vec), g_res$support,
                            pmin(support_vec, g_res$support, na.rm = TRUE))
      n_valid_vec <- n_valid_vec + g_res$n_valid
    }
    any_uncertain <- apply(group_type_mat, 1, function(r) any(is.na(r)))
    all_super <- apply(group_type_mat, 1, function(r) all(!is.na(r) & r == super_label))
    consensus_type <- rep(typ_label, nrow(se))
    consensus_type[all_super] <- super_label
    consensus_type[any_uncertain] <- NA_character_
    rowData(se)[["Breadth_Domain_Type"]] <- consensus_type
    rowData(se)[["Breadth_Replicate_Support"]] <- support_vec
    rowData(se)[["Breadth_N_Broad_Replicates"]] <-
      rowSums(domain_evidence == super_label, na.rm = TRUE)
    if (verbose) message(sprintf("Global Breadth consensus: %d Breadth-Super.",
                                 sum(consensus_type == super_label, na.rm = TRUE)))

  } else if (mode == "per_group") {
    if (!group_var %in% colnames(meta)) stop("group_var not found in colData.")
    groups <- unique(meta[[group_var]])
    breadth_group_prov <- list()
    for (g in groups) {
      idx <- which(meta[[group_var]] == g)
      g_res <- .replicate_support_call(domain_evidence[, idx, drop = FALSE],
                                       feature, support_rule = support_rule,
                                       min_replicate_support = min_replicate_support,
                                       min_valid_replicates = min_valid_replicates)
      rowData(se)[[paste0("Breadth_Call__", g)]] <- g_res$group_type
      rowData(se)[[paste0("Breadth_Support__", g)]] <- g_res$support
      rowData(se)[[paste0("Breadth_N_Broad_Replicates__", g)]] <-
        rowSums(domain_evidence[, idx, drop = FALSE] == super_label, na.rm = TRUE)
      # expose the domain x replicate Broad evidence matrix for
      # get_replicate_calls() (Breadth replicate audit trail must
      # be symmetric with Intensity, design 2026-08-10).
      breadth_group_prov[[g]] <- list(
        replicate_call_matrix = {
          m <- domain_evidence[, idx, drop = FALSE]
          colnames(m) <- colnames(se)[idx]
          rownames(m) <- rownames(se)
          m
        },
        support_rule = support_rule,
        min_valid_replicates = min_valid_replicates,
        n_super = sum(g_res$group_type == super_label, na.rm = TRUE)
      )
      if (verbose) {
        message(sprintf("  [Breadth] %s: %d Breadth-Super (rule '%s')",
                        g, sum(g_res$group_type == super_label, na.rm = TRUE),
                        support_rule))
      }
    }

  } else { # per_sample
    for (s in colnames(se)) {
      rowData(se)[[paste0("Breadth_Call__", s)]] <- domain_evidence[, s]
    }
  }

  # ---- 4. provenance -------------------------------------------------------
  if (is.null(S4Vectors::metadata(se)$superdomain_calls))
    S4Vectors::metadata(se)$superdomain_calls <- list()
  S4Vectors::metadata(se)$superdomain_calls[["Breadth"]] <- list(
    feature = "Breadth",
    mode = mode,
    group_var = if (mode %in% c("per_group", "global_consensus")) group_var else NULL,
    method = method,
    quantile_cutoff = quantile_cutoff,
    min_quality = min_quality,
    tie_policy = tie_policy,
    support_rule = support_rule,
    min_replicate_support = min_replicate_support,
    min_valid_replicates = min_valid_replicates,
    n_bootstrap = n_bootstrap,
    seed = seed,
    min_peak_overlap_fraction = min_peak_overlap_fraction,
    valid_chroms = valid_chroms,
    replicates = prov_replicates,
    groups = if (exists("breadth_group_prov", inherits = FALSE)) breadth_group_prov else NULL,
    n_total = nrow(se),
    final_n_super = if (mode == "global_consensus") {
      sum(consensus_type == super_label, na.rm = TRUE)
    } else NULL,
    calling_paradigm = "peak-level native PeakWidth, unique mapping, replicate aggregation")
  if (length(peak_tables) > 0) {
    S4Vectors::metadata(se)$breadth_peak_calls <- do.call(rbind, peak_tables)
    S4Vectors::metadata(se)$breadth_peak_mapping <- do.call(rbind, mapping_prov)
  }

  se
}


# Basic fixed QC for the genome-wide eligible native peak set used by
# Breadth-Super calling (design: cutoff estimated on the eligible set, NOT on
# domain-restricted peaks). Filters: valid width > 0 and, when provided,
# allowed chromosomes.
.eligible_native_peaks <- function(peaks, valid_chroms = NULL) {
  ok <- GenomicRanges::width(peaks) > 0 &
    !is.na(GenomicRanges::seqnames(peaks))
  if (!is.null(valid_chroms)) {
    ok <- ok & as.character(GenomicRanges::seqnames(peaks)) %in% valid_chroms
  }
  peaks[ok]
}


