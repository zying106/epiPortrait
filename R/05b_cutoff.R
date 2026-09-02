
#' Find the Hockey-Stick Inflection Point
#'
#' Locates the inflection point of a ranked (ascending) value curve used to
#' separate super-domains from typical domains. Two heuristics are supported:
#'
#' \itemize{
#'   \item \code{"elbow"} (default): the point with the maximum perpendicular
#'         distance to the line connecting the first and last ranked points.
#'         Classic knee/elbow detection; cheap and robust when signals are
#'         log-transformed.
#'   \item \code{"tangent"}: a ROSE-inspired tangent-optimization inflection
#'         (Whyte et al., 2013, Cell 153:307-319). Slides a line whose slope is
#'         \code{(max - min)/n} and minimizes the number of points below it.
#'         Note that this is a geometric variant of ROSE's
#'         \code{calculate_cutoff()} step — not a bit-exact reproduction — and
#'         implements only an inflection step, not the full ROSE pipeline
#'         (stitching, input subtraction). epiPortrait does not log-transform
#'         the input unless \code{log_transform} is requested.
#' }
#'
#' A method-agnostic \code{quality_score} (max perpendicular distance relative
#' to the value range) is returned. When it falls below \code{min_quality},
#' \code{call_status} is \code{"no_call"} and the inflection is not reliable.
#'
#' @param x Numeric vector of feature values (e.g. per-domain intensity).
#' @param method Character. \code{"elbow"} (default) or \code{"tangent"}.
#' @param min_quality Numeric. Minimum \code{quality_score} for a reliable call.
#'   The score is dimensionless (0-1, computed on \eqn{[0,1]}-normalized coordinates)
#'   and measures curve prominence: ~0 for a linear curve, ~0.3-0.7 for real
#'   hockey-stick profiles. Default 0.1 rejects near-linear curves while
#'   retaining genuine inflection points.
#' @param verbose Logical. Print warnings on low-quality calls (default FALSE).
#'
#' @return A list with elements:
#'   \item{inflection_idx}{Index of the inflection point in the sorted vector.}
#'   \item{cutoff_value}{Feature value at the inflection point.}
#'   \item{quality_score}{Right-tail curve-prominence quality (0-1). A
#'   left-loaded / inverse hockey-stick curve has quality 0 and is not called.}
#'   \item{method}{The method used.}
#'   \item{call_status}{\code{"called"} or \code{"no_call"}.}
#'   \item{reason}{Reason for \code{"no_call"}, if any.}
#' @examples
#' x <- c(rep(1, 50), seq(1, 100, length.out = 50))
#' find_hockey_inflection(x, method = "elbow")
#' find_hockey_inflection(x, method = "tangent")
#' @export
find_hockey_inflection <- function(x, method = c("elbow", "tangent"),
                                   min_quality = 0.1, verbose = FALSE) {
  method <- match.arg(method)
  if (length(min_quality) != 1L || !is.numeric(min_quality) ||
      !is.finite(min_quality) || min_quality < 0 || min_quality > 1) {
    stop("min_quality must be a finite number in [0, 1].", call. = FALSE)
  }
  x <- x[is.finite(x)]
  n <- length(x)

  no_call <- function(reason) {
    if (verbose) warning(reason, call. = FALSE)
    list(inflection_idx = NA_integer_, cutoff_value = NA_real_,
         quality_score = NA_real_, method = method, call_status = "no_call",
         reason = reason)
  }

  if (n < 3) return(no_call("Fewer than 3 finite values. No reliable inflection."))
  if (diff(range(x)) <= .Machine$double.eps)
    return(no_call("Constant feature values. No reliable inflection."))

  xs <- sort(x)
  if (method == "tangent") {
    xs_clamp <- xs
    xs_clamp[xs_clamp < 0] <- 0   # ROSE clamps negative (control-subtracted) signals to 0
    slope <- (max(xs_clamp) - min(xs_clamp)) / n
    if (slope <= .Machine$double.eps)
      return(no_call("Zero slope in tangent fit. No reliable inflection."))
    numPts_below <- function(z) {
      z <- pmin(pmax(z, 1), n)
      yPt <- xs_clamp[z]           # R truncates non-integer index like ROSE's myVector[x]
      b <- yPt - (slope * z)
      sum(xs_clamp <= (seq_len(n) * slope + b))
    }
    opt <- stats::optimize(numPts_below, lower = 1, upper = n)
    inflection_idx <- as.integer(floor(opt$minimum))
    inflection_idx <- min(max(inflection_idx, 1L), n)
  } else {
    # elbow: SIGNED perpendicular distance below the first-last line.
    # Both coordinates are normalized to [0,1] (P1-7) so the geometry is
    # dimensionless and independent of the feature scale / sample size.
    # The first point maps to (0,0) and the last to (1,1); a perfectly linear
    # curve then coincides with the diagonal y = x (distance 0).
    xr <- (seq_len(n) - 1) / max(n - 1, 1)              # rank in [0,1]
    yr <- (xs - xs[1]) / max(xs[n] - xs[1], 1e-10)      # signal in [0,1]
    # Line through the two endpoints; on the normalized square the endpoints
    # are (0,0) and (1,1) so slope = 1 exactly.
    slope <- (yr[n] - yr[1]) / max(xr[n] - xr[1], 1e-10)
    denom <- sqrt(slope^2 + 1)
    # A valid ascending hockey stick has its middle ranks BELOW the endpoint
    # line (xr - yr > 0) and then rises in the right tail. Do not use abs():
    # it makes the mirror-image, left-loaded curve look equally prominent and
    # can classify almost the entire universe as Super.
    perp <- (slope * xr - yr) / denom
    inflection_idx <- which.max(perp)
  }

  # Method-agnostic curve-prominence quality, computed on normalized coords so
  # min_quality is a dimensionless, scale-invariant threshold (P1-7).
  xr0 <- (seq_len(n) - 1) / max(n - 1, 1)
  yr0 <- (xs - xs[1]) / max(xs[n] - xs[1], 1e-10)
  slope0 <- (yr0[n] - yr0[1]) / max(xr0[n] - xr0[1], 1e-10)
  perp0 <- (slope0 * xr0 - yr0) / sqrt(slope0^2 + 1)
  # Endpoints are exactly zero. An inverse/left-loaded curve has no positive
  # right-tail prominence, hence quality 0 rather than a large false score.
  quality_score <- max(0, max(perp0))

  cutoff_value <- xs[inflection_idx]
  call_status <- "called"
  reason <- NA_character_
  if (is.finite(quality_score) && quality_score < min_quality) {
    call_status <- "no_call"
    reason <- sprintf("quality=%.4f < %.4f. No reliable inflection.", quality_score, min_quality)
    if (verbose) warning(reason, call. = FALSE)
  }

  list(inflection_idx = inflection_idx, cutoff_value = cutoff_value,
       quality_score = quality_score, method = method,
       call_status = call_status, reason = reason)
}


# ---- Internal: single-value-vector super-domain calling -------------------
.call_super_domains_on_vector <- function(val_vector, feature, quantile_cutoff,
                                           log_transform, verbose,
                                           n_bootstrap = NULL, seed = NULL,
                                           method = "elbow", min_quality = 0.1,
                                           tie_policy = "strict") {

  # Return a clean no-call record for degenerate input rather than stopping:
  # Bioconductor packages should not error on legitimate boundary data
  # (constant values, all-NA). The caller records no_call and the combined
  # taxonomy propagates it as "Uncertain".
  .no_call <- function(reason) {
    list(Domain_Type = rep(NA_character_, length(val_vector)),
         Rank = rep(NA_real_, length(val_vector)),
         Value_Used = rep(NA_real_, length(val_vector)),
         cutoff_value = NA_real_,
         inflection_idx = NA_integer_,
         inflection_method = NA_character_,
         quality_score = NA_real_,
         call_status = "no_call",
         reason = reason,
         n_total = length(val_vector),
         n_super = NA_integer_,
         cutoff_stability_interval = NULL,
         bootstrap_success_rate = NULL)
  }

  if (!any(is.finite(val_vector)) || all(is.na(val_vector)))
    return(.no_call("All values are NA or non-finite."))
  if (diff(range(val_vector, na.rm = TRUE)) <= .Machine$double.eps)
    return(.no_call("Constant feature values. No reliable inflection."))

  # Resolve log-transform ONCE (P0-2). The user may pass NULL (auto), TRUE,
  # or FALSE. We resolve to a concrete value here so that bootstrap recursion
  # uses the SAME effective scale as the main analysis. Passing NULL into the
  # recursion would re-resolve to the auto rule, which could differ from an
  # explicit user choice (e.g. log_transform = FALSE would silently become
  # TRUE for Intensity inside the bootstrap).
  effective_log_transform <- if (is.null(log_transform)) {
    feature %in% c("Intensity", "SignalDispersion")
  } else {
    isTRUE(log_transform)
  }

  # Keep the RAW (untransformed) values for bootstrap resampling. Sampling
  # must be done on the raw scale, then the same single transform applied
  # inside the recursive call — otherwise a log-transformed vector would be
  # log-transformed a second time (double-log), corrupting the cutoff CI.
  raw_vector <- val_vector
  if (effective_log_transform) val_vector <- log10(pmax(val_vector, 0) + 1)

  pid <- names(val_vector)
  if (is.null(pid)) pid <- as.character(seq_along(val_vector))
  rank_df <- data.frame(Peak_ID = pid,
                         Value = val_vector, stringsAsFactors = FALSE)
  rank_df <- rank_df[order(rank_df$Value, rank_df$Peak_ID), ]
  rank_df$Cumulative_Rank <- seq_len(nrow(rank_df))

  n <- nrow(rank_df)
  call_status <- "called"
  inflection_idx <- NA_integer_
  inflection_method <- NA_character_
  quality_score <- NA_real_
  cutoff_value <- NA_real_

  if (!is.null(quantile_cutoff)) {
    cutoff_value <- stats::quantile(rank_df$Value, probs = quantile_cutoff, na.rm = TRUE)
    inflection_idx <- which(rank_df$Value >= cutoff_value)[1]
    if (is.na(inflection_idx)) inflection_idx <- n
    inflection_method <- sprintf("quantile_%.2f", quantile_cutoff)
  } else {
    inflect <- find_hockey_inflection(rank_df$Value, method = method,
                                      min_quality = min_quality, verbose = verbose)
    inflection_idx <- inflect$inflection_idx
    quality_score <- inflect$quality_score
    call_status <- inflect$call_status
    inflection_method <- if (method == "tangent") "tangent" else "elbow_distance"
    if (call_status == "no_call" && verbose)
      message(sprintf("  [%s] %s", feature, inflect$reason))
  }

  if (call_status == "called") {
    cutoff_value <- rank_df$Value[inflection_idx]
    inflection_idx <- min(inflection_idx, n)
    # P0-D: tie policy. ROSE selects super-enhancers strictly ABOVE the
    # cutoff (> cutoff). "strict" (default) follows ROSE and avoids inflating
    # the super set when many discrete values (e.g. integer widths) tie exactly
    # at the cutoff. "inclusive" (>= cutoff) is the legacy behaviour.
    if (tie_policy == "strict") {
      rank_df$Domain_Type <- ifelse(rank_df$Value > cutoff_value,
                                    paste0(feature, "_Super_Element"),
                                    paste0(feature, "_Typical"))
    } else {
      rank_df$Domain_Type <- ifelse(rank_df$Value >= cutoff_value,
                                    paste0(feature, "_Super_Element"),
                                    paste0(feature, "_Typical"))
    }
    if (effective_log_transform) cutoff_value <- 10^cutoff_value - 1
  } else {
    rank_df$Domain_Type <- NA_character_
    cutoff_value <- NA_real_
  }

  map_idx <- match(names(val_vector), rank_df$Peak_ID)

  cutoff_stability_interval <- NULL
  bootstrap_success_rate <- NULL
  if (!is.null(n_bootstrap) && n_bootstrap > 0 && call_status == "called") {
    boot_cutoffs <- .with_opt_seed(seed, {
      cutoffs <- numeric(n_bootstrap)
      n_peaks <- length(raw_vector)
      for (b in seq_len(n_bootstrap)) {
        boot_idx <- sample(n_peaks, replace = TRUE)
        boot_vec <- raw_vector[boot_idx]
        names(boot_vec) <- seq_len(n_peaks)
        # Recurse from the RAW vector using the SAME effective transform scale
        # as the main analysis (P0-2): FALSE stays FALSE, TRUE stays TRUE, and
        # AUTO is resolved once here. Never re-resolve inside the recursion.
        boot_res <- tryCatch(
          .call_super_domains_on_vector(boot_vec, feature, quantile_cutoff,
                                        effective_log_transform,
                                        verbose = FALSE, method = method,
                                        min_quality = min_quality,
                                        tie_policy = tie_policy),
          error = function(e) list(cutoff_value = NA_real_)
        )
        cutoffs[b] <- boot_res$cutoff_value
      }
      cutoffs
    })
    # The bootstrap reflects resampling stability over the fixed candidate
    # universe, not a classical iid confidence interval; report it as a
    # stability interval and record how many resamples yielded a valid cutoff
    # (review #9).
    finite_idx <- is.finite(boot_cutoffs)
    cutoff_stability_interval <- if (any(finite_idx)) {
      stats::quantile(boot_cutoffs[finite_idx],
                      probs = c(0.025, 0.975), na.rm = TRUE)
    } else NULL
    bootstrap_success_rate <- mean(finite_idx)
  }

  list(
    Domain_Type = rank_df$Domain_Type[map_idx],
    Rank = rank_df$Cumulative_Rank[map_idx],
    Value_Used = rank_df$Value[map_idx],
    cutoff_value = cutoff_value,
    inflection_idx = inflection_idx,
    inflection_method = inflection_method,
    quality_score = quality_score,
    call_status = call_status,
    effective_log_transform = effective_log_transform,
    n_total = n,
    n_super = sum(grepl("_Super", rank_df$Domain_Type), na.rm = TRUE),
    cutoff_stability_interval = cutoff_stability_interval,
    bootstrap_success_rate = bootstrap_success_rate
  )
}


# Compute a replicate-support group call from a matrix of per-replicate
# domain types (rows = domains, cols = replicates). Implements P0-3:
#
#   * support_rule = "majority": need > n/2 super replicates
#     (floor(n/2)+1); n=2 needs 2, n=3 needs 2.
#   * support_rule = "all": every replicate must be super.
#   * support_rule = "fraction": fraction >= min_replicate_support.
#
#   * no-call (NA) replicates are NOT counted as super, and they are NOT
#     silently dropped from the denominator in a way that inflates support.
#     If the number of valid (non-NA) replicates is < min_valid_replicates,
#     the group call is NA (Uncertain), never Super.
#
# Returns list(group_type, support, n_valid, n_super).
.replicate_support_call <- function(type_mat, feature,
                                     support_rule = c("majority", "all", "fraction"),
                                     min_replicate_support = 0.5,
                                     min_valid_replicates = NULL) {
  support_rule <- match.arg(support_rule)
  n_super <- rowSums(type_mat == paste0(feature, "_Super_Element"), na.rm = TRUE)
  n_valid <- rowSums(!is.na(type_mat))
  n_reps <- ncol(type_mat)
  if (n_reps < 1L) stop("type_mat must contain at least one replicate column.")
  if (!is.null(min_valid_replicates) &&
      (length(min_valid_replicates) != 1L ||
       !is.numeric(min_valid_replicates) ||
       !is.finite(min_valid_replicates) ||
       min_valid_replicates < 1 ||
       min_valid_replicates != floor(min_valid_replicates))) {
    stop("min_valid_replicates must be NULL or a positive integer.")
  }

  required <- switch(support_rule,
    majority = floor(n_reps / 2) + 1L,
    all = n_reps,
    fraction = ceiling(min_replicate_support * n_reps)
  )

  # E: if min_valid_replicates is NULL, require at least the number of valid
  # replicates the support rule needs. This prevents "Super + no_call" in an
  # n=2 group from being called Typical (it is actually Uncertain), because a
  # single valid replicate cannot establish 2-replicate reproducibility.
  if (is.null(min_valid_replicates)) min_valid_replicates <- required

  # Support as fraction of ALL replicates (not just valid ones): a no-call
  # replicate must not boost support (P0-3).
  support <- n_super / max(n_reps, 1)
  enough_valid <- n_valid >= min_valid_replicates

  group_type <- ifelse(enough_valid & n_super >= required,
                       paste0(feature, "_Super_Element"),
                       ifelse(n_valid == 0, NA_character_,
                              ifelse(enough_valid, paste0(feature, "_Typical"), NA_character_)))
  list(group_type = group_type, support = support,
       n_valid = n_valid, n_super = n_super)
}


# Classify per-replicate signal vectors against a FIXED cutoff, returning a
# matrix of domain types (rows = domains, cols = replicates). Used by
# compare_superdomains(cutoff_scope = "reference"/"pooled") to keep the
# transition replicate-aware (P0-C). NA values stay NA (never coerced to
# "Typical"), and the tie policy (strict ">" vs inclusive ">=") is inherited
# from the primary call (freeze review 2026-08-11).
.classify_vs_cutoff <- function(mat, cutoff, feature,
                                tie_policy = c("strict", "inclusive")) {
  tie_policy <- match.arg(tie_policy)
  super <- paste0(feature, "_Super_Element")
  typ <- paste0(feature, "_Typical")
  res <- apply(mat, 2, function(v) {
    hit <- if (tie_policy == "inclusive") v >= cutoff else v > cutoff
    ifelse(is.finite(v) & hit, super,
           ifelse(is.finite(v), typ, NA_character_))
  })
  if (is.null(dim(res))) res <- matrix(res, nrow = nrow(mat))
  res
}
