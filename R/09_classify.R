#' Classify Shape Shift Events into Biologically Interpretable Categories
#'
#' @description Translates the abstract 4D geometric changes detected by
#' \code{test_global_shape_shift} into named, interpretable categories:
#' \emph{Concentration}, \emph{Flattening}, \emph{Polarity Shift},
#' \emph{Global Gain}, \emph{Global Loss}, or \emph{Complex}. This bridges the
#' gap between multivariate statistics and biological narrative.
#'
#' @param se A \code{SummarizedExperiment} object from \code{build_portrait_matrix()}.
#' @param shift_res A \code{data.frame} returned by \code{test_global_shape_shift()}.
#' @param group_var Character. Column in \code{colData(se)} with the grouping variable.
#' @param target_group Character. The treatment group.
#' @param ref_group Character. The reference/control group.
#' @param fc_threshold Numeric. Log2 fold-change threshold for calling a
#'   dimension as changed (default: 0.5, i.e. ~1.4-fold).
#' @param skewness_threshold Numeric. Absolute difference threshold for
#'   Skewness shift (default: 0.3).
#'
#' @return The input \code{shift_res} data.frame augmented with columns:
#'   \code{Intensity_FC}, \code{Height_FC}, \code{Skewness_Delta},
#'   \code{Shape_Class}, and \code{Class_Label} (human-readable description).
#'
#' @import SummarizedExperiment
#' @importFrom stats setNames
#' @export
classify_shape_shift <- function(se, shift_res,
                                  group_var = "Condition",
                                  target_group, ref_group,
                                  fc_threshold = 0.5,
                                  skewness_threshold = 0.3) {

  stopifnot(is(se, "SummarizedExperiment"))
  stopifnot(is.data.frame(shift_res))
  stopifnot("Peak_ID" %in% colnames(shift_res))

  meta <- as.data.frame(colData(se))
  idx_target <- which(meta[[group_var]] == target_group)
  idx_ref    <- which(meta[[group_var]] == ref_group)

  if (length(idx_target) == 0) stop(sprintf("Group '%s' not found.", target_group))
  if (length(idx_ref) == 0)    stop(sprintf("Group '%s' not found.", ref_group))

  peak_ids <- shift_res$Peak_ID
  if (!all(peak_ids %in% rownames(se))) {
    missing <- setdiff(peak_ids, rownames(se))
    stop(sprintf("%d peak(s) in shift_res not found in se rownames.", length(missing)))
  }

  # Extract per-dimension group means -------------------------------------------------
  i_mat <- assay(se, "Intensity")
  h_mat <- assay(se, "Height")
  s_mat <- assay(se, "Skewness")

  intensity_fc <- vapply(peak_ids, function(pid) {
    tgt <- mean(i_mat[pid, idx_target], na.rm = TRUE)
    ref <- mean(i_mat[pid, idx_ref],    na.rm = TRUE)
    log2((tgt + 1) / (ref + 1))
  }, numeric(1))

  height_fc <- vapply(peak_ids, function(pid) {
    tgt <- mean(h_mat[pid, idx_target], na.rm = TRUE)
    ref <- mean(h_mat[pid, idx_ref],    na.rm = TRUE)
    log2((tgt + 1) / (ref + 1))
  }, numeric(1))

  skewness_delta <- vapply(peak_ids, function(pid) {
    mean(s_mat[pid, idx_target], na.rm = TRUE) -
    mean(s_mat[pid, idx_ref],    na.rm = TRUE)
  }, numeric(1))

  # Classification logic --------------------------------------------------------------
  classify_one <- function(ifc, hfc, sdelta) {
    i_up   <- ifc  >  fc_threshold
    i_down <- ifc  < -fc_threshold
    h_up   <- hfc  >  fc_threshold
    h_down <- hfc  < -fc_threshold
    s_big  <- abs(sdelta) > skewness_threshold

    # Concentration: Height rises, Intensity does not keep pace
    if (h_up && !i_up && !s_big) {
      return(c("Concentration", "Signal concentrates at summit (Height increases without Intensity gain)"))
    }

    # Flattening: Height drops, Intensity does not collapse proportionally
    if (h_down && !i_down && !s_big) {
      return(c("Flattening", "Signal spreads evenly across domain (Height decreases while Intensity stable)"))
    }

    # Polarity Shift: Skewness dominates, other dims stable
    if (s_big && !i_up && !i_down && !h_up && !h_down) {
      direction <- if (sdelta > 0) "right-skewed" else "left-skewed"
      return(c("Polarity Shift",
               sprintf("Signal redistributes within domain (%s)", direction)))
    }

    # Global Gain: Intensity up with concordant/stable Height, no skewness shift
    if (i_up && (h_up || !h_down) && !s_big) {
      return(c("Global Gain", "Coordinated increase in signal abundance"))
    }

    # Global Loss: Intensity down with concordant/stable Height, no skewness shift
    if (i_down && (h_down || !h_up) && !s_big) {
      return(c("Global Loss", "Coordinated decrease in signal abundance"))
    }

    # Complex: multiple dimensions shifting in discordant directions
    n_changed <- sum(i_up, i_down, h_up, h_down, s_big)
    if (n_changed >= 2) {
      return(c("Complex", "Multiple geometric dimensions shift discordantly"))
    }

    c("Stable", "No significant geometric change detected")
  }

  class_mat <- vapply(seq_along(peak_ids), function(j) {
    classify_one(intensity_fc[j], height_fc[j], skewness_delta[j])
  }, character(2))

  shift_res$Intensity_FC  <- intensity_fc
  shift_res$Height_FC     <- height_fc
  shift_res$Skewness_Delta <- skewness_delta
  shift_res$Shape_Class   <- class_mat[1, ]
  shift_res$Class_Label   <- class_mat[2, ]

  # Summary message -------------------------------------------------------------------
  class_counts <- table(shift_res$Shape_Class)
  msg <- paste(sprintf("%s: %d", names(class_counts), class_counts), collapse = ", ")
  message("Shape Shift classification complete: ", msg)

  return(shift_res)
}
