# ---- Academic visualization suite (epiPortrait_object_visualization) --------
# Follows the academic visualization design: every plot is a view of the
# SummarizedExperiment object; plotting never re-computes biological calls.

# Publication palette: same biological class always the same colour. Keys are
# the DISPLAY labels, so the palette adapts to mark-aware terminology
# (e.g. H3K27me3: Extended-Domain / Dual-Extreme) while colours stay fixed.
.epi_publication_palette <- function(mark = NULL) {
  cols <- c("Typical"          = "#B8B8B8",
            "Intensity-Super"  = "#D55E00",
            "Breadth-Super"    = "#0072B2",
            "Dual-Super"       = "#009E73",
            "Uncertain"        = "#666666")
  if (!is.null(mark)) {
    preset <- tryCatch(get_mark_preset(mark), error = function(e) NULL)
    if (!is.null(preset)) {
      # map by NAME, not by position: the canonical colour is bound to the
      # canonical class key, and the value is the mark-aware display label.
      # A positional assignment names(cols) <- preset$class_labels would
      # scramble the colours when class_labels are not in the same order as
      # cols (e.g. marks whose labels start with the super classes).
      lab <- stats::setNames(preset$class_labels,
                             c("Intensity-Super", "Breadth-Super",
                               "Dual-Super", "Typical", "Uncertain"))
      names(cols) <- unname(lab[names(cols)])
    }
  }
  cols
}

.epi_status_palette <- function() {
  c("Super" = "#D55E00", "Typical" = "#B8B8B8", "Uncertain" = "#666666")
}

# Blend a colour towards white. white_fraction = 0 keeps the full colour,
# 1 returns pure white. Used for support-level gradients on the
# replicate-support bars (full saturation only at maximum support).
.epi_blend_with_white <- function(col, white_fraction) {
  m <- grDevices::col2rgb(col)[, 1]
  out <- m + white_fraction * (255 - m)
  grDevices::rgb(out[1], out[2], out[3], maxColorValue = 255)
}

# Display labels of the three super classes in canonical Intensity/Breadth/Dual
# order (mark-aware). The publication palette is ordered with "Typical" first,
# so the super classes must NOT be taken as names(pal)[1:3].
.epi_super_class_labels <- function(mark = NULL) {
  preset <- NULL
  if (!is.null(mark)) {
    preset <- tryCatch(get_mark_preset(mark), error = function(e) NULL)
  }
  unname(.epi_class_display(c("Intensity-Super", "Breadth-Super", "Dual-Super"),
                            preset))
}

.epi_group_palette <- function(groups) {
  grp <- unique(as.character(groups))
  n <- length(grp)
  cols <- if (n <= 3) {
    c("#D55E00", "#0072B2", "#009E73")
  } else {
    grDevices::colorRampPalette(c("#D55E00", "#0072B2", "#009E73"))(n)
  }
  stats::setNames(cols, grp)
}

.epi_theme_publication <- function(base_size = 10.5) {
  # Nature-style publication theme: no background grid lines, thin black axis
  # lines, Arial sans-serif, compact sizes, frameless legend (nature-figure
  # skill R quick-start).
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(color = "grey35", size = base_size - 1),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.spacing.x = ggplot2::unit(0.15, "cm"),
      legend.spacing.y = ggplot2::unit(0.05, "cm"),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      strip.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Wrap a colour/fill legend into `nrow` rows so long legends never overflow
# the plot width when positioned at the top (ggplot2 does not auto-wrap top
# legends). `n_items` estimates the number of keys; a sensible nrow is chosen
# automatically when nrow is NULL.
.epi_wrap_legend <- function(scale = "colour", nrow = NULL, n_items = NULL) {
  if (is.null(nrow)) {
    nrow <- if (is.null(n_items)) 1 else ceiling(n_items / 4)
    nrow <- max(nrow, 1)
  }
  ggplot2::guides(
    setNames(
      list(ggplot2::guide_legend(nrow = nrow, byrow = TRUE)),
      scale
    )
  )
}

# Display class: render legacy "Width-Super" as "Breadth-Super" for backward
# compatibility with pre-v1.0 objects, and apply mark-aware display labels
# when a preset is provided (e.g. H3K27me3 -> Extended-Domain / Dual-Extreme).
# Internal class values are never changed (heterochromatin plan 2026-08-10).
.epi_display_class <- function(x, mark = NULL) {
  if (!is.null(mark)) {
    preset <- tryCatch(get_mark_preset(mark), error = function(e) NULL)
    if (!is.null(preset)) return(.epi_class_display(x, preset))
  }
  ifelse(!is.na(x) & x == "Width-Super", "Breadth-Super", x)
}

.epi_resolve_feature <- function(se, feature) {
  avail <- assayNames(se)
  if (feature %in% avail) return(feature)
  alt <- .resolve_assay(se, feature)
  if (alt %in% avail) return(alt)
  stop(sprintf("Feature '%s' not found in assays.", feature))
}

.epi_feature_vector <- function(se, feature, group = NULL, group_var = "Condition") {
  feat <- .epi_resolve_feature(se, feature)
  v <- rowMeans(assay(se, feat), na.rm = TRUE)
  if (!is.null(group)) {
    meta <- as.data.frame(colData(se))
    if (!group_var %in% colnames(meta)) {
      stop(sprintf("group_var '%s' not found in colData.", group_var))
    }
    if (!group %in% unique(meta[[group_var]])) {
      stop(sprintf("group '%s' not found in colData.", group))
    }
    idx <- which(meta[[group_var]] == group)
    if (length(idx) > 0) v <- rowMeans(assay(se, feat)[, idx, drop = FALSE], na.rm = TRUE)
  }
  v
}

.epi_transform_value <- function(x, transform = c("none", "log10p1")) {
  transform <- match.arg(transform)
  if (transform == "log10p1") return(log10(pmax(x, 0) + 1))
  x
}

.epi_class_col <- function(se, group = NULL) {
  if (is.null(group)) {
    "Combined_Domain_Class"
  } else {
    paste0("Combined_Class__", group)
  }
}


#' Plot Intensity x Breadth Architecture Landscape
#'
#' @description Signature figure: each domain as a point in the
#' Intensity x Breadth plane, coloured by combined architecture class.
#' This is the visual that distinguishes epiPortrait from single-intensity
#' callers (Intensity-extreme vs breadth-extreme vs dual).
#'
#' @param se A SummarizedExperiment with combined class columns
#'   (\code{combine_superdomain_calls}).
#' @param group Character or NULL. Condition group whose combined class to
#'   colour by (NULL uses the consensus \code{Combined_Domain_Class}).
#' @param transform Character. "log10p1" (default) or "none" (display only;
#'   classification is never recomputed).
#' @param label_n Integer. Number of representative domains to label per
#'   super class (for \code{label_mode = "per_class"}) or number of
#'   top-intensity domains (for \code{label_mode = "top_intensity"}).
#'   Default 0 = no domain labels (the architecture scatter is the message;
#'   labels are opt-in via \code{label_n > 0}).
#' @param label_col Character. Column in rowData used for labels (default
#'   "SYMBOL" if present, else Domain_ID).
#' @param label_mode Character. \code{"per_class"} (default): label the most
#'   representative (prototype) domains of each super class — top intensity
#'   within Intensity-Super, top breadth within Breadth-Super, top combined
#'   rank within Dual-Super — so each architecture class is named in the
#'   figure. \code{"top_intensity"}: label the globally top-intensity domains
#'   (legacy behaviour).
#' @param show_counts Logical. Append the per-class domain count to the legend
#'   labels, e.g. "Intensity-Super (n=607)" (default TRUE).
#' @param x_feature Character. Assay for the x axis (default "Intensity").
#' @param y_feature Character. Assay for the y axis (default
#'   "NativeMaxPeakWidth", the domain-level native breadth summary; falls back
#'   to "NativeOccupiedWidth" then "SignalDispersion" when unavailable).
#' @param mark Character or NULL. Histone mark preset for display terminology
#'   (e.g. "H3K27me3" renders Extended-Domain / Dual-Extreme labels). NULL uses
#'   the canonical Intensity-Super / Breadth-Super / Dual-Super labels.
#' @param group_var Character. colData column used for grouping when
#'   \code{group} is provided (default "Condition"). Must match the group_var
#'   used by \code{call_super_domains(mode = "per_group")} /
#'   \code{combine_superdomain_calls()} for the requested group.
#' @return A ggplot.
#' @import ggplot2
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' se <- call_super_domains(se, feature = "Breadth",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' se <- combine_superdomain_calls(se, group_var = "Condition")
#' plot_domain_landscape(se, group = "Control")
#' @export
plot_domain_landscape <- function(se, group = NULL,
                                  transform = c("log10p1", "none"),
                                  label_n = 0, label_col = NULL,
                                  label_mode = c("per_class", "top_intensity"),
                                  x_feature = "Intensity",
                                  y_feature = "NativeMaxPeakWidth",
                                  mark = NULL,
                                  group_var = "Condition",
                                  show_counts = TRUE) {
  transform <- match.arg(transform)
  label_mode <- match.arg(label_mode)
  rd <- as.data.frame(rowData(se), optional = TRUE)
  cc <- .epi_class_col(se, group)
  if (!cc %in% colnames(rd)) {
    stop("Run combine_superdomain_calls() first (with group_var for per-group classes).")
  }
  cl <- .epi_display_class(rd[[cc]], mark = mark)
  # y-axis: NativeMaxPeakWidth preferred, then NativeOccupiedWidth, then
  # SignalDispersion (design: Breadth-Super calls come from peak-level native
  # width, so the landscape y-axis uses the domain-level native breadth summary).
  if (!y_feature %in% assayNames(se)) {
    for (alt in c("NativeMaxPeakWidth", "NativeOccupiedWidth", "SignalDispersion")) {
      if (alt %in% assayNames(se)) { y_feature <- alt; break }
    }
  }
  df <- data.frame(
    Domain_ID = rownames(se),
    Intensity = .epi_transform_value(.epi_feature_vector(se, x_feature, group, group_var), transform),
    Breadth = .epi_transform_value(.epi_feature_vector(se, y_feature, group, group_var), transform),
    Class = factor(cl, levels = names(.epi_publication_palette(mark)))
  )
  if (is.null(label_col)) {
    df$Label <- if ("SYMBOL" %in% colnames(rd)) rd$SYMBOL else rownames(se)
  } else if (label_col %in% colnames(rd)) {
    df$Label <- rd[[label_col]]
  } else {
    df$Label <- rownames(se)
  }

  # Per-class domain counts appended to the legend labels (e.g. "n=607").
  pal <- .epi_publication_palette(mark)
  if (show_counts) {
    cnt <- table(factor(cl, levels = names(pal)))
    leg_labels <- sprintf("%s (n=%d)", names(cnt), as.vector(cnt))
  } else {
    leg_labels <- names(pal)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Intensity, y = Breadth,
                                        color = Class)) +
    ggplot2::geom_point(size = 1.3, alpha = 0.7) +
    ggplot2::scale_color_manual(values = pal,
                                name = "Architecture class",
                                labels = leg_labels) +
    .epi_theme_publication() +
    ggplot2::labs(
      title = "Domain Architecture Landscape",
      x = if (transform == "log10p1") paste0(x_feature, " (log10)") else x_feature,
      y = if (transform == "log10p1") paste0(y_feature, " (log10)") else y_feature)

  if (label_n > 0) {
    if (label_mode == "top_intensity") {
      top <- df[order(-df$Intensity), ][seq_len(min(label_n, nrow(df))), ]
    } else {
      # per-class prototype labeling: pick the most representative (top-ranked)
      # domains of each super class, so every architecture class is named.
      super_names <- .epi_super_class_labels(mark)
      frac_rank <- function(x) {
        r <- rank(x, na.last = "keep")
        m <- max(r, na.rm = TRUE)
        if (!is.finite(m) || m == 0) rep(0, length(r)) else r / m
      }
      top <- do.call(rbind, lapply(super_names, function(cl) {
        sub <- df[df$Class == cl, , drop = FALSE]
        if (nrow(sub) == 0) return(NULL)
        pos <- which(cl == super_names)[1]
        score <- switch(pos,
          frac_rank(sub$Intensity),
          frac_rank(sub$Breadth),
          (frac_rank(sub$Intensity) + frac_rank(sub$Breadth)) / 2)
        sub$Score <- score
        sub[order(-sub$Score), ][seq_len(min(label_n, nrow(sub))), , drop = FALSE]
      }))
      if (is.null(top) || nrow(top) == 0) top <- data.frame()
    }
    if (nrow(top) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = top, ggplot2::aes(label = Label),
        size = 3, max.overlaps = 30, box.padding = 0.4,
        point.padding = 0.2, color = "black", segment.color = "grey50")
    }
  }
  p + .epi_wrap_legend("colour", n_items = length(pal))
}


#' Plot Combined Architecture Transition Matrix
#'
#' @description Heatmap of the combined architecture class transition between
#' two conditions (Control rows -> Treatment columns). Displays row percentage
#' and domain counts. This is a RELATIVE architecture-state transition (per
#' condition-group cutoffs) and must not be read as absolute gain/loss.
#'
#' @param se A SummarizedExperiment with \code{Combined_Class__<group>} columns.
#' @param ref_group Character. Reference condition (rows).
#' @param target_group Character. Treatment / target condition (columns).
#' @param normalize Character. "row_percent" (default), "count", or "none".
#' @param include_uncertain Logical. Include Uncertain class (default TRUE).
#' @param mark Character or NULL. Histone mark preset for display terminology.
#' @param show_axis_labels Logical. Show the class names on the heatmap
#'   margins (default TRUE). The axis lines and ticks are always removed, so
#'   the plot reads as a clean heatmap that keeps the row/column class
#'   identity; set \code{FALSE} for a label-free heatmap.
#' @param ref_role Character or NULL. Optional semantic role label appended to
#'   the reference group name in the subtitle (e.g. "control"), so the plot
#'   states unambiguously which group is the reference. Default NULL.
#' @param target_role Character or NULL. Optional semantic role label appended
#'   to the target group name in the subtitle (e.g. "tumor"). Default NULL.
#' @return A ggplot.
#' @import ggplot2
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' se <- call_super_domains(se, feature = "Breadth",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' se <- combine_superdomain_calls(se, group_var = "Condition")
#' plot_class_transition(se, ref_group = "Control", target_group = "Treatment")
#' @export
plot_class_transition <- function(se, ref_group, target_group,
                                  normalize = c("row_percent", "count", "none"),
                                  include_uncertain = TRUE,
                                  mark = NULL,
                                  ref_role = NULL, target_role = NULL,
                                  show_axis_labels = TRUE) {
  normalize <- match.arg(normalize)
  rc <- paste0("Combined_Class__", ref_group)
  tc <- paste0("Combined_Class__", target_group)
  rd <- as.data.frame(rowData(se), optional = TRUE)
  if (!rc %in% colnames(rd) || !tc %in% colnames(rd)) {
    stop("Run combine_superdomain_calls(se, group_var = ...) first.")
  }
  r <- .epi_display_class(rd[[rc]], mark = mark)
  t <- .epi_display_class(rd[[tc]], mark = mark)
  tab <- table(r, t)
  if (!include_uncertain) {
    tab <- tab[!rownames(tab) %in% "Uncertain", , drop = FALSE]
    tab <- tab[, !colnames(tab) %in% "Uncertain", drop = FALSE]
  }
  if (normalize == "row_percent") {
    tab_p <- sweep(tab, 1, rowSums(tab), "/") * 100
    value <- as.vector(tab_p)
    label <- sprintf("%.1f%%\n(n=%d)", value, as.vector(tab))
  } else {
    value <- as.vector(tab)
    label <- sprintf("%d", value)
  }
  df <- data.frame(
    From = rep(rownames(tab), ncol(tab)),
    To = rep(colnames(tab), each = nrow(tab)),
    Value = value, Label = label, stringsAsFactors = FALSE
  )
  df$From <- factor(df$From, levels = rownames(tab))
  df$To <- factor(df$To, levels = colnames(tab))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = To, y = From, fill = Value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = Label), size = 3, color = "black") +
    ggplot2::scale_fill_gradient2(low = "#F7FBFF", high = "#08519C", mid = "#6BAED6",
                                  midpoint = if (normalize == "row_percent") 50 else NA,
                                  name = if (normalize == "row_percent")
                                    "Row %" else "Domain count") +
    .epi_theme_publication() +
    ggplot2::labs(
      title = "Architecture Class Transition",
      subtitle = .epi_transition_subtitle(ref_group, target_group,
                                          ref_role, target_role),
      x = if (is.null(target_role)) sprintf("%s (target)", target_group)
          else sprintf("%s (%s, target)", target_group, target_role),
      y = if (is.null(ref_role)) sprintf("%s (reference)", ref_group)
          else sprintf("%s (%s, reference)", ref_group, ref_role)) +
    ggplot2::theme(axis.text.x = if (show_axis_labels)
                     ggplot2::element_text(angle = 45, hjust = 1)
                   else ggplot2::element_blank(),
                   axis.text.y = if (show_axis_labels)
                     ggplot2::element_text() else ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank())
  p
}


# Subtitle for the class-transition heatmap. Optionally annotates the semantic
# role of each group (e.g. control / tumor) so the arrow direction is
# unambiguous when group names do not express it.
.epi_transition_subtitle <- function(ref_group, target_group,
                                     ref_role = NULL, target_role = NULL) {
  if (is.null(ref_role) && is.null(target_role)) {
    return(sprintf("%s -> %s (relative, per-group cutoffs)",
                   ref_group, target_group))
  }
  ref_lab <- if (is.null(ref_role)) ref_group else sprintf("%s (%s)", ref_group, ref_role)
  tgt_lab <- if (is.null(target_role)) target_group else sprintf("%s (%s)", target_group, target_role)
  sprintf("%s -> %s (relative, per-group cutoffs)", ref_lab, tgt_lab)
}

#' Plot Replicate-Support Distribution
#'
#' @description Bar plot of the fraction of domains with each replicate-support
#' level (0/k, 1/k, ..., k/k), stacked by the final group call. Demonstrates
#' replicate-aware calling and the effect of the support rule.
#'
#' @param se A SummarizedExperiment after \code{call_super_domains} per_group.
#' @param feature Character. Feature (default "Intensity").
#' @param group Character. Condition group (default first group).
#' @param group_var Character. colData column for groups.
#' @param min_support Integer. Drop support levels with fewer than this many
#'   supporting replicates. The 0/k bin (domains never super in any replicate)
#'   usually dominates the bar and is grey; set \code{min_support = 1} (default)
#'   to zoom into the replicate-aware distribution.
#' @param show_counts Logical. Draw the domain count above each bar
#'   (default TRUE).
#' @return A ggplot.
#' @import ggplot2
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' plot_replicate_support(se, feature = "Intensity", group = "Control")
#' @export
plot_replicate_support <- function(se, feature = "Intensity", group = NULL,
                                    group_var = "Condition",
                                    min_support = 1, show_counts = TRUE) {
  # Breadth support lives in rowData (Breadth_Support__<group>), not in an
  # assay, so bypass .epi_resolve_feature() for it (heterochromatin review).
  feat <- if (identical(feature, "Breadth")) "Breadth" else .epi_resolve_feature(se, feature)
  rd <- as.data.frame(rowData(se), optional = TRUE)
  meta <- as.data.frame(colData(se))
  if (!group_var %in% colnames(meta)) stop("group_var not found in colData.")
  groups <- unique(meta[[group_var]])
  if (is.null(group)) group <- groups[1]
  if (!group %in% groups) stop(sprintf("group '%s' not found.", group))
  n_rep <- sum(meta[[group_var]] == group)

  support_col <- sprintf("%s_Support__%s", feat, group)
  call_col <- sprintf("%s_Call__%s", feat, group)
  if (!support_col %in% colnames(rd)) {
    stop("Run call_super_domains(mode='per_group') first.")
  }
  sup <- rd[[support_col]]
  n_super <- round(sup * n_rep)
  level <- sprintf("%d/%d", n_super, n_rep)
  status <- ifelse(is.na(rd[[call_col]]), "Uncertain",
                   ifelse(grepl("_Super", rd[[call_col]]), "Super", "Typical"))
  df <- data.frame(Level = level, Status = status, stringsAsFactors = FALSE)
  df$Level <- factor(df$Level, levels = sprintf("%d/%d", 0:n_rep, n_rep))

  all_levels <- sprintf("%d/%d", 0:n_rep, n_rep)
  keep_levels <- all_levels[as.integer(sub("/.*", "", all_levels)) >= min_support]
  # Defensive: never return an empty plot if every domain is in the 0/k bin.
  if (length(keep_levels) == 0 || !any(df$Level %in% keep_levels)) {
    keep_levels <- all_levels
    if (min_support > 0) {
      warning("No domains exceed min_support; showing all support levels.",
              call. = FALSE)
    }
  }
  df <- df[df$Level %in% keep_levels, , drop = FALSE]
  df$Level <- factor(df$Level, levels = keep_levels)

  # Feature-specific Super hue, consistent with the landscape palette
  # (Intensity-Super = orange, Breadth-Super = blue, Dual-Super = green).
  super_hue <- switch(feat,
                      "Breadth" = "#0072B2",
                      "Width"   = "#0072B2",
                      "Intensity" = "#D55E00",
                      "SignalDispersion" = "#009E73",
                      "#D55E00")

  # Per-domain fill key = Status x support level; the Super colour is
  # desaturated towards white as support decreases (full colour at k/k).
  # Build the fill keys in a DETERMINISTIC order (Status first, then support
  # level) so the manual scale and the legend are never scrambled by row order.
  df$SupportK <- as.integer(sub("/.*", "", as.character(df$Level)))
  df$FillKey <- paste(df$Status, df$Level)
  key_order <- c(
    sprintf("Super %s", keep_levels),
    sprintf("Typical %s", keep_levels),
    sprintf("Uncertain %s", keep_levels)
  )
  key_order <- key_order[key_order %in% unique(df$FillKey)]
  # MUST be a factor whose levels exactly match the manual-scale names,
  # otherwise ggplot treats the fill as an unordered discrete scale and the
  # colour mapping / legend both break ("No shared levels found").
  df$FillKey <- factor(df$FillKey, levels = key_order)
  fill_vals <- vapply(seq_len(nrow(df)), function(i) {
    st <- df$Status[i]
    if (st == "Typical")  return(.epi_status_palette()[["Typical"]])
    if (st == "Uncertain") return(.epi_status_palette()[["Uncertain"]])
    # Super: blend full hue -> white as support fraction drops below k/k.
    frac <- df$SupportK[i] / n_rep
    .epi_blend_with_white(super_hue, 1 - frac)
  }, character(1))
  fill_map <- stats::setNames(fill_vals, as.character(df$FillKey))
  # Collapse duplicate keys (same Status x Level -> same colour) preserving the
  # deterministic key_order for the legend.
  fill_map <- fill_map[!duplicated(as.character(df$FillKey))]
  fill_map <- fill_map[key_order]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Level, fill = FillKey)) +
    ggplot2::geom_bar(position = "stack") +
    ggplot2::scale_fill_manual(values = fill_map, breaks = names(fill_map),
                               name = "Group call") +
    .epi_theme_publication()

  if (show_counts) {
    cnt <- as.data.frame(table(df$Level), stringsAsFactors = FALSE)
    colnames(cnt) <- c("Level", "N")
    cnt$Level <- factor(cnt$Level, levels = keep_levels)
    p <- p + ggplot2::geom_text(
      data = cnt, ggplot2::aes(x = Level, y = N, label = N),
      size = 2.6, vjust = -0.6, inherit.aes = FALSE)
  }

  p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      title = sprintf("Replicate Support (%s)", feat),
      subtitle = sprintf("%s, n = %d replicates%s", group, n_rep,
                         if (min_support > 0)
                           sprintf(" (support >= %d/%d shown)", min_support, n_rep) else ""),
      x = "Support level", y = "Number of domains",
      fill = "Group call")
  p + .epi_wrap_legend("fill", n_items = length(fill_map))
}


#' Plot Domain Class Composition
#'
#' @description 100% stacked bar of the combined architecture class
#' composition per condition group.
#'
#' @param se A SummarizedExperiment with \code{Combined_Class__<group>} columns.
#' @param group_var Character. colData column for groups.
#' @param mark Character or NULL. Histone mark preset for display terminology.
#' @param focus Character. \code{"super"} (default): restrict the 100% stacked
#'   bar to the super classes (Intensity/Breadth/Dual) so the composition is
#'   readable when non-super domains dominate the total, and annotate the super
#'   share of all domains above each bar. \code{"all"}: full composition
#'   including Typical / Uncertain.
#' @param show_counts Logical. Draw absolute domain counts inside each segment
#'   (default TRUE).
#' @return A ggplot.
#' @import ggplot2
#' @importFrom stats ave
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' se <- call_super_domains(se, feature = "Breadth",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' se <- combine_superdomain_calls(se, group_var = "Condition")
#' plot_domain_class_composition(se)
#' @export
plot_domain_class_composition <- function(se, group_var = "Condition",
                                          mark = NULL,
                                          focus = c("super", "all"),
                                          show_counts = TRUE) {
  focus <- match.arg(focus)
  meta <- as.data.frame(colData(se))
  if (!group_var %in% colnames(meta)) stop("group_var not found in colData.")
  groups <- unique(meta[[group_var]])
  rd <- as.data.frame(rowData(se), optional = TRUE)
  cl_cols <- paste0("Combined_Class__", groups)
  if (!all(cl_cols %in% colnames(rd))) {
    stop("Run combine_superdomain_calls(se, group_var = ...) first.")
  }
  pal <- .epi_publication_palette(mark)
  class_names <- names(pal)
  full <- do.call(rbind, lapply(groups, function(g) {
    cl <- .epi_display_class(rd[[paste0("Combined_Class__", g)]], mark = mark)
    tab <- table(factor(cl, levels = class_names))
    data.frame(Group = g, Class = names(tab), N = as.vector(tab), stringsAsFactors = FALSE)
  }))
  full$Class <- factor(full$Class, levels = class_names)
  full$Group <- factor(full$Group, levels = groups)

  super_classes <- .epi_super_class_labels(mark)
  df <- if (focus == "super") full[full$Class %in% super_classes, , drop = FALSE] else full
  renorm <- function(n) {
    s <- sum(n)
    if (is.na(s) || s == 0) rep(0, length(n)) else n / s
  }
  df$Prop <- ave(df$N, df$Group, FUN = renorm)
  df$LabelY <- ave(df$Prop, df$Group, FUN = function(p) cumsum(p) - p / 2)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Prop, fill = Class)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.7) +
    ggplot2::scale_y_continuous(labels = .epi_percent_labels) +
    ggplot2::scale_fill_manual(values = pal) +
    .epi_theme_publication()

  if (show_counts) {
    lbl <- df[df$N > 0 & (focus == "super" | df$Prop >= 0.02), , drop = FALSE]
    if (nrow(lbl) > 0) {
      p <- p + ggplot2::geom_text(data = lbl,
                                  ggplot2::aes(y = LabelY, label = N),
                                  size = 2.6, color = "black")
    }
  }
  if (focus == "super") {
    n_super <- tapply(full$N[full$Class %in% super_classes],
                      full$Group[full$Class %in% super_classes], sum)
    n_all <- tapply(full$N, full$Group, sum)
    frac <- ifelse(is.na(n_all) | n_all == 0, 0, n_super / n_all)
    ann <- data.frame(Group = factor(names(frac), levels = groups),
                      Frac = 100 * frac, stringsAsFactors = FALSE)
    p <- p + ggplot2::geom_text(data = ann,
                                ggplot2::aes(x = Group, y = 1.06,
                                             label = sprintf("%.1f%%", Frac)),
                                size = 2.6, color = "grey25", inherit.aes = FALSE) +
      ggplot2::coord_cartesian(ylim = c(0, 1.22))
  }

  p + ggplot2::labs(
    title = "Domain Class Composition",
    subtitle = if (focus == "super")
      "super classes only; % above bar = super share of all domains (Typical/Uncertain excluded)"
    else "100% stacked composition by condition group",
    x = group_var, y = "Proportion of domains", fill = "Class") +
    .epi_wrap_legend("fill", n_items = length(pal))
}

.epi_percent_labels <- function(x) paste0(round(100 * x), "%")


#' Plot Raw Domain Feature Profile
#'
#' @description Quantitative per-domain profile across conditions: biological
#'   replicates as points and group means as larger open points, faceted by
#'   feature. A dashed line connects the group means so the across-condition
#'   direction of change is read directly. When combined architecture classes
#'   are present (\code{combine_superdomain_calls}), the domain's class in each
#'   condition is annotated under the axis, linking the quantitative values to
#'   its architectural identity. Replaces the qualitative radar as the
#'   quantitative view.
#'
#' @param se A SummarizedExperiment.
#' @param peak_id Character. Domain ID.
#' @param group_var Character. colData column for groups.
#' @param features Character vector of assays to show (default the three
#'   dynamic assays).
#' @param transform Character. "log10p1" or "none".
#' @param annotate_class Logical. Annotate the domain's combined architecture
#'   class under each condition group (default TRUE; silently skipped when the
#'   combined-class columns are absent).
#' @param mark Character or NULL. Histone mark preset for class display
#'   terminology.
#' @return A ggplot.
#' @import ggplot2
#' @examples
#' data(example_se)
#' plot_domain_feature_profile(example_se, peak_id = rownames(example_se)[1])
#' @export
plot_domain_feature_profile <- function(se, peak_id, group_var = "Condition",
                                        features = NULL,
                                        transform = c("log10p1", "none"),
                                        annotate_class = TRUE,
                                        mark = NULL) {
  transform <- match.arg(transform)
  if (!peak_id %in% rownames(se)) stop("peak_id not found.")
  meta <- as.data.frame(colData(se))
  if (!group_var %in% colnames(meta)) stop("group_var not found in colData.")
  groups <- unique(as.character(meta[[group_var]]))
  if (is.null(features)) features <- intersect(
    c("Intensity", "SignalDispersion", "NativeMaxPeakWidth"), assayNames(se))
  features <- vapply(features, .epi_resolve_feature, character(1), se = se)

  rows <- lapply(features, function(f) {
    v <- assay(se, f)[peak_id, ]
    data.frame(Feature = f, Group = meta[[group_var]],
               Value = .epi_transform_value(v, transform),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)

  # Combined architecture class per condition (annotated under the axis when
  # the object carries combine_superdomain_calls() output).
  class_lab <- NULL
  if (annotate_class) {
    rd <- as.data.frame(rowData(se), optional = TRUE)
    cl_cols <- paste0("Combined_Class__", groups)
    if (all(cl_cols %in% colnames(rd))) {
      cls <- vapply(cl_cols, function(cc) {
        .epi_display_class(rd[[cc]][match(peak_id, rownames(se))], mark = mark)
      }, character(1))
      class_lab <- stats::setNames(as.character(cls), groups)
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Value, color = Group)) +
    ggplot2::stat_summary(fun = mean, geom = "line",
                          aes(group = 1), linewidth = 0.5,
                          linetype = "dashed", color = "grey45", inherit.aes = TRUE) +
    ggplot2::geom_point(size = 2, alpha = 0.7, position = ggplot2::position_jitter(0.05)) +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 4, shape = 21,
                          fill = "white", stroke = 1.2) +
    ggplot2::facet_wrap(~Feature, scales = "free_y") +
    ggplot2::scale_color_manual(values = .epi_group_palette(df$Group)) +
    .epi_theme_publication()

  if (!is.null(class_lab)) {
    p <- p + ggplot2::scale_x_discrete(labels = function(grp) {
      vapply(as.character(grp), function(g) {
        if (is.na(class_lab[[g]])) g else paste0(g, "\n", class_lab[[g]])
      }, character(1))
    })
  }

  p <- p + ggplot2::labs(
    title = sprintf("Domain Feature Profile: %s", peak_id),
    subtitle = sprintf(
      "replicates (points) and group means (open circles); dashed line = mean trend%s%s",
      if (!is.null(class_lab)) "; axis class label = combined architecture call per condition" else "",
      if (transform == "log10p1") " (log10, display only)" else ""),
    x = group_var,
    y = if (transform == "log10p1") "Value (log10)" else "Value",
    color = group_var) +
    ggplot2::theme(legend.position = "none")
  p
}


#' Plot the Cause Composition of Uncertain Calls
#'
#' @description Bar chart of why domains are \code{Uncertain} for a feature's
#'   per-group call (see \code{get_uncertain_cause}): no native peak in the
#'   condition / insufficient valid replicates / inflection no-call. Helps
#'   audit whether the \code{Uncertain} fraction reflects condition-specific
#'   loci rather than analysis artifacts.
#'
#' @param se A SummarizedExperiment after \code{call_super_domains(mode =
#'   "per_group")}.
#' @param feature Character. Feature (default "Breadth").
#' @param group Character. Condition group.
#' @param group_var Character or NULL. Grouping column (resolved from stored
#'   provenance when NULL).
#' @param show_counts Logical. Label bars with counts (default TRUE).
#' @return A ggplot.
#' @import ggplot2
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Breadth",
#'                          mode = "per_group", group_var = "Condition",
#'                          verbose = FALSE)
#' plot_uncertain_cause(se, feature = "Breadth", group = "Control")
#' @export
plot_uncertain_cause <- function(se, feature = "Breadth", group,
                                 group_var = NULL, show_counts = TRUE) {
  uc <- get_uncertain_cause(se, feature = feature, group = group,
                            group_var = group_var)
  df <- uc[!is.na(uc$Cause), , drop = FALSE]
  tab <- as.data.frame(table(Cause = df$Cause), stringsAsFactors = FALSE)
  tab <- tab[order(tab$Freq, decreasing = TRUE), , drop = FALSE]
  tab$Cause <- factor(tab$Cause, levels = rev(tab$Cause))

  p <- ggplot2::ggplot(tab, ggplot2::aes(x = Cause, y = Freq)) +
    ggplot2::geom_col(fill = "#D55E00", alpha = 0.85, width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    .epi_theme_publication() +
    ggplot2::labs(
      title = sprintf("Uncertain Cause (%s, %s)", feature, group),
      subtitle = sprintf("%d uncertain domains (%.1f%% of %d)",
                         sum(tab$Freq), 100 * sum(tab$Freq) / nrow(se), nrow(se)),
      x = NULL, y = "Number of domains")
  if (show_counts) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Freq),
                                hjust = -0.2, size = 3.2)
  }
  p
}


#' Save an epiPortrait Figure in Vector or Raster Format
#'
#' @param p A ggplot.
#' @param file Output path (extension determines format: .pdf, .svg, .png,
#'   .tiff).
#' @param width,height Numeric in inches.
#' @param dpi Numeric. For raster output (default 600).
#' @param ... Passed to the underlying ggsave / ragg device.
#' @return Invisibly the file path.
#' @examples
#' data(example_se)
#' p <- plot_portrait_pca(example_se, feature = "Intensity")
#' save_epiportrait_figure(p, file = tempfile(fileext = ".png"))
#' @export
save_epiportrait_figure <- function(p, file, width = 3.5, height = 3.3,
                                    dpi = 600, ...) {
  ext <- tolower(tools::file_ext(file))
  if (!ext %in% c("pdf", "svg", "png", "tiff", "tif")) {
    stop("file extension must be pdf, svg, png or tiff.")
  }
  if (ext %in% c("png", "tiff", "tif")) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      devfun <- if (ext == "png") ragg::agg_png else ragg::agg_tiff
      ggplot2::ggsave(file, plot = p, width = width, height = height,
                      units = "in", dpi = dpi, device = devfun, ...)
    } else {
      ggplot2::ggsave(file, plot = p, width = width, height = height, dpi = dpi, ...)
    }
  } else {
    ggplot2::ggsave(file, plot = p, width = width, height = height, ...)
  }
  invisible(file)
}
