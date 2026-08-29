# Rank-axis breaks for hockey-stick plots: target ~5 ticks (4-6) so the
# x-axis is neither sparse (the original complaint) nor cluttered on large
# domain universes. Step is a "nice" number (1/2/2.5/5/10 x 10^k) sized to the
# range; among steps that give 4-6 breaks the densest is preferred.
.epi_rank_breaks <- function(limits) {
  lo <- limits[1]; hi <- limits[2]
  span <- hi - lo
  if (!is.finite(span) || span <= 0) return(pretty(limits, n = 5))
  target <- span / 5
  e <- 10^floor(log10(target))
  if (!is.finite(e) || e <= 0) return(pretty(limits, n = 5))
  nices <- c(1, 2, 2.5, 5, 10) * e
  cand <- lapply(nices, function(s) {
    start <- ceiling(lo / s) * s
    brk <- seq(start, hi, by = s)
    brk[brk >= lo & brk <= hi]
  })
  ns <- vapply(cand, length, integer(1))
  in_range <- which(ns >= 4 & ns <= 6)
  pick <- if (length(in_range) > 0) in_range[1] else which.min(abs(ns - 5))
  cand[[pick]]
}


#' Plot Hockey Stick Curve for Super-Domain Identification
#'
#' @description Generates a ranking plot for super-domain identification. Domains are ranked by
#' their feature value; super-domains are highlighted and the cutoff (and, when
#' available, the replicate-cutoff stability band) is drawn from the object's
#' stored provenance. In per-group mode the cutoff is shown as a descriptive
#' median replicate cutoff, never as the unique biological decision boundary.
#'
#' @param se A SummarizedExperiment object after running \code{call_super_domains}.
#' @param feature Character. Which feature to plot ("Intensity",
#'   "SignalDispersion", or "Breadth"; legacy aliases accepted). For
#'   \code{"Breadth"} the plot shows the replicate-specific native PeakWidth
#'   rank curve from \code{metadata(se)$breadth_peak_calls}, and \code{group}
#'   must be a sample (replicate) name.
#' @param label_genes Logical. Label top super-domains (default FALSE; the
#'   ranking curve and cutoff are the message, per-domain labels are opt-in).
#' @param top_n_label Integer. Number of top domains to label (default: 10).
#' @param group Character or NULL. For signal features: a condition group
#'   (per-group call/rank columns \code{<feature>_Call__<group>}, descriptive
#'   median replicate cutoff), OR a single sample name when the object was
#'   called with \code{mode = "per_sample"} — the sample's own ranked
#'   distribution is plotted with its per-sample cutoff. For
#'   \code{feature = "Breadth"}: the replicate whose native PeakWidth
#'   distribution to plot (required). NULL (default) uses the consensus
#'   \code{<feature>_Domain_Type} for signal features.
#' @param group_var Character. colData column used for groups when
#'   \code{group} is provided (default "Condition").
#' @param show_cutoff Logical. Draw the cutoff as a dashed horizontal guide
#'   (default TRUE).
#' @param show_cutoff_band Logical. Draw a light horizontal band for the
#'   replicate cutoff stability interval (default TRUE).
#'
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats setNames median quantile
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity",
#'                          method = "tangent", log_transform = FALSE)
#' plot_hockey_stick(se, feature = "Intensity", label_genes = FALSE)
#' @export
plot_hockey_stick <- function(se, feature = "Intensity", label_genes = FALSE,
                              top_n_label = 10, group = NULL,
                              group_var = "Condition",
                              show_cutoff = TRUE, show_cutoff_band = TRUE) {

  # ---- Breadth branch (P1-1) ------------------------------------------------
  # Breadth-Super is a PEAK-LEVEL call: there is no shared-domain breadth rank
  # curve. We plot the replicate-specific native PeakWidth rank curve directly
  # from metadata(se)$breadth_peak_calls. The caller must specify a sample
  # (replicate); plotting a pooled/domain curve would be misleading.
  if (feature == "Breadth") {
    if (is.null(group) || !group %in% colnames(se)) {
      stop("For feature = 'Breadth', pass the replicate name via `group` ",
           "(e.g. group = 'Control_1'). Breadth-Super is a peak-level call and ",
           "has no shared-domain rank curve.")
    }
    pc <- S4Vectors::metadata(se)$breadth_peak_calls
    if (is.null(pc)) {
      stop("No Breadth peak-level calls found. Run call_super_domains(se, ",
           "feature = 'Breadth') first.")
    }
    pc_s <- pc[pc$SampleID == group, , drop = FALSE]
    if (nrow(pc_s) == 0) {
      stop(sprintf("No peak-level breadth calls for sample '%s'.", group))
    }
    cutoff <- unique(pc_s$PeakWidthCutoff)
    cutoff <- if (length(cutoff) == 1 && is.finite(cutoff)) cutoff else NA_real_
    ord <- order(pc_s$PeakWidth, pc_s$NativePeakID)
    plot_df <- data.frame(
      Peak_ID = pc_s$NativePeakID[ord],
      Value = pc_s$PeakWidth[ord],
      Rank = seq_len(nrow(pc_s)),
      Type = ifelse(pc_s$PeakBroadCall[ord] == "Broad",
                    "Breadth_Super_Element", "Breadth_Typical"),
      stringsAsFactors = FALSE
    )
    super_type <- "Breadth_Super_Element"
    typical_type <- "Breadth_Typical"
    display_label <- "Breadth-Super"
    plot_df$Type_Label <- ifelse(plot_df$Type == super_type, display_label,
                                 "Typical")
    plot_df$Type_Label <- factor(plot_df$Type_Label,
                                 levels = c(display_label, "Typical"))
    plot_df$SYMBOL <- plot_df$Peak_ID
    pal <- stats::setNames(c("#0072B2", "#B8B8B8"),
                           c(display_label, "Typical"))

    p <- ggplot(plot_df, aes(x = Rank, y = Value, color = Type_Label)) +
      geom_line(color = "grey80", linewidth = 0.5) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_manual(values = pal, name = "Class") +
      scale_x_continuous(breaks = .epi_rank_breaks) +
      .epi_theme_publication(base_size = 11) +
      labs(title = sprintf("Hockey Stick Plot: %s (%s)", feature, group),
           subtitle = sprintf("Native PeakWidth distribution - %s  -  %d broad (%.1f%% of %d peaks)",
                              display_label,
                              sum(!is.na(plot_df$Type) & plot_df$Type == super_type),
                              100 * sum(!is.na(plot_df$Type) & plot_df$Type == super_type) / nrow(plot_df),
                              nrow(plot_df)),
           x = "Peak Width Rank", y = "Peak Width (bp)")
    if (show_cutoff && is.finite(cutoff)) {
      p <- p + geom_hline(yintercept = cutoff, linetype = "dashed",
                          color = "grey25", linewidth = 0.6) +
        annotate("text", x = 0, y = cutoff, hjust = 0, vjust = -0.4,
                 label = sprintf("cutoff = %.0f bp", cutoff),
                 size = 3, color = "grey25")
    }
    return(p)
  }

  feature <- .resolve_assay(se, feature)
  rd <- as.data.frame(rowData(se), optional = TRUE)

  # Determine which call columns to use (consensus vs per-group).
  if (is.null(group)) {
    col_type <- paste0(feature, "_Domain_Type")
    col_rank <- paste0(feature, "_Rank")
  } else {
    col_type <- paste0(feature, "_Call__", group)
    col_rank <- paste0(feature, "_Rank__", group)
  }
  if (!col_type %in% colnames(rd)) {
    # P2-fix: after a per_group call there is no consensus _Domain_Type column;
    # the previous error ("Run call_super_domains first") was misleading because
    # the user HAD called it. Point directly to the available per-group columns.
    if (is.null(group)) {
      pg_cols <- grep(sprintf("^%s_Call__", feature), colnames(rd), value = TRUE)
      if (length(pg_cols) > 0) {
        avail_groups <- sub(sprintf("^%s_Call__", feature), "", pg_cols)
        stop(paste0(
          "No consensus column '", col_type, "' in the object. The calls ",
          "were made with mode='per_group'; pass `group` to plot a ",
          "condition group (available: ",
          paste(avail_groups, collapse = ", "), ")."))
      }
    }
    stop(sprintf("Run call_super_domains(se, feature='%s'%s) first.",
                 feature, if (is.null(group)) "" else
                   if (group %in% colnames(se))
                     sprintf(", mode='per_sample' (sample '%s')", group)
                   else sprintf(", mode='per_group' (group '%s')", group)))
  }

  # Extract value vector (feature may be an assay or a static rowData column).
  # `group` may be a condition group (per_group/consensus) OR a single sample
  # (per_sample mode): when it matches a sample name, plot that sample's own
  # ranked distribution.
  if (feature %in% assayNames(se)) {
    if (is.null(group)) {
      val_vector <- rowMeans(assay(se, feature), na.rm = TRUE)
    } else if (group %in% colnames(se)) {
      val_vector <- assay(se, feature)[, group]
    } else {
      meta <- as.data.frame(colData(se))
      idx <- which(meta[[group_var]] == group)
      if (length(idx) == 0) stop(sprintf("group '%s' not found in colData.", group))
      val_vector <- rowMeans(assay(se, feature)[, idx, drop = FALSE], na.rm = TRUE)
    }
  } else if (feature %in% colnames(rd)) {
    val_vector <- as.numeric(rd[[feature]])
  } else {
    stop(sprintf("Feature '%s' not found in assays or rowData.", feature))
  }

  plot_df <- data.frame(
    Peak_ID = rownames(se),
    Value = val_vector,
    Rank = rd[[col_rank]],
    Type = rd[[col_type]],
    stringsAsFactors = FALSE
  )
  plot_df$SYMBOL <- if ("SYMBOL" %in% colnames(rd)) rd$SYMBOL else rownames(se)

  # Caller labels: Super / Typical / NA (no-call -> Uncertain).
  # SignalDispersion is a secondary architecture descriptor, NOT a canonical
  # Super axis (v1.0 design): label its generic caller output as a plain
  # "Dispersion" extreme rather than a taxonomy class.
  super_type <- paste0(feature, "_Super_Element")
  typical_type <- paste0(feature, "_Typical")
  display_label <- switch(feature,
                          "Intensity"      = "Intensity-Super",
                          "Breadth"        = "Breadth-Super",
                          "Width"          = "Breadth-Super",
                          "SignalDispersion" = "Dispersion-Extreme",
                          "Super-Domain")
  plot_df$Type_Label <- ifelse(is.na(plot_df$Type), "Uncertain",
                               ifelse(plot_df$Type == super_type, display_label,
                                      ifelse(plot_df$Type == typical_type, "Typical", "Uncertain")))
  plot_df$Type_Label <- factor(plot_df$Type_Label,
                               levels = c(display_label, "Typical", "Uncertain"))
  # Super point colour: match the package publication palette so every plot
  # assigns the same biological class the same colour (Intensity = orange,
  # Breadth/Width = blue, SignalDispersion = green).
  super_col <- switch(feature,
                      "Intensity"      = "#D55E00",
                      "Breadth"        = "#0072B2",
                      "Width"          = "#0072B2",
                      "SignalDispersion" = "#009E73",
                      "#D55E00")
  pal <- c(stats::setNames(c(super_col, "#B8B8B8", "#666666"),
                           c(display_label, "Typical", "Uncertain")))

  # Cutoff guide line stays neutral grey for both features.
  cutoff_col <- "grey25"

  # Super-domain count for the subtitle annotation (plotted, not recomputed).
  n_super <- sum(!is.na(plot_df$Type) & plot_df$Type == super_type)

  p <- ggplot(plot_df, aes(x = Rank, y = Value, color = Type_Label)) +
    geom_line(color = "grey80", linewidth = 0.5) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = pal, name = "Class") +
    scale_x_continuous(breaks = .epi_rank_breaks) +
    .epi_theme_publication(base_size = 11) +
    labs(title = sprintf("Hockey Stick Plot: %s", feature),
         subtitle = if (is.null(group))
           sprintf("Identification of %s  -  %d super (%.1f%% of %d domains)",
                   display_label, n_super, 100 * n_super / nrow(plot_df),
                   nrow(plot_df))
         else if (group %in% colnames(se))
           sprintf("%s  - %s (single replicate)  -  %d super (%.1f%%)",
                   display_label, group, n_super, 100 * n_super / nrow(plot_df))
         else sprintf("%s  - %s (per-group replicate-aware)  -  %d super (%.1f%%)",
                      display_label, group, n_super,
                      100 * n_super / nrow(plot_df)),
         x = paste(feature, "Rank"),
         y = paste(feature, "Value"))

  # Cutoff guide and replicate-cutoff stability band, read from the object's
  # stored provenance (never recomputed here).
  if (show_cutoff || show_cutoff_band) {
    prov <- get_call_provenance(se, feature)
    cutoff <- NULL
    band <- NULL
    if (!is.null(prov)) {
      if (is.null(group)) {
        cutoff <- prov$cutoff
        band <- prov$cutoff_stability_interval
      } else if (!is.null(prov$replicates) && !is.null(prov$replicates[[group]])) {
        # per_sample mode: draw the per-sample cutoff / stability interval.
        cutoff <- prov$replicates[[group]]$cutoff
        band <- prov$replicates[[group]]$cutoff_stability_interval
      } else if (!is.null(prov$groups) && !is.null(prov$groups[[group]])) {
        rep_cut <- vapply(prov$groups[[group]]$replicate_calls,
                          function(r) r$cutoff, numeric(1))
        rep_cut <- rep_cut[is.finite(rep_cut)]
        if (length(rep_cut) > 0) {
          cutoff <- stats::median(rep_cut)
          band <- stats::quantile(rep_cut, c(0.25, 0.75), na.rm = TRUE)
        }
      }
    }
    if (show_cutoff_band && !is.null(band) && all(is.finite(band))) {
      p <- p + annotate("rect",
                        xmin = -Inf, xmax = Inf,
                        ymin = band[1], ymax = band[2],
                        fill = "grey70", alpha = 0.15)
    }
    if (show_cutoff && !is.null(cutoff) && is.finite(cutoff)) {
      p <- p + geom_hline(yintercept = cutoff, linetype = "dashed",
                          color = cutoff_col, linewidth = 0.6)
      if (is.null(group)) {
        p <- p + annotate("text", x = 0, y = cutoff, hjust = 0, vjust = -0.4,
                          label = sprintf("cutoff = %.1f", cutoff),
                          size = 3, color = cutoff_col)
      } else {
        p <- p + annotate("text", x = 0, y = cutoff, hjust = 0, vjust = -0.4,
                          label = sprintf("median replicate cutoff (descriptive) = %.1f", cutoff),
                          size = 3, color = cutoff_col)
      }
    }
  }

  # Label top-ranked super-domains
  if (label_genes) {
    top_peaks <- plot_df[!is.na(plot_df$Type) & plot_df$Type == super_type, ]
    top_peaks <- top_peaks[order(top_peaks$Value, decreasing = TRUE), ]
    top_peaks <- head(top_peaks, top_n_label)

    if (nrow(top_peaks) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = top_peaks,
        aes(label = SYMBOL),
        size = 3.5,
        box.padding = 0.5,
        point.padding = 0.3,
        fontface = "bold.italic",
        color = "black",
        segment.color = "grey50"
      )
    }
  }

  return(p)
}


#' Plot Physical Coverage Profile for a Specific Peak (Track View)
#'
#' @description Generates a track-like visualization of the raw BigWig coverage
#' for a specific domain across all samples, for direct visual validation.
#'
#' @param se A SummarizedExperiment object from build_portrait_matrix().
#' @param peak_id Character. The ID of the peak to plot.
#' @param group_var Character. The grouping variable in colData(se).
#' @param extend_bp Numeric. Number of base pairs to extend the view around the peak (default: 1000).
#' @param smooth_window Numeric. Window size for smoothing the coverage (default: 1, no smoothing).
#' @param free_y Logical. Use independent y-scales per facet panel (default
#'   FALSE = fixed y, so cross-sample intensity differences are visually
#'   preserved rather than auto-scaled away).
#' @param show_native_occupancy Logical. Overlay the per-sample
#'   NativeOccupiedWidth (sum of reduced native peak widths inside the domain)
#'   as horizontal bands (default FALSE). The band is centred on the candidate
#'   domain and is a per-sample geometry view, not a new biological call.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @import rtracklayer
#' @import SummarizedExperiment
#' @examples
#' if (.Platform$OS.type != "windows") {
#'   extdata <- system.file("extdata", package = "epiPortrait")
#'   samples <- data.frame(
#'     SampleID  = c("C1", "C2", "T1", "T2"),
#'     Condition = c("Control", "Control", "Treatment", "Treatment"),
#'     bw_path   = file.path(extdata, c("C1.bw", "C2.bw", "T1.bw", "T2.bw")))
#'   se_tiny <- build_portrait_matrix(
#'     samples, consensus_peaks = rtracklayer::import(file.path(extdata, "peaks.bed")),
#'     workers = 1)
#'   plot_peak_track(se_tiny, peak_id = rownames(se_tiny)[1],
#'                   group_var = "Condition")
#' }
#' @export
plot_peak_track <- function(se, peak_id, group_var = "Condition", extend_bp = 1000,
                            smooth_window = 1, free_y = FALSE,
                            show_native_occupancy = FALSE) {
  
  if (!peak_id %in% rownames(se)) stop("Peak ID not found.")
  
  # 1. Extract region coordinates with flanking extension
  region <- rowRanges(se)[peak_id]
  plot_region <- region
  start(plot_region) <- max(1, start(region) - extend_bp)
  end(plot_region) <- end(region) + extend_bp
  
  # 2. Prepare sample metadata
  meta <- as.data.frame(colData(se))
  if (!"bw_path" %in% colnames(meta)) stop("bw_path column missing in colData(se).")
  
  message(sprintf("Loading BigWig signals for %s...", peak_id))
  
  # 3. Extract BigWig signals for all samples
  track_list <- lapply(seq_len(nrow(meta)), function(i) {
    bw_path <- meta$bw_path[i]
    sample_id <- meta$SampleID[i]
    group <- meta[[group_var]][i]
    
    # Read only the region of interest
    cvg <- rtracklayer::import(bw_path, which = plot_region, as = "NumericList")[[1]]
    
    # Construct plot data frame
    df <- data.frame(
      Position = seq(start(plot_region), end(plot_region)),
      Score = as.numeric(cvg),
      Sample = sample_id,
      Group = group
    )
    
    # Optional rolling-mean smoothing; backfill boundary NAs with original values.
    # Smoothing only affects the DISPLAY, never any stored feature.
    if (smooth_window > 1) {
      smoothed <- stats::filter(df$Score, rep(1 / smooth_window, smooth_window), sides = 2)
      na_idx <- which(is.na(smoothed))
      smoothed[na_idx] <- df$Score[na_idx]
      df$Score <- as.numeric(smoothed)
    }
    return(df)
  })
  
  plot_df <- do.call(rbind, track_list)
  
  # 4. Native occupancy span per sample (from the assay, centred on the domain)
  occ_bands <- NULL
  if (show_native_occupancy) {
    if ("NativeOccupiedWidth" %in% assayNames(se)) {
      occ <- as.numeric(assay(se, "NativeOccupiedWidth")[peak_id, ])
      dom_start <- start(region)
      dom_end <- end(region)
      dom_mid <- (dom_start + dom_end) / 2
      half <- pmax(occ, 1) / 2
      occ_bands <- data.frame(
        Sample = meta$SampleID,
        Occ_Start = pmax(dom_mid - half, dom_start),
        Occ_End = pmin(dom_mid + half, dom_end),
        stringsAsFactors = FALSE
      )
    } else {
      warning("NativeOccupiedWidth assay not present; skipping native-occupancy overlay.")
    }
  }
  
  # 5. Render coverage track
  p <- ggplot(plot_df, aes(x = Position, y = Score)) +
    # candidate domain background (very light grey) spanning all panels
    annotate("rect",
             xmin = start(region), xmax = end(region),
             ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5)
  if (!is.null(occ_bands)) {
    # per-sample NativeOccupiedWidth span (light green band)
    p <- p + geom_rect(
      data = occ_bands,
      aes(xmin = Occ_Start, xmax = Occ_End, ymin = -Inf, ymax = Inf),
      fill = "#009E73", alpha = 0.15, inherit.aes = FALSE)
  }
  p <- p +
    geom_area(aes(fill = Group), alpha = 0.8) +
    facet_wrap(~Sample, ncol = 1, scales = if (free_y) "free_y" else "fixed") +
    theme_classic(base_size = 12, base_family = "sans") +
    # P2-fix: scale_fill_brewer("Set1") capped at 9 groups; use the package's
    # own condition palette (ramped beyond 3 groups) for consistency.
    scale_fill_manual(values = .epi_group_palette(plot_df$Group)) +
    labs(title = sprintf("Coverage Track: %s", peak_id),
         subtitle = sprintf("Region: %s:%d-%d%s",
                            seqnames(region), start(region), end(region),
                            if (show_native_occupancy) "  - green band = NativeOccupiedWidth span" else ""),
         x = "Genomic Coordinates", y = "Signal Intensity (BigWig)") +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.text = element_text(colour = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) +
    .epi_wrap_legend("fill", n_items = length(unique(plot_df$SampleID)))
    
  return(p)
}
