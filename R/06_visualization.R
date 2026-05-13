#' Plot Hockey Stick Curve for Super-Domain Identification
#'
#' @description Generates a ROSE-style ranking plot. Peaks are ranked by their feature 
#' value (Intensity, Width, or Height). Super-domains above the inflection point 
#' are highlighted.
#'
#' @param se A SummarizedExperiment object after running \code{call_super_domains}.
#' @param feature Character. Which feature to plot ("Intensity", "Width", or "Height").
#' @param label_genes Logical. If SYMBOL column exists in rowData, label top super-domains?
#' @param top_n_label Integer. Number of top domains to label (default: 10).
#'
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats setNames
#' @export
plot_hockey_stick <- function(se, feature = "Intensity", label_genes = TRUE, top_n_label = 10) {
  
  col_type <- paste0(feature, "_Domain_Type")
  col_rank <- paste0(feature, "_Rank")
  
  if (!col_type %in% colnames(rowData(se))) {
    stop(sprintf("Please run call_super_domains(se, feature='%s') first.", feature))
  }
  
  # 提取绘图数据
  mat <- assay(se, feature)
  val_vector <- rowMeans(mat, na.rm = TRUE)
  
  plot_df <- data.frame(
    Peak_ID = rownames(se),
    Value = val_vector,
    Rank = rowData(se)[[col_rank]],
    Type = rowData(se)[[col_type]]
  )
  
  # 如果有基因注释，则加入
  if ("SYMBOL" %in% colnames(rowData(se))) {
    plot_df$SYMBOL <- rowData(se)$SYMBOL
  } else {
    plot_df$SYMBOL <- plot_df$Peak_ID
  }
  
  # 确定颜色
  target_label <- switch(feature,
                         "Intensity" = "Super-Element",
                         "Width"     = "Broad-Domain",
                         "Height"    = "Steep-Peak",
                         "Super-Domain")
  
  p <- ggplot(plot_df, aes(x = Rank, y = Value, color = Type)) +
    geom_line(color = "grey80", size = 0.5) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = stats::setNames(c("#E41A1C", "#377EB8"), c(target_label, "Typical"))) +
    theme_classic(base_size = 14) +
    labs(title = sprintf("Hockey Stick Plot: %s", feature),
         subtitle = paste0("Identification of ", target_label, "s"),
         x = paste(feature, "Rank"),
         y = paste(feature, "Value (Mean across samples)")) +
    theme(legend.position = "bottom")
  
  # Label top-ranked Super-Element / Broad-Domain / Steep-Peak
  if (label_genes) {
    top_peaks <- plot_df[plot_df$Type == target_label, ]
    top_peaks <- top_peaks[order(top_peaks$Value, decreasing = TRUE), ]
    top_peaks <- head(top_peaks, top_n_label)
    
    p <- p + ggrepel::geom_text_repel(
      data = top_peaks,
      aes(label = SYMBOL),
      size = 4,
      box.padding = 0.5,
      point.padding = 0.3,
      fontface = "bold.italic",
      color = "black",
      segment.color = "grey50"
    )
  }
  
  return(p)
}


#' Plot Physical Coverage Profile for a Specific Peak (Track View)
#'
#' @description Generates a track-like visualization of the raw BigWig coverage for a specific 
#' peak across all samples. This allows for direct visual validation of the "Shape Shift" 
#' detected by epiPortrait.
#'
#' @param se A SummarizedExperiment object from build_portrait_matrix().
#' @param peak_id Character. The ID of the peak to plot.
#' @param group_var Character. The grouping variable in colData(se).
#' @param extend_bp Numeric. Number of base pairs to extend the view around the peak (default: 1000).
#' @param smooth_window Numeric. Window size for smoothing the coverage (default: 1, no smoothing).
#'
#' @return A ggplot object.
#' @import ggplot2
#' @import rtracklayer
#' @import SummarizedExperiment
#' @export
plot_peak_track <- function(se, peak_id, group_var = "Condition", extend_bp = 1000, smooth_window = 1) {
  
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
    
    # 仅读取感兴趣的区域以提高速度
    cvg <- rtracklayer::import(bw_path, which = plot_region, as = "NumericList")[[1]]
    
    # 构造绘图数据框
    df <- data.frame(
      Position = seq(start(plot_region), end(plot_region)),
      Score = as.numeric(cvg),
      Sample = sample_id,
      Group = group
    )
    
    # Optional rolling-mean smoothing; backfill boundary NAs with original values
    if (smooth_window > 1) {
      smoothed <- stats::filter(df$Score, rep(1 / smooth_window, smooth_window), sides = 2)
      na_idx <- which(is.na(smoothed))
      smoothed[na_idx] <- df$Score[na_idx]
      df$Score <- as.numeric(smoothed)
    }
    return(df)
  })
  
  plot_df <- do.call(rbind, track_list)
  
  # 4. Render coverage track
  p <- ggplot(plot_df, aes(x = Position, y = Score, fill = Group)) +
    geom_area(alpha = 0.8) +
    facet_wrap(~Sample, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 14) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = sprintf("Coverage Track: %s", peak_id),
         subtitle = sprintf("Region: %s:%d-%d", seqnames(region), start(region), end(region)),
         x = "Genomic Coordinates", y = "Signal Intensity (BigWig)") +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
    
  return(p)
}


#' Plot 4D Epigenetic Portrait for a Specific Peak
#'
#' @description Visualizes the multi-dimensional conformation of a specific peak across experimental
#' conditions using a radar/spider chart approach. All 4 dimensions are Min-Max scaled (0 to 1)
#' for visual comparability.
#'
#' @param se A SummarizedExperiment object from build_portrait_matrix().
#' @param peak_id Character. The ID of the peak to plot (rownames of se).
#' @param group_var Character. The grouping variable in colData(se).
#' @return A ggplot object.
#' @import ggplot2
#' @import SummarizedExperiment
#' @export
plot_peak_portrait <- function(se, peak_id, group_var = "Condition") {

  if (!peak_id %in% rownames(se)) stop("Peak ID not found in the SummarizedExperiment object.")

  meta <- as.data.frame(colData(se))
  if (!group_var %in% colnames(meta)) stop(sprintf("Group variable '%s' not found in colData.", group_var))
  groups <- meta[[group_var]]

  # Optimisation 1: vapply + tapply for fast group aggregation, avoiding slow data.frame ops
  features <- c("Intensity", "Height", "Width", "Skewness")
  group_levels <- unique(groups)

  # Extract mean vector per group; produce Group × Feature matrix directly
  agg_mat <- vapply(features, function(f) {
    vals <- assay(se, f)[peak_id, ]
    tapply(vals, groups, mean, na.rm = TRUE)
  }, numeric(length(group_levels)))

  # Optimisation 2: column-wise Min-Max scaling to [0,1]
  scaled_mat <- apply(agg_mat, 2, function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[1] == rng[2]) return(rep(0.5, length(x))) # guard against zero-variance NaN
    (x - rng[1]) / (rng[2] - rng[1])
  })

  # Melt to long format for ggplot
  plot_df <- expand.grid(Group = group_levels, Feature = features, stringsAsFactors = FALSE)
  plot_df$ScaledValue <- as.vector(scaled_mat)

  # 确保特征绘制顺序的固定性
  plot_df$Feature <- factor(plot_df$Feature, levels = features)

  # 绘制极具极客美学的高级雷达图
  p <- ggplot(plot_df, aes(x = Feature, y = ScaledValue, group = Group, color = Group, fill = Group)) +
    geom_polygon(alpha = 0.25, size = 1.2) +
    geom_point(size = 3.5, shape = 21, stroke = 1.5, fill = "white") +
    coord_polar() +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(face = "bold", color = "black", size = 13),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    ) +
    labs(title = sprintf("epiPortrait: %s", peak_id),
         subtitle = "4D Conformational Shift (Min-Max Scaled)") +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1")

  return(p)
}


#' Plot Multi-dimensional Portrait Volcano (Global Landscape)
#'
#' @description A next-generation volcano plot visualizing 4 dimensions simultaneously:
#' X-axis (Intensity LogFC), Y-axis (-log10 P-value), Point Size (Absolute Width of Domain),
#' and Point Color (Skewness shift).
#'
#' @param res_df A \code{data.frame} generated by \code{test_differential_feature(se, feature = "Intensity")}, 
#' containing 'logFC' and 'adj.P.Val' columns.
#' @param se The original \code{SummarizedExperiment} object from \code{build_portrait_matrix()}.
#' @param group_var Character. The grouping variable in \code{colData(se)}.
#' @param target_group Character. The treatment group (numerator for FC).
#' @param ref_group Character. The control group (denominator for FC).
#' @param p_cutoff Numeric. adj.P.Val cutoff for significance labeling (default: 0.05).
#' @param fc_cutoff Numeric. Log2FC cutoff for significance labeling (default: 1).
#' @return A ggplot object.
#' @import ggplot2
#' @import ggrepel
#' @import SummarizedExperiment
#' @export
plot_portrait_volcano <- function(res_df, se, group_var = "Condition", target_group, ref_group,
                                  p_cutoff = 0.05, fc_cutoff = 1) {

  meta <- as.data.frame(colData(se))
  idx_target <- which(meta[[group_var]] == target_group)
  idx_ref <- which(meta[[group_var]] == ref_group)

  if (length(idx_target) == 0 || length(idx_ref) == 0) {
    stop("Target or Reference group not found in colData.")
  }

  # Width is constant across samples — map as absolute physical size (log10-compressed)
  # rather than computing an artificial logFC
  w_mat <- assay(se, "Width")
  s_mat <- assay(se, "Skewness")

  # Ensure res_df rownames match se rownames (supports subsets)
  peak_ids <- rownames(res_df)
  if (is.null(peak_ids)) {
    stop("res_df must have rownames matching the peak IDs in the SummarizedExperiment object.")
  }
  if (!any(peak_ids %in% rownames(se))) {
    stop("None of the res_df rownames match the peak IDs in the SE object.")
  }
  if (!all(peak_ids %in% rownames(se))) {
    warning(sprintf("%d row(s) in res_df not found in SE object and will be excluded.",
                    sum(!peak_ids %in% rownames(se))))
    res_df <- res_df[peak_ids %in% rownames(se), , drop = FALSE]
    peak_ids <- rownames(res_df)
  }

  # 计算辅助维度
  abs_width <- log10(w_mat[peak_ids, 1] + 1) # 取第一列即可，因为宽度固定
  skewness_shift <- rowMeans(s_mat[peak_ids, idx_target, drop=FALSE], na.rm = TRUE) -
    rowMeans(s_mat[peak_ids, idx_ref, drop=FALSE], na.rm = TRUE)

  # 合并绘图数据
  plot_df <- res_df
  plot_df$Log10_Width <- abs_width
  plot_df$Skewness_Shift <- skewness_shift
  plot_df$LogP <- -log10(plot_df$adj.P.Val)

  # 显著性分级
  plot_df$Significance <- "Not Sig"
  plot_df$Significance[plot_df$adj.P.Val < p_cutoff & plot_df$logFC > fc_cutoff] <- "Up"
  plot_df$Significance[plot_df$adj.P.Val < p_cutoff & plot_df$logFC < -fc_cutoff] <- "Down"
  plot_df$Significance <- factor(plot_df$Significance, levels = c("Up", "Down", "Not Sig"))

  # 构建极其精美的多维火山图
  p <- ggplot(plot_df, aes(x = logFC, y = LogP)) +
    geom_point(aes(size = Log10_Width, color = Skewness_Shift, alpha = Significance)) +
    scale_alpha_manual(values = c("Up" = 0.9, "Down" = 0.9, "Not Sig" = 0.15)) +
    scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0,
                          name = "Skewness\nShift (Target-Ref)") +
    scale_size_continuous(range = c(0.5, 6), name = "Domain Size\n(Log10 bp)") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey30", size = 0.5) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey30", size = 0.5) +
    theme_classic(base_size = 14) +
    labs(x = "Intensity Log2 Fold Change",
         y = "-Log10(adj.P.Val)",
         title = "epiPortrait Global Landscape",
         subtitle = sprintf("Conformational shift: %s vs %s", target_group, ref_group)) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40"),
      legend.position = "right",
      legend.background = element_rect(fill = "transparent", color = NA)
    )

  # Smart non-overlapping gene labels (if SYMBOL column exists)
  if ("SYMBOL" %in% colnames(plot_df)) {
    top_genes <- subset(plot_df, Significance != "Not Sig")
    # Order by P-value; label at most the top 15 genes
    if (nrow(top_genes) > 15) top_genes <- top_genes[order(top_genes$adj.P.Val)[1:15], ]

    if (nrow(top_genes) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = SYMBOL),
        size = 4,
        max.overlaps = 30,
        box.padding = 0.8,
        point.padding = 0.3,
        fontface = "bold.italic",
        color = "black",
        segment.color = "grey50"
      )
    }
  }

  return(p)
}


#' Plot Group-Level Coverage Overlay for Shape Shift Validation
#'
#' @description Overlays mean coverage profiles of two experimental groups on
#' the same panel with translucent ribbons (mean \eqn{\pm} SEM). This is the
#' most intuitive way to visually confirm a Shape Shift event: instead of
#' faceting by sample, the two group averages are directly superimposed so
#' differences in peak geometry are immediately apparent.
#'
#' @param se A \code{SummarizedExperiment} object from \code{build_portrait_matrix()}.
#' @param peak_id Character. The peak ID to plot (must match \code{rownames(se)}).
#' @param group_var Character. Column in \code{colData(se)} with the grouping variable.
#' @param group1 Character. Identifier for the first group (e.g., reference).
#' @param group2 Character. Identifier for the second group (e.g., treatment).
#' @param extend_bp Numeric. Base pairs to extend the view beyond the peak boundaries.
#' @param se_ribbon Logical. If \code{TRUE}, plot mean \eqn{\pm} SEM ribbons.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @import rtracklayer
#' @import SummarizedExperiment
#' @export
plot_shift_comparison <- function(se, peak_id, group_var = "Condition",
                                   group1, group2,
                                   extend_bp = 1000,
                                   se_ribbon = TRUE) {

  if (!peak_id %in% rownames(se)) stop("Peak ID not found.")
  meta <- as.data.frame(colData(se))
  if (!"bw_path" %in% colnames(meta)) stop("bw_path column missing in colData(se).")
  if (!group_var %in% colnames(meta)) stop(sprintf("'%s' not found in colData.", group_var))

  region <- rowRanges(se)[peak_id]
  plot_region <- region
  start(plot_region) <- max(1, start(region) - extend_bp)
  end(plot_region) <- end(region) + extend_bp

  .load_group_track <- function(idx_vec, group_label) {
    track_parts <- lapply(idx_vec, function(i) {
      bw_path <- meta$bw_path[i]
      cvg <- rtracklayer::import(bw_path, which = plot_region,
                                 as = "NumericList")[[1]]
      data.frame(
        Position = seq(start(plot_region), end(plot_region)),
        Score    = as.numeric(cvg),
        Replicate = meta$SampleID[i]
      )
    })
    df <- do.call(rbind, track_parts)

    pos_levels <- unique(df$Position)
    agg <- vapply(pos_levels, function(p) {
      vals <- df$Score[df$Position == p]
      c(mean = mean(vals, na.rm = TRUE),
        sem  = if (length(vals) > 1) stats::sd(vals, na.rm = TRUE) / sqrt(length(vals)) else 0)
    }, numeric(2))

    data.frame(
      Position = pos_levels,
      Mean = agg["mean", ],
      SEM  = agg["sem", ],
      Group = group_label
    )
  }

  idx1 <- which(meta[[group_var]] == group1)
  idx2 <- which(meta[[group_var]] == group2)
  if (length(idx1) == 0) stop(sprintf("Group '%s' has no samples.", group1))
  if (length(idx2) == 0) stop(sprintf("Group '%s' has no samples.", group2))

  message(sprintf("Loading BigWig signals for %s (%s vs %s)...",
                  peak_id, group1, group2))
  df1 <- .load_group_track(idx1, group1)
  df2 <- .load_group_track(idx2, group2)
  plot_df <- rbind(df1, df2)

  region_label <- sprintf("%s:%d-%d",
                          as.character(seqnames(region)),
                          start(region), end(region))

  p <- ggplot(plot_df, aes(x = Position, y = Mean, fill = Group, color = Group)) +
    geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                alpha = 0.25, color = NA, show.legend = !se_ribbon) +
    geom_line(linewidth = 1.2) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("Shape Shift Validation: %s", peak_id),
      subtitle = region_label,
      x = "Genomic Coordinates",
      y = "Signal Intensity (BigWig)",
      fill = group_var, color = group_var
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  return(p)
}


#' Manhattan Plot of 4D Conformational Outliers
#'
#' @description Plots the genomic distribution of conformational outlier scores
#' from \code{detect_conformational_outliers()}. Each point is a genomic domain;
#' the y-axis shows -log10(adj.P.Val) or Mahalanobis D². Chromosomes alternate
#' in color, and significant outliers (adj.P.Val < fdr_cutoff) are highlighted
#' with larger, more opaque points.
#'
#' @param outlier_res A \code{data.frame} from \code{detect_conformational_outliers()}.
#' @param fdr_cutoff Numeric. FDR threshold for highlighting (default: 0.05).
#' @param y_axis Character. \code{"logP"} for -log10(P-value) or \code{"D2"}
#'   for Mahalanobis distance. Default \code{"logP"}.
#' @param label_n Integer. Number of top outliers to label (default: 10).
#'   Set to 0 to suppress labels.
#' @param label_col Character or NULL. Optional column name for labeling
#'   (e.g., \code{"SYMBOL"} if gene-annotated). Uses Peak_ID if NULL or
#'   if the named column is absent.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats qchisq setNames
#' @export
plot_outlier_manhattan <- function(outlier_res,
                                    fdr_cutoff = 0.05,
                                    y_axis = c("logP", "D2"),
                                    label_n = 10,
                                    label_col = NULL) {

  y_axis <- match.arg(y_axis)

  req_cols <- c("seqnames", "start", "Peak_ID", "Mahalanobis_D2", "adj.P.Val")
  missing <- setdiff(req_cols, colnames(outlier_res))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns in outlier_res: %s",
                 paste(missing, collapse = ", ")))
  }

  dd <- outlier_res
  dd$chr <- as.character(dd$seqnames)

  if (y_axis == "logP") {
    dd$Y <- -log10(dd$adj.P.Val)
    y_label <- expression(-log[10](adj.P.Value))
    thresh_y <- -log10(fdr_cutoff)
  } else {
    dd$Y <- dd$Mahalanobis_D2
    y_label <- expression(Mahalanobis~D^2)
    thresh_y <- stats::qchisq(1 - fdr_cutoff, df = 4)
  }

  dd <- dd[!is.na(dd$Y) & is.finite(dd$Y), , drop = FALSE]
  if (nrow(dd) == 0) stop("No finite y-values to plot.")

  # ---- Cumulative genomic coordinates ------------------------------------
  chr_levels <- unique(dd$chr)
  chr_nums <- suppressWarnings(as.integer(chr_levels))
  chr_levels <- if (any(is.na(chr_nums))) sort(chr_levels) else chr_levels[order(chr_nums)]

  chr_offsets <- list()
  cumulative <- 0
  for (ch in chr_levels) {
    chr_offsets[[ch]] <- cumulative
    dd_ch <- dd[dd$chr == ch, , drop = FALSE]
    cumulative <- cumulative + max(dd_ch$start, na.rm = TRUE)
  }

  dd$GenomicPos <- dd$start + vapply(dd$chr, function(ch) chr_offsets[[ch]], numeric(1))

  chr_midpoints <- vapply(chr_levels, function(ch) {
    rows <- dd[dd$chr == ch, , drop = FALSE]
    if (nrow(rows) == 0) return(chr_offsets[[ch]])
    chr_offsets[[ch]] + median(rows$start, na.rm = TRUE)
  }, numeric(1))

  dd$chr_color <- factor((match(dd$chr, chr_levels) - 1) %% 2)
  dd$Is_Sig <- !is.na(dd$adj.P.Val) & dd$adj.P.Val < fdr_cutoff

  # ---- Labels -------------------------------------------------------------
  label_rows <- NULL
  has_gene_label <- !is.null(label_col) && label_col %in% colnames(dd)
  if (label_n > 0) {
    sig_rows <- dd[dd$Is_Sig, , drop = FALSE]
    if (nrow(sig_rows) > 0) {
      if (has_gene_label) {
        sig_rows <- sig_rows[!is.na(sig_rows[[label_col]]) &
                             sig_rows[[label_col]] != "", , drop = FALSE]
      }
      n_label <- min(label_n, nrow(sig_rows))
      if (n_label > 0) {
        label_rows <- sig_rows[order(sig_rows$adj.P.Val), ][seq_len(n_label), ]
      }
    }
  }

  # ---- Render --------------------------------------------------------------
  p <- ggplot(dd, aes(x = .data$GenomicPos, y = .data$Y)) +
    geom_point(
      aes(color = .data$chr_color, alpha = .data$Is_Sig, size = .data$Is_Sig)
    ) +
    scale_color_manual(
      values = stats::setNames(c("#4A5568", "#8899AA"), c("0", "1")),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c("FALSE" = 0.25, "TRUE" = 0.85),
      guide = "none"
    ) +
    scale_size_manual(
      values = c("FALSE" = 0.7, "TRUE" = 1.6),
      guide = "none"
    ) +
    geom_hline(
      yintercept = thresh_y, linetype = "dashed",
      color = "#E87722", linewidth = 0.6
    ) +
    scale_x_continuous(
      breaks = chr_midpoints,
      labels = gsub("^chr", "", chr_levels),
      expand = c(0.02, 0.02)
    ) +
    scale_y_continuous(expand = c(0, 0.08)) +
    theme_classic(base_size = 14) +
    labs(
      title = "Conformational Outlier Landscape",
      subtitle = sprintf("%d domains tested, %d significant (FDR < %.2f)",
                         nrow(dd), sum(dd$Is_Sig), fdr_cutoff),
      x = "Chromosome",
      y = y_label
    ) +
    theme(
      axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "none"
    )

  # ---- Gene labels ---------------------------------------------------------
  if (!is.null(label_rows) && nrow(label_rows) > 0) {
    label_text <- if (has_gene_label) label_rows[[label_col]] else label_rows$Peak_ID
    p <- p + ggrepel::geom_text_repel(
      data = label_rows,
      aes(x = .data$GenomicPos, y = .data$Y, label = label_text),
      size = 3.5,
      max.overlaps = 30,
      box.padding = 0.6,
      point.padding = 0.3,
      fontface = "bold.italic",
      color = "#1B3A5C",
      segment.color = "grey50",
      segment.size = 0.3
    )
  }

  return(p)
}
