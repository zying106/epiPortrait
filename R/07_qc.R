#' Exploratory PCA Plot for epiPortrait Features
#'
#' @description Performs Principal Component Analysis (PCA) on a specific
#'   quantitative domain feature
#' (e.g., Intensity or SignalDispersion) to visualize sample clustering and batch effects.
#'
#' @param se A SummarizedExperiment object.
#' @param feature Character. The assay to use for PCA (default: "Intensity").
#' @param group_var Character. The grouping variable for coloring points.
#' @param n_top Numeric. Use only the top N most variable peaks (default: 500).
#' @param features Character vector. Additional assays to include in a
#'   combined-feature PCA (e.g. \code{c("Intensity", "SignalDispersion")}). The
#'   top-variable domains of each feature are concatenated into the PCA
#'   variable set. Default NULL uses \code{feature} only.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom stats prcomp var
#' @examples
#' data(example_se)
#' plot_portrait_pca(example_se, feature = "Intensity", group_var = "Condition")
#' @export
plot_portrait_pca <- function(se, feature = "Intensity", group_var = "Condition",
                              n_top = 500, features = NULL) {
  # For combined-feature PCA, features = c("Intensity",
  # "SignalDispersion") concatenates the top-variable domains of each feature
  # into one variable set before prcomp.
  if (is.null(features)) {
    features <- feature
  } else {
    features <- c(feature, setdiff(features, feature))
  }
  features <- vapply(features, .epi_resolve_feature, character(1), se = se)

  meta <- as.data.frame(colData(se))
  # Single-sample PCA is undefined: fail with a clear message before any
  # variance selection (a 1-column matrix would otherwise crash inside var()).
  if (ncol(se) < 2) {
    stop("PCA requires at least 2 samples. This object has only ", ncol(se),
         " column(s); a single sample cannot define a PC space.")
  }
  mats <- lapply(features, function(f) {
    m <- assay(se, f)
    # log10(1+x) stabilises heavy-tailed signal for PCA
    m <- log10(pmax(m, 0) + 1)
    rv <- apply(m, 1, stats::var)
    keep <- is.finite(rv)
    select <- order(rv[keep], decreasing = TRUE)[seq_len(min(n_top, sum(keep)))]
    nm <- rownames(m)[keep][select]
    m[nm, , drop = FALSE]
  })
  # row-bind with unique composite rownames to avoid collisions
  combined <- do.call(rbind, lapply(seq_along(mats), function(i) {
    rownames(mats[[i]]) <- paste0(features[i], "::", rownames(mats[[i]]))
    mats[[i]]
  }))
  # drop any remaining non-finite rows
  good <- stats::complete.cases(combined)
  combined <- combined[good, , drop = FALSE]
  if (nrow(combined) < 2) stop("Too few complete variables for PCA.")
  if (ncol(combined) < 2) {
    stop("PCA requires at least 2 samples. This object has only 1 column ",
         "(a single sample cannot define a PC space).")
  }

  pca <- stats::prcomp(t(combined), scale. = TRUE)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  plot_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Group = meta[[group_var]],
    SampleID = meta$SampleID
  )

  p <- ggplot(plot_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(aes(label = SampleID), size = 3,
                             max.overlaps = 30, box.padding = 0.4) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    .epi_theme_publication(base_size = 11) +
    scale_color_manual(values = .epi_group_palette(meta[[group_var]])) +
    labs(title = sprintf("PCA: %s", paste(features, collapse = " + "))) +
    .epi_wrap_legend("colour", n_items = length(unique(meta[[group_var]])))

  return(p)
}

#' Sample Correlation Heatmap for epiPortrait Features
#'
#' @description Computes and visualizes the Spearman (default, robust to
#' heavy-tailed epigenomic signal) or Pearson correlation between samples.
#'
#' @param se A SummarizedExperiment object.
#' @param feature Character. The assay to use (default: "Intensity").
#' @param method Character. "spearman" (default) or "pearson".
#' @param limits Numeric of length 2 or NULL. Colour-scale limits for the
#'   correlation. Default \code{NULL} follows the observed data range so that
#'   within/between-group differences remain visible even when all samples are
#'   highly correlated (typical for biological replicates); when the observed
#'   range does not span 0, the colour midpoint is the range centre (see
#'   \code{dynamic_mid}). Pass an explicit fixed scale (e.g. \code{c(-1, 1)})
#'   to keep a statistically anchored, cross-dataset-comparable scale with
#'   0 = unrelated as the semantic midpoint.
#' @param dynamic_mid Logical. When \code{limits = NULL}, whether to centre the
#'   colour scale on the observed data midpoint (TRUE, default) instead of the
#'   fixed 0 anchor. Centring on the data maximises visual contrast between
#'   samples; centring on 0 keeps the "unrelated" semantics but compresses the
#'   palette when all correlations are high.
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom stats cor
#' @examples
#' data(example_se)
#' plot_portrait_correlation(example_se, feature = "Intensity")
#' plot_portrait_correlation(example_se, feature = "Intensity", limits = c(-1, 1))
#' @export
plot_portrait_correlation <- function(se, feature = "Intensity",
                                      method = c("spearman", "pearson"),
                                      limits = NULL,
                                      dynamic_mid = TRUE) {
  method <- match.arg(method)
  mat <- assay(se, feature)
  cor_mat <- stats::cor(mat, method = method)

  mid <- 0
  if (is.null(limits)) {
    # Dynamic scale: follow the observed data range so high-correlation samples
    # (replicates, same-treatment tumours) still show differences. Clamp the
    # range to at most [-1, 1].
    limits <- range(cor_mat, na.rm = TRUE)
    limits[1] <- max(-1, min(limits[1], 1))
    limits[2] <- max(-1, min(limits[2], 1))
    if (dynamic_mid && limits[1] >= 0 && limits[2] > 0) {
      # Range does not span 0: centre the palette on the data midpoint so the
      # colour gradient spans the full observed range instead of being squashed
      # against the high end.
      mid <- (limits[1] + limits[2]) / 2
    } else if (limits[1] < 0 && limits[2] > 0) {
      mid <- 0
    } else {
      mid <- (limits[1] + limits[2]) / 2
    }
  }
  if (length(limits) != 2 || !is.numeric(limits) ||
      any(!is.finite(limits)) || limits[1] >= limits[2]) {
    stop("limits must be NULL or a numeric vector c(low, high) with low < high.")
  }
  if (limits[1] < -1 || limits[2] > 1) {
    warning("limits outside [-1, 1] are allowed only for the dynamic (NULL) ",
            "scale; fixed limits should stay within the correlation range.",
            call. = FALSE)
  }

  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Sample1", "Sample2", "Correlation")
  cor_df$Sample1 <- factor(cor_df$Sample1, levels = colnames(mat))
  cor_df$Sample2 <- factor(cor_df$Sample2, levels = rev(colnames(mat)))

  p <- ggplot(cor_df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3.5) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = mid, limits = limits, oob = scales::squish
    ) +
    coord_fixed() +
    theme_classic(base_size = 11, base_family = "sans") +
    labs(
      title = sprintf("Sample Correlation: %s", feature),
      x = "", y = ""
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  return(p)
}
