#' Exploratory PCA Plot for epiPortrait Features
#'
#' @description Performs Principal Component Analysis (PCA) on a specific geometric feature
#' (e.g., Intensity or Height) to visualize sample clustering and batch effects.
#'
#' @param se A SummarizedExperiment object.
#' @param feature Character. The assay to use for PCA (default: "Intensity").
#' @param group_var Character. The grouping variable for coloring points.
#' @param n_top Numeric. Use only the top N most variable peaks (default: 500).
#'
#' @return A ggplot object.
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom stats prcomp var
#' @export
plot_portrait_pca <- function(se, feature = "Intensity", group_var = "Condition", n_top = 500) {
  
  mat <- assay(se, feature)
  meta <- as.data.frame(colData(se))
  
  # Select peaks with highest variance
  rv <- apply(mat, 1, stats::var)
  select <- order(rv, decreasing = TRUE)[seq_len(min(n_top, nrow(mat)))]
  
  pca <- stats::prcomp(t(mat[select, ]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  
  plot_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Group = meta[[group_var]],
    SampleID = meta$SampleID
  )
  
  p <- ggplot(plot_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 4) +
    geom_text(aes(label = SampleID), vjust = 1.5, size = 3) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_bw(base_size = 14) +
    scale_color_brewer(palette = "Set1") +
    labs(title = sprintf("PCA Plot: %s", feature))
    
  return(p)
}

#' Sample Correlation Heatmap for epiPortrait Features
#'
#' @description Computes and visualizes the Pearson correlation between samples based
#' on their geometric portraits. Returns a ggplot object consistent with the rest of the
#' epiPortrait visualization suite.
#'
#' @param se A SummarizedExperiment object.
#' @param feature Character. The assay to use (default: "Intensity").
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom stats cor
#' @export
plot_portrait_correlation <- function(se, feature = "Intensity") {

  mat <- assay(se, feature)
  cor_mat <- stats::cor(mat, method = "pearson")

  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Sample1", "Sample2", "Correlation")
  cor_df$Sample1 <- factor(cor_df$Sample1, levels = colnames(mat))
  cor_df$Sample2 <- factor(cor_df$Sample2, levels = rev(colnames(mat)))

  p <- ggplot(cor_df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3.5) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1)
    ) +
    coord_fixed() +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("Sample Correlation: %s", feature),
      x = "", y = ""
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  return(p)
}
