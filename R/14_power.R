#' Simulate Statistical Power for Shape Shift Detection
#'
#' @description Estimates the statistical power of epiPortrait's MANOVA-based
#' shape shift test across a range of sample sizes and effect sizes. Generates
#' synthetic 4D portrait data with known ground-truth effects and computes the
#' proportion of peaks detected at a given alpha. Use this before designing an
#' experiment to determine the minimum number of replicates needed.
#'
#' @param n_range Integer vector. Sample sizes per group to test (e.g., 3:10).
#' @param eta2_range Numeric vector. Target partial eta-squared values to test
#'   (e.g., \code{c(0.1, 0.2, 0.4, 0.6)}). These are Cohen-style benchmarks:
#'   0.01 = small, 0.06 = medium, 0.14 = large (for multivariate context).
#' @param n_dims Integer. Number of dimensions carrying the true effect
#'   (1 = only Intensity shifts, 4 = all dimensions shift). Default: 2.
#' @param n_sim Integer. Number of simulation replicates per condition
#'   (default: 200). Higher values give smoother power curves.
#' @param alpha Numeric. Significance threshold (default: 0.05).
#' @param workers Integer. Parallel workers for simulation (default: 1).
#'
#' @return A data.frame with columns \code{n}, \code{eta2}, \code{n_dims},
#'   \code{power}, \code{power_se} suitable for plotting power curves.
#'
#' @import BiocParallel
#' @importFrom stats rnorm manova
#' @export
simulate_power <- function(n_range = 3:10,
                            eta2_range = c(0.06, 0.14, 0.30, 0.50),
                            n_dims = 2,
                            n_sim = 200,
                            alpha = 0.05,
                            workers = 1) {

  if (workers > 1) {
    param <- if (.Platform$OS.type == "unix") BiocParallel::MulticoreParam(workers)
             else BiocParallel::SnowParam(workers)
  } else {
    param <- BiocParallel::SerialParam()
  }

  # Build grid of simulation conditions
  grid <- expand.grid(
    n     = n_range,
    eta2  = eta2_range,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$eta2, grid$n), ]

  message(sprintf("Simulating power for %d conditions x %d replicates = %d MANOVA runs...",
                  nrow(grid), n_sim, nrow(grid) * n_sim))

  # Pre-compute the mean-shift vector for each eta2 value
  # For k dimensions with equal effect: eta2_pillai ~ (d^2 * k) / (d^2 * k + 4*k)
  # Simplified calibration: d = sqrt(4 * eta2 / (1 - eta2))
  shift_cache <- vapply(eta2_range, function(e2) {
    d <- sqrt(4 * e2 / max(1 - e2, 1e-8))
    # Build a k=4 mean vector where only n_dims dimensions carry the effect
    shift <- rep(0, 4)
    if (n_dims > 0) shift[1:n_dims] <- d / sqrt(n_dims)
    shift
  }, numeric(4))
  colnames(shift_cache) <- as.character(eta2_range)

  res_list <- BiocParallel::bplapply(seq_len(nrow(grid)), function(idx) {

    n     <- grid$n[idx]
    eta2  <- grid$eta2[idx]
    shift <- shift_cache[, as.character(eta2)]

    sig_count <- 0
    for (sim_i in seq_len(n_sim)) {
      # Generate control: 4D MVN(0, I)
      ctrl <- matrix(rnorm(n * 4, mean = 0, sd = 1), nrow = n, ncol = 4)

      # Generate treatment: 4D MVN(shift, I)
      treat <- matrix(rnorm(n * 4, mean = shift, sd = 1), nrow = n, ncol = 4)

      Y <- rbind(ctrl, treat)
      colnames(Y) <- c("W", "I", "H", "S")
      groups <- factor(rep(c("Control", "Treatment"), each = n))

      fit <- tryCatch(
        stats::manova(Y ~ groups),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      res <- tryCatch(
        summary(fit, test = "Pillai"),
        error = function(e) NULL
      )
      if (is.null(res)) next

      pval <- res$stats["groups", "Pr(>F)"]
      if (!is.na(pval) && pval < alpha) sig_count <- sig_count + 1
    }

    power <- sig_count / n_sim
    power_se <- sqrt(power * (1 - power) / n_sim)

    data.frame(
      n        = n,
      eta2     = eta2,
      n_dims   = n_dims,
      power    = power,
      power_se = power_se,
      stringsAsFactors = FALSE
    )

  }, BPPARAM = param)

  result <- do.call(rbind, res_list)
  rownames(result) <- NULL

  message(sprintf("Power simulation complete. Max power: %.2f at n=%d, eta2=%.2f.",
                  max(result$power), result$n[which.max(result$power)],
                  result$eta2[which.max(result$power)]))

  return(result)
}


#' Plot Power Curve from Simulation Results
#'
#' @param power_df A data.frame returned by \code{simulate_power()}.
#'
#' @return A \code{ggplot} object showing power curves with error bars.
#' @import ggplot2
#' @export
plot_power_curve <- function(power_df) {

  power_df$eta2_label <- factor(
    sprintf("eta^2 == %.2f", power_df$eta2),
    levels = sprintf("eta^2 == %.2f", sort(unique(power_df$eta2)))
  )

  p <- ggplot(power_df, aes(x = n, y = power, color = eta2_label)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = pmax(0, power - 1.96 * power_se),
                      ymax = pmin(1, power + 1.96 * power_se)),
                  width = 0.2, alpha = 0.4) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50") +
    annotate("text", x = min(power_df$n), y = 0.82, label = "80% power",
             hjust = 0, size = 3.5, color = "grey50") +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(breaks = unique(power_df$n)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Statistical Power for 4D Shape Shift Detection",
      subtitle = sprintf("%d dimensions with true effect, MANOVA-Pillai",
                         power_df$n_dims[1]),
      x = "Sample Size per Group",
      y = "Power",
      color = expression(partial~eta^2)
    ) +
    theme_minimal(base_size = 14)

  return(p)
}
