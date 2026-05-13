#' Render epiPortrait Parameterized Report
#'
#' @description Renders a publication-ready HTML report from the epiPortrait
#' parameterized R Markdown template. Modules are automatically hidden when their
#' corresponding optional parameters are left at `NULL`.
#'
#' @param sample_sheet A `data.frame` or path to a CSV file with columns
#'   `SampleID`, `Condition`, and `bw_path`.
#' @param consensus_peaks A `GRanges` object or path to a BED file.
#' @param group_var Column name in sample_sheet for the grouping variable (default: "Condition").
#' @param target_group Treatment group identifier.
#' @param ref_group Reference/control group identifier.
#' @param output_file Path for the output HTML report.
#' @param stitch_distance Optional. Maximum bp between peaks to stitch (default: `NULL` — skip stitching).
#' @param normalization_method Normalization method: `"TotalSignal"`, `"TMM"`, `"Quantile"`, `"Z-score"`, or `"None"`.
#' @param bg_quantile Background quantile gate for robust skewness (default: 0.1).
#' @param workers Number of parallel workers (default: 1).
#' @param on_disk Use HDF5-backed on-disk storage (default: `FALSE`).
#' @param p_cutoff Adjusted P-value cutoff for significance (default: 0.05).
#' @param fc_cutoff Log2 fold-change cutoff for volcano plot (default: 1).
#' @param annotation_genome Optional genome for ChIPseeker annotation (`"hg38"`, `"hg19"`, `"mm10"`).
#'   If `NULL`, the annotation section is hidden.
#' @param super_domain_feature Optional feature for ROSE-style ranking (`"Intensity"`, `"Width"`, `"Height"`).
#'   If `NULL`, the super-domain section is hidden.
#' @param report_title Title for the report.
#' @param author Author name for the report.
#' @param ... Additional arguments passed to `rmarkdown::render()`.
#'
#' @return The path to the rendered HTML report (invisibly).
#'
#' @export
render_report <- function(sample_sheet,
                          consensus_peaks,
                          group_var = "Condition",
                          target_group,
                          ref_group,
                          output_file = "epiPortrait_report.html",
                          stitch_distance = NULL,
                          normalization_method = "TotalSignal",
                          bg_quantile = 0.1,
                          workers = 1,
                          on_disk = FALSE,
                          p_cutoff = 0.05,
                          fc_cutoff = 1,
                          annotation_genome = NULL,
                          super_domain_feature = NULL,
                          report_title = "epiPortrait Analysis Report",
                          author = "epiPortrait",
                          ...) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required to render reports. Install it with install.packages('rmarkdown').")
  }

  template_rmd <- system.file(
    "rmarkdown", "templates", "epiPortrait-report", "skeleton", "skeleton.Rmd",
    package = "epiPortrait"
  )

  if (template_rmd == "") {
    stop("Report template not found. Re-install epiPortrait or check the installation.")
  }

  output_dir <- dirname(normalizePath(output_file, mustWork = FALSE))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  params_list <- list(
    sample_sheet          = sample_sheet,
    consensus_peaks       = consensus_peaks,
    group_var             = group_var,
    target_group          = target_group,
    ref_group             = ref_group,
    stitch_distance       = stitch_distance,
    normalization_method  = normalization_method,
    bg_quantile           = bg_quantile,
    workers               = workers,
    on_disk               = on_disk,
    p_cutoff              = p_cutoff,
    fc_cutoff             = fc_cutoff,
    annotation_genome     = annotation_genome,
    super_domain_feature  = super_domain_feature,
    report_title          = report_title,
    author                = author,
    output_dir            = output_dir
  )

  rmarkdown::render(
    input = template_rmd,
    output_file = output_file,
    params = params_list,
    envir = new.env(),
    ...
  )

  invisible(output_file)
}
