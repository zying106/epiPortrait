#' Get Mark-Specific Preset
#'
#' @param mark "H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3", etc.
#' @return List of preset parameters, including:
#'   \itemize{
#'     \item \code{label}: human-readable domain label
#'     \item \code{stitch_distance}: default stitching gap in bp. 0 / NULL is a
#'           SENTINEL meaning "do not automatically re-stitch at the
#'           preset-workflow level" — treat it as "skip \code{stitch_epi_peaks()}"
#'           (e.g. broad repressive marks such as H3K27me3 / H3K9me3, whose
#'           domains are already spatially clustered by the upstream broad-domain
#'           caller, and unknown/generic marks). Note that
#'           \code{stitch_epi_peaks(x, 0)} itself still merges
#'           overlapping/adjacent ranges; the sentinel must be handled by the
#'           caller (see the vignette preset workflow).
#'     \item \code{primary_feature}: ranking feature
#'     \item \code{mark_class}: "active" (H3K27ac, H3K4me3, H3K4me1) or
#'           "broad_repressive" (H3K27me3, H3K9me3)
#'     \item \code{taxonomy_style}: "super_domain" (active marks) or
#'           "repressive_remodeling" (broad repressive marks). This controls
#'           DISPLAY terminology only; the internal canonical classes
#'           (Intensity-Super / Breadth-Super / Dual-Super) are unchanged.
#'     \item \code{class_labels}: display aliases for the canonical taxonomy.
#'   }
#' @examples
#' get_mark_preset("H3K27ac")
#' get_mark_preset("H3K27me3")
#' @export
get_mark_preset <- function(mark) {
  active <- list(mark_class = "active", taxonomy_style = "super_domain",
                 class_labels = c("Intensity-Super", "Breadth-Super",
                                  "Dual-Super", "Typical", "Uncertain"))
  repressive <- list(mark_class = "broad_repressive",
                     taxonomy_style = "repressive_remodeling",
                     class_labels = c("Intensity-Extreme", "Extended-Domain",
                                      "Dual-Extreme", "Typical", "Uncertain"))
  p <- switch(mark,
    H3K27ac = c(list(label = "active super-enhancer-like domain",
                     stitch_distance = 12500L,
                     exclude_promoter = TRUE,
                     primary_feature = "Intensity",
                     normalize = "None"), active),
    H3K4me3 = c(list(label = "broad promoter-associated domain",
                     stitch_distance = 5000L,
                     exclude_promoter = FALSE,
                     primary_feature = "Breadth",
                     normalize = "None"), active),
    H3K27me3 = c(list(label = "Polycomb-repressive broad domain",
                      stitch_distance = 0L,
                      exclude_promoter = FALSE,
                      primary_feature = "Breadth",
                      normalize = "None"), repressive),
    H3K9me3 = c(list(label = "heterochromatic broad domain",
                     stitch_distance = 0L,
                     exclude_promoter = FALSE,
                     primary_feature = "Breadth",
                     normalize = "None"), repressive),
    H3K4me1 = c(list(label = "enhancer-associated domain",
                     stitch_distance = 12500L,
                     exclude_promoter = TRUE,
                     primary_feature = "Intensity",
                     normalize = "None"), active),
    ATAC = c(list(label = "accessible chromatin domain",
                  stitch_distance = 10000L,
                  exclude_promoter = FALSE,
                  primary_feature = "Intensity",
                  normalize = "None"), active),
    {
      # Unknown / unregistered mark: do NOT silently treat it as active.
      # Broad-repressive marks that are not yet preset (e.g. H4K20me3) would
      # otherwise receive the wrong taxonomy and stitching defaults.
      # stitch_distance = 0L = "do not automatically re-stitch at the
      # preset-workflow level" (a safe sentinel for unknown marks); users who
      # know their mark should stitch can call stitch_epi_peaks() explicitly.
      warning(sprintf(
        "Unknown mark '%s': using a generic preset (no mark-specific stitching ",
        mark, "or terminology). Register the mark via get_mark_preset() if you ",
        "need mark-aware behaviour."), call. = FALSE)
      c(list(label = "generic epigenomic domain",
             stitch_distance = 0L,
             exclude_promoter = FALSE,
             primary_feature = "Intensity",
             normalize = "None"),
        list(mark_class = "generic", taxonomy_style = "generic",
             class_labels = c("Intensity-Super", "Breadth-Super",
                              "Dual-Super", "Typical", "Uncertain")))
    }
  )
  p$mark <- mark
  p
}


# Map a canonical epiPortrait class to the mark-appropriate DISPLAY label.
# The internal class values are never changed; only the surface terminology
# adapts (active marks keep Intensity-Super / Breadth-Super / Dual-Super;
# broad repressive marks use Intensity-Extreme / Extended-Domain / Dual-Extreme).
.epi_class_display <- function(class_value, preset = NULL) {
  if (is.null(preset)) {
    map <- c("Intensity-Super" = "Intensity-Super",
             "Breadth-Super"   = "Breadth-Super",
             "Dual-Super"      = "Dual-Super",
             "Typical"         = "Typical",
             "Uncertain"       = "Uncertain")
  } else {
    map <- stats::setNames(preset$class_labels,
                           c("Intensity-Super", "Breadth-Super", "Dual-Super",
                             "Typical", "Uncertain"))
  }
  out <- class_value
  for (nm in names(map)) out[!is.na(out) & out == nm] <- unname(map[nm])
  out
}


# Internal: resolve assay name with backward compatibility.
# Canonical assay names (v1.0): Intensity, SignalDispersion,
# NativeMaxPeakWidth, NativeOccupiedWidth, NativePeakCount. Legacy aliases
# (TotalIntensity, Width) are mapped to the canonical ones.
.resolve_assay <- function(se, name) {
  avail <- assayNames(se)
  avail_row <- colnames(rowData(se))
  if (name %in% avail || name %in% avail_row) return(name)
  map <- c("TotalIntensity" = "Intensity", "Intensity" = "TotalIntensity",
           "IntervalWidth" = "Width", "Width" = "IntervalWidth")
  alt <- map[name]
  if (!is.na(alt) && (alt %in% avail || alt %in% avail_row)) return(alt)
  name
}
