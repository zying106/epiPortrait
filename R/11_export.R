#' Export Shape-Shifted Regions to IGV Session
#'
#' @description Generates an IGV session XML with BigWig tracks for all samples
#' and a BED annotation track of shape-shifted regions colored by Shape_Class.
#' Also writes the accompanying BED file. Double-click the XML to open
#' everything in IGV for direct visual validation.
#'
#' @param shift_res A data.frame from \code{classify_shape_shift()} containing
#'   classified shape-shift results with genomic coordinates.
#' @param se A \code{SummarizedExperiment} object from \code{build_portrait_matrix()}.
#' @param genome Character. Reference genome (e.g., \code{"hg38"}, \code{"hg19"},
#'   \code{"mm10"}).
#' @param session_file Character. Path for the output IGV session XML.
#' @param bed_file Character. Path for the output annotation BED. Defaults to
#'   the session file stem with \code{_annotation.bed} suffix.
#' @param fdr_cutoff Numeric. Only include peaks with adj.P.Val below this
#'   threshold (default: 0.05).
#' @param locus Character or NULL. Initial view locus (e.g., \code{"chr1:1000-2000"}).
#'   If \code{NULL} (default), zooms to the top shifted peak.
#' @param label_var Character. Column in \code{colData(se)} used for track group
#'   labels (default: \code{"Condition"}).
#'
#' @return Invisibly, the path to the session XML file.
#' @export
export_igv_session <- function(shift_res, se, genome = "hg38",
                                session_file = "epiPortrait_session.xml",
                                bed_file = NULL,
                                fdr_cutoff = 0.05,
                                locus = NULL,
                                label_var = "Condition") {

  required_cols <- c("seqnames", "start", "end", "Peak_ID",
                     "adj.P.Val", "Shape_Shift_Score", "Shape_Class")
  missing_cols <- setdiff(required_cols, colnames(shift_res))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in shift_res: %s.",
         paste(missing_cols, collapse = ", ")))
  }

  if (is.null(bed_file)) {
    bed_file <- sub("\\.xml$", "_annotation.bed", session_file)
  }

  # ---- 1. Write annotated BED file ------------------------------------------
  sig <- shift_res[shift_res$adj.P.Val <= fdr_cutoff, , drop = FALSE]
  if (nrow(sig) == 0) {
    warning("No peaks pass the FDR cutoff. No session file written.")
    return(invisible(NULL))
  }

  # Color palette for shape classes
  class_colors <- c(
    Concentration  = "205,92,92",
    Flattening     = "70,130,180",
    Polarity_Shift = "218,165,32",
    `Global Gain`  = "60,179,113",
    `Global Loss`  = "147,112,219",
    Complex        = "128,128,128",
    Stable         = "192,192,192"
  )

  bed_lines <- vapply(seq_len(nrow(sig)), function(i) {
    cl <- sig$Shape_Class[i]
    color <- if (cl %in% names(class_colors)) class_colors[cl] else "128,128,128"
    paste(
      sig$seqnames[i],
      as.integer(sig$start[i]) - 1,
      as.integer(sig$end[i]),
      sprintf("%s|%s|SC=%.2f|FDR=%.2e",
              sig$Peak_ID[i], cl, sig$Shape_Shift_Score[i], sig$adj.P.Val[i]),
      round(sig$Shape_Shift_Score[i] * 1000),
      ".",
      as.integer(sig$start[i]) - 1,
      as.integer(sig$end[i]),
      color,
      sep = "\t"
    )
  }, character(1))

  writeLines(bed_lines, bed_file)
  message(sprintf("Wrote %d annotated regions to %s", nrow(sig), bed_file))

  # ---- 2. Collect BigWig paths from SE colData ------------------------------
  meta <- as.data.frame(colData(se))
  if (!"bw_path" %in% colnames(meta)) {
    stop("colData(se) must contain a 'bw_path' column.")
  }
  # Resolve to absolute paths
  bw_paths <- normalizePath(meta$bw_path, mustWork = FALSE)
  if (!label_var %in% colnames(meta)) {
    warning(sprintf("Label variable '%s' not found in colData(se). Using SampleID only.", label_var))
    bw_labels <- meta$SampleID
  } else {
    bw_labels <- paste0(meta$SampleID, " (", meta[[label_var]], ")")
  }

  # ---- 3. Determine initial locus -------------------------------------------
  if (is.null(locus)) {
    top <- sig[which.min(sig$adj.P.Val), ]
    locus <- sprintf("%s:%d-%d", top$seqnames[1],
                     as.integer(top$start[1]), as.integer(top$end[1]))
  }

  # ---- 4. Write IGV session XML ---------------------------------------------
  xml_lines <- c(
    '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
    sprintf('<Session genome="%s" hasGeneTrack="true" hasSequenceTrack="true"',
            genome),
    sprintf('  locus="%s" version="8">', locus),
    '  <Resources>'
  )
  for (i in seq_along(bw_paths)) {
    xml_lines <- c(xml_lines,
      sprintf('    <Resource path="%s" label="%s"/>', bw_paths[i], bw_labels[i]))
  }
  xml_lines <- c(xml_lines,
    sprintf('    <Resource path="%s" label="epiPortrait Shape-Shifted Regions"/>',
            normalizePath(bed_file, mustWork = FALSE)),
    '  </Resources>',
    '  <Panel name="DataPanel">'
  )
  # BigWig tracks
  for (i in seq_along(bw_paths)) {
    xml_lines <- c(xml_lines,
      sprintf('    <Track altColor="175,175,175" autoScale="true" color="0,0,178"',
              bw_labels[i]),
      sprintf('           displayMode="COLLAPSED" featureVisibilityWindow="-1"',
              bw_labels[i]),
      sprintf('           fontSize="10" id="%s" name="%s"',
              bw_labels[i], bw_labels[i]),
      sprintf('           showDataRange="true" visible="true"/>'))
  }
  # Annotation track
  xml_lines <- c(xml_lines,
    '    <Track altColor="175,175,175" autoScale="false" color="0,0,178"',
    '           displayMode="EXPANDED" featureVisibilityWindow="-1"',
    '           fontSize="10" id="epiPortrait_annotations"',
    '           name="epiPortrait Shape-Shifted Regions"',
    '           renderer="BASIC_FEATURE" showDataRange="true" visible="true"/>',
    '  </Panel>',
    '  <PanelLayout>',
    '    <Panel name="DataPanel" height="600"/>',
    '  </PanelLayout>',
    '</Session>'
  )

  writeLines(xml_lines, session_file)
  message(sprintf("Wrote IGV session to %s (locus: %s)", session_file, locus))

  invisible(session_file)
}
