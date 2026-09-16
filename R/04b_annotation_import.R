# Import standardized domain-gene annotations produced outside epiPortrait.
#
# External programs are intentionally not parsed here. Users convert their
# output to the documented long-table contract, and this module validates that
# contract before rebuilding the same three annotation views used by the native
# annotator.

.canonical_annotation_relations <- function() {
  c("nearest_tss", "promoter_overlap", "gene_body_overlap",
    "fully_contained", "bedpe_promoter_contact")
}

.annotation_link_columns <- function() {
  c("domain_id", "gene_id", "gene_symbol", "relation_type",
    "distance_to_tss_bp", "overlap_bp", "domain_overlap_fraction",
    "feature_overlap_fraction", "bedpe_record_id", "evidence_source",
    "contact_score", "annotation_source")
}

.standardize_annotation_links <- function(x, source) {
  if (!is.data.frame(x)) {
    stop("annotations must be a data.frame or a path to a tab-delimited file.",
         call. = FALSE)
  }
  required <- c("domain_id", "gene_id", "relation_type")
  missing_required <- setdiff(required, colnames(x))
  if (length(missing_required) > 0L) {
    stop("annotations is missing required column(s): ",
         paste(missing_required, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(x) == 0L) {
    stop("annotations must contain at least one domain-gene relationship.",
         call. = FALSE)
  }

  for (nm in required) {
    x[[nm]] <- as.character(x[[nm]])
    if (any(is.na(x[[nm]]) | !nzchar(x[[nm]]))) {
      stop("annotations$", nm, " must contain non-missing, non-empty values.",
           call. = FALSE)
    }
  }
  allowed <- .canonical_annotation_relations()
  invalid_relations <- setdiff(unique(x$relation_type), allowed)
  if (length(invalid_relations) > 0L) {
    stop("Unsupported relation_type value(s): ",
         paste(invalid_relations, collapse = ", "), ". Convert external labels ",
         "to one of: ", paste(allowed, collapse = ", "), ".", call. = FALSE)
  }

  defaults <- list(
    gene_symbol = NA_character_,
    distance_to_tss_bp = NA_real_,
    overlap_bp = NA_real_,
    domain_overlap_fraction = NA_real_,
    feature_overlap_fraction = NA_real_,
    bedpe_record_id = NA_character_,
    contact_score = NA_real_)
  for (nm in names(defaults)) {
    if (!nm %in% colnames(x)) x[[nm]] <- defaults[[nm]]
  }
  x$gene_symbol <- as.character(x$gene_symbol)
  x$bedpe_record_id <- as.character(x$bedpe_record_id)
  numeric_columns <- c("distance_to_tss_bp", "overlap_bp",
                       "domain_overlap_fraction", "feature_overlap_fraction",
                       "contact_score")
  for (nm in numeric_columns) {
    old <- x[[nm]]
    if (is.numeric(old)) {
      converted <- as.numeric(old)
    } else {
      old_character <- as.character(old)
      missing_value <- is.na(old) | !nzchar(old_character) |
        old_character == "NA"
      converted <- rep(NA_real_, length(old_character))
      if (any(!missing_value)) {
        parsed <- utils::type.convert(old_character[!missing_value],
                                      as.is = TRUE)
        if (!is.numeric(parsed)) {
          stop("annotations$", nm, " must be numeric or NA.",
               call. = FALSE)
        }
        converted[!missing_value] <- as.numeric(parsed)
      }
    }
    x[[nm]] <- converted
    if (any(!is.na(converted) & !is.finite(converted))) {
      stop("annotations$", nm, " must contain finite values or NA.",
           call. = FALSE)
    }
  }
  if (any(x$overlap_bp < 0, na.rm = TRUE)) {
    stop("annotations$overlap_bp must be non-negative.", call. = FALSE)
  }
  for (nm in c("domain_overlap_fraction", "feature_overlap_fraction")) {
    if (any(x[[nm]] < 0 | x[[nm]] > 1, na.rm = TRUE)) {
      stop("annotations$", nm, " must be between 0 and 1.", call. = FALSE)
    }
  }

  # evidence_source is a semantic category used internally; annotation_source
  # records the external program or workflow named by the user.
  x$evidence_source <- ifelse(
    x$relation_type == "bedpe_promoter_contact", "bedpe", "external")
  x$annotation_source <- source
  # bedpe_record_id is the contact identity used for support deduplication.
  # It must be a stable, user-supplied value: silently synthesizing one per call
  # would make repeated mode = "append" imports reuse the same ids for different
  # records, undercounting contact support. Fail loudly instead.
  is_bedpe <- x$relation_type == "bedpe_promoter_contact"
  missing_record <- is_bedpe &
    (is.na(x$bedpe_record_id) | !nzchar(x$bedpe_record_id))
  if (any(missing_record)) {
    stop("annotations$bedpe_record_id is required for every ",
         "bedpe_promoter_contact row; ", sum(missing_record),
         " row(s) are missing it. Provide a stable per-record identifier ",
         "(it is the contact support deduplication key).", call. = FALSE)
  }

  # Keep user-supplied audit columns after the stable core schema.
  core <- .annotation_link_columns()
  x <- x[, c(core, setdiff(colnames(x), core)), drop = FALSE]
  rownames(x) <- NULL
  x
}

.align_annotation_link_columns <- function(x, columns) {
  for (nm in setdiff(columns, colnames(x))) x[[nm]] <- NA
  x[, columns, drop = FALSE]
}

#' Import Standardized External Domain-Gene Annotations
#'
#' @description Imports domain-gene relationships generated outside
#'   epiPortrait after the caller has converted them to the package's standard
#'   long-table contract. The function does not parse software-specific output
#'   or guess relation types.
#'
#' @param se A SummarizedExperiment with non-empty, unique row names used as
#'   domain identifiers.
#' @param annotations A data.frame, or a path to a tab-delimited file with a
#'   header. Required columns are \code{domain_id}, \code{gene_id} and
#'   \code{relation_type}. \code{bedpe_record_id} is additionally REQUIRED for
#'   every \code{bedpe_promoter_contact} row: it is the stable per-record
#'   identity used to deduplicate contact support, so it must be supplied by the
#'   caller and is never generated. Optional standard columns are
#'   \code{gene_symbol}, \code{distance_to_tss_bp}, \code{overlap_bp},
#'   \code{domain_overlap_fraction}, \code{feature_overlap_fraction} and
#'   \code{contact_score}. Additional columns are retained for auditing.
#' @param source A non-empty label identifying the external program or workflow
#'   (for example, \code{"ABC_loop_links"} or
#'   \code{"custom_loop_workflow"}).
#' @param mode \code{"replace"} (default) replaces the active domain-gene link
#'   table. \code{"append"} adds the imported relationships to an existing
#'   annotation; use append only for evidence that should coexist with the
#'   existing links.
#' @param nearest_tss_cutoff_bp Distance below which an imported
#'   \code{nearest_tss} relationship receives proximal evidence tier 2.
#'
#' @details \code{relation_type} must already use one of epiPortrait's canonical
#'   values: \code{nearest_tss}, \code{promoter_overlap},
#'   \code{gene_body_overlap}, \code{fully_contained}, or
#'   \code{bedpe_promoter_contact}. Software-specific labels should be retained
#'   in an additional column such as \code{external_relation}.
#'
#'   Every \code{bedpe_promoter_contact} row must carry a stable
#'   \code{bedpe_record_id}. The value is the contact identity used to collapse
#'   duplicate evidence and to count supporting records
#'   (\code{bedpe_support_count}, \code{bedpe_contact_score_n}); a missing value
#'   is an error rather than an auto-generated placeholder, because generated ids
#'   restart per import and would collide across repeated \code{mode = "append"}
#'   calls, undercounting contact support.
#'
#'   Importing annotations changes the candidate-gene universe. Stored
#'   enrichment results are therefore removed, with a warning, while expression
#'   summaries are retained.
#'
#' @return The updated SummarizedExperiment with rebuilt raw, pair-level and
#'   per-domain annotation tables and import provenance.
#' @examples
#' data(example_se)
#' links <- data.frame(
#'   domain_id = rownames(example_se)[1:2],
#'   gene_id = c("101", "102"),
#'   gene_symbol = c("GENE1", "GENE2"),
#'   relation_type = rep("bedpe_promoter_contact", 2),
#'   bedpe_record_id = c("loop_1", "loop_2"),
#'   abc_score = c(0.04, 0.03)
#' )
#' se2 <- import_domain_annotations(example_se, links,
#'                                  source = "ABC_loop_links")
#' S4Vectors::metadata(se2)$annotation_summary
#' @export
import_domain_annotations <- function(
    se, annotations, source, mode = c("replace", "append"),
    nearest_tss_cutoff_bp = getOption("epiPortrait.nearest_tss_cutoff_bp",
                                      10000)) {
  mode <- match.arg(mode)
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment.", call. = FALSE)
  }
  domain_ids <- rownames(se)
  if (is.null(domain_ids) || any(is.na(domain_ids) | !nzchar(domain_ids)) ||
      anyDuplicated(domain_ids)) {
    stop("se must have non-missing, non-empty, unique row names.",
         call. = FALSE)
  }
  if (!is.character(source) || length(source) != 1L || is.na(source) ||
      !nzchar(source)) {
    stop("source must be one non-empty character string.", call. = FALSE)
  }
  if (length(nearest_tss_cutoff_bp) != 1L ||
      !is.numeric(nearest_tss_cutoff_bp) ||
      !is.finite(nearest_tss_cutoff_bp) || nearest_tss_cutoff_bp < 0) {
    stop("nearest_tss_cutoff_bp must be a finite non-negative number.",
         call. = FALSE)
  }
  if (is.character(annotations) && length(annotations) == 1L) {
    if (!file.exists(annotations)) {
      stop("Annotation file not found: ", annotations, call. = FALSE)
    }
    annotations <- utils::read.delim(
      annotations, header = TRUE, stringsAsFactors = FALSE,
      check.names = FALSE)
  }
  imported <- .standardize_annotation_links(annotations, source = source)
  unknown <- setdiff(unique(imported$domain_id), domain_ids)
  if (length(unknown) > 0L) {
    stop(length(unknown), " imported domain_id value(s) are not present in se: ",
         paste(utils::head(unknown, 5L), collapse = ", "), ".",
         call. = FALSE)
  }
  duplicated_rows <- duplicated(imported)
  n_duplicates <- sum(duplicated_rows)
  if (n_duplicates > 0L) imported <- imported[!duplicated_rows, , drop = FALSE]

  old_links <- S4Vectors::metadata(se)$domain_gene_links
  if (mode == "append" && (is.null(old_links) || nrow(old_links) == 0L)) {
    stop("mode = 'append' requires existing domain-gene links. Use mode = ",
         "'replace' for an unannotated object.", call. = FALSE)
  }
  if (mode == "append") {
    if (!"annotation_source" %in% colnames(old_links)) {
      old_links$annotation_source <- ifelse(
        old_links$evidence_source == "bedpe", "epiPortrait:BEDPE",
        "epiPortrait:native")
    }
    columns <- union(colnames(old_links), colnames(imported))
    old_links <- .align_annotation_link_columns(old_links, columns)
    imported <- .align_annotation_link_columns(imported, columns)
    links <- rbind(old_links, imported)
    links <- links[!duplicated(links), , drop = FALSE]
  } else {
    links <- imported
    S4Vectors::metadata(se)$bedpe_provenance <- NULL
  }
  rownames(links) <- NULL

  se <- .rebuild_annotation_views(
    se, links, nearest_tss_cutoff_bp = nearest_tss_cutoff_bp)
  entry <- list(
    source = source,
    mode = mode,
    imported_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    n_input_rows = nrow(annotations),
    n_imported_rows = nrow(imported),
    n_exact_duplicates_removed = n_duplicates,
    relation_counts = as.list(table(imported$relation_type)),
    nearest_tss_cutoff_bp = nearest_tss_cutoff_bp)
  history <- S4Vectors::metadata(se)$annotation_import_provenance
  if (is.null(history)) history <- list()
  S4Vectors::metadata(se)$annotation_import_provenance <- c(history, list(entry))

  provenance <- S4Vectors::metadata(se)$annotation_provenance
  if (is.null(provenance)) provenance <- list()
  provenance$annotation_mode <- paste0("external_", mode)
  provenance$external_import_sources <- unique(c(
    provenance$external_import_sources, source))
  provenance$active_gene_links <- if (mode == "replace") {
    paste0("external standardized links: ", source)
  } else {
    "native and externally imported standardized links"
  }
  provenance$nearest_tss_cutoff_bp <- nearest_tss_cutoff_bp
  S4Vectors::metadata(se)$annotation_provenance <- provenance

  invalidated <- character()
  for (nm in c("enrichment", "enrichment_comparison")) {
    if (!is.null(S4Vectors::metadata(se)[[nm]])) {
      S4Vectors::metadata(se)[[nm]] <- NULL
      invalidated <- c(invalidated, nm)
    }
  }
  if (length(invalidated) > 0L) {
    warning("Imported annotations changed the candidate-gene universe; removed ",
            "stored ", paste(invalidated, collapse = " and "), " results.",
            call. = FALSE)
  }
  se
}
