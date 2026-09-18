#' epiPortrait: Replicate-Aware Epigenomic Domain Profiling
#'
#' @description
#' epiPortrait provides replicate-aware quantitative phenotyping of epigenomic
#' domains across biological samples and conditions. Each shared genomic domain
#' is characterized by three measurement layers: **Intensity** (integrated
#' signal magnitude), **native peak breadth geometry** (per-replicate native
#' peak width), and **SignalDispersion** (continuous within-domain signal
#' architecture). Intensity and Breadth are the canonical classification axes;
#' super-domain states (Intensity-Super, Breadth-Super, Dual-Super) are called
#' with replicate-level evidence and explicit support rules. An orthogonal
#' Breadth evidence layer distinguishes operational peak-call absence from a
#' technical or assignment no-call, while SignalDispersion is retained as a
#' secondary architecture descriptor. The package supports condition
#' transitions, continuous domain-width remodeling,
#' signal quality control, genome-aware domain annotation, optional BEDPE
#' contact evidence, expression-aware candidate-gene prioritization,
#' domain-aware functional interpretation (GO ORA / GSEA against a
#' domain-derived testable universe), and result visualization. Primary use
#' case: histone-mark ChIP-seq / CUT&Tag; continuous-track workflows also apply
#' to ATAC-seq-derived domains.
#'
#' @section Intensity ranking scale:
#' For heavy-tailed active-mark Intensity such as H3K27ac, the default
#' \code{log10(x + 1)} ranking can place the inflection cutoff too low and
#' over-call Intensity-Super. Ranking on the raw scale
#' (\code{log_transform = FALSE}) gives Intensity-Super counts and
#' ROSE-style concordance much closer to a stitched-enhancer (ROSE) reference;
#' \code{method = "elbow"} and \code{"tangent"} are near-equivalent on the raw
#' scale, and native-peak Breadth (for example H3K4me3 broadPeak) is insensitive
#' to the cutoff method.
#'
#' @section Core workflow:
#' \enumerate{
#'   \item \code{\link{get_consensus_peaks}} — Build consensus domain universe
#'   \item \code{\link{stitch_epi_peaks}} — Stitch proximal peaks into macro-domains
#'   \item \code{\link{build_portrait_matrix}} — Extract portrait matrix from BigWigs
#'   \item \code{\link{normalize_portrait}} — Normalize Intensity while preserving SignalDispersion
#'   \item \code{\link{call_super_domains}} — Intensity / peak-level Breadth calling
#'   \item \code{\link{get_breadth_evidence}} — Audit peak presence versus no-call
#'   \item \code{\link{annotate_epi_domains}} — Domain-aware annotation & candidate genes
#' }
#'
#' @section Visualization:
#' \itemize{
#'   \item \code{\link{plot_portrait_pca}} — PCA of quantitative domain features
#'   \item \code{\link{plot_portrait_correlation}} — Sample correlation heatmap
#'   \item \code{\link{plot_peak_track}} — Raw BigWig coverage track
#'   \item \code{\link{plot_hockey_stick}} — Super-domain ranking plot
#' }
#'
#' @section Supplementary:
#' \itemize{
#'   \item \code{\link{combine_superdomain_calls}} — Combine Intensity/Breadth super calls
#'   \item \code{\link{compare_superdomains}} — Condition-aware super-domain transitions
#'   \item \code{\link{filter_promoter_peaks}} — Exclude promoter-proximal peaks
#'   \item \code{\link{filter_blacklist}} — Filter blacklisted regions
#' }
#'
#' @section Functional interpretation:
#' \itemize{
#'   \item \code{\link{get_domain_gene_universe}} — Domain-derived testable gene universe
#'   \item \code{\link{enrich_epi_domains}} — ORA / GSEA on domain phenotypes and remodeling
#'   \item \code{\link{rank_epi_genes}} — Continuous domain score to signed gene ranking
#'   \item \code{\link{compare_epi_enrichment}} — Shared-universe comparison across groups
#'   \item \code{\link{plot_epi_enrichment}} — Comparative effect-size heatmap
#' }
#'
#' @docType package
#' @name epiPortrait
#' @importFrom methods is
#' @importFrom utils head
#' @examples
#' data(example_se)
#' se <- call_super_domains(example_se, feature = "Intensity", verbose = FALSE)
#' table(SummarizedExperiment::rowData(se)$Intensity_Domain_Type)
"_PACKAGE"

utils::globalVariables(c(
  "Rank", "Value", "Type", "Type_Label", "SYMBOL",
  "Feature", "ScaledValue", "Group",
  "Score", "Sample1", "Sample2", "Correlation",
  "PC1", "PC2", "SampleID",
  "Mean", "SEM",
  "Intensity", "Breadth", "Class", "Label",
  "From", "To", "Level", "Status", "N",
  "Occ_Start", "Occ_End",
  "Cause", "Freq", "Prop", "LabelY", "Frac", "FillKey",
  "logFC", "negLog10P",
  "group", "term_name", "effect", "label"
))
