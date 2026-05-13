#' Example 4D Epigenetic Portrait Dataset
#'
#' @description A simulated \code{SummarizedExperiment} containing 500 genomic
#' domains with 4D geometric features (Width, Intensity, Height, Skewness)
#' measured across 6 samples (3 Control, 3 Treatment). The dataset includes
#' 5 types of ground-truth shape shifts embedded by simulation, enabling
#' method validation against a known answer.
#'
#' @format A \code{SummarizedExperiment} object with:
#' \describe{
#'   \item{assays}{4 matrices (Width, Intensity, Height, Skewness), each 500 x 6}
#'   \item{rowRanges}{GRanges with mcols including \code{True_Class} (factor with
#'     levels: Concentration, Flattening, Polarity_Shift, Global_Gain, Global_Loss)}
#'   \item{colData}{SampleID, Condition (Control/Treatment), bw_path}
#' }
#'
#' @usage data("example_se")
#'
#' @examples
#' data(example_se)
#' example_se
#' table(SummarizedExperiment::mcols(
#'   SummarizedExperiment::rowRanges(example_se))$True_Class)
"example_se"
