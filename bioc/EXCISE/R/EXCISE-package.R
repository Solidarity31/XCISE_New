#' EXCISE: X-chromosome inactivation calling from single-cell RNA-seq
#'
#' EXCISE infers which X chromosome is inactive in each cell from
#' allele-specific single-cell RNA-seq reads (STAR WASP tags) at heterozygous
#' SNVs. It is an enhanced version of XCISE (Henning et al., 2024). The main
#' function is \code{\link{runExcise}}; results are returned as an
#' \code{\link{ExciseResult-class}} object and can be added to a
#' \code{SingleCellExperiment} with \code{\link{addXciCalls}}.
#'
#' @keywords internal
"_PACKAGE"
