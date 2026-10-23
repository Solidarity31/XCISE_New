#' Add EXCISE calls to a SingleCellExperiment
#'
#' Matches the cell barcodes of an \code{\link{ExciseResult-class}} to the
#' columns of a \code{SummarizedExperiment} (for example a
#' \code{SingleCellExperiment}) and adds the calls to its \code{colData}.
#'
#' @param x A \code{SummarizedExperiment} or \code{SingleCellExperiment}.
#' @param result An \code{\link{ExciseResult-class}} object.
#' @param barcodeColumn Name of the \code{colData} column holding the cell
#'     barcodes. If \code{NULL}, \code{colnames(x)} are used.
#' @param prefix Prefix of the added columns.
#'
#' @return \code{x} with three new \code{colData} columns:
#'     \code{<prefix>call} (factor; \code{NA} for cells without EXCISE
#'     evidence), \code{<prefix>hap1} and \code{<prefix>hap2} (UMIs supporting
#'     haplotype X1 and X2).
#'
#' @examples
#' bam <- system.file("extdata", "sim_chrX.bam", package = "EXCISE")
#' vcf <- system.file("extdata", "sim_chrX_hets.vcf", package = "EXCISE")
#' res <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 1)
#' bc <- xciCalls(res)$barcode
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = matrix(0L, 2, length(bc),
#'         dimnames = list(c("g1", "g2"), bc))))
#' sce <- addXciCalls(sce, res)
#' table(sce$xci_call, useNA = "ifany")
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @export
addXciCalls <- function(x, result, barcodeColumn = NULL, prefix = "xci_") {
    stopifnot(is(x, "SummarizedExperiment"), is(result, "ExciseResult"))
    cd <- colData(x)
    bc <- if (is.null(barcodeColumn)) colnames(x) else cd[[barcodeColumn]]
    if (is.null(bc))
        stop("no barcodes: set colnames(x) or give 'barcodeColumn'")
    calls <- xciCalls(result)
    m <- match(as.character(bc), calls$barcode)
    if (all(is.na(m)))
        warning("no barcode of 'x' was found in the EXCISE result")
    cd[[paste0(prefix, "call")]] <- calls$call[m]
    cd[[paste0(prefix, "hap1")]] <- calls$hap1[m]
    cd[[paste0(prefix, "hap2")]] <- calls$hap2[m]
    colData(x) <- cd
    x
}
