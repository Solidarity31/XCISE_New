#' The ExciseResult class
#'
#' Results of \code{\link{runExcise}}: one X-inactivation call per cell, the
#' phased heterozygous SNVs, and a run summary.
#'
#' @slot calls A \code{\link[S4Vectors]{DataFrame}} with one row per cell
#'     barcode: \code{barcode}, \code{hap1} and \code{hap2} (UMIs supporting
#'     haplotype X1 and X2) and \code{call}, a factor with levels \code{X1},
#'     \code{X2}, \code{Both}, \code{Low_coverage} and \code{Unknown}.
#' @slot snvs A \code{\link[GenomicRanges]{GRanges}} with one range per SNV
#'     that passed the filters. Metadata columns: \code{id}, \code{ref},
#'     \code{alt}, \code{x1Allele} (\code{"Ref"}, \code{"Alt"} or \code{"Unk"}
#'     for unassigned SNVs, e.g. escape candidates), \code{refUmis},
#'     \code{altUmis}, \code{x1Umis}, \code{x2Umis} and \code{monoallelic}.
#' @slot summary A list: best score, UMI counts, discordancy, SNV counts,
#'     median minor-allele fraction, cell counts per call and whitelist
#'     statistics.
#' @slot profile \code{NULL}, or (with \code{returnProfile = TRUE}) a
#'     \code{DataFrame} with one row per cell and SNV with at least one UMI.
#' @slot parameters A list of the parameters used.
#'
#' @return An \code{ExciseResult} object (constructed by
#'     \code{\link{runExcise}}).
#'
#' @seealso \code{\link{runExcise}}, \code{\link{xciCalls}}
#' @examples
#' bam <- system.file("extdata", "sim_chrX.bam", package = "EXCISE")
#' vcf <- system.file("extdata", "sim_chrX_hets.vcf", package = "EXCISE")
#' res <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 1)
#' res
#' @aliases ExciseResult
#' @importFrom methods new is
#' @exportClass ExciseResult
setClass("ExciseResult",
    representation(
        calls = "DataFrame",
        snvs = "GRanges",
        summary = "list",
        profile = "ANY",
        parameters = "list"
    )
)

#' Accessors for ExciseResult objects
#'
#' @param x An \code{\link{ExciseResult-class}} object.
#'
#' @return \code{xciCalls}: per-cell calls (\code{DataFrame}).
#'     \code{xciSnvs}: phased SNVs (\code{GRanges}).
#'     \code{xciSummary}: run summary (list).
#'     \code{xciProfile}: per cell and SNV UMI counts (\code{DataFrame}), or
#'     \code{NULL} if \code{runExcise} was called with
#'     \code{returnProfile = FALSE}.
#'
#' @examples
#' bam <- system.file("extdata", "sim_chrX.bam", package = "EXCISE")
#' vcf <- system.file("extdata", "sim_chrX_hets.vcf", package = "EXCISE")
#' res <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 1)
#' table(xciCalls(res)$call)
#' xciSnvs(res)
#' xciSummary(res)$bestScore
#' @name ExciseResult-accessors
NULL

#' @rdname ExciseResult-accessors
#' @export
xciCalls <- function(x) {
    stopifnot(is(x, "ExciseResult"))
    x@calls
}

#' @rdname ExciseResult-accessors
#' @export
xciSnvs <- function(x) {
    stopifnot(is(x, "ExciseResult"))
    x@snvs
}

#' @rdname ExciseResult-accessors
#' @export
xciSummary <- function(x) {
    stopifnot(is(x, "ExciseResult"))
    x@summary
}

#' @rdname ExciseResult-accessors
#' @export
xciProfile <- function(x) {
    stopifnot(is(x, "ExciseResult"))
    x@profile
}

#' @rdname ExciseResult-class
#' @param object An \code{ExciseResult} object.
#' @importFrom methods show
#' @export
setMethod("show", "ExciseResult", function(object) {
    s <- object@summary
    n <- s$cellCounts
    msg <- c(
        "ExciseResult",
        sprintf(paste0("  cells: %d  (X1 %d, X2 %d, Both %d, ",
            "Low_coverage %d, Unknown %d)"),
            nrow(object@calls), n[["X1"]], n[["X2"]], n[["Both"]],
            n[["Low_coverage"]], n[["Unknown"]]),
        sprintf("  SNVs:  %d  (phased %d, unassigned %d, monoallelic %d)",
            s$nSnvs, s$nSnvs - s$nNonInformative, s$nNonInformative,
            s$nMonoallelic),
        sprintf("  best score: %s   discordancy: %s",
            format(s$bestScore),
            if (is.na(s$discordancy)) "NA" else
                formatC(s$discordancy, digits = 4, format = "f")),
        "  accessors: xciCalls(), xciSnvs(), xciSummary(), xciProfile()"
    )
    cat(msg, sep = "\n")
})
