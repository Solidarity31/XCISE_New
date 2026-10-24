#' Call X-chromosome inactivation states in single cells
#'
#' Runs EXCISE 1.5 on WASP-tagged single-cell RNA-seq alignments. The
#' heterozygous SNVs on one chromosome (usually X) are phased into two
#' haplotypes, X1 and X2, by maximising the number of concordant minus a
#' penalty times the number of discordant UMIs per cell; each cell is then
#' called \code{X1}, \code{X2}, \code{Both} (biallelic), \code{Low_coverage}
#' or \code{Unknown}. The implementation is a C++ port of
#' \code{excise1.5.pl} and gives the same results for the same optimum.
#'
#' @param bam Character vector: path(s) to coordinate-sorted, indexed BAM
#'     files aligned by STAR in WASP mode (tags \code{vA}, \code{vG},
#'     \code{vW}), with cell barcodes in \code{CB} (or \code{RG}) and UMIs in
#'     \code{UB} (optional).
#' @param vcf Path to a VCF (plain or gzipped) with the heterozygous SNVs.
#' @param chromosome Chromosome to analyse, as named in the BAM and VCF.
#' @param tries Number of independent optimiser restarts; the best-scoring
#'     one is kept.
#' @param minUmis Minimum number of UMIs per SNV.
#' @param minMaf Minimum minor-allele fraction per SNV (0 to 0.5).
#' @param penalty Penalty per discordant UMI in the score (at least 1).
#' @param pretrain Seed the tries from a co-segregation graph of SNVs. This
#'     is recommended; \code{FALSE} matches the \code{excise1.5.pl} default.
#' @param jitter Per-try perturbation of the pretrained seed (0 to 1).
#' @param preCap,preMinWeight,preTieSkip Pretraining settings: maximum SNV
#'     pairs per cell, minimum absolute edge weight, and whether tied allele
#'     calls are skipped.
#' @param phase3 Try to assign SNVs left unassigned by the optimiser
#'     (Phase 3).
#' @param phase3Ratio Minimum concordance ratio for a Phase 3 assignment.
#' @param whitelist Optional path to a barcode whitelist (one barcode per
#'     line, first whitespace-delimited column; may be gzipped, e.g. a 10x
#'     \code{barcodes.tsv.gz}).
#' @param whitelistStripGem Strip a trailing \code{-1}, \code{-2}, ... from
#'     barcodes before matching the whitelist.
#' @param skipSoftClipped Ignore reads with soft-clipped bases.
#' @param monostretch Ignore reads containing a mononucleotide run of at
#'     least this length.
#' @param workers Number of threads for the optimiser tries.
#' @param bamThreads Number of threads for reading BAM files.
#' @param seed Integer seed. If \code{NULL}, a seed is drawn from R's random
#'     number generator, so \code{set.seed()} makes runs reproducible.
#' @param outPrefix If not \code{NULL}, also write the five
#'     \code{excise1.5.pl} output files
#'     \code{<outPrefix>_chr<chromosome>_XCISE*}
#'     (identical in format to the Perl script).
#' @param returnProfile Return the per cell and SNV UMI counts
#'     (\code{\link{xciProfile}}); can be large.
#' @param verbose Report progress with \code{message()}.
#'
#' @return An \code{\link{ExciseResult-class}} object.
#'
#' @details Arguments map to \code{excise1.5.pl} flags: \code{tries} (-t),
#'     \code{minUmis} (-u), \code{minMaf} (-m), \code{penalty} (-p),
#'     \code{pretrain} (-pretrain), \code{jitter} (-jitter), \code{phase3}
#'     (-ph), \code{phase3Ratio} (-ph3_ratio), \code{whitelist} (-wl),
#'     \code{whitelistStripGem} (-wl_strip_gem), \code{skipSoftClipped} (-sm),
#'     \code{monostretch} (-ms), \code{workers} (-j), \code{bamThreads}
#'     (-samthreads), \code{seed} (-seed). Out-of-range values are reset as in
#'     the Perl script (e.g. \code{minUmis <= 0} becomes 10).
#'
#'     Labels X1 and X2 are arbitrary; they are swapped so that there are at
#'     least as many X2 as X1 cells.
#'
#' @references XCISE: Henning RH, Rust TM, Dijksterhuis K, Eggen BJL,
#'     Guryev V. bioRxiv (2024). \doi{10.1101/2024.08.29.610317}
#'
#' @examples
#' bam <- system.file("extdata", "sim_chrX.bam", package = "EXCISE")
#' vcf <- system.file("extdata", "sim_chrX_hets.vcf", package = "EXCISE")
#' res <- runExcise(bam, vcf, tries = 10, pretrain = TRUE, seed = 42)
#' res
#' table(xciCalls(res)$call)
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom stats runif
#' @useDynLib EXCISE, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
runExcise <- function(bam, vcf, chromosome = "X", tries = 100L, minUmis = 4,
    minMaf = 0, penalty = 5, pretrain = FALSE, jitter = 0.25,
    preCap = 5000L, preMinWeight = 3, preTieSkip = TRUE, phase3 = FALSE,
    phase3Ratio = 0.9, whitelist = NULL, whitelistStripGem = FALSE,
    skipSoftClipped = FALSE, monostretch = 15L, workers = 1L,
    bamThreads = workers, seed = NULL, outPrefix = NULL,
    returnProfile = FALSE, verbose = FALSE) {
    .checkFiles(bam, "bam")
    .checkFiles(vcf, "vcf")
    .checkIndex(bam)
    if (!is.null(whitelist)) .checkFiles(whitelist, "whitelist")
    if (is.null(seed)) seed <- floor(runif(1L, 0, .Machine$integer.max))
    opts <- list(
        bam = normalizePath(bam), vcf = normalizePath(vcf),
        chromosome = .scalarChar(chromosome, "chromosome"),
        whitelist = if (is.null(whitelist)) "" else normalizePath(whitelist),
        whitelistStripGem = .flag(whitelistStripGem, "whitelistStripGem"),
        tries = .count(tries, "tries"),
        minUmis = .number(minUmis, "minUmis"),
        minMaf = .number(minMaf, "minMaf"),
        penalty = .number(penalty, "penalty"),
        skipSoftClipped = .flag(skipSoftClipped, "skipSoftClipped"),
        monostretch = .count(monostretch, "monostretch"),
        pretrain = .flag(pretrain, "pretrain"),
        preCap = .count(preCap, "preCap"),
        preMinWeight = .number(preMinWeight, "preMinWeight"),
        preTieSkip = .flag(preTieSkip, "preTieSkip"),
        jitter = .number(jitter, "jitter"),
        phase3 = .flag(phase3, "phase3"),
        phase3Ratio = .number(phase3Ratio, "phase3Ratio"),
        workers = .count(workers, "workers"),
        bamThreads = .count(bamThreads, "bamThreads"),
        seed = .number(seed, "seed"),
        outPrefix = if (is.null(outPrefix)) "" else
            .scalarChar(outPrefix, "outPrefix"),
        returnProfile = .flag(returnProfile, "returnProfile")
    )
    raw <- .exciseRun(opts, .flag(verbose, "verbose"))
    .asExciseResult(raw, opts)
}

.asExciseResult <- function(raw, opts) {
    levels <- c("X1", "X2", "Both", "Low_coverage", "Unknown")
    cells <- DataFrame(
        barcode = raw$cells$barcode,
        hap1 = as.integer(raw$cells$hap1),
        hap2 = as.integer(raw$cells$hap2),
        call = factor(raw$cells$call, levels = levels)
    )
    sv <- raw$snvs
    dir <- sv$dir
    x1 <- ifelse(dir == 1L, sv$refUmis, ifelse(dir == -1L, sv$altUmis, 0))
    x2 <- ifelse(dir == -1L, sv$refUmis, ifelse(dir == 1L, sv$altUmis, 0))
    snvs <- GRanges(seqnames = sv$chrom,
        ranges = IRanges(start = sv$pos, width = 1L))
    mcols(snvs) <- DataFrame(
        id = sv$id, ref = sv$ref, alt = sv$alt,
        x1Allele = c("Alt", "Unk", "Ref")[dir + 2L],
        refUmis = as.integer(sv$refUmis), altUmis = as.integer(sv$altUmis),
        x1Umis = as.integer(x1), x2Umis = as.integer(x2),
        monoallelic = as.logical(sv$mono)
    )
    s <- raw$summary
    tot <- s$concordant + s$discordant
    s$discordancy <- if (tot > 0) s$discordant / tot else NA_real_
    profile <- NULL
    if (isTRUE(opts$returnProfile)) {
        p <- raw$profile
        profile <- DataFrame(
            barcode = cells$barcode[p$cell + 1L],
            pos = sv$pos[p$snv + 1L],
            refUmis = p$refUmis, altUmis = p$altUmis,
            x1Allele = mcols(snvs)$x1Allele[p$snv + 1L],
            call = cells$call[p$cell + 1L]
        )
    }
    new("ExciseResult", calls = cells, snvs = snvs, summary = s,
        profile = profile, parameters = opts)
}

.checkFiles <- function(x, what) {
    if (!is.character(x) || length(x) == 0L || anyNA(x))
        stop("'", what, "' must be a character vector of file paths")
    missing <- x[!file.exists(x)]
    if (length(missing))
        stop(what, " file not found: ", paste(missing, collapse = ", "))
    invisible(TRUE)
}

.checkIndex <- function(bam) {
    idx <- vapply(bam, function(b) {
        any(file.exists(c(paste0(b, ".bai"), paste0(b, ".csi"),
            sub("\\.bam$", ".bai", b))))
    }, logical(1))
    if (!all(idx))
        stop("BAM index (.bai or .csi) not found for: ",
            paste(bam[!idx], collapse = ", "),
            "; create it with Rsamtools::indexBam() or 'samtools index'")
    invisible(TRUE)
}

.scalarChar <- function(x, what) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
        stop("'", what, "' must be a single non-empty string")
    x
}

.flag <- function(x, what) {
    if (!is.logical(x) || length(x) != 1L || is.na(x))
        stop("'", what, "' must be TRUE or FALSE")
    x
}

.number <- function(x, what) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x))
        stop("'", what, "' must be a single number")
    as.numeric(x)
}

.count <- function(x, what) {
    x <- .number(x, what)
    if (x != round(x)) stop("'", what, "' must be a whole number")
    as.integer(x)
}
