bam <- system.file("extdata", "sim_chrX.bam", package = "EXCISE")
vcf <- system.file("extdata", "sim_chrX_hets.vcf", package = "EXCISE")
truthFile <- system.file("extdata", "sim_chrX_truth.tsv", package = "EXCISE")

test_that("runExcise returns a well-formed ExciseResult", {
    res <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 1)
    expect_s4_class(res, "ExciseResult")
    calls <- xciCalls(res)
    expect_s4_class(calls, "DataFrame")
    expect_identical(nrow(calls), 300L)
    expect_identical(levels(calls$call),
        c("X1", "X2", "Both", "Low_coverage", "Unknown"))
    expect_false(is.unsorted(calls$barcode))
    snvs <- xciSnvs(res)
    expect_s4_class(snvs, "GRanges")
    expect_true(all(mcols(snvs)$x1Allele %in% c("Ref", "Alt", "Unk")))
    s <- xciSummary(res)
    expect_identical(s$nSnvs, length(snvs))
    expect_equal(sum(s$cellCounts), 300)
    expect_true(s$cellCounts[["X2"]] >= s$cellCounts[["X1"]])
    expect_null(xciProfile(res))
    expect_output(show(res), "ExciseResult")
})

test_that("results are reproducible", {
    a <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 7,
        workers = 2)
    b <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 7,
        workers = 1)
    expect_identical(xciCalls(a), xciCalls(b))
    set.seed(3)
    c1 <- runExcise(bam, vcf, tries = 4)
    set.seed(3)
    c2 <- runExcise(bam, vcf, tries = 4)
    expect_identical(xciCalls(c1), xciCalls(c2))
})

test_that("calls agree with the simulated truth", {
    res <- runExcise(bam, vcf, tries = 20, pretrain = TRUE, seed = 42)
    truth <- read.delim(truthFile, header = FALSE,
        col.names = c("barcode", "state"))
    calls <- xciCalls(res)
    pred <- as.character(calls$call)[match(truth$barcode, calls$barcode)]
    swapped <- c(X1 = "X2", X2 = "X1", Both = "Both")[pred]
    swapped[is.na(swapped)] <- pred[is.na(swapped)]
    acc <- max(mean(pred == truth$state), mean(swapped == truth$state))
    expect_gt(acc, 0.9)
})

test_that("outPrefix writes the excise1.5.pl files", {
    out <- file.path(tempdir(), "excise_test")
    res <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 1,
        outPrefix = out)
    files <- paste0(out, "_chrX_XCISE", c("_bc2xci.txt", "_summary.txt",
        ".vcf", "_haplotypes.tsv", "_cell_snv_profile.tsv"))
    expect_true(all(file.exists(files)))
    b <- read.delim(files[1], header = FALSE,
        col.names = c("barcode", "hap1", "hap2", "call"))
    calls <- xciCalls(res)
    expect_identical(b$barcode, calls$barcode)
    expect_identical(b$call, as.character(calls$call))
    expect_identical(as.integer(b$hap1), calls$hap1)
})

test_that("returnProfile gives per cell x SNV counts", {
    res <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 1,
        returnProfile = TRUE)
    p <- xciProfile(res)
    expect_s4_class(p, "DataFrame")
    expect_true(all(p$barcode %in% xciCalls(res)$barcode))
    expect_true(all(p$refUmis + p$altUmis > 0))
})

test_that("bad inputs give informative errors", {
    expect_error(runExcise("missing.bam", vcf), "not found")
    expect_error(runExcise(bam, "missing.vcf"), "not found")
    expect_error(runExcise(bam, vcf, tries = 2.5), "whole number")
    expect_error(runExcise(bam, vcf, pretrain = "yes"), "TRUE or FALSE")
    noidx <- file.path(tempdir(), "noindex.bam")
    file.copy(bam, noidx, overwrite = TRUE)
    expect_error(runExcise(noidx, vcf), "index")
})

test_that("addXciCalls adds calls to colData", {
    res <- runExcise(bam, vcf, tries = 4, pretrain = TRUE, seed = 1)
    bc <- c(xciCalls(res)$barcode[1:5], "NOT_A_CELL")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(0L, 2, 6,
            dimnames = list(c("g1", "g2"), bc))))
    se <- addXciCalls(se, res)
    expect_identical(as.character(se$xci_call[1:5]),
        as.character(xciCalls(res)$call[1:5]))
    expect_true(is.na(se$xci_call[6]))
})
