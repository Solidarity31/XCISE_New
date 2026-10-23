# EXCISE R/Bioconductor package

`EXCISE/` is the R package for Bioconductor submission. It runs the C++ port of
`excise1.5.pl` (`cpp/excise_core.cpp`), linked to htslib through Rhtslib.

```r
library(EXCISE)
bam <- system.file("extdata", "sim_chrX.bam", package = "EXCISE")
vcf <- system.file("extdata", "sim_chrX_hets.vcf", package = "EXCISE")
res <- runExcise(bam, vcf, tries = 100, pretrain = TRUE, workers = 4, seed = 1)
xciCalls(res)      # per-cell calls (DataFrame)
xciSnvs(res)       # phased SNVs (GRanges)
sce <- addXciCalls(sce, res)   # add to a SingleCellExperiment
```

`runExcise(..., outPrefix = "sample")` also writes the five `excise1.5.pl`
output files. With the same seed and settings these are byte-identical to
`cpp/excise15` (tested on two benchmark datasets, 2,000 and 16,000 cells).

## Development

```bash
make -C cpp sync-r               # copy the C++ core into EXCISE/src after changing it
make -C cpp check-sync           # verify the two copies are identical
Rscript -e 'Rcpp::compileAttributes("bioc/EXCISE"); roxygen2::roxygenise("bioc/EXCISE")'
rm -f bioc/EXCISE/src/*.o bioc/EXCISE/src/*.so   # roxygen leaves -O0 objects behind
R CMD build bioc/EXCISE && R CMD check EXCISE_0.99.0.tar.gz
Rscript -e 'BiocCheck::BiocCheck("EXCISE_0.99.0.tar.gz")'
```

Checks in this environment (R 4.3.3, Bioconductor 3.18; no internet, so
BiocCheck's online checks were skipped):

- `R CMD check`: 0 errors, 0 warnings, 2 notes. Installed size is 6.9 MB because
  htslib is linked statically. The other note is GNU make, which every package that
  links Rhtslib needs.
- `BiocCheck`: 0 errors, 0 warnings, 3 notes: maintainer ORCID, licence, and
  indentation in roxygen-generated `\usage` sections.

## Before submitting

1. **Licence (blocking).** Neither EXCISE nor XCISE has a licence yet.
   Bioconductor requires an open-source licence. The code derives from XCISE, so
   agree a licence with the XCISE authors, then replace `EXCISE/LICENSE` and the
   `License:` field (e.g. `Artistic-2.0`, `GPL-3`, `MIT + file LICENSE`).
2. **Authors.** The XCISE authors are listed as contributors (`ctb`) because the
   code derives from `xcise.pl`. Confirm the roles and full names with them. Add
   your ORCID: `comment = c(ORCID = "...")`.
3. **Citation.** Add `inst/CITATION` once the XCISE (and EXCISE) paper titles
   are final.
4. **Own repository.** Bioconductor needs the package at the top level of its
   own Git repository. Copy `bioc/EXCISE/` into a new repository, e.g.
   `Solidarity31/EXCISE-bioc`, or make this repository package-only.
5. **Current Bioconductor.** Re-run `R CMD check` and `BiocCheck` with internet
   access, under the current release and devel R/Bioconductor.
6. **Submit.** Open an issue at https://github.com/Bioconductor/Contributions
   and follow the reviewer checklist.
