# EXCISE 0.99.0

* Initial Bioconductor submission.
* `runExcise()`: C++ port of `excise1.5.pl` (htslib via Rhtslib), returning an
  `ExciseResult` with per-cell calls (`DataFrame`), phased SNVs (`GRanges`)
  and a run summary; optionally writes the five `excise1.5.pl` output files.
* `addXciCalls()`: adds the calls to a `SingleCellExperiment`.
