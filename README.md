# XCISE_New

**XCISE_New**, referred to in this repository as **EXCISE** (*An Enhanced XCISE pipeline*), is an extended pipeline for inferring **X chromosome inactivation (XCI)** states from allele specific single cell sequencing data. The project is designed for workflows that use **heterozygous SNVs**, **allele specific read or UMI evidence**, and **cell barcodes** to assign each cell to an XCI related state.

Importantly, **our work is established on the original XCISE framework**. XCISE_New keeps the same overall biological goal and core logic of allele specific XCI inference, while adding practical improvements in computation, initialization, tuning flexibility, and downstream evaluation.

## Built on the original XCISE

The starting point of this work is the original **XCISE** pipeline. XCISE_New is not meant to replace the conceptual foundation of XCISE, but to extend it. In particular, this repository preserves the central idea that XCI status can be inferred from informative allele specific signals on the X chromosome at the single cell level.

On top of that foundation, XCISE_New introduces a more flexible and research friendly implementation that is easier to tune, benchmark, and compare across datasets.

## What XCISE_New aims to do

XCISE_New is developed to:

1. extract informative allele specific signal from BAM and VCF based inputs,
2. infer cell level XCI states from heterozygous X chromosomal SNVs,
3. improve computational efficiency relative to earlier workflows,
4. support more stable initialization and search behavior during inference, and
5. facilitate downstream comparison, annotation, and validation.

In practice, the pipeline is intended to help classify cells into biologically meaningful categories such as cells favoring one X haplotype, the other haplotype, mixed or ambiguous states, and low coverage or uninformative cells.

## Main improvements in this repository

Compared with the original XCISE style workflow, this repository adds several practical extensions.

### 1. Improved computational efficiency

The repository includes support for parallel execution and more flexible runtime control. This is especially useful when processing larger BAM files or running repeated optimization or search procedures.

### 2. More flexible initialization

The pipeline includes optional **pretraining or seeding style controls**, which are intended to stabilize inference and improve starting configurations before deeper optimization.

### 3. Better reproducibility

The repository exposes options for **random seed control**, making results easier to reproduce and compare across runs.

### 4. More tunable inference behavior

The enhanced scripts expose more parameters for controlling search depth, penalties, and block level updates. This makes the method more suitable for methodological development and benchmarking.

### 5. Stronger downstream utilities

This repository also includes helper scripts for:

- comparing two `*_bc2xci.txt` outputs,
- annotating VCF results with gene annotation information, and
- supporting downstream result interpretation.

## Repository structure

The current repository contains the following main components:

- `excise.pl`  
  Main enhanced pipeline script.

- `xcise_mt.pl`  
  Multi try or multi process version of the XCI inference workflow.

- `xcise_mt.R`  
  R side support for the inference workflow.

- `compare_bc2xci.R`  
  Utility script to compare two `*_bc2xci.txt` outputs.

- `annotate_vcf_with_gtf.R` and `annotation.R`  
  Scripts for annotation and downstream interpretation.

- `Comparison between two outputs.R`  
  Additional comparison or benchmarking helper script.

## Input data

XCISE_New is designed around the same type of core inputs used by XCISE style workflows:

1. **BAM files** containing aligned reads with cell barcode and, when available, UMI information.
2. **VCF files** containing heterozygous SNVs, especially on the X chromosome.
3. Optional metadata or annotation resources for downstream interpretation.

Because XCI inference depends on informative allele specific evidence, input quality remains critical. In particular, mapping quality, barcode quality, UMI evidence, and the number of informative SNVs per cell can all affect final cell assignments.

## Output interpretation

A main output of the pipeline is the `*_bc2xci.txt` style result table, which records barcode level summaries and cell assignments. Depending on the evidence available for each cell, assignments may include classes such as:

- `X1`
- `X2`
- `Both`
- `Low_coverage`
- `Unknown`

These outputs can then be compared across runs, parameter settings, or preprocessing strategies using the comparison scripts included in the repository.

## Example usage

A typical run follows the same high level style as the original XCISE workflow, while allowing additional enhanced options:

```bash
perl excise.pl \
  -o sample_name \
  -s heterozygous_snvs.vcf \
  -b input.bam \
  -r X \
  -j 4 \
  -samthreads 4 \
  -seed 1 \
  -pretrain
```

For output comparison, the repository also includes:

```bash
Rscript compare_bc2xci.R file1_bc2xci.txt file2_bc2xci.txt
```

## Why this project matters

X chromosome inactivation is a central biological process in many developmental and disease related settings. Reliable cell level inference of XCI status can help characterize heterogeneity across cells, study escape from inactivation, and improve the interpretation of X linked regulatory behavior.

By building directly on the original XCISE framework and extending it with improved efficiency, tunability, and downstream support, XCISE_New aims to provide a stronger platform for both methodological development and applied biological analysis.

## Acknowledgement

**XCISE_New is established on the original XCISE framework**, and this repository should be understood as an enhanced continuation of that line of work rather than a disconnected method. We acknowledge the conceptual and methodological foundation provided by XCISE.
