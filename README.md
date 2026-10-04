# EXCISE — *An Enhanced XCISE Pipeline*

<p align="center">
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"/></a>
  <img src="https://img.shields.io/badge/language-Perl%20%7C%20R-blue.svg" alt="Language"/>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg" alt="Platform"/>
  <img src="https://img.shields.io/badge/input-10x%20Genomics%20scRNA--seq-green.svg" alt="Input"/>
  <img src="https://img.shields.io/badge/WASP-bias%20corrected-orange.svg" alt="WASP"/>
</p>

<p align="center">
  <b>Cell-level X chromosome inactivation inference from allele-specific single-cell RNA-seq data</b><br/>
  Built on the original <a href="https://github.com/Vityay/XCISE">XCISE</a> framework · WASP-corrected · MCMC-optimized · Parallelized
</p>

---

## Overview

**EXCISE** extends the [XCISE](https://github.com/Vityay/XCISE) framework for inferring X chromosome inactivation (XCI) states at single-cell resolution. Starting from allele-specific BAM files and heterozygous SNV calls on the X chromosome, EXCISE assigns each cell barcode to one of five XCI categories: `X1`, `X2`, `Both`, `Low_coverage`, or `Unknown`.

XCISE (*XCI calling from Single cell Expression data*) was developed by Henning, Rust, Dijksterhuis, Eggen & Guryev ([bioRxiv 2024](https://doi.org/10.1101/2024.08.29.610317); code: [Vityay/XCISE](https://github.com/Vityay/XCISE)). `xcise_mt.pl` and `xcise_mcmc.pl` are direct modifications of the original `xcise.pl`, and `excise.pl` is derived from them.

Key advances over the original XCISE include **MCMC-based haplotype optimization** with simulated annealing, **observation-based co-segregation seeding**, **parallel multi-try execution**, and a suite of downstream utilities for annotation and inter-run comparison.

---

## Features

| Feature | XCISE (original) | EXCISE |
|---|:---:|:---:|
| WASP mapping-bias correction | ✅ | ✅ |
| Greedy haplotype optimization | ✅ | ✅ |
| MCMC + simulated annealing | ❌ | ✅ |
| Observation-based co-segregation seed | ❌ | ✅ |
| Multi-SNV block proposals | ❌ | ✅ |
| Parallel tries (`-j`) | ❌ | ✅ |
| Unsupervised pretraining seed | ❌ | ✅ |
| Co-segregation graph export (TSV + DOT) | ❌ | ✅ |
| Reproducible RNG (`-seed`) | ❌ | ✅ |
| VCF annotation with GTF | ❌ | ✅ |
| Inter-run comparison + Cohen's κ | ❌ | ✅ |

---

## Repository Structure

```
EXCISE/
├── excise.pl                    # Main pipeline (MCMC + pretraining)
├── xcise_mt.pl                  # Greedy multi-try variant
├── xcise_mcmc.pl                # MCMC + co-segregation graph variant
├── compare_bc2xci.R             # Compare two bc2xci outputs (kappa, Bowker test)
├── Comparison between two outputs.R  # Ad hoc pairwise comparison of two bc2xci files
├── annotation.R                 # Annotate XCISE VCF output with GTF gene models
└── annotate_vcf_with_gtf.R      # Auto-download GTF (GENCODE v46 / Ensembl 113)
```

---

## Requirements

| Tool / Package | Version |
|---|---|
| Perl | ≥ 5.26 |
| `Parallel::ForkManager` | any |
| `File::Temp` | any |
| `List::Util` | any |
| `samtools` | ≥ 1.15 |
| R | ≥ 4.1 |
| R `data.table` | ≥ 1.14 |
| R `rtracklayer` *(optional)* | ≥ 1.56 |

Install Perl dependencies:

```bash
cpan Parallel::ForkManager File::Temp List::Util
```

Install R dependencies:

```r
install.packages("data.table")
BiocManager::install("rtracklayer")   # optional but recommended
```

---

## Input Requirements

| Input | Description |
|---|---|
| BAM file(s) | WASP-corrected, indexed; must carry `CB` / `UB` tags and `vG` / `vA` / `vW` custom tags |
| VCF file | Heterozygous SNVs on the X chromosome (e.g., from GATK or BCFtools) |
| GTF file | Gene annotation for VCF annotation step (GENCODE v46 auto-downloaded if omitted) |

> **Note:** WASP re-mapping is required upstream. Reads without the `vW:i:1` tag are silently skipped.
>
> EXCISE does not re-implement the upstream preprocessing. Follow the step-by-step guide in the [original XCISE README](https://github.com/Vityay/XCISE#readme) (STAR/STARsolo alignment → BCFtools variant calling → common-SNV filtering → `vcf4wasp.pl` → STAR in WASP mode with `--waspOutputMode SAMtag --outSAMattributes vA vG ...`) to produce compatible BAM and VCF inputs.

---

## Usage

### Main pipeline (`excise.pl`)

```bash
perl excise.pl \
  -o  SAMPLE_NAME         \  # output prefix
  -s  hets_chrX.vcf       \  # heterozygous SNV VCF
  -b  sample_WASP.bam     \  # WASP-corrected BAM (repeatable)
  -r  X                   \  # chromosome (default: X)
  -j  8                   \  # parallel tries
  -t  100                 \  # total tries
  -cs 10000               \  # MCMC steps per try
  -T0 1.0 -Tmin 0.01      \  # simulated annealing schedule
  -blocksize 3            \  # SNVs flipped per proposal
  -seed 42                \  # reproducible RNG
  -pretrain                  # observation-based initialization
```

### Greedy multi-try variant (`xcise_mt.pl`)

```bash
perl xcise_mt.pl \
  -o SAMPLE -s hets_chrX.vcf -b sample_WASP.bam \
  -j 8 -t 100 -samthreads 4 -seed 42 -pretrain
```

### MCMC + co-segregation graph (`xcise_mcmc.pl`)

```bash
perl xcise_mcmc.pl \
  -o SAMPLE -s hets_chrX.vcf -b sample_WASP.bam \
  -j 8 -cs 4000 -blocksize 3 -init observed \
  -writegraph 1 -co_min_pairs 2 -samthreads 8
```

---

## Full Parameter Reference

### Core parameters

| Flag | Default | Description |
|---|---|---|
| `-o <prefix>` | — | Output file prefix |
| `-s <vcf>` | — | Input VCF with heterozygous SNVs |
| `-b <bam>` | — | Input BAM (can be repeated for multiple files) |
| `-r <chrom>` | `X` | Chromosome to analyse |
| `-t <int>` | `100` | Number of random restarts (tries) |
| `-u <int>` | `10` | Min UMIs per SNV to retain |
| `-m <float>` | `0` | Min minor allele frequency per SNV |
| `-p <int>` | `5` | Discordant UMI penalty |

### Parallelism & I/O

| Flag | Default | Description |
|---|---|---|
| `-j <int>` | `1` | Parallel processes (one per try) |
| `-samthreads <int>` | `1` | Threads passed to `samtools view -@` |
| `-seed <int>` | — | Base RNG seed for reproducibility |

### MCMC / simulated annealing

| Flag | Default | Description |
|---|---|---|
| `-cs <int>` | `10000` | MCMC chain steps per try |
| `-T0 <float>` | `1.0` | Initial temperature |
| `-Tmin <float>` | `0.01` | Final temperature |
| `-blocksize <int>` | `1` | SNVs flipped per proposal |
| `-blocktype <str>` | `random` | Proposal type: `random` or `window` |

### Pretraining (unsupervised seeding)

| Flag | Default | Description |
|---|---|---|
| `-pretrain` | off | Enable graph-based initialization |
| `-pre_cap <int>` | `5000` | Max SNV-pair edges sampled per cell |
| `-pre_wmin <int>` | `3` | Min edge weight for BFS propagation |
| `-pre_tieskip` | on | Skip allele calls with 1:1 tie |

### `xcise_mcmc.pl`-specific

| Flag | Default | Description |
|---|---|---|
| `-init <str>` | `observed` | Init mode: `observed` or `random` |
| `-co_min_pairs <int>` | `2` | Min co-observed molecules to trust a graph edge |
| `-init_jitter <float>` | `0.0` | Fraction of seeded SNVs to randomly perturb |
| `-writegraph <0\|1>` | `1` | Write co-segregation graph TSV + DOT files |
| `-anchor_sign <+1\|-1>` | `+1` | Anchor sign for connected components |

---

## Outputs

| File | Description |
|---|---|
| `*_XCISE_bc2xci.txt` | Per-cell XCI assignment (barcode, hap1 UMIs, hap2 UMIs, label) |
| `*_XCISE.vcf` | Phased SNV VCF with X1-allele and allele-depth annotations |
| `*_XCISE_summary.txt` | Run statistics: best score, concordancy rate, cell counts |
| `*_coGraph_nodes.tsv` | Co-segregation graph node table *(xcise_mcmc.pl only)* |
| `*_coGraph_edges.tsv` | Co-segregation graph edge table *(xcise_mcmc.pl only)* |
| `*_coGraph.dot` | Graphviz DOT file for graph visualization *(xcise_mcmc.pl only)* |
| `*.annot.tsv` | VCF with gene name, biotype, and nearest-gene annotations |

### Cell classification thresholds

| Label | Criterion |
|---|---|
| `X1` | hap1 ≥ 2 UMIs **and** hap1 / (hap1 + hap2) ≥ 0.9 |
| `X2` | hap2 ≥ 2 UMIs **and** hap2 / (hap1 + hap2) ≥ 0.9 |
| `Both` | both hap1 > 0 and hap2 > 0 (not meeting X1/X2 threshold) |
| `Low_coverage` | exactly one haplotype has evidence, but count < 2 |
| `Unknown` | no allele-specific UMIs observed |

---

## Downstream Utilities

### Compare two runs

```bash
Rscript compare_bc2xci.R  run1_chrX_XCISE_bc2xci.txt  run2_chrX_XCISE_bc2xci.txt
```

Outputs: overlap statistics, contingency table (LC/Unknown dropped), **Cohen's κ**, and a **Bowker symmetry test** (or McNemar for 2×2 tables).

### Annotate the output VCF

```bash
# Auto-download GENCODE v46 and annotate
Rscript annotate_vcf_with_gtf.R \
  --vcf=SAMPLE_chrX_XCISE.vcf \
  --gtf=gencode_v46_GRCh38 \
  --out=SAMPLE_chrX_XCISE.annot.tsv

# Use a local GTF
Rscript annotation.R \
  --vcf=SAMPLE_chrX_XCISE.vcf \
  --fasta=/path/to/genome.fa \
  --gtf=/path/to/annotation.gtf \
  --out=SAMPLE_chrX_XCISE.annot.tsv
```

Outputs: per-SNV gene name, biotype, and nearest upstream/downstream gene with distance in bp.

---

## Recommended Run Settings

| Dataset size | Recommended flags |
|---|---|
| Small (< 1,000 cells, < 50 SNVs) | `-t 50 -cs 5000 -j 4` |
| Medium (1,000–5,000 cells) | `-t 100 -cs 10000 -j 8 -pretrain` |
| Large (> 5,000 cells) | `-t 100 -cs 20000 -j 16 -pretrain -blocksize 3` |

---

## Citation

EXCISE is a derivative of XCISE. If you use EXCISE in your work, **please cite the original XCISE paper**:

> Henning RH, Rust TM, Dijksterhuis K, Eggen BJL, Guryev V. **Single-cell X-chromosome inactivation analysis links biased chimerism to differential gene expression and epigenetic erosion.** *bioRxiv* (2024). doi: [10.1101/2024.08.29.610317](https://doi.org/10.1101/2024.08.29.610317)
>
> Original code: <https://github.com/Vityay/XCISE>

and acknowledge this repository:

> **EXCISE:** Zhuang Z. *EXCISE — An Enhanced XCISE Pipeline.* GitHub: <https://github.com/Solidarity31/XCISE_New>

<details>
<summary>BibTeX</summary>

```bibtex
@article{henning2024xcise,
  title   = {Single-cell X-chromosome inactivation analysis links biased chimerism
             to differential gene expression and epigenetic erosion},
  author  = {Henning, Robert H. and Rust, Thomas M. and Dijksterhuis, Kasper and
             Eggen, Bart J. L. and Guryev, Victor},
  journal = {bioRxiv},
  year    = {2024},
  doi     = {10.1101/2024.08.29.610317}
}

@article{wang2025femxpress,
  title   = {FemXpress: Systematic Analysis of X Chromosome Inactivation
             Heterogeneity in Female Single-Cell RNA-Seq Samples},
  author  = {Wang, Xin and Ma, Yingke and Li, Fan and Cui, Wentao and Pan, Tianshi and
             Wang, Siqi and Ma, Sinan and Shan, Qingtong and Liu, Chao and Wang, Yukai and
             Zhang, Ying and Zhou, Yuanchun and Li, Wei and Wang, Pengfei and
             Zhou, Qi and Feng, Guihai},
  journal = {Advanced Science},
  volume  = {12},
  number  = {35},
  year    = {2025},
  doi     = {10.1002/advs.202504754}
}
```

</details>

---

## Related Tools

- **[XCISE](https://github.com/Vityay/XCISE)** — the original pipeline that EXCISE extends (Henning *et al.*, bioRxiv 2024).
- **FemXpress** — an independent tool that uses X-linked SNPs to group cells in female scRNA-seq data by the parental origin of the inactivated X chromosome, without requiring parental genotypes, and additionally identifies genes that escape X inactivation. Wang X, *et al.* **FemXpress: Systematic Analysis of X Chromosome Inactivation Heterogeneity in Female Single-Cell RNA-Seq Samples.** *Advanced Science* 12(35) (2025). doi: [10.1002/advs.202504754](https://doi.org/10.1002/advs.202504754)

EXCISE and FemXpress address the same biological question (per-cell XCI state from allele-specific scRNA-seq) with different algorithms; running both on the same sample and comparing per-cell assignments (e.g. with `compare_bc2xci.R` after reformatting) is a useful orthogonal check. If you use FemXpress, please cite its paper.

---

## Acknowledgements

EXCISE is built directly on the conceptual and algorithmic foundation of the original **[XCISE](https://github.com/Vityay/XCISE)** framework by Henning, Rust, Dijksterhuis, Eggen & Guryev. All core biological logic — allele-specific UMI phasing, the concordance-minus-discordance scoring function, and the X1/X2 classification scheme — originates with XCISE. This repository extends that foundation with improved optimization, initialization, and downstream analysis capabilities.

---

<p align="center">
  <sub>For questions or issues, please open a <a href="../../issues">GitHub Issue</a>.</sub>
</p>
