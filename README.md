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

XCISE (*XCI calling from Single cell Expression data*) was developed by Henning, Rust, Dijksterhuis, Eggen & Guryev ([bioRxiv 2024](https://doi.org/10.1101/2024.08.29.610317); code: [Vityay/XCISE](https://github.com/Vityay/XCISE)). `xcise_mt.pl` and `xcise_mcmc.pl` are direct modifications of the original `xcise.pl`; `excise.pl` and `excise1.5.pl` are derived from them.

Key advances over the original XCISE include **MCMC-based haplotype optimization** with simulated annealing, **observation-based co-segregation seeding**, **parallel multi-try execution**, and a suite of downstream utilities for annotation and inter-run comparison.

### What's new in EXCISE v1.5 (`excise1.5.pl`)

`excise1.5.pl` is the latest script. It builds on the greedy multi-try optimizer of `xcise_mt.pl` (it does **not** use MCMC) and adds explicit handling of SNVs that cannot be phased:

- **Escape-candidate SNVs.** A biallelic SNV can stay unassigned (`dir = 0`, reported as `X1A=Unk`) when neither orientation improves the score, so putative escape genes do not distort cell calls. Unassigned SNVs get a symmetric test of both orientations at every pass.
- **Monoallelic SNVs.** SNVs where only one allele is observed are flagged (`MONO`), are never moved to 0 by the optimizer, and are re-oriented after optimization by a majority vote of cells classified from biallelic SNVs (Phase 2).
- **Optional Phase 3 (`-ph`).** Remaining unassigned biallelic SNVs are rescued only if one orientation is supported with concordance ≥ `-ph3_ratio` (default 0.90).
- **Uncertainty-scaled jitter.** With `-pretrain`, each try perturbs the seeded orientation of every SNV with probability `jitter × 1/(1 + w)`, where `w` is the weakest co-segregation edge on its seeding path, so weakly supported SNVs are explored more.
- **Cell-barcode whitelist** (`-wl`, plain or gzipped, optional `-wl_strip_gem`).
- **Richer outputs.** A valid VCF with `##INFO` headers and per-haplotype UMI counts, a haplotype table, and a per-cell × per-SNV UMI profile.
- **Defaults.** `-samthreads` defaults to `-j`, and the minimum UMIs per SNV (`-u`) defaults to **4** (10 in the other scripts).

---

## Features

| Feature | XCISE (original) | EXCISE (`excise.pl`, `xcise_mt.pl`, `xcise_mcmc.pl`) | EXCISE v1.5 (`excise1.5.pl`) |
|---|:---:|:---:|:---:|
| WASP mapping-bias correction | ✅ | ✅ | ✅ |
| Greedy haplotype optimization | ✅ | ✅ | ✅ |
| MCMC + simulated annealing | ❌ | ✅ | ❌ |
| Multi-SNV block proposals | ❌ | ✅ | ❌ |
| Observation-based co-segregation seed (`-pretrain`) | ❌ | ✅ | ✅ |
| Uncertainty-scaled per-try jitter | ❌ | ❌ | ✅ |
| Parallel tries (`-j`) | ❌ | ✅ | ✅ |
| Reproducible RNG (`-seed`) | ❌ | ✅ | ✅ |
| Co-segregation graph export (TSV + DOT) | ❌ | ✅ *(xcise_mcmc.pl)* | ❌ |
| Symmetric ±1 test for unassigned (escape-candidate) SNVs | ❌ | ❌ | ✅ |
| Monoallelic SNV flagging + cell-vote correction | ❌ | ❌ | ✅ |
| Optional Phase-3 rescue of unassigned SNVs (`-ph`) | ❌ | ❌ | ✅ |
| Cell-barcode whitelist (`-wl`) | ❌ | ❌ | ✅ |
| Haplotype table + per-cell × per-SNV profile | ❌ | ❌ | ✅ |
| VCF annotation with GTF *(R utilities)* | ❌ | ✅ | ✅ |
| Inter-run comparison + Cohen's κ *(R utilities)* | ❌ | ✅ | ✅ |

---

## Repository Structure

```
EXCISE/
├── excise1.5.pl                 # EXCISE v1.5: greedy multi-try + escape/monoallelic handling (latest)
├── excise15.py                  # Fast Python/numba port of excise1.5.pl (same flags and outputs)
├── cpp/excise15.cpp, Makefile   # Fastest: C++/htslib port of excise1.5.pl (same flags and outputs)
├── simulations/                 # XCI-skew simulations, Perl/Python/C++ benchmark, figures
├── extract_bc2xci_overlap.R     # Cells shared by two bc2xci files with each run's label
├── tests/simulate_xci_data.py   # Simulated WASP BAM + VCF with known truth (testing/benchmarking)
├── excise.pl                    # MCMC + simulated annealing + pretraining
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
| `gzip` | any *(needed only for gzipped VCF / whitelist input)* |
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
| VCF file | Heterozygous SNVs on the X chromosome (e.g., from GATK or BCFtools); may be gzipped |
| Barcode whitelist *(optional, v1.5)* | One barcode per line (first whitespace-delimited column is used, so a 10x `barcodes.tsv[.gz]` works) |
| GTF file | Gene annotation for VCF annotation step (GENCODE v46 auto-downloaded if omitted) |

> **Note:** WASP re-mapping is required upstream. Reads without the `vW:i:1` tag are silently skipped. Reads without a `CB` tag fall back to the `RG` tag as the cell ID (plate-based / Smart-seq data); reads without `UB` get a synthetic UMI built from barcode, chromosome, position, flag, CIGAR and template length.
>
> EXCISE does not re-implement the upstream preprocessing. Follow the step-by-step guide in the [original XCISE README](https://github.com/Vityay/XCISE#readme) (STAR/STARsolo alignment → BCFtools variant calling → common-SNV filtering → `vcf4wasp.pl` → STAR in WASP mode with `--waspOutputMode SAMtag --outSAMattributes vA vG ...`) to produce compatible BAM and VCF inputs.

---

## Usage

### EXCISE v1.5 (`excise1.5.pl`)

```bash
perl excise1.5.pl \
  -o SAMPLE_NAME \
  -s hets_chrX.vcf.gz \
  -b sample_WASP.bam \
  -r X \
  -t 100 -j 8 \
  -seed 42 \
  -pretrain -jitter 0.25 \
  -wl barcodes.tsv.gz -wl_strip_gem \
  -ph -ph3_ratio 0.9
```

`-o` output prefix · `-s` heterozygous SNV VCF (plain or gzipped) · `-b` WASP-corrected BAM(s), several may follow one `-b` · `-r` chromosome · `-t`/`-j` tries, run 8 at a time (samtools also gets 8 threads) · `-seed` reproducible RNG · `-pretrain -jitter` co-segregation seeding with uncertainty-scaled jitter · `-wl -wl_strip_gem` *optional*, keep only whitelisted cells · `-ph -ph3_ratio` *optional*, Phase-3 rescue of unassigned SNVs.

The pipeline runs these stages in order:

1. **Read evidence.** It collects allele-specific UMIs per SNV and per cell from `samtools view` (skipping secondary alignments, reads without `vW:i:1`, reads with homopolymer runs ≥ `-ms`, and soft-clipped reads if `-sm` is set). It then drops SNVs with fewer than `-u` UMIs or with allele frequency outside [`-m`, 1 − `-m`].
2. **Pretraining** *(`-pretrain`)*. Each cell makes a majority allele call per SNV. Every pair of SNVs called in the same cell then votes ±1 on a co-segregation graph (capped at `-pre_cap` pairs per cell). BFS from the highest-degree SNV propagates orientations along edges with |weight| ≥ `-pre_wmin`. SNVs that are not connected stay unassigned.
3. **Parallel greedy tries.** Each try starts from the (jittered) seed, or from a random orientation without `-pretrain`, and hill-climbs the XCISE score *concordant − p × discordant* over orientations {−1, 0, +1}. The optimizer never moves a monoallelic SNV to 0. Ties prefer 0, and when a pass stalls the search makes one random probe before stopping. The highest-scoring try is kept.
4. **Phase 2.** Monoallelic SNVs are re-oriented by the majority of covering cells already called X1 or X2 from biallelic SNVs.
5. **Phase 3** *(`-ph`)*. Two passes over the remaining unassigned biallelic SNVs. The first picks the orientation that yields more cells with ≥ 90% haplotype purity. The second votes using cells classified from phased SNVs. An orientation is assigned only if its support (the dominant-haplotype UMI fraction in the first pass, the concordance ratio in the second) is ≥ `-ph3_ratio`.
6. **Output.** Haplotype labels are swapped if needed so that `#X2 cells ≥ #X1 cells`, and all output files are written.

> `excise1.5.pl` does not accept the MCMC flags (`-cs`, `-T0`, `-Tmin`, `-blocksize`, `-blocktype`) or the `xcise_mcmc.pl` graph flags. Any unknown flag aborts with `Unexpected/incomplete parameter`.

### Fast Python port (`excise15.py`)

`excise15.py` runs the same algorithm as `excise1.5.pl`, with the same flags and the same output files, but much faster:

```bash
pip install numpy numba pysam
python3 excise15.py -o SAMPLE_NAME -s hets_chrX.vcf.gz -b sample_WASP.bam \
  -r X -t 100 -j 8 -seed 42 -pretrain -ph
```

The main speedup is algorithmic. After each single-SNV move, `excise1.5.pl` rescores the run by looping over **every cell**. The score is a sum of per-cell terms, so `excise15.py` re-evaluates only the cells that cover the moved SNV. The resulting score is identical, but each move is far cheaper. The hill-climber, pretraining graph and BFS are compiled with numba. BAM evidence is packed once into flat integer arrays, so the inner loops no longer split `"umi\tcb"` strings. BAM input is read with `pysam`, so no `samtools` binary is needed, but the BAM must be indexed.

Validation on simulated data (`tests/simulate_xci_data.py`):

- When both scripts reach the same optimum, all five output files are byte-identical to `excise1.5.pl`. This holds with and without `-pretrain`, and with `-ph`.
- Given the same optimizer result, Phase 2, both Phase-3 passes and all output writers match exactly, including on sparse data with many monoallelic SNVs.
- On 1,000–16,000 cells (2,000 SNVs, `-t 4 -j 4 -pretrain`, 3 datasets per size), it was 2.5× to 13.7× faster than Perl. The optimizer took 0.04–0.06 s, so most of the remaining time is BAM reading and numba compilation. See `simulations/`.
- The compiled numba kernels are cached on disk (`__pycache__/`). The first run in a new location compiles them, and later runs skip that step, saving about 1.3 s per run.

Results are not bit-identical to the Perl script for the same `-seed`, because the RNGs differ and Perl iterates hashes in random order. On hard, sparse data, individual tries can therefore reach different local optima, as they do between two Perl runs with different seeds.

### Fast C++ port (`cpp/excise15`)

`cpp/excise15` is the same program again in C++17 with htslib, with the same flags and the same five output files. It is the fastest option:

```bash
make -C cpp                    # needs htslib + zlib (apt install libhts-dev, or conda install -c bioconda htslib)
cpp/excise15 -o SAMPLE_NAME -s hets_chrX.vcf.gz -b sample_WASP.bam \
  -r X -t 100 -j 8 -seed 42 -pretrain -ph
```

The Python port spends almost all of its time outside the optimizer: decoding BAM records one at a time through pysam, plus numba compilation at start-up. The C++ port removes both. htslib decodes the BAM natively, and the binary is compiled ahead of time. The chromosome is also split into chunks that are decompressed and parsed in parallel on `-samthreads` threads (default `-j`). Each read belongs to the chunk holding its start, and chunks are merged in genomic order, so the result is the same as a single sequential pass. The UMI collapsing (as in `%umis` / `%cb2allele`) uses one sort over a flat array, and the co-segregation votes are summed by sorting rather than through hash tables. The optimizer uses the same delta scoring as the Python port, and tries run on `-j` threads. For a non-system htslib, use `make -C cpp HTSLIB_PREFIX=$CONDA_PREFIX`.

Validation against `excise1.5.pl`:

- On the 15 benchmark datasets, all five output files are byte-identical in 12 runs. In the remaining 3 runs (1,000 cells, where the search is hardest), every cell gets the same `X1`/`X2`/`Both` label as in Perl. The haplotype counts differ slightly because the search reached a different local optimum.
- Given Perl's optimizer result (`make -C cpp excise15_test`, `INJECT=<file>`), Phase 2, both Phase-3 passes and all writers are byte-identical on sparse data. This holds with and without `-ph`, at several `-u` and `-p` settings.
- The input paths also match Perl byte for byte: a gzipped whitelist with `-wl_strip_gem`, a gzipped VCF, reads without `UB` (soft UMIs), `-sm`, `-ms`, and several BAMs.
- Speed: 174–199× faster than Perl at every size from 1,000 to 16,000 cells, which is 14–75× faster than the Python port (geometric means per size). At 16,000 cells (0.96 M reads), it took 0.95–0.99 s, compared with 181–185 s for Perl and 13–14 s for Python. On a BAM shaped like real 10x data (10 M chrX reads, 90% without WASP tags), BAM reading took 1.5 s, against 5.0 s single-threaded and 2.0 s for `samtools view -@4 -c`. The output was byte-identical.

### MCMC pipeline (`excise.pl`)

```bash
perl excise.pl \
  -o SAMPLE_NAME \
  -s hets_chrX.vcf \
  -b sample_WASP.bam \
  -r X \
  -j 8 -t 100 \
  -cs 10000 \
  -T0 1.0 -Tmin 0.01 \
  -blocksize 3 \
  -seed 42 \
  -pretrain
```

`-cs` MCMC steps per try · `-T0`/`-Tmin` simulated-annealing schedule · `-blocksize` SNVs flipped per proposal · `-pretrain` observation-based initialization; other flags as above.

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

### `excise1.5.pl`-specific

`excise1.5.pl` accepts the core, parallelism and pretraining flags above, with these differences and additions:

| Flag | Default | Description |
|---|---|---|
| `-u <int>` | **`4`** | Min UMIs per SNV (10 in the other scripts) |
| `-samthreads <int>` | **= `-j`** | Threads for `samtools view -@` |
| `-x` | off | Randomly scramble alleles (negative control) |
| `-i` | off | Adopt the best score from an existing `*_XCISE_summary.txt` if it is higher |
| `-sm` | off | Skip soft-clipped alignments |
| `-ms <int>` | `15` | Skip reads containing a mononucleotide run of this length |
| `-wl <file>` | — | Barcode whitelist (plain text/TSV, or gzipped; gzip is detected by extension or magic bytes) |
| `-wl_strip_gem` | off | Strip a trailing `-<digits>` (e.g. `-1`) from both BAM and whitelist barcodes before matching; requires `-wl` |
| `-jitter <float>` | `0.25` | Scales per-try flip probability `jitter × 1/(1 + bottleneck edge weight)`; clamped to [0, 1]. Without `-pretrain`, the starting orientation is random instead |
| `-ph` | off | Enable Phase-3 assignment of unassigned biallelic SNVs |
| `-ph3_ratio <float>` | `0.90` | Min support ratio for Phase 3 to assign an orientation |

> Tip: if your BAM barcodes carry the 10x GEM suffix (`AAAC…-1`) but the whitelist does not (or vice versa), use `-wl_strip_gem`. Otherwise no cells match and the run aborts with `Zero SNVs after BAM reading`.

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
| `*_XCISE_summary.txt` | Run statistics: best score, concordancy rate, cell counts *(v1.5 adds SNV counts — informative / non-informative / monoallelic — and whitelist statistics)* |
| `*_XCISE_haplotypes.tsv` | `CHROM POS ID REF ALT X1_allele X2_allele mono`, one row per phased SNV; unassigned SNVs are omitted *(excise1.5.pl only)* |
| `*_XCISE_cell_snv_profile.tsv` | One row per (cell, SNV) with ≥ 1 UMI: `cell_barcode pos rs ref alt snv_dir x1a ref_umis alt_umis total_umis hap_x1_umis hap_x2_umis cell_xci mono` *(excise1.5.pl only)* |
| `*_coGraph_nodes.tsv` | Co-segregation graph node table *(xcise_mcmc.pl only)* |
| `*_coGraph_edges.tsv` | Co-segregation graph edge table *(xcise_mcmc.pl only)* |
| `*_coGraph.dot` | Graphviz DOT file for graph visualization *(xcise_mcmc.pl only)* |
| `*.annot.tsv` | VCF with gene name, biotype, and nearest-gene annotations |

All file names carry the prefix `<prefix>_chr<chrom>_`, for example `SAMPLE_chrX_XCISE_bc2xci.txt`. `excise1.5.pl` keeps the `_XCISE_` infix so its outputs drop into the same downstream tools.

#### VCF INFO fields written by `excise1.5.pl`

| Tag | Meaning |
|---|---|
| `X1A` | Allele on the X1 haplotype: `Ref`, `Alt`, or `Unk` (unassigned / escape candidate) |
| `AD` | UMI depth as `ref,alt` |
| `X1`, `X2` | UMIs supporting the X1 / X2 haplotype (0 for `Unk` SNVs) |
| `DP` | `X1 + X2`; also written in the `QUAL` column |
| `MONO` | Flag: monoallelic SNV, oriented by cell vote |

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

# Use a local GTF; merge the excise1.5 haplotype table, add XCI status, summary and plots
Rscript annotation.R \
  --vcf=SAMPLE_chrX_XCISE.vcf \
  --hap=SAMPLE_chrX_XCISE_haplotypes.tsv \
  --fasta=/path/to/genome.fa \
  --gtf=/path/to/annotation.gtf \
  --xci_ref=xci_reference.tsv \
  --summary=TRUE --plot=TRUE \
  --out=SAMPLE_chrX_XCISE.annot.tsv
```

`annotation.R` outputs per-SNV gene name, biotype, XCI status (with `--xci_ref`), and the nearest upstream/downstream gene with distance in bp. `--hap` merges the haplotype TSV (X1/X2 alleles) with the VCF (`X1A` tag and `AD` UMI counts), cross-checks the two phasings and writes `*.haplotype_merged.tsv` with per-SNV X1/X2 UMIs. `--summary` writes a readable `*.summary.txt` and a `*.gene_summary.tsv`. `--plot` draws X1/X2 UMIs along chrX, the per-SNV UMI distribution on each haplotype, and the per-SNV X1 fraction. Without `--hap`, the phase comes from the VCF `X1A` tag, so the original XCISE VCF also works. The GTF is read through `grep`, which keeps only the `gene` rows before parsing (about 2% of a GENCODE GTF). This makes the import roughly 9–13× faster than parsing the whole file, and gives the same gene table. On systems without `grep` (or `gzip`, for `.gtf.gz`), it falls back to rtracklayer or a full `fread`.

---

## Recommended Run Settings

| Dataset size | Recommended flags |
|---|---|
| Small (< 1,000 cells, < 50 SNVs) | `-t 50 -cs 5000 -j 4` |
| Medium (1,000–5,000 cells) | `-t 100 -cs 10000 -j 8 -pretrain` |
| Large (> 5,000 cells) | `-t 100 -cs 20000 -j 16 -pretrain -blocksize 3` |

These settings apply to the MCMC scripts. For `excise1.5.pl`, a reasonable starting point is `-t 100 -j <cores> -seed 42 -pretrain`. Add `-ph` only if you want unassigned SNVs rescued, and add `-u 10` to match the stricter SNV filter of the other scripts.

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
