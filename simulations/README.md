# EXCISE 1.5 simulations

Simulation and benchmarking for `excise15.py`, the Python/numba port of `excise1.5.pl`.
Everything in this directory comes from **simulated data only**. No real
biological datasets are read or written here.

```
simulations/
├── run_xci_skew_sweep.py     # XCI-skew experiment: simulate -> excise15.py -> score
├── run_speed_benchmark.py    # Perl vs Python runtime on increasing cell numbers
├── make_figures.R            # all figures (ggplot2); reads results/ only
├── results/                  # numerical results (TSV / CSV) behind every figure
├── figures/                  # PDF + 300 dpi PNG
├── sweep.log, speed.log      # run logs
└── data/                     # simulated BAMs + per-run outputs (git-ignored, regenerable)
```

## Reproduce

```bash
pip install numpy numba pysam          # plus samtools and perl Parallel::ForkManager
python3 simulations/run_xci_skew_sweep.py --reps 10 --jobs 4
python3 simulations/run_speed_benchmark.py --jobs 4
Rscript simulations/make_figures.R     # ggplot2, dplyr, tidyr, readr, patchwork, scales
```

Both runners reuse existing datasets and runs, and never overwrite them. Run the two
experiments one after the other, not at the same time, so that the runtimes are not
measured under CPU contention.

## Designs

**XCI-skew sweep** (`tests/simulate_xci_data.py` with `--assign exact`). Only the
X1:X2 ratio changes between conditions:

| Parameter | Value |
|---|---|
| X1:X2 among **monoallelic** cells | 50:50, 55:45, 70:30, 90:10, 95:5, 99:1 |
| Both cells, as a share of **all** cells | 5% (100 of 2,000) |
| Cells / SNVs | 2,000 / 1,000 on a simulated 5 Mb chrX |
| Coverage | 20–60 reads per cell (mean 40), each cell expresses 3–25 SNVs |
| Allelic noise | 2% of reads show the inactive allele |
| Escape SNVs | 5% (biallelic in every cell) |
| Replicates | 10 per ratio; replicate *k* uses seed 100+*k* for every ratio |
| EXCISE | `excise15.py -t 100 -j 4 -seed <seed> -pretrain` |

At 99:1, the 1,900 monoallelic cells are 1,881 X1 and **19 X2**, so the minority
class is 0.95% of all cells.

**Speed benchmark.** 2,000 SNVs, ~60 reads per cell, and 1,000–16,000 cells × 3
seeds. Both tools run with `-t 4 -j 4 -pretrain`. The 8,000-cell, seed 11 dataset is the
dataset of the originally reported benchmark, which is reused.

## Metrics (`tests/evaluate_xci_sim.py`)

EXCISE's X1/X2 names are arbitrary, so each run gets **one global orientation**,
chosen from the SNV haplotypes. Every SNV and cell metric uses that same mapping.

- **Haplotype accuracy**: correctly oriented / phased non-escape SNVs.
- **Unassigned SNVs**: SNVs that pass the filters but are left `Unk`, reported separately for
  non-escape SNVs and for escape SNVs (where `Unk` is the desired outcome).
- **Cell accuracy**: correct calls / all cells. `Low_coverage` and `Unknown` count as
  errors. Mean per-class recall over X1, X2 and Both is also given.
- **Minority-class recall and precision**: for the true X2 cells.
- **Informative cell recovery**: cells called X1, X2 or Both / all cells, and true
  X1/X2 cells called X1 or X2.
- **Runtime**: wall-clock time of the whole `excise15.py` run.
