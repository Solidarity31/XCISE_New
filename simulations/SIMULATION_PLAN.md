# EXCISE simulation plan

This plan covers what has been simulated so far, what is still needed for a preprint, and the design of each test.
All tests use `tests/simulate_xci_data.py` (extended where noted). They are scored with `tests/evaluate_xci_sim.py`,
which allows one global X1/X2 label switch per run.

## Shared design

| Item | Value |
|---|---|
| Truth per cell | `X1` (haplotype A active), `X2` (haplotype B active), `Both` (biallelic) |
| Truth per SNV | which allele is on haplotype A; escape flag |
| EXCISE settings | defaults unless the test varies them: `-t 100 -pretrain -u 4 -p 5` |
| Replicates | ≥ 5 independent seeds per condition (10 when the minority class is small) |
| Implementation | `cpp/excise15`, which matches `excise1.5.pl`; spot-check each new test with the Perl script |

**Metrics, reported for every test:**

- **Phasing:** haplotype accuracy, and the fraction of SNVs left unassigned (non-escape and escape separately).
- **Calls overall:** cell accuracy, with `Low_coverage`/`Unknown` counted as errors.
- **Per class:** recall and precision for each class, and wrong-X calls (the most serious error).
- **Recovery:** the fraction of informative cells recovered.
- **Cost:** runtime and peak memory.

Never report overall accuracy alone. When one class is rare, always report the actual number of cells in it.

## Completed

| # | Test | Design | Main result | Where |
|---|---|---|---|---|
| C1 | Speed: original XCISE vs EXCISE 1.5 (Perl, Python, C++) | 1,000–16,000 cells × 3 seeds; 2,000 SNVs; `-t 4` | C++ is 1,849–5,812× faster than the original, with the same calls | Fig. 1, `results/speed_benchmark.tsv` |
| C2 | Default tries | `-t 100`, 1,000 cells | Original 56 min, C++ 0.15 s, same accuracy (0.996) | `results/default_tries_benchmark.tsv` |
| C3 | Equivalence of implementations | Byte-identical outputs; injected optimizer result; whitelist, gz, soft UMI, `-sm`, `-ms`, multi-BAM | Identical whenever the same optimum is reached | `README.md` |
| C4 | XCI skew | X1:X2 = 50:50 … 99:1 among monoallelic cells; 5% Both; 10 independent replicates (+ paired control) | No effect at standard coverage; no wrong-X calls; minority recall ≥ 0.968 | Figs. 2–5 |
| C5 | `-u` × `-p` sweep | 7 × 7 grid × 4 regimes (standard, low coverage, high noise, sparse SNVs) × 5 replicates | Defaults are near-optimal except under high noise; `-p` extremes harm results | Fig. 6 |

## Still needed

Priority 1 is needed for a methods paper, priority 2 strengthens it, and priority 3 is optional.

| # | Priority | Question | Factors and levels (others at standard) | Reps | New simulator features |
|---|---|---|---|---|---|
| N1 | 1 | When does EXCISE break? Coverage × noise | reads/cell {5, 10, 20, 40, 80} × noise {0, 2, 5, 10, 15%} | 5 | none |
| N2 | 1 | How many cells and SNVs are needed to phase? | cells {100, 250, 500, 1,000, 5,000, 20,000} × SNVs {50, 100, 300, 1,000, 3,000} | 5 | none |
| N3 | 1 | Does skew matter when data are hard? | X1:X2 {50:50, 90:10, 99:1} × regime {standard, low coverage, high noise, sparse} | 10 | none |
| N4 | 1 | Realistic expression | Gene-clustered SNVs (1–10 per gene); negative-binomial gene expression; reads spanning 2+ SNVs; 1–5 reads per UMI | 5 | gene model, multi-SNV reads, PCR duplicates |
| N5 | 1 | Escape | escape fraction {0, 5, 15, 30%} × escape level (inactive-X expression {10, 30, 50%} of active) | 5 | partial escape |
| N6 | 1 | Biallelic cells and doublets | Both fraction {0, 5, 20, 40%}; doublets formed by pooling the reads of one X1 and one X2 cell | 5 | doublet model |
| N7 | 2 | VCF errors | false heterozygous SNVs (truly homozygous) {0, 5, 20%}; missing true hets {0, 20%} | 5 | genotype errors |
| N8 | 2 | Contamination | ambient RNA / barcode swapping: {0, 2, 5, 10%} of reads from other cells | 5 | read mixing |
| N9 | 2 | Residual reference bias after WASP | REF over-representation {1.0, 1.1, 1.3} | 5 | allele bias |
| N10 | 2 | Optimizer settings | `-t` {1, 4, 10, 25, 100} × `-pretrain` {on, off}; `-jitter` {0, 0.25, 0.5}; `-ph` {off, on}; `-m` {0, 0.05, 0.1}; seed stability | 5 | none |
| N11 | 2 | The 90% purity threshold | Post-hoc relabelling at {0.8, 0.85, 0.9, 0.95} on N1 runs (no re-run needed; uses hap1/hap2 counts) | — | none |
| N12 | 2 | Comparison with other tools | FemXpress, XC-ID (and XIPHER, if available) on N1/N3/N4 datasets | 5 | input converters |
| N13 | 3 | Semi-real ground truth: male mixture | Pool cells from two male donors (one X each, different haplotypes) into a pseudo-female dataset. Donor identity = true X1/X2. Vary mixing ratio (50:50 … 95:5) | per available donors | real BAMs + donor genotypes |
| N14 | 3 | Memory and runtime at real scale | 10,000–100,000 cells; 10–50 M chrX reads; peak RSS + wall time | 3 | large-BAM generator (as in the realistic-BAM test) |

### Notes on the most informative designs

- **N1 and N2 are the core robustness results.** Present each as a heatmap of cell accuracy and of wrong-X calls (the
  Fig. 6 layout). Mark where accuracy falls below 0.9 and where wrong-X calls exceed 1%.
- **N3 completes C4.** The skew test showed no effect only because the standard regime is easy. Repeat it under the
  hard regimes, where a 99:1 minority (19 cells) is most at risk.
- **N4 is required before claiming real-data behaviour.** The current simulator draws SNVs uniformly, gives one SNV per
  read, and rarely duplicates UMIs. Real 10x data have SNVs clustered in a few highly expressed genes, which changes
  how SNVs are co-observed and therefore phasing.
- **N13 is the strongest validation**, because it uses real reads, real noise and real expression with known truth.
  It needs two male samples from the same platform and their genotypes, or genotypes called from the samples.
- **N11 costs nothing:** `*_bc2xci.txt` already stores the haplotype UMI counts, so other thresholds can be applied
  after the run.

### Suggested order

N1 → N2 → N3 (simulator unchanged) → N4 → N5/N6 → N10/N11 → N7–N9 → N12 → N13/N14.
