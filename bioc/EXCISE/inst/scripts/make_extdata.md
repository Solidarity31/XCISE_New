# Example data in inst/extdata

The files are simulated; they contain no real sequencing data.

- `sim_chrX.bam`, `sim_chrX.bam.bai`: 300 cells with WASP-style tags
  (`vA`, `vG`, `vW`), `CB` and `UB`, on a simulated 5 Mb chromosome `X`.
  X1:X2 = 40:60 among monoallelic cells, 5% biallelic (`Both`) cells,
  2% allelic noise, 5% escape SNVs.
- `sim_chrX_hets.vcf`: the 200 heterozygous SNVs.
- `sim_chrX_truth.tsv`: true state of each cell (barcode, X1/X2/Both).

They were generated with `tests/simulate_xci_data.py` from the EXCISE
repository (https://github.com/Solidarity31/EXCISE):

```
python3 tests/simulate_xci_data.py out 300 200 30 2026 --x1-frac 0.4 --both-frac 0.05
```

and renamed (`sim.bam` -> `sim_chrX.bam`, `sim.vcf` -> `sim_chrX_hets.vcf`,
`truth.tsv` -> `sim_chrX_truth.tsv`).
