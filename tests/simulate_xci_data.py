#!/usr/bin/env python3
"""
Simulate a small WASP-tagged scRNA-seq BAM + het-SNV VCF on chrX with known
per-cell XCI truth, for testing / benchmarking excise1.5.pl vs excise15.py.

Usage:
  python3 tests/simulate_xci_data.py <outdir> [n_cells] [n_snvs] [reads_per_cell] [seed]
      [--x1-frac F] [--both-frac F] [--assign exact|bernoulli]
      [--noise F] [--escape-frac F] [--snvs-per-cell MIN MAX]

Cell states:
  X1   : haplotype A is the active X (cell expresses hapA alleles)
  X2   : haplotype B is the active X
  Both : biallelic cell (each read picks an allele at random)

Proportions are given on two different bases:
  --both-frac  fraction of ALL cells that are Both
  --x1-frac    fraction of the MONOALLELIC (X1 + X2) cells that are X1,
               i.e. the X1:X2 ratio. --x1-frac 0.99 means X1:X2 = 99:1
               among monoallelic cells; with --both-frac 0.05 that is
               94.05% X1, 0.95% X2 and 5% Both of all cells.
  --assign exact      exact counts: n_both = round(N * both_frac),
                      n_x1 = round((N - n_both) * x1_frac), states shuffled
  --assign bernoulli  each cell drawn independently with those probabilities

Without --x1-frac / --both-frac the legacy model is used (each cell
independently X1 40% / X2 55% / Both 5%), and output is identical to
the previous version of this script for the same positional arguments.

Writes <outdir>/sim.bam(.bai), sim.vcf, truth.tsv (barcode, state),
truth_snvs.tsv (pos, hapA_allele [1=REF, 2=ALT], escape) and sim_params.tsv.
Requires samtools on PATH.
"""
import argparse
import os
import random
import subprocess
import sys


def parse():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("n_cells", nargs="?", type=int, default=500)
    ap.add_argument("n_snvs", nargs="?", type=int, default=300)
    ap.add_argument("reads_per_cell", nargs="?", type=int, default=40)
    ap.add_argument("seed", nargs="?", type=int, default=1)
    ap.add_argument("--x1-frac", type=float, default=None,
                    help="X1 share among monoallelic (X1+X2) cells")
    ap.add_argument("--both-frac", type=float, default=None,
                    help="Both share among all cells")
    ap.add_argument("--assign", choices=("exact", "bernoulli"), default="exact")
    ap.add_argument("--noise", type=float, default=0.02,
                    help="per-read probability of observing the inactive allele")
    ap.add_argument("--escape-frac", type=float, default=None,
                    help="fraction of SNVs that escape XCI (default: n_snvs // 20)")
    ap.add_argument("--snvs-per-cell", type=int, nargs=2, default=(3, 25),
                    metavar=("MIN", "MAX"),
                    help="each cell expresses a random number of SNVs in [MIN, MAX]")
    a = ap.parse_args()
    if a.x1_frac is not None and not 0 <= a.x1_frac <= 1:
        ap.error("--x1-frac must be in [0, 1]")
    if a.both_frac is not None and not 0 <= a.both_frac <= 1:
        ap.error("--both-frac must be in [0, 1]")
    return a


def assign_states(a, n_cells):
    """State list for configured proportions (None = legacy per-cell draw)."""
    if a.x1_frac is None and a.both_frac is None:
        return None
    x1 = 0.5 if a.x1_frac is None else a.x1_frac
    both = 0.05 if a.both_frac is None else a.both_frac
    # separate RNG so the read-generation stream stays independent of states
    srnd = random.Random(a.seed * 7919 + 17)
    if a.assign == "exact":
        n_both = int(round(n_cells * both))
        n_x1 = int(round((n_cells - n_both) * x1))
        states = ["Both"] * n_both + ["X1"] * n_x1 + ["X2"] * (n_cells - n_both - n_x1)
        srnd.shuffle(states)
        return states
    states = []
    for _ in range(n_cells):
        if srnd.random() < both:
            states.append("Both")
        else:
            states.append("X1" if srnd.random() < x1 else "X2")
    return states


def main():
    a = parse()
    outdir, n_cells, n_snvs = a.outdir, a.n_cells, a.n_snvs
    reads_per_cell, seed = a.reads_per_cell, a.seed
    rnd = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)

    chrom_len = 5_000_000
    positions = sorted(rnd.sample(range(1000, chrom_len - 1000), n_snvs))
    hapA = {p: rnd.choice((1, 2)) for p in positions}       # allele on haplotype A
    n_escape = (max(1, n_snvs // 20) if a.escape_frac is None
                else int(round(n_snvs * a.escape_frac)))
    escape = set(rnd.sample(positions, n_escape))           # biallelic everywhere
    noise = a.noise
    smin, smax = a.snvs_per_cell
    preset = assign_states(a, n_cells)

    with open(os.path.join(outdir, "sim.vcf"), "w") as fv:
        fv.write("##fileformat=VCFv4.2\n")
        fv.write("##contig=<ID=X,length=%d>\n" % chrom_len)
        fv.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for p in positions:
            fv.write("X\t%d\trs%d\tA\tG\t.\tPASS\t.\tGT\t0|1\n" % (p, p))

    truth = []
    sam = os.path.join(outdir, "sim.sam")
    bases = "ACGT"
    with open(sam, "w") as fs:
        fs.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:X\tLN:%d\n" % chrom_len)
        rid = 0
        for c in range(n_cells):
            cb = "".join(rnd.choice(bases) for _ in range(16)) + "-1"
            if preset is None:
                r = rnd.random()
                state = "X1" if r < 0.4 else ("X2" if r < 0.95 else "Both")
            else:
                state = preset[c]
            truth.append((cb, state))
            # each cell expresses a subset of SNV-bearing genes
            cell_snvs = rnd.sample(positions, min(len(positions), rnd.randint(smin, smax)))
            for _ in range(rnd.randint(reads_per_cell // 2, reads_per_cell * 3 // 2)):
                p = rnd.choice(cell_snvs)
                if p in escape or state == "Both":
                    allele = rnd.choice((1, 2))
                else:
                    al = hapA[p] if state == "X1" else 3 - hapA[p]
                    allele = al if rnd.random() > noise else 3 - al
                umi = "".join(rnd.choice(bases) for _ in range(10))
                start = p - rnd.randint(1, 40)
                seq = "".join(rnd.choice(bases) for _ in range(50))
                rid += 1
                fs.write("r%d\t0\tX\t%d\t60\t50M\t*\t0\t0\t%s\t%s\t"
                         "CB:Z:%s\tUB:Z:%s\tvA:B:c,%d\tvG:B:i,%d\tvW:i:1\n"
                         % (rid, start, seq, "I" * 50, cb, umi, allele, p - 1))
    bam = os.path.join(outdir, "sim.bam")
    subprocess.check_call(["samtools", "sort", "-o", bam, sam])
    subprocess.check_call(["samtools", "index", bam])
    os.remove(sam)
    with open(os.path.join(outdir, "truth.tsv"), "w") as ft:
        for cb, s in truth:
            ft.write("%s\t%s\n" % (cb, s))
    with open(os.path.join(outdir, "truth_snvs.tsv"), "w") as ft:
        ft.write("pos\thapA_allele\tescape\n")
        for p in positions:
            ft.write("%d\t%d\t%d\n" % (p, hapA[p], int(p in escape)))
    counts = {s: sum(1 for _, t in truth if t == s) for s in ("X1", "X2", "Both")}
    with open(os.path.join(outdir, "sim_params.tsv"), "w") as fp:
        params = [("n_cells", n_cells), ("n_snvs", n_snvs),
                  ("reads_per_cell", reads_per_cell), ("seed", seed),
                  ("model", "legacy" if preset is None else a.assign),
                  ("x1_frac_of_mono", a.x1_frac), ("both_frac_of_all", a.both_frac),
                  ("noise", noise), ("n_escape", n_escape),
                  ("snvs_per_cell_min", smin), ("snvs_per_cell_max", smax),
                  ("n_reads", rid), ("n_X1", counts["X1"]), ("n_X2", counts["X2"]),
                  ("n_Both", counts["Both"])]
        for k, v in params:
            fp.write("%s\t%s\n" % (k, "NA" if v is None else v))
    print("wrote", bam, "with", rid, "reads,", n_cells, "cells,", n_snvs, "SNVs;",
          "X1/X2/Both = %d/%d/%d" % (counts["X1"], counts["X2"], counts["Both"]))


if __name__ == "__main__":
    main()
