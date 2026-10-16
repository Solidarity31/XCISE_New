#!/usr/bin/env python3
"""
Simulate a small WASP-tagged scRNA-seq BAM + het-SNV VCF on chrX with known
per-cell XCI truth, for testing / benchmarking excise1.5.pl vs excise15.py.

Usage:
  python3 tests/simulate_xci_data.py <outdir> [n_cells] [n_snvs] [reads_per_cell] [seed]

Writes <outdir>/sim.bam(.bai), <outdir>/sim.vcf and <outdir>/truth.tsv
(barcode, true_state with X1 = haplotype A active, X2 = haplotype B active,
Both = biallelic cell). Requires samtools on PATH.
"""
import os
import random
import subprocess
import sys


def main():
    outdir = sys.argv[1]
    n_cells = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    n_snvs = int(sys.argv[3]) if len(sys.argv) > 3 else 300
    reads_per_cell = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    rnd = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)

    chrom_len = 5_000_000
    positions = sorted(rnd.sample(range(1000, chrom_len - 1000), n_snvs))
    hapA = {p: rnd.choice((1, 2)) for p in positions}       # allele on haplotype A
    escape = set(rnd.sample(positions, max(1, n_snvs // 20)))  # biallelic everywhere
    noise = 0.02

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
            r = rnd.random()
            state = "X1" if r < 0.4 else ("X2" if r < 0.95 else "Both")
            truth.append((cb, state))
            # each cell expresses a subset of SNV-bearing genes
            cell_snvs = rnd.sample(positions, min(len(positions), rnd.randint(3, 25)))
            for _ in range(rnd.randint(reads_per_cell // 2, reads_per_cell * 3 // 2)):
                p = rnd.choice(cell_snvs)
                if p in escape or state == "Both":
                    allele = rnd.choice((1, 2))
                else:
                    a = hapA[p] if state == "X1" else 3 - hapA[p]
                    allele = a if rnd.random() > noise else 3 - a
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
    print("wrote", bam, "with", rid, "reads,", n_cells, "cells,", n_snvs, "SNVs")


if __name__ == "__main__":
    main()
