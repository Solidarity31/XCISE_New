#!/usr/bin/env python3
"""
Runtime benchmark: excise1.5.pl (Perl) vs excise15.py (Python/numba)
on simulated datasets of increasing cell number.

Usage (from the repository root):
  python3 simulations/run_speed_benchmark.py [--cells 1000 2000 4000 8000 16000]
                                             [--seeds 11 12 13] [--jobs 4]

Datasets: 2000 SNVs, 60 reads/cell, legacy simulator state model
(X1 40% / X2 55% / Both 5%), seed = replicate seed. Both tools run with
-t 4 -j <jobs> -seed <seed> -pretrain, one after the other (never
concurrently), and wall-clock time is recorded for the full run including
BAM reading; the Python time includes its one-off numba JIT compilation.
Whether the two bc2xci outputs are byte-identical is recorded per run.

Existing datasets / runs are reused, never overwritten.
Output: simulations/results/speed_benchmark.tsv
"""
import argparse
import filecmp
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def timed(cmd, cwd, log):
    t0 = time.time()
    with open(log, "w") as lf:
        subprocess.check_call(cmd, cwd=cwd, stderr=lf)
    return time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, nargs="+", default=[1000, 2000, 4000, 8000, 16000])
    ap.add_argument("--seeds", type=int, nargs="+", default=[11, 12, 13])
    ap.add_argument("--snvs", type=int, default=2000)
    ap.add_argument("--reads", type=int, default=60)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--reuse", default=None,
                    help="existing dataset dir to use for (8000 cells, seed 11), "
                         "e.g. the dataset of the first reported benchmark")
    a = ap.parse_args()
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    out = os.path.join(HERE, "results", "speed_benchmark.tsv")
    done = {}
    if os.path.exists(out):
        with open(out) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for l in fh:
                r = dict(zip(hdr, l.rstrip("\n").split("\t")))
                done[(int(r["n_cells"]), int(r["seed"]))] = r
    rows = []
    for n in a.cells:
        for s in a.seeds:
            if (n, s) in done:
                rows.append(done[(n, s)]); continue
            d = os.path.join(HERE, "data", "speed", "c%d_s%d" % (n, s))
            os.makedirs(d, exist_ok=True)
            src = "simulated"
            if not os.path.exists(os.path.join(d, "sim.bam")):
                if a.reuse and n == 8000 and s == 11 and os.path.exists(os.path.join(a.reuse, "sim.bam")):
                    for f in ("sim.bam", "sim.bam.bai", "sim.vcf", "truth.tsv"):
                        os.symlink(os.path.join(a.reuse, f), os.path.join(d, f))
                    src = "reused:" + a.reuse
                else:
                    subprocess.check_call(
                        [sys.executable, os.path.join(ROOT, "tests", "simulate_xci_data.py"),
                         d, str(n), str(a.snvs), str(a.reads), str(s)], stdout=subprocess.DEVNULL)
            elif os.path.islink(os.path.join(d, "sim.bam")):
                src = "reused:" + os.path.dirname(os.path.realpath(os.path.join(d, "sim.bam")))
            nreads = int(subprocess.check_output(["samtools", "view", "-c", os.path.join(d, "sim.bam")]))
            common = ["-o", "b", "-s", "../sim.vcf", "-b", "../sim.bam", "-t", "4",
                      "-j", str(a.jobs), "-seed", str(s), "-pretrain"]
            for sub in ("perl", "python"):
                os.makedirs(os.path.join(d, sub), exist_ok=True)
            tpy = timed([sys.executable, os.path.join(ROOT, "excise15.py")] + common,
                        os.path.join(d, "python"), os.path.join(d, "python", "log"))
            tpl = timed(["perl", os.path.join(ROOT, "excise1.5.pl")] + common,
                        os.path.join(d, "perl"), os.path.join(d, "perl", "log"))
            opt = re.search(r"Optimiser wall time: ([\d.]+)",
                            open(os.path.join(d, "python", "log")).read())
            same = filecmp.cmp(os.path.join(d, "perl", "b_chrX_XCISE_bc2xci.txt"),
                               os.path.join(d, "python", "b_chrX_XCISE_bc2xci.txt"), shallow=False)
            r = dict(n_cells=n, n_snvs=a.snvs, seed=s, n_reads=nreads, perl_s="%.2f" % tpl,
                     python_s="%.2f" % tpy, python_optimiser_s=opt.group(1) if opt else "NA",
                     speedup="%.2f" % (tpl / tpy), bc2xci_identical=int(same), dataset=src)
            rows.append(r)
            print("cells=%d seed=%d  perl=%.1fs python=%.1fs  speedup=%.1fx identical=%s"
                  % (n, s, tpl, tpy, tpl / tpy, same), flush=True)
            with open(out, "w") as fh:  # rewrite after every run (resumable)
                keys = list(rows[0].keys())
                fh.write("\t".join(keys) + "\n")
                for x in rows:
                    fh.write("\t".join(str(x[k]) for k in keys) + "\n")


if __name__ == "__main__":
    main()
