#!/usr/bin/env python3
"""
Runtime at the default number of tries (-t 100): original XCISE xcise.pl vs
excise1.5.pl (Perl) vs excise15.py (Python) vs cpp/excise15 (C++).

Usage (from the repository root):
  python3 simulations/run_default_tries_benchmark.py --xcise /path/to/xcise.pl \
      [--datasets c1000_s11 ...] [--tries 100] [--jobs 4]

Datasets are the speed-benchmark datasets in simulations/data/speed/.
  original : perl xcise.pl -t <tries> -u 4      (no -j / -seed / -pretrain;
                                                 tries run one after another)
  EXCISE   : -t <tries> -j <jobs> -seed <seed> -pretrain   (Perl, Python, C++)
-u 4 gives the original the same SNV filter as the EXCISE defaults. Tools run
one at a time (never concurrently), fastest first; finished runs are reused.
Output: simulations/results/default_tries_benchmark.tsv
"""
import argparse
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from run_speed_benchmark import cell_accuracy  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xcise", required=True, help="path to the original XCISE xcise.pl")
    ap.add_argument("--datasets", nargs="+", default=["c1000_s11"])
    ap.add_argument("--tries", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    out = os.path.join(HERE, "results", "default_tries_benchmark.tsv")
    done = {}
    if os.path.exists(out):
        with open(out) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for l in fh:
                r = dict(zip(hdr, l.rstrip("\n").split("\t")))
                done[(r["dataset"], r["tool"], r["tries"])] = r
    rows = list(done.values())
    for ds in a.datasets:
        d = os.path.join(HERE, "data", "speed", ds)
        seed = ds.split("_s")[-1]
        common = ["-o", "b", "-s", "../../sim.vcf", "-b", "../../sim.bam", "-t", str(a.tries)]
        excise = common + ["-j", str(a.jobs), "-seed", seed, "-pretrain"]
        tools = [
            ("cpp", [os.path.join(ROOT, "cpp", "excise15")] + excise),
            ("python", [sys.executable, os.path.join(ROOT, "excise15.py")] + excise),
            ("perl1.5", ["perl", os.path.join(ROOT, "excise1.5.pl")] + excise),
            ("original_xcise", ["perl", os.path.abspath(a.xcise)] + common + ["-u", "4"]),
        ]
        truth = dict(l.split()[:2] for l in open(os.path.join(d, "truth.tsv")))
        for name, cmd in tools:
            if (ds, name, str(a.tries)) in done:
                continue
            run = os.path.join(d, "t%d" % a.tries, name)
            os.makedirs(run, exist_ok=True)
            t0 = time.time()
            with open(os.path.join(run, "log"), "w") as lf:
                subprocess.check_call(cmd, cwd=run, stderr=lf)
            wall = time.time() - t0
            with open(os.path.join(run, "b_chrX_XCISE_summary.txt")) as fh:
                best = re.match(r"Best score\s+:\s+(\S+)", fh.readline()).group(1)
            lab = dict((l.split("\t")[0], l.rstrip("\n").split("\t")[3])
                       for l in open(os.path.join(run, "b_chrX_XCISE_bc2xci.txt")))
            r = dict(dataset=ds, n_cells=ds.split("_")[0][1:], seed=seed, tries=a.tries, tool=name,
                     wall_s="%.2f" % wall, best_score=best, cell_acc="%.4f" % cell_accuracy(truth, lab))
            rows.append(r)
            print("%s %-15s t=%d  %.1f s  best=%s  cell_acc=%s" % (ds, name, a.tries, wall, best, r["cell_acc"]),
                  flush=True)
            keys = list(rows[0].keys())
            with open(out, "w") as fh:
                fh.write("\t".join(keys) + "\n")
                for x in rows:
                    fh.write("\t".join(str(x[k]) for k in keys) + "\n")


if __name__ == "__main__":
    main()
