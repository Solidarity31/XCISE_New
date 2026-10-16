#!/usr/bin/env python3
"""
XCI-skew sweep: simulate datasets at several X1:X2 ratios, run excise15.py
on each, and score the result against the simulated truth.

Usage (from the repository root):
  python3 simulations/run_xci_skew_sweep.py [--reps 10] [--jobs 4]

Fixed across conditions (only the X1:X2 ratio changes):
  2000 cells, 1000 SNVs, 40 reads/cell (uniform 20-60), 3-25 expressed
  SNVs per cell, allelic noise 0.02, 5% escape SNVs, Both = 5% of ALL cells,
  exact state counts. Replicate k uses simulation seed 100+k for every
  ratio (common random numbers) and the same EXCISE -seed.
EXCISE settings: -t 100 -j <jobs> -seed <seed> -pretrain (README defaults).

Existing simulations / runs are reused, never overwritten.
Outputs:
  simulations/data/skew/x1_<pct>/rep<k>/        simulated data + EXCISE run
  simulations/results/skew_replicates.tsv       one row per run
  simulations/results/skew_snv_level.tsv.gz     per-SNV truth vs inference
"""
import argparse
import gzip
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tests"))
from evaluate_xci_sim import evaluate, snv_table  # noqa: E402

RATIOS = [(50, 50), (55, 45), (70, 30), (90, 10), (95, 5), (99, 1)]
SIM = dict(n_cells=2000, n_snvs=1000, reads_per_cell=40, both_frac=0.05,
           noise=0.02, escape_frac=0.05, spc=(3, 25))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--tries", type=int, default=100)
    a = ap.parse_args()
    res_dir = os.path.join(HERE, "results")
    os.makedirs(res_dir, exist_ok=True)
    rows, snv_rows = [], []
    for x1, x2 in RATIOS:
        for rep in range(1, a.reps + 1):
            seed = 100 + rep
            d = os.path.join(HERE, "data", "skew", "x1_%02d" % x1, "rep%d" % rep)
            if not os.path.exists(os.path.join(d, "sim.bam")):
                subprocess.check_call(
                    [sys.executable, os.path.join(ROOT, "tests", "simulate_xci_data.py"), d,
                     str(SIM["n_cells"]), str(SIM["n_snvs"]), str(SIM["reads_per_cell"]), str(seed),
                     "--x1-frac", str(x1 / 100.0), "--both-frac", str(SIM["both_frac"]),
                     "--assign", "exact", "--noise", str(SIM["noise"]),
                     "--escape-frac", str(SIM["escape_frac"]),
                     "--snvs-per-cell", str(SIM["spc"][0]), str(SIM["spc"][1])],
                    stdout=subprocess.DEVNULL)
            run = os.path.join(d, "excise15")
            prefix = os.path.join(run, "run_chrX_XCISE")
            log = os.path.join(run, "excise15.log")
            timing = os.path.join(run, "wall_seconds.txt")
            if not os.path.exists(prefix + "_summary.txt"):
                os.makedirs(run, exist_ok=True)
                t0 = time.time()
                with open(log, "w") as lf:
                    subprocess.check_call(
                        [sys.executable, os.path.join(ROOT, "excise15.py"), "-o", "run",
                         "-s", "../sim.vcf", "-b", "../sim.bam", "-t", str(a.tries),
                         "-j", str(a.jobs), "-seed", str(seed), "-pretrain"],
                        cwd=run, stderr=lf)
                with open(timing, "w") as tf:
                    tf.write("%.3f\n" % (time.time() - t0))
            wall = float(open(timing).read())
            opt = re.search(r"Optimiser wall time: ([\d.]+)", open(log).read())
            m = evaluate(d, prefix)
            row = dict(ratio="%d:%d" % (x1, x2), x1_pct_of_mono=x1, rep=rep, seed=seed,
                       wall_s=wall, optimiser_s=float(opt.group(1)) if opt else float("nan"))
            row.update(m)
            row["minority_pct_of_all"] = 100.0 * m["n_minority"] / m["n_cells"]
            rows.append(row)
            for r in snv_table(d, prefix, bool(m["label_switched"])):
                snv_rows.append((row["ratio"], rep) + r)
            print("%-6s rep%-2d  hap_acc=%.4f cell_acc=%.4f minority %d/%d recall=%.3f  wall=%.1fs"
                  % (row["ratio"], rep, m["hap_acc"], m["cell_acc"], m["minority_correct"],
                     m["n_minority"], m["minority_recall"], wall), flush=True)

    keys = list(rows[0].keys())
    with open(os.path.join(res_dir, "skew_replicates.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(("%.6g" % r[k]) if isinstance(r[k], float) else str(r[k])
                               for k in keys) + "\n")
    with gzip.open(os.path.join(res_dir, "skew_snv_level.tsv.gz"), "wt") as fh:
        fh.write("ratio\trep\tpos\ttrue_x1_allele\tescape\tinferred_x1_allele\tmono\tdepth\tstatus\n")
        for r in snv_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


if __name__ == "__main__":
    main()
