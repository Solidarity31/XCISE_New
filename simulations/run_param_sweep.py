#!/usr/bin/env python3
"""
Parameter sweep for EXCISE 1.5: min UMIs per SNV (-u) x discordant penalty (-p).

Usage (from the repository root, after `make -C cpp`):
  python3 simulations/run_param_sweep.py [--reps 5] [--jobs 4]

The runs use cpp/excise15, which implements excise1.5.pl exactly (byte-identical
outputs to the Perl script when both reach the same optimum); a few grid points
are re-run with excise1.5.pl itself as a check (--perl-check).

Regimes (2,000 cells; 5% Both; X1:X2 = 50:50 among monoallelic cells; 5% escape SNVs):
  standard      1,000 SNVs, 40 reads/cell, 2% noise   (= skew_indep 50:50 datasets, reused)
  low_coverage  1,000 SNVs, 10 reads/cell, 2% noise
  high_noise    1,000 SNVs, 40 reads/cell, 10% noise
  sparse_snvs   3,000 SNVs, 15 reads/cell, 2% noise
Grid: -u {1,2,4,6,10,15,20} x -p {1,2,3,5,8,12,20}; all other flags fixed:
  -t 100 -j <jobs> -seed <dataset seed> -pretrain
(excise1.5.pl resets -u <= 0 to 10 and -p < 1 to 5, so the grid starts at 1.)

Outputs:
  simulations/data/params/<regime>/rep<k>/            simulated data (+ runs/)
  simulations/results/param_sweep_replicates.tsv      one row per run
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tests"))
from evaluate_xci_sim import evaluate  # noqa: E402

U_GRID = [1, 2, 4, 6, 10, 15, 20]
P_GRID = [1, 2, 3, 5, 8, 12, 20]
REGIMES = {
    # name: (n_snvs, reads_per_cell, noise, seed base)
    "standard": (1000, 40, 0.02, None),        # reuse skew_indep x1_50
    "low_coverage": (1000, 10, 0.02, 20000),
    "high_noise": (1000, 40, 0.10, 30000),
    "sparse_snvs": (3000, 15, 0.02, 40000),
}


def dataset(regime, rep):
    n_snvs, rpc, noise, base = REGIMES[regime]
    if base is None:
        d = os.path.join(HERE, "data", "skew_indep", "x1_50", "rep%d" % rep)
        return d, 10000 + 100 * 50 + rep
    seed = base + rep
    d = os.path.join(HERE, "data", "params", regime, "rep%d" % rep)
    if not os.path.exists(os.path.join(d, "sim.bam")):
        subprocess.check_call(
            [sys.executable, os.path.join(ROOT, "tests", "simulate_xci_data.py"), d, "2000", str(n_snvs),
             str(rpc), str(seed), "--x1-frac", "0.5", "--both-frac", "0.05", "--assign", "exact",
             "--noise", str(noise), "--escape-frac", "0.05"], stdout=subprocess.DEVNULL)
    return d, seed


def wrong_x(sim_dir, prefix, switched):
    """True X1/X2 cells called as the opposite X (after the global label switch)."""
    truth = dict(l.split()[:2] for l in open(os.path.join(sim_dir, "truth.tsv")))
    n = 0
    for l in open(prefix + "_bc2xci.txt"):
        cb, _, _, lab = l.rstrip("\n").split("\t")
        if switched and lab in ("X1", "X2"):
            lab = "X2" if lab == "X1" else "X1"
        if truth.get(cb) in ("X1", "X2") and lab in ("X1", "X2") and lab != truth[cb]:
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--tries", type=int, default=100)
    a = ap.parse_args()
    exe = os.path.join(ROOT, "cpp", "excise15")
    if not os.path.exists(exe):
        sys.exit("build the C++ port first: make -C cpp")
    out = os.path.join(HERE, "results", "param_sweep_replicates.tsv")
    done = set()
    rows = []
    if os.path.exists(out):
        with open(out) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for l in fh:
                r = dict(zip(hdr, l.rstrip("\n").split("\t")))
                rows.append(r)
                done.add((r["regime"], r["rep"], r["min_umis"], r["penalty"]))
    keys = None
    for regime in REGIMES:
        for rep in range(1, a.reps + 1):
            d, seed = dataset(regime, rep)
            for u in U_GRID:
                for p in P_GRID:
                    if (regime, str(rep), str(u), str(p)) in done:
                        continue
                    run = os.path.join(d, "params", "u%d_p%d" % (u, p))
                    os.makedirs(run, exist_ok=True)
                    t0 = time.time()
                    with open(os.path.join(run, "log"), "w") as lf:
                        subprocess.check_call(
                            [exe, "-o", "run", "-s", os.path.join(d, "sim.vcf"), "-b", os.path.join(d, "sim.bam"),
                             "-t", str(a.tries), "-j", str(a.jobs), "-seed", str(seed), "-pretrain",
                             "-u", str(u), "-p", str(p)], cwd=run, stderr=lf)
                    wall = time.time() - t0
                    prefix = os.path.join(run, "run_chrX_XCISE")
                    m = evaluate(d, prefix)
                    r = dict(regime=regime, rep=rep, seed=seed, min_umis=u, penalty=p, wall_s="%.2f" % wall)
                    r.update({k: ("%.6g" % v) if isinstance(v, float) else v for k, v in m.items()})
                    r["wrong_x_cells"] = wrong_x(d, prefix, bool(m["label_switched"]))
                    rows.append(r)
                    keys = keys or list(r.keys())
            print("%s rep%d done (%d runs total)" % (regime, rep, len(rows)), flush=True)
            keys = list(rows[-1].keys())
            with open(out, "w") as fh:
                fh.write("\t".join(keys) + "\n")
                for x in rows:
                    fh.write("\t".join(str(x.get(k, "NA")) for k in keys) + "\n")


if __name__ == "__main__":
    main()
