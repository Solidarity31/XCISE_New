#!/usr/bin/env python3
"""
Runtime benchmark: excise1.5.pl (Perl) vs excise15.py (Python/numba) vs
cpp/excise15 (C++/htslib) on simulated datasets of increasing cell number.

Usage (from the repository root):
  python3 simulations/run_speed_benchmark.py [--cells 1000 2000 4000 8000 16000]
                                             [--seeds 11 12 13] [--jobs 4]

Datasets: 2000 SNVs, 60 reads/cell, legacy simulator state model
(X1 40% / X2 55% / Both 5%), seed = replicate seed. Both tools run with
-t 4 -j <jobs> -seed <seed> -pretrain, one after the other (never
concurrently), and wall-clock time is recorded for the full run including
BAM reading; the Python time includes its one-off numba JIT compilation.
Whether the two bc2xci outputs are byte-identical is recorded per run.
If cpp/excise15 has been built (make -C cpp), the C++ port is timed on the
same datasets with the same flags, after the Perl/Python runs; rows that
already have a C++ time are not re-run.
With --xcise <path to the original XCISE xcise.pl>, the original is timed
too, with -t 4 -u 4 (same tries and SNV filter as the EXCISE runs). It has
no -j, -seed or -pretrain: its tries run one after another on one core.
Cell accuracy against the simulated truth is recorded for every tool
(one global X1/X2 label switch allowed, since the names are arbitrary).

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
    ap.add_argument("--xcise", default=None, help="path to the original XCISE xcise.pl")
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
            write_rows(out, rows)
    # C++ port, timed separately (never concurrently with anything else)
    cpp = os.path.join(ROOT, "cpp", "excise15")
    if os.path.exists(cpp):
        for r in rows:
            if r.get("cpp_s", "NA") not in ("NA", ""):
                continue
            d = os.path.join(HERE, "data", "speed", "c%s_s%s" % (r["n_cells"], r["seed"]))
            os.makedirs(os.path.join(d, "cpp"), exist_ok=True)
            tc = timed([cpp, "-o", "b", "-s", "../sim.vcf", "-b", "../sim.bam", "-t", "4", "-j", str(a.jobs),
                        "-seed", str(r["seed"]), "-pretrain"], os.path.join(d, "cpp"), os.path.join(d, "cpp", "log"))
            opt = re.search(r"Optimiser wall time: ([\d.]+)", open(os.path.join(d, "cpp", "log")).read())
            r["cpp_s"] = "%.2f" % tc
            r["cpp_optimiser_s"] = opt.group(1) if opt else "NA"
            r["cpp_speedup_vs_perl"] = "%.2f" % (float(r["perl_s"]) / tc)
            r["cpp_speedup_vs_python"] = "%.2f" % (float(r["python_s"]) / tc)
            r["cpp_bc2xci_identical_perl"] = int(filecmp.cmp(
                os.path.join(d, "perl", "b_chrX_XCISE_bc2xci.txt"),
                os.path.join(d, "cpp", "b_chrX_XCISE_bc2xci.txt"), shallow=False))
            print("cells=%s seed=%s  cpp=%.2fs  (perl %ss, python %ss)  identical_to_perl=%d"
                  % (r["n_cells"], r["seed"], tc, r["perl_s"], r["python_s"], r["cpp_bc2xci_identical_perl"]),
                  flush=True)
            write_rows(out, rows)
    # original XCISE (single process), timed separately
    if a.xcise:
        for r in rows:
            if r.get("xcise_s", "NA") not in ("NA", ""):
                continue
            d = os.path.join(HERE, "data", "speed", "c%s_s%s" % (r["n_cells"], r["seed"]))
            os.makedirs(os.path.join(d, "xcise_orig"), exist_ok=True)
            tx = timed(["perl", os.path.abspath(a.xcise), "-o", "b", "-s", "../sim.vcf", "-b", "../sim.bam",
                        "-t", "4", "-u", "4"], os.path.join(d, "xcise_orig"), os.path.join(d, "xcise_orig", "log"))
            r["xcise_s"] = "%.2f" % tx
            r["perl15_speedup_vs_xcise"] = "%.2f" % (tx / float(r["perl_s"]))
            if r.get("cpp_s", "NA") not in ("NA", ""):
                r["cpp_speedup_vs_xcise"] = "%.2f" % (tx / float(r["cpp_s"]))
            print("cells=%s seed=%s  original xcise=%.1fs  (perl1.5 %ss, cpp %ss)"
                  % (r["n_cells"], r["seed"], tx, r["perl_s"], r.get("cpp_s")), flush=True)
            write_rows(out, rows)
    # Output agreement beyond byte identity: both tools are stochastic, so on
    # sparse data their best tries can reach different local optima.
    for r in rows:
        d = os.path.join(HERE, "data", "speed", "c%s_s%s" % (r["n_cells"], r["seed"]))
        best, lab = {}, {}
        subs = [x for x in ("perl", "python", "cpp", "xcise_orig")
                if os.path.exists(os.path.join(d, x, "b_chrX_XCISE_summary.txt"))]
        for sub in subs:
            with open(os.path.join(d, sub, "b_chrX_XCISE_summary.txt")) as fh:
                best[sub] = re.match(r"Best score\s+:\s+(\S+)", fh.readline()).group(1)
            with open(os.path.join(d, sub, "b_chrX_XCISE_bc2xci.txt")) as fh:
                lab[sub] = dict((l.split("\t")[0], l.rstrip("\n").split("\t")[3]) for l in fh)
        r["perl_best_score"], r["python_best_score"] = best["perl"], best["python"]
        r["cells_label_differs"] = sum(1 for cb in lab["perl"] if lab["python"].get(cb) != lab["perl"][cb])
        if "cpp" in best:
            r["cpp_best_score"] = best["cpp"]
            r["cpp_cells_label_differs_perl"] = sum(1 for cb in lab["perl"] if lab["cpp"].get(cb) != lab["perl"][cb])
        if "xcise_orig" in best:
            r["xcise_best_score"] = best["xcise_orig"]
            r["xcise_cells_label_differs_perl"] = sum(1 for cb in lab["perl"] if lab["xcise_orig"].get(cb) != lab["perl"][cb])
        truth = os.path.join(d, "truth.tsv")
        if os.path.exists(truth):
            t = dict(l.split()[:2] for l in open(truth))
            for sub in subs:
                r[{"perl": "perl", "python": "python", "cpp": "cpp", "xcise_orig": "xcise"}[sub] + "_cell_acc"] = \
                    "%.4f" % cell_accuracy(t, lab[sub])
    write_rows(out, rows)


def cell_accuracy(truth, called):
    """Fraction of true cells called correctly, allowing one global X1/X2 swap."""
    sw = {"X1": "X2", "X2": "X1"}
    a = sum(1 for cb, s in truth.items() if called.get(cb) == s)
    b = sum(1 for cb, s in truth.items() if sw.get(called.get(cb), called.get(cb)) == s)
    return max(a, b) / len(truth)


def write_rows(out, rows):
    keys = list(rows[0].keys())
    for x in rows:
        keys += [k for k in x if k not in keys]
    with open(out, "w") as fh:  # rewritten after every run (resumable)
        fh.write("\t".join(keys) + "\n")
        for x in rows:
            fh.write("\t".join(str(x.get(k, "NA")) for k in keys) + "\n")


if __name__ == "__main__":
    main()
