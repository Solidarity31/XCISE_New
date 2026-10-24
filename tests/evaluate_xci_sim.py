#!/usr/bin/env python3
"""
Score an EXCISE run on a simulated dataset from tests/simulate_xci_data.py.

Usage:
  python3 tests/evaluate_xci_sim.py <sim_dir> <excise_prefix>
    e.g. <excise_prefix> = runs/rep1/out_chrX_XCISE  (the run's file prefix)

EXCISE's X1/X2 names are arbitrary (it relabels so that #X2 >= #X1), so one
global label switch is chosen from the SNV haplotypes: if most phased SNVs
put the hapB allele on EXCISE-X1, EXCISE-X1 is mapped to truth-X2 (and
vice versa) for every SNV and cell. The same mapping is used for all metrics.

Metrics (cells = all simulated cells; truth states X1/X2/Both):
  hap_acc                fraction of phased (X1A != Unk) non-escape SNVs whose
                         orientation is correct after the global switch
  snv_unassigned         fraction of SNVs passing EXCISE filters left Unk
  snv_unassigned_nonesc  same, non-escape SNVs only
  escape_unassigned      fraction of escape SNVs left Unk (ideal: 1)
  cell_acc               fraction of all cells whose call equals truth
                         (Low_coverage / Unknown count as wrong)
  cell_acc_informative   same, among cells called X1/X2/Both
  informative_recovery   fraction of all cells called X1, X2 or Both
  mono_recovery          fraction of true X1/X2 cells called X1 or X2
  recall_X1/X2/Both      per-class recall; macro_recall = their mean
  minority_*             the smaller true monoallelic class (X2 on ties)
"""
import os
import sys


def read_kv(path):
    with open(path) as fh:
        return dict(l.rstrip("\n").split("\t", 1) for l in fh if l.strip())


def evaluate(sim_dir, prefix):
    truth = dict(l.split() for l in open(os.path.join(sim_dir, "truth.tsv")))
    snv_truth = {}
    with open(os.path.join(sim_dir, "truth_snvs.tsv")) as fh:
        next(fh)
        for l in fh:
            p, a, e = l.split()
            snv_truth[int(p)] = (int(a), e == "1")

    # ---- SNV haplotypes from the EXCISE VCF ----
    n_pass = n_unk = n_ne = n_ne_unk = n_esc = n_esc_unk = 0
    agree = phased = 0
    with open(prefix + ".vcf") as fh:
        for l in fh:
            if l.startswith("#"):
                continue
            f = l.rstrip("\n").split("\t")
            pos = int(f[1])
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            x1a = info["X1A"]
            hapA, esc = snv_truth[pos]
            n_pass += 1
            unk = x1a == "Unk"
            n_unk += unk
            if esc:
                n_esc += 1; n_esc_unk += unk
                continue
            n_ne += 1; n_ne_unk += unk
            if not unk:
                phased += 1
                agree += (1 if x1a == "Ref" else 2) == hapA
    switched = phased > 0 and agree < phased / 2.0
    hap_acc = (max(agree, phased - agree) / phased) if phased else float("nan")

    # ---- cells ----
    pred = {}
    with open(prefix + "_bc2xci.txt") as fh:
        for l in fh:
            cb, _, _, lab = l.rstrip("\n").split("\t")
            if switched and lab in ("X1", "X2"):
                lab = "X2" if lab == "X1" else "X1"
            pred[cb] = lab
    cls = ("X1", "X2", "Both")
    n_true = {c: sum(1 for s in truth.values() if s == c) for c in cls}
    n_pred = {c: sum(1 for cb in truth if pred.get(cb) == c) for c in cls}
    hit = {c: sum(1 for cb, s in truth.items() if s == c and pred.get(cb) == c) for c in cls}
    n = len(truth)
    informative = [cb for cb in truth if pred.get(cb) in cls]
    mono_true = [cb for cb, s in truth.items() if s in ("X1", "X2")]
    rec = {c: (hit[c] / n_true[c] if n_true[c] else float("nan")) for c in cls}
    minority = "X1" if n_true["X1"] < n_true["X2"] else "X2"
    majority = "X2" if minority == "X1" else "X1"
    nan = float("nan")
    return {
        "n_cells": n, "n_X1": n_true["X1"], "n_X2": n_true["X2"], "n_Both": n_true["Both"],
        "minority_class": minority, "n_minority": n_true[minority],
        "minority_called": n_pred[minority], "minority_correct": hit[minority],
        "minority_recall": rec[minority],
        "minority_precision": (hit[minority] / n_pred[minority]) if n_pred[minority] else nan,
        "majority_recall": rec[majority],
        "recall_X1": rec["X1"], "recall_X2": rec["X2"], "recall_Both": rec["Both"],
        "macro_recall": sum(rec[c] for c in cls if n_true[c]) / sum(1 for c in cls if n_true[c]),
        "cell_acc": sum(1 for cb, s in truth.items() if pred.get(cb) == s) / n,
        "cell_acc_informative": (sum(1 for cb in informative if pred[cb] == truth[cb]) / len(informative))
                                if informative else nan,
        "informative_recovery": len(informative) / n,
        "mono_recovery": (sum(1 for cb in mono_true if pred.get(cb) in ("X1", "X2")) / len(mono_true))
                         if mono_true else nan,
        "est_minority_share_of_mono": (n_pred[minority] / (n_pred["X1"] + n_pred["X2"]))
                                      if (n_pred["X1"] + n_pred["X2"]) else nan,
        "hap_acc": hap_acc, "snvs_phased_nonesc": phased, "label_switched": int(switched),
        "snvs_pass_filter": n_pass,
        "snv_unassigned": n_unk / n_pass if n_pass else nan,
        "snv_unassigned_nonesc": n_ne_unk / n_ne if n_ne else nan,
        "escape_unassigned": n_esc_unk / n_esc if n_esc else nan,
    }


def snv_table(sim_dir, prefix, switched):
    """Per-SNV truth vs inference after the global label switch.

    status: correct / incorrect / unassigned (X1A=Unk) / filtered (SNV not
    in the EXCISE VCF, i.e. removed by the -u / -m filters).
    hapA_on_X1 (truth, in EXCISE orientation): 1 if the hapA allele is the
    allele EXCISE should place on X1 after the label switch.
    """
    calls = {}
    with open(prefix + ".vcf") as fh:
        for l in fh:
            if l.startswith("#"):
                continue
            f = l.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            calls[int(f[1])] = (info["X1A"], "MONO" in f[7].split(";"), int(f[5]))
    rows = []
    with open(os.path.join(sim_dir, "truth_snvs.tsv")) as fh:
        next(fh)
        for l in fh:
            p, a, e = l.split()
            p, a, esc = int(p), int(a), e == "1"
            # allele (1=REF, 2=ALT) on truth haplotype X1 (= hapA)
            true_x1 = a
            if p not in calls:
                rows.append((p, true_x1, int(esc), "NA", 0, 0, "filtered"))
                continue
            x1a, mono, dp = calls[p]
            if x1a == "Unk":
                inf = "NA"; st = "unassigned"
            else:
                inf = 1 if x1a == "Ref" else 2
                if switched:
                    inf = 3 - inf
                st = "escape_phased" if esc else ("correct" if inf == true_x1 else "incorrect")
            if esc and x1a == "Unk":
                st = "escape_unassigned"
            rows.append((p, true_x1, int(esc), inf, int(mono), dp, st))
    return rows


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    for k, v in evaluate(sys.argv[1], sys.argv[2]).items():
        print("%s\t%s" % (k, ("%.4f" % v) if isinstance(v, float) else v))
