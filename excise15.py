#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
excise15.py -- Python/numba port of excise1.5.pl (EXCISE v1.5).

Derived from xcise.pl of XCISE (https://github.com/Vityay/XCISE).
Please cite: Henning RH, Rust TM, Dijksterhuis K, Eggen BJL, Guryev V.
  bioRxiv (2024). doi:10.1101/2024.08.29.610317

Same algorithm, flags and output files as excise1.5.pl; it is faster because:
  1. Incremental (delta) scoring. excise1.5.pl rescans every cell to rescore
     after each single-SNV move; here only the cells covering that SNV are
     re-evaluated. The score is a sum of per-cell terms, so the result is
     identical, but each move costs O(cells at SNV) instead of O(all cells).
  2. The hill-climber, pretraining graph and BFS are compiled with numba.
  3. BAM/UMI data are packed into flat integer arrays (CSR) once, instead of
     nested hashes of "umi\\tcb" strings that are re-split on every move.
  4. Tries run in parallel worker processes (-j), as in the Perl script.

Outputs (identical names/formats to excise1.5.pl):
  <prefix>_chr<chrom>_XCISE_bc2xci.txt
  <prefix>_chr<chrom>_XCISE_summary.txt
  <prefix>_chr<chrom>_XCISE.vcf
  <prefix>_chr<chrom>_XCISE_haplotypes.tsv
  <prefix>_chr<chrom>_XCISE_cell_snv_profile.tsv

Requirements: python >= 3.8, numpy, numba, pysam
  pip install numpy numba pysam

Run results are not bit-identical to the Perl script for the same -seed:
the RNGs differ, and Perl iterates hashes in random order. On well-supported
data both converge to the same phasing.
"""

import gzip
import multiprocessing as mp
import os
import re
import sys
import time

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover
    sys.exit("numba is required: pip install numba")


def warn(*args):
    sys.stderr.write("".join(str(a) for a in args))
    sys.stderr.flush()


def perl_num(x):
    """Format a number the way Perl's default stringification does."""
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    x = float(x)
    if x == int(x) and abs(x) < 1e15:
        return str(int(x))
    return "%.15g" % x


# =============================================================================
# CLI (same flags and quirks as excise1.5.pl)
# =============================================================================
USAGE = ("Usage: python3 " + os.path.basename(sys.argv[0]) +
         " -o <prefix> -s <vcf> -b <bam> [...] [-r <chrom>] [-t <tries>]"
         " [-u <min_umis>] [-m <min_maf>] [-p <penalty>]"
         " [-j <jobs>] [-samthreads <N>] [-seed <int>]"
         " [-wl <barcodes.txt|tsv[.gz]>] [-wl_strip_gem]"
         " [-pretrain] [-pre_cap <N>] [-pre_wmin <W>] [-pre_tieskip] [-jitter <f>]"
         " [-ph] [-ph3_ratio <f>]\n")


def parse_args(argv):
    o = dict(sample=None, chromosome=None, snv_file=None, tries=100,
             min_maf=0.0, min_umis=4, penalty=5.0, scramble=False,
             improve_existing=False, no_softmasking=False, monostretch=15,
             jobs=1, samthreads=-1, seed=None, whitelist=None,
             wl_strip_gem=False, pretrain=False, pre_cap=5000,
             pre_wmin=3.0, pre_tieskip=True, jitter=0.25, phase3=False,
             ph3_ratio=0.90, bams=[])
    a = argv
    n = len(a)
    i = 0
    has_val = lambda: i < n - 1
    while i < n:
        f = a[i]
        if f == "-b":
            while i < n - 1 and not a[i + 1].startswith("-"):
                if not os.path.exists(a[i + 1]):
                    sys.exit("BAM file " + a[i + 1] + " does not exist")
                if a[i + 1] not in o["bams"]:
                    o["bams"].append(a[i + 1])
                i += 1
            i += 1
        elif f == "-o" and has_val(): o["sample"] = a[i + 1]; i += 2
        elif f == "-s" and has_val(): o["snv_file"] = a[i + 1]; i += 2
        elif f == "-r" and has_val(): o["chromosome"] = a[i + 1]; i += 2
        elif f == "-t" and has_val(): o["tries"] = int(a[i + 1]); i += 2
        elif f == "-u" and has_val(): o["min_umis"] = float(a[i + 1]); i += 2
        elif f == "-m" and has_val(): o["min_maf"] = float(a[i + 1]); i += 2
        elif f == "-p" and has_val(): o["penalty"] = float(a[i + 1]); i += 2
        elif f == "-x": o["scramble"] = True; i += 1
        elif f == "-i": o["improve_existing"] = True; i += 1
        elif f == "-sm": o["no_softmasking"] = True; i += 1
        elif f == "-ms" and has_val() and re.match(r"^(\d+)", a[i + 1]):
            o["monostretch"] = int(re.match(r"^(\d+)", a[i + 1]).group(1)); i += 2
        elif f == "-j" and has_val(): o["jobs"] = int(a[i + 1]); i += 2
        elif f == "-samthreads" and has_val(): o["samthreads"] = int(a[i + 1]); i += 2
        elif f == "-seed" and has_val(): o["seed"] = int(a[i + 1]); i += 2
        elif f == "-wl" and has_val():
            o["whitelist"] = a[i + 1]
            if not os.path.exists(o["whitelist"]):
                sys.exit("Barcode whitelist " + o["whitelist"] + " does not exist")
            i += 2
        elif f == "-wl_strip_gem": o["wl_strip_gem"] = True; i += 1
        elif f == "-pretrain": o["pretrain"] = True; i += 1
        elif f == "-pre_cap" and has_val(): o["pre_cap"] = int(a[i + 1]); i += 2
        elif f == "-pre_wmin" and has_val(): o["pre_wmin"] = float(a[i + 1]); i += 2
        elif f == "-pre_tieskip": o["pre_tieskip"] = True; i += 1
        elif f == "-jitter" and has_val(): o["jitter"] = float(a[i + 1]); i += 2
        elif f == "-ph": o["phase3"] = True; i += 1
        elif f == "-ph3_ratio" and has_val(): o["ph3_ratio"] = float(a[i + 1]); i += 2
        else:
            sys.exit("Unexpected/incomplete parameter: " + f + "\n")

    if not o["sample"]: sys.exit("No sample name\n" + USAGE)
    if not o["snv_file"]: sys.exit("No VCF file\n" + USAGE)
    if not o["bams"]: sys.exit("No BAM file\n" + USAGE)
    if o["wl_strip_gem"] and o["whitelist"] is None:
        sys.exit("-wl_strip_gem requires -wl <file>\n" + USAGE)

    if not o["chromosome"]: o["chromosome"] = "X"
    if o["tries"] <= 0: o["tries"] = 100
    if o["min_umis"] <= 0: o["min_umis"] = 10
    if o["min_maf"] < 0 or o["min_maf"] > 0.5: o["min_maf"] = 0.0
    if o["penalty"] < 1: o["penalty"] = 5.0
    if o["jobs"] < 1: o["jobs"] = 1
    if o["pre_cap"] < 100: o["pre_cap"] = 100
    if o["pre_wmin"] < 0: o["pre_wmin"] = 0.0
    o["jitter"] = min(max(o["jitter"], 0.0), 1.0)
    if o["samthreads"] < 0: o["samthreads"] = o["jobs"]
    if o["samthreads"] < 1: o["samthreads"] = 1
    return o


# =============================================================================
# INPUT
# =============================================================================
def open_maybe_gzip(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if path.lower().endswith(".gz") or magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


def load_whitelist(path, strip_gem):
    allowed = set()
    nonempty = 0
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r").lstrip()
            if line == "":
                continue
            bc = line.split(None, 1)[0]
            if not bc:
                continue
            nonempty += 1
            if strip_gem:
                bc = re.sub(r"-\d+$", "", bc)
            if not bc:
                continue
            allowed.add(bc)
    if not allowed:
        sys.exit("Zero barcodes loaded from whitelist %s\n" % path)
    return allowed, nonempty


def load_snvs(path, chromosome):
    snp = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            a = line.rstrip("\n").split("\t")
            if a[0] != chromosome:
                continue
            snp[int(a[1])] = [a[2], a[3], a[4]]
    return snp


def read_bams(o, snp_info, allowed, rng):
    """Parse WASP-tagged BAMs.

    Returns obs: dict (pos, cb_id, umi) -> [n_reads_allele1, n_reads_allele2,
    last_allele], plus the list of cell barcodes. This holds the same
    information as the two Perl hashes:
      %umis{pos}{allele}{"umi\\tcb"}  -> obs[...][allele-1] (> 0 means present)
      %cb2allele{cb}{pos}{umi}       -> obs[...][2] (last allele seen wins)
    """
    import pysam

    chrom = o["chromosome"]
    mono_re = re.compile("G{%d}|A{%d}|T{%d}|C{%d}" % ((o["monostretch"],) * 4))
    strip_gem = o["wl_strip_gem"]
    use_wl = allowed is not None
    obs = {}
    cb_index = {}
    cb_list = []
    st = dict(inf_reads=0, inf_alleles=0, hard=0, soft=0,
              wl_cand=0, wl_kept=0, wl_rej=0)
    wl_seen, wl_matched = set(), set()

    for path in sorted(o["bams"]):
        warn("    Reading %s ...\n" % path)
        bam = pysam.AlignmentFile(path, "rb", threads=o["samthreads"])
        if chrom not in bam.references:
            warn("    [warn] %s not in BAM header of %s; skipped\n" % (chrom, path))
            bam.close()
            continue
        for r in bam.fetch(chrom):
            if r.is_secondary:  # samtools view -F 256
                continue
            if not r.has_tag("vG") or not r.has_tag("vW") or r.get_tag("vW") != 1:
                continue
            seq = r.query_sequence
            if seq and mono_re.search(seq):
                continue
            cig = r.cigarstring or "*"
            if o["no_softmasking"] and "S" in cig:
                continue
            vG = r.get_tag("vG")
            vA = r.get_tag("vA") if r.has_tag("vA") else ()
            if r.has_tag("CB"):
                cb = r.get_tag("CB")
            elif r.has_tag("RG"):
                cb = r.get_tag("RG")
            else:
                sys.exit("No CB/RG tag in read:\n" + r.to_string())
            if cb == "-":
                continue
            if use_wl:
                st["wl_cand"] += 1
                m = re.sub(r"-\d+$", "", cb) if strip_gem else cb
                wl_seen.add(m)
                if m not in allowed:
                    st["wl_rej"] += 1
                    continue
                st["wl_kept"] += 1
                wl_matched.add(m)
            if r.has_tag("UB"):
                umi = r.get_tag("UB")
                if umi == "-":
                    continue
                st["hard"] += 1
            else:
                umi = "_".join((cb, r.reference_name, str(r.reference_start + 1),
                                str(r.flag), cig, str(r.template_length)))
                st["soft"] += 1
            st["inf_reads"] += 1
            cid = None
            for k in range(len(vG)):
                pos = vG[k] + 1
                if pos not in snp_info:
                    continue
                if k >= len(vA):
                    continue
                allele = vA[k]
                if allele != 1 and allele != 2:
                    continue
                if o["scramble"]:
                    allele = 1 + int(rng.integers(0, 2))
                st["inf_alleles"] += 1
                if cid is None:
                    cid = cb_index.get(cb)
                    if cid is None:
                        cid = len(cb_list)
                        cb_index[cb] = cid
                        cb_list.append(cb)
                key = (pos, cid, umi)
                rec = obs.get(key)
                if rec is None:
                    rec = [0, 0, 0]
                    obs[key] = rec
                rec[allele - 1] += 1
                rec[2] = allele
        bam.close()
        warn("    Reads with allelic info: %d  Alleles: %d  UMIs: %d/%d (hard/soft)\n"
             % (st["inf_reads"], st["inf_alleles"], st["hard"], st["soft"]))
    return obs, cb_list, st, wl_seen, wl_matched


def build_tables(obs, n_cells):
    """Collapse per-UMI observations to per-(SNV position, cell) counts.

    For every (pos, cell) pair we keep:
      n1, n2 : # distinct UMIs observed with allele 1 / 2   (%umis based;
               drives the optimiser score and VCF AD)
      c1, c2 : # UMIs whose last observed allele is 1 / 2   (%cb2allele based;
               drives bc2xci labels, pretraining and the profile)
      r1, r2 : read count of the first UMI seen with allele 1 / 2
               (used by Phase-3b, which in Perl reads one arbitrary UMI's count)
    """
    m = len(obs)
    pos = np.empty(m, np.int64)
    cid = np.empty(m, np.int64)
    a1 = np.empty(m, np.int64)
    a2 = np.empty(m, np.int64)
    last = np.empty(m, np.int8)
    for i, (k, v) in enumerate(obs.items()):
        pos[i] = k[0]; cid[i] = k[1]
        a1[i] = v[0]; a2[i] = v[1]; last[i] = v[2]
    # dict iteration is insertion order: the first record per (pos, cell)
    # in this order is the first UMI seen for that cell at that SNV.
    pair = pos * n_cells + cid
    order = np.argsort(pair, kind="stable")
    pair_s = pair[order]
    uniq, start = np.unique(pair_s, return_index=True)
    g = np.repeat(np.arange(uniq.size), np.diff(np.append(start, pair_s.size)))
    has1 = (a1[order] > 0).astype(np.int64)
    has2 = (a2[order] > 0).astype(np.int64)
    n1 = np.bincount(g, weights=has1, minlength=uniq.size).astype(np.int64)
    n2 = np.bincount(g, weights=has2, minlength=uniq.size).astype(np.int64)
    c1 = np.bincount(g, weights=(last[order] == 1), minlength=uniq.size).astype(np.int64)
    c2 = np.bincount(g, weights=(last[order] == 2), minlength=uniq.size).astype(np.int64)
    # first UMI (in insertion order) with allele 1 / 2 for each group
    r1 = np.zeros(uniq.size, np.int64)
    r2 = np.zeros(uniq.size, np.int64)
    a1s, a2s = a1[order], a2[order]
    idx1 = np.nonzero(a1s > 0)[0]
    g1, f1 = np.unique(g[idx1], return_index=True)
    r1[g1] = a1s[idx1[f1]]
    idx2 = np.nonzero(a2s > 0)[0]
    g2, f2 = np.unique(g[idx2], return_index=True)
    r2[g2] = a2s[idx2[f2]]
    return dict(pos=uniq // n_cells, cell=uniq % n_cells,
                n1=n1, n2=n2, c1=c1, c2=c2, r1=r1, r2=r2)


# =============================================================================
# NUMBA KERNELS
# =============================================================================
@njit(cache=False)
def _cell_term(h1, h2, pen):
    mx = h1 if h1 > h2 else h2
    mn = h2 if h1 > h2 else h1
    if mx == 0:
        return 0.0
    return (mx - 1) - pen * mn


@njit(cache=False)
def _move_delta(s, old, new, ptr, ecell, en1, en2, h1, h2, pen):
    d = 0.0
    for k in range(ptr[s], ptr[s + 1]):
        c = ecell[k]
        a = h1[c]; b = h2[c]
        na = a; nb = b
        if old == 1:
            na -= en1[k]; nb -= en2[k]
        elif old == -1:
            na -= en2[k]; nb -= en1[k]
        if new == 1:
            na += en1[k]; nb += en2[k]
        elif new == -1:
            na += en2[k]; nb += en1[k]
        d += _cell_term(na, nb, pen) - _cell_term(a, b, pen)
    return d


@njit(cache=False)
def _move_apply(s, old, new, ptr, ecell, en1, en2, h1, h2):
    for k in range(ptr[s], ptr[s + 1]):
        c = ecell[k]
        if old == 1:
            h1[c] -= en1[k]; h2[c] -= en2[k]
        elif old == -1:
            h1[c] -= en2[k]; h2[c] -= en1[k]
        if new == 1:
            h1[c] += en1[k]; h2[c] += en2[k]
        elif new == -1:
            h1[c] += en2[k]; h2[c] += en1[k]


@njit(cache=False)
def hill_climb(order, dirs, mono, ptr, ecell, en1, en2, n_cells, pen, seed,
               pass_imps, pass_score):
    """Greedy hill-climber of excise1.5.pl, with delta scoring.

    order[ele] = SNV index at slot ele, dirs[ele] = its direction (-1/0/+1).
    Modifies dirs in place. Returns (initial_score, final_score, n_passes).
    """
    np.random.seed(seed)
    h1 = np.zeros(n_cells, np.int64)
    h2 = np.zeros(n_cells, np.int64)
    n = order.size
    for ele in range(n):
        if dirs[ele] != 0:
            _move_apply(order[ele], 0, dirs[ele], ptr, ecell, en1, en2, h1, h2)
    score = 0.0
    for c in range(n_cells):
        score += _cell_term(h1[c], h2[c], pen)
    init_score = score

    npass = 0
    imps = 1
    last_imp = 0
    while imps > 0:
        imps = 0
        for ele in range(n):
            s = order[ele]
            old = dirs[ele]
            if old == 0:
                # dir=0 SNV: evaluate both +-1, keep the better if it improves
                dm = _move_delta(s, 0, -1, ptr, ecell, en1, en2, h1, h2, pen)
                dp = _move_delta(s, 0, 1, ptr, ecell, en1, en2, h1, h2, pen)
                if dp >= dm:
                    best_new = 1; bd = dp
                else:
                    best_new = -1; bd = dm
                if bd > 0:
                    _move_apply(s, 0, best_new, ptr, ecell, en1, en2, h1, h2)
                    dirs[ele] = best_new
                    score += bd
                    imps += 1; last_imp = ele
            else:
                for new in (-1, 0, 1):
                    if new == old:
                        continue
                    if new == 0 and mono[s]:
                        continue
                    d = _move_delta(s, old, new, ptr, ecell, en1, en2, h1, h2, pen)
                    # accept if strictly better, or tie and new == 0
                    if d > 0 or (d == 0 and new == 0):
                        _move_apply(s, old, new, ptr, ecell, en1, en2, h1, h2)
                        imps += 1; last_imp = ele
                        dirs[ele] = new
                        old = new
                        score += d
            # Early exit with one random probe (same as excise1.5.pl)
            if imps == 0 and last_imp - ele == 1:
                rnd = np.random.randint(0, n)
                if rnd != ele and rnd != ele - 1:
                    rs = order[rnd]
                    rold = dirs[rnd]
                    for rnew in (-1, 0, 1):
                        if rnew == rold:
                            continue
                        if rnew == 0 and mono[rs]:
                            continue
                        d = _move_delta(rs, rold, rnew, ptr, ecell, en1, en2, h1, h2, pen)
                        if d > 0 or (d == 0 and rnew == 0):
                            _move_apply(rs, rold, rnew, ptr, ecell, en1, en2, h1, h2)
                            imps += 1; last_imp = rnd
                            dirs[rnd] = rnew
                            score += d
                            break
                if imps == 0:
                    break
        if npass < pass_imps.size:
            pass_imps[npass] = imps
            pass_score[npass] = score
        npass += 1
    return init_score, score, npass


@njit(cache=False)
def _pair_votes(cptr, csnv, cal, n_snv, cap, cell_lo, cell_hi, out_off, keys, sgn, seed):
    """Co-segregation votes for cells [cell_lo, cell_hi) (excise1.5 pretrain)."""
    np.random.seed(seed)
    for c in range(cell_lo, cell_hi):
        lo = cptr[c]
        L = cptr[c + 1] - lo
        if L < 2:
            continue
        o = out_off[c - cell_lo]
        allp = L * (L - 1) // 2
        if allp <= cap:
            for i in range(L):
                for j in range(i + 1, L):
                    si = csnv[lo + i]; sj = csnv[lo + j]
                    if si > sj:
                        si, sj = sj, si
                    keys[o] = si * n_snv + sj
                    sgn[o] = 1 if cal[lo + i] == cal[lo + j] else -1
                    o += 1
        else:
            seen = set()
            seen.add(np.int64(-1))
            tries = 0
            got = 0
            while got < cap and tries < 5 * cap:
                tries += 1
                i = np.random.randint(0, L)
                j = np.random.randint(0, L)
                if i == j:
                    continue
                if i > j:
                    i, j = j, i
                kk = np.int64(i) * L + j
                if kk in seen:
                    continue
                seen.add(kk)
                si = csnv[lo + i]; sj = csnv[lo + j]
                if si > sj:
                    si, sj = sj, si
                keys[o] = si * n_snv + sj
                sgn[o] = 1 if cal[lo + i] == cal[lo + j] else -1
                o += 1
                got += 1
            # unused slots keep key -1 and are dropped by the caller


@njit(cache=False)
def _bfs_seed(cand, aptr, anbr, aw, n_snv, wmin):
    """BFS seeding with bottleneck weights (_bfs_seed_with_uncertainty)."""
    dirs = np.zeros(n_snv, np.int8)
    seen = np.zeros(n_snv, np.bool_)
    bott = np.zeros(n_snv, np.float64)
    q = np.empty(n_snv, np.int64)
    for start in cand:
        if seen[start]:
            continue
        if aptr[start + 1] == aptr[start]:
            continue
        seen[start] = True
        dirs[start] = 1
        bott[start] = 1e9
        qh = 0; qt = 0
        q[qt] = start; qt += 1
        while qh < qt:
            u = q[qh]; qh += 1
            for k in range(aptr[u], aptr[u + 1]):
                v = anbr[k]
                if seen[v]:
                    continue
                w = aw[k]
                if abs(w) < wmin:
                    continue
                seen[v] = True
                dirs[v] = dirs[u] if w >= 0 else -dirs[u]
                bv = abs(w)
                bott[v] = bv if bv < bott[u] else bott[u]
                q[qt] = v; qt += 1
    return dirs, seen, bott


@njit(cache=False)
def _phase3a(cand, ptr, ecell, en1, en2, h1, h2, ratio):
    """Phase-3 directionality assignment against fixed cb2phase_best."""
    out = np.zeros(cand.size, np.int8)
    for t in range(cand.size):
        s = cand[t]
        if ptr[s + 1] == ptr[s]:
            continue
        nd = np.zeros(2, np.int64)
        sh1 = np.zeros(2, np.float64)
        sh2 = np.zeros(2, np.float64)
        for ti in range(2):
            trial = 1 if ti == 0 else -1
            for k in range(ptr[s], ptr[s + 1]):
                c = ecell[k]
                if trial == 1:
                    a = h1[c] + en1[k]; b = h2[c] + en2[k]
                else:
                    a = h1[c] + en2[k]; b = h2[c] + en1[k]
                tot = a + b
                if tot <= 0:
                    continue
                sh1[ti] += a; sh2[ti] += b
                if a / tot >= 0.9 or b / tot >= 0.9:
                    nd[ti] += 1
        if nd[0] > nd[1]:
            wi = 0
        elif nd[1] > nd[0]:
            wi = 1
        else:
            dpos = sh1[0] if sh1[0] > sh2[0] else sh2[0]
            dneg = sh1[1] if sh1[1] > sh2[1] else sh2[1]
            wi = 0 if dpos >= dneg else 1
        tot = sh1[wi] + sh2[wi]
        frac = 0.0
        if tot > 0:
            frac = (sh1[wi] if sh1[wi] > sh2[wi] else sh2[wi]) / tot
        if frac >= ratio:
            out[t] = 1 if wi == 0 else -1
    return out


# =============================================================================
# HELPERS ON THE CSR TABLES
# =============================================================================
def hap_totals(dirs, snv_of_entry, cell, x1, x2, n_cells):
    """Per-cell haplotype UMI totals given SNV directions (vectorised)."""
    d = dirs[snv_of_entry]
    w1 = np.where(d == 1, x1, np.where(d == -1, x2, 0))
    w2 = np.where(d == 1, x2, np.where(d == -1, x1, 0))
    h1 = np.bincount(cell, weights=w1, minlength=n_cells).astype(np.int64)
    h2 = np.bincount(cell, weights=w2, minlength=n_cells).astype(np.int64)
    return h1, h2


def cell_hap_from(h1, h2):
    tot = h1 + h2
    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = np.where(tot > 0, h1 / np.maximum(tot, 1), 0.0)
        f2 = np.where(tot > 0, h2 / np.maximum(tot, 1), 0.0)
    return np.where(tot == 0, 0, np.where(f1 >= 0.9, 1, np.where(f2 >= 0.9, 2, 0)))


def classify(h1, h2):
    tot = h1 + h2
    safe = np.maximum(tot, 1)
    lab = np.full(h1.size, 4, np.int8)  # Low_coverage
    both = (h1 > 0) & (h2 > 0)
    lab[both] = 3
    lab[(h2 >= 2) & (h2 / safe >= 0.9)] = 2
    lab[(h1 >= 2) & (h1 / safe >= 0.9)] = 1
    lab[tot == 0] = 0
    return lab


LABELS = np.array(["Unknown", "X1", "X2", "Both", "Low_coverage"], dtype=object)


def score_parts(h1, h2, pen):
    mx = np.maximum(h1, h2)
    mn = np.minimum(h1, h2)
    has = mx > 0
    total = int((h1 + h2).sum())
    disc = int(mn[has].sum())
    conc = int((mx[has] - 1).sum())
    return total, disc, conc


# =============================================================================
# PARALLEL TRY WORKER (globals are inherited by fork)
# =============================================================================
G = {}


def run_try(t):
    o = G["o"]
    n = G["n_snv"]
    seed = (o["seed"] if o["seed"] is not None else G["time_seed"]) + t * 1337
    rng = np.random.default_rng(seed)
    warn("Try #%d, initialising ...\n" % t)
    local_pos = np.arange(n, dtype=np.int64)
    if o["pretrain"]:
        local_dir = G["pre_seed"].copy()
        u = G["unc"] * o["jitter"]
        flip = rng.random(n) < u
        rnd_sign = np.where(rng.random(n) < 0.5, 1, -1).astype(np.int8)
        local_dir = np.where(flip, np.where(local_dir == 0, rnd_sign, -local_dir),
                             local_dir).astype(np.int8)
        for i in range(n - 1, 0, -1):  # Fisher-Yates, dirs move with positions
            j = int(rng.integers(0, i + 1))
            if i != j:
                local_pos[i], local_pos[j] = local_pos[j], local_pos[i]
                local_dir[i], local_dir[j] = local_dir[j], local_dir[i]
    else:
        local_dir = np.zeros(n, np.int8)
        mono = G["mono"]
        for i in range(n):
            a = int(rng.integers(0, n)); b = int(rng.integers(0, n))
            if a != b:
                local_pos[a], local_pos[b] = local_pos[b], local_pos[a]
            if mono[local_pos[i]]:
                local_dir[i] = 1 if rng.random() < 0.5 else -1
            else:
                local_dir[i] = int(rng.integers(0, 3)) - 1

    pi = np.zeros(10000, np.int64)
    ps = np.zeros(10000, np.float64)
    T = G["opt"]
    init, best, npass = hill_climb(local_pos, local_dir, G["mono"], T["ptr"],
                                   T["cell"], T["n1"], T["n2"], G["n_cells"],
                                   float(o["penalty"]), seed % (2 ** 31), pi, ps)
    warn("    SNVs: %d  Initial score: %s\n" % (n, perl_num(init)))
    for k in range(min(npass, pi.size)):
        warn("    Pass: %d  Improvements: %d  Score: %s\n" % (k + 1, pi[k], perl_num(ps[k])))
    warn("    Finished Try #%d  Final score: %s\n" % (t, perl_num(best)))
    dirs = np.zeros(n, np.int8)
    dirs[local_pos] = local_dir  # realign to canonical SNV order
    return t, best, dirs


# =============================================================================
# MAIN
# =============================================================================
def main(argv):
    t_start = time.time()
    warn("Reading parameters ...\n")
    o = parse_args(argv)
    rng = np.random.default_rng(o["seed"])
    chrom = o["chromosome"]
    warn("Sample             : %s\n" % o["sample"])
    warn("VCF file           : %s\n" % o["snv_file"])
    warn("BAM files          : %d (%s)\n" % (len(o["bams"]), ",".join(o["bams"])))
    warn("Chromosome         : %s\n" % chrom)
    warn("Tries              : %d\n" % o["tries"])
    warn("Min UMIs/SNV       : %s\n" % perl_num(o["min_umis"]))
    warn("Min MAF            : %s\n" % perl_num(o["min_maf"]))
    warn("Discordant penalty : %s\n" % perl_num(o["penalty"]))
    warn("Parallel jobs      : %d\n" % o["jobs"])
    warn("samtools threads   : %d\n" % o["samthreads"])
    warn("Barcode whitelist  : %s\n" % (o["whitelist"] or "disabled"))
    if o["whitelist"]:
        warn("Strip GEM suffix   : %s\n" % ("yes" if o["wl_strip_gem"] else "no"))
    warn("Pretrain           : %s\n" % ("yes" if o["pretrain"] else "no"))
    warn("Jitter scale       : %s\n" % perl_num(o["jitter"]))
    warn("Phase 3            : %s\n" % ("enabled" if o["phase3"] else "disabled"))

    allowed = None
    wl_nonempty = 0
    if o["whitelist"]:
        warn("Reading barcode whitelist ...\n")
        allowed, wl_nonempty = load_whitelist(o["whitelist"], o["wl_strip_gem"])
        warn("Whitelist barcodes loaded: %d unique from %d non-empty records\n"
             % (len(allowed), wl_nonempty))

    warn("Reading SNVs ...\n")
    snp_info = load_snvs(o["snv_file"], chrom)
    warn("%d SNVs loaded from %s\n" % (len(snp_info), o["snv_file"]))
    if not snp_info:
        sys.exit("Zero SNVs, run terminated\n")

    warn("Reading BAM files ...\n")
    obs, cb_list, st, wl_seen, wl_matched = read_bams(o, snp_info, allowed, rng)
    n_cells = len(cb_list)
    if allowed is not None:
        warn("Whitelist candidate barcodes observed in BAM: %d\n" % len(wl_seen))
        warn("Whitelist barcodes matched in BAM: %d / %d loaded\n" % (len(wl_matched), len(allowed)))
        warn("Whitelist CB-tagged candidate reads: checked=%d kept=%d rejected=%d\n"
             % (st["wl_cand"], st["wl_kept"], st["wl_rej"]))
        warn("Whitelist barcodes with usable chr%s allelic evidence: %d / %d loaded\n"
             % (chrom, n_cells, len(allowed)))
    tb = build_tables(obs, max(n_cells, 1))
    del obs
    all_pos = np.unique(tb["pos"])
    warn("%d SNV positions covered in %d cells\n" % (all_pos.size, n_cells))
    if all_pos.size == 0:
        sys.exit("Zero SNVs after BAM reading\n")
    if n_cells == 0:
        sys.exit("Zero cells after BAM reading\n")

    # ---------------- FILTER SNVs ----------------
    warn("Checking AFs after WASP genotyping ...\n")
    sidx = np.searchsorted(all_pos, tb["pos"])
    N1 = np.bincount(sidx, weights=tb["n1"], minlength=all_pos.size).astype(np.int64)
    N2 = np.bincount(sidx, weights=tb["n2"], minlength=all_pos.size).astype(np.int64)
    tot = N1 + N2
    low_umis_m = tot < o["min_umis"]
    mn = np.minimum(N1, N2)
    ratios = np.sort(np.where(tot > 0, 100.0 * mn / np.maximum(tot, 1), 0.0)[~low_umis_m])
    af = np.where(tot > 0, N1 / np.maximum(tot, 1), 0.0)
    low_maf_m = (~low_umis_m) & ((af < o["min_maf"]) | (af > 1 - o["min_maf"]))
    keep = ~low_umis_m & ~low_maf_m
    phased_pos = all_pos[keep]
    if phased_pos.size == 0:
        sys.exit("No SNVs to phase\n")
    pN1, pN2 = N1[keep], N2[keep]
    mono = (pN1 == 0) | (pN2 == 0)
    n_snv = phased_pos.size
    warn("Excluded %d SNVs < %s UMIs; %d SNVs with MAF < %s\n"
         % (int(low_umis_m.sum()), perl_num(o["min_umis"]), int(low_maf_m.sum()), perl_num(o["min_maf"])))
    warn("Monoallelic SNVs in pool: %d / %d\n" % (int(mono.sum()), n_snv))
    if ratios.size % 2:
        median_maf = ratios[ratios.size // 2]
    else:
        median_maf = (ratios[ratios.size // 2 - 1] + ratios[ratios.size // 2]) / 2
    warn("Median MAF: %s\n" % perl_num(median_maf))

    # entries restricted to phased SNVs, sorted by SNV then cell (CSR by SNV)
    snv_map = np.full(all_pos.size, -1, np.int64)
    snv_map[np.nonzero(keep)[0]] = np.arange(n_snv)
    e_snv = snv_map[sidx]
    m = e_snv >= 0
    order = np.lexsort((tb["cell"][m], e_snv[m]))
    E = {k: tb[k][m][order] for k in ("cell", "n1", "n2", "c1", "c2", "r1", "r2")}
    E["snv"] = e_snv[m][order]
    ptr = np.zeros(n_snv + 1, np.int64)
    np.cumsum(np.bincount(E["snv"], minlength=n_snv), out=ptr[1:])
    opt = dict(ptr=ptr, cell=E["cell"], n1=E["n1"], n2=E["n2"])
    del tb

    # ---------------- PRETRAINING ----------------
    pre_seed = np.zeros(n_snv, np.int8)
    unc = np.ones(n_snv, np.float64)
    if o["pretrain"]:
        warn("Running improved pretraining ...\n")
        warn("  [pretrain] Building cell majority calls ...\n")
        c1, c2 = E["c1"], E["c2"]
        ok = (c1 + c2) > 0
        if o["pre_tieskip"]:
            ok &= c1 != c2
        call = np.where(c1 >= c2, 1, 2).astype(np.int8)
        corder = np.lexsort((E["snv"][ok], E["cell"][ok]))
        ccell = E["cell"][ok][corder]
        csnv = E["snv"][ok][corder]
        cal = call[ok][corder]
        cptr = np.zeros(n_cells + 1, np.int64)
        np.cumsum(np.bincount(ccell, minlength=n_cells), out=cptr[1:])

        warn("  [pretrain] Building co-segregation graph ...\n")
        L = np.diff(cptr)
        npairs = np.where(L >= 2, np.minimum(L * (L - 1) // 2, o["pre_cap"]), 0)
        acc_k, acc_w, acc_d = [], [], np.zeros(n_snv, np.float64)
        batch_lim = 20_000_000
        c0 = 0
        while c0 < n_cells:
            cs = np.cumsum(npairs[c0:])
            c1_ = c0 + max(1, int(np.searchsorted(cs, batch_lim, side="right")))
            c1_ = min(c1_, n_cells)
            nb = npairs[c0:c1_]
            off = np.zeros(nb.size, np.int64)
            off[1:] = np.cumsum(nb)[:-1]
            tot_b = int(nb.sum())
            if tot_b:
                keys = np.full(tot_b, -1, np.int64)
                sg = np.zeros(tot_b, np.int64)
                _pair_votes(cptr, csnv, cal, n_snv, o["pre_cap"], c0, c1_, off,
                            keys, sg, int(rng.integers(0, 2 ** 31)))
                v = keys >= 0
                keys, sg = keys[v], sg[v]
                acc_d += np.bincount(keys // n_snv, minlength=n_snv)
                acc_d += np.bincount(keys % n_snv, minlength=n_snv)
                uk, inv = np.unique(keys, return_inverse=True)
                acc_k.append(uk)
                acc_w.append(np.bincount(inv, weights=sg))
                if len(acc_k) > 1:
                    kk = np.concatenate(acc_k); ww = np.concatenate(acc_w)
                    uk, inv = np.unique(kk, return_inverse=True)
                    acc_k, acc_w = [uk], [np.bincount(inv, weights=ww)]
            c0 = c1_
        if acc_k:
            ukeys, uw = acc_k[0], acc_w[0]
        else:
            ukeys, uw = np.zeros(0, np.int64), np.zeros(0)
        ei, ej = ukeys // n_snv, ukeys % n_snv
        src = np.concatenate([ei, ej]); dst = np.concatenate([ej, ei]); ww = np.concatenate([uw, uw])
        ao = np.lexsort((dst, src))
        anbr, aw = dst[ao], ww[ao]
        aptr = np.zeros(n_snv + 1, np.int64)
        np.cumsum(np.bincount(src, minlength=n_snv), out=aptr[1:])

        warn("  [pretrain] BFS seeding (wmin=%s) ...\n" % perl_num(o["pre_wmin"]))
        cand = np.argsort(-acc_d, kind="stable").astype(np.int64)
        pre_seed, assigned_m, bott = _bfs_seed(cand, aptr, anbr, aw, n_snv, float(o["pre_wmin"]))
        bott[~assigned_m] = 0.0
        unc = 1.0 / (1.0 + bott)
        edges = int(anbr.size)
        strong = int((np.abs(aw) >= o["pre_wmin"]).sum())
        assigned = int((pre_seed != 0).sum())
        warn("  [pretrain] nodes=%d  edges~%d  strong-edges~%d  assigned=%d (%.1f%%)\n"
             % (n_snv, edges, strong, assigned, 100.0 * assigned / n_snv if n_snv else 0))

    # ---------------- PARALLEL TRIES ----------------
    warn("Starting XCI calling\n")
    G.update(o=o, n_snv=n_snv, n_cells=n_cells, mono=mono, pre_seed=pre_seed,
             unc=unc, opt=opt, time_seed=int(time.time()))
    # compile kernels once in the parent so forked workers inherit them
    hill_climb(np.zeros(1, np.int64), np.zeros(1, np.int8), np.zeros(1, np.bool_),
               np.zeros(2, np.int64), np.zeros(0, np.int64), np.zeros(0, np.int64),
               np.zeros(0, np.int64), 1, 1.0, 0, np.zeros(1, np.int64), np.zeros(1))
    t_opt = time.time()
    results = {}
    if o["jobs"] > 1 and o["tries"] > 1:
        ctx = mp.get_context("fork")
        with ctx.Pool(o["jobs"]) as pool:
            for t, sc, d in pool.imap_unordered(run_try, range(1, o["tries"] + 1)):
                results[t] = (sc, d)
    else:
        for t in range(1, o["tries"] + 1):
            _, sc, d = run_try(t)
            results[t] = (sc, d)
    warn("    Optimiser wall time: %.2f s\n" % (time.time() - t_opt))

    global_best = None
    best_dir = None
    for t in range(1, o["tries"] + 1):
        sc, d = results[t]
        if global_best is None or sc > global_best:
            global_best, best_dir = sc, d.copy()

    summ_path = "%s_chr%s_XCISE_summary.txt" % (o["sample"], chrom)
    if o["improve_existing"] and os.path.exists(summ_path) and os.path.getsize(summ_path) > 0:
        with open(summ_path) as fh:
            line = fh.readline()
        mt = re.match(r"^Best\sscore\s+:\s+(\d+)", line or "")
        if mt and (not global_best or int(mt.group(1)) > global_best):
            global_best = int(mt.group(1))
            warn("    Adopted best score from previous run: %d\n" % global_best)

    dirs = best_dir
    pen = float(o["penalty"])
    esnv, ecell = E["snv"], E["cell"]

    # ---------------- PHASE 2: monoallelic SNVs by cell vote ----------------
    warn("    Phase-2 cell-vote correction for monoallelic SNVs ...\n")
    bi_dirs = np.where(mono, 0, dirs).astype(np.int8)
    ch = cell_hap_from(*hap_totals(bi_dirs, esnv, ecell, E["c1"], E["c2"], n_cells))
    obs_is2 = pN2 > 0
    e_obs = np.where(obs_is2[esnv], E["n2"], E["n1"]) > 0
    eh = ch[ecell]
    v1 = np.bincount(esnv, weights=(e_obs & (eh == 1)), minlength=n_snv)
    v2 = np.bincount(esnv, weights=(e_obs & (eh == 2)), minlength=n_snv)
    tgt = mono & ((v1 + v2) > 0)
    obs_allele1 = ~obs_is2
    correct = np.where(v1 >= v2, np.where(obs_allele1, 1, -1), np.where(obs_allele1, -1, 1))
    changed = tgt & (correct != dirs)
    warn("    Mono SNVs corrected: %d  already correct: %d\n"
         % (int(changed.sum()), int((tgt & ~changed).sum())))
    dirs = np.where(tgt, correct, dirs).astype(np.int8)

    # ---------------- PHASE 3 (opt-in) ----------------
    if o["phase3"]:
        warn("    Phase-3 directionality assignment for dir=0 biallelic SNVs ...\n")
        h1, h2 = hap_totals(dirs, esnv, ecell, E["n1"], E["n2"], n_cells)
        cand = np.nonzero((dirs == 0) & ~mono)[0].astype(np.int64)
        res = _phase3a(cand, ptr, ecell, E["n1"], E["n2"], h1, h2, float(o["ph3_ratio"]))
        has_cov = np.diff(ptr)[cand] > 0
        warn("    Phase-3: assigned=%d  stayed Unk=%d\n"
             % (int((res != 0).sum()), int(((res == 0) & has_cov).sum())))
        dirs[cand] = np.where(res != 0, res, dirs[cand])

        warn("    Phase-3 cell-vote assignment for dir=0 biallelic SNVs ...\n")
        ch = cell_hap_from(*hap_totals(dirs, esnv, ecell, E["c1"], E["c2"], n_cells))
        eh = ch[ecell]
        r1 = np.where(E["n1"] > 0, E["r1"], 0)
        r2 = np.where(E["n2"] > 0, E["r2"], 0)
        bc = lambda w: np.bincount(esnv, weights=w, minlength=n_snv)
        conc_pos = bc(np.where(eh == 1, r1, 0)) + bc(np.where(eh == 2, r2, 0))
        conc_neg = bc(np.where(eh == 1, r2, 0)) + bc(np.where(eh == 2, r1, 0))
        tot_ev = conc_pos + conc_neg
        c3 = (dirs == 0) & ~mono & (tot_ev > 0)
        win = np.where(conc_pos >= conc_neg, 1, -1)
        wconc = np.where(conc_pos >= conc_neg, conc_pos, conc_neg)
        ok = c3 & (wconc / np.maximum(tot_ev, 1) >= o["ph3_ratio"])
        warn("    Phase-3: assigned=%d  stayed Unk=%d\n" % (int(ok.sum()), int((c3 & ~ok).sum())))
        dirs = np.where(ok, win, dirs).astype(np.int8)
    else:
        warn("    Phase 3 disabled (use -ph to enable); dir=0 SNVs remain unresolved.\n")

    # ---------------- OUTPUTS (flip so that X2 >= X1) ----------------
    while True:
        h1, h2 = hap_totals(dirs, esnv, ecell, E["c1"], E["c2"], n_cells)
        lab = classify(h1, h2)
        counts = np.bincount(lab, minlength=5)
        t0, t1, t2, t3, t4 = (int(x) for x in counts)
        warn("    X1/X2/Both/LowC/Unk: %d / %d / %d / %d / %d\n" % (t1, t2, t3, t4, t0))
        if t1 > t2:
            warn("    Flipping X1<->X2 so that X2 >= X1 ...\n")
            dirs = (-dirs).astype(np.int8)
            continue
        break

    prefix = "%s_chr%s_XCISE" % (o["sample"], chrom)
    # Python str ordering == Perl's default sort for byte/ASCII barcodes
    cb_rank_order = np.array(sorted(range(n_cells), key=lambda i: cb_list[i]), dtype=np.int64)

    warn("    Outputting barcode-to-XCI ...\n")
    with open(prefix + "_bc2xci.txt", "w") as fh:
        fh.write("".join("%s\t%d\t%d\t%s\n" % (cb_list[c], h1[c], h2[c], LABELS[lab[c]])
                         for c in cb_rank_order))

    hb1, hb2 = hap_totals(dirs, esnv, ecell, E["n1"], E["n2"], n_cells)
    total_umi, disc, conc = score_parts(hb1, hb2, pen)
    ref1 = int((dirs == 1).sum()); alt1 = int((dirs == -1).sum()); unk1 = int((dirs == 0).sum())
    n_mono = int(mono.sum())
    grand = t0 + t1 + t2 + t3 + t4
    warn("    Total XCI-informative UMIs   : %d\n" % total_umi)
    warn("    SNVs informative/non-informative/mono : %d/%d/%d\n" % (ref1 + alt1 - n_mono, unk1, n_mono))
    pct = lambda x: 100.0 * x / grand if grand else 0.0
    with open(prefix + "_summary.txt", "w") as fs:
        fs.write("Best score  : %s\n" % perl_num(global_best))
        fs.write("Total UMIs  : %d\n" % total_umi)
        fs.write("Concordant  : %d\n" % conc)
        fs.write("Discordant  : %d\n" % disc)
        fs.write("Discordancy : %s\n" % (perl_num(disc / (conc + disc)) if conc + disc else "N/A"))
        fs.write("Total SNVs        : %d\n" % n_snv)
        fs.write("Informative SNVs  : %d\n" % (ref1 + alt1 - n_mono))
        fs.write("Non-informative   : %d\n" % unk1)
        fs.write("Monoallelic SNVs  : %d\n" % n_mono)
        fs.write("Median MAF        : %3.2f %%\n" % median_maf)
        fs.write("X1 cells    : %5d( %3.2f %% )\n" % (t1, pct(t1)))
        fs.write("X2 cells    : %5d( %3.2f %% )\n" % (t2, pct(t2)))
        fs.write("Both X      : %5d( %3.2f %% )\n" % (t3, pct(t3)))
        fs.write("LowCoverage : %5d( %3.2f %% )\n" % (t4, pct(t4)))
        fs.write("Unknown     : %5d( %3.2f %% )\n" % (t0, pct(t0)))
        if o["whitelist"]:
            fs.write("Whitelist file     : %s\n" % o["whitelist"])
            fs.write("Whitelist strip GEM: %s\n" % ("yes" if o["wl_strip_gem"] else "no"))
            fs.write("Whitelist loaded   : %d\n" % len(allowed))
            fs.write("Whitelist BAM match: %d\n" % len(wl_matched))
            fs.write("Whitelist evidence  : %d\n" % n_cells)
            fs.write("Whitelist reads     : checked=%d kept=%d rejected=%d\n"
                     % (st["wl_cand"], st["wl_kept"], st["wl_rej"]))
        else:
            fs.write("Barcode whitelist  : disabled\n")

    warn("    Outputting VCF ...\n")
    with open(prefix + ".vcf", "w") as fv:
        fv.write("##fileformat=VCFv4.2\n")
        fv.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]) + "\n")
        fv.write('##INFO=<ID=X1A,Number=1,Type=String,Description="X1 haplotype allele: Ref, Alt, or Unk (escape candidate)">\n')
        fv.write('##INFO=<ID=MONO,Number=0,Type=Flag,Description="Monoallelic SNV: direction by cell-vote not optimizer">\n')
        fv.write('##INFO=<ID=AD,Number=2,Type=Integer,Description="Allele depth: ref_umi,alt_umi">\n')
        fv.write('##INFO=<ID=X1,Number=1,Type=Integer,Description="UMI count supporting X1 haplotype">\n')
        fv.write('##INFO=<ID=X2,Number=1,Type=Integer,Description="UMI count supporting X2 haplotype">\n')
        fv.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total UMI depth (X1 + X2)">\n')
        out = []
        for i in range(n_snv):
            p = int(phased_pos[i]); d = int(dirs[i]); n1 = int(pN1[i]); n2 = int(pN2[i])
            x1a = "Unk" if d == 0 else ("Ref" if d == 1 else "Alt")
            x1u = n1 if d == 1 else (n2 if d == -1 else 0)
            x2u = n1 if d == -1 else (n2 if d == 1 else 0)
            dp = x1u + x2u
            rs, ref, alt = snp_info[p]
            out.append("%s\t%d\t%s\t%s\t%s\t%d\tPASS\tX1A=%s;AD=%d,%d;X1=%d;X2=%d;DP=%d%s\n"
                       % (chrom, p, rs, ref, alt, dp, x1a, n1, n2, x1u, x2u, dp,
                          ";MONO" if mono[i] else ""))
        fv.write("".join(out))

    warn("    Outputting haplotype sequence table ...\n")
    with open(prefix + "_haplotypes.tsv", "w") as fh:
        fh.write("CHROM\tPOS\tID\tREF\tALT\tX1_allele\tX2_allele\tmono\n")
        out = []
        for i in range(n_snv):
            d = int(dirs[i])
            if d == 0:
                continue
            p = int(phased_pos[i]); rs, ref, alt = snp_info[p]
            x1, x2 = (ref, alt) if d == 1 else (alt, ref)
            out.append("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\n" % (chrom, p, rs, ref, alt, x1, x2, int(mono[i])))
        fh.write("".join(out))

    warn("    Outputting per-cell SNV profile ...\n")
    rank = np.empty(n_cells, np.int64)
    rank[cb_rank_order] = np.arange(n_cells)
    po = np.lexsort((esnv, rank[ecell]))
    with open(prefix + "_cell_snv_profile.tsv", "w") as fp:
        fp.write("\t".join(["cell_barcode", "pos", "rs", "ref", "alt", "snv_dir", "x1a",
                            "ref_umis", "alt_umis", "total_umis", "hap_x1_umis",
                            "hap_x2_umis", "cell_xci", "mono"]) + "\n")
        x1a_s = {1: "Ref", -1: "Alt", 0: "Unk"}
        buf = []
        for k in po:
            c = int(ecell[k]); s = int(esnv[k])
            r_, a_ = int(E["c1"][k]), int(E["c2"][k])
            if r_ + a_ <= 0:
                continue
            d = int(dirs[s]); p = int(phased_pos[s]); rs, ref, alt = snp_info[p]
            hx1, hx2 = (r_, a_) if d == 1 else ((a_, r_) if d == -1 else (0, 0))
            buf.append("%s\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\n"
                       % (cb_list[c], p, rs, ref, alt, d, x1a_s[d], r_, a_, r_ + a_,
                          hx1, hx2, LABELS[lab[c]], int(mono[s])))
            if len(buf) >= 100000:
                fp.write("".join(buf)); buf = []
        fp.write("".join(buf))
    warn("    Per-cell SNV profile written: %s_cell_snv_profile.tsv\n" % prefix)
    warn("    Total wall time: %.2f s\n" % (time.time() - t_start))


if __name__ == "__main__":
    main(sys.argv[1:])
