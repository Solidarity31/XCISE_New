// main.cpp -- command-line interface of the C++ port of excise1.5.pl.
//
// Derived from xcise.pl of XCISE (https://github.com/Vityay/XCISE).
// Please cite: Henning RH, Rust TM, Dijksterhuis K, Eggen BJL, Guryev V.
//   bioRxiv (2024). doi:10.1101/2024.08.29.610317
//
// Build:  make -C cpp            (needs htslib + zlib, C++17)
// Usage:  cpp/excise15 -o <prefix> -s <vcf> -b <bam> [...]   (flags as excise1.5.pl)

#include "excise_core.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

using std::string;
using std::vector;

[[noreturn]] static void die(const string& m) {
    fputs(m.c_str(), stderr);
    if (m.empty() || m.back() != '\n') fputc('\n', stderr);
    exit(1);
}

static bool is_number_prefix(const char* s, int* out) {
    if (!isdigit((unsigned char)s[0])) return false;
    *out = atoi(s);
    return true;
}

// ---------------------------------------------------------------------------
// options (same flags, defaults and quirks as excise1.5.pl)
// ---------------------------------------------------------------------------


static const char* USAGE =
    "Usage: excise15 -o <prefix> -s <vcf> -b <bam> [...] [-r <chrom>] [-t <tries>]"
    " [-u <min_umis>] [-m <min_maf>] [-p <penalty>]"
    " [-j <jobs>] [-samthreads <N>] [-seed <int>]"
    " [-wl <barcodes.txt|tsv[.gz]>] [-wl_strip_gem]"
    " [-pretrain] [-pre_cap <N>] [-pre_wmin <W>] [-pre_tieskip] [-jitter <f>]"
    " [-ph] [-ph3_ratio <f>]\n";

static excise::Opt parse_args(int argc, char** argv) {
    excise::Opt o;
    vector<string> a(argv + 1, argv + argc);
    size_t n = a.size(), i = 0;
    auto hv = [&]() { return i + 1 < n; };
    auto fexists = [](const string& p) { FILE* f = fopen(p.c_str(), "r"); if (f) fclose(f); return f != nullptr; };
    while (i < n) {
        const string& f = a[i];
        int ms;
        if (f == "-b") {
            while (i + 1 < n && a[i + 1][0] != '-') {
                if (!fexists(a[i + 1])) die("BAM file " + a[i + 1] + " does not exist");
                if (std::find(o.bams.begin(), o.bams.end(), a[i + 1]) == o.bams.end()) o.bams.push_back(a[i + 1]);
                i++;
            }
            i++;
        }
        else if (f == "-o" && hv()) { o.sample = a[i + 1]; i += 2; }
        else if (f == "-s" && hv()) { o.snv_file = a[i + 1]; i += 2; }
        else if (f == "-r" && hv()) { o.chromosome = a[i + 1]; i += 2; }
        else if (f == "-t" && hv()) { o.tries = atoi(a[i + 1].c_str()); i += 2; }
        else if (f == "-u" && hv()) { o.min_umis = atof(a[i + 1].c_str()); i += 2; }
        else if (f == "-m" && hv()) { o.min_maf = atof(a[i + 1].c_str()); i += 2; }
        else if (f == "-p" && hv()) { o.penalty = atof(a[i + 1].c_str()); i += 2; }
        else if (f == "-x") { o.scramble = true; i++; }
        else if (f == "-i") { o.improve_existing = true; i++; }
        else if (f == "-sm") { o.no_softmasking = true; i++; }
        else if (f == "-ms" && hv() && is_number_prefix(a[i + 1].c_str(), &ms)) { o.monostretch = ms; i += 2; }
        else if (f == "-j" && hv()) { o.jobs = atoi(a[i + 1].c_str()); i += 2; }
        else if (f == "-samthreads" && hv()) { o.samthreads = atoi(a[i + 1].c_str()); i += 2; }
        else if (f == "-seed" && hv()) { o.seed = atoll(a[i + 1].c_str()); o.has_seed = true; i += 2; }
        else if (f == "-wl" && hv()) {
            o.whitelist = a[i + 1];
            if (!fexists(o.whitelist)) die("Barcode whitelist " + o.whitelist + " does not exist");
            i += 2;
        }
        else if (f == "-wl_strip_gem") { o.wl_strip_gem = true; i++; }
        else if (f == "-pretrain") { o.pretrain = true; i++; }
        else if (f == "-pre_cap" && hv()) { o.pre_cap = atoi(a[i + 1].c_str()); i += 2; }
        else if (f == "-pre_wmin" && hv()) { o.pre_wmin = atof(a[i + 1].c_str()); i += 2; }
        else if (f == "-pre_tieskip") { o.pre_tieskip = true; i++; }
        else if (f == "-jitter" && hv()) { o.jitter = atof(a[i + 1].c_str()); i += 2; }
        else if (f == "-ph") { o.phase3 = true; i++; }
        else if (f == "-ph3_ratio" && hv()) { o.ph3_ratio = atof(a[i + 1].c_str()); i += 2; }
        else die("Unexpected/incomplete parameter: " + f + "\n");
    }
    if (o.sample.empty()) die(string("No sample name\n") + USAGE);
    if (o.snv_file.empty()) die(string("No VCF file\n") + USAGE);
    if (o.bams.empty()) die(string("No BAM file\n") + USAGE);
    if (o.wl_strip_gem && o.whitelist.empty()) die(string("-wl_strip_gem requires -wl <file>\n") + USAGE);
    excise::normalize_options(o);
    return o;
}

int main(int argc, char** argv) {
    fputs("Reading parameters ...\n", stderr);
    excise::Opt o = parse_args(argc, argv);
    excise::Results res;
    try {
        excise::run_excise(o, res, [](const char* m) { fputs(m, stderr); });
    } catch (const std::exception& e) {
        die(e.what());
    }
    return 0;
}
