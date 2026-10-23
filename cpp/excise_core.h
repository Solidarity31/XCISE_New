// excise_core.h -- EXCISE 1.5 (C++ port of excise1.5.pl) as a library.
//
// Derived from xcise.pl of XCISE (https://github.com/Vityay/XCISE).
// Please cite: Henning RH, Rust TM, Dijksterhuis K, Eggen BJL, Guryev V.
//   bioRxiv (2024). doi:10.1101/2024.08.29.610317
//
// Used by the command-line tool (main.cpp) and by the R/Bioconductor package.
// run_excise() never exits or writes to the console: errors are thrown as
// std::runtime_error and messages go to the LogFn callback, which is only ever
// called from the calling thread.

#ifndef EXCISE_CORE_H
#define EXCISE_CORE_H

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace excise {

typedef std::function<void(const char*)> LogFn;

// Options; field names follow the excise1.5.pl flags (in brackets).
struct Opt {
    std::string sample;            // [-o] output prefix (files: <sample>_chr<chrom>_XCISE*)
    std::string chromosome;        // [-r] default X
    std::string snv_file;          // [-s] heterozygous SNV VCF (plain or gzipped)
    std::string whitelist;         // [-wl] optional barcode whitelist
    int tries = 100;               // [-t]
    double min_maf = 0;            // [-m]
    double min_umis = 4;           // [-u]
    double penalty = 5;            // [-p]
    bool scramble = false;         // [-x]
    bool improve_existing = false; // [-i]
    bool no_softmasking = false;   // [-sm]
    int monostretch = 15;          // [-ms]
    int jobs = 1;                  // [-j]
    int samthreads = -1;           // [-samthreads] (-1: same as jobs)
    bool has_seed = false;         // [-seed] given?
    int64_t seed = 0;              // [-seed]
    bool wl_strip_gem = false;     // [-wl_strip_gem]
    bool pretrain = false;         // [-pretrain]
    bool pre_tieskip = true;       // [-pre_tieskip]
    bool phase3 = false;           // [-ph]
    int pre_cap = 5000;            // [-pre_cap]
    double pre_wmin = 3;           // [-pre_wmin]
    double jitter = 0.25;          // [-jitter]
    double ph3_ratio = 0.90;       // [-ph3_ratio]
    std::vector<std::string> bams; // [-b]
    // library-only switches
    bool write_files = true;       // write the five excise1.5.pl output files
    bool return_profile = false;   // fill Results::prof_* (per cell x SNV UMIs)
};

// Applies excise1.5.pl's defaults and range rules (e.g. -u <= 0 becomes 10).
// Throws std::runtime_error for missing required inputs.
void normalize_options(Opt& o);

// Everything the output files contain, as vectors.
struct Results {
    // cells, sorted by barcode (= *_bc2xci.txt)
    std::vector<std::string> barcode;
    std::vector<int64_t> hap1, hap2;
    std::vector<int> label;        // 0 Unknown, 1 X1, 2 X2, 3 Both, 4 Low_coverage
    // SNVs that passed the filters, by position (= *_XCISE.vcf)
    std::vector<std::string> chrom, id, ref, alt;
    std::vector<int64_t> pos, ref_umis, alt_umis;
    std::vector<int> dir;          // 1: REF on X1, -1: ALT on X1, 0: unassigned
    std::vector<int> mono;
    // summary (= *_summary.txt)
    double best_score = 0, median_maf = 0;
    int64_t total_umis = 0, concordant = 0, discordant = 0;
    int n_snvs = 0, n_informative = 0, n_noninformative = 0, n_monoallelic = 0;
    int64_t cell_counts[5] = {0, 0, 0, 0, 0};  // indexed like label
    bool whitelist_used = false;
    int64_t wl_loaded = 0, wl_matched = 0, wl_reads_checked = 0, wl_reads_kept = 0, wl_reads_rejected = 0;
    // optional per cell x SNV UMI counts (indices into barcode / pos)
    std::vector<int> prof_cell, prof_snv, prof_ref, prof_alt;
};

// Runs EXCISE 1.5. Returns 0; throws std::runtime_error on bad input.
int run_excise(const Opt& o, Results& res, LogFn log);

extern const char* const LABELS[5];

}  // namespace excise

#endif
