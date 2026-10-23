// excise15.cpp -- C++/htslib port of excise1.5.pl (EXCISE v1.5).
//
// Derived from xcise.pl of XCISE (https://github.com/Vityay/XCISE).
// Please cite: Henning RH, Rust TM, Dijksterhuis K, Eggen BJL, Guryev V.
//   bioRxiv (2024). doi:10.1101/2024.08.29.610317
//
// Same flags, algorithm and output files as excise1.5.pl and excise15.py.
// Compared with excise15.py it removes the remaining non-optimiser costs:
// BAM records are parsed natively with htslib (no per-read Python calls),
// there is no JIT compilation at start-up, and outputs are written through
// large buffers. The optimiser uses the same delta scoring as excise15.py:
// a move of SNV s only re-evaluates the cells covering s, which gives the
// same score as a full rescore because the score is a sum of per-cell terms.
//
// This file is the algorithm as a library (excise::run_excise). The command-line
// tool is main.cpp; the R/Bioconductor package wraps the same file.

#include "excise_core.h"

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <zlib.h>

#include <algorithm>
#include <cstdarg>
#include <tuple>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <string_view>
#include <map>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <exception>
#include <stdexcept>

using std::string;
using excise::Opt;
using excise::Results;
using std::vector;
typedef int64_t i64;

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------
[[noreturn]] static void die(const string& m) { throw std::runtime_error(m); }

// Log messages go through a callback (stderr for the CLI, R's message stream for
// the R package). Only the calling thread logs; worker threads never do.
static excise::LogFn g_log;
static void LOGS(const char* s) { if (g_log) g_log(s); }
static void LOGS(const string& s) { LOGS(s.c_str()); }
#define LOGF(...) do { char _lb[4096]; snprintf(_lb, sizeof _lb, __VA_ARGS__); LOGS(_lb); } while (0)

static double now_s() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Perl-style number stringification (integers without decimals, else %.15g)
static string perl_num(double x) {
    char b[64];
    if (std::isfinite(x) && x == std::floor(x) && std::fabs(x) < 1e15) snprintf(b, sizeof b, "%lld", (long long)x);
    else snprintf(b, sizeof b, "%.15g", x);
    return b;
}

static string strip_gem(const string& bc) {  // s/-\d+$//
    size_t i = bc.size();
    while (i > 0 && isdigit((unsigned char)bc[i - 1])) i--;
    if (i < bc.size() && i > 0 && bc[i - 1] == '-') return bc.substr(0, i - 1);
    return bc;
}

// reads a (possibly gzipped) text file line by line; gzopen handles both
struct LineReader {
    gzFile f;
    vector<char> buf;
    explicit LineReader(const string& p) : f(gzopen(p.c_str(), "rb")), buf(1 << 16) {
        if (!f) die("Cannot open " + p);
        gzbuffer(f, 1 << 20);
    }
    ~LineReader() { if (f) gzclose(f); }
    bool next(string& line) {
        line.clear();
        while (true) {
            if (!gzgets(f, buf.data(), (int)buf.size())) return !line.empty();
            line += buf.data();
            if (!line.empty() && line.back() == '\n') { line.pop_back(); return true; }
        }
    }
};

// buffered output file
struct Out {
    FILE* f;
    string b;
    explicit Out(const string& p) : f(fopen(p.c_str(), "w")) {
        if (!f) die("Cannot write " + p);
        b.reserve(1 << 22);
    }
    void put(const char* s, size_t n) { b.append(s, n); if (b.size() > (1 << 22)) flush(); }
    void put(const string& s) { put(s.data(), s.size()); }
    void printf(const char* fmt, ...) __attribute__((format(printf, 2, 3))) {
        char tmp[4096];
        va_list ap; va_start(ap, fmt);
        int n = vsnprintf(tmp, sizeof tmp, fmt, ap);
        va_end(ap);
        put(tmp, (size_t)std::min(n, (int)sizeof tmp - 1));
    }
    void flush() { fwrite(b.data(), 1, b.size(), f); b.clear(); }
    ~Out() { flush(); fclose(f); }
};

// ---------------------------------------------------------------------------
// input
// ---------------------------------------------------------------------------
struct SnvInfo { string rs, ref, alt; };

// one allele observation (read x SNV), in the order seen. The UMI is kept as
// bytes in an arena plus a 64-bit hash: sorting by (pos, cell, hash) and then
// by the bytes groups identical UMIs exactly, without a global UMI hash table.
struct Rec { i64 pos; uint64_t uhash; int64_t uoff, seq; int32_t cid; uint32_t ulen; int8_t allele; };
// one (pos, cell, UMI) after collapsing reads; first = first time it was seen
struct Obs { i64 pos; int32_t cid; int64_t first; int32_t a1, a2; int8_t last; };

struct BamStats { i64 inf_reads = 0, inf_alleles = 0, hard = 0, soft = 0, wl_cand = 0, wl_kept = 0, wl_rej = 0; };

static bool has_monostretch(const bam1_t* b, int n) {
    const uint8_t* s = bam_get_seq(b);
    int len = b->core.l_qseq, run = 0;
    char prev = 0;
    for (int i = 0; i < len; i++) {
        char c = seq_nt16_str[bam_seqi(s, i)];
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            run = (c == prev) ? run + 1 : 1;
            if (run >= n) return true;
        } else run = 0;
        prev = c;
    }
    return false;
}

static string cigar_string(const bam1_t* b) {
    if (b->core.n_cigar == 0) return "*";
    string s;
    const uint32_t* c = bam_get_cigar(b);
    for (uint32_t i = 0; i < b->core.n_cigar; i++) {
        s += std::to_string(bam_cigar_oplen(c[i]));
        s += bam_cigar_opchr(c[i]);
    }
    return s;
}

// ---------------------------------------------------------------------------
// optimiser (identical rules to excise1.5.pl / excise15.py)
// ---------------------------------------------------------------------------
struct Opt_tables {
    vector<i64> ptr;          // per SNV entry range
    vector<int32_t> cell, n1, n2;
};

static inline double cell_term(i64 h1, i64 h2, double pen) {
    i64 mx = h1 > h2 ? h1 : h2, mn = h1 > h2 ? h2 : h1;
    if (mx == 0) return 0.0;
    return (double)(mx - 1) - pen * (double)mn;
}

static inline double move_delta(const Opt_tables& T, int s, int old_, int new_, const vector<i64>& h1,
                                const vector<i64>& h2, double pen) {
    double d = 0;
    for (i64 k = T.ptr[s]; k < T.ptr[s + 1]; k++) {
        int c = T.cell[k];
        i64 a = h1[c], b = h2[c], na = a, nb = b;
        if (old_ == 1) { na -= T.n1[k]; nb -= T.n2[k]; }
        else if (old_ == -1) { na -= T.n2[k]; nb -= T.n1[k]; }
        if (new_ == 1) { na += T.n1[k]; nb += T.n2[k]; }
        else if (new_ == -1) { na += T.n2[k]; nb += T.n1[k]; }
        d += cell_term(na, nb, pen) - cell_term(a, b, pen);
    }
    return d;
}

static inline void move_apply(const Opt_tables& T, int s, int old_, int new_, vector<i64>& h1, vector<i64>& h2) {
    for (i64 k = T.ptr[s]; k < T.ptr[s + 1]; k++) {
        int c = T.cell[k];
        if (old_ == 1) { h1[c] -= T.n1[k]; h2[c] -= T.n2[k]; }
        else if (old_ == -1) { h1[c] -= T.n2[k]; h2[c] -= T.n1[k]; }
        if (new_ == 1) { h1[c] += T.n1[k]; h2[c] += T.n2[k]; }
        else if (new_ == -1) { h1[c] += T.n2[k]; h2[c] += T.n1[k]; }
    }
}

struct TryResult { double init = 0, best = 0; vector<int8_t> dirs; vector<std::pair<int, double>> passes; string log; };

static TryResult hill_climb(const Opt_tables& T, const vector<char>& mono, int n_cells, double pen,
                            vector<int>& order, vector<int8_t>& dirs, std::mt19937_64& rng) {
    TryResult R;
    vector<i64> h1(n_cells, 0), h2(n_cells, 0);
    int n = (int)order.size();
    for (int e = 0; e < n; e++) if (dirs[e] != 0) move_apply(T, order[e], 0, dirs[e], h1, h2);
    double score = 0;
    for (int c = 0; c < n_cells; c++) score += cell_term(h1[c], h2[c], pen);
    R.init = score;
    std::uniform_int_distribution<int> U(0, std::max(n - 1, 0));
    int imps = 1, last_imp = 0;
    static const int cands[3] = {-1, 0, 1};
    while (imps > 0) {
        imps = 0;
        for (int ele = 0; ele < n; ele++) {
            int s = order[ele], old_ = dirs[ele];
            if (old_ == 0) {
                double dm = move_delta(T, s, 0, -1, h1, h2, pen);
                double dp = move_delta(T, s, 0, 1, h1, h2, pen);
                int bn; double bd;
                if (dp >= dm) { bn = 1; bd = dp; } else { bn = -1; bd = dm; }
                if (bd > 0) {
                    move_apply(T, s, 0, bn, h1, h2);
                    dirs[ele] = (int8_t)bn; score += bd; imps++; last_imp = ele;
                }
            } else {
                for (int nw : cands) {
                    if (nw == old_) continue;
                    if (nw == 0 && mono[s]) continue;
                    double d = move_delta(T, s, old_, nw, h1, h2, pen);
                    if (d > 0 || (d == 0 && nw == 0)) {
                        move_apply(T, s, old_, nw, h1, h2);
                        imps++; last_imp = ele; dirs[ele] = (int8_t)nw; old_ = nw; score += d;
                    }
                }
            }
            // early exit with one random probe (same as excise1.5.pl)
            if (imps == 0 && last_imp - ele == 1) {
                int rnd = U(rng);
                if (rnd != ele && rnd != ele - 1) {
                    int rs = order[rnd], rold = dirs[rnd];
                    for (int rn : cands) {
                        if (rn == rold) continue;
                        if (rn == 0 && mono[rs]) continue;
                        double d = move_delta(T, rs, rold, rn, h1, h2, pen);
                        if (d > 0 || (d == 0 && rn == 0)) {
                            move_apply(T, rs, rold, rn, h1, h2);
                            imps++; last_imp = rnd; dirs[rnd] = (int8_t)rn; score += d;
                            break;
                        }
                    }
                }
                if (imps == 0) break;
            }
        }
        R.passes.push_back({imps, score});
    }
    R.best = score;
    return R;
}

// ---------------------------------------------------------------------------
// helpers on the per-SNV tables
// ---------------------------------------------------------------------------
static void hap_totals(const vector<int8_t>& dirs, const vector<int32_t>& esnv, const vector<int32_t>& ecell,
                       const vector<int32_t>& x1, const vector<int32_t>& x2, int n_cells,
                       vector<i64>& h1, vector<i64>& h2) {
    h1.assign(n_cells, 0); h2.assign(n_cells, 0);
    for (size_t k = 0; k < esnv.size(); k++) {
        int d = dirs[esnv[k]];
        if (d == 1) { h1[ecell[k]] += x1[k]; h2[ecell[k]] += x2[k]; }
        else if (d == -1) { h1[ecell[k]] += x2[k]; h2[ecell[k]] += x1[k]; }
    }
}

static vector<int8_t> cell_hap_from(const vector<i64>& h1, const vector<i64>& h2) {
    vector<int8_t> r(h1.size(), 0);
    for (size_t c = 0; c < h1.size(); c++) {
        i64 t = h1[c] + h2[c];
        if (t == 0) continue;
        if ((double)h1[c] / t >= 0.9) r[c] = 1;
        else if ((double)h2[c] / t >= 0.9) r[c] = 2;
    }
    return r;
}

static int classify1(i64 h1, i64 h2) {  // 0 Unknown, 1 X1, 2 X2, 3 Both, 4 Low_coverage
    i64 t = h1 + h2;
    if (t == 0) return 0;
    if (h1 >= 2 && (double)h1 / t >= 0.9) return 1;
    if (h2 >= 2 && (double)h2 / t >= 0.9) return 2;
    if (h1 > 0 && h2 > 0) return 3;
    return 4;
}
namespace excise {
const char* const LABELS[5] = {"Unknown", "X1", "X2", "Both", "Low_coverage"};

void normalize_options(Opt& o) {
    if (o.snv_file.empty()) throw std::runtime_error("No VCF file");
    if (o.bams.empty()) throw std::runtime_error("No BAM file");
    if (o.wl_strip_gem && o.whitelist.empty()) throw std::runtime_error("-wl_strip_gem requires -wl <file>");
    if (o.chromosome.empty()) o.chromosome = "X";
    if (o.tries <= 0) o.tries = 100;
    if (o.min_umis <= 0) o.min_umis = 10;
    if (o.min_maf < 0 || o.min_maf > 0.5) o.min_maf = 0;
    if (o.penalty < 1) o.penalty = 5;
    if (o.jobs < 1) o.jobs = 1;
    if (o.pre_cap < 100) o.pre_cap = 100;
    if (o.pre_wmin < 0) o.pre_wmin = 0;
    o.jitter = std::min(std::max(o.jitter, 0.0), 1.0);
    if (o.samthreads < 0) o.samthreads = o.jobs;
    if (o.samthreads < 1) o.samthreads = 1;
}
}  // namespace excise
using excise::LABELS;

// ---------------------------------------------------------------------------
// library entry point (the CLI is in main.cpp, the R interface in the R package)
// ---------------------------------------------------------------------------
namespace excise {

int run_excise(const Opt& o, Results& res, LogFn log) {
    g_log = log;
    res = Results();
    if (o.write_files && o.sample.empty()) die("No sample name (output prefix)");
    double t_start = now_s();
    std::mt19937_64 rng(o.has_seed ? (uint64_t)o.seed : std::random_device{}());
    const string& chrom = o.chromosome;
    LOGF("Sample             : %s\n", o.sample.c_str());
    LOGF("VCF file           : %s\n", o.snv_file.c_str());
    {
        string j; for (size_t i = 0; i < o.bams.size(); i++) j += (i ? "," : "") + o.bams[i];
        LOGF("BAM files          : %zu (%s)\n", o.bams.size(), j.c_str());
    }
    LOGF("Chromosome         : %s\n", chrom.c_str());
    LOGF("Tries              : %d\n", o.tries);
    LOGF("Min UMIs/SNV       : %s\n", perl_num(o.min_umis).c_str());
    LOGF("Min MAF            : %s\n", perl_num(o.min_maf).c_str());
    LOGF("Discordant penalty : %s\n", perl_num(o.penalty).c_str());
    LOGF("Parallel jobs      : %d\n", o.jobs);
    LOGF("samtools threads   : %d\n", o.samthreads);
    LOGF("Barcode whitelist  : %s\n", o.whitelist.empty() ? "disabled" : o.whitelist.c_str());
    if (!o.whitelist.empty()) LOGF("Strip GEM suffix   : %s\n", o.wl_strip_gem ? "yes" : "no");
    LOGF("Pretrain           : %s\n", o.pretrain ? "yes" : "no");
    LOGF("Jitter scale       : %s\n", perl_num(o.jitter).c_str());
    LOGF("Phase 3            : %s\n", o.phase3 ? "enabled" : "disabled");

    // ---- whitelist ----
    std::unordered_set<string> allowed;
    i64 wl_nonempty = 0;
    bool use_wl = !o.whitelist.empty();
    if (use_wl) {
        LOGS("Reading barcode whitelist ...\n");
        LineReader lr(o.whitelist);
        string line;
        while (lr.next(line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            size_t p = 0;
            while (p < line.size() && isspace((unsigned char)line[p])) p++;
            if (p == line.size()) continue;
            size_t q = p;
            while (q < line.size() && !isspace((unsigned char)line[q])) q++;
            string bc = line.substr(p, q - p);
            wl_nonempty++;
            if (o.wl_strip_gem) bc = strip_gem(bc);
            if (bc.empty()) continue;
            allowed.insert(bc);
        }
        if (allowed.empty()) die("Zero barcodes loaded from whitelist " + o.whitelist + "\n");
        LOGF("Whitelist barcodes loaded: %zu unique from %lld non-empty records\n", allowed.size(),
                (long long)wl_nonempty);
    }

    // ---- SNVs ----
    LOGS("Reading SNVs ...\n");
    std::unordered_map<i64, SnvInfo> snp;
    {
        LineReader lr(o.snv_file);
        string line;
        while (lr.next(line)) {
            if (line.empty() || line[0] == '#') continue;
            vector<string> f;
            size_t p = 0;
            for (int k = 0; k < 5; k++) {
                size_t q = line.find('\t', p);
                f.push_back(line.substr(p, q == string::npos ? string::npos : q - p));
                if (q == string::npos) break;
                p = q + 1;
            }
            if (f.size() < 5 || f[0] != chrom) continue;
            snp[atoll(f[1].c_str())] = SnvInfo{f[2], f[3], f[4]};
        }
    }
    LOGF("%zu SNVs loaded from %s\n", snp.size(), o.snv_file.c_str());
    if (snp.empty()) die("Zero SNVs, run terminated\n");

    // ---- BAMs ----
    // Each BAM is split into genomic chunks parsed in parallel (own file handle
    // per worker). A read belongs to the chunk holding its start position, and
    // chunk results are concatenated in genomic order, so the observations,
    // their order and all counters equal a single sequential pass.
    LOGS("Reading BAM files ...\n");
    vector<Rec> recs;
    // cell-barcode keys point into deque storage (stable addresses)
    std::unordered_map<std::string_view, int32_t> cb_index;
    string umi_arena;
    std::deque<string> cb_list;
    BamStats st;
    std::unordered_set<string> wl_seen, wl_matched;
    struct Chunk {
        hts_pos_t beg = 0, end = 0;
        vector<Rec> recs;                 // cid = local cell id, uoff into arena, seq local
        string arena;
        vector<string> cbs;               // local cell ids, first-seen order
        BamStats st;
        std::unordered_set<string> wl_seen, wl_matched;
    };
    auto parse_chunk = [&](const string& path, Chunk& C) {
        samFile* fp = sam_open(path.c_str(), "r");
        if (!fp) die("Cannot open BAM " + path);
        sam_hdr_t* hdr = sam_hdr_read(fp);
        hts_idx_t* idx = sam_index_load(fp, path.c_str());
        if (!hdr || !idx) die("Cannot read header/index of " + path + " (run samtools index)");
        int tid = sam_hdr_name2tid(hdr, chrom.c_str());
        hts_itr_t* itr = sam_itr_queryi(idx, tid, C.beg, C.end);
        std::unordered_map<std::string_view, int32_t> lcb;
        std::deque<string> lcb_store;
        bam1_t* b = bam_init1();
        while (itr && sam_itr_next(fp, itr, b) >= 0) {
            if (b->core.pos < C.beg) continue;          // started in the previous chunk
            if (b->core.flag & BAM_FSECONDARY) continue;
            uint8_t* vG = bam_aux_get(b, "vG");
            if (!vG) continue;
            uint8_t* vW = bam_aux_get(b, "vW");
            if (!vW || bam_aux2i(vW) != 1) continue;
            if (vG[0] != 'B') continue;
            if (o.monostretch > 0 && has_monostretch(b, o.monostretch)) continue;
            if (o.no_softmasking) {
                const uint32_t* c = bam_get_cigar(b);
                bool soft = false;
                for (uint32_t i = 0; i < b->core.n_cigar; i++) if (bam_cigar_op(c[i]) == BAM_CSOFT_CLIP) soft = true;
                if (soft) continue;
            }
            uint8_t* vA = bam_aux_get(b, "vA");
            uint32_t nvG = bam_auxB_len(vG), nvA = (vA && vA[0] == 'B') ? bam_auxB_len(vA) : 0;
            const char* cbp = nullptr;
            uint8_t* t = bam_aux_get(b, "CB");
            if (t) cbp = bam_aux2Z(t);
            else if ((t = bam_aux_get(b, "RG"))) cbp = bam_aux2Z(t);
            else die(string("No CB/RG tag in read: ") + bam_get_qname(b));
            if (!cbp) continue;
            std::string_view cb(cbp);
            if (cb == "-") continue;
            if (use_wl) {
                C.st.wl_cand++;
                string m = o.wl_strip_gem ? strip_gem(string(cb)) : string(cb);
                C.wl_seen.insert(m);
                if (!allowed.count(m)) { C.st.wl_rej++; continue; }
                C.st.wl_kept++;
                C.wl_matched.insert(m);
            }
            std::string_view umi;
            string soft_umi;
            uint8_t* ub = bam_aux_get(b, "UB");
            if (ub) {
                const char* u = bam_aux2Z(ub);
                umi = u ? u : "";
                if (umi == "-") continue;
                C.st.hard++;
            } else {
                soft_umi = string(cb) + "_" + sam_hdr_tid2name(hdr, b->core.tid) + "_" + std::to_string(b->core.pos + 1) + "_" +
                      std::to_string(b->core.flag) + "_" + cigar_string(b) + "_" + std::to_string(b->core.isize);
                umi = soft_umi;
                C.st.soft++;
            }
            C.st.inf_reads++;
            int32_t cid = -1;
            int64_t uoff = -1;
            uint64_t uhash = 0;
            for (uint32_t k = 0; k < nvG; k++) {
                i64 pos = bam_auxB2i(vG, k) + 1;
                if (!snp.count(pos)) continue;
                if (k >= nvA) continue;
                int allele = (int)bam_auxB2i(vA, k);
                if (allele != 1 && allele != 2) continue;
                if (o.scramble) allele = 1 + (int)(rng() & 1);   // single chunk when -x
                C.st.inf_alleles++;
                if (cid < 0) {
                    auto it = lcb.find(cb);
                    if (it == lcb.end()) {
                        cid = (int32_t)C.cbs.size(); C.cbs.emplace_back(cb);
                        lcb_store.emplace_back(cb); lcb.emplace(lcb_store.back(), cid);
                    } else cid = it->second;
                    uoff = (int64_t)C.arena.size();
                    C.arena.append(umi.data(), umi.size());
                    uhash = std::hash<std::string_view>{}(umi);
                }
                C.recs.push_back(Rec{pos, uhash, uoff, (int64_t)C.recs.size(), cid, (uint32_t)umi.size(), (int8_t)allele});
            }
        }
        bam_destroy1(b); if (itr) hts_itr_destroy(itr); hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
    };
    vector<string> bams = o.bams;
    std::sort(bams.begin(), bams.end());
    for (const string& path : bams) {
        LOGF("    Reading %s ...\n", path.c_str());
        hts_pos_t clen;
        {
            samFile* fp = sam_open(path.c_str(), "r");
            if (!fp) die("Cannot open BAM " + path);
            sam_hdr_t* hdr = sam_hdr_read(fp);
            if (!hdr) die("Cannot read header of " + path);
            int tid = sam_hdr_name2tid(hdr, chrom.c_str());
            if (tid < 0) {
                LOGF("    [warn] %s not in BAM header of %s; skipped\n", chrom.c_str(), path.c_str());
                sam_hdr_destroy(hdr); sam_close(fp); continue;
            }
            hts_idx_t* idx = sam_index_load(fp, path.c_str());
            if (!idx) die("Cannot load index for " + path + " (run samtools index)");
            clen = sam_hdr_tid2len(hdr, tid);
            hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
        }
        // more chunks than threads for load balance; one chunk when -x (RNG order)
        int nthreads = o.scramble ? 1 : std::max(1, o.samthreads);
        int nchunks = o.scramble ? 1 : nthreads * 8;
        vector<Chunk> chunks(nchunks);
        for (int k = 0; k < nchunks; k++) {
            chunks[k].beg = clen * k / nchunks;
            chunks[k].end = (k == nchunks - 1) ? HTS_POS_MAX : clen * (k + 1) / nchunks;
        }
        std::atomic<int> next_chunk(0);
        std::exception_ptr err;
        std::mutex err_mu;
        vector<std::thread> pool;
        for (int w = 0; w < nthreads; w++)
            pool.emplace_back([&]() {
                try {
                    for (int k; (k = next_chunk.fetch_add(1)) < nchunks;) parse_chunk(path, chunks[k]);
                } catch (...) {
                    std::lock_guard<std::mutex> g(err_mu);
                    if (!err) err = std::current_exception();
                    next_chunk.store(nchunks);
                }
            });
        for (auto& th : pool) th.join();
        if (err) std::rethrow_exception(err);
        // merge in genomic order: global cell ids in first-seen order
        for (Chunk& C : chunks) {
            vector<int32_t> gid(C.cbs.size());
            for (size_t l = 0; l < C.cbs.size(); l++) {
                auto it = cb_index.find(C.cbs[l]);
                if (it == cb_index.end()) {
                    gid[l] = (int32_t)cb_list.size(); cb_list.emplace_back(C.cbs[l]);
                    cb_index.emplace(cb_list.back(), gid[l]);
                } else gid[l] = it->second;
            }
            int64_t aoff = (int64_t)umi_arena.size(), soff = (int64_t)recs.size();
            umi_arena += C.arena;
            for (Rec r : C.recs) { r.cid = gid[r.cid]; r.uoff += aoff; r.seq += soff; recs.push_back(r); }
            st.inf_reads += C.st.inf_reads; st.inf_alleles += C.st.inf_alleles;
            st.hard += C.st.hard; st.soft += C.st.soft;
            st.wl_cand += C.st.wl_cand; st.wl_kept += C.st.wl_kept; st.wl_rej += C.st.wl_rej;
            wl_seen.insert(C.wl_seen.begin(), C.wl_seen.end());
            wl_matched.insert(C.wl_matched.begin(), C.wl_matched.end());
            vector<Rec>().swap(C.recs); string().swap(C.arena);
        }
        LOGF("    Reads with allelic info: %lld  Alleles: %lld  UMIs: %lld/%lld (hard/soft)\n",
                (long long)st.inf_reads, (long long)st.inf_alleles, (long long)st.hard, (long long)st.soft);
    }
    LOGF("    [time] BAM reading: %.2f s\n", now_s() - t_start);

    int n_cells = (int)cb_list.size();
    if (use_wl) {
        LOGF("Whitelist candidate barcodes observed in BAM: %zu\n", wl_seen.size());
        LOGF("Whitelist barcodes matched in BAM: %zu / %zu loaded\n", wl_matched.size(), allowed.size());
        LOGF("Whitelist CB-tagged candidate reads: checked=%lld kept=%lld rejected=%lld\n",
                (long long)st.wl_cand, (long long)st.wl_kept, (long long)st.wl_rej);
        LOGF("Whitelist barcodes with usable chr%s allelic evidence: %d / %zu loaded\n", chrom.c_str(),
                n_cells, allowed.size());
    }

    // ---- collapse reads to (pos, cell, UMI) observations ----
    auto umi_of = [&](const Rec& r) { return std::string_view(umi_arena.data() + r.uoff, r.ulen); };
    std::sort(recs.begin(), recs.end(), [&](const Rec& x, const Rec& y) {
        if (x.pos != y.pos) return x.pos < y.pos;
        if (x.cid != y.cid) return x.cid < y.cid;
        if (x.uhash != y.uhash) return x.uhash < y.uhash;
        int c = umi_of(x).compare(umi_of(y));
        if (c != 0) return c < 0;
        return x.seq < y.seq;
    });
    auto same_umi = [&](const Rec& x, const Rec& y) {
        return x.pos == y.pos && x.cid == y.cid && x.uhash == y.uhash && umi_of(x) == umi_of(y);
    };
    vector<Obs> obs;
    for (size_t i = 0; i < recs.size();) {
        Obs ob{recs[i].pos, recs[i].cid, recs[i].seq, 0, 0, 0};
        size_t j = i;
        for (; j < recs.size() && same_umi(recs[j], recs[i]); j++) {
            if (recs[j].allele == 1) ob.a1++; else ob.a2++;
            ob.last = recs[j].allele;  // last allele seen wins (as %cb2allele)
        }
        obs.push_back(ob);
        i = j;
    }
    vector<Rec>().swap(recs);
    string().swap(umi_arena);
    // per (pos, cell) groups, UMIs in the order they were first seen
    std::sort(obs.begin(), obs.end(), [](const Obs& x, const Obs& y) {
        if (x.pos != y.pos) return x.pos < y.pos;
        if (x.cid != y.cid) return x.cid < y.cid;
        return x.first < y.first;
    });
    struct Grp { i64 pos; int32_t cell, n1, n2, c1, c2, r1, r2; };
    vector<Grp> grp;
    for (size_t i = 0; i < obs.size();) {
        size_t j = i;
        const Obs& f = obs[i];
        Grp g{f.pos, f.cid, 0, 0, 0, 0, 0, 0};
        bool got1 = false, got2 = false;
        while (j < obs.size() && obs[j].pos == f.pos && obs[j].cid == f.cid) {
            const Obs& r = obs[j];
            if (r.a1 > 0) { g.n1++; if (!got1) { g.r1 = r.a1; got1 = true; } }
            if (r.a2 > 0) { g.n2++; if (!got2) { g.r2 = r.a2; got2 = true; } }
            if (r.last == 1) g.c1++; else if (r.last == 2) g.c2++;
            j++;
        }
        grp.push_back(g);
        i = j;
    }
    vector<Obs>().swap(obs);

    vector<i64> all_pos;
    for (const Grp& g : grp) if (all_pos.empty() || all_pos.back() != g.pos) all_pos.push_back(g.pos);
    LOGF("%zu SNV positions covered in %d cells\n", all_pos.size(), n_cells);
    if (all_pos.empty()) die("Zero SNVs after BAM reading\n");
    if (n_cells == 0) die("Zero cells after BAM reading\n");

    // ---- SNV filters ----
    LOGS("Checking AFs after WASP genotyping ...\n");
    size_t NP = all_pos.size();
    vector<i64> N1(NP, 0), N2(NP, 0);
    vector<int32_t> gsnv(grp.size());
    {
        size_t s = 0;
        for (size_t k = 0; k < grp.size(); k++) {
            while (all_pos[s] != grp[k].pos) s++;
            gsnv[k] = (int32_t)s;
            N1[s] += grp[k].n1; N2[s] += grp[k].n2;
        }
    }
    vector<int32_t> snv_map(NP, -1);
    vector<i64> phased_pos, pN1, pN2;
    vector<double> ratios;
    i64 low_umis = 0, low_maf = 0;
    for (size_t s = 0; s < NP; s++) {
        i64 tot = N1[s] + N2[s];
        if ((double)tot < o.min_umis) { low_umis++; continue; }
        i64 mn = std::min(N1[s], N2[s]);
        ratios.push_back(tot > 0 ? 100.0 * mn / tot : 0.0);
        double af = tot > 0 ? (double)N1[s] / tot : 0.0;
        if (af < o.min_maf || af > 1 - o.min_maf) { low_maf++; continue; }
        snv_map[s] = (int32_t)phased_pos.size();
        phased_pos.push_back(all_pos[s]); pN1.push_back(N1[s]); pN2.push_back(N2[s]);
    }
    int n_snv = (int)phased_pos.size();
    if (n_snv == 0) die("No SNVs to phase\n");
    vector<char> mono(n_snv);
    int n_mono = 0;
    for (int s = 0; s < n_snv; s++) { mono[s] = (pN1[s] == 0 || pN2[s] == 0); n_mono += mono[s]; }
    LOGF("Excluded %lld SNVs < %s UMIs; %lld SNVs with MAF < %s\n", (long long)low_umis,
            perl_num(o.min_umis).c_str(), (long long)low_maf, perl_num(o.min_maf).c_str());
    LOGF("Monoallelic SNVs in pool: %d / %d\n", n_mono, n_snv);
    std::sort(ratios.begin(), ratios.end());
    double median_maf = 0;
    if (!ratios.empty())
        median_maf = ratios.size() % 2 ? ratios[ratios.size() / 2]
                                       : (ratios[ratios.size() / 2 - 1] + ratios[ratios.size() / 2]) / 2;
    LOGF("Median MAF: %s\n", perl_num(median_maf).c_str());

    // entries for phased SNVs, sorted by SNV then cell (groups already are)
    vector<int32_t> esnv, ecell, en1, en2, ec1, ec2, er1, er2;
    for (size_t k = 0; k < grp.size(); k++) {
        int s = snv_map[gsnv[k]];
        if (s < 0) continue;
        esnv.push_back(s); ecell.push_back(grp[k].cell);
        en1.push_back(grp[k].n1); en2.push_back(grp[k].n2);
        ec1.push_back(grp[k].c1); ec2.push_back(grp[k].c2);
        er1.push_back(grp[k].r1); er2.push_back(grp[k].r2);
    }
    vector<Grp>().swap(grp);
    Opt_tables T;
    T.ptr.assign(n_snv + 1, 0);
    for (int s : esnv) T.ptr[s + 1]++;
    for (int s = 0; s < n_snv; s++) T.ptr[s + 1] += T.ptr[s];
    T.cell = ecell; T.n1 = en1; T.n2 = en2;

    LOGF("    [time] tables built: %.2f s\n", now_s() - t_start);
    // ---- pretraining ----
    vector<int8_t> pre_seed(n_snv, 0);
    vector<double> unc(n_snv, 1.0);
    if (o.pretrain) {
        LOGS("Running improved pretraining ...\n");
        LOGS("  [pretrain] Building cell majority calls ...\n");
        // per-cell (snv, call) lists
        vector<vector<std::pair<int32_t, int8_t>>> calls(n_cells);
        for (size_t k = 0; k < esnv.size(); k++) {
            int c1 = ec1[k], c2 = ec2[k];
            if (c1 + c2 <= 0) continue;
            if (o.pre_tieskip && c1 == c2) continue;
            calls[ecell[k]].push_back({esnv[k], (int8_t)(c1 >= c2 ? 1 : 2)});
        }
        LOGS("  [pretrain] Building co-segregation graph ...\n");
        vector<double> deg(n_snv, 0);
        vector<std::pair<uint64_t, int8_t>> votes;  // key = i*n_snv + j (i < j), sign
        for (int c = 0; c < n_cells; c++) {
            auto& v = calls[c];
            i64 L = (i64)v.size();
            if (L < 2) continue;
            auto vote = [&](int i, int j) {
                int si = v[i].first, sj = v[j].first;
                if (si > sj) std::swap(si, sj);
                votes.push_back({(uint64_t)si * n_snv + sj, (int8_t)((v[i].second == v[j].second) ? 1 : -1)});
                deg[si] += 1; deg[sj] += 1;
            };
            if (L * (L - 1) / 2 <= o.pre_cap) {
                for (int i = 0; i < L; i++) for (int j = i + 1; j < L; j++) vote(i, j);
            } else {
                std::unordered_set<i64> seen;
                std::uniform_int_distribution<int> U(0, (int)L - 1);
                i64 tr = 0, got = 0;
                while (got < o.pre_cap && tr < 5LL * o.pre_cap) {
                    tr++;
                    int i = U(rng), j = U(rng);
                    if (i == j) continue;
                    if (i > j) std::swap(i, j);
                    if (!seen.insert((i64)i * L + j).second) continue;
                    vote(i, j);
                    got++;
                }
            }
        }
        vector<vector<std::pair<int32_t, int8_t>>>().swap(calls);
        // sum votes per SNV pair
        std::sort(votes.begin(), votes.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
        vector<uint64_t> ekey;
        vector<double> ew;
        for (size_t k = 0; k < votes.size(); k++) {
            if (ekey.empty() || ekey.back() != votes[k].first) { ekey.push_back(votes[k].first); ew.push_back(0); }
            ew.back() += votes[k].second;
        }
        vector<std::pair<uint64_t, int8_t>>().swap(votes);
        // adjacency (both directions). Edges are visited in (i, j) key order, so
        // every neighbour list comes out sorted by neighbour index.
        vector<i64> aptr(n_snv + 1, 0);
        for (uint64_t k : ekey) { aptr[k / n_snv + 1]++; aptr[k % n_snv + 1]++; }
        for (int s = 0; s < n_snv; s++) aptr[s + 1] += aptr[s];
        vector<int32_t> anbr(aptr[n_snv]);
        vector<double> aw(aptr[n_snv]);
        {
            vector<i64> fill(aptr.begin(), aptr.end() - 1);
            for (size_t k = 0; k < ekey.size(); k++) {
                int i = (int)(ekey[k] / n_snv), j = (int)(ekey[k] % n_snv);
                anbr[fill[i]] = j; aw[fill[i]++] = ew[k];
                anbr[fill[j]] = i; aw[fill[j]++] = ew[k];
            }
        }
        LOGF("  [pretrain] BFS seeding (wmin=%s) ...\n", perl_num(o.pre_wmin).c_str());
        vector<int> cand(n_snv);
        std::iota(cand.begin(), cand.end(), 0);
        std::stable_sort(cand.begin(), cand.end(), [&](int x, int y) { return deg[x] > deg[y]; });
        vector<char> seen(n_snv, 0);
        vector<double> bott(n_snv, 0);
        vector<int> q(n_snv);
        for (int start : cand) {
            if (seen[start] || aptr[start + 1] == aptr[start]) continue;
            seen[start] = 1; pre_seed[start] = 1; bott[start] = 1e9;
            int qh = 0, qt = 0;
            q[qt++] = start;
            while (qh < qt) {
                int u = q[qh++];
                for (i64 k = aptr[u]; k < aptr[u + 1]; k++) {
                    int v = anbr[k];
                    if (seen[v]) continue;
                    double w = aw[k];
                    if (std::fabs(w) < o.pre_wmin) continue;
                    seen[v] = 1;
                    pre_seed[v] = (int8_t)(w >= 0 ? pre_seed[u] : -pre_seed[u]);
                    bott[v] = std::min(std::fabs(w), bott[u]);
                    q[qt++] = v;
                }
            }
        }
        i64 strong = 0;
        for (double w : aw) strong += std::fabs(w) >= o.pre_wmin;
        int assigned = 0;
        for (int s = 0; s < n_snv; s++) {
            if (!seen[s]) bott[s] = 0;
            unc[s] = 1.0 / (1.0 + bott[s]);
            assigned += pre_seed[s] != 0;
        }
        LOGF("  [pretrain] nodes=%d  edges~%zu  strong-edges~%lld  assigned=%d (%.1f%%)\n", n_snv,
                anbr.size(), (long long)strong, assigned, n_snv ? 100.0 * assigned / n_snv : 0.0);
    }

    // ---- parallel tries ----
    LOGF("    [time] pretraining done: %.2f s\n", now_s() - t_start);
    LOGS("Starting XCI calling\n");
    double t_opt = now_s();
    i64 time_seed = (i64)time(nullptr);
    vector<TryResult> results(o.tries + 1);
    std::atomic<int> next_try(1);
    auto worker = [&]() {
        while (true) {
            int t = next_try.fetch_add(1);
            if (t > o.tries) return;
            uint64_t seed = (uint64_t)((o.has_seed ? o.seed : time_seed) + (i64)t * 1337);
            std::mt19937_64 r(seed);
            std::uniform_real_distribution<double> U01(0.0, 1.0);
            vector<int> lpos(n_snv);
            std::iota(lpos.begin(), lpos.end(), 0);
            vector<int8_t> ldir(n_snv, 0);
            if (o.pretrain) {
                for (int i = 0; i < n_snv; i++) {
                    int8_t d = pre_seed[i];
                    bool flip = U01(r) < unc[i] * o.jitter;
                    int8_t rs = U01(r) < 0.5 ? 1 : -1;
                    ldir[i] = flip ? (d == 0 ? rs : (int8_t)-d) : d;
                }
                for (int i = n_snv - 1; i > 0; i--) {  // Fisher-Yates, dirs move with positions
                    int j = std::uniform_int_distribution<int>(0, i)(r);
                    if (i != j) { std::swap(lpos[i], lpos[j]); std::swap(ldir[i], ldir[j]); }
                }
            } else {
                std::uniform_int_distribution<int> UN(0, n_snv - 1), U3(0, 2);
                for (int i = 0; i < n_snv; i++) {
                    int a = UN(r), b = UN(r);
                    if (a != b) std::swap(lpos[a], lpos[b]);
                    if (mono[lpos[i]]) ldir[i] = U01(r) < 0.5 ? 1 : -1;
                    else ldir[i] = (int8_t)(U3(r) - 1);
                }
            }
            TryResult R = hill_climb(T, mono, n_cells, o.penalty, lpos, ldir, r);
            R.dirs.assign(n_snv, 0);
            for (int i = 0; i < n_snv; i++) R.dirs[lpos[i]] = ldir[i];
            {
                char lb[256];
                snprintf(lb, sizeof lb, "Try #%d, initialising ...\n    SNVs: %d  Initial score: %s\n", t, n_snv,
                         perl_num(R.init).c_str());
                R.log += lb;
                for (size_t k = 0; k < R.passes.size(); k++) {
                    snprintf(lb, sizeof lb, "    Pass: %zu  Improvements: %d  Score: %s\n", k + 1, R.passes[k].first,
                             perl_num(R.passes[k].second).c_str());
                    R.log += lb;
                }
                snprintf(lb, sizeof lb, "    Finished Try #%d  Final score: %s\n", t, perl_num(R.best).c_str());
                R.log += lb;
            }
            results[t] = std::move(R);
        }
    };
    {
        vector<std::thread> pool;
        int nt = std::min(o.jobs, o.tries);
        for (int k = 0; k < nt; k++) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
    }
    for (int t = 1; t <= o.tries; t++) LOGS(results[t].log);
    LOGF("    Optimiser wall time: %.2f s\n", now_s() - t_opt);

    double global_best = 0;
    bool have_best = false;
    vector<int8_t> dirs;
    for (int t = 1; t <= o.tries; t++)
        if (!have_best || results[t].best > global_best) { global_best = results[t].best; dirs = results[t].dirs; have_best = true; }
    results.clear();

    string prefix = o.sample + "_chr" + chrom + "_XCISE";
    if (o.improve_existing) {
        FILE* f = fopen((prefix + "_summary.txt").c_str(), "r");
        if (f) {
            char line[512];
            long long prev;
            if (fgets(line, sizeof line, f) && sscanf(line, "Best score : %lld", &prev) == 1 && (global_best == 0 || prev > global_best)) {
                global_best = (double)prev;
                LOGF("    Adopted best score from previous run: %lld\n", prev);
            }
            fclose(f);
        }
    }

#ifdef EXCISE_TEST_INJECT
    // test builds only: replace the optimiser result with directions from a file
    // (first line: best score; then "pos<TAB>dir"), to compare post-optimiser stages
    if (const char* inj = getenv("INJECT")) {
        FILE* f = fopen(inj, "r");
        long long bs = 0;
        if (fscanf(f, "%lld", &bs) != 1) die("bad INJECT file");
        global_best = (double)bs;
        std::unordered_map<i64, int> m; long long p; int d;
        while (fscanf(f, "%lld %d", &p, &d) == 2) m[p] = d;
        fclose(f);
        for (int s = 0; s < n_snv; s++) dirs[s] = (int8_t)m.at(phased_pos[s]);
    }
#endif

    vector<i64> h1, h2;
    // ---- Phase 2: monoallelic SNVs by cell vote ----
    LOGS("    Phase-2 cell-vote correction for monoallelic SNVs ...\n");
    {
        vector<int8_t> bi(n_snv);
        for (int s = 0; s < n_snv; s++) bi[s] = mono[s] ? 0 : dirs[s];
        hap_totals(bi, esnv, ecell, ec1, ec2, n_cells, h1, h2);
        vector<int8_t> ch = cell_hap_from(h1, h2);
        vector<i64> v1(n_snv, 0), v2(n_snv, 0);
        for (size_t k = 0; k < esnv.size(); k++) {
            int s = esnv[k];
            bool obs2 = pN2[s] > 0;
            if ((obs2 ? en2[k] : en1[k]) <= 0) continue;
            if (ch[ecell[k]] == 1) v1[s]++; else if (ch[ecell[k]] == 2) v2[s]++;
        }
        int changed = 0, kept = 0;
        for (int s = 0; s < n_snv; s++) {
            if (!mono[s] || v1[s] + v2[s] == 0) continue;
            bool obs1 = !(pN2[s] > 0);
            int correct = v1[s] >= v2[s] ? (obs1 ? 1 : -1) : (obs1 ? -1 : 1);
            if (correct != dirs[s]) { dirs[s] = (int8_t)correct; changed++; } else kept++;
        }
        LOGF("    Mono SNVs corrected: %d  already correct: %d\n", changed, kept);
    }

    // ---- Phase 3 (opt-in) ----
    if (o.phase3) {
        LOGS("    Phase-3 directionality assignment for dir=0 biallelic SNVs ...\n");
        hap_totals(dirs, esnv, ecell, en1, en2, n_cells, h1, h2);
        vector<int8_t> res(n_snv, 0);
        int assigned = 0, stayed = 0;
        for (int s = 0; s < n_snv; s++) {
            if (dirs[s] != 0 || mono[s]) continue;
            if (T.ptr[s + 1] == T.ptr[s]) continue;
            i64 nd[2] = {0, 0};
            double sh1[2] = {0, 0}, sh2[2] = {0, 0};
            for (int ti = 0; ti < 2; ti++) {
                for (i64 k = T.ptr[s]; k < T.ptr[s + 1]; k++) {
                    int c = ecell[k];
                    i64 a = ti == 0 ? h1[c] + en1[k] : h1[c] + en2[k];
                    i64 b = ti == 0 ? h2[c] + en2[k] : h2[c] + en1[k];
                    i64 tot = a + b;
                    if (tot <= 0) continue;
                    sh1[ti] += a; sh2[ti] += b;
                    if ((double)a / tot >= 0.9 || (double)b / tot >= 0.9) nd[ti]++;
                }
            }
            int wi;
            if (nd[0] > nd[1]) wi = 0;
            else if (nd[1] > nd[0]) wi = 1;
            else {
                double dpos = std::max(sh1[0], sh2[0]), dneg = std::max(sh1[1], sh2[1]);
                wi = dpos >= dneg ? 0 : 1;
            }
            double tot = sh1[wi] + sh2[wi];
            double frac = tot > 0 ? std::max(sh1[wi], sh2[wi]) / tot : 0.0;
            if (frac >= o.ph3_ratio) { res[s] = (int8_t)(wi == 0 ? 1 : -1); assigned++; }
            else stayed++;
        }
        LOGF("    Phase-3: assigned=%d  stayed Unk=%d\n", assigned, stayed);
        for (int s = 0; s < n_snv; s++) if (res[s] != 0) dirs[s] = res[s];

        LOGS("    Phase-3 cell-vote assignment for dir=0 biallelic SNVs ...\n");
        hap_totals(dirs, esnv, ecell, ec1, ec2, n_cells, h1, h2);
        vector<int8_t> ch = cell_hap_from(h1, h2);
        vector<double> cpos(n_snv, 0), cneg(n_snv, 0);
        for (size_t k = 0; k < esnv.size(); k++) {
            int s = esnv[k], hp = ch[ecell[k]];
            double r1 = en1[k] > 0 ? er1[k] : 0, r2 = en2[k] > 0 ? er2[k] : 0;
            if (hp == 1) { cpos[s] += r1; cneg[s] += r2; }
            else if (hp == 2) { cpos[s] += r2; cneg[s] += r1; }
        }
        assigned = 0; stayed = 0;
        for (int s = 0; s < n_snv; s++) {
            double tev = cpos[s] + cneg[s];
            if (dirs[s] != 0 || mono[s] || tev <= 0) continue;
            int win = cpos[s] >= cneg[s] ? 1 : -1;
            double wc = cpos[s] >= cneg[s] ? cpos[s] : cneg[s];
            if (wc / std::max(tev, 1.0) >= o.ph3_ratio) { dirs[s] = (int8_t)win; assigned++; }
            else stayed++;
        }
        LOGF("    Phase-3: assigned=%d  stayed Unk=%d\n", assigned, stayed);
    } else {
        LOGS("    Phase 3 disabled (use -ph to enable); dir=0 SNVs remain unresolved.\n");
    }

    // ---- outputs (flip so that X2 >= X1) ----
    vector<int8_t> lab(n_cells);
    i64 cnt[5];
    while (true) {
        hap_totals(dirs, esnv, ecell, ec1, ec2, n_cells, h1, h2);
        std::fill(cnt, cnt + 5, 0);
        for (int c = 0; c < n_cells; c++) { lab[c] = (int8_t)classify1(h1[c], h2[c]); cnt[lab[c]]++; }
        LOGF("    X1/X2/Both/LowC/Unk: %lld / %lld / %lld / %lld / %lld\n", (long long)cnt[1],
                (long long)cnt[2], (long long)cnt[3], (long long)cnt[4], (long long)cnt[0]);
        if (cnt[1] > cnt[2]) {
            LOGS("    Flipping X1<->X2 so that X2 >= X1 ...\n");
            for (auto& d : dirs) d = (int8_t)-d;
            continue;
        }
        break;
    }

    LOGF("    [time] phases 2/3 done: %.2f s\n", now_s() - t_start);
    vector<int> cb_order(n_cells);
    std::iota(cb_order.begin(), cb_order.end(), 0);
    std::sort(cb_order.begin(), cb_order.end(), [&](int x, int y) { return cb_list[x] < cb_list[y]; });

    // ---- results returned to the caller ----
    for (int c : cb_order) {
        res.barcode.push_back(cb_list[c]); res.hap1.push_back(h1[c]); res.hap2.push_back(h2[c]);
        res.label.push_back(lab[c]);
    }
    for (int s = 0; s < n_snv; s++) {
        const SnvInfo& si = snp.at(phased_pos[s]);
        res.chrom.push_back(chrom); res.pos.push_back(phased_pos[s]); res.id.push_back(si.rs);
        res.ref.push_back(si.ref); res.alt.push_back(si.alt); res.dir.push_back(dirs[s]);
        res.ref_umis.push_back(pN1[s]); res.alt_umis.push_back(pN2[s]); res.mono.push_back(mono[s]);
    }
    res.median_maf = median_maf;
    for (int k = 0; k < 5; k++) res.cell_counts[k] = cnt[k];
    if (o.return_profile) {
        vector<int> rank(n_cells);
        for (int i = 0; i < n_cells; i++) rank[cb_order[i]] = i;
        for (size_t k = 0; k < esnv.size(); k++) {
            if (ec1[k] + ec2[k] <= 0) continue;
            res.prof_cell.push_back(rank[ecell[k]]); res.prof_snv.push_back(esnv[k]);
            res.prof_ref.push_back(ec1[k]); res.prof_alt.push_back(ec2[k]);
        }
    }
    if (o.write_files) LOGS("    Outputting barcode-to-XCI ...\n");
    if (o.write_files) {
        Out f(prefix + "_bc2xci.txt");
        for (int c : cb_order) f.printf("%s\t%lld\t%lld\t%s\n", cb_list[c].c_str(), (long long)h1[c], (long long)h2[c], LABELS[lab[c]]);
    }

    vector<i64> b1, b2;
    hap_totals(dirs, esnv, ecell, en1, en2, n_cells, b1, b2);
    i64 total_umi = 0, disc = 0, conc = 0;
    for (int c = 0; c < n_cells; c++) {
        total_umi += b1[c] + b2[c];
        i64 mx = std::max(b1[c], b2[c]), mn = std::min(b1[c], b2[c]);
        if (mx > 0) { disc += mn; conc += mx - 1; }
    }
    int ref1 = 0, alt1 = 0, unk1 = 0;
    for (int8_t d : dirs) { ref1 += d == 1; alt1 += d == -1; unk1 += d == 0; }
    i64 grand = cnt[0] + cnt[1] + cnt[2] + cnt[3] + cnt[4];
    LOGF("    Total XCI-informative UMIs   : %lld\n", (long long)total_umi);
    LOGF("    SNVs informative/non-informative/mono : %d/%d/%d\n", ref1 + alt1 - n_mono, unk1, n_mono);
    auto pct = [&](i64 x) { return grand ? 100.0 * x / grand : 0.0; };
    res.best_score = global_best; res.total_umis = total_umi; res.concordant = conc; res.discordant = disc;
    res.n_snvs = n_snv; res.n_informative = ref1 + alt1 - n_mono; res.n_noninformative = unk1;
    res.n_monoallelic = n_mono; res.whitelist_used = use_wl;
    if (use_wl) {
        res.wl_loaded = (i64)allowed.size(); res.wl_matched = (i64)wl_matched.size();
        res.wl_reads_checked = st.wl_cand; res.wl_reads_kept = st.wl_kept; res.wl_reads_rejected = st.wl_rej;
    }
    if (o.write_files) {
        Out f(prefix + "_summary.txt");
        f.printf("Best score  : %s\n", perl_num(global_best).c_str());
        f.printf("Total UMIs  : %lld\n", (long long)total_umi);
        f.printf("Concordant  : %lld\n", (long long)conc);
        f.printf("Discordant  : %lld\n", (long long)disc);
        f.printf("Discordancy : %s\n", conc + disc ? perl_num((double)disc / (conc + disc)).c_str() : "N/A");
        f.printf("Total SNVs        : %d\n", n_snv);
        f.printf("Informative SNVs  : %d\n", ref1 + alt1 - n_mono);
        f.printf("Non-informative   : %d\n", unk1);
        f.printf("Monoallelic SNVs  : %d\n", n_mono);
        f.printf("Median MAF        : %3.2f %%\n", median_maf);
        f.printf("X1 cells    : %5lld( %3.2f %% )\n", (long long)cnt[1], pct(cnt[1]));
        f.printf("X2 cells    : %5lld( %3.2f %% )\n", (long long)cnt[2], pct(cnt[2]));
        f.printf("Both X      : %5lld( %3.2f %% )\n", (long long)cnt[3], pct(cnt[3]));
        f.printf("LowCoverage : %5lld( %3.2f %% )\n", (long long)cnt[4], pct(cnt[4]));
        f.printf("Unknown     : %5lld( %3.2f %% )\n", (long long)cnt[0], pct(cnt[0]));
        if (use_wl) {
            f.printf("Whitelist file     : %s\n", o.whitelist.c_str());
            f.printf("Whitelist strip GEM: %s\n", o.wl_strip_gem ? "yes" : "no");
            f.printf("Whitelist loaded   : %zu\n", allowed.size());
            f.printf("Whitelist BAM match: %zu\n", wl_matched.size());
            f.printf("Whitelist evidence  : %d\n", n_cells);
            f.printf("Whitelist reads     : checked=%lld kept=%lld rejected=%lld\n", (long long)st.wl_cand,
                     (long long)st.wl_kept, (long long)st.wl_rej);
        } else {
            f.printf("Barcode whitelist  : disabled\n");
        }
    }

    if (o.write_files) LOGS("    Outputting VCF ...\n");
    if (o.write_files) {
        Out f(prefix + ".vcf");
        f.put("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        f.put("##INFO=<ID=X1A,Number=1,Type=String,Description=\"X1 haplotype allele: Ref, Alt, or Unk (escape candidate)\">\n");
        f.put("##INFO=<ID=MONO,Number=0,Type=Flag,Description=\"Monoallelic SNV: direction by cell-vote not optimizer\">\n");
        f.put("##INFO=<ID=AD,Number=2,Type=Integer,Description=\"Allele depth: ref_umi,alt_umi\">\n");
        f.put("##INFO=<ID=X1,Number=1,Type=Integer,Description=\"UMI count supporting X1 haplotype\">\n");
        f.put("##INFO=<ID=X2,Number=1,Type=Integer,Description=\"UMI count supporting X2 haplotype\">\n");
        f.put("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total UMI depth (X1 + X2)\">\n");
        for (int s = 0; s < n_snv; s++) {
            int d = dirs[s];
            i64 n1 = pN1[s], n2 = pN2[s];
            i64 x1u = d == 1 ? n1 : (d == -1 ? n2 : 0), x2u = d == -1 ? n1 : (d == 1 ? n2 : 0);
            const SnvInfo& si = snp.at(phased_pos[s]);
            f.printf("%s\t%lld\t%s\t%s\t%s\t%lld\tPASS\tX1A=%s;AD=%lld,%lld;X1=%lld;X2=%lld;DP=%lld%s\n", chrom.c_str(),
                     (long long)phased_pos[s], si.rs.c_str(), si.ref.c_str(), si.alt.c_str(), (long long)(x1u + x2u),
                     d == 0 ? "Unk" : (d == 1 ? "Ref" : "Alt"), (long long)n1, (long long)n2, (long long)x1u,
                     (long long)x2u, (long long)(x1u + x2u), mono[s] ? ";MONO" : "");
        }
    }

    if (o.write_files) LOGS("    Outputting haplotype sequence table ...\n");
    if (o.write_files) {
        Out f(prefix + "_haplotypes.tsv");
        f.put("CHROM\tPOS\tID\tREF\tALT\tX1_allele\tX2_allele\tmono\n");
        for (int s = 0; s < n_snv; s++) {
            int d = dirs[s];
            if (d == 0) continue;
            const SnvInfo& si = snp.at(phased_pos[s]);
            const string& x1 = d == 1 ? si.ref : si.alt;
            const string& x2 = d == 1 ? si.alt : si.ref;
            f.printf("%s\t%lld\t%s\t%s\t%s\t%s\t%s\t%d\n", chrom.c_str(), (long long)phased_pos[s], si.rs.c_str(),
                     si.ref.c_str(), si.alt.c_str(), x1.c_str(), x2.c_str(), (int)mono[s]);
        }
    }

    if (o.write_files) LOGS("    Outputting per-cell SNV profile ...\n");
    if (o.write_files) {
        vector<int> rank(n_cells);
        for (int i = 0; i < n_cells; i++) rank[cb_order[i]] = i;
        vector<size_t> po(esnv.size());
        std::iota(po.begin(), po.end(), 0);
        std::sort(po.begin(), po.end(), [&](size_t x, size_t y) {
            if (rank[ecell[x]] != rank[ecell[y]]) return rank[ecell[x]] < rank[ecell[y]];
            return esnv[x] < esnv[y];
        });
        Out f(prefix + "_cell_snv_profile.tsv");
        f.put("cell_barcode\tpos\trs\tref\talt\tsnv_dir\tx1a\tref_umis\talt_umis\ttotal_umis\thap_x1_umis\thap_x2_umis\tcell_xci\tmono\n");
        for (size_t k : po) {
            int r = ec1[k], a = ec2[k];
            if (r + a <= 0) continue;
            int s = esnv[k], c = ecell[k], d = dirs[s];
            int hx1 = d == 1 ? r : (d == -1 ? a : 0), hx2 = d == 1 ? a : (d == -1 ? r : 0);
            const SnvInfo& si = snp.at(phased_pos[s]);
            f.printf("%s\t%lld\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\n", cb_list[c].c_str(),
                     (long long)phased_pos[s], si.rs.c_str(), si.ref.c_str(), si.alt.c_str(), d,
                     d == 0 ? "Unk" : (d == 1 ? "Ref" : "Alt"), r, a, r + a, hx1, hx2, LABELS[lab[c]], (int)mono[s]);
        }
    }
    if (o.write_files) LOGF("    Per-cell SNV profile written: %s_cell_snv_profile.tsv\n", prefix.c_str());
    LOGF("    Total wall time: %.2f s\n", now_s() - t_start);
    g_log = nullptr;
    return 0;
}

}  // namespace excise
