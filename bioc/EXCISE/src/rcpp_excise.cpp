// R interface to the EXCISE 1.5 C++ core (excise_core.cpp).
//
// The core never calls R from its worker threads: messages reach this file
// only through the log callback, which runs on the calling (R) thread.

#include <Rcpp.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "excise_core.h"

// [[Rcpp::export(.exciseRun)]]
Rcpp::List exciseRun(Rcpp::List opts, bool verbose) {
    excise::Opt o;
    o.bams = Rcpp::as<std::vector<std::string>>(opts["bam"]);
    o.snv_file = Rcpp::as<std::string>(opts["vcf"]);
    o.chromosome = Rcpp::as<std::string>(opts["chromosome"]);
    o.whitelist = Rcpp::as<std::string>(opts["whitelist"]);
    o.wl_strip_gem = Rcpp::as<bool>(opts["whitelistStripGem"]);
    o.tries = Rcpp::as<int>(opts["tries"]);
    o.min_umis = Rcpp::as<double>(opts["minUmis"]);
    o.min_maf = Rcpp::as<double>(opts["minMaf"]);
    o.penalty = Rcpp::as<double>(opts["penalty"]);
    o.no_softmasking = Rcpp::as<bool>(opts["skipSoftClipped"]);
    o.monostretch = Rcpp::as<int>(opts["monostretch"]);
    o.pretrain = Rcpp::as<bool>(opts["pretrain"]);
    o.pre_cap = Rcpp::as<int>(opts["preCap"]);
    o.pre_wmin = Rcpp::as<double>(opts["preMinWeight"]);
    o.pre_tieskip = Rcpp::as<bool>(opts["preTieSkip"]);
    o.jitter = Rcpp::as<double>(opts["jitter"]);
    o.phase3 = Rcpp::as<bool>(opts["phase3"]);
    o.ph3_ratio = Rcpp::as<double>(opts["phase3Ratio"]);
    o.jobs = Rcpp::as<int>(opts["workers"]);
    o.samthreads = Rcpp::as<int>(opts["bamThreads"]);
    o.has_seed = true;
    o.seed = static_cast<int64_t>(Rcpp::as<double>(opts["seed"]));
    o.sample = Rcpp::as<std::string>(opts["outPrefix"]);
    o.write_files = !o.sample.empty();
    o.return_profile = Rcpp::as<bool>(opts["returnProfile"]);

    excise::Results r;
    Rcpp::Function message("message");
    excise::LogFn log;
    if (verbose) {
        log = [&message](const char* s) {
            std::string m(s);
            if (!m.empty() && m.back() == '\n') m.pop_back();
            if (!m.empty()) message(m);
        };
    }
    try {
        excise::normalize_options(o);
        excise::run_excise(o, r, log);
    } catch (const std::exception& e) {
        Rcpp::stop(e.what());
    }

    Rcpp::CharacterVector labels(r.label.size());
    for (size_t i = 0; i < r.label.size(); i++) labels[i] = excise::LABELS[r.label[i]];
    Rcpp::List cells = Rcpp::List::create(
        Rcpp::Named("barcode") = r.barcode,
        Rcpp::Named("hap1") = Rcpp::NumericVector(r.hap1.begin(), r.hap1.end()),
        Rcpp::Named("hap2") = Rcpp::NumericVector(r.hap2.begin(), r.hap2.end()),
        Rcpp::Named("call") = labels);
    Rcpp::List snvs = Rcpp::List::create(
        Rcpp::Named("chrom") = r.chrom,
        Rcpp::Named("pos") = Rcpp::NumericVector(r.pos.begin(), r.pos.end()),
        Rcpp::Named("id") = r.id,
        Rcpp::Named("ref") = r.ref,
        Rcpp::Named("alt") = r.alt,
        Rcpp::Named("dir") = r.dir,
        Rcpp::Named("refUmis") = Rcpp::NumericVector(r.ref_umis.begin(), r.ref_umis.end()),
        Rcpp::Named("altUmis") = Rcpp::NumericVector(r.alt_umis.begin(), r.alt_umis.end()),
        Rcpp::Named("mono") = r.mono);
    Rcpp::NumericVector counts(5);
    for (int k = 0; k < 5; k++) counts[k] = static_cast<double>(r.cell_counts[k]);
    counts.names() = Rcpp::CharacterVector::create("Unknown", "X1", "X2", "Both", "Low_coverage");
    Rcpp::List summary = Rcpp::List::create(
        Rcpp::Named("bestScore") = r.best_score,
        Rcpp::Named("totalUmis") = static_cast<double>(r.total_umis),
        Rcpp::Named("concordant") = static_cast<double>(r.concordant),
        Rcpp::Named("discordant") = static_cast<double>(r.discordant),
        Rcpp::Named("nSnvs") = r.n_snvs,
        Rcpp::Named("nInformative") = r.n_informative,
        Rcpp::Named("nNonInformative") = r.n_noninformative,
        Rcpp::Named("nMonoallelic") = r.n_monoallelic,
        Rcpp::Named("medianMaf") = r.median_maf,
        Rcpp::Named("cellCounts") = counts,
        Rcpp::Named("whitelistUsed") = r.whitelist_used,
        Rcpp::Named("whitelistLoaded") = static_cast<double>(r.wl_loaded),
        Rcpp::Named("whitelistMatched") = static_cast<double>(r.wl_matched),
        Rcpp::Named("whitelistReadsChecked") = static_cast<double>(r.wl_reads_checked),
        Rcpp::Named("whitelistReadsKept") = static_cast<double>(r.wl_reads_kept),
        Rcpp::Named("whitelistReadsRejected") = static_cast<double>(r.wl_reads_rejected));
    Rcpp::List profile = Rcpp::List::create(
        Rcpp::Named("cell") = r.prof_cell, Rcpp::Named("snv") = r.prof_snv,
        Rcpp::Named("refUmis") = r.prof_ref, Rcpp::Named("altUmis") = r.prof_alt);
    return Rcpp::List::create(Rcpp::Named("cells") = cells, Rcpp::Named("snvs") = snvs,
                              Rcpp::Named("summary") = summary, Rcpp::Named("profile") = profile);
}
