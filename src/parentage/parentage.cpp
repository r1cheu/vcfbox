#include "parentage/parentage.h"

#include <Eigen/Core>
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "common/line_writer.h"
#include "common/lines.h"
#include "common/progress.h"
#include "common/strings.h"
#include "parentage/bitmatrix.h"
#include "parentage/likelihood.h"
#include "parentage/pileup.h"
#include "parentage/report.h"
#include "parentage/sites.h"

namespace vcfbox
{
namespace
{
struct LoadedMatrices
{
    std::vector<Site> sites;
    Eigen::MatrixXd m0;
    Eigen::MatrixXd m1;
    Eigen::MatrixXd p0;
    Eigen::MatrixXd p1;
    std::vector<std::string> maternal_names;
    std::vector<std::string> paternal_names;
};

LoadedMatrices load_matrices(const std::string& prefix)
{
    LoadedMatrices loaded;
    loaded.sites = load_sites(prefix + ".sites.bed");
    loaded.m0 = load_bitmatrix(prefix + ".M0.bin");
    loaded.m1 = load_bitmatrix(prefix + ".M1.bin");
    loaded.p0 = load_bitmatrix(prefix + ".P0.bin");
    loaded.p1 = load_bitmatrix(prefix + ".P1.bin");
    loaded.maternal_names = read_lines(prefix + ".maternal.tsv");
    loaded.paternal_names = read_lines(prefix + ".paternal.tsv");

    const auto n = static_cast<Eigen::Index>(loaded.sites.size());
    if (loaded.m0.rows() != n || loaded.m1.rows() != n || loaded.p0.rows() != n
        || loaded.p1.rows() != n)
    {
        throw std::runtime_error(
            "Matrix row count does not match sites count in prefix: " + prefix);
    }
    if (loaded.m0.cols()
            != static_cast<Eigen::Index>(loaded.maternal_names.size())
        || loaded.p0.cols()
               != static_cast<Eigen::Index>(loaded.paternal_names.size()))
    {
        throw std::runtime_error(
            "Matrix column count does not match parent name count in prefix: "
            + prefix);
    }
    return loaded;
}

void write_top_k(
    LineWriter& out,
    const std::string& sample,
    const SampleReport& report,
    std::span<const std::string> maternal_names,
    std::span<const std::string> paternal_names)
{
    const double best_ll = report.best.loglik;
    for (size_t i = 0; i < report.top_k.size(); ++i)
    {
        const auto& cand = report.top_k[i];
        out.write_fields(
            sample,
            i + 1,
            maternal_names[cand.m],
            paternal_names[cand.f],
            format_number(cand.loglik, 8),
            format_number(cand.posterior, 6),
            format_number(cand.loglik - best_ll, 8));
    }
}

void write_summary_row(
    LineWriter& out,
    const std::string& sample,
    const SampleReport& report,
    std::span<const std::string> maternal_names,
    std::span<const std::string> paternal_names,
    double threshold)
{
    out.write_fields(
        sample,
        maternal_names[report.best.m],
        paternal_names[report.best.f],
        format_number(report.best.loglik, 8),
        format_number(report.best.posterior, 6),
        format_number(report.m_marg_post, 6),
        format_number(report.f_marg_post, 6),
        to_string(classify_call(report, threshold)));
}

void write_raw_block(
    LineWriter& out,
    const std::string& sample,
    const Eigen::MatrixXd& ll,
    std::span<const std::string> maternal_names,
    std::span<const std::string> paternal_names)
{
    for (Eigen::Index r = 0; r < ll.rows(); ++r)
    {
        for (Eigen::Index c = 0; c < ll.cols(); ++c)
        {
            out.write_fields(
                sample,
                maternal_names[r],
                paternal_names[c],
                format_number(ll(r, c), 8));
        }
    }
}

void emit_reports(
    const ParentageTestOptions& options,
    std::span<const std::string> bam_paths,
    const LoadedMatrices& loaded,
    std::span<const SampleReport> reports,
    std::span<const Eigen::MatrixXd> ll_per_bam)
{
    const bool emit_raw = !options.raw_path.empty();
    const std::span<const std::string> mat_names = loaded.maternal_names;
    const std::span<const std::string> pat_names = loaded.paternal_names;

    LineWriter top(options.output_path);
    top.write_fields(
        "sample", "rank", "mother", "father",
        "loglik", "posterior", "delta_best");
    LineWriter summary(options.summary_path);
    summary.write_fields(
        "sample", "mother", "father", "loglik",
        "pair_post", "m_marg_post", "f_marg_post", "call");
    std::optional<LineWriter> raw;
    if (emit_raw)
    {
        raw.emplace(options.raw_path);
        raw->write_fields("sample", "mother", "father", "loglik");
    }

    for (size_t i = 0; i < bam_paths.size(); ++i)
    {
        const auto sample
            = std::filesystem::path(bam_paths[i]).stem().string();
        write_top_k(top, sample, reports[i], mat_names, pat_names);
        write_summary_row(
            summary,
            sample,
            reports[i],
            mat_names,
            pat_names,
            options.threshold);
        if (emit_raw)
        {
            write_raw_block(
                *raw, sample, ll_per_bam[i], mat_names, pat_names);
        }
    }
}
}  // namespace

void test_parentage(const ParentageTestOptions& options)
{
    const auto loaded = load_matrices(options.prefix);
    const auto bam_paths = read_lines(options.bam_list_path);

    const PileupOptions pileup_opts{options.min_mapq, options.min_baseq};
    const int n_threads = std::max(1, options.threads);
    const bool emit_raw = !options.raw_path.empty();

    std::vector<Eigen::MatrixXd> ll_per_bam;
    if (emit_raw)
    {
        ll_per_bam.resize(bam_paths.size());
    }
    std::vector<SampleReport> reports(bam_paths.size());
    std::atomic<size_t> next_bam{0};
    size_t processed = 0;
    auto counter = create_counter("Scoring BAMs", processed, "bam/s");
    counter->show();

    auto worker = [&] {
        while (true)
        {
            const size_t i = next_bam.fetch_add(1);
            if (i >= bam_paths.size())
            {
                break;
            }
            const auto counts
                = count_alleles(bam_paths[i], loaded.sites, pileup_opts);
            auto ll = compute_pair_loglik(
                counts,
                loaded.m0,
                loaded.m1,
                loaded.p0,
                loaded.p1,
                options.error_rate);
            reports[i] = build_report(ll, options.top_k);
            if (emit_raw)
            {
                ll_per_bam[i] = std::move(ll);
            }
            ++processed;
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(n_threads - 1);
    for (int t = 0; t < n_threads - 1; ++t)
    {
        pool.emplace_back(worker);
    }
    worker();
    for (auto& t : pool)
    {
        t.join();
    }
    counter->done();

    emit_reports(options, bam_paths, loaded, reports, ll_per_bam);
}
}  // namespace vcfbox
