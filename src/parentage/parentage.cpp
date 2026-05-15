#include "parentage/parentage.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "common/line_writer.h"
#include "common/lines.h"
#include "common/progress.h"
#include "parentage/bitmatrix.h"
#include "parentage/likelihood.h"
#include "parentage/pileup.h"
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

std::string format_loglik(double ll)
{
    std::array<char, 32> buf{};
    const auto [ptr, ec] = std::to_chars(
        buf.data(), buf.data() + buf.size(), ll, std::chars_format::general, 8);
    if (ec != std::errc{})
    {
        throw std::runtime_error("Failed to format loglik value");
    }
    return std::string(buf.data(), ptr);
}
}  // namespace

void test_parentage(const ParentageTestOptions& options)
{
    const auto loaded = load_matrices(options.prefix);
    const auto bam_paths = read_lines(options.bam_list_path);

    const PileupOptions pileup_opts{options.min_mapq, options.min_baseq};
    const int n_threads = std::max(1, options.threads);

    std::vector<Eigen::MatrixXd> ll_per_bam(bam_paths.size());
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
            ll_per_bam[i] = compute_pair_loglik(
                counts,
                loaded.m0,
                loaded.m1,
                loaded.p0,
                loaded.p1,
                options.error_rate);
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

    LineWriter out(options.output_path);
    out.write_line("sample\tmaternal\tpaternal\tloglik");
    for (size_t i = 0; i < bam_paths.size(); ++i)
    {
        const auto sample
            = std::filesystem::path(bam_paths[i]).stem().string();
        const auto& ll = ll_per_bam[i];
        for (Eigen::Index r = 0; r < ll.rows(); ++r)
        {
            for (Eigen::Index c = 0; c < ll.cols(); ++c)
            {
                out.write_line(
                    sample + "\t" + loaded.maternal_names[r] + "\t"
                    + loaded.paternal_names[c] + "\t"
                    + format_loglik(ll(r, c)));
            }
        }
    }
}
}  // namespace vcfbox
