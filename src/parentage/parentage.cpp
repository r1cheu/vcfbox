#include "parentage/parentage.h"

#include <Eigen/Core>
#include <array>
#include <charconv>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
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
    const auto loaded = load_matrices(options.matrix_prefix);
    const auto bam_paths = read_lines(options.bam_list_path);

    const PileupOptions pileup_opts{options.min_mapq, options.min_baseq};

    LineWriter out(options.output_path);
    out.write_line("sample\tmaternal\tpaternal\tloglik");

    size_t processed = 0;
    auto counter = create_counter("Scoring BAMs", processed, "bam/s");
    counter->show();

    for (const auto& bam_path : bam_paths)
    {
        const auto sample = std::filesystem::path(bam_path).stem().string();
        const auto counts = count_alleles(bam_path, loaded.sites, pileup_opts);
        const auto ll = compute_pair_loglik(
            counts,
            loaded.m0,
            loaded.m1,
            loaded.p0,
            loaded.p1,
            options.error_rate);

        for (Eigen::Index i = 0; i < ll.rows(); ++i)
        {
            for (Eigen::Index j = 0; j < ll.cols(); ++j)
            {
                out.write_line(
                    sample + "\t" + loaded.maternal_names[i] + "\t"
                    + loaded.paternal_names[j] + "\t"
                    + format_loglik(ll(i, j)));
            }
        }
        ++processed;
    }
    counter->done();
}
}  // namespace vcfbox
