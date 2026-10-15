#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <format>
#include <string>
#include <vector>

#include "file_fixture.h"
#include "parentage/counts.h"
#include "parentage/counts_io.h"
#include "parentage/parentage.h"
#include "sam_fixture.h"
#include "vcf_fixture.h"

using Catch::Matchers::ContainsSubstring;
using vcfbox::parentage::read_counts;
using vcfbox::parentage::run;
using vcfbox::parentage::RunOptions;
using vcfbox::parentage::SampleCounts;
using vcfbox::test::FileFixture;
using vcfbox::test::PanelVcf;
using vcfbox::test::repeated_sites;
using vcfbox::test::SamFile;
using vcfbox::test::SamRead;

namespace
{
// One read of `n` bases from chr01:1 carrying allele0 (A) everywhere.
auto allele0_sam(int n) -> std::string
{
    SamFile sam;
    sam.reads
        = {{.cigar = std::format("{}M", n),
            .seq = std::string(static_cast<size_t>(n), 'A'),
            .qual = std::string(static_cast<size_t>(n), 'I')}};
    return sam.text();
}

auto n0_of(const SampleCounts& counts) -> std::vector<uint8_t>
{
    return {counts.n0().begin(), counts.n0().end()};
}
}  // namespace

TEST_CASE(
    "run processes every sample on the worker threads",
    "[parentage][run]")
{
    FileFixture files;
    PanelVcf panel;
    panel.sites = repeated_sites("chr01", 3, {"0/0", "1/1", "0/0"});
    RunOptions options;
    options.vcf_path = files.create_text_file(panel.text(), ".vcf");
    options.counts_dir = files.missing_path("counts");
    options.threads = 2;
    const std::filesystem::path sam_a
        = files.create_text_file(allele0_sam(3), ".sam");
    const std::filesystem::path sam_b
        = files.create_text_file(allele0_sam(2), ".sam");

    SECTION("writes the counts file of every sample")
    {
        options.bams_path = files.create_text_file(
            std::format("A\t{}\nB\t{}\n", sam_a.string(), sam_b.string()));

        run(options);

        const SampleCounts a = read_counts(
            (std::filesystem::path(options.counts_dir) / "A.counts.gz")
                .string(),
            3);
        const SampleCounts b = read_counts(
            (std::filesystem::path(options.counts_dir) / "B.counts.gz")
                .string(),
            3);
        REQUIRE(n0_of(a) == std::vector<uint8_t>{1, 1, 1});
        REQUIRE(n0_of(b) == std::vector<uint8_t>{1, 1, 0});
    }
    SECTION("a failing sample fails the run")
    {
        options.bams_path = files.create_text_file(
            std::format(
                "A\t{}\nB\t{}\n",
                sam_a.string(),
                files.missing_path(".sam").string()));

        REQUIRE_THROWS_WITH(run(options), ContainsSubstring(".sam"));
    }
    SECTION("zero threads is rejected")
    {
        options.bams_path
            = files.create_text_file(std::format("A\t{}\n", sam_a.string()));
        options.threads = 0;

        REQUIRE_THROWS_WITH(run(options), ContainsSubstring("--threads"));
    }
}
