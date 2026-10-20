#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
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

auto file_text(const std::string& path) -> std::string
{
    std::ifstream in(path);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
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
    options.output_prefix = files.missing_path("").string();
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
    SECTION("no output prefix is rejected")
    {
        options.bams_path
            = files.create_text_file(std::format("A\t{}\n", sam_a.string()));
        options.output_prefix.clear();

        REQUIRE_THROWS_WITH(run(options), ContainsSubstring("--output"));
    }
    SECTION("zero threads is rejected")
    {
        options.bams_path
            = files.create_text_file(std::format("A\t{}\n", sam_a.string()));
        options.threads = 0;

        REQUIRE_THROWS_WITH(run(options), ContainsSubstring("--threads"));
    }
}

TEST_CASE(
    "run writes a report row per sample with the expected pair",
    "[parentage][run]")
{
    // P1 and P3 are identical, so every sample is best explained by P1
    // selfing, with (P1, P3) tied second: AMBIGUOUS, with the pair leading.
    FileFixture files;
    PanelVcf panel;
    panel.sites = repeated_sites("chr01", 3, {"0/0", "1/1", "0/0"});
    RunOptions options;
    options.vcf_path = files.create_text_file(panel.text(), ".vcf");
    options.bams_path = files.create_text_file(
        std::format(
            "A\t{}\nB\t{}\n",
            files.create_text_file(allele0_sam(3), ".sam").string(),
            files.create_text_file(allele0_sam(3), ".sam").string()));
    options.expected_path = files.create_text_file(
        "sample\tmother\tfather\nA\tP3\tP1\nB\tP1\tQ9\n");
    options.output_prefix = files.missing_path("").string();

    run(options);

    const std::string ambiguous
        = file_text(options.output_prefix + ".ambiguous.tsv");
    const std::string summary
        = file_text(options.output_prefix + ".summary.txt");
    CAPTURE(ambiguous, summary);
    SECTION("writes an AMBIGUOUS sample with the contest it fails")
    {
        REQUIRE_THAT(
            ambiguous,
            ContainsSubstring("\nA\t1.000\tmodel\tP1/P1\tP1\t0.5\tPASS\t"));
    }
    SECTION("compares a present expected pair with the best pair")
    {
        REQUIRE_THAT(ambiguous, ContainsSubstring("\tP3\tP1\tMISMATCH\t0.0\n"));
    }
    SECTION("marks an expected parent missing from the panel ABSENT")
    {
        REQUIRE_THAT(ambiguous, ContainsSubstring("\tP1\tQ9\tABSENT\tNA\n"));
    }
    SECTION("writes only the header of a status no sample has")
    {
        REQUIRE(
            file_text(options.output_prefix + ".two_parent.tsv")
            == "sample\tdepth_mean\ta\tb\tis_selfing\thom_mismatch_rate\t"
               "het_alt_fraction\tfit\texpected_a\texpected_b\t"
               "expected_status\texpected_lod_vs_best\n");
    }
    SECTION("lists the top pairs of every sample")
    {
        REQUIRE_THAT(
            file_text(options.output_prefix + ".topk.tsv"),
            ContainsSubstring("\nB\t3\tP3\tP3\t3\t0\t0\t0\t"));
    }
    SECTION("counts statuses, TWO_PARENT fits and expected statuses")
    {
        REQUIRE_THAT(summary, ContainsSubstring("samples\t2\n"));
        REQUIRE_THAT(summary, ContainsSubstring("\nAMBIGUOUS\t2\n"));
        REQUIRE_THAT(summary, ContainsSubstring("\nPASS\t0\n"));
        REQUIRE_THAT(
            summary, ContainsSubstring("\nAMBIGUOUS\t0\t1\t1\t0\n"));
    }
}
