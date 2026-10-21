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
using vcfbox::test::repeated_reads;
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

// One allele0 read of 3 bases from chr01:1 followed by `n_unmapped` unmapped
// reads.
auto partly_unmapped_sam(int n_unmapped) -> std::string
{
    SamFile sam;
    sam.reads = repeated_reads(
        {.qname = "u", .flag = 0x4, .cigar = "*", .seq = "A", .qual = "I"},
        n_unmapped);
    sam.reads.insert(
        sam.reads.begin(), {.cigar = "3M", .seq = "AAA", .qual = "III"});
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
    SECTION("writes the panel and every sample summary to the log")
    {
        options.bams_path = files.create_text_file(
            std::format("A\t{}\nB\t{}\n", sam_a.string(), sam_b.string()));

        run(options);

        const std::string log
            = file_text(std::format("{}.log", options.output_prefix));
        REQUIRE_THAT(log, ContainsSubstring("Panel: 3 parents, 3 sites"));
        REQUIRE_THAT(
            log,
            ContainsSubstring("Sample A (pileup): 1 of 1 reads used"));
        REQUIRE_THAT(
            log,
            ContainsSubstring("Sample B (pileup): 1 of 1 reads used"));
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
    // selfing, with (P1, P3) tied second: ambiguous, with the pair leading.
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

    const std::string table = file_text(options.output_prefix + ".tsv");
    const std::string summary
        = file_text(options.output_prefix + ".summary.txt");
    CAPTURE(table, summary);
    SECTION("writes the header")
    {
        REQUIRE(table.starts_with(
            "sample\tcoverage\tbest\tbest_a\tbest_b\texpected_a\t"
            "expected_b\tnull\tsingle\tbest_pair\tsecond_pair\texpected\t"
            "pair_hom_mismatch_rate\tpair_het_alt_fraction\tpair_fit\t"
            "single_allele_absent_rate\tcontest\texpected_pair_rank\n"));
    }
    SECTION("writes the call, the leading parents and the expected parents")
    {
        REQUIRE_THAT(
            table, ContainsSubstring("\nA\t1.000\tambiguous\tP1\tP1\tP3\tP1\t"));
    }
    SECTION("writes the log-likelihoods of every model and the expected pair")
    {
        REQUIRE_THAT(
            table, ContainsSubstring("\tP3\tP1\t-1.2\t-0.6\t-0.1\t-0.1\t-0.1\t"));
    }
    SECTION("writes the diagnostics, the contest and the expected rank")
    {
        REQUIRE_THAT(
            table, ContainsSubstring("\t0.00000\tNA\tPASS\tNA\tmodel\t2\n"));
    }
    SECTION("an expected parent missing from the panel is ABSENT")
    {
        REQUIRE_THAT(
            table,
            ContainsSubstring("\tP1\tQ9\t-1.2\t-0.6\t-0.1\t-0.1\tNA\t"));
        REQUIRE_THAT(table, ContainsSubstring("\tmodel\tABSENT\n"));
    }
    SECTION("lists the top pairs of every sample")
    {
        REQUIRE_THAT(
            file_text(options.output_prefix + ".topk.tsv"),
            ContainsSubstring("\nB\t3\tP3\tP3\t3\t0\t0\t0\t"));
    }
    SECTION("counts every best, the fits of pair calls and the expected ranks")
    {
        REQUIRE_THAT(summary, ContainsSubstring("samples\t2\n"));
        REQUIRE_THAT(summary, ContainsSubstring("\nambiguous\t2\n"));
        REQUIRE_THAT(summary, ContainsSubstring("\nPASS\t0\n"));
        REQUIRE_THAT(
            summary,
            ContainsSubstring("best\\expected_pair_rank\t1\t2+\tABSENT\tNA\n"));
        REQUIRE_THAT(summary, ContainsSubstring("\nambiguous\t0\t1\t1\t0\n"));
    }
}

TEST_CASE(
    "run drops a sample with more than max_unmapped_rate of reads unmapped",
    "[parentage][run]")
{
    // A has 1 of 2 reads unmapped, at the rate; B has 2 of 3, above it.
    FileFixture files;
    PanelVcf panel;
    panel.sites = repeated_sites("chr01", 3, {"0/0", "1/1", "0/0"});
    RunOptions options;
    options.vcf_path = files.create_text_file(panel.text(), ".vcf");
    options.bams_path = files.create_text_file(
        std::format(
            "A\t{}\nB\t{}\n",
            files.create_text_file(partly_unmapped_sam(1), ".sam").string(),
            files.create_text_file(partly_unmapped_sam(2), ".sam").string()));
    options.output_prefix = files.missing_path("").string();
    options.max_unmapped_rate = 0.5;

    run(options);

    SECTION("scores the sample at the rate and leaves the other out")
    {
        const std::string topk
            = file_text(options.output_prefix + ".topk.tsv");
        CAPTURE(topk);
        REQUIRE_THAT(topk, ContainsSubstring("\nA\t"));
        REQUIRE_THAT(topk, !ContainsSubstring("\nB\t"));
    }
    SECTION("writes the dropped sample as a row without scores")
    {
        REQUIRE_THAT(
            file_text(options.output_prefix + ".tsv"),
            ContainsSubstring(
                "\nB\t1.000\tdropped\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"
                "NA\tNA\tNA\tNA\tNA\tNA\n"));
    }
    SECTION("logs the dropped sample with its unmapped reads")
    {
        REQUIRE_THAT(
            file_text(std::format("{}.log", options.output_prefix)),
            ContainsSubstring("Sample B dropped: 2 of 3 reads unmapped"));
    }
    SECTION("counts the dropped sample in the summary")
    {
        const std::string summary
            = file_text(options.output_prefix + ".summary.txt");
        REQUIRE_THAT(summary, ContainsSubstring("samples\t2\n"));
        REQUIRE_THAT(summary, ContainsSubstring("\ndropped\t1\n"));
    }
}
