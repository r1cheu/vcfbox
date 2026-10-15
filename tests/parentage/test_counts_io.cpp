#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <string>
#include <vector>

#include "file_fixture.h"
#include "parentage/counts.h"
#include "parentage/counts_io.h"

using Catch::Matchers::ContainsSubstring;
using vcfbox::parentage::PileupStats;
using vcfbox::parentage::read_counts;
using vcfbox::parentage::SampleCounts;
using vcfbox::parentage::write_counts;
using vcfbox::test::FileFixture;

TEST_CASE("counts file round-trips every field", "[parentage][counts_io]")
{
    FileFixture files;
    const std::string path = files.missing_path(".counts.gz");
    SampleCounts counts("SH1", 5);
    counts.add0(0);
    counts.add1(1);
    counts.add1(1);
    counts.add_other(4);
    counts.stats().reads_total = 10;
    counts.stats().reads_used = 7;
    counts.stats().reads_skipped_flag = 1;
    counts.stats().reads_skipped_improper = 6;
    counts.stats().reads_skipped_mapq = 2;
    counts.stats().reads_off_panel = 3;
    counts.stats().bases_skipped_baseq = 4;

    write_counts(counts, path);
    const SampleCounts back = read_counts(path, 5);

    REQUIRE(back.sample() == "SH1");
    REQUIRE(back.n_sites() == 5);
    REQUIRE(
        std::vector<uint8_t>(back.n0().begin(), back.n0().end())
        == std::vector<uint8_t>{1, 0, 0, 0, 0});
    REQUIRE(
        std::vector<uint8_t>(back.n1().begin(), back.n1().end())
        == std::vector<uint8_t>{0, 2, 0, 0, 0});
    REQUIRE(
        std::vector<uint8_t>(back.n_other().begin(), back.n_other().end())
        == std::vector<uint8_t>{0, 0, 0, 0, 1});
    REQUIRE(back.stats().reads_total == 10);
    REQUIRE(back.stats().reads_used == 7);
    REQUIRE(back.stats().reads_skipped_flag == 1);
    REQUIRE(back.stats().reads_skipped_improper == 6);
    REQUIRE(back.stats().reads_skipped_mapq == 2);
    REQUIRE(back.stats().reads_off_panel == 3);
    REQUIRE(back.stats().bases_skipped_baseq == 4);
}

TEST_CASE(
    "read_counts rejects files that do not fit the panel",
    "[parentage][counts_io]")
{
    FileFixture files;

    SECTION("site count differs")
    {
        const std::string path = files.missing_path(".counts.gz");
        write_counts(SampleCounts("SH1", 5), path);
        REQUIRE_THROWS_WITH(
            read_counts(path, 4),
            ContainsSubstring("holds 5 sites but the panel has 4"));
    }
    SECTION("not a counts file")
    {
        REQUIRE_THROWS_WITH(
            read_counts(files.create_text_file("hello world"), 5),
            ContainsSubstring("is not a vcfbox counts file"));
    }
    SECTION("truncated after the magic")
    {
        REQUIRE_THROWS_WITH(
            read_counts(files.create_text_file("VBXCNT01"), 5),
            ContainsSubstring("is truncated"));
    }
    SECTION("file does not exist")
    {
        REQUIRE_THROWS_WITH(
            read_counts(files.missing_path(".counts.gz"), 5),
            ContainsSubstring("Failed to open counts file"));
    }
}

TEST_CASE(
    "SampleCounts rejects arrays of different sizes",
    "[parentage][counts_io]")
{
    REQUIRE_THROWS_WITH(
        SampleCounts("S", {0, 0}, {0}, {0, 0}, PileupStats{}),
        ContainsSubstring("count arrays differ in size (2, 1, 2)"));
}
