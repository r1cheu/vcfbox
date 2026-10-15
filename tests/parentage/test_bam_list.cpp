#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <vector>

#include "file_fixture.h"
#include "parentage/parentage.h"

using Catch::Matchers::ContainsSubstring;
using vcfbox::parentage::BamEntry;
using vcfbox::parentage::read_bam_list;
using vcfbox::test::FileFixture;

TEST_CASE(
    "read_bam_list takes a sample name or derives it from the path",
    "[parentage][bam_list]")
{
    FileFixture files;

    SECTION("two columns")
    {
        const std::vector<BamEntry> entries
            = read_bam_list(files.create_text_file("SH1\t/x/a.bam\n"));
        REQUIRE(entries.size() == 1);
        REQUIRE(entries[0].sample == "SH1");
        REQUIRE(entries[0].path == "/x/a.bam");
    }
    SECTION("one column strips the .bam extension")
    {
        const std::vector<BamEntry> entries
            = read_bam_list(files.create_text_file("/x/a.sorted.bam\n"));
        REQUIRE(entries[0].sample == "a.sorted");
        REQUIRE(entries[0].path == "/x/a.sorted.bam");
    }
    SECTION("one column without .bam keeps the file name")
    {
        const std::vector<BamEntry> entries
            = read_bam_list(files.create_text_file("/x/a.cram\n"));
        REQUIRE(entries[0].sample == "a.cram");
    }
    SECTION("blank lines and comments are skipped")
    {
        const std::vector<BamEntry> entries = read_bam_list(
            files.create_text_file("# list\n\nA\t/a.bam\n\nB\t/b.bam\n"));
        REQUIRE(entries.size() == 2);
        REQUIRE(entries[1].sample == "B");
    }
}

TEST_CASE("read_bam_list rejects malformed lists", "[parentage][bam_list]")
{
    FileFixture files;

    SECTION("three columns")
    {
        REQUIRE_THROWS_WITH(
            read_bam_list(files.create_text_file("A\t/a.bam\textra\n")),
            ContainsSubstring(
                ":1: expected `sample<TAB>path` or `path`, found 3 fields"));
    }
    SECTION("sample listed twice")
    {
        REQUIRE_THROWS_WITH(
            read_bam_list(files.create_text_file("A\t/a.bam\nA\t/b.bam\n")),
            ContainsSubstring(":2: sample A is listed twice"));
    }
}
