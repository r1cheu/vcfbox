#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <stdexcept>
#include <vector>

#include "parentage/sites.h"

using Catch::Matchers::ContainsSubstring;
using vcfbox::parentage::ContigInfo;
using vcfbox::parentage::decode_base;
using vcfbox::parentage::encode_base;
using vcfbox::parentage::match_contigs;
using vcfbox::parentage::no_base;
using vcfbox::parentage::SiteTable;

TEST_CASE(
    "encode_base maps ACGT in either case and rejects the rest",
    "[parentage][sites]")
{
    REQUIRE(encode_base('A') == 0);
    REQUIRE(encode_base('c') == 1);
    REQUIRE(encode_base('G') == 2);
    REQUIRE(encode_base('t') == 3);
    REQUIRE(encode_base('N') == no_base);
    REQUIRE(encode_base('-') == no_base);
    REQUIRE(decode_base(2) == 'G');
    REQUIRE(decode_base(no_base) == 'N');
}

TEST_CASE(
    "SiteTable keeps sites grouped by contig with 0-based positions",
    "[parentage][sites]")
{
    SiteTable table;
    table.add_contig({"chr01", 1000});
    table.add_site(100, encode_base('A'), encode_base('G'), 0.25F);
    table.add_site(204, encode_base('C'), encode_base('T'), 0.5F);
    table.add_contig({"chr02", 500});
    table.add_site(6, encode_base('G'), encode_base('A'), 0.0F);
    table.finish();

    REQUIRE(table.size() == 3);
    REQUIRE(table.contigs().size() == 2);
    REQUIRE(table.contigs()[0].info.name == "chr01");
    REQUIRE(table.contigs()[0].info.length == 1000);
    REQUIRE(table.contigs()[0].first_site == 0);
    REQUIRE(table.contigs()[0].n_sites == 2);
    REQUIRE(table.contigs()[1].first_site == 2);
    REQUIRE(table.contigs()[1].n_sites == 1);
    REQUIRE(table.pos(2) == 6);
    REQUIRE(table.allele0(1) == encode_base('C'));
    REQUIRE(table.allele1(1) == encode_base('T'));
    REQUIRE(table.freq1(0) == 0.25F);
    REQUIRE(table.contig_positions(0).size() == 2);
    REQUIRE(table.contig_positions(0)[1] == 204);
    REQUIRE(table.describe(0) == "chr01:101 A>G");
}

TEST_CASE(
    "SiteTable drops contigs that received no sites",
    "[parentage][sites]")
{
    SiteTable table;
    table.add_contig({"empty_head", 10});
    table.add_contig({"chr01", 1000});
    table.add_site(0, 0, 1, 0.0F);
    table.add_contig({"empty_tail", 10});
    table.finish();

    REQUIRE(table.contigs().size() == 1);
    REQUIRE(table.contigs()[0].info.name == "chr01");
}

TEST_CASE(
    "SiteTable rejects input that breaks its ordering invariant",
    "[parentage][sites]")
{
    SiteTable table;

    SECTION("site before any contig")
    {
        REQUIRE_THROWS_AS(table.add_site(0, 0, 1, 0.0F), std::logic_error);
    }
    SECTION("contig name used twice")
    {
        table.add_contig({"chr01", 1000});
        table.add_site(0, 0, 1, 0.0F);
        table.add_contig({"chr02", 500});
        table.add_site(0, 0, 1, 0.0F);
        REQUIRE_THROWS_WITH(
            table.add_contig({"chr01", 1000}),
            ContainsSubstring("Contig chr01 appears twice"));
    }
    SECTION("position equal to the previous one")
    {
        table.add_contig({"chr01", 1000});
        table.add_site(100, 0, 1, 0.0F);
        REQUIRE_THROWS_WITH(
            table.add_site(100, 0, 1, 0.0F),
            ContainsSubstring(
                "chr01:101: panel positions must be strictly increasing"));
    }
    SECTION("position below the previous one")
    {
        table.add_contig({"chr01", 1000});
        table.add_site(100, 0, 1, 0.0F);
        REQUIRE_THROWS_WITH(
            table.add_site(99, 0, 1, 0.0F),
            ContainsSubstring(
                "chr01:100: panel positions must be strictly increasing"));
    }
    SECTION("position beyond int32")
    {
        table.add_contig({"chr01", 1000});
        REQUIRE_THROWS_WITH(
            table.add_site(int64_t{1} << 31U, 0, 1, 0.0F),
            ContainsSubstring("exceeds the supported range"));
    }
    SECTION("finish with no sites")
    {
        table.add_contig({"chr01", 1000});
        REQUIRE_THROWS_WITH(
            table.finish(),
            ContainsSubstring("Panel contains no usable sites"));
    }
}

TEST_CASE(
    "match_contigs maps panel contigs to BAM tids by name",
    "[parentage][sites]")
{
    SiteTable table;
    table.add_contig({"chr01", 1000});
    table.add_site(0, 0, 1, 0.0F);
    table.add_contig({"chr02", 500});
    table.add_site(0, 0, 1, 0.0F);
    table.finish();

    SECTION("same names in a different order")
    {
        const std::vector<ContigInfo> bam{
            {"chr02", 500}, {"chrM", 16}, {"chr01", 1000}};
        REQUIRE(match_contigs(table, bam, "x.bam") == std::vector<int>{2, 0});
    }
    SECTION("contig absent from the BAM")
    {
        const std::vector<ContigInfo> bam{{"chr01", 1000}};
        REQUIRE_THROWS_WITH(
            match_contigs(table, bam, "x.bam"),
            ContainsSubstring("Panel contig chr02 is not present in x.bam"));
    }
    SECTION("contig length differs")
    {
        const std::vector<ContigInfo> bam{{"chr01", 1000}, {"chr02", 501}};
        REQUIRE_THROWS_WITH(
            match_contigs(table, bam, "x.bam"),
            ContainsSubstring(
                "Contig chr02 has length 500 in the panel but 501 in x.bam"));
    }
}
