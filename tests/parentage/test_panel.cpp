#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include "file_fixture.h"
#include "parentage/bitrows.h"
#include "parentage/panel.h"
#include "parentage/sites.h"
#include "vcf_fixture.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using vcfbox::parentage::BitRows;
using vcfbox::parentage::encode_base;
using vcfbox::parentage::load_panel;
using vcfbox::parentage::Panel;
using vcfbox::parentage::PanelOptions;
using vcfbox::parentage::SiteTable;
using vcfbox::test::FileFixture;
using vcfbox::test::PanelVcf;
using vcfbox::test::repeated_sites;

namespace
{
auto load(FileFixture& files, const PanelVcf& vcf, PanelOptions options = {})
    -> Panel
{
    return load_panel(files.create_text_file(vcf.text(), ".vcf"), options);
}
}  // namespace

TEST_CASE(
    "load_panel fills sites, names and genotype bits from a VCF",
    "[parentage][panel]")
{
    FileFixture files;
    PanelVcf vcf;
    vcf.sites = {
        {"chr01", 101, "A", "G", {"1/1", "0/0", "0/0"}},
        {"chr01", 205, "C", "T", {"0/0", "1/1", "1/1"}},
        {"chr02", 7, "G", "A", {"0/0", "0/0", "1/1"}},
    };

    const Panel panel = load(files, vcf);

    REQUIRE(panel.n_parents() == 3);
    REQUIRE(panel.names()[1] == "P2");
    REQUIRE(panel.n_sites() == 3);
    REQUIRE(panel.n_excluded_sites() == 0);
    REQUIRE(panel.near_duplicates().empty());

    const SiteTable& sites = panel.sites();
    REQUIRE(sites.contigs().size() == 2);
    REQUIRE(sites.contigs()[0].info.name == "chr01");
    REQUIRE(sites.contigs()[0].info.length == 1000);
    REQUIRE(sites.contigs()[0].n_sites == 2);
    REQUIRE(sites.contigs()[1].info.name == "chr02");
    REQUIRE(sites.contigs()[1].first_site == 2);
    REQUIRE(sites.pos(0) == 100);
    REQUIRE(sites.pos(2) == 6);
    REQUIRE(sites.allele0(1) == encode_base('C'));
    REQUIRE(sites.allele1(1) == encode_base('T'));
    REQUIRE_THAT(sites.freq1(0), WithinAbs(1.0 / 3, 1e-6));
    REQUIRE_THAT(sites.freq1(1), WithinAbs(2.0 / 3, 1e-6));
    REQUIRE_THAT(sites.freq1(2), WithinAbs(1.0 / 3, 1e-6));

    const BitRows& geno = panel.geno();
    REQUIRE(geno.rows() == 3);
    REQUIRE(geno.n_bits() == 3);
    REQUIRE(geno.get(0, 0));
    REQUIRE_FALSE(geno.get(0, 1));
    REQUIRE_FALSE(geno.get(1, 0));
    REQUIRE(geno.get(1, 1));
    REQUIRE(geno.get(2, 1));
    REQUIRE(geno.get(2, 2));
    REQUIRE_FALSE(geno.get(1, 2));
}

TEST_CASE(
    "load_panel places bits correctly across 64-site word boundaries",
    "[parentage][panel]")
{
    FileFixture files;
    PanelVcf vcf;

    SECTION("exactly one word")
    {
        vcf.sites = repeated_sites("chr01", 64, {"1/1", "0/0", "0/0"});
        const Panel panel = load(files, vcf);
        REQUIRE(panel.n_sites() == 64);
        REQUIRE(panel.geno().words_per_row() == 1);
        REQUIRE(panel.geno().get(0, 63));
        REQUIRE_FALSE(panel.geno().get(1, 63));
    }
    SECTION("one site into the second word")
    {
        vcf.sites = repeated_sites("chr01", 65, {"1/1", "0/0", "0/0"});
        const Panel panel = load(files, vcf);
        REQUIRE(panel.geno().words_per_row() == 2);
        REQUIRE(panel.geno().get(0, 63));
        REQUIRE(panel.geno().get(0, 64));
        REQUIRE_FALSE(panel.geno().get(1, 64));
        REQUIRE(panel.geno().row(0)[1] == 1);
    }
    SECTION("three words")
    {
        vcf.sites = repeated_sites("chr01", 130, {"0/0", "1/1", "0/0"});
        const Panel panel = load(files, vcf);
        REQUIRE(panel.geno().words_per_row() == 3);
        REQUIRE(panel.geno().get(1, 127));
        REQUIRE(panel.geno().get(1, 128));
        REQUIRE(panel.geno().get(1, 129));
        REQUIRE(panel.geno().row(1)[2] == 3);
        REQUIRE(panel.geno().row(0)[2] == 0);
    }
}

TEST_CASE(
    "load_panel rejects genotypes that are not homozygous diploid calls",
    "[parentage][panel]")
{
    FileFixture files;
    PanelVcf vcf;

    SECTION("missing")
    {
        vcf.sites = {{"chr01", 205, "A", "G", {"0/0", "./.", "0/0"}}};
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring("chr01:205 P2: missing genotype"));
    }
    SECTION("heterozygous")
    {
        vcf.sites = {{"chr01", 205, "A", "G", {"0/0", "0/0", "0/1"}}};
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring("chr01:205 P3: heterozygous genotype 0/1"));
    }
    SECTION("haploid on one sample")
    {
        vcf.sites = {{"chr01", 205, "A", "G", {"1", "0/0", "0/0"}}};
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring("chr01:205 P1: haploid genotype 1"));
    }
    SECTION("haploid on every sample")
    {
        vcf.sites = {{"chr01", 205, "A", "G", {"1", "0", "0"}}};
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring(
                "chr01:205: expected diploid GT for every parent"));
    }
}

TEST_CASE(
    "load_panel rejects sites that are not bi-allelic ACGT SNPs",
    "[parentage][panel]")
{
    FileFixture files;
    PanelVcf vcf;

    SECTION("three alleles")
    {
        vcf.sites = {{"chr01", 205, "A", "G,T", {"0/0", "1/1", "2/2"}}};
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring(
                "chr01:205: panel sites must be bi-allelic, found 3 "
                "allele(s)"));
    }
    SECTION("indel")
    {
        vcf.sites = {{"chr01", 205, "AT", "A", {"0/0", "1/1", "0/0"}}};
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring(
                "chr01:205: panel sites must be SNPs with A/C/G/T alleles, "
                "found AT"));
    }
    SECTION("N allele")
    {
        vcf.sites = {{"chr01", 205, "A", "N", {"0/0", "1/1", "0/0"}}};
        REQUIRE_THROWS_WITH(load(files, vcf), ContainsSubstring("found N"));
    }
}

TEST_CASE("load_panel rejects sites out of order", "[parentage][panel]")
{
    FileFixture files;
    PanelVcf vcf;

    SECTION("position decreases")
    {
        vcf.sites = {
            {"chr01", 205, "A", "G", {"0/0", "0/0", "0/0"}},
            {"chr01", 101, "A", "G", {"0/0", "0/0", "0/0"}},
        };
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring(
                "chr01:101: panel positions must be strictly increasing"));
    }
    SECTION("position repeats")
    {
        vcf.sites = {
            {"chr01", 205, "A", "G", {"0/0", "0/0", "0/0"}},
            {"chr01", 205, "C", "T", {"0/0", "0/0", "0/0"}},
        };
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring(
                "chr01:205: panel positions must be strictly increasing"));
    }
    SECTION("contig interleaved")
    {
        vcf.sites = {
            {"chr01", 1, "A", "G", {"0/0", "0/0", "0/0"}},
            {"chr02", 1, "A", "G", {"0/0", "0/0", "0/0"}},
            {"chr01", 2, "A", "G", {"0/0", "0/0", "0/0"}},
        };
        REQUIRE_THROWS_WITH(
            load(files, vcf), ContainsSubstring("Contig chr01 appears twice"));
    }
}

TEST_CASE(
    "load_panel rejects files that cannot describe a panel",
    "[parentage][panel]")
{
    FileFixture files;
    PanelVcf vcf;

    SECTION("contig without a length")
    {
        vcf.contigs = {{"chr01", 0}};
        vcf.sites = {{"chr01", 1, "A", "G", {"0/0", "0/0", "0/0"}}};
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring(
                "Panel VCF header has no length for contig chr01"));
    }
    SECTION("no samples")
    {
        vcf.samples = {};
        REQUIRE_THROWS_WITH(
            load(files, vcf), ContainsSubstring("Panel VCF has no samples"));
    }
    SECTION("no sites")
    {
        REQUIRE_THROWS_WITH(
            load(files, vcf),
            ContainsSubstring("Panel contains no usable sites"));
    }
    SECTION("file does not exist")
    {
        REQUIRE_THROWS_WITH(
            load_panel(files.missing_path(".vcf"), {}),
            ContainsSubstring("Failed to open panel VCF"));
    }
}

TEST_CASE(
    "load_panel drops sites listed in the exclusion file",
    "[parentage][panel]")
{
    FileFixture files;
    PanelVcf vcf;
    vcf.sites = {
        {"chr01", 101, "A", "G", {"1/1", "0/0", "0/0"}},
        {"chr01", 205, "C", "T", {"0/0", "1/1", "1/1"}},
        {"chr02", 7, "G", "A", {"0/0", "0/0", "0/0"}},
    };
    PanelOptions options;

    SECTION("one site")
    {
        options.exclude_sites_path = files.create_text_file("chr01\t101\n");
        const Panel panel = load(files, vcf, options);
        REQUIRE(panel.n_sites() == 2);
        REQUIRE(panel.n_excluded_sites() == 1);
        REQUIRE(panel.sites().pos(0) == 204);
        REQUIRE(panel.geno().get(1, 0));
    }
    SECTION("every site of a contig removes the contig")
    {
        options.exclude_sites_path
            = files.create_text_file("# comment\nchr02\t7\n");
        const Panel panel = load(files, vcf, options);
        REQUIRE(panel.n_sites() == 2);
        REQUIRE(panel.sites().contigs().size() == 1);
        REQUIRE(panel.sites().contigs()[0].info.name == "chr01");
    }
    SECTION("an excluded site is not validated")
    {
        vcf.sites[1].gts = {"0/0", "0/1", "./."};
        options.exclude_sites_path = files.create_text_file("chr01\t205\n");
        const Panel panel = load(files, vcf, options);
        REQUIRE(panel.n_sites() == 2);
    }
    SECTION("malformed line")
    {
        options.exclude_sites_path = files.create_text_file("chr01 101\n");
        REQUIRE_THROWS_WITH(
            load(files, vcf, options),
            ContainsSubstring(":1: expected `contig<TAB>pos`"));
    }
}

TEST_CASE(
    "load_panel reports near-duplicate parents without merging them",
    "[parentage][panel]")
{
    FileFixture files;
    PanelVcf vcf;
    // P1 == P2; P3 differs from both at the last two of four sites.
    vcf.sites = {
        {"chr01", 1, "A", "G", {"1/1", "1/1", "1/1"}},
        {"chr01", 2, "A", "G", {"1/1", "1/1", "1/1"}},
        {"chr01", 3, "A", "G", {"0/0", "0/0", "1/1"}},
        {"chr01", 4, "A", "G", {"0/0", "0/0", "1/1"}},
    };
    PanelOptions options;

    SECTION("default threshold reports only the identical pair")
    {
        const Panel panel = load(files, vcf, options);
        REQUIRE(panel.n_parents() == 3);
        REQUIRE(panel.near_duplicates().size() == 1);
        REQUIRE(panel.near_duplicates()[0].parent_a == 0);
        REQUIRE(panel.near_duplicates()[0].parent_b == 1);
        REQUIRE(panel.near_duplicates()[0].n_discordant == 0);
    }
    SECTION("threshold at exactly the discordance includes the pair")
    {
        options.duplicate_max_discordance = 0.5;  // 2 of 4 sites
        const Panel panel = load(files, vcf, options);
        REQUIRE(panel.near_duplicates().size() == 3);
        REQUIRE(panel.near_duplicates()[1].parent_b == 2);
        REQUIRE(panel.near_duplicates()[1].n_discordant == 2);
    }
    SECTION("threshold just below the discordance excludes the pair")
    {
        options.duplicate_max_discordance = 0.49;
        const Panel panel = load(files, vcf, options);
        REQUIRE(panel.near_duplicates().size() == 1);
    }
}

TEST_CASE(
    "Panel rejects containers that disagree on their dimensions",
    "[parentage][panel]")
{
    SiteTable sites;
    sites.add_contig({"chr01", 10});
    sites.add_site(0, 0, 1, 0.0F);
    sites.finish();

    SECTION("names versus rows")
    {
        REQUIRE_THROWS_WITH(
            Panel(sites, {"a", "b"}, BitRows(1, 1), 0.001, 0),
            ContainsSubstring("2 parent names but 1 genotype rows"));
    }
    SECTION("sites versus columns")
    {
        REQUIRE_THROWS_WITH(
            Panel(sites, {"a"}, BitRows(1, 3), 0.001, 0),
            ContainsSubstring("1 sites but 3 genotype columns"));
    }
}
