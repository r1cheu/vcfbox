#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "combine/combine_genotypes.h"
#include "file_fixture.h"
#include "vcf_fixture.h"

using Catch::Matchers::ContainsSubstring;
using vcfbox::combine_genotypes;
using vcfbox::test::FileFixture;
using vcfbox::test::PanelVcf;

namespace
{
auto file_text(const std::filesystem::path& path) -> std::string
{
    std::ifstream in(path);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}
}  // namespace

TEST_CASE(
    "combine_genotypes writes one homozygous-pair genotype per sample pair",
    "[combine]")
{
    FileFixture files;
    PanelVcf vcf;
    vcf.sites = {{"chr01", 5, "A", "G", {"0/0", "1/1", "0/1"}}};
    const auto vcf_path = files.create_text_file(vcf.text(), ".vcf");
    const auto out_path = files.missing_path(".vcf");

    SECTION("without old samples, only the pair column is written")
    {
        combine_genotypes(
            vcf_path.string(), {{"P1", "P2"}}, false, out_path.string());

        const std::string out = file_text(out_path);
        CAPTURE(out);
        REQUIRE_THAT(out, ContainsSubstring("FORMAT\tP1~P2\n"));
        REQUIRE_THAT(
            out, ContainsSubstring("chr01\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1\n"));
    }

    SECTION("with old samples, heterozygous old genotypes become missing")
    {
        combine_genotypes(
            vcf_path.string(), {{"P1", "P2"}}, true, out_path.string());

        const std::string out = file_text(out_path);
        CAPTURE(out);
        REQUIRE_THAT(out, ContainsSubstring("FORMAT\tP1\tP2\tP3\tP1~P2\n"));
        REQUIRE_THAT(
            out,
            ContainsSubstring(
                "chr01\t5\t.\tA\tG\t.\t.\t.\tGT\t0/0\t1/1\t./.\t0/1\n"));
    }
}

TEST_CASE(
    "combine_genotypes rejects a record whose genotypes are not diploid",
    "[combine]")
{
    FileFixture files;
    const bool keep_old_samples = GENERATE(false, true);
    const auto out_path = files.missing_path(".vcf");
    PanelVcf vcf;

    SECTION("haploid record")
    {
        vcf.sites = {
            {"chr01", 3, "A", "G", {"0/0", "1/1", "0/0"}},
            {"chr02", 7, "C", "T", {"0", "1", "0"}},
        };
        const auto vcf_path = files.create_text_file(vcf.text(), ".vcf");

        REQUIRE_THROWS_WITH(
            combine_genotypes(
                vcf_path.string(),
                {{"P1", "P2"}},
                keep_old_samples,
                out_path.string()),
            ContainsSubstring("chr02:7 (ploidy 1)"));
    }

    SECTION("triploid record")
    {
        vcf.sites = {
            {"chr01", 3, "A", "G", {"0/0", "1/1", "0/0"}},
            {"chr02", 7, "C", "T", {"0/0/1", "1/1/1", "0/0/0"}},
        };
        const auto vcf_path = files.create_text_file(vcf.text(), ".vcf");

        REQUIRE_THROWS_WITH(
            combine_genotypes(
                vcf_path.string(),
                {{"P1", "P2"}},
                keep_old_samples,
                out_path.string()),
            ContainsSubstring("chr02:7 (ploidy 3)"));
    }
}
