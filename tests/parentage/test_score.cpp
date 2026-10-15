#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <string>
#include <vector>

#include "parentage/bitrows.h"
#include "parentage/counts.h"
#include "parentage/panel.h"
#include "parentage/score.h"
#include "parentage/sites.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using vcfbox::parentage::BitRows;
using vcfbox::parentage::Panel;
using vcfbox::parentage::SampleCounts;
using vcfbox::parentage::SampleScores;
using vcfbox::parentage::score_sample;
using vcfbox::parentage::ScoreOptions;
using vcfbox::parentage::SiteTable;

namespace
{
// `n_sites` A>G sites on chr01 with the given allele1 frequencies.
auto make_sites(const std::vector<float>& freq1) -> SiteTable
{
    SiteTable sites;
    sites.add_contig({"chr01", 1000});
    for (size_t s = 0; s < freq1.size(); ++s)
    {
        sites.add_site(static_cast<int64_t>(s), 0, 2, freq1[s]);
    }
    sites.finish();
    return sites;
}

// Parent rows given as strings of '0' / '1' / 'h' (heterozygous) per site.
auto make_panel(
    const std::vector<float>& freq1,
    const std::vector<std::string>& rows) -> Panel
{
    BitRows geno(rows.size(), freq1.size());
    BitRows het(rows.size(), freq1.size());
    for (size_t r = 0; r < rows.size(); ++r)
    {
        for (size_t s = 0; s < rows[r].size(); ++s)
        {
            if (rows[r][s] == '1')
            {
                geno.set(r, s);
            }
            if (rows[r][s] == 'h')
            {
                het.set(r, s);
            }
        }
    }
    std::vector<std::string> names;
    for (size_t r = 0; r < rows.size(); ++r)
    {
        names.push_back(std::format("P{}", r + 1));
    }
    return Panel(
        make_sites(freq1), names, std::move(geno), std::move(het), 0.0, 0);
}

auto make_counts(const std::vector<uint8_t>& n0, const std::vector<uint8_t>& n1)
    -> SampleCounts
{
    SampleCounts counts("S", n0.size());
    for (size_t s = 0; s < n0.size(); ++s)
    {
        for (uint8_t i = 0; i < n0[s]; ++i)
        {
            counts.add0(s);
        }
        for (uint8_t i = 0; i < n1[s]; ++i)
        {
            counts.add1(s);
        }
    }
    return counts;
}

// Three parents, six sites, two of them uncovered:
//   site      0  1  2  3  4  5
//   P1        0  0  0  0  0  0
//   P2        1  1  0  0  0  0
//   P3        1  1  1  1  0  0
//   n0        2  0  1  0  3  0
//   n1        0  1  1  0  0  0
auto three_parent_scores(ScoreOptions options = {}) -> SampleScores
{
    const Panel panel = make_panel(
        {2.0F / 3, 2.0F / 3, 1.0F / 3, 1.0F / 3, 0, 0},
        {"000000", "110000", "111100"});
    return score_sample(
        panel, make_counts({2, 0, 1, 0, 3, 0}, {0, 1, 1, 0, 0, 0}), options);
}
}  // namespace

TEST_CASE("score_sample splits reads by parent pair", "[parentage][score]")
{
    const SampleScores scores = three_parent_scores();

    REQUIRE(scores.n_covered == 4);
    REQUIRE(scores.pairs.size() == 6);

    SECTION("ranks pairs by score, ties by index")
    {
        REQUIRE(scores.pairs[0].a == 0);
        REQUIRE(scores.pairs[0].b == 2);
        REQUIRE(scores.pairs[0].score == 2049);
        REQUIRE(scores.pairs[1].a == 0);
        REQUIRE(scores.pairs[1].b == 1);
        REQUIRE(scores.pairs[1].score == -1180);
        REQUIRE(scores.pairs[2].a == 0);
        REQUIRE(scores.pairs[2].b == 0);
        REQUIRE(scores.pairs[2].score == -3726);
        REQUIRE(scores.pairs[3].a == 1);
        REQUIRE(scores.pairs[3].b == 2);
        REQUIRE(scores.pairs[3].score == -5092);
        REQUIRE(scores.pairs[4].a == 1);
        REQUIRE(scores.pairs[4].b == 1);
        REQUIRE(scores.pairs[4].score == -8321);
        REQUIRE(scores.pairs[5].a == 2);
        REQUIRE(scores.pairs[5].b == 2);
        REQUIRE(scores.pairs[5].score == -8321);
    }
    SECTION("counts of a mixed pair")
    {
        const auto& p13 = scores.pairs[0];
        REQUIRE(p13.hom_conc == 3);
        REQUIRE(p13.hom_disc == 0);
        REQUIRE(p13.het_0 == 3);
        REQUIRE(p13.het_1 == 2);
    }
    SECTION("counts of a selfing pair have no heterozygous sites")
    {
        const auto& p22 = scores.pairs[4];
        REQUIRE(p22.hom_conc == 5);
        REQUIRE(p22.hom_disc == 3);
        REQUIRE(p22.het_0 == 0);
        REQUIRE(p22.het_1 == 0);
        REQUIRE(std::isnan(scores.top[4].het_alt_fraction));
    }
    SECTION("float statistics of the top pairs")
    {
        REQUIRE(scores.top.size() == 6);
        REQUIRE_THAT(scores.top[0].lod_vs_null, WithinAbs(0.60753, 1e-4));
        REQUIRE_THAT(scores.top[1].lod_vs_null, WithinAbs(-2.62140, 1e-4));
        REQUIRE_THAT(scores.top[0].hom_mismatch_rate, WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(scores.top[0].het_alt_fraction, WithinAbs(0.4, 1e-12));
        REQUIRE_THAT(scores.top[1].hom_mismatch_rate, WithinAbs(0.2, 1e-12));
        REQUIRE_THAT(scores.lod_vs_second, WithinAbs(3.229, 1e-12));
    }
    SECTION("single-parent model")
    {
        REQUIRE(scores.singles.size() == 3);
        REQUIRE(scores.singles[0].parent == 0);
        REQUIRE_THAT(scores.singles[0].lod_vs_null, WithinAbs(-0.13059, 1e-4));
        REQUIRE(scores.singles[1].parent == 2);
        REQUIRE_THAT(scores.singles[1].lod_vs_null, WithinAbs(-0.37335, 1e-4));
        REQUIRE(scores.singles[2].parent == 1);
        REQUIRE_THAT(scores.singles[2].lod_vs_null, WithinAbs(-1.00988, 1e-4));
        // Sites with two reads: 0, 2, 4. P3 carries allele1 at site 0, where
        // both reads are allele0.
        REQUIRE_THAT(scores.singles[0].allele_absent_rate, WithinAbs(0, 1e-12));
        REQUIRE_THAT(
            scores.singles[1].allele_absent_rate, WithinAbs(1.0 / 3, 1e-12));
    }
}

// As above, but P2 is heterozygous at site 2 (allele1 dosage 3 of 6).
auto het_parent_scores() -> SampleScores
{
    const Panel panel = make_panel(
        {4.0F / 6, 4.0F / 6, 3.0F / 6, 2.0F / 6, 0, 0},
        {"000000", "11h000", "111100"});
    return score_sample(
        panel, make_counts({2, 0, 1, 0, 3, 0}, {0, 1, 1, 0, 0, 0}), {});
}

TEST_CASE(
    "score_sample leaves heterozygous parent sites out of the counts",
    "[parentage][score]")
{
    const SampleScores scores = het_parent_scores();

    SECTION("the site counts for neither pair of the parent")
    {
        // (P1,P2) and (P1,P3) tie at 2049; (P1,P2) first by index.
        const auto& p12 = scores.pairs[0];
        REQUIRE(p12.a == 0);
        REQUIRE(p12.b == 1);
        REQUIRE(p12.hom_conc == 3);
        REQUIRE(p12.hom_disc == 0);
        REQUIRE(p12.het_0 == 2);
        REQUIRE(p12.het_1 == 1);
        REQUIRE(p12.score == 2049);
        const auto& p22 = scores.pairs[3];
        REQUIRE(p22.a == 1);
        REQUIRE(p22.b == 1);
        REQUIRE(p22.hom_conc == 4);
        REQUIRE(p22.hom_disc == 2);
    }
    SECTION("lod_vs_null adds the exact mixture at the uncertain site")
    {
        REQUIRE_THAT(scores.top[0].lod_vs_null, WithinAbs(-0.15509, 1e-4));
        REQUIRE_THAT(scores.top[1].lod_vs_null, WithinAbs(0.49922, 1e-4));
    }
    SECTION("single-parent model marginalises the heterozygous site")
    {
        REQUIRE(scores.singles[0].parent == 0);
        REQUIRE_THAT(scores.singles[0].lod_vs_null, WithinAbs(0.12919, 1e-4));
        REQUIRE(scores.singles[1].parent == 1);
        REQUIRE_THAT(scores.singles[1].lod_vs_null, WithinAbs(-0.75011, 1e-4));
        // Two-read sites 0, 2, 4; site 2 is heterozygous for P2 and dropped.
        REQUIRE_THAT(
            scores.singles[1].allele_absent_rate, WithinAbs(0.5, 1e-12));
    }
}

TEST_CASE(
    "score_sample uses the error rate as read weights",
    "[parentage][score]")
{
    const Panel panel = make_panel({0}, {"0"});

    SECTION("a concordant read at e = 0.01 scores +683")
    {
        const SampleScores scores
            = score_sample(panel, make_counts({1}, {0}), {});
        REQUIRE(scores.pairs[0].score == 683);
    }
    SECTION("a discordant read at e = 0.01 scores -3912")
    {
        const SampleScores scores
            = score_sample(panel, make_counts({0}, {1}), {});
        REQUIRE(scores.pairs[0].score == -3912);
    }
    SECTION("e = 0.1 gives +588 / -1609")
    {
        ScoreOptions options;
        options.error_rate = 0.1;
        REQUIRE(
            score_sample(panel, make_counts({1}, {0}), options).pairs[0].score
            == 588);
        REQUIRE(
            score_sample(panel, make_counts({0}, {1}), options).pairs[0].score
            == -1609);
    }
    SECTION("error rate outside (0, 0.5) is rejected")
    {
        ScoreOptions options;
        options.error_rate = 0.5;
        REQUIRE_THROWS_WITH(
            score_sample(panel, make_counts({1}, {0}), options),
            ContainsSubstring("error rate must be in (0, 0.5)"));
    }
}

TEST_CASE(
    "score_sample handles word boundaries and saturation",
    "[parentage][score]")
{
    SECTION("a read at site 64 lands in the second word")
    {
        const Panel panel
            = make_panel(std::vector<float>(65, 0), {std::string(65, '0')});
        std::vector<uint8_t> n0(65, 0);
        n0[64] = 1;
        const SampleScores scores = score_sample(
            panel, make_counts(n0, std::vector<uint8_t>(65, 0)), {});
        REQUIRE(scores.pairs[0].hom_conc == 1);
        REQUIRE(scores.pairs[0].hom_disc == 0);
        REQUIRE(scores.n_covered == 1);
    }
    SECTION("a word whose only counts are even is not skipped")
    {
        const Panel panel = make_panel({0}, {"0"});
        const SampleScores scores
            = score_sample(panel, make_counts({2}, {2}), {});
        REQUIRE(scores.pairs[0].hom_conc == 2);
        REQUIRE(scores.pairs[0].hom_disc == 2);
    }
    SECTION("a saturated count of 255 uses all eight planes")
    {
        const Panel panel = make_panel({0}, {"0"});
        const SampleScores scores
            = score_sample(panel, make_counts({255}, {0}), {});
        REQUIRE(scores.pairs[0].hom_conc == 255);
    }
    SECTION("a single pair has no second")
    {
        const Panel panel = make_panel({0}, {"0"});
        const SampleScores scores
            = score_sample(panel, make_counts({1}, {0}), {});
        REQUIRE(std::isnan(scores.lod_vs_second));
    }
    SECTION("top is capped by the option")
    {
        ScoreOptions options;
        options.top = 2;
        REQUIRE(three_parent_scores(options).top.size() == 2);
    }
}

TEST_CASE(
    "score_sample pairs only the best single parents",
    "[parentage][score]")
{
    // Single-parent order is P1, P3, P2 (see the three-parent case).
    ScoreOptions options;
    options.candidates = 2;

    const SampleScores scores = three_parent_scores(options);

    SECTION("keeps every pair of the candidates, selfing included")
    {
        REQUIRE(scores.pairs.size() == 3);
        REQUIRE(scores.pairs[0].a == 0);
        REQUIRE(scores.pairs[0].b == 2);
        REQUIRE(scores.pairs[0].score == 2049);
        REQUIRE(scores.pairs[1].a == 0);
        REQUIRE(scores.pairs[1].b == 0);
        REQUIRE(scores.pairs[2].a == 2);
        REQUIRE(scores.pairs[2].b == 2);
    }
    SECTION("the second best is taken among the candidates")
    {
        REQUIRE_THAT(scores.lod_vs_second, WithinAbs(5.775, 1e-9));
    }
    SECTION("every single parent is still scored")
    {
        REQUIRE(scores.singles.size() == 3);
    }
    SECTION("zero means every parent")
    {
        options.candidates = 0;
        REQUIRE(three_parent_scores(options).pairs.size() == 6);
    }
    SECTION("more candidates than parents means every parent")
    {
        options.candidates = 50;
        REQUIRE(three_parent_scores(options).pairs.size() == 6);
    }
}

TEST_CASE("score_sample rejects counts of another panel", "[parentage][score]")
{
    const Panel panel = make_panel({0, 0}, {"00"});
    REQUIRE_THROWS_WITH(
        score_sample(panel, make_counts({1}, {0}), {}),
        ContainsSubstring("has 1 sites but the panel has 2"));
}
