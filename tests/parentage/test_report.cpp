#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "file_fixture.h"
#include "parentage/report.h"
#include "parentage/score.h"

using Catch::Matchers::ContainsSubstring;
using vcfbox::parentage::Ambiguity;
using vcfbox::parentage::assign_fit;
using vcfbox::parentage::assign_status;
using vcfbox::parentage::CallOptions;
using vcfbox::parentage::Contest;
using vcfbox::parentage::ExpectedParents;
using vcfbox::parentage::find_ambiguity;
using vcfbox::parentage::Fit;
using vcfbox::parentage::model_margin;
using vcfbox::parentage::PairReport;
using vcfbox::parentage::PairStats;
using vcfbox::parentage::read_expected_pairs;
using vcfbox::parentage::SampleScores;
using vcfbox::parentage::single_lod_vs_second;
using vcfbox::parentage::SingleStats;
using vcfbox::parentage::Status;
using vcfbox::parentage::to_string;
using vcfbox::test::FileFixture;

namespace
{
// Best pair P0 x P1 with mismatch 0.01, lod_vs_null 1000, lod_vs_second 100;
// singles P0 at 500 and P2 at 400: every margin is clear.
auto passing_scores() -> SampleScores
{
    SampleScores scores;
    PairStats best;
    best.a = 0;
    best.b = 1;
    best.hom_conc = 19800;
    best.hom_disc = 200;
    best.het_0 = 5000;
    best.het_1 = 5000;
    scores.top = {PairReport{best, 1000, 0.01, 0.5}};
    scores.lod_vs_second = 100;
    scores.singles = {SingleStats{0, 500, 0.02}, SingleStats{2, 400, 0.02}};
    return scores;
}
}  // namespace

TEST_CASE(
    "model_margin is the leading model's lod over the best other model",
    "[parentage][report]")
{
    SampleScores scores = passing_scores();

    SECTION("the pair leading is measured against the single")
    {
        REQUIRE(model_margin(scores) == 500);
    }
    SECTION("the single leading is measured against the pair")
    {
        scores.top[0].lod_vs_null = 100;
        REQUIRE(model_margin(scores) == 400);
    }
    SECTION("no parent leading is measured against the better of the two")
    {
        scores.top[0].lod_vs_null = -100;
        scores.singles[0].lod_vs_null = -200;
        REQUIRE(model_margin(scores) == 100);
    }
}

TEST_CASE(
    "single_lod_vs_second is the best single over the second",
    "[parentage][report]")
{
    SampleScores scores = passing_scores();

    SECTION("two singles give their difference")
    {
        REQUIRE(single_lod_vs_second(scores) == 100);
    }
    SECTION("a single parent gives NaN")
    {
        scores.singles.pop_back();
        REQUIRE(std::isnan(single_lod_vs_second(scores)));
    }
}

TEST_CASE(
    "find_ambiguity names the first margin below min_lod",
    "[parentage][report]")
{
    const CallOptions options;
    SampleScores scores = passing_scores();

    SECTION("every margin clear finds nothing")
    {
        REQUIRE_FALSE(find_ambiguity(scores, options).has_value());
    }
    SECTION("a pair close to the single fails the model contest")
    {
        scores.singles[0].lod_vs_null = 970;
        const std::optional<Ambiguity> found = find_ambiguity(scores, options);
        REQUIRE(found.has_value());
        REQUIRE(found->contest == Contest::model);
        REQUIRE(found->margin == 30);
    }
    SECTION("the model contest is checked before the pair contest")
    {
        scores.singles[0].lod_vs_null = 970;
        scores.lod_vs_second = 10;
        REQUIRE(find_ambiguity(scores, options)->contest == Contest::model);
    }
    SECTION("a close second pair fails the pair contest")
    {
        scores.lod_vs_second = 10;
        const std::optional<Ambiguity> found = find_ambiguity(scores, options);
        REQUIRE(found.has_value());
        REQUIRE(found->contest == Contest::pair);
        REQUIRE(found->margin == 10);
    }
    SECTION("a close second single fails the single contest")
    {
        scores.top[0].lod_vs_null = 100;
        scores.singles[1].lod_vs_null = 480;
        const std::optional<Ambiguity> found = find_ambiguity(scores, options);
        REQUIRE(found.has_value());
        REQUIRE(found->contest == Contest::single);
        REQUIRE(found->margin == 20);
    }
    SECTION("a close second single is ignored while the pair leads")
    {
        scores.singles[1].lod_vs_null = 480;
        REQUIRE_FALSE(find_ambiguity(scores, options).has_value());
    }
    SECTION("contests are written in lower case")
    {
        REQUIRE(to_string(Contest::single) == "single");
    }
}

TEST_CASE(
    "assign_status names the leading model when every margin is clear",
    "[parentage][report]")
{
    const CallOptions options;
    SampleScores scores = passing_scores();

    SECTION("a pair clear of the single and of the second pair is TWO_PARENT")
    {
        REQUIRE(assign_status(scores, options) == Status::two_parent);
    }
    SECTION("a selfing best pair is TWO_PARENT")
    {
        scores.top[0].stats.b = 0;
        REQUIRE(assign_status(scores, options) == Status::two_parent);
    }
    SECTION("a pair less than min_lod above the single is AMBIGUOUS")
    {
        scores.singles[0].lod_vs_null = 950.1;
        REQUIRE(assign_status(scores, options) == Status::ambiguous);
    }
    SECTION("a pair exactly min_lod above the single is TWO_PARENT")
    {
        scores.singles[0].lod_vs_null = 950;
        REQUIRE(assign_status(scores, options) == Status::two_parent);
    }
    SECTION("lod_vs_second below min_lod is AMBIGUOUS")
    {
        scores.lod_vs_second = 49.9;
        REQUIRE(assign_status(scores, options) == Status::ambiguous);
    }
    SECTION("a single pair (NaN lod_vs_second) is TWO_PARENT")
    {
        scores.lod_vs_second = std::numeric_limits<double>::quiet_NaN();
        REQUIRE(assign_status(scores, options) == Status::two_parent);
    }
    SECTION("a single clear of the pair and of the second single is ONE_PARENT")
    {
        scores.top[0].lod_vs_null = 100;
        REQUIRE(assign_status(scores, options) == Status::one_parent);
    }
    SECTION("a leading single outside the best pair is ONE_PARENT")
    {
        scores.top[0].lod_vs_null = 100;
        scores.singles = {SingleStats{2, 500, 0.02}, SingleStats{0, 400, 0.02}};
        REQUIRE(assign_status(scores, options) == Status::one_parent);
    }
    SECTION("a single less than min_lod above the second single is AMBIGUOUS")
    {
        scores.top[0].lod_vs_null = 100;
        scores.singles[1].lod_vs_null = 450.1;
        REQUIRE(assign_status(scores, options) == Status::ambiguous);
    }
    SECTION("a single parent in the panel (NaN second) is ONE_PARENT")
    {
        scores.top[0].lod_vs_null = 100;
        scores.singles.pop_back();
        REQUIRE(assign_status(scores, options) == Status::one_parent);
    }
    SECTION("pair and single min_lod below the null is ZERO_PARENT")
    {
        scores.top[0].lod_vs_null = -50;
        scores.singles[0].lod_vs_null = -60;
        REQUIRE(assign_status(scores, options) == Status::zero_parent);
    }
    SECTION("a model less than min_lod from the null is AMBIGUOUS")
    {
        scores.top[0].lod_vs_null = 30;
        scores.singles[0].lod_vs_null = 10;
        REQUIRE(assign_status(scores, options) == Status::ambiguous);
    }
    SECTION("statuses are written in upper case")
    {
        REQUIRE(to_string(Status::zero_parent) == "ZERO_PARENT");
    }
}

TEST_CASE(
    "assign_fit judges the best pair only when the pair model leads",
    "[parentage][report]")
{
    const CallOptions options;
    SampleScores scores = passing_scores();

    SECTION("hom_mismatch at the maximum is PASS")
    {
        scores.top[0].hom_mismatch_rate = 0.03;
        REQUIRE(assign_fit(scores, options) == Fit::pass);
    }
    SECTION("hom_mismatch above the maximum is FAIL")
    {
        scores.top[0].hom_mismatch_rate = 0.031;
        REQUIRE(assign_fit(scores, options) == Fit::fail);
    }
    SECTION("a leading pair within min_lod of the single is still judged")
    {
        scores.singles[0].lod_vs_null = 990;
        scores.top[0].hom_mismatch_rate = 0.031;
        REQUIRE(assign_fit(scores, options) == Fit::fail);
    }
    SECTION("a leading single is NA")
    {
        scores.top[0].lod_vs_null = 100;
        scores.top[0].hom_mismatch_rate = 0.5;
        REQUIRE(assign_fit(scores, options) == Fit::na);
    }
    SECTION("a NaN hom_mismatch is NA")
    {
        scores.top[0].hom_mismatch_rate
            = std::numeric_limits<double>::quiet_NaN();
        REQUIRE(assign_fit(scores, options) == Fit::na);
    }
    SECTION("fits are written as PASS, FAIL and NA")
    {
        REQUIRE(to_string(Fit::pass) == "PASS");
        REQUIRE(to_string(Fit::fail) == "FAIL");
        REQUIRE(to_string(Fit::na) == "NA");
    }
}

TEST_CASE(
    "read_expected_pairs resolves parent names against the panel",
    "[parentage][report]")
{
    FileFixture files;
    const std::vector<std::string> names{"P1", "P2", "P3"};

    SECTION("pairs, selfing and absent parents")
    {
        const auto expected = read_expected_pairs(
            files.create_text_file(
                "sample\tmother\tfather\n"
                "S1\tP3\tP1\n"
                "\n"
                "# comment\n"
                "S2\tP2\tP2\n"
                "S3\tP1\tQ9\n"),
            names);

        REQUIRE(expected.size() == 3);
        const ExpectedParents& s1 = expected.at("S1");
        REQUIRE(s1.a == "P3");
        REQUIRE(s1.b == "P1");
        REQUIRE(s1.pair->a == 2);
        REQUIRE(s1.pair->b == 0);
        REQUIRE(expected.at("S2").pair->a == 1);
        REQUIRE(expected.at("S2").pair->b == 1);
        REQUIRE(expected.at("S3").b == "Q9");
        REQUIRE_FALSE(expected.at("S3").pair);
    }
    SECTION("a line without three fields is rejected")
    {
        REQUIRE_THROWS_WITH(
            read_expected_pairs(files.create_text_file("S1\tP1\n"), names),
            ContainsSubstring(":1: expected `sample<TAB>mother<TAB>father`"));
    }
    SECTION("a repeated sample is rejected")
    {
        REQUIRE_THROWS_WITH(
            read_expected_pairs(
                files.create_text_file("S1\tP1\tP2\nS1\tP1\tP3\n"), names),
            ContainsSubstring(":2: sample S1 is listed twice"));
    }
}
