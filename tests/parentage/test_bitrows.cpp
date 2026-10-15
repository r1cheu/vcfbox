#include <catch2/catch_test_macros.hpp>

#include "parentage/bitrows.h"

using vcfbox::parentage::BitRows;
using vcfbox::parentage::popcount_xor;

TEST_CASE(
    "BitRows rounds the row width up to whole words",
    "[parentage][bitrows]")
{
    SECTION("exact multiple of 64")
    {
        const BitRows rows(2, 128);
        REQUIRE(rows.rows() == 2);
        REQUIRE(rows.n_bits() == 128);
        REQUIRE(rows.words_per_row() == 2);
        REQUIRE(rows.row(1).size() == 2);
    }
    SECTION("one bit past a word boundary")
    {
        const BitRows rows(1, 65);
        REQUIRE(rows.words_per_row() == 2);
    }
    SECTION("empty")
    {
        const BitRows rows;
        REQUIRE(rows.rows() == 0);
        REQUIRE(rows.words_per_row() == 0);
    }
}

TEST_CASE("BitRows set and get address a single bit", "[parentage][bitrows]")
{
    BitRows rows(2, 100);

    rows.set(1, 70);

    REQUIRE(rows.get(1, 70));
    REQUIRE_FALSE(rows.get(1, 69));
    REQUIRE_FALSE(rows.get(1, 71));
    REQUIRE_FALSE(rows.get(0, 70));
    REQUIRE(rows.row(1)[1] == (uint64_t{1} << 6U));
    REQUIRE(rows.row(1)[0] == 0);
}

TEST_CASE(
    "popcount_xor counts differing bits across all words",
    "[parentage][bitrows]")
{
    BitRows rows(2, 128);
    rows.set(0, 0);
    rows.set(0, 64);
    rows.set(1, 64);
    rows.set(1, 127);

    REQUIRE(popcount_xor(rows.row(0), rows.row(1)) == 2);
    REQUIRE(popcount_xor(rows.row(0), rows.row(0)) == 0);
}
