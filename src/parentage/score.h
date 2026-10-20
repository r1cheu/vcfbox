#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "parentage/counts.h"
#include "parentage/panel.h"

namespace vcfbox::parentage
{
struct ScoreOptions
{
    double error_rate = 0.01;  // per-read probability of the wrong allele
    size_t top = 10;           // best pairs kept in SampleScores::top
    // Parents with the best single-parent lod among which the pairs are
    // formed; 0 pairs every parent.
    size_t candidates = 20;
};

// A sample's reads split by what parents `a` and `b` carry at each site.
// Sites where either parent is heterozygous or missing are left out of every
// count.
struct PairStats
{
    uint32_t a = 0;  // a <= b; a == b is selfing
    uint32_t b = 0;
    uint64_t hom_conc = 0;  // both parents carry the read's allele
    uint64_t hom_disc = 0;  // both parents carry the other allele
    uint64_t het_0 = 0;     // parents differ; read is allele0
    uint64_t het_1 = 0;     // parents differ; read is allele1
};

// Float statistics of a pair. Rates are NaN when their
// denominator is 0; `het_alt_fraction` is therefore NaN for selfing.
struct PairReport
{
    PairStats stats;
    // ln L(pair) - ln L(random panel pair) over all covered sites; at sites
    // where a parent is heterozygous or missing the F1 genotype is
    // marginalised over the two gametes.
    double lod_vs_null = 0;
    double hom_mismatch_rate = 0;  // hom_disc / (hom_conc + hom_disc)
    double het_alt_fraction = 0;   // het_1 / (het_0 + het_1)
};

// One parent with the other one marginalised over the panel allele
// frequencies (and its own heterozygous sites over the two gametes; its
// missing sites score as the null).
struct SingleStats
{
    uint32_t parent = 0;
    double lod_vs_null = 0;
    // Over the sites with at least two reads where this parent is
    // homozygous: fraction where every read carries the other allele; NaN
    // when there is no such site.
    double allele_absent_rate = 0;
};

// The expected pair of a sample and its rank among every scored pair
// (1 = best).
struct ExpectedScore
{
    PairReport report;
    size_t rank = 0;
};

struct SampleScores
{
    size_t n_covered = 0;  // sites with at least one allele0/allele1 read
    // The `ScoreOptions::top` best of the pairs a <= b of the candidate
    // parents (plus the expected parents, when given), best lod_vs_null
    // first; ties broken by (a, b).
    std::vector<PairReport> top;
    // lod_vs_null of the best pair minus that of the second; NaN with a
    // single pair.
    double lod_vs_second = 0;
    // Every parent, best first; ties broken by index.
    std::vector<SingleStats> singles;
    std::optional<ExpectedScore> expected;  // set when an expected pair is given
};

// Two parent indices; order does not matter.
struct ParentPair
{
    uint32_t a = 0;
    uint32_t b = 0;
};

// Scores every single parent and the pairs of the candidate parents against
// the sample's counts. The parents of `expected` are always candidates. Throws when the counts and
// the panel differ in site number.
[[nodiscard]] auto score_sample(
    const Panel& panel,
    const SampleCounts& counts,
    const ScoreOptions& options,
    std::optional<ParentPair> expected = std::nullopt) -> SampleScores;
}  // namespace vcfbox::parentage
