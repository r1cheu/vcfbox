#include "parentage/score.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <initializer_list>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "parentage/bitrows.h"
#include "parentage/sites.h"

namespace vcfbox::parentage
{
namespace
{
constexpr size_t n_planes = 8;  // counts are uint8
constexpr double nan = std::numeric_limits<double>::quiet_NaN();

// Per-read log-likelihoods: `conc` for a read matching a homozygous F1,
// `disc` for the other allele, `het` for any read at a heterozygous F1.
struct ReadWeights
{
    double conc;
    double disc;
    double het;
    int64_t score_conc;  // round(1000 * (conc - het))
    int64_t score_disc;  // round(1000 * (disc - het))
};

auto make_weights(double error_rate) -> ReadWeights
{
    if (error_rate <= 0 || error_rate >= 0.5)
    {
        throw std::runtime_error(
            std::format("error rate must be in (0, 0.5), got {}", error_rate));
    }
    const double conc = std::log1p(-error_rate);
    const double disc = std::log(error_rate);
    const double het = std::log(0.5);
    return {
        conc,
        disc,
        het,
        std::llround(1000 * (conc - het)),
        std::llround(1000 * (disc - het))};
}

// Log-likelihoods of one site's reads under the three F1 genotypes.
struct SiteLikelihoods
{
    double l00;
    double l01;
    double l11;
};

auto site_likelihoods(uint8_t n0, uint8_t n1, const ReadWeights& w)
    -> SiteLikelihoods
{
    return {
        n0 * w.conc + n1 * w.disc,
        (n0 + n1) * w.het,
        n1 * w.conc + n0 * w.disc};
}

auto log_sum_exp(std::initializer_list<double> terms) -> double
{
    const double max = std::ranges::max(terms);
    double sum = 0;
    for (const double t : terms)
    {
        sum += std::exp(t - max);
    }
    return max + std::log(sum);
}

auto log_or_neg_inf(double x) -> double
{
    return x > 0 ? std::log(x) : -std::numeric_limits<double>::infinity();
}

// Bit k of every n0 / n1 count, laid out like a panel row, so that a count
// weighted by a parent-pair mask is a popcount per plane.
struct CountPlanes
{
    std::array<std::vector<uint64_t>, n_planes> p0;
    std::array<std::vector<uint64_t>, n_planes> p1;
    std::vector<uint32_t> covered_words;  // words holding a covered site
};

auto build_planes(const SampleCounts& counts, size_t words) -> CountPlanes
{
    CountPlanes planes;
    for (size_t k = 0; k < n_planes; ++k)
    {
        planes.p0[k].assign(words, 0);
        planes.p1[k].assign(words, 0);
    }
    const std::span<const uint8_t> n0 = counts.n0();
    const std::span<const uint8_t> n1 = counts.n1();
    uint32_t last_covered_word = UINT32_MAX;
    for (size_t s = 0; s < n0.size(); ++s)
    {
        if (n0[s] + n1[s] == 0)
        {
            continue;
        }
        const auto word = static_cast<uint32_t>(s / 64);
        if (word != last_covered_word)
        {
            planes.covered_words.push_back(word);
            last_covered_word = word;
        }
        const uint64_t bit = uint64_t{1} << (s % 64);
        for (size_t k = 0; k < n_planes; ++k)
        {
            if (((n0[s] >> k) & 1) != 0)
            {
                planes.p0[k][word] |= bit;
            }
            if (((n1[s] >> k) & 1) != 0)
            {
                planes.p1[k][word] |= bit;
            }
        }
    }
    return planes;
}

auto popcount(uint64_t x) -> uint64_t
{
    return static_cast<uint64_t>(std::popcount(x));
}

// Genotype and heterozygosity rows of one parent.
struct ParentRows
{
    std::span<const uint64_t> geno;
    std::span<const uint64_t> het;
};

auto parent_rows(const Panel& panel, size_t parent) -> ParentRows
{
    return {panel.geno().row(parent), panel.het().row(parent)};
}

// The four read splits of one pair over the sites where neither parent is
// heterozygous. Padding bits past the last site are set in `both0`, which is
// harmless because the planes are zero there.
auto count_pair(
    const ParentRows& pi,
    const ParentRows& pj,
    const CountPlanes& planes) -> PairStats
{
    PairStats stats;
    for (const uint32_t w : planes.covered_words)
    {
        const uint64_t called = ~(pi.het[w] | pj.het[w]);
        const uint64_t both0 = ~(pi.geno[w] | pj.geno[w]) & called;
        const uint64_t both1 = pi.geno[w] & pj.geno[w] & called;
        const uint64_t mixed = (pi.geno[w] ^ pj.geno[w]) & called;
        for (size_t k = 0; k < n_planes; ++k)
        {
            const uint64_t p0 = planes.p0[k][w];
            const uint64_t p1 = planes.p1[k][w];
            if ((p0 | p1) == 0)
            {
                continue;
            }
            stats.hom_conc += (popcount(both0 & p0) + popcount(both1 & p1))
                              << k;
            stats.hom_disc += (popcount(both0 & p1) + popcount(both1 & p0))
                              << k;
            stats.het_0 += popcount(mixed & p0) << k;
            stats.het_1 += popcount(mixed & p1) << k;
        }
    }
    return stats;
}

// Scores every pair of `parents` (ascending indices), including selfing.
auto score_pairs(
    const Panel& panel,
    std::span<const uint32_t> parents,
    const CountPlanes& planes,
    const ReadWeights& weights) -> std::vector<PairStats>
{
    const size_t k = parents.size();
    std::vector<PairStats> pairs;
    pairs.reserve(k * (k + 1) / 2);
    for (size_t i = 0; i < k; ++i)
    {
        const ParentRows pi = parent_rows(panel, parents[i]);
        for (size_t j = i; j < k; ++j)
        {
            PairStats stats
                = count_pair(pi, parent_rows(panel, parents[j]), planes);
            stats.a = parents[i];
            stats.b = parents[j];
            stats.score
                = weights.score_conc * static_cast<int64_t>(stats.hom_conc)
                  + weights.score_disc * static_cast<int64_t>(stats.hom_disc);
            pairs.push_back(stats);
        }
    }
    std::ranges::sort(
        pairs,
        [](const PairStats& x, const PairStats& y)
        {
            if (x.score != y.score)
            {
                return x.score > y.score;
            }
            return x.a != y.a ? x.a < y.a : x.b < y.b;
        });
    return pairs;
}

// Per-site terms of the null and single-parent models (design §3.4), summed
// where the parent genotype does not matter and kept per site where it does.
struct SiteTerms
{
    double null_sum = 0;            // Σ ln L(random panel pair)
    double m0_sum = 0;              // Σ [ln L(parent carries allele0) - L_null]
    std::vector<uint32_t> covered;  // sites with reads
    std::vector<double> m1_minus_m0;    // per covered site
    std::vector<double> mhet_minus_m0;  // per covered site; parent het
    std::vector<uint64_t> two_reads;    // bit set: n0 + n1 >= 2
    std::vector<uint64_t> absent0;      // bit set: n0 == 0 and n1 >= 2
    std::vector<uint64_t> absent1;      // bit set: n1 == 0 and n0 >= 2
};

auto build_site_terms(
    const SiteTable& sites,
    const SampleCounts& counts,
    const ReadWeights& w,
    size_t words) -> SiteTerms
{
    SiteTerms terms;
    terms.two_reads.assign(words, 0);
    terms.absent0.assign(words, 0);
    terms.absent1.assign(words, 0);
    const std::span<const uint8_t> n0 = counts.n0();
    const std::span<const uint8_t> n1 = counts.n1();
    for (size_t s = 0; s < n0.size(); ++s)
    {
        const unsigned reads = n0[s] + n1[s];
        if (reads == 0)
        {
            continue;
        }
        terms.covered.push_back(static_cast<uint32_t>(s));
        const auto [l00, l01, l11] = site_likelihoods(n0[s], n1[s], w);
        const double q = sites.freq1(s);
        const double p = 1 - q;
        const double l_null = log_sum_exp(
            {log_or_neg_inf(p * p) + l00,
             log_or_neg_inf(2 * p * q) + l01,
             log_or_neg_inf(q * q) + l11});
        const double m0
            = log_sum_exp({log_or_neg_inf(p) + l00, log_or_neg_inf(q) + l01})
              - l_null;
        const double m1
            = log_sum_exp({log_or_neg_inf(p) + l01, log_or_neg_inf(q) + l11})
              - l_null;
        const double mhet = log_sum_exp({m0, m1}) + w.het;
        terms.null_sum += l_null;
        terms.m0_sum += m0;
        terms.m1_minus_m0.push_back(m1 - m0);
        terms.mhet_minus_m0.push_back(mhet - m0);
        if (reads >= 2)
        {
            const uint64_t bit = uint64_t{1} << (s % 64);
            terms.two_reads[s / 64] |= bit;
            if (n0[s] == 0)
            {
                terms.absent0[s / 64] |= bit;
            }
            if (n1[s] == 0)
            {
                terms.absent1[s / 64] |= bit;
            }
        }
    }
    return terms;
}

auto ratio(uint64_t num, uint64_t den) -> double
{
    return den == 0 ? nan : static_cast<double>(num) / static_cast<double>(den);
}

// Σ ln L(pair) over the covered sites where a parent is heterozygous: the
// F1 is homozygous or heterozygous with probability 1/2 each (1/4, 1/2, 1/4
// over the three genotypes when both parents are).
auto uncertain_log_likelihood(
    const ParentRows& pi,
    const ParentRows& pj,
    const CountPlanes& planes,
    const SampleCounts& counts,
    const ReadWeights& w) -> double
{
    const std::span<const uint8_t> n0 = counts.n0();
    const std::span<const uint8_t> n1 = counts.n1();
    const double ln_half = w.het;
    double sum = 0;
    for (const uint32_t word : planes.covered_words)
    {
        uint64_t uncertain = pi.het[word] | pj.het[word];
        while (uncertain != 0)
        {
            const int bit = std::countr_zero(uncertain);
            uncertain &= uncertain - 1;
            const size_t s
                = (static_cast<size_t>(word) * 64) + static_cast<size_t>(bit);
            if (n0[s] + n1[s] == 0)
            {
                continue;
            }
            const auto [l00, l01, l11] = site_likelihoods(n0[s], n1[s], w);
            const uint64_t mask = uint64_t{1} << bit;
            if ((pi.het[word] & pj.het[word] & mask) != 0)
            {
                sum += log_sum_exp(
                    {l00 + (2 * ln_half), l01 + ln_half, l11 + (2 * ln_half)});
            }
            else if (((pi.geno[word] | pj.geno[word]) & mask) != 0)
            {
                sum += log_sum_exp({l01, l11}) + ln_half;
            }
            else
            {
                sum += log_sum_exp({l00, l01}) + ln_half;
            }
        }
    }
    return sum;
}

auto make_report(
    const Panel& panel,
    const PairStats& stats,
    const CountPlanes& planes,
    const SampleCounts& counts,
    const ReadWeights& w,
    double null_sum) -> PairReport
{
    // Over the sites where both parents are homozygous, Σ_s ln L(pair)
    // collapses to the four counts because every read contributes one of
    // three constants.
    const double log_l
        = w.conc * static_cast<double>(stats.hom_conc)
          + w.disc * static_cast<double>(stats.hom_disc)
          + w.het * static_cast<double>(stats.het_0 + stats.het_1)
          + uncertain_log_likelihood(
              parent_rows(panel, stats.a),
              parent_rows(panel, stats.b),
              planes,
              counts,
              w);
    return {
        stats,
        log_l - null_sum,
        ratio(stats.hom_disc, stats.hom_conc + stats.hom_disc),
        ratio(stats.het_1, stats.het_0 + stats.het_1)};
}

auto score_singles(const Panel& panel, const SiteTerms& terms)
    -> std::vector<SingleStats>
{
    std::vector<SingleStats> singles;
    singles.reserve(panel.n_parents());
    for (size_t i = 0; i < panel.n_parents(); ++i)
    {
        const ParentRows parent = parent_rows(panel, i);
        double lod = terms.m0_sum;
        for (size_t c = 0; c < terms.covered.size(); ++c)
        {
            const uint32_t s = terms.covered[c];
            const uint64_t mask = uint64_t{1} << (s % 64);
            if ((parent.het[s / 64] & mask) != 0)
            {
                lod += terms.mhet_minus_m0[c];
            }
            else if ((parent.geno[s / 64] & mask) != 0)
            {
                lod += terms.m1_minus_m0[c];
            }
        }
        uint64_t absent = 0;
        uint64_t two_reads = 0;
        for (size_t w = 0; w < parent.geno.size(); ++w)
        {
            const uint64_t called = ~parent.het[w];
            absent += popcount(~parent.geno[w] & called & terms.absent0[w])
                      + popcount(parent.geno[w] & called & terms.absent1[w]);
            two_reads += popcount(called & terms.two_reads[w]);
        }
        singles.push_back(
            {static_cast<uint32_t>(i), lod, ratio(absent, two_reads)});
    }
    std::ranges::sort(
        singles,
        [](const SingleStats& x, const SingleStats& y)
        {
            if (x.lod_vs_null != y.lod_vs_null)
            {
                return x.lod_vs_null > y.lod_vs_null;
            }
            return x.parent < y.parent;
        });
    return singles;
}

// The `n` best single parents by index (every parent when n == 0).
auto candidate_parents(std::span<const SingleStats> singles, size_t n)
    -> std::vector<uint32_t>
{
    if (n == 0 || n > singles.size())
    {
        n = singles.size();
    }
    std::vector<uint32_t> parents;
    parents.reserve(n);
    for (const SingleStats& single : singles.first(n))
    {
        parents.push_back(single.parent);
    }
    std::ranges::sort(parents);
    return parents;
}
}  // namespace

auto score_sample(
    const Panel& panel,
    const SampleCounts& counts,
    const ScoreOptions& options) -> SampleScores
{
    if (counts.n_sites() != panel.n_sites())
    {
        throw std::logic_error(
            std::format(
                "sample {} has {} sites but the panel has {}",
                counts.sample(),
                counts.n_sites(),
                panel.n_sites()));
    }
    const ReadWeights weights = make_weights(options.error_rate);
    const size_t words = panel.geno().words_per_row();
    const CountPlanes planes = build_planes(counts, words);
    const SiteTerms terms
        = build_site_terms(panel.sites(), counts, weights, words);

    SampleScores scores;
    scores.n_covered = terms.covered.size();
    scores.singles = score_singles(panel, terms);
    scores.pairs = score_pairs(
        panel,
        candidate_parents(scores.singles, options.candidates),
        planes,
        weights);
    const size_t n_top = std::min(options.top, scores.pairs.size());
    scores.top.reserve(n_top);
    for (size_t r = 0; r < n_top; ++r)
    {
        scores.top.push_back(make_report(
            panel, scores.pairs[r], planes, counts, weights, terms.null_sum));
    }
    scores.lod_vs_second
        = scores.pairs.size() < 2
              ? nan
              : static_cast<double>(
                    scores.pairs[0].score - scores.pairs[1].score)
                    / 1000;
    return scores;
}
}  // namespace vcfbox::parentage
