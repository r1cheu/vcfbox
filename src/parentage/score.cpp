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
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
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
};

auto make_weights(double error_rate) -> ReadWeights
{
    if (error_rate <= 0 || error_rate >= 0.5)
    {
        throw std::runtime_error(
            std::format("error rate must be in (0, 0.5), got {}", error_rate));
    }
    return {std::log1p(-error_rate), std::log(error_rate), std::log(0.5)};
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

// ln L(null) of one site: both gametes drawn at the panel frequency.
auto null_log_likelihood(const SiteLikelihoods& l, double q) -> double
{
    const double p = 1 - q;
    return log_sum_exp(
        {log_or_neg_inf(p * p) + l.l00,
         log_or_neg_inf(2 * p * q) + l.l01,
         log_or_neg_inf(q * q) + l.l11});
}

// Maps a site's lod vs the null, x, to ln((1 - γ) e^x + γ): with
// probability γ the F1 genotype is drawn as under the null instead of from
// the model (the random genotype replacement model of Kalinowski, Taper &
// Marshall 2007, Mol Ecol 16:1099).
class SiteMixture
{
   public:
    explicit SiteMixture(double site_error_rate)
        : gamma_(site_error_rate),
          log_keep_(std::log1p(-site_error_rate)),
          log_replace_(std::log(site_error_rate))
    {
        if (site_error_rate < 0 || site_error_rate >= 1)
        {
            throw std::runtime_error(
                std::format(
                    "site error rate must be in [0, 1), got {}",
                    site_error_rate));
        }
    }

    auto operator()(double x) const -> double
    {
        return gamma_ == 0 ? x : log_sum_exp({log_keep_ + x, log_replace_});
    }

   private:
    double gamma_;
    double log_keep_;
    double log_replace_;
};

// Bit k of every n0 / n1 count, laid out like a panel row, so that a count
// weighted by a parent-pair mask is a popcount per plane.
struct CountPlanes
{
    std::array<std::vector<uint64_t>, n_planes> p0;
    std::array<std::vector<uint64_t>, n_planes> p1;
    std::vector<uint32_t> covered_words;  // words holding a covered site
    std::vector<uint64_t> covered;  // covered sites of each covered word
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
            planes.covered.push_back(0);
            last_covered_word = word;
        }
        const uint64_t bit = uint64_t{1} << (s % 64);
        planes.covered.back() |= bit;
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

// The four read splits of one pair over the sites where both parents are
// homozygous. Padding bits past the last site are set in `both0`, which is
// harmless because the planes are zero there.
auto count_pair(
    const ParentRows& pi,
    const ParentRows& pj,
    const CountPlanes& planes) -> PairStats
{
    PairStats stats;
    for (const uint32_t w : planes.covered_words)
    {
        const uint64_t called = both_homozygous(pi, pj, w);
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

// Per-site lods vs the null after the site mixture (design §3.4): of the
// single-parent model, summed where the parent genotype does not matter and
// kept per site where it does, and of a fixed F1 genotype at every site.
struct SiteTerms
{
    double null_sum = 0;            // Σ ln L(random panel pair)
    double m0_sum = 0;              // Σ lod(parent carries allele0)
    // lod of F1 genotype 00 / 01 / 11 at site s at [3 s] / [3 s + 1] /
    // [3 s + 2]; 0 where the site has no reads.
    std::vector<float> fixed;
    std::vector<uint32_t> covered;  // sites with reads
    std::vector<double> m0;             // per covered site
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
    const SiteMixture& mix,
    size_t words) -> SiteTerms
{
    SiteTerms terms;
    terms.fixed.assign(3 * sites.size(), 0.0F);
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
        const SiteLikelihoods l = site_likelihoods(n0[s], n1[s], w);
        const double q = sites.freq1(s);
        const double p = 1 - q;
        const double l_null = null_log_likelihood(l, q);
        terms.fixed[3 * s] = static_cast<float>(mix(l.l00 - l_null));
        terms.fixed[(3 * s) + 1] = static_cast<float>(mix(l.l01 - l_null));
        terms.fixed[(3 * s) + 2] = static_cast<float>(mix(l.l11 - l_null));
        const double raw0 = log_sum_exp(
                                {log_or_neg_inf(p) + l.l00,
                                 log_or_neg_inf(q) + l.l01})
                            - l_null;
        const double raw1 = log_sum_exp(
                                {log_or_neg_inf(p) + l.l01,
                                 log_or_neg_inf(q) + l.l11})
                            - l_null;
        const double m0 = mix(raw0);
        const double m1 = mix(raw1);
        const double mhet = mix(log_sum_exp({raw0, raw1}) + w.het);
        terms.null_sum += l_null;
        terms.m0_sum += m0;
        terms.m0.push_back(m0);
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

// Probability that a parent passes on allele1 at bit `mask` of word `w`: 0
// or 1 when homozygous, 1/2 when heterozygous, the panel frequency `freq1`
// when missing.
auto gamete_allele1(
    const ParentRows& parent,
    size_t w,
    uint64_t mask,
    double freq1) -> double
{
    if ((parent.missing[w] & mask) != 0)
    {
        return freq1;
    }
    if ((parent.het[w] & mask) != 0)
    {
        return 0.5;
    }
    return (parent.geno[w] & mask) != 0 ? 1 : 0;
}

// Σ lod(pair) over the covered sites where both parents are homozygous, each
// site taking the lod of the F1 genotype they fix.
auto homozygous_log_odds(
    const ParentRows& pi,
    const ParentRows& pj,
    const CountPlanes& planes,
    const SiteTerms& terms) -> double
{
    double sum = 0;
    const auto add = [&](uint32_t word, uint64_t sites, size_t genotype)
    {
        while (sites != 0)
        {
            const auto bit = static_cast<size_t>(std::countr_zero(sites));
            sites &= sites - 1;
            const size_t s = (static_cast<size_t>(word) * 64) + bit;
            sum += terms.fixed[(3 * s) + genotype];
        }
    };
    for (const auto [i, w] : std::views::enumerate(planes.covered_words))
    {
        const uint64_t called = both_homozygous(pi, pj, w)
                                & planes.covered[static_cast<size_t>(i)];
        add(w, ~(pi.geno[w] | pj.geno[w]) & called, 0);
        add(w, (pi.geno[w] ^ pj.geno[w]) & called, 1);
        add(w, pi.geno[w] & pj.geno[w] & called, 2);
    }
    return sum;
}

// Σ lod(pair) over the covered sites where a parent is heterozygous or
// missing: the F1 genotype is marginalised over the two gametes.
auto uncertain_log_odds(
    const ParentRows& pi,
    const ParentRows& pj,
    const CountPlanes& planes,
    const SampleCounts& counts,
    const SiteTable& sites,
    const ReadWeights& w,
    const SiteMixture& mix) -> double
{
    const std::span<const uint8_t> n0 = counts.n0();
    const std::span<const uint8_t> n1 = counts.n1();
    double sum = 0;
    for (const auto [i, word] : std::views::enumerate(planes.covered_words))
    {
        uint64_t uncertain = ~both_homozygous(pi, pj, word)
                             & planes.covered[static_cast<size_t>(i)];
        while (uncertain != 0)
        {
            const int bit = std::countr_zero(uncertain);
            uncertain &= uncertain - 1;
            const size_t s
                = (static_cast<size_t>(word) * 64) + static_cast<size_t>(bit);
            const SiteLikelihoods l = site_likelihoods(n0[s], n1[s], w);
            const uint64_t mask = uint64_t{1} << bit;
            const double q = sites.freq1(s);
            const double gi = gamete_allele1(pi, word, mask, q);
            const double gj = gamete_allele1(pj, word, mask, q);
            const double log_l = log_sum_exp(
                {log_or_neg_inf((1 - gi) * (1 - gj)) + l.l00,
                 log_or_neg_inf((gi * (1 - gj)) + ((1 - gi) * gj)) + l.l01,
                 log_or_neg_inf(gi * gj) + l.l11});
            sum += mix(log_l - null_log_likelihood(l, q));
        }
    }
    return sum;
}

auto score_pair(
    const Panel& panel,
    uint32_t a,
    uint32_t b,
    const CountPlanes& planes,
    const SampleCounts& counts,
    const SiteTerms& terms,
    const ReadWeights& w,
    const SiteMixture& mix) -> PairReport
{
    const ParentRows pa = panel.parent_rows(a);
    const ParentRows pb = panel.parent_rows(b);
    PairStats stats = count_pair(pa, pb, planes);
    stats.a = a;
    stats.b = b;
    const double lod
        = homozygous_log_odds(pa, pb, planes, terms)
          + uncertain_log_odds(pa, pb, planes, counts, panel.sites(), w, mix);
    return {
        stats,
        lod,
        ratio(stats.hom_disc, stats.hom_conc + stats.hom_disc),
        ratio(stats.het_1, stats.het_0 + stats.het_1)};
}

// Scores every pair of `parents` (ascending indices), including selfing;
// best lod_vs_null first, ties broken by (a, b).
auto score_pairs(
    const Panel& panel,
    std::span<const uint32_t> parents,
    const CountPlanes& planes,
    const SampleCounts& counts,
    const SiteTerms& terms,
    const ReadWeights& w,
    const SiteMixture& mix) -> std::vector<PairReport>
{
    const size_t k = parents.size();
    std::vector<PairReport> pairs;
    pairs.reserve(k * (k + 1) / 2);
    for (size_t i = 0; i < k; ++i)
    {
        for (size_t j = i; j < k; ++j)
        {
            pairs.push_back(score_pair(
                panel,
                parents[i],
                parents[j],
                planes,
                counts,
                terms,
                w,
                mix));
        }
    }
    std::ranges::sort(
        pairs,
        [](const PairReport& x, const PairReport& y)
        {
            if (x.lod_vs_null != y.lod_vs_null)
            {
                return x.lod_vs_null > y.lod_vs_null;
            }
            return x.stats.a != y.stats.a ? x.stats.a < y.stats.a
                                          : x.stats.b < y.stats.b;
        });
    return pairs;
}

auto score_singles(const Panel& panel, const SiteTerms& terms)
    -> std::vector<SingleStats>
{
    std::vector<SingleStats> singles;
    singles.reserve(panel.n_parents());
    for (size_t i = 0; i < panel.n_parents(); ++i)
    {
        const ParentRows parent = panel.parent_rows(i);
        double lod = terms.m0_sum;
        for (size_t c = 0; c < terms.covered.size(); ++c)
        {
            const uint32_t s = terms.covered[c];
            const uint64_t mask = uint64_t{1} << (s % 64);
            if ((parent.missing[s / 64] & mask) != 0)
            {
                lod -= terms.m0[c];
            }
            else if ((parent.het[s / 64] & mask) != 0)
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
            const uint64_t called = ~(parent.het[w] | parent.missing[w]);
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
    const ScoreOptions& options,
    std::optional<ParentPair> expected) -> SampleScores
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
    const SiteMixture mix(options.site_error_rate);
    const size_t words = panel.geno().words_per_row();
    const CountPlanes planes = build_planes(counts, words);
    const SiteTerms terms
        = build_site_terms(panel.sites(), counts, weights, mix, words);

    SampleScores scores;
    scores.n_covered = terms.covered.size();
    scores.null_log_likelihood = terms.null_sum;
    scores.singles = score_singles(panel, terms);
    std::vector<uint32_t> parents
        = candidate_parents(scores.singles, options.candidates);
    if (expected)
    {
        parents.push_back(expected->a);
        parents.push_back(expected->b);
        std::ranges::sort(parents);
        const auto [first, last] = std::ranges::unique(parents);
        parents.erase(first, last);
    }
    std::vector<PairReport> pairs = score_pairs(
        panel, parents, planes, counts, terms, weights, mix);
    scores.lod_vs_second = pairs.size() < 2
                               ? nan
                               : pairs[0].lod_vs_null - pairs[1].lod_vs_null;
    if (expected)
    {
        const uint32_t a = std::min(expected->a, expected->b);
        const uint32_t b = std::max(expected->a, expected->b);
        const auto found = std::ranges::find_if(
            pairs,
            [&](const PairReport& p)
            { return p.stats.a == a && p.stats.b == b; });
        scores.expected = ExpectedScore{
            *found, static_cast<size_t>(found - pairs.begin()) + 1};
    }
    pairs.resize(std::min(options.top, pairs.size()));
    scores.top = std::move(pairs);
    return scores;
}
}  // namespace vcfbox::parentage
