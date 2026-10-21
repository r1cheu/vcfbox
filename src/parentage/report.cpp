#include "parentage/report.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/lines.h"
#include "common/strings.h"
#include "parentage/counts.h"
#include "parentage/score.h"

namespace vcfbox::parentage
{
namespace
{
constexpr std::string_view na = "NA";
constexpr std::string_view dropped = "dropped";
constexpr std::string_view absent = "ABSENT";
constexpr std::string_view rank_first = "1";
constexpr std::string_view rank_lower = "2+";
// Classes of expected_pair_rank in the summary.
constexpr std::array<std::string_view, 4> rank_classes{
    rank_first, rank_lower, absent, na};

// `value` with `precision` decimals, or NA when it is NaN.
auto fixed(double value, int precision) -> std::string
{
    return std::isnan(value) ? std::string(na)
                             : std::format("{:.{}f}", value, precision);
}

auto open_output(const std::string& path) -> std::ofstream
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error(std::format("cannot create {}", path));
    }
    return out;
}

enum class Model
{
    none,
    single,
    pair,
};

// The model with the highest lod_vs_null; ties go to the pair, then the
// single parent.
auto leading_model(const SampleScores& scores) -> Model
{
    const double pair = scores.top.front().lod_vs_null;
    const double single = scores.singles.front().lod_vs_null;
    if (pair >= single && pair >= 0)
    {
        return Model::pair;
    }
    return single >= 0 ? Model::single : Model::none;
}

constexpr std::string_view table_header
    = "sample\tcoverage\tbest\tbest_a\tbest_b\texpected_a\texpected_b\t"
      "null\tsingle\tbest_pair\tsecond_pair\texpected\t"
      "pair_hom_mismatch_rate\tpair_het_alt_fraction\tpair_fit\t"
      "single_allele_absent_rate\tcontest\texpected_pair_rank";

// Columns after `expected_b`.
constexpr size_t n_scored_columns = 11;

auto mean_depth(const SampleCounts& counts) -> double
{
    return static_cast<double>(counts.bases_used())
           / static_cast<double>(counts.n_sites());
}

// The parents of the leading model: the best pair, the best single parent
// and NA, or NA twice for the null.
auto leading_parents(
    const SampleScores& scores,
    std::span<const std::string> names) -> std::pair<std::string, std::string>
{
    switch (leading_model(scores))
    {
        case Model::pair:
            return {
                names[scores.top.front().stats.a],
                names[scores.top.front().stats.b]};
        case Model::single:
            return {names[scores.singles.front().parent], std::string(na)};
        case Model::none:
            return {std::string(na), std::string(na)};
    }
    throw std::logic_error("unknown model");
}

// The expected pair's rank among the scored pairs, ABSENT when a parent is
// not in the panel, NA without an expected pair.
auto expected_rank(
    const SampleScores& scores,
    const ExpectedParents* expected) -> std::string
{
    if (expected == nullptr)
    {
        return std::string(na);
    }
    return scores.expected ? std::to_string(scores.expected->rank)
                           : std::string(absent);
}

auto rank_class(const SampleScores& scores, const ExpectedParents* expected)
    -> std::string_view
{
    if (expected == nullptr)
    {
        return na;
    }
    if (!scores.expected)
    {
        return absent;
    }
    return scores.expected->rank == 1 ? rank_first : rank_lower;
}

// `sample<TAB>coverage<TAB>best<TAB>best_a<TAB>best_b<TAB>expected_a<TAB>expected_b`
auto write_leading_columns(
    std::ostream& out,
    const SampleCounts& counts,
    std::string_view best,
    const std::pair<std::string, std::string>& parents,
    const ExpectedParents* expected) -> void
{
    out << std::format(
        "{}\t{:.3f}\t{}\t{}\t{}\t{}\t{}",
        counts.sample(),
        mean_depth(counts),
        best,
        parents.first,
        parents.second,
        expected == nullptr ? na : std::string_view(expected->a),
        expected == nullptr ? na : std::string_view(expected->b));
}
}  // namespace

auto to_string(Status status) -> std::string_view
{
    switch (status)
    {
        case Status::null:
            return "null";
        case Status::single:
            return "single";
        case Status::pair:
            return "pair";
        case Status::ambiguous:
            return "ambiguous";
    }
    throw std::logic_error("unknown status");
}

auto to_string(Fit fit) -> std::string_view
{
    switch (fit)
    {
        case Fit::pass:
            return "PASS";
        case Fit::fail:
            return "FAIL";
        case Fit::na:
            return na;
    }
    throw std::logic_error("unknown fit");
}

auto model_margin(const SampleScores& scores) -> double
{
    const double pair = scores.top.front().lod_vs_null;
    const double single = scores.singles.front().lod_vs_null;
    switch (leading_model(scores))
    {
        case Model::pair:
            return pair - std::max(single, 0.0);
        case Model::single:
            return single - std::max(pair, 0.0);
        case Model::none:
            return -std::max(pair, single);
    }
    throw std::logic_error("unknown model");
}

auto single_lod_vs_second(const SampleScores& scores) -> double
{
    if (scores.singles.size() < 2)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return scores.singles[0].lod_vs_null - scores.singles[1].lod_vs_null;
}

auto to_string(Contest contest) -> std::string_view
{
    switch (contest)
    {
        case Contest::model:
            return "model";
        case Contest::pair:
            return "pair";
        case Contest::single:
            return "single";
    }
    throw std::logic_error("unknown contest");
}

auto find_ambiguity(const SampleScores& scores, const CallOptions& options)
    -> std::optional<Ambiguity>
{
    const double between = model_margin(scores);
    if (between < options.min_lod)
    {
        return Ambiguity{Contest::model, between};
    }
    const Model leader = leading_model(scores);
    if (leader == Model::pair && scores.lod_vs_second < options.min_lod)
    {
        return Ambiguity{Contest::pair, scores.lod_vs_second};
    }
    const double single = single_lod_vs_second(scores);
    if (leader == Model::single && single < options.min_lod)
    {
        return Ambiguity{Contest::single, single};
    }
    return std::nullopt;
}

auto assign_status(const SampleScores& scores, const CallOptions& options)
    -> Status
{
    if (find_ambiguity(scores, options))
    {
        return Status::ambiguous;
    }
    switch (leading_model(scores))
    {
        case Model::pair:
            return Status::pair;
        case Model::single:
            return Status::single;
        case Model::none:
            return Status::null;
    }
    throw std::logic_error("unknown model");
}

auto assign_fit(const SampleScores& scores, const CallOptions& options)
    -> Fit
{
    const double rate = scores.top.front().hom_mismatch_rate;
    if (std::isnan(rate))
    {
        return Fit::na;
    }
    return rate > options.max_hom_mismatch ? Fit::fail : Fit::pass;
}

auto read_expected_pairs(
    const std::string& path,
    std::span<const std::string> parent_names)
    -> std::unordered_map<std::string, ExpectedParents>
{
    std::unordered_map<std::string_view, uint32_t> index;
    for (size_t i = 0; i < parent_names.size(); ++i)
    {
        index.emplace(parent_names[i], static_cast<uint32_t>(i));
    }
    std::unordered_map<std::string, ExpectedParents> expected;
    size_t line_no = 0;
    for (const auto& line : read_lines(path))
    {
        ++line_no;
        if (line.empty() || line[0] == '#'
            || (line_no == 1 && line.starts_with("sample")))
        {
            continue;
        }
        const auto fields = split(line, '\t');
        if (fields.size() != 3)
        {
            throw std::runtime_error(
                std::format(
                    "{}:{}: expected `sample<TAB>mother<TAB>father`, found {} "
                    "fields",
                    path,
                    line_no,
                    fields.size()));
        }
        ExpectedParents parents{
            std::string(fields[1]), std::string(fields[2]), std::nullopt};
        const auto a = index.find(fields[1]);
        const auto b = index.find(fields[2]);
        if (a != index.end() && b != index.end())
        {
            parents.pair = ParentPair{a->second, b->second};
        }
        if (!expected.emplace(std::string(fields[0]), std::move(parents))
                 .second)
        {
            throw std::runtime_error(
                std::format(
                    "{}:{}: sample {} is listed twice",
                    path,
                    line_no,
                    fields[0]));
        }
    }
    return expected;
}

ReportWriter::ReportWriter(
    std::string prefix,
    std::span<const std::string> parent_names,
    bool with_expected,
    const CallOptions& options,
    std::ofstream table,
    std::ofstream topk)
    : prefix_(std::move(prefix)),
      names_(parent_names),
      with_expected_(with_expected),
      options_(options),
      table_(std::move(table)),
      topk_(std::move(topk))
{
}

auto open_report(
    const std::string& prefix,
    std::span<const std::string> parent_names,
    bool with_expected,
    const CallOptions& options) -> ReportWriter
{
    std::ofstream table = open_output(std::format("{}.tsv", prefix));
    table << table_header << '\n' << std::flush;
    std::ofstream topk = open_output(std::format("{}.topk.tsv", prefix));
    topk << "sample\trank\ta\tb\thom_conc\thom_disc\thet_0\thet_1\t"
            "lod_vs_null\thom_mismatch_rate\thet_alt_fraction\n"
         << std::flush;
    return {
        prefix,
        parent_names,
        with_expected,
        options,
        std::move(table),
        std::move(topk)};
}

auto ReportWriter::add(
    const SampleCounts& counts,
    const SampleScores& scores,
    const ExpectedParents* expected) -> void
{
    const PairReport& best = scores.top.front();
    const SingleStats& single = scores.singles.front();
    const Status call = assign_status(scores, options_);
    const Fit fit = assign_fit(scores, options_);
    const std::optional<Ambiguity> ambiguity
        = find_ambiguity(scores, options_);
    const double null = scores.null_log_likelihood;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    write_leading_columns(
        table_,
        counts,
        to_string(call),
        leading_parents(scores, names_),
        expected);
    table_ << std::format(
        "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        fixed(null, 1),
        fixed(null + single.lod_vs_null, 1),
        fixed(null + best.lod_vs_null, 1),
        fixed(scores.top.size() < 2 ? nan : null + scores.top[1].lod_vs_null, 1),
        fixed(
            scores.expected ? null + scores.expected->report.lod_vs_null : nan,
            1),
        fixed(best.hom_mismatch_rate, 5),
        fixed(best.het_alt_fraction, 4),
        to_string(fit),
        fixed(single.allele_absent_rate, 5),
        ambiguity ? to_string(ambiguity->contest) : na,
        expected_rank(scores, expected));
    table_ << std::flush;

    for (size_t rank = 0; rank < scores.top.size(); ++rank)
    {
        const PairReport& report = scores.top[rank];
        const PairStats& st = report.stats;
        topk_ << std::format(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            counts.sample(),
            rank + 1,
            names_[st.a],
            names_[st.b],
            st.hom_conc,
            st.hom_disc,
            st.het_0,
            st.het_1,
            fixed(report.lod_vs_null, 1),
            fixed(report.hom_mismatch_rate, 5),
            fixed(report.het_alt_fraction, 4));
    }
    topk_ << std::flush;

    ++tally_[{to_string(call), rank_class(scores, expected)}];
    if (call == Status::pair)
    {
        ++fits_[to_string(fit)];
    }
    check_written();
}

auto ReportWriter::add_dropped(
    const SampleCounts& counts,
    const ExpectedParents* expected) -> void
{
    write_leading_columns(
        table_, counts, dropped, {std::string(na), std::string(na)}, expected);
    for (size_t i = 0; i < n_scored_columns; ++i)
    {
        table_ << '\t' << na;
    }
    table_ << '\n' << std::flush;
    ++tally_[{dropped, na}];
    check_written();
}

auto ReportWriter::check_written() -> void
{
    if (!table_ || !topk_)
    {
        throw std::runtime_error(
            std::format("cannot write the report of {}", prefix_));
    }
}

auto ReportWriter::finish() -> void
{
    std::map<std::string_view, size_t> by_best;
    size_t n_samples = 0;
    for (const auto& [key, n] : tally_)
    {
        by_best[key.first] += n;
        n_samples += n;
    }
    std::vector<std::string_view> bests;
    for (const Status call : all_statuses)
    {
        bests.push_back(to_string(call));
    }
    bests.push_back(dropped);

    const std::string path = std::format("{}.summary.txt", prefix_);
    std::ofstream out = open_output(path);
    out << std::format("samples\t{}\n\nbest\tsamples\n", n_samples);
    for (const std::string_view best : bests)
    {
        out << std::format("{}\t{}\n", best, by_best[best]);
    }
    out << "\npair_fit\tsamples\n";
    for (const Fit fit : all_fits)
    {
        out << std::format("{}\t{}\n", to_string(fit), fits_[to_string(fit)]);
    }
    if (with_expected_)
    {
        out << "\nbest\\expected_pair_rank";
        for (const std::string_view rank : rank_classes)
        {
            out << '\t' << rank;
        }
        out << '\n';
        for (const std::string_view best : bests)
        {
            out << best;
            for (const std::string_view rank : rank_classes)
            {
                out << '\t' << tally_[{best, rank}];
            }
            out << '\n';
        }
    }
    out.flush();
    if (!out)
    {
        throw std::runtime_error(std::format("cannot write {}", path));
    }
}
}  // namespace vcfbox::parentage
