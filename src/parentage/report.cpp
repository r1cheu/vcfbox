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

#include "common/lines.h"
#include "common/strings.h"
#include "parentage/counts.h"
#include "parentage/score.h"

namespace vcfbox::parentage
{
namespace
{
constexpr std::string_view na = "NA";
constexpr std::string_view match = "MATCH";
constexpr std::string_view mismatch = "MISMATCH";
constexpr std::string_view absent = "ABSENT";
constexpr std::array<std::string_view, 3> expected_statuses{
    match, mismatch, absent};

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

// MATCH when the expected pair ranks first, MISMATCH otherwise, ABSENT when
// a parent is not in the panel.
auto expected_status(const SampleScores& scores) -> std::string_view
{
    if (!scores.expected)
    {
        return absent;
    }
    return scores.expected->rank == 1 ? match : mismatch;
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

// The model with the highest lod_vs_null after the leader.
auto runner_up_model(const SampleScores& scores) -> Model
{
    const double pair = scores.top.front().lod_vs_null;
    const double single = scores.singles.front().lod_vs_null;
    switch (leading_model(scores))
    {
        case Model::pair:
            return single >= 0 ? Model::single : Model::none;
        case Model::single:
            return pair >= 0 ? Model::pair : Model::none;
        case Model::none:
            return pair >= single ? Model::pair : Model::single;
    }
    throw std::logic_error("unknown model");
}

// `a/b` for a pair, the name for a single parent.
auto pair_name(const PairStats& pair, std::span<const std::string> names)
    -> std::string
{
    return std::format("{}/{}", names[pair.a], names[pair.b]);
}

// The best pair, the best single parent, or NONE for the null model.
auto model_name(
    const SampleScores& scores,
    Model model,
    std::span<const std::string> names) -> std::string
{
    switch (model)
    {
        case Model::pair:
            return pair_name(scores.top.front().stats, names);
        case Model::single:
            return names[scores.singles.front().parent];
        case Model::none:
            return "NONE";
    }
    throw std::logic_error("unknown model");
}

// The two candidates `ambiguity` could not separate, best first.
auto contenders(
    const SampleScores& scores,
    const Ambiguity& ambiguity,
    std::span<const std::string> names) -> std::pair<std::string, std::string>
{
    switch (ambiguity.contest)
    {
        case Contest::model:
            return {
                model_name(scores, leading_model(scores), names),
                model_name(scores, runner_up_model(scores), names)};
        case Contest::pair:
            return {
                pair_name(scores.top[0].stats, names),
                pair_name(scores.top[1].stats, names)};
        case Contest::single:
            return {
                names[scores.singles[0].parent],
                names[scores.singles[1].parent]};
    }
    throw std::logic_error("unknown contest");
}

constexpr std::array<std::string_view, all_statuses.size()> table_headers{
    "sample\tdepth_mean",
    "sample\tdepth_mean\tparent",
    "sample\tdepth_mean\ta\tb\tis_selfing\thom_mismatch_rate\t"
    "het_alt_fraction\tfit",
    "sample\tdepth_mean\tcontest\tfirst\tsecond\tmargin\tfit",
};

constexpr std::array<std::string_view, all_statuses.size()> table_names{
    "zero_parent", "one_parent", "two_parent", "ambiguous"};
}  // namespace

auto to_string(Status status) -> std::string_view
{
    switch (status)
    {
        case Status::zero_parent:
            return "ZERO_PARENT";
        case Status::one_parent:
            return "ONE_PARENT";
        case Status::two_parent:
            return "TWO_PARENT";
        case Status::ambiguous:
            return "AMBIGUOUS";
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
            return Status::two_parent;
        case Model::single:
            return Status::one_parent;
        case Model::none:
            return Status::zero_parent;
    }
    throw std::logic_error("unknown model");
}

auto assign_fit(const SampleScores& scores, const CallOptions& options)
    -> Fit
{
    const double rate = scores.top.front().hom_mismatch_rate;
    if (leading_model(scores) != Model::pair || std::isnan(rate))
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
    std::array<std::ofstream, all_statuses.size()> tables,
    std::ofstream topk)
    : prefix_(std::move(prefix)),
      names_(parent_names),
      with_expected_(with_expected),
      options_(options),
      tables_(std::move(tables)),
      topk_(std::move(topk))
{
}

auto open_report(
    const std::string& prefix,
    std::span<const std::string> parent_names,
    bool with_expected,
    const CallOptions& options) -> ReportWriter
{
    std::array<std::ofstream, all_statuses.size()> tables;
    for (size_t i = 0; i < tables.size(); ++i)
    {
        std::ofstream& table = tables[i];
        table = open_output(std::format("{}.{}.tsv", prefix, table_names[i]));
        table << table_headers[i];
        if (with_expected)
        {
            table << "\texpected_a\texpected_b\texpected_status\t"
                     "expected_lod_vs_best";
        }
        table << '\n' << std::flush;
    }
    std::ofstream topk = open_output(std::format("{}.topk.tsv", prefix));
    topk << "sample\trank\ta\tb\thom_conc\thom_disc\thet_0\thet_1\t"
            "lod_vs_null\thom_mismatch_rate\thet_alt_fraction\n"
         << std::flush;
    return {
        prefix,
        parent_names,
        with_expected,
        options,
        std::move(tables),
        std::move(topk)};
}

auto ReportWriter::write_expected(
    std::ostream& out,
    const SampleScores& scores,
    const ExpectedParents* expected) -> std::string_view
{
    if (!with_expected_)
    {
        return na;
    }
    if (expected == nullptr)
    {
        out << std::format("\t{0}\t{0}\t{0}\t{0}", na);
        return na;
    }
    const std::string_view status = expected_status(scores);
    const double lod_vs_best = scores.expected
                                   ? scores.expected->report.lod_vs_null
                                         - scores.top.front().lod_vs_null
                                   : std::numeric_limits<double>::quiet_NaN();
    out << std::format(
        "\t{}\t{}\t{}\t{}",
        expected->a,
        expected->b,
        status,
        fixed(lod_vs_best, 1));
    return status;
}

auto ReportWriter::add(
    const SampleCounts& counts,
    const SampleScores& scores,
    const ExpectedParents* expected) -> void
{
    const std::string& sample = counts.sample();
    const PairReport& best = scores.top.front();
    const Status call = assign_status(scores, options_);
    const Fit fit = assign_fit(scores, options_);
    std::ofstream& table = tables_[static_cast<size_t>(call)];

    table << std::format(
        "{}\t{:.3f}",
        sample,
        static_cast<double>(counts.bases_used())
            / static_cast<double>(counts.n_sites()));
    switch (call)
    {
        case Status::zero_parent:
            break;
        case Status::one_parent:
            table << std::format(
                "\t{}", names_[scores.singles.front().parent]);
            break;
        case Status::two_parent:
            table << std::format(
                "\t{}\t{}\t{}\t{}\t{}\t{}",
                names_[best.stats.a],
                names_[best.stats.b],
                best.stats.a == best.stats.b ? 1 : 0,
                fixed(best.hom_mismatch_rate, 5),
                fixed(best.het_alt_fraction, 4),
                to_string(fit));
            break;
        case Status::ambiguous:
        {
            const Ambiguity ambiguity = *find_ambiguity(scores, options_);
            const auto [first, second] = contenders(scores, ambiguity, names_);
            table << std::format(
                "\t{}\t{}\t{}\t{}\t{}",
                to_string(ambiguity.contest),
                first,
                second,
                fixed(ambiguity.margin, 1),
                to_string(fit));
            break;
        }
    }
    const std::string_view status = write_expected(table, scores, expected);
    table << '\n' << std::flush;

    for (size_t rank = 0; rank < scores.top.size(); ++rank)
    {
        const PairReport& report = scores.top[rank];
        const PairStats& st = report.stats;
        topk_ << std::format(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            sample,
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

    ++n_samples_;
    ++tally_[{to_string(call), status}];
    if (call == Status::two_parent)
    {
        ++fits_[to_string(fit)];
    }
    if (!table || !topk_)
    {
        throw std::runtime_error(
            std::format("cannot write the report of {}", prefix_));
    }
}

auto ReportWriter::finish() -> void
{
    std::map<std::string_view, size_t> by_call;
    std::map<std::string_view, size_t> by_expected;
    for (const auto& [key, n] : tally_)
    {
        by_call[key.first] += n;
        by_expected[key.second] += n;
    }
    const std::string path = std::format("{}.summary.txt", prefix_);
    std::ofstream out = open_output(path);
    out << std::format(
        "samples\t{}\ndropped\t{}\n\nstatus\tsamples\n",
        n_samples_,
        n_dropped_);
    for (const Status call : all_statuses)
    {
        out << std::format(
            "{}\t{}\n", to_string(call), by_call[to_string(call)]);
    }
    out << "\ntwo_parent_fit\tsamples\n";
    for (const Fit fit : all_fits)
    {
        out << std::format("{}\t{}\n", to_string(fit), fits_[to_string(fit)]);
    }
    if (with_expected_)
    {
        out << "\nexpected_status\tsamples\n";
        for (const std::string_view status : expected_statuses)
        {
            out << std::format("{}\t{}\n", status, by_expected[status]);
        }
        out << std::format("{}\t{}\n", na, by_expected[na]);
        out << "\nstatus";
        for (const std::string_view status : expected_statuses)
        {
            out << '\t' << status;
        }
        out << '\t' << na << '\n';
        for (const Status call : all_statuses)
        {
            out << to_string(call);
            for (const std::string_view status : expected_statuses)
            {
                out << '\t' << tally_[{to_string(call), status}];
            }
            out << '\t' << tally_[{to_string(call), na}] << '\n';
        }
    }
    out.flush();
    if (!out)
    {
        throw std::runtime_error(std::format("cannot write {}", path));
    }
}
}  // namespace vcfbox::parentage
