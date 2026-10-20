#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "parentage/counts.h"
#include "parentage/score.h"

namespace vcfbox::parentage
{
struct CallOptions
{
    double min_lod = 50;             // every margin a call needs
    double max_hom_mismatch = 0.03;  // best-pair hom_mismatch_rate for fit
};

// How many panel parents the sample's reads identify.
enum class Status
{
    zero_parent,
    one_parent,
    two_parent,
    ambiguous,
};

inline constexpr std::array<Status, 4> all_statuses{
    Status::zero_parent,
    Status::one_parent,
    Status::two_parent,
    Status::ambiguous};

// Whether the best pair explains the reads; `na` when the pair model does
// not lead or its hom_mismatch_rate is NaN.
enum class Fit
{
    pass,
    fail,
    na,
};

inline constexpr std::array<Fit, 3> all_fits{Fit::pass, Fit::fail, Fit::na};

[[nodiscard]] auto to_string(Status status) -> std::string_view;
[[nodiscard]] auto to_string(Fit fit) -> std::string_view;

// The lod_vs_null of the leading model minus the best of the other two,
// over the models: no panel parent (0), the best single parent and the best
// pair.
[[nodiscard]] auto model_margin(const SampleScores& scores) -> double;

// The best single parent's lod_vs_null minus the second's; NaN with a single
// parent.
[[nodiscard]] auto single_lod_vs_second(const SampleScores& scores) -> double;

// The margin an AMBIGUOUS call fails: between the leading model and the
// runner-up model, or within the leading model between the best pair and the
// second (pair) or the best single parent and the second (single).
enum class Contest
{
    model,
    pair,
    single,
};

[[nodiscard]] auto to_string(Contest contest) -> std::string_view;

struct Ambiguity
{
    Contest contest;
    double margin;
};

// The first margin below min_lod: model_margin, then lod_vs_second when the
// pair leads or single_lod_vs_second when the single leads. Empty when every
// margin is clear; a NaN margin counts as clear.
[[nodiscard]] auto find_ambiguity(
    const SampleScores& scores,
    const CallOptions& options) -> std::optional<Ambiguity>;

// ambiguous when find_ambiguity finds a margin, otherwise the leading model.
[[nodiscard]] auto assign_status(
    const SampleScores& scores,
    const CallOptions& options) -> Status;

// pass when the best pair's hom_mismatch_rate <= max_hom_mismatch, fail
// above it, na when the pair model does not lead or the rate is NaN.
[[nodiscard]] auto assign_fit(
    const SampleScores& scores,
    const CallOptions& options) -> Fit;

// A sample's expected parents as named in the expected-pairs file.
struct ExpectedParents
{
    std::string a;
    std::string b;
    std::optional<ParentPair> pair;  // unset when a name is not in the panel
};

// Reads `sample<TAB>mother<TAB>father` lines (mother == father for selfing)
// and resolves the names against `parent_names`. A first line starting with
// `sample` is a header; blank lines and `#` comments are skipped. Throws on a
// line without three fields or a repeated sample.
[[nodiscard]] auto read_expected_pairs(
    const std::string& path,
    std::span<const std::string> parent_names)
    -> std::unordered_map<std::string, ExpectedParents>;

// Writes one table per status, `<prefix>.<status>.tsv` (one row per sample,
// with the columns that status needs), and `<prefix>.topk.tsv` (the top
// pairs of every sample) as samples are added, and `<prefix>.summary.txt` on
// `finish()`. Rows are flushed as they are written. Not thread-safe.
class ReportWriter
{
   public:
    auto add(
        const SampleCounts& counts,
        const SampleScores& scores,
        const ExpectedParents* expected) -> void;

    // Counts a sample left out of the tables, reported in the summary.
    auto add_dropped() -> void { ++n_dropped_; }

    // Writes the dropped count, the status counts, the fit counts of TWO_PARENT samples, and the
    // expected-status counts when expected pairs were given. Throws when a file cannot be written.
    auto finish() -> void;

   private:
    friend auto open_report(
        const std::string& prefix,
        std::span<const std::string> parent_names,
        bool with_expected,
        const CallOptions& options) -> ReportWriter;

    ReportWriter(
        std::string prefix,
        std::span<const std::string> parent_names,
        bool with_expected,
        const CallOptions& options,
        std::array<std::ofstream, all_statuses.size()> tables,
        std::ofstream topk);

    auto write_expected(
        std::ostream& out,
        const SampleScores& scores,
        const ExpectedParents* expected) -> std::string_view;

    std::string prefix_;
    std::span<const std::string> names_;
    bool with_expected_;
    CallOptions options_;
    std::array<std::ofstream, all_statuses.size()> tables_;  // by Status
    std::ofstream topk_;
    size_t n_samples_ = 0;  // scored
    size_t n_dropped_ = 0;
    // (status, expected_status) -> samples
    std::map<std::pair<std::string_view, std::string_view>, size_t> tally_;
    std::map<std::string_view, size_t> fits_;  // of TWO_PARENT samples
};

// Opens the per-sample outputs and writes their headers; `with_expected`
// adds the expected-pair columns. `parent_names` must outlive the writer.
// Throws when a file cannot be created.
[[nodiscard]] auto open_report(
    const std::string& prefix,
    std::span<const std::string> parent_names,
    bool with_expected,
    const CallOptions& options) -> ReportWriter;
}  // namespace vcfbox::parentage
