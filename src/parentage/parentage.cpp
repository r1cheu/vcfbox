#include "parentage/parentage.h"

#include <atomic>
#include <exception>
#include <filesystem>
#include <format>
#include <iostream>
#include <mutex>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/lines.h"
#include "common/strings.h"
#include "parentage/counts.h"
#include "parentage/counts_io.h"
#include "parentage/panel.h"
#include "parentage/pileup.h"
#include "parentage/report.h"
#include "parentage/score.h"

namespace fs = std::filesystem;

namespace vcfbox::parentage
{
namespace
{
auto print_panel_summary(const Panel& panel) -> void
{
    std::cerr << std::format(
        "Panel: {} parents, {} sites on {} contigs{}, {} heterozygous and {} "
        "missing calls\n",
        panel.n_parents(),
        panel.n_sites(),
        panel.sites().contigs().size(),
        panel.n_excluded_sites() > 0
            ? std::format(" ({} excluded)", panel.n_excluded_sites())
            : "",
        panel.n_het_calls(),
        panel.n_missing_calls());

    for (const auto& dup : panel.near_duplicates())
    {
        std::cerr << std::format(
            "Warning: parents {} and {} differ at only {} of {} sites where "
            "both are homozygous; "
            "samples derived from either will be AMBIGUOUS\n",
            panel.names()[dup.parent_a],
            panel.names()[dup.parent_b],
            dup.n_discordant,
            dup.n_compared);
    }
}

// Counts of one sample and where they came from.
struct ObtainedCounts
{
    SampleCounts counts;
    std::string origin;  // "pileup" or "counts file"
};

auto print_sample_summary(const SampleCounts& counts, const std::string& origin)
    -> void
{
    const PileupStats& stats = counts.stats();
    std::cerr << std::format(
        "Sample {} ({}): {} of {} reads used ({} flag, {} improper pair, {} "
        "MAPQ, {} off-panel skipped); {} of {} sites covered, mean depth "
        "{:.3f}; {} bases below baseQ\n",
        counts.sample(),
        origin,
        stats.reads_used,
        stats.reads_total,
        stats.reads_skipped_flag,
        stats.reads_skipped_improper,
        stats.reads_skipped_mapq,
        stats.reads_off_panel,
        counts.sites_covered(),
        counts.n_sites(),
        static_cast<double>(counts.bases_used())
            / static_cast<double>(counts.n_sites()),
        stats.bases_skipped_baseq);
}

auto sample_name_from_path(const std::string& path) -> std::string
{
    const fs::path p(path);
    return p.extension() == ".bam" ? p.stem().string() : p.filename().string();
}

// Counts from `counts_dir` when a file for the sample exists there, else a
// fresh pileup, written to `counts_dir` when one is given.
auto obtain_counts(
    const BamEntry& bam,
    const SiteTable& sites,
    const RunOptions& options,
    const PileupOptions& pileup_options) -> ObtainedCounts
{
    if (options.counts_dir.empty())
    {
        return {
            pileup_bam(bam.path, bam.sample, sites, pileup_options), "pileup"};
    }
    const fs::path file = fs::path(options.counts_dir)
                          / std::format("{}.counts.gz", bam.sample);
    if (fs::exists(file))
    {
        return {read_counts(file.string(), sites.size()), "counts file"};
    }
    SampleCounts counts
        = pileup_bam(bam.path, bam.sample, sites, pileup_options);
    write_counts(counts, file.string());
    return {std::move(counts), "pileup"};
}

// Processes the samples on `options.threads` workers, each taking the next
// unprocessed sample, and adds them to `report` in completion order. The
// first exception thrown by a worker stops the others after their current
// sample and is rethrown here.
auto process_samples(
    const Panel& panel,
    std::span<const BamEntry> bams,
    const std::unordered_map<std::string, ExpectedParents>& expected,
    const RunOptions& options,
    const PileupOptions& pileup_options,
    const ScoreOptions& score_options,
    ReportWriter& report) -> void
{
    std::atomic<size_t> next = 0;
    std::atomic<bool> failed = false;
    std::mutex output_mutex;
    std::vector<std::exception_ptr> errors(options.threads);
    const auto worker = [&](size_t thread)
    {
        try
        {
            for (size_t i = next++; i < bams.size() && !failed; i = next++)
            {
                const ObtainedCounts obtained = obtain_counts(
                    bams[i], panel.sites(), options, pileup_options);
                const auto found = expected.find(bams[i].sample);
                const ExpectedParents* parents
                    = found == expected.end() ? nullptr : &found->second;
                const SampleScores scores = score_sample(
                    panel,
                    obtained.counts,
                    score_options,
                    parents == nullptr ? std::nullopt : parents->pair);
                const std::lock_guard lock(output_mutex);
                print_sample_summary(obtained.counts, obtained.origin);
                report.add(obtained.counts, scores, parents);
            }
        }
        catch (...)
        {
            errors[thread] = std::current_exception();
            failed = true;
        }
    };
    {
        std::vector<std::jthread> workers;
        workers.reserve(options.threads);
        for (size_t t = 0; t < options.threads; ++t)
        {
            workers.emplace_back(worker, t);
        }
    }
    for (const std::exception_ptr& error : errors)
    {
        if (error)
        {
            std::rethrow_exception(error);
        }
    }
}
}  // namespace

auto read_bam_list(const std::string& path) -> std::vector<BamEntry>
{
    std::vector<BamEntry> entries;
    std::unordered_set<std::string> seen;
    size_t line_no = 0;
    for (const auto& line : read_lines(path))
    {
        ++line_no;
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        const auto fields = split(line, '\t');
        if (fields.size() > 2)
        {
            throw std::runtime_error(
                std::format(
                    "{}:{}: expected `sample<TAB>path` or `path`, found {} "
                    "fields",
                    path,
                    line_no,
                    fields.size()));
        }
        BamEntry entry;
        entry.path = std::string(fields.back());
        entry.sample = fields.size() == 2 ? std::string(fields[0])
                                          : sample_name_from_path(entry.path);
        if (!seen.insert(entry.sample).second)
        {
            throw std::runtime_error(
                std::format(
                    "{}:{}: sample {} is listed twice",
                    path,
                    line_no,
                    entry.sample));
        }
        entries.push_back(std::move(entry));
    }
    return entries;
}

auto run(const RunOptions& options) -> void
{
    PanelOptions panel_options;
    panel_options.exclude_sites_path = options.exclude_sites_path;
    panel_options.duplicate_max_discordance = options.duplicate_max_discordance;

    const Panel panel = load_panel(options.vcf_path, panel_options);
    print_panel_summary(panel);

    if (options.bams_path.empty())
    {
        return;
    }
    if (options.top < 2)
    {
        throw std::runtime_error("--top must be at least 2");
    }
    if (options.threads == 0)
    {
        throw std::runtime_error("--threads must be at least 1");
    }
    if (options.output_prefix.empty())
    {
        throw std::runtime_error("--output is required with --bams");
    }
    const std::vector<BamEntry> bams = read_bam_list(options.bams_path);
    const std::unordered_map<std::string, ExpectedParents> expected
        = options.expected_path.empty()
              ? std::unordered_map<std::string, ExpectedParents>{}
              : read_expected_pairs(options.expected_path, panel.names());
    if (!options.counts_dir.empty())
    {
        fs::create_directories(options.counts_dir);
    }
    PileupOptions pileup_options;
    pileup_options.min_mapq = options.min_mapq;
    pileup_options.min_baseq = options.min_baseq;
    pileup_options.dedup_overlaps = !options.keep_overlaps;
    pileup_options.progress = options.threads == 1;
    ScoreOptions score_options;
    score_options.error_rate = options.error_rate;
    score_options.top = options.top;
    score_options.candidates = options.candidates;
    ReportWriter report = open_report(
        options.output_prefix,
        panel.names(),
        !options.expected_path.empty(),
        options.call);
    process_samples(
        panel,
        bams,
        expected,
        options,
        pileup_options,
        score_options,
        report);
    report.finish();
}
}  // namespace vcfbox::parentage
