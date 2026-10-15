#include "parentage/parentage.h"

#include <filesystem>
#include <format>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "common/lines.h"
#include "common/strings.h"
#include "parentage/counts.h"
#include "parentage/counts_io.h"
#include "parentage/panel.h"
#include "parentage/pileup.h"
#include "parentage/score.h"

namespace fs = std::filesystem;

namespace vcfbox::parentage
{
namespace
{
auto print_panel_summary(const Panel& panel) -> void
{
    std::cerr << std::format(
        "Panel: {} parents, {} sites on {} contigs{}, {} heterozygous calls\n",
        panel.n_parents(),
        panel.n_sites(),
        panel.sites().contigs().size(),
        panel.n_excluded_sites() > 0
            ? std::format(" ({} excluded)", panel.n_excluded_sites())
            : "",
        panel.n_het_calls());

    for (const auto& dup : panel.near_duplicates())
    {
        std::cerr << std::format(
            "Warning: parents {} and {} differ at only {} of {} sites; "
            "samples derived from either will be AMBIGUOUS\n",
            panel.names()[dup.parent_a],
            panel.names()[dup.parent_b],
            dup.n_discordant,
            panel.n_sites());
    }
}

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

// Temporary stderr output; the TSV report replaces it.
auto print_scores(
    const Panel& panel,
    const std::string& sample,
    const SampleScores& scores) -> void
{
    const std::span<const std::string> names = panel.names();
    const PairReport& best = scores.top.front();
    const SingleStats& single = scores.singles.front();
    std::cerr << std::format(
        "Sample {}: {} covered sites; best {} x {} (lod_vs_null {:.1f}, "
        "lod_vs_second {:.1f}, hom_mismatch {:.4f}, het_alt {:.3f}); best "
        "single {} (lod_vs_null {:.1f}, allele_absent {:.4f})\n",
        sample,
        scores.n_covered,
        names[best.stats.a],
        names[best.stats.b],
        best.lod_vs_null,
        scores.lod_vs_second,
        best.hom_mismatch_rate,
        best.het_alt_fraction,
        names[single.parent],
        single.lod_vs_null,
        single.allele_absent_rate);
    std::cerr << "  rank\ta\tb\tscore\thom_conc\thom_disc\thet_0\thet_1\t"
                 "lod_vs_null\thom_mismatch\thet_alt\n";
    for (size_t rank = 0; rank < scores.top.size(); ++rank)
    {
        const PairReport& report = scores.top[rank];
        const PairStats& st = report.stats;
        std::cerr << std::format(
            "  {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.4f}\t{:.3f}\n",
            rank + 1,
            names[st.a],
            names[st.b],
            st.score,
            st.hom_conc,
            st.hom_disc,
            st.het_0,
            st.het_1,
            report.lod_vs_null,
            report.hom_mismatch_rate,
            report.het_alt_fraction);
    }
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
    const PileupOptions& pileup_options) -> SampleCounts
{
    if (options.counts_dir.empty())
    {
        SampleCounts counts
            = pileup_bam(bam.path, bam.sample, sites, pileup_options);
        print_sample_summary(counts, "pileup");
        return counts;
    }
    fs::create_directories(options.counts_dir);
    const fs::path file = fs::path(options.counts_dir)
                          / std::format("{}.counts.gz", bam.sample);
    if (fs::exists(file))
    {
        SampleCounts counts = read_counts(file.string(), sites.size());
        print_sample_summary(counts, "counts file");
        return counts;
    }
    SampleCounts counts
        = pileup_bam(bam.path, bam.sample, sites, pileup_options);
    write_counts(counts, file.string());
    print_sample_summary(counts, "pileup");
    return counts;
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
    if (options.top == 0)
    {
        throw std::runtime_error("--top must be at least 1");
    }
    const std::vector<BamEntry> bams = read_bam_list(options.bams_path);
    PileupOptions pileup_options;
    pileup_options.min_mapq = options.min_mapq;
    pileup_options.min_baseq = options.min_baseq;
    pileup_options.dedup_overlaps = !options.keep_overlaps;
    ScoreOptions score_options;
    score_options.error_rate = options.error_rate;
    score_options.top = options.top;
    score_options.candidates = options.candidates;
    for (const auto& bam : bams)
    {
        const SampleCounts counts
            = obtain_counts(bam, panel.sites(), options, pileup_options);
        print_scores(
            panel, bam.sample, score_sample(panel, counts, score_options));
    }
}
}  // namespace vcfbox::parentage
