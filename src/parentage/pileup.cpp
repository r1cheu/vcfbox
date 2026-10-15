#include "parentage/pileup.h"

#include <cstdint>
#include <exception>
#include <format>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "common/progress.h"
#include "hts/hts_raii.h"

extern "C"
{
#include <htslib/sam.h>
}

namespace vcfbox::parentage
{
namespace
{
constexpr uint16_t skip_flags
    = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP;
constexpr int mapq_unavailable = 255;
constexpr int max_depth = 8000;  // samtools mpileup default (-d)

struct OpenedBam
{
    hts_file_ptr file;
    sam_hdr_ptr hdr;
};

auto open_bam(const std::string& bam_path) -> OpenedBam
{
    hts_file_ptr file(sam_open(bam_path.c_str(), "r"));
    if (!file)
    {
        throw std::runtime_error(
            std::format("Failed to open BAM: {}", bam_path));
    }
    sam_hdr_ptr hdr(sam_hdr_read(file.get()));
    if (!hdr)
    {
        throw std::runtime_error(
            std::format("Failed to read BAM header: {}", bam_path));
    }
    return {std::move(file), std::move(hdr)};
}

// Panel contig index for every BAM tid, or -1 when the panel has no sites on
// that contig.
auto panel_contig_of_tid(
    const SiteTable& sites,
    const sam_hdr_t* hdr,
    const std::string& bam_path) -> std::vector<int>
{
    const int n_targets = sam_hdr_nref(hdr);
    std::vector<ContigInfo> bam_contigs;
    bam_contigs.reserve(static_cast<size_t>(n_targets));
    for (int tid = 0; tid < n_targets; ++tid)
    {
        bam_contigs.push_back(
            {sam_hdr_tid2name(hdr, tid),
             static_cast<int64_t>(sam_hdr_tid2len(hdr, tid))});
    }
    const std::vector<int> tids = match_contigs(sites, bam_contigs, bam_path);

    std::vector<int> contig_of_tid(static_cast<size_t>(n_targets), -1);
    for (size_t contig = 0; contig < tids.size(); ++contig)
    {
        contig_of_tid[static_cast<size_t>(tids[contig])]
            = static_cast<int>(contig);
    }
    return contig_of_tid;
}

// State of the read callback that feeds htslib's pileup engine. Reads are
// filtered here, before the engine sees them, exactly as samtools mpileup
// does; errors are parked in `error` because the callback is called from C.
struct ReadSource
{
    htsFile* file = nullptr;
    sam_hdr_t* hdr = nullptr;
    const std::string* bam_path = nullptr;
    const PileupOptions* options = nullptr;
    std::span<const int> contig_of_tid;
    PileupStats* stats = nullptr;
    size_t* processed = nullptr;
    int prev_tid = -1;
    int64_t prev_pos = -1;
    std::exception_ptr error;
};

auto check_sorted(ReadSource& src, const bam1_t* rec) -> void
{
    const int tid = rec->core.tid;
    const int64_t pos = rec->core.pos;
    if (tid < src.prev_tid || (tid == src.prev_tid && pos < src.prev_pos))
    {
        throw std::runtime_error(
            std::format(
                "{}: BAM is not coordinate sorted ({}:{} follows {}:{})",
                *src.bam_path,
                sam_hdr_tid2name(src.hdr, tid),
                pos + 1,
                sam_hdr_tid2name(src.hdr, src.prev_tid),
                src.prev_pos + 1));
    }
    src.prev_tid = tid;
    src.prev_pos = pos;
}

// Whether a read takes part in the pileup; updates the skip counters.
auto passes_read_filters(ReadSource& src, const bam1_t* rec) -> bool
{
    const auto flag = rec->core.flag;
    PileupStats& stats = *src.stats;
    if ((flag & skip_flags) != 0)
    {
        ++stats.reads_skipped_flag;
        return false;
    }
    check_sorted(src, rec);
    if ((flag & BAM_FPAIRED) != 0 && (flag & BAM_FPROPER_PAIR) == 0)
    {
        ++stats.reads_skipped_improper;
        return false;
    }
    if (rec->core.qual < src.options->min_mapq
        || rec->core.qual == mapq_unavailable)
    {
        ++stats.reads_skipped_mapq;
        return false;
    }
    if (src.contig_of_tid[static_cast<size_t>(rec->core.tid)] < 0)
    {
        ++stats.reads_off_panel;
        return false;
    }
    ++stats.reads_used;
    return true;
}

// bam_plp_auto_f: next read that passes the filters; -1 at end of file, a
// smaller value on error (parked in `src.error`).
auto next_read(void* data, bam1_t* rec) -> int
{
    auto& src = *static_cast<ReadSource*>(data);
    try
    {
        while (true)
        {
            const int ret = sam_read1(src.file, src.hdr, rec);
            if (ret < -1)
            {
                throw std::runtime_error(
                    std::format(
                        "{}: error while reading record {}",
                        *src.bam_path,
                        *src.processed + 1));
            }
            if (ret < 0)
            {
                return ret;
            }
            ++*src.processed;
            ++src.stats->reads_total;
            if (passes_read_filters(src, rec))
            {
                return ret;
            }
        }
    }
    catch (...)
    {
        src.error = std::current_exception();
        return -2;
    }
}

auto count_column(
    std::span<const bam_pileup1_t> column,
    size_t site,
    const SiteTable& sites,
    const PileupOptions& options,
    SampleCounts& counts) -> void
{
    for (const bam_pileup1_t& p : column)
    {
        if (p.is_del != 0 || p.is_refskip != 0)
        {
            continue;
        }
        const auto qpos = static_cast<size_t>(p.qpos);
        if (bam_get_qual(p.b)[qpos] < options.min_baseq)
        {
            ++counts.stats().bases_skipped_baseq;
            continue;
        }
        const uint8_t code
            = encode_base(seq_nt16_str[bam_seqi(bam_get_seq(p.b), qpos)]);
        if (code == sites.allele0(site))
        {
            counts.add0(site);
        }
        else if (code == sites.allele1(site))
        {
            counts.add1(site);
        }
        else if (code != no_base)
        {
            counts.add_other(site);
        }
    }
}
}  // namespace

auto pileup_bam(
    const std::string& bam_path,
    const std::string& sample,
    const SiteTable& sites,
    const PileupOptions& options) -> SampleCounts
{
    auto [file, hdr] = open_bam(bam_path);
    const std::vector<int> contig_of_tid
        = panel_contig_of_tid(sites, hdr.get(), bam_path);

    SampleCounts counts(sample, sites.size());
    size_t processed = 0;
    auto counter
        = create_counter(std::format("Pileup {}", sample), processed, "read/s");
    if (options.progress)
    {
        counter->show();
    }

    ReadSource src;
    src.file = file.get();
    src.hdr = hdr.get();
    src.bam_path = &bam_path;
    src.options = &options;
    src.contig_of_tid = contig_of_tid;
    src.stats = &counts.stats();
    src.processed = &processed;
    void* data = &src;

    bam_mplp_ptr iter(bam_mplp_init(1, next_read, &data));
    if (options.dedup_overlaps)
    {
        bam_mplp_init_overlaps(iter.get());
    }
    bam_mplp_set_maxcnt(iter.get(), max_depth);

    const std::span<const int32_t> positions = sites.positions();
    int current_tid = -1;
    const ContigRange* contig = nullptr;
    size_t cursor = 0;
    size_t contig_end = 0;

    int tid = 0;
    hts_pos_t pos = 0;
    int n_plp = 0;
    const bam_pileup1_t* plp = nullptr;
    while (bam_mplp64_auto(iter.get(), &tid, &pos, &n_plp, &plp) > 0)
    {
        if (tid != current_tid)
        {
            current_tid = tid;
            const int c = contig_of_tid[static_cast<size_t>(tid)];
            contig = c < 0 ? nullptr : &sites.contigs()[static_cast<size_t>(c)];
            cursor = contig == nullptr ? 0 : contig->first_site;
            contig_end = contig == nullptr ? 0 : cursor + contig->n_sites;
        }
        while (cursor < contig_end && positions[cursor] < pos)
        {
            ++cursor;
        }
        if (cursor < contig_end && positions[cursor] == pos)
        {
            count_column(
                std::span(plp, static_cast<size_t>(n_plp)),
                cursor,
                sites,
                options,
                counts);
        }
    }
    if (src.error)
    {
        std::rethrow_exception(src.error);
    }
    if (options.progress)
    {
        counter->done();
    }
    return counts;
}
}  // namespace vcfbox::parentage
