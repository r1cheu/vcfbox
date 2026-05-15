#include "parentage/pileup.h"

#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "hts/hts_raii.h"

namespace vcfbox
{
namespace
{
constexpr uint16_t kSkipMask
    = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP;

int32_t ref_to_query_pos(const bam1_t* rec, int64_t target)
{
    const uint32_t* cigar = bam_get_cigar(rec);
    const uint32_t n_cigar = rec->core.n_cigar;
    int64_t ref_pos = rec->core.pos;
    int32_t q_pos = 0;
    for (uint32_t i = 0; i < n_cigar; ++i)
    {
        const auto op = bam_cigar_op(cigar[i]);
        const auto len = bam_cigar_oplen(cigar[i]);
        const int type = bam_cigar_type(op);
        const bool consume_query = (type & 1) != 0;
        const bool consume_ref = (type & 2) != 0;
        if (consume_ref && consume_query)
        {
            const int64_t end = ref_pos + len;
            if (target >= ref_pos && target < end)
            {
                return q_pos + static_cast<int32_t>(target - ref_pos);
            }
        }
        else if (consume_ref)
        {
            const int64_t end = ref_pos + len;
            if (target >= ref_pos && target < end)
            {
                return -1;
            }
        }
        if (consume_ref)
        {
            ref_pos += static_cast<int64_t>(len);
        }
        if (consume_query)
        {
            q_pos += static_cast<int32_t>(len);
        }
        if (ref_pos > target)
        {
            break;
        }
    }
    return -1;
}

char to_upper(char c)
{
    return static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
}

void tally_read(
    const bam1_t* rec,
    std::span<const Site> sites,
    size_t& window_lo,
    size_t window_hi,
    const PileupOptions& opts,
    std::span<AlleleCount> counts)
{
    if ((rec->core.flag & kSkipMask) != 0 || rec->core.qual < opts.min_mapq)
    {
        return;
    }
    const int64_t read_start = rec->core.pos;
    const int64_t read_end = bam_endpos(rec);
    while (window_lo < window_hi && sites[window_lo].pos < read_start)
    {
        ++window_lo;
    }
    for (size_t s = window_lo; s < window_hi && sites[s].pos < read_end; ++s)
    {
        const int32_t qpos = ref_to_query_pos(rec, sites[s].pos);
        if (qpos < 0)
        {
            continue;
        }
        if (bam_get_qual(rec)[qpos] < opts.min_baseq)
        {
            continue;
        }
        const char base = seq_nt16_str[bam_seqi(bam_get_seq(rec), qpos)];
        if (base == to_upper(sites[s].ref))
        {
            ++counts[s].n_ref;
        }
        else if (base == to_upper(sites[s].alt))
        {
            ++counts[s].n_alt;
        }
    }
}

void scan_chrom(
    htsFile* file,
    hts_idx_t* idx,
    int tid,
    std::span<const Site> sites,
    size_t chrom_lo,
    size_t chrom_hi,
    const PileupOptions& opts,
    bam1_t* rec,
    std::span<AlleleCount> counts)
{
    const int64_t lo = sites[chrom_lo].pos;
    const int64_t hi = sites[chrom_hi - 1].pos + 1;
    hts_itr_ptr it(sam_itr_queryi(idx, tid, lo, hi));
    if (!it)
    {
        return;
    }
    size_t window_lo = chrom_lo;
    while (sam_itr_next(file, it.get(), rec) >= 0)
    {
        tally_read(rec, sites, window_lo, chrom_hi, opts, counts);
    }
}
}  // namespace

std::vector<AlleleCount> count_alleles(
    const std::string& bam_path,
    std::span<const Site> sites,
    const PileupOptions& opts)
{
    hts_file_ptr file(sam_open(bam_path.c_str(), "r"));
    if (!file)
    {
        throw std::runtime_error("Cannot open BAM: " + bam_path);
    }
    sam_hdr_ptr hdr(sam_hdr_read(file.get()));
    if (!hdr)
    {
        throw std::runtime_error("Cannot read BAM header: " + bam_path);
    }
    hts_idx_ptr idx(sam_index_load(file.get(), bam_path.c_str()));
    if (!idx)
    {
        throw std::runtime_error(
            "Cannot load BAM index for: " + bam_path
            + " (run `samtools index`)");
    }

    bam1_ptr rec(bam_init1());
    std::vector<AlleleCount> counts(sites.size(), AlleleCount{0, 0});

    size_t i = 0;
    while (i < sites.size())
    {
        size_t j = i + 1;
        while (j < sites.size() && sites[j].chrom == sites[i].chrom)
        {
            ++j;
        }
        const int tid = sam_hdr_name2tid(hdr.get(), sites[i].chrom.c_str());
        if (tid >= 0)
        {
            scan_chrom(
                file.get(),
                idx.get(),
                tid,
                sites,
                i,
                j,
                opts,
                rec.get(),
                counts);
        }
        i = j;
    }
    return counts;
}
}  // namespace vcfbox
