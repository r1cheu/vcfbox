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
constexpr uint16_t kSkipMask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY
                               | BAM_FQCFAIL | BAM_FDUP;

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

void count_site(
    htsFile* file,
    hts_idx_t* idx,
    int tid,
    const Site& site,
    const PileupOptions& opts,
    bam1_t* rec,
    AlleleCount& out)
{
    hts_itr_ptr it(sam_itr_queryi(idx, tid, site.pos, site.pos + 1));
    if (!it)
    {
        return;
    }
    const char ref_upper = static_cast<char>(
        std::toupper(static_cast<unsigned char>(site.ref)));
    const char alt_upper = static_cast<char>(
        std::toupper(static_cast<unsigned char>(site.alt)));
    while (sam_itr_next(file, it.get(), rec) >= 0)
    {
        if ((rec->core.flag & kSkipMask) != 0)
        {
            continue;
        }
        if (rec->core.qual < opts.min_mapq)
        {
            continue;
        }
        const int32_t qpos = ref_to_query_pos(rec, site.pos);
        if (qpos < 0)
        {
            continue;
        }
        const uint8_t baseq = bam_get_qual(rec)[qpos];
        if (baseq < opts.min_baseq)
        {
            continue;
        }
        const uint8_t enc = bam_seqi(bam_get_seq(rec), qpos);
        const char base = seq_nt16_str[enc];
        if (base == ref_upper)
        {
            ++out.n_ref;
        }
        else if (base == alt_upper)
        {
            ++out.n_alt;
        }
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

    std::string last_chrom;
    int last_tid = -1;
    for (size_t i = 0; i < sites.size(); ++i)
    {
        const Site& site = sites[i];
        if (site.chrom != last_chrom)
        {
            last_tid = sam_hdr_name2tid(hdr.get(), site.chrom.c_str());
            last_chrom = site.chrom;
        }
        if (last_tid < 0)
        {
            continue;
        }
        count_site(
            file.get(), idx.get(), last_tid, site, opts, rec.get(), counts[i]);
    }
    return counts;
}
}  // namespace vcfbox
