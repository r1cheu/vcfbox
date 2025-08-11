#include "utils.h"

#include <set>
#include <string>

#include "barkeep.h"
#include "vcf_raii.h"

extern "C"
{
#include <htslib/vcf.h>
}
namespace bk = barkeep;
namespace vcfbox
{
std::string parse_mode(std::string_view file_path)
{
    auto ext = file_path.substr(file_path.find_last_of('.') + 1);
    if (ext == "bam")
    {
        return "wb";
    }
    if (ext == "sam")
    {
        return "w";
    }
    if (ext == "cram")
    {
        return "wc";
    }
    if (ext == "bcf")
    {
        return "wb";
    }
    if (ext == "vcf")
    {
        return "w";
    }
    if (ext == "gz")
    {
        auto base = file_path.substr(0, file_path.find_last_of('.'));
        auto base_ext = base.substr(base.find_last_of('.') + 1);
        if (base_ext == "vcf")
        {
            return "wz";
        }
        if (base_ext == "bcf")
        {
            return "wb";
        }
    }
    return "w";
}

size_t count_records(std::string_view vcf_path)
{
    size_t rec_count = 0;
    auto counter = bk::Counter(
        &rec_count,
        {
            .message = "Counting SNPs",
            .speed = 1.,
            .speed_unit = "snp/s",
        });

    HtsFile vcf_file(bcf_open(vcf_path.data(), "r"));
    if (!vcf_file)
    {
        throw std::runtime_error(
            "Could not open VCF file: " + std::string(vcf_path));
    }
    BcfHdr header(bcf_hdr_read(vcf_file.get()));
    if (!header)
    {
        throw std::runtime_error(
            "Could not read VCF header from: " + std::string(vcf_path));
    }
    BcfRec rec(bcf_init());
    if (!rec)
    {
        throw std::runtime_error("Failed to initialize VCF record.");
    }

    while (bcf_read(vcf_file.get(), header.get(), rec.get()) == 0)
    {
        rec_count++;
    }
    return rec_count;
}

}  // namespace vcfbox

namespace detail
{
namespace bk = barkeep;
std::shared_ptr<barkeep::CompositeDisplay> create_progress(
    size_t total,
    size_t& progress_counters)
{
    auto anim = bk::Animation(
        {.style = bk::Strings{"⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"},
         .interval = 0.08,
         .show = false});

    bk::BarParts custom_bar_style;
    custom_bar_style.left = "[";
    custom_bar_style.right = "]";
    custom_bar_style.fill = {"\033[1;33m━\033[0m"};
    custom_bar_style.empty = {"─"};

    auto pbar = bk::ProgressBar(
        &progress_counters,
        {.total = total,
         .format = "Adding SNPs {bar} {value}/{total} ({speed:.1f}/s)",
         .speed = 0.1,
         .style = custom_bar_style,
         .show = false});

    return bk::Composite({anim, pbar}, " ");
}

std::shared_ptr<barkeep::CompositeDisplay> create_counter(
    const std::string& message,
    size_t& progress_counters)
{
    auto anim = bk::Animation(
        {.style = bk::Strings{"⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"},
         .interval = 0.08,
         .show = false});

    auto pbar = bk::Counter(
        &progress_counters,
        {
            .message = message,
            .speed = 1.,
            .speed_unit = "snp/s",
            .show = false,
        });

    return bk::Composite({anim, pbar}, " ");
}
void check_sample_consistence(
    std::string_view vcf_path,
    const std::vector<SamplePair>& sample_pairs)
{
    HtsFile vcf_file(bcf_open(vcf_path.data(), "r"));
    if (!vcf_file)
    {
        throw std::runtime_error(
            "Could not open VCF file: " + std::string(vcf_path));
    }

    BcfHdr header(bcf_hdr_read(vcf_file.get()));

    if (!header)
    {
        throw std::runtime_error(
            "Could not read VCF header from: " + std::string(vcf_path));
    }

    char** sample_names = header->samples;
    size_t num_samples = bcf_hdr_nsamples(header);

    std::set<std::string> all_samples;

    for (size_t i = 0; i < num_samples; ++i)
    {
        all_samples.insert(sample_names[i]);
    }

    std::vector<std::string> missing_samples;

    for (const auto& pair : sample_pairs)
    {
        if (!all_samples.contains(pair.first))
        {
            missing_samples.push_back(pair.first);
        }
        if (!all_samples.contains(pair.second))
        {
            missing_samples.push_back(pair.second);
        }
    }

    if (!missing_samples.empty())
    {
        std::string missing_str;
        for (const auto& s : missing_samples)
        {
            missing_str += s + ", ";
        }
        throw std::runtime_error(
            "Samples not found in VCF: " + missing_str
            + "make sure sample list is correct and matches VCF file.");
    }
}

bcf_hdr_t* init_bcf_head(
    bcf_hdr_t* header,
    const std::vector<SamplePair>& sample_pairs,
    bool keep_old_samples)
{
    bcf_hdr_t* output_header = bcf_hdr_init("w");

    for (int i = 0; i < header->nhrec; ++i)
    {
        if (header->hrec[i]->type == BCF_HL_CTG)
        {
            bcf_hdr_add_hrec(output_header, bcf_hrec_dup(header->hrec[i]));
        }
    }

    bcf_hrec_t* gt_hrec
        = bcf_hdr_get_hrec(header, BCF_HL_FMT, "ID", "GT", nullptr);
    if (gt_hrec != nullptr)
    {
        bcf_hdr_add_hrec(output_header, bcf_hrec_dup(gt_hrec));
    }
    else
    {
        bcf_hdr_append(
            output_header,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    }

    if (keep_old_samples)
    {
        for (int i = 0; i < bcf_hdr_nsamples(header); ++i)
        {
            bcf_hdr_add_sample(output_header, header->samples[i]);
        }
    }

    for (const auto& pair : sample_pairs)
    {
        bcf_hdr_add_sample(
            output_header, (pair.first + "_" + pair.second).c_str());
    }

    bcf_hdr_add_sample(output_header, nullptr);  // 更新样本列表

    return output_header;
}

void copy_rec_info(
    bcf_hdr_t* header,
    bcf_hdr_t* output_header,
    bcf1_t* in_rec,
    bcf1_t* out_rec)
{
    bcf_clear(out_rec);
    out_rec->rid
        = bcf_hdr_name2id(output_header, bcf_hdr_id2name(header, in_rec->rid));
    out_rec->pos = in_rec->pos;
    bcf_update_id(output_header, out_rec, in_rec->d.id);
    bcf_update_alleles(
        output_header,
        out_rec,
        const_cast<const char**>(in_rec->d.allele),
        in_rec->n_allele);

    out_rec->qual = in_rec->qual;
}

std::vector<int32_t> concat_gt(
    const std::vector<std::pair<std::string, std::string>>& sample_pairs,
    const std::unordered_map<std::string, int>& sample_to_idx,
    const int32_t* gt_arr,
    bool keep_old_samples,
    int n_gt)
{
    std::vector<int32_t> out_gts;
    if (keep_old_samples)
    {
        out_gts.reserve(n_gt + (sample_pairs.size() * 2));
        for (int i = 0; i < n_gt / 2; ++i)
        {
            int32_t gt0 = gt_arr[i * 2];
            int32_t gt1 = gt_arr[i * 2 + 1];
            if (bcf_gt_is_missing(gt0)
                || bcf_gt_allele(gt0) != bcf_gt_allele(gt1))
            {
                out_gts.push_back(bcf_gt_missing);
                out_gts.push_back(bcf_gt_missing);
            }
            else
            {
                out_gts.push_back(gt0);
                out_gts.push_back(gt1);
            }
        }
    }
    else
    {
        out_gts.reserve(sample_pairs.size() * 2);
    }
    for (const auto& pair : sample_pairs)
    {
        int f_idx = sample_to_idx.at(pair.first);
        int m_idx = sample_to_idx.at(pair.second);

        int32_t f_gt0 = gt_arr[f_idx * 2];
        int32_t f_gt1 = gt_arr[f_idx * 2 + 1];
        int32_t m_gt0 = gt_arr[m_idx * 2];
        int32_t m_gt1 = gt_arr[m_idx * 2 + 1];

        if (bcf_gt_is_missing(f_gt0) || bcf_gt_is_missing(f_gt1)
            || bcf_gt_is_missing(m_gt0) || bcf_gt_is_missing(m_gt1)
            || bcf_gt_allele(f_gt0) != bcf_gt_allele(f_gt1)
            || bcf_gt_allele(m_gt0) != bcf_gt_allele(m_gt1))
        {
            out_gts.push_back(bcf_gt_missing);
            out_gts.push_back(bcf_gt_missing);
        }
        else
        {
            out_gts.push_back(bcf_gt_unphased(bcf_gt_allele(f_gt0)));
            out_gts.push_back(bcf_gt_unphased(bcf_gt_allele(m_gt0)));
        }
    }
    return out_gts;
}

}  // namespace detail
