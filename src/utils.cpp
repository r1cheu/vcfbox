#include "utils.h"
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

}  // namespace detail
