#include <string>
#include "barkeep.h"

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

namespace detail
{
namespace bk = barkeep;
auto create_progress(size_t total, size_t& progress_counters)
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
