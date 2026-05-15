#include "common/progress.h"

namespace bk = barkeep;

namespace vcfbox
{
std::shared_ptr<barkeep::CompositeDisplay> create_counter(
    const std::string& message,
    size_t& progress_counters)
{
    auto anim = bk::Animation(
        {.style = bk::Strings{"|", "/", "-", "\\"},
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
}  // namespace vcfbox
