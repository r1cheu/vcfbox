#include "common/strings.h"

namespace vcfbox
{
std::string join(std::span<const std::string> parts, std::string_view sep)
{
    std::string out;
    for (size_t i = 0; i < parts.size(); ++i)
    {
        if (i != 0)
        {
            out.append(sep);
        }
        out.append(parts[i]);
    }
    return out;
}
}  // namespace vcfbox
