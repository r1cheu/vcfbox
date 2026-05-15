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

std::vector<std::string_view> split(std::string_view s, char delim)
{
    std::vector<std::string_view> out;
    while (true)
    {
        const auto pos = s.find(delim);
        if (pos == std::string_view::npos)
        {
            out.push_back(s);
            return out;
        }
        out.push_back(s.substr(0, pos));
        s.remove_prefix(pos + 1);
    }
}
}  // namespace vcfbox
