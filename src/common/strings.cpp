#include "common/strings.h"

#include <array>
#include <charconv>
#include <stdexcept>

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

std::string format_number(double value, int precision)
{
    std::array<char, 32> buf{};
    const auto [ptr, ec] = std::to_chars(
        buf.data(),
        buf.data() + buf.size(),
        value,
        std::chars_format::general,
        precision);
    if (ec != std::errc{})
    {
        throw std::runtime_error("Failed to format number");
    }
    return std::string(buf.data(), ptr);
}
}  // namespace vcfbox
