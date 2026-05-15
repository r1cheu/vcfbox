#pragma once

#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace vcfbox
{
std::string join(std::span<const std::string> parts, std::string_view sep);

std::vector<std::string_view> split(std::string_view s, char delim);

std::string format_number(double value, int precision);
}  // namespace vcfbox
