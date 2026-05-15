#pragma once

#include <span>
#include <string>
#include <string_view>

namespace vcfbox
{
std::string join(std::span<const std::string> parts, std::string_view sep);
}  // namespace vcfbox
