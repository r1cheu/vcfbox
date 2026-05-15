#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include <barkeep.h>

namespace vcfbox
{
std::shared_ptr<barkeep::CompositeDisplay> create_counter(
    const std::string& message,
    size_t& progress_counters,
    const std::string& speed_unit = "snp/s");
}  // namespace vcfbox
