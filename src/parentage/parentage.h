#pragma once

#include <string>

namespace vcfbox::parentage
{
struct RunOptions
{
    std::string vcf_path;
    std::string exclude_sites_path;
    double duplicate_max_discordance = 0.001;
};

// Loads the panel and prints a summary. Scoring stages are added on top of
// this entry point as they are implemented.
auto run(const RunOptions& options) -> void;
}  // namespace vcfbox::parentage
