#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace vcfbox
{
struct Site
{
    std::string chrom;
    int64_t pos;
    char ref;
    char alt;
};

std::vector<Site> load_sites(const std::string& path);
}  // namespace vcfbox
