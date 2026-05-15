#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "parentage/sites.h"

namespace vcfbox
{
struct AlleleCount
{
    uint32_t n_ref;
    uint32_t n_alt;
};

struct PileupOptions
{
    int min_mapq = 20;
    int min_baseq = 13;
};

std::vector<AlleleCount> count_alleles(
    const std::string& bam_path,
    std::span<const Site> sites,
    const PileupOptions& opts);
}  // namespace vcfbox
