#pragma once

#include <string>

#include "parentage/counts.h"
#include "parentage/sites.h"

namespace vcfbox::parentage
{
struct PileupOptions
{
    int min_mapq = 20;
    int min_baseq = 20;
    // Count a base once when both mates of a pair cover the same site.
    bool dedup_overlaps = true;
};

// Scans the whole BAM sequentially (no index) and counts, at every panel
// site, the reads carrying allele0, allele1 or a third base. Throws when the
// file cannot be read, is not coordinate sorted, or its contigs do not match
// the panel (`match_contigs`).
[[nodiscard]] auto pileup_bam(
    const std::string& bam_path,
    const std::string& sample,
    const SiteTable& sites,
    const PileupOptions& options) -> SampleCounts;
}  // namespace vcfbox::parentage
