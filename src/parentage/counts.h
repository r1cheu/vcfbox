#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace vcfbox::parentage
{
struct PileupStats
{
    uint64_t reads_total = 0;         // records in the BAM
    uint64_t reads_used = 0;          // passed every read-level filter
    uint64_t reads_skipped_flag = 0;  // unmapped/secondary/supplementary/...
    uint64_t reads_skipped_improper = 0;  // paired but not a proper pair
    uint64_t reads_skipped_mapq = 0;      // below --min-mapq or 255
    uint64_t reads_off_panel = 0;         // contig the panel has no sites on
    uint64_t bases_skipped_baseq = 0;     // panel-site bases below --min-baseq
                                          // (includes mate-overlap duplicates,
                                          // whose quality is set to 0)
};

inline auto saturating_increment(uint8_t& n) -> void
{
    if (n < UINT8_MAX)
    {
        ++n;
    }
}

// Per-site read counts of one sample, indexed like the SiteTable. Counts
// saturate at 255.
class SampleCounts
{
   public:
    SampleCounts(std::string sample, size_t n_sites);

    // Adopts arrays read back from a counts file. Throws when their sizes
    // differ.
    SampleCounts(
        std::string sample,
        std::vector<uint8_t> n0,
        std::vector<uint8_t> n1,
        std::vector<uint8_t> n_other,
        PileupStats stats);

    auto sample() const -> const std::string& { return sample_; }
    auto n_sites() const -> size_t { return n0_.size(); }
    auto n0() const -> std::span<const uint8_t> { return n0_; }
    auto n1() const -> std::span<const uint8_t> { return n1_; }
    auto n_other() const -> std::span<const uint8_t> { return n_other_; }
    auto stats() const -> const PileupStats& { return stats_; }
    auto stats() -> PileupStats& { return stats_; }

    auto add0(size_t s) -> void { saturating_increment(n0_[s]); }
    auto add1(size_t s) -> void { saturating_increment(n1_[s]); }
    auto add_other(size_t s) -> void { saturating_increment(n_other_[s]); }

    // Sites with at least one allele0 or allele1 read.
    auto sites_covered() const -> size_t;
    // Sum of allele0 and allele1 reads over all sites.
    auto bases_used() const -> uint64_t;

   private:
    std::string sample_;
    std::vector<uint8_t> n0_;
    std::vector<uint8_t> n1_;
    std::vector<uint8_t> n_other_;
    PileupStats stats_;
};
}  // namespace vcfbox::parentage
