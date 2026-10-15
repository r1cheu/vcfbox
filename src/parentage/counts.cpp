#include "parentage/counts.h"

#include <format>
#include <stdexcept>
#include <utility>

namespace vcfbox::parentage
{
SampleCounts::SampleCounts(std::string sample, size_t n_sites)
    : sample_(std::move(sample)),
      n0_(n_sites, 0),
      n1_(n_sites, 0),
      n_other_(n_sites, 0)
{
}

SampleCounts::SampleCounts(
    std::string sample,
    std::vector<uint8_t> n0,
    std::vector<uint8_t> n1,
    std::vector<uint8_t> n_other,
    PileupStats stats)
    : sample_(std::move(sample)),
      n0_(std::move(n0)),
      n1_(std::move(n1)),
      n_other_(std::move(n_other)),
      stats_(stats)
{
    if (n0_.size() != n1_.size() || n0_.size() != n_other_.size())
    {
        throw std::logic_error(
            std::format(
                "SampleCounts {}: count arrays differ in size ({}, {}, {})",
                sample_,
                n0_.size(),
                n1_.size(),
                n_other_.size()));
    }
}

auto SampleCounts::sites_covered() const -> size_t
{
    size_t n = 0;
    for (size_t s = 0; s < n0_.size(); ++s)
    {
        n += (n0_[s] + n1_[s] > 0) ? 1 : 0;
    }
    return n;
}

auto SampleCounts::bases_used() const -> uint64_t
{
    uint64_t n = 0;
    for (size_t s = 0; s < n0_.size(); ++s)
    {
        n += n0_[s] + n1_[s];
    }
    return n;
}
}  // namespace vcfbox::parentage
