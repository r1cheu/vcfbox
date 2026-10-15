#include "parentage/sites.h"

#include <format>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace vcfbox::parentage
{
auto SiteTable::add_contig(ContigInfo info) -> void
{
    drop_empty_tail();
    for (const auto& c : contigs_)
    {
        if (c.info.name == info.name)
        {
            throw std::runtime_error(
                std::format(
                    "Contig {} appears twice; panel sites must be grouped by "
                    "contig",
                    info.name));
        }
    }
    contigs_.push_back({std::move(info), pos_.size(), 0});
}

auto SiteTable::add_site(
    int64_t pos,
    uint8_t allele0,
    uint8_t allele1,
    float freq1) -> void
{
    if (contigs_.empty())
    {
        throw std::logic_error("SiteTable::add_site called before add_contig");
    }
    ContigRange& contig = contigs_.back();
    if (pos > std::numeric_limits<int32_t>::max())
    {
        throw std::runtime_error(
            std::format(
                "{}:{}: position exceeds the supported range",
                contig.info.name,
                pos + 1));
    }
    if (contig.n_sites > 0 && pos <= pos_.back())
    {
        throw std::runtime_error(
            std::format(
                "{}:{}: panel positions must be strictly increasing within a "
                "contig",
                contig.info.name,
                pos + 1));
    }
    pos_.push_back(static_cast<int32_t>(pos));
    allele_.push_back(static_cast<uint8_t>((allele0 << 4) | allele1));
    freq1_.push_back(freq1);
    ++contig.n_sites;
}

auto SiteTable::finish() -> void
{
    drop_empty_tail();
    if (pos_.empty())
    {
        throw std::runtime_error("Panel contains no usable sites");
    }
}

auto SiteTable::drop_empty_tail() -> void
{
    if (!contigs_.empty() && contigs_.back().n_sites == 0)
    {
        contigs_.pop_back();
    }
}

auto SiteTable::contig_positions(size_t contig_index) const
    -> std::span<const int32_t>
{
    const auto& c = contigs_.at(contig_index);
    return std::span<const int32_t>(pos_).subspan(c.first_site, c.n_sites);
}

auto SiteTable::describe(size_t s) const -> std::string
{
    for (const auto& c : contigs_)
    {
        if (s >= c.first_site && s < c.first_site + c.n_sites)
        {
            return std::format(
                "{}:{} {}>{}",
                c.info.name,
                pos_[s] + 1,
                decode_base(allele0(s)),
                decode_base(allele1(s)));
        }
    }
    return std::format("site#{}", s);
}

auto match_contigs(
    const SiteTable& sites,
    std::span<const ContigInfo> bam_contigs,
    const std::string& bam_label) -> std::vector<int>
{
    std::unordered_map<std::string_view, int> by_name;
    by_name.reserve(bam_contigs.size());
    for (size_t i = 0; i < bam_contigs.size(); ++i)
    {
        by_name.emplace(bam_contigs[i].name, static_cast<int>(i));
    }

    std::vector<int> tids;
    tids.reserve(sites.contigs().size());
    for (const auto& c : sites.contigs())
    {
        const auto it = by_name.find(c.info.name);
        if (it == by_name.end())
        {
            throw std::runtime_error(
                std::format(
                    "Panel contig {} is not present in {}; the panel and the "
                    "BAM "
                    "must use the same reference",
                    c.info.name,
                    bam_label));
        }
        const auto& bc = bam_contigs[static_cast<size_t>(it->second)];
        if (bc.length != c.info.length)
        {
            throw std::runtime_error(
                std::format(
                    "Contig {} has length {} in the panel but {} in {}; the "
                    "panel "
                    "and the BAM must use the same reference",
                    c.info.name,
                    c.info.length,
                    bc.length,
                    bam_label));
        }
        tids.push_back(it->second);
    }
    return tids;
}
}  // namespace vcfbox::parentage
