#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "parentage/bitrows.h"
#include "parentage/sites.h"

namespace vcfbox::parentage
{
struct PanelOptions
{
    // Optional file of `contig<TAB>pos` (1-based) sites to drop at load.
    std::string exclude_sites_path;
    // Parent pairs whose genotypes differ on at most this fraction of sites
    // are reported as near-duplicates. They are never merged.
    double duplicate_max_discordance = 0.001;
};

// A pair of parents whose genotypes differ on `n_discordant` sites.
struct NearDuplicate
{
    uint32_t parent_a = 0;
    uint32_t parent_b = 0;
    size_t n_discordant = 0;
};

// Homozygous, fully called parent panel. Row k of `geno()` is parent k; bit s
// is 1 when the parent carries allele1 at site s and 0 for allele0.
class Panel
{
   public:
    // Takes ownership of the three parallel containers and records every
    // parent pair differing on at most `duplicate_max_discordance` of the
    // sites. Throws when `names`, `geno` and `sites` disagree on the number
    // of parents or sites.
    Panel(
        SiteTable sites,
        std::vector<std::string> names,
        BitRows geno,
        double duplicate_max_discordance,
        size_t n_excluded_sites);

    auto n_parents() const -> size_t { return names_.size(); }
    auto n_sites() const -> size_t { return sites_.size(); }
    auto n_excluded_sites() const -> size_t { return n_excluded_sites_; }

    auto sites() const -> const SiteTable& { return sites_; }
    auto names() const -> std::span<const std::string> { return names_; }
    auto geno() const -> const BitRows& { return geno_; }
    auto near_duplicates() const -> std::span<const NearDuplicate>
    {
        return near_duplicates_;
    }

   private:
    SiteTable sites_;
    std::vector<std::string> names_;
    BitRows geno_;
    std::vector<NearDuplicate> near_duplicates_;
    size_t n_excluded_sites_ = 0;
};

// Loads a bi-allelic SNP panel. Every genotype must be homozygous and
// called; anything else (missing, heterozygous, multi-allelic, indel,
// unsorted or duplicated positions) is an error.
[[nodiscard]] auto load_panel(
    const std::string& vcf_path,
    const PanelOptions& options) -> Panel;
}  // namespace vcfbox::parentage
