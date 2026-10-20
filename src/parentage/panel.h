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

// Parent panel. Row k of `geno()` is parent k; bit s is 1 when the parent
// carries allele1 at site s and 0 for allele0. Bit s of `het()` is 1 where
// the parent is heterozygous (residual heterozygosity of an inbred line) and
// bit s of `missing()` where it has no call; its `geno` bit is then 0 and
// carries no information.
class Panel
{
   public:
    // Takes ownership of the parallel containers and records every parent
    // pair differing on at most `duplicate_max_discordance` of the sites
    // where both are homozygous. Throws when `names`, `geno`, `het`,
    // `missing` and `sites` disagree on the number of parents or sites.
    Panel(
        SiteTable sites,
        std::vector<std::string> names,
        BitRows geno,
        BitRows het,
        BitRows missing,
        double duplicate_max_discordance,
        size_t n_excluded_sites);

    auto n_parents() const -> size_t { return names_.size(); }
    auto n_sites() const -> size_t { return sites_.size(); }
    auto n_excluded_sites() const -> size_t { return n_excluded_sites_; }
    // Heterozygous calls over all parents and sites.
    auto n_het_calls() const -> size_t { return n_het_calls_; }
    // Missing calls over all parents and sites.
    auto n_missing_calls() const -> size_t { return n_missing_calls_; }

    auto sites() const -> const SiteTable& { return sites_; }
    auto names() const -> std::span<const std::string> { return names_; }
    auto geno() const -> const BitRows& { return geno_; }
    auto het() const -> const BitRows& { return het_; }
    auto missing() const -> const BitRows& { return missing_; }
    auto parent_rows(size_t parent) const -> ParentRows
    {
        return {geno_.row(parent), het_.row(parent), missing_.row(parent)};
    }
    auto near_duplicates() const -> std::span<const NearDuplicate>
    {
        return near_duplicates_;
    }

   private:
    SiteTable sites_;
    std::vector<std::string> names_;
    BitRows geno_;
    BitRows het_;
    BitRows missing_;
    std::vector<NearDuplicate> near_duplicates_;
    size_t n_excluded_sites_ = 0;
    size_t n_het_calls_ = 0;
    size_t n_missing_calls_ = 0;
};

// Loads a bi-allelic SNP panel. Every genotype must be diploid; a genotype
// with a missing allele (`./.`, `./1`) is a missing call. Haploid calls,
// multi-allelic sites, indels, unsorted or duplicated positions and sites
// where no parent is called are errors. `freq1` of a site is the allele1
// dosage over the called parents.
[[nodiscard]] auto load_panel(
    const std::string& vcf_path,
    const PanelOptions& options) -> Panel;
}  // namespace vcfbox::parentage
