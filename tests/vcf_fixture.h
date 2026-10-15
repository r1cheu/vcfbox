#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace vcfbox::test
{
struct VcfContig
{
    std::string name;
    int64_t length = 0;  // 0 omits `length=` from the header line
};

struct VcfSite
{
    std::string contig;
    int64_t pos = 0;  // 1-based, as written in the file
    std::string ref;
    std::string alt;
    std::vector<std::string> gts;  // one per sample, e.g. "0/0"
};

// Builds the text of a GT-only VCF. Defaults give two contigs and three
// samples so a test only lists the sites it cares about.
struct PanelVcf
{
    std::vector<VcfContig> contigs{{"chr01", 1000}, {"chr02", 500}};
    std::vector<std::string> samples{"P1", "P2", "P3"};
    std::vector<VcfSite> sites;

    auto text() const -> std::string;
};

// `n` A>G sites on `contig` at positions 1..n, all with genotypes `gts`.
auto repeated_sites(
    const std::string& contig,
    int n,
    const std::vector<std::string>& gts) -> std::vector<VcfSite>;
}  // namespace vcfbox::test
