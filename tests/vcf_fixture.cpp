#include "vcf_fixture.h"

#include <format>
#include <iterator>

namespace vcfbox::test
{
auto PanelVcf::text() const -> std::string
{
    std::string out;
    auto append = std::back_inserter(out);
    std::format_to(append, "##fileformat=VCFv4.2\n");
    for (const auto& c : contigs)
    {
        if (c.length > 0)
        {
            std::format_to(
                append, "##contig=<ID={},length={}>\n", c.name, c.length);
        }
        else
        {
            std::format_to(append, "##contig=<ID={}>\n", c.name);
        }
    }
    std::format_to(
        append,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    std::format_to(append, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if (!samples.empty())
    {
        std::format_to(append, "\tFORMAT");
    }
    for (const auto& s : samples)
    {
        std::format_to(append, "\t{}", s);
    }
    std::format_to(append, "\n");
    for (const auto& site : sites)
    {
        std::format_to(
            append,
            "{}\t{}\t.\t{}\t{}\t.\t.\t.",
            site.contig,
            site.pos,
            site.ref,
            site.alt);
        if (!samples.empty())
        {
            std::format_to(append, "\tGT");
        }
        for (const auto& gt : site.gts)
        {
            std::format_to(append, "\t{}", gt);
        }
        std::format_to(append, "\n");
    }
    return out;
}

auto repeated_sites(
    const std::string& contig,
    int n,
    const std::vector<std::string>& gts) -> std::vector<VcfSite>
{
    std::vector<VcfSite> sites;
    sites.reserve(static_cast<size_t>(n));
    for (int pos = 1; pos <= n; ++pos)
    {
        sites.push_back({contig, pos, "A", "G", gts});
    }
    return sites;
}
}  // namespace vcfbox::test
