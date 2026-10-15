#include "sam_fixture.h"

#include <format>
#include <iterator>

namespace vcfbox::test
{
auto SamFile::text() const -> std::string
{
    std::string out;
    auto append = std::back_inserter(out);
    std::format_to(append, "@HD\tVN:1.6\tSO:coordinate\n");
    for (const auto& c : contigs)
    {
        std::format_to(append, "@SQ\tSN:{}\tLN:{}\n", c.name, c.length);
    }
    for (const auto& r : reads)
    {
        std::format_to(
            append,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.qname,
            r.flag,
            r.contig,
            r.pos,
            r.mapq,
            r.cigar,
            r.mate_contig,
            r.mpos,
            r.tlen,
            r.seq,
            r.qual);
    }
    return out;
}

auto repeated_reads(const SamRead& read, int n) -> std::vector<SamRead>
{
    std::vector<SamRead> reads(static_cast<size_t>(n), read);
    for (int i = 0; i < n; ++i)
    {
        reads[static_cast<size_t>(i)].qname
            = std::format("{}_{}", read.qname, i);
    }
    return reads;
}
}  // namespace vcfbox::test
