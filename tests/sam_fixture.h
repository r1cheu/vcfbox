#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace vcfbox::test
{
struct SamContig
{
    std::string name;
    int64_t length = 0;
};

// One SAM record. Fields follow the SAM column order so designated
// initializers read like a SAM line.
struct SamRead
{
    std::string qname = "r";
    int flag = 0;
    std::string contig = "chr01";
    int64_t pos = 1;  // 1-based
    int mapq = 60;
    std::string cigar;
    std::string mate_contig = "*";
    int64_t mpos = 0;  // 1-based
    int64_t tlen = 0;
    std::string seq;
    std::string qual;
};

// Builds the text of a coordinate-sorted SAM. Defaults give the same two
// contigs as `PanelVcf`.
struct SamFile
{
    std::vector<SamContig> contigs{{"chr01", 1000}, {"chr02", 500}};
    std::vector<SamRead> reads;

    auto text() const -> std::string;
};

// `n` copies of `read`, with distinct query names.
auto repeated_reads(const SamRead& read, int n) -> std::vector<SamRead>;
}  // namespace vcfbox::test
