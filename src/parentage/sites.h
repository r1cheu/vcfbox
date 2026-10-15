#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace vcfbox::parentage
{
// Base codes: A=0 C=1 G=2 T=3; `no_base` for anything else.
inline constexpr uint8_t no_base = 0xF;

namespace detail
{
consteval auto make_base_table() -> std::array<uint8_t, 256>
{
    std::array<uint8_t, 256> table{};
    table.fill(no_base);
    table['A'] = table['a'] = 0;
    table['C'] = table['c'] = 1;
    table['G'] = table['g'] = 2;
    table['T'] = table['t'] = 3;
    return table;
}
inline constexpr auto base_table = make_base_table();
}  // namespace detail

constexpr auto encode_base(char base) -> uint8_t
{
    return detail::base_table[static_cast<unsigned char>(base)];
}

constexpr auto decode_base(uint8_t code) -> char
{
    constexpr std::string_view bases = "ACGT";
    return code < bases.size() ? bases[code] : 'N';
}

struct ContigInfo
{
    std::string name;
    int64_t length = 0;
};

// A contig's slice of the flat site arrays: sites
// [first_site, first_site + n_sites) belong to `info`.
struct ContigRange
{
    ContigInfo info;
    size_t first_site = 0;
    size_t n_sites = 0;
};

// The single source of site order for the whole module: panel rows, sample
// counts and site reports are all indexed by position in this table.
// Sites are grouped by contig in insertion order and strictly increasing in
// position within a contig; `add_contig` / `add_site` enforce both.
class SiteTable
{
   public:
    // Starts a new contig; following `add_site` calls belong to it. A
    // preceding contig that received no sites is dropped. Throws on a
    // duplicate contig name.
    auto add_contig(ContigInfo info) -> void;

    // Appends a site to the current contig. Throws when no contig has been
    // started, when `pos` does not increase, or when it exceeds int32.
    auto add_site(int64_t pos, uint8_t allele0, uint8_t allele1, float freq1)
        -> void;

    // Drops a trailing empty contig. Throws when the table holds no sites.
    auto finish() -> void;

    auto size() const -> size_t { return pos_.size(); }
    auto contigs() const -> std::span<const ContigRange> { return contigs_; }
    auto positions() const -> std::span<const int32_t> { return pos_; }
    auto freqs1() const -> std::span<const float> { return freq1_; }

    auto pos(size_t s) const -> int32_t { return pos_[s]; }
    auto allele0(size_t s) const -> uint8_t { return allele_[s] >> 4; }
    auto allele1(size_t s) const -> uint8_t { return allele_[s] & 0xF; }
    auto freq1(size_t s) const -> float { return freq1_[s]; }

    auto contig_positions(size_t contig_index) const
        -> std::span<const int32_t>;

    // "chr01:12345 A>G" for messages.
    auto describe(size_t s) const -> std::string;

   private:
    auto drop_empty_tail() -> void;

    std::vector<ContigRange> contigs_;
    std::vector<int32_t> pos_;     // 0-based
    std::vector<uint8_t> allele_;  // high nibble allele0, low nibble allele1
    std::vector<float> freq1_;     // frequency of allele1 among parents
};

// Maps every panel contig to the matching BAM contig (by name) and checks
// that the lengths agree. Throws when a contig is absent or has a different
// length. Returns the BAM tid for each panel contig.
auto match_contigs(
    const SiteTable& sites,
    std::span<const ContigInfo> bam_contigs,
    const std::string& bam_label) -> std::vector<int>;
}  // namespace vcfbox::parentage
