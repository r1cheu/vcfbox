#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace vcfbox::parentage
{
// Dense bit matrix with `rows` rows of `n_bits` bits each. Rows are stored
// contiguously as 64-bit words so a row can be fed to AND/XOR + popcount
// kernels directly.
class BitRows
{
   public:
    BitRows() = default;
    BitRows(size_t rows, size_t n_bits)
        : rows_(rows),
          n_bits_(n_bits),
          words_per_row_((n_bits + 63) / 64),
          bits_(rows * words_per_row_, 0)
    {
    }

    auto rows() const -> size_t { return rows_; }
    auto n_bits() const -> size_t { return n_bits_; }
    auto words_per_row() const -> size_t { return words_per_row_; }

    auto row(size_t r) -> std::span<uint64_t>
    {
        return std::span(bits_).subspan(r * words_per_row_, words_per_row_);
    }
    auto row(size_t r) const -> std::span<const uint64_t>
    {
        return std::span(bits_).subspan(r * words_per_row_, words_per_row_);
    }

    auto set(size_t r, size_t bit) -> void
    {
        row(r)[bit / 64] |= uint64_t{1} << (bit % 64);
    }
    auto get(size_t r, size_t bit) const -> bool
    {
        return ((row(r)[bit / 64] >> (bit % 64)) & 1) != 0;
    }

   private:
    size_t rows_ = 0;
    size_t n_bits_ = 0;
    size_t words_per_row_ = 0;
    std::vector<uint64_t> bits_;
};

inline auto popcount_rows(std::span<const uint64_t> a) -> size_t
{
    size_t n = 0;
    for (const uint64_t word : a)
    {
        n += static_cast<size_t>(std::popcount(word));
    }
    return n;
}

// Bit planes of one panel parent: `geno` is 1 for allele1, `het` for a
// heterozygous call and `missing` for no call; `geno` is 0 wherever `het` or
// `missing` is set.
struct ParentRows
{
    std::span<const uint64_t> geno;
    std::span<const uint64_t> het;
    std::span<const uint64_t> missing;
};

// Bits where both parents are homozygous.
inline auto both_homozygous(const ParentRows& a, const ParentRows& b, size_t w)
    -> uint64_t
{
    return ~(a.het[w] | a.missing[w] | b.het[w] | b.missing[w]);
}

// Sites where both parents are homozygous and carry different alleles.
inline auto popcount_discordant(const ParentRows& a, const ParentRows& b)
    -> size_t
{
    size_t n = 0;
    for (size_t w = 0; w < a.geno.size(); ++w)
    {
        n += static_cast<size_t>(
            std::popcount((a.geno[w] ^ b.geno[w]) & both_homozygous(a, b, w)));
    }
    return n;
}
}  // namespace vcfbox::parentage
