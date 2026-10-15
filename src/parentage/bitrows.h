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

// Sites where the two genotype rows differ and neither parent is
// heterozygous.
inline auto popcount_discordant(
    std::span<const uint64_t> geno_a,
    std::span<const uint64_t> geno_b,
    std::span<const uint64_t> het_a,
    std::span<const uint64_t> het_b) -> size_t
{
    size_t n = 0;
    for (size_t w = 0; w < geno_a.size(); ++w)
    {
        n += static_cast<size_t>(
            std::popcount((geno_a[w] ^ geno_b[w]) & ~(het_a[w] | het_b[w])));
    }
    return n;
}
}  // namespace vcfbox::parentage
