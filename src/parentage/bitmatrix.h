#pragma once

#include <Eigen/Core>
#include <array>
#include <cstdint>
#include <span>
#include <string>

namespace vcfbox
{
inline constexpr std::array<char, 8>
    kBitmatrixMagic{'V', 'B', 'X', 'B', 'I', 'T', '0', '1'};

void write_bitmatrix(
    const std::string& path,
    std::span<const uint8_t> bits,
    uint64_t rows,
    uint64_t cols);

Eigen::MatrixXf load_bitmatrix(const std::string& path);

}  // namespace vcfbox
