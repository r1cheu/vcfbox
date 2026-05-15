#include "parentage/bitmatrix.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace vcfbox
{
void write_bitmatrix(
    const std::string& path,
    std::span<const uint8_t> bits,
    uint64_t rows,
    uint64_t cols)
{
    if (bits.size() != rows * cols)
    {
        throw std::runtime_error("write_bitmatrix size mismatch for: " + path);
    }

    std::ofstream out(path, std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open output file: " + path);
    }

    out.write(kBitmatrixMagic.data(), kBitmatrixMagic.size());
    out.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(cols));

    const uint64_t stride = (rows + 7) / 8;
    std::vector<uint8_t> col_bytes(stride);
    for (uint64_t col = 0; col < cols; ++col)
    {
        std::fill(col_bytes.begin(), col_bytes.end(), uint8_t{0});
        for (uint64_t row = 0; row < rows; ++row)
        {
            if (bits[(row * cols) + col] != 0)
            {
                col_bytes[row / 8] |= static_cast<uint8_t>(1U << (row % 8));
            }
        }
        out.write(
            reinterpret_cast<const char*>(col_bytes.data()),
            static_cast<std::streamsize>(stride));
    }
    if (!out)
    {
        throw std::runtime_error("Failed writing bitmatrix: " + path);
    }
}

Eigen::MatrixXd load_bitmatrix(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
    {
        throw std::runtime_error("Cannot open bitmatrix file: " + path);
    }

    std::array<char, 8> magic{};
    in.read(magic.data(), magic.size());
    if (!in || magic != kBitmatrixMagic)
    {
        throw std::runtime_error("Bad bitmatrix magic in: " + path);
    }

    uint64_t rows = 0;
    uint64_t cols = 0;
    in.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    in.read(reinterpret_cast<char*>(&cols), sizeof(cols));
    if (!in)
    {
        throw std::runtime_error("Bad bitmatrix header in: " + path);
    }

    Eigen::MatrixXd mat(
        static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
    const uint64_t stride = (rows + 7) / 8;
    std::vector<uint8_t> col_bytes(stride);
    for (uint64_t col = 0; col < cols; ++col)
    {
        in.read(
            reinterpret_cast<char*>(col_bytes.data()),
            static_cast<std::streamsize>(stride));
        if (!in)
        {
            throw std::runtime_error("Truncated bitmatrix in: " + path);
        }
        for (uint64_t row = 0; row < rows; ++row)
        {
            const uint32_t byte = col_bytes[row / 8];
            const uint32_t bit = (byte >> (row % 8)) & 1U;
            mat(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col))
                = static_cast<double>(bit);
        }
    }
    return mat;
}
}  // namespace vcfbox
