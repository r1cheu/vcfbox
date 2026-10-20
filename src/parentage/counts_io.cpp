#include "parentage/counts_io.h"

#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <format>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "hts/hts_raii.h"

extern "C"
{
#include <htslib/bgzf.h>
#include <htslib/hts_endian.h>
}

namespace vcfbox::parentage
{
namespace
{
constexpr std::string_view magic = "VBXCNT01";
constexpr size_t n_stat_fields = 7;

auto stats_to_array(const PileupStats& stats)
    -> std::array<uint64_t, n_stat_fields>
{
    return {
        stats.reads_total,
        stats.reads_used,
        stats.reads_skipped_flag,
        stats.reads_skipped_improper,
        stats.reads_skipped_mapq,
        stats.reads_off_panel,
        stats.bases_skipped_baseq};
}

auto stats_from_array(const std::array<uint64_t, n_stat_fields>& fields)
    -> PileupStats
{
    PileupStats stats;
    stats.reads_total = fields[0];
    stats.reads_used = fields[1];
    stats.reads_skipped_flag = fields[2];
    stats.reads_skipped_improper = fields[3];
    stats.reads_skipped_mapq = fields[4];
    stats.reads_off_panel = fields[5];
    stats.bases_skipped_baseq = fields[6];
    return stats;
}

class Writer
{
   public:
    explicit Writer(const std::string& path)
        : path_(path), file_(bgzf_open(path.c_str(), "w"))
    {
        if (!file_)
        {
            throw std::runtime_error(
                std::format("Failed to create counts file {}", path));
        }
    }

    auto bytes(std::span<const uint8_t> data) -> void
    {
        if (bgzf_write(file_.get(), data.data(), data.size())
            != static_cast<ssize_t>(data.size()))
        {
            throw std::runtime_error(
                std::format("Failed to write counts file {}", path_));
        }
    }

    auto u64(uint64_t value) -> void
    {
        std::array<uint8_t, 8> buf{};
        u64_to_le(value, buf.data());
        bytes(buf);
    }

    auto text(std::string_view s) -> void
    {
        bytes(std::span(reinterpret_cast<const uint8_t*>(s.data()), s.size()));
    }

    // Flushes and closes the file; throws when that fails.
    auto close() -> void
    {
        if (bgzf_close(file_.release()) != 0)
        {
            throw std::runtime_error(
                std::format("Failed to write counts file {}", path_));
        }
    }

   private:
    std::string path_;
    bgzf_ptr file_;
};

class Reader
{
   public:
    explicit Reader(const std::string& path)
        : path_(path), file_(bgzf_open(path.c_str(), "r"))
    {
        if (!file_)
        {
            throw std::runtime_error(
                std::format("Failed to open counts file {}", path));
        }
    }

    auto bytes(std::span<uint8_t> data) -> void
    {
        if (bgzf_read(file_.get(), data.data(), data.size())
            != static_cast<ssize_t>(data.size()))
        {
            throw std::runtime_error(
                std::format("Counts file {} is truncated", path_));
        }
    }

    auto u64() -> uint64_t
    {
        std::array<uint8_t, 8> buf{};
        bytes(buf);
        return le_to_u64(buf.data());
    }

    auto text(size_t n) -> std::string
    {
        std::string s(n, '\0');
        bytes(std::span(reinterpret_cast<uint8_t*>(s.data()), s.size()));
        return s;
    }

    auto array(size_t n) -> std::vector<uint8_t>
    {
        std::vector<uint8_t> v(n);
        bytes(v);
        return v;
    }

   private:
    std::string path_;
    bgzf_ptr file_;
};
}  // namespace

auto write_counts(const SampleCounts& counts, const std::string& path) -> void
{
    const std::string partial = std::format("{}.partial", path);
    Writer out(partial);
    out.text(magic);
    out.u64(counts.sample().size());
    out.text(counts.sample());
    out.u64(counts.n_sites());
    for (const uint64_t field : stats_to_array(counts.stats()))
    {
        out.u64(field);
    }
    out.bytes(counts.n0());
    out.bytes(counts.n1());
    out.bytes(counts.n_other());
    out.close();
    std::filesystem::rename(partial, path);
}

auto read_counts(const std::string& path, size_t expected_n_sites)
    -> SampleCounts
{
    Reader in(path);
    if (in.text(magic.size()) != magic)
    {
        throw std::runtime_error(
            std::format("{} is not a vcfbox counts file", path));
    }
    std::string sample = in.text(in.u64());
    const uint64_t n_sites = in.u64();
    if (n_sites != expected_n_sites)
    {
        throw std::runtime_error(
            std::format(
                "{} holds {} sites but the panel has {}; it was made with a "
                "different panel",
                path,
                n_sites,
                expected_n_sites));
    }
    std::array<uint64_t, n_stat_fields> fields{};
    for (uint64_t& field : fields)
    {
        field = in.u64();
    }
    std::vector<uint8_t> n0 = in.array(n_sites);
    std::vector<uint8_t> n1 = in.array(n_sites);
    std::vector<uint8_t> n_other = in.array(n_sites);
    return {
        std::move(sample),
        std::move(n0),
        std::move(n1),
        std::move(n_other),
        stats_from_array(fields)};
}
}  // namespace vcfbox::parentage
