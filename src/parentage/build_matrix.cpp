#include "parentage/parentage.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <regex>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "common/progress.h"
#include "common/strings.h"
#include "hts/hts_raii.h"

extern "C"
{
#include <htslib/vcf.h>
}

namespace
{
constexpr std::array<char, 8>
    kBitmatrixMagic{'V', 'B', 'X', 'B', 'I', 'T', '0', '1'};

std::vector<std::regex> compile_patterns(std::span<const std::string> patterns)
{
    std::vector<std::regex> compiled;
    compiled.reserve(patterns.size());
    for (const auto& p : patterns)
    {
        try
        {
            compiled.emplace_back(p);
        }
        catch (const std::regex_error& e)
        {
            throw std::runtime_error(
                "Invalid regex pattern '" + p + "': " + e.what());
        }
    }
    return compiled;
}

bool matches_any(
    const std::string& sample,
    std::span<const std::regex> patterns)
{
    return std::ranges::any_of(
        patterns,
        [&](const std::regex& r) { return std::regex_match(sample, r); });
}

struct SampleGroups
{
    std::vector<int> maternal_idx;
    std::vector<int> paternal_idx;
    std::vector<std::string> maternal_names;
    std::vector<std::string> paternal_names;
};

SampleGroups partition_samples(
    bcf_hdr_t* header,
    std::span<const std::string> maternal_patterns,
    std::span<const std::string> paternal_patterns)
{
    auto maternal_regex = compile_patterns(maternal_patterns);
    auto paternal_regex = compile_patterns(paternal_patterns);

    SampleGroups g;
    std::vector<std::string> conflicts;
    const int n = bcf_hdr_nsamples(header);
    for (int i = 0; i < n; ++i)
    {
        std::string name = header->samples[i];
        const bool in_m = matches_any(name, maternal_regex);
        const bool in_p = matches_any(name, paternal_regex);
        if (in_m && in_p)
        {
            conflicts.push_back(name);
            continue;
        }
        if (in_m)
        {
            g.maternal_idx.push_back(i);
            g.maternal_names.push_back(std::move(name));
        }
        else if (in_p)
        {
            g.paternal_idx.push_back(i);
            g.paternal_names.push_back(std::move(name));
        }
    }
    if (!conflicts.empty())
    {
        throw std::runtime_error(
            "Samples matched by both maternal and paternal patterns: "
            + vcfbox::join(conflicts, ", "));
    }
    if (g.maternal_idx.empty())
    {
        throw std::runtime_error(
            "No samples matched maternal patterns: "
            + vcfbox::join(maternal_patterns, ", "));
    }
    if (g.paternal_idx.empty())
    {
        throw std::runtime_error(
            "No samples matched paternal patterns: "
            + vcfbox::join(paternal_patterns, ", "));
    }
    return g;
}

enum class HomoState : uint8_t
{
    Other = 0,
    Ref = 1,
    Alt = 2,
};

HomoState classify_homo(int32_t allele0, int32_t allele1)
{
    if (bcf_gt_is_missing(allele0) || bcf_gt_is_missing(allele1))
    {
        return HomoState::Other;
    }
    const int v0 = bcf_gt_allele(allele0);
    if (v0 != bcf_gt_allele(allele1))
    {
        return HomoState::Other;
    }
    if (v0 == 0)
    {
        return HomoState::Ref;
    }
    if (v0 == 1)
    {
        return HomoState::Alt;
    }
    return HomoState::Other;
}

void append_states(
    std::span<const int32_t> gt_arr,
    std::span<const int> indices,
    std::vector<HomoState>& out)
{
    out.reserve(out.size() + indices.size());
    for (int idx : indices)
    {
        out.push_back(classify_homo(gt_arr[(idx * 2)], gt_arr[(idx * 2) + 1]));
    }
}

void write_bitmatrix_eq(
    const std::string& path,
    std::span<const HomoState> states,
    uint64_t rows,
    uint64_t cols,
    HomoState target)
{
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
            if (states[(row * cols) + col] == target)
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

void write_lines(const std::string& path, std::span<const std::string> lines)
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("Cannot open output file: " + path);
    }
    for (const auto& l : lines)
    {
        out << l << '\n';
    }
    if (!out)
    {
        throw std::runtime_error("Failed writing: " + path);
    }
}
}  // namespace

namespace vcfbox
{
void build_parentage_matrices(const ParentageMatrixOptions& options)
{
    HtsFile vcf_file(bcf_open(options.parents_path.c_str(), "r"));
    if (!vcf_file)
    {
        throw std::runtime_error(
            "Could not open VCF file: " + options.parents_path);
    }
    BcfHdr header(bcf_hdr_read(vcf_file.get()));
    if (!header)
    {
        throw std::runtime_error(
            "Could not read VCF header from: " + options.parents_path);
    }

    const auto groups = partition_samples(
        header.get(), options.maternal_patterns, options.paternal_patterns);

    std::vector<HomoState> maternal_states;
    std::vector<HomoState> paternal_states;
    std::vector<std::string> site_lines;

    BcfRec rec(bcf_init());
    Genotypes gt;

    size_t processed = 0;
    auto counter = create_counter("Indicating sites", processed);
    counter->show();

    while (bcf_read(vcf_file.get(), header.get(), rec.get()) == 0)
    {
        processed++;
        bcf_unpack(rec.get(), BCF_UN_ALL);
        if (rec->n_allele > 2)
        {
            throw std::runtime_error(
                "Multi-allelic site at "
                + std::string(bcf_hdr_id2name(header.get(), rec->rid)) + ":"
                + std::to_string(rec->pos + 1)
                + " (filter to bi-allelic before build-matrix)");
        }
        if (bcf_get_genotypes(header.get(), rec.get(), &gt.p_, &gt.n_) <= 0)
        {
            throw std::runtime_error(
                "Missing GT field at "
                + std::string(bcf_hdr_id2name(header.get(), rec->rid)) + ":"
                + std::to_string(rec->pos + 1));
        }

        const std::span gt_span{gt.p_, static_cast<size_t>(gt.n_)};
        append_states(gt_span, groups.maternal_idx, maternal_states);
        append_states(gt_span, groups.paternal_idx, paternal_states);

        const char* chrom = bcf_hdr_id2name(header.get(), rec->rid);
        const char* ref = rec->n_allele > 0 ? rec->d.allele[0] : ".";
        const char* alt = rec->n_allele > 1 ? rec->d.allele[1] : ".";
        site_lines.push_back(
            std::string(chrom) + "\t" + std::to_string(rec->pos) + "\t"
            + std::to_string(rec->pos + 1) + "\t" + ref + "\t" + alt);
    }
    counter->done();

    const uint64_t rows = site_lines.size();
    const uint64_t n_mat = groups.maternal_idx.size();
    const uint64_t n_pat = groups.paternal_idx.size();
    write_bitmatrix_eq(
        options.matrix_prefix + ".M0.bin",
        maternal_states,
        rows,
        n_mat,
        HomoState::Ref);
    write_bitmatrix_eq(
        options.matrix_prefix + ".M1.bin",
        maternal_states,
        rows,
        n_mat,
        HomoState::Alt);
    write_bitmatrix_eq(
        options.matrix_prefix + ".P0.bin",
        paternal_states,
        rows,
        n_pat,
        HomoState::Ref);
    write_bitmatrix_eq(
        options.matrix_prefix + ".P1.bin",
        paternal_states,
        rows,
        n_pat,
        HomoState::Alt);
    write_lines(options.matrix_prefix + ".sites.bed", site_lines);
    write_lines(options.matrix_prefix + ".maternal.tsv", groups.maternal_names);
    write_lines(options.matrix_prefix + ".paternal.tsv", groups.paternal_names);
}
}  // namespace vcfbox
