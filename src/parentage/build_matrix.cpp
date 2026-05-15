#include "parentage/parentage.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <regex>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "common/line_writer.h"
#include "common/progress.h"
#include "common/strings.h"
#include "hts/hts_raii.h"
#include "parentage/bitmatrix.h"

extern "C"
{
#include <htslib/vcf.h>
}

namespace
{
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

struct HomoBits
{
    uint8_t is_ref;
    uint8_t is_alt;
};

HomoBits classify_homo_bits(int32_t allele0, int32_t allele1)
{
    if (bcf_gt_is_missing(allele0) || bcf_gt_is_missing(allele1))
    {
        return {0, 0};
    }
    const int v0 = bcf_gt_allele(allele0);
    if (v0 != bcf_gt_allele(allele1))
    {
        return {0, 0};
    }
    if (v0 == 0)
    {
        return {1, 0};
    }
    if (v0 == 1)
    {
        return {0, 1};
    }
    return {0, 0};
}

void append_homo_bits(
    std::span<const int32_t> gt_arr,
    std::span<const int> indices,
    std::vector<uint8_t>& ref_bits,
    std::vector<uint8_t>& alt_bits)
{
    for (int idx : indices)
    {
        const auto bits
            = classify_homo_bits(gt_arr[(idx * 2)], gt_arr[(idx * 2) + 1]);
        ref_bits.push_back(bits.is_ref);
        alt_bits.push_back(bits.is_alt);
    }
}

}  // namespace

namespace vcfbox
{
void build_parentage_matrices(const ParentageMatrixOptions& options)
{
    hts_file_ptr vcf_file(bcf_open(options.parents_path.c_str(), "r"));
    if (!vcf_file)
    {
        throw std::runtime_error(
            "Could not open VCF file: " + options.parents_path);
    }
    bcf_hdr_ptr header(bcf_hdr_read(vcf_file.get()));
    if (!header)
    {
        throw std::runtime_error(
            "Could not read VCF header from: " + options.parents_path);
    }

    const auto groups = partition_samples(
        header.get(), options.maternal_patterns, options.paternal_patterns);

    std::vector<uint8_t> m0_bits;
    std::vector<uint8_t> m1_bits;
    std::vector<uint8_t> p0_bits;
    std::vector<uint8_t> p1_bits;

    LineWriter bed(options.matrix_prefix + ".sites.bed");

    bcf1_ptr rec(bcf_init());
    genotype_buffer gt;

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
        if (read_genotypes(header.get(), rec.get(), gt) <= 0)
        {
            throw std::runtime_error(
                "Missing GT field at "
                + std::string(bcf_hdr_id2name(header.get(), rec->rid)) + ":"
                + std::to_string(rec->pos + 1));
        }

        const std::span gt_span{
            gt.data.get(), static_cast<size_t>(gt.capacity)};
        append_homo_bits(gt_span, groups.maternal_idx, m0_bits, m1_bits);
        append_homo_bits(gt_span, groups.paternal_idx, p0_bits, p1_bits);

        const char* chrom = bcf_hdr_id2name(header.get(), rec->rid);
        const char* ref = rec->n_allele > 0 ? rec->d.allele[0] : ".";
        const char* alt = rec->n_allele > 1 ? rec->d.allele[1] : ".";
        bed.write_line(
            std::string(chrom) + "\t" + std::to_string(rec->pos) + "\t"
            + std::to_string(rec->pos + 1) + "\t" + ref + "\t" + alt);
    }
    counter->done();

    const uint64_t rows = processed;
    const uint64_t n_mat = groups.maternal_idx.size();
    const uint64_t n_pat = groups.paternal_idx.size();
    write_bitmatrix(options.matrix_prefix + ".M0.bin", m0_bits, rows, n_mat);
    write_bitmatrix(options.matrix_prefix + ".M1.bin", m1_bits, rows, n_mat);
    write_bitmatrix(options.matrix_prefix + ".P0.bin", p0_bits, rows, n_pat);
    write_bitmatrix(options.matrix_prefix + ".P1.bin", p1_bits, rows, n_pat);

    LineWriter mat_tsv(options.matrix_prefix + ".maternal.tsv");
    for (const auto& name : groups.maternal_names)
    {
        mat_tsv.write_line(name);
    }
    LineWriter pat_tsv(options.matrix_prefix + ".paternal.tsv");
    for (const auto& name : groups.paternal_names)
    {
        pat_tsv.write_line(name);
    }
}
}  // namespace vcfbox
