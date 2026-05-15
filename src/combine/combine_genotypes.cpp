#include "combine/combine_genotypes.h"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "common/path.h"
#include "common/progress.h"
#include "hts/hts_raii.h"

extern "C"
{
#include <htslib/vcf.h>
}

namespace
{
using SampleIndex = std::unordered_map<std::string, int>;

SampleIndex build_sample_index(bcf_hdr_t* header)
{
    SampleIndex index;
    const int n = bcf_hdr_nsamples(header);
    index.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        index.emplace(header->samples[i], i);
    }
    return index;
}

void check_pairs_in_header(
    const std::vector<vcfbox::SamplePair>& sample_pairs,
    const SampleIndex& sample_to_idx)
{
    std::vector<std::string> missing;
    for (const auto& pair : sample_pairs)
    {
        if (!sample_to_idx.contains(pair.first))
        {
            missing.push_back(pair.first);
        }
        if (!sample_to_idx.contains(pair.second))
        {
            missing.push_back(pair.second);
        }
    }
    if (missing.empty())
    {
        return;
    }
    std::string list;
    for (size_t i = 0; i < missing.size(); ++i)
    {
        if (i != 0)
        {
            list += ", ";
        }
        list += missing[i];
    }
    throw std::runtime_error(
        "Samples not found in VCF: " + list
        + ". Make sure sample list matches VCF file.");
}

bcf_hdr_t* init_bcf_header(
    bcf_hdr_t* header,
    const std::vector<vcfbox::SamplePair>& sample_pairs,
    bool keep_old_samples)
{
    bcf_hdr_t* output_header = bcf_hdr_init("w");

    for (int i = 0; i < header->nhrec; ++i)
    {
        if (header->hrec[i]->type == BCF_HL_CTG)
        {
            bcf_hdr_add_hrec(output_header, bcf_hrec_dup(header->hrec[i]));
        }
    }

    bcf_hrec_t* gt_hrec
        = bcf_hdr_get_hrec(header, BCF_HL_FMT, "ID", "GT", nullptr);
    if (gt_hrec != nullptr)
    {
        bcf_hdr_add_hrec(output_header, bcf_hrec_dup(gt_hrec));
    }
    else
    {
        bcf_hdr_append(
            output_header,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    }

    if (keep_old_samples)
    {
        for (int i = 0; i < bcf_hdr_nsamples(header); ++i)
        {
            bcf_hdr_add_sample(output_header, header->samples[i]);
        }
    }

    for (const auto& pair : sample_pairs)
    {
        bcf_hdr_add_sample(
            output_header, (pair.first + "~" + pair.second).c_str());
    }

    bcf_hdr_add_sample(output_header, nullptr);

    return output_header;
}

void copy_rec_info(
    bcf_hdr_t* header,
    bcf_hdr_t* output_header,
    bcf1_t* in_rec,
    bcf1_t* out_rec)
{
    bcf_clear(out_rec);
    out_rec->rid
        = bcf_hdr_name2id(output_header, bcf_hdr_id2name(header, in_rec->rid));
    out_rec->pos = in_rec->pos;
    bcf_update_id(output_header, out_rec, in_rec->d.id);
    bcf_update_alleles(
        output_header,
        out_rec,
        const_cast<const char**>(in_rec->d.allele),
        in_rec->n_allele);

    out_rec->qual = in_rec->qual;
}

std::vector<int32_t> concat_gt(
    const std::vector<vcfbox::SamplePair>& sample_pairs,
    const SampleIndex& sample_to_idx,
    const int32_t* gt_arr,
    bool keep_old_samples,
    int n_gt)
{
    std::vector<int32_t> out_gts;
    if (keep_old_samples)
    {
        out_gts.reserve(n_gt + (sample_pairs.size() * 2));
        for (int i = 0; i < n_gt / 2; ++i)
        {
            int32_t gt0 = gt_arr[i * 2];
            int32_t gt1 = gt_arr[i * 2 + 1];
            if (bcf_gt_is_missing(gt0)
                || bcf_gt_allele(gt0) != bcf_gt_allele(gt1))
            {
                out_gts.push_back(bcf_gt_missing);
                out_gts.push_back(bcf_gt_missing);
            }
            else
            {
                out_gts.push_back(gt0);
                out_gts.push_back(gt1);
            }
        }
    }
    else
    {
        out_gts.reserve(sample_pairs.size() * 2);
    }
    for (const auto& pair : sample_pairs)
    {
        int f_idx = sample_to_idx.at(pair.first);
        int m_idx = sample_to_idx.at(pair.second);

        int32_t f_gt0 = gt_arr[f_idx * 2];
        int32_t f_gt1 = gt_arr[f_idx * 2 + 1];
        int32_t m_gt0 = gt_arr[m_idx * 2];
        int32_t m_gt1 = gt_arr[m_idx * 2 + 1];

        if (bcf_gt_is_missing(f_gt0) || bcf_gt_is_missing(f_gt1)
            || bcf_gt_is_missing(m_gt0) || bcf_gt_is_missing(m_gt1)
            || bcf_gt_allele(f_gt0) != bcf_gt_allele(f_gt1)
            || bcf_gt_allele(m_gt0) != bcf_gt_allele(m_gt1))
        {
            out_gts.push_back(bcf_gt_missing);
            out_gts.push_back(bcf_gt_missing);
        }
        else
        {
            out_gts.push_back(bcf_gt_unphased(bcf_gt_allele(f_gt0)));
            out_gts.push_back(bcf_gt_unphased(bcf_gt_allele(m_gt0)));
        }
    }
    return out_gts;
}
}  // namespace

namespace vcfbox
{
std::vector<SamplePair> parse_sample_pairs(const std::string& file_path)
{
    std::vector<SamplePair> pairs;
    std::ifstream file(file_path);

    if (!file)
    {
        throw std::runtime_error("Cannot open sample pairs file: " + file_path);
    }
    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        std::istringstream iss(line);
        std::string s1;
        std::string s2;
        if (iss >> s1 >> s2)
        {
            pairs.emplace_back(s1, s2);
        }
    }
    return pairs;
}

void combine_genotypes(
    const std::string& vcf_path,
    const std::vector<SamplePair>& sample_pairs,
    bool keep_old_samples,
    const std::string& out_path)
{
    HtsFile vcf_file(bcf_open(vcf_path.c_str(), "r"));
    if (!vcf_file)
    {
        throw std::runtime_error("Could not open VCF file: " + vcf_path);
    }
    BcfHdr header(bcf_hdr_read(vcf_file.get()));
    if (!header)
    {
        throw std::runtime_error("Could not read VCF header from: " + vcf_path);
    }

    auto sample_to_idx = build_sample_index(header.get());
    check_pairs_in_header(sample_pairs, sample_to_idx);

    const auto mode = parse_mode(out_path);
    HtsFile output_file(hts_open(out_path.c_str(), mode.c_str()));
    if (!output_file)
    {
        throw std::runtime_error("Could not open output file: " + out_path);
    }

    BcfHdr output_header(
        init_bcf_header(header.get(), sample_pairs, keep_old_samples));

    if (bcf_hdr_write(output_file.get(), output_header.get()) != 0)
    {
        throw std::runtime_error("Failed to write output header");
    }

    BcfRec in_rec(bcf_init());
    BcfRec out_rec(bcf_init());
    Genotypes gt;

    size_t processed_snp = 0;
    auto counter = create_counter("Adding SNPs", processed_snp);
    counter->show();
    while (bcf_read(vcf_file.get(), header.get(), in_rec.get()) == 0)
    {
        processed_snp++;
        bcf_unpack(in_rec.get(), BCF_UN_ALL);
        if (in_rec->n_allele > 2)
        {
            continue;
        }
        if (bcf_get_genotypes(header.get(), in_rec.get(), &gt.p_, &gt.n_) <= 0)
        {
            continue;
        }

        auto out_gts = concat_gt(
            sample_pairs, sample_to_idx, gt.p_, keep_old_samples, gt.n_);
        copy_rec_info(
            header.get(), output_header.get(), in_rec.get(), out_rec.get());
        bcf_update_genotypes(
            output_header.get(), out_rec.get(), out_gts.data(), out_gts.size());

        if (bcf_write(output_file.get(), output_header.get(), out_rec.get())
            != 0)
        {
            throw std::runtime_error("Failed to write VCF record");
        }
    }
    counter->done();
}
}  // namespace vcfbox
