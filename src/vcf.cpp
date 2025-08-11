#include "vcf.h"

#include <cstddef>
#include <cstdlib>
#include <format>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "utils.h"
#include "vcf_raii.h"
namespace bk = barkeep;
namespace vcfbox
{
std::vector<detail::SamplePair> parse_sample_pairs(const std::string& file_path)
{
    std::vector<detail::SamplePair> pairs;
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
    const std::vector<detail::SamplePair>& sample_pairs,
    bool keep_old_samples,
    const std::string& out_path,
    const std::string& mode)
{
    detail::check_sample_consistence(vcf_path, sample_pairs);
    size_t n_lines = vcfbox::count_records(vcf_path);

    HtsFile vcf_file(bcf_open(vcf_path.c_str(), "r"));
    BcfHdr header(bcf_hdr_read(vcf_file.get()));
    std::unordered_map<std::string, int> sample_to_idx;
    for (int i = 0; i < bcf_hdr_nsamples(header); ++i)
    {
        sample_to_idx[header->samples[i]] = i;
    }

    HtsFile output_file(hts_open(out_path.c_str(), mode.c_str()));
    if (!output_file)
    {
        throw std::runtime_error("Could not open output file: " + out_path);
    }

    BcfHdr output_header(
        detail::init_bcf_head(header.get(), sample_pairs, keep_old_samples));

    if (bcf_hdr_write(output_file.get(), output_header.get()) != 0)
    {
        throw std::runtime_error("Failed to write output header");
    }

    BcfRec in_rec(bcf_init());
    BcfRec out_rec(bcf_init());
    Genotypes gt;

    size_t processd_snp = 0;
    auto bar = detail::create_progress(n_lines, processd_snp);
    bar->show();
    while (bcf_read(vcf_file.get(), header.get(), in_rec.get()) == 0)
    {
        processd_snp++;
        bcf_unpack(in_rec.get(), BCF_UN_ALL);
        if (in_rec->n_allele > 2)
        {
            continue;
        }
        if (bcf_get_genotypes(header.get(), in_rec.get(), &gt.p_, &gt.n_) <= 0)
        {
            continue;
        }

        auto out_gts = detail::concat_gt(
            sample_pairs, sample_to_idx, gt.p_, keep_old_samples, gt.n_);
        detail::copy_rec_info(
            header.get(), output_header.get(), in_rec.get(), out_rec.get());
        bcf_update_genotypes(
            output_header.get(), out_rec.get(), out_gts.data(), out_gts.size());

        if (bcf_write(output_file.get(), output_header.get(), out_rec.get())
            != 0)
        {
            throw std::runtime_error("Failed to write VCF record");
        }
    }
    bar->done();
}

void to_hapmap(const std::string& vcf_path, const std::string& out_path)
{
    HtsFile vcf_file(bcf_open(vcf_path.c_str(), "r"));
    BcfHdr header(bcf_hdr_read(vcf_file.get()));
    if (!header)
    {
        throw std::runtime_error("Failed to read VCF header");
    }

    std::ofstream stream(out_path);
    if (!stream)
    {
        throw std::runtime_error("Failed to open output file: " + out_path);
    }

    stream << "rs\talleles\tchrom\tpos\tstrand\t"
              "assembly\tcenter\tprotLSID\tassayLSID\t"
              "panel\tQCcode\t";

    for (int i = 0; i < bcf_hdr_nsamples(header); ++i)
    {
        stream << header->samples[i] << "\t";
    }
    stream << "\n";
    BcfRec in_rec(bcf_init());
    Genotypes gt;
    size_t processd_snp = 0;
    auto counter
        = detail::create_counter("Converting to HapMap format", processd_snp);

    counter->show();
    while (bcf_read(vcf_file.get(), header.get(), in_rec.get()) == 0)
    {
        processd_snp++;
        bcf_unpack(in_rec.get(), BCF_UN_ALL);
        if (in_rec->n_allele > 2)
        {
            continue;
        }

        std::string ref = in_rec->d.allele[0];
        std::string alt = in_rec->d.allele[1];
        std::string chrom = std::format("chr{:02d}", in_rec->rid + 1);
        int32_t pos = in_rec->pos + 1;
        std::string rs = std::format("{}_{:d}_{}_{}", chrom, pos, ref, alt);

        stream << rs << "\t" << ref << "/" << alt << "\t" << chrom << "\t"
               << pos << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t";

        if (bcf_get_genotypes(header.get(), in_rec.get(), &gt.p_, &gt.n_) <= 0)
        {
            continue;
        }

        for (int i = 0; i < gt.n_ / 2; ++i)
        {
            int32_t gt0 = gt.p_[i * 2];
            int32_t gt1 = gt.p_[i * 2 + 1];
            if (bcf_gt_is_missing(gt0) || bcf_gt_is_missing(gt1))
            {
                stream << "NN\t";
            }
            else if (gt0 == gt1)
            {
                if (bcf_gt_allele(gt0) == 0)
                {
                    stream << ref << ref << "\t";
                }
                else
                {
                    stream << alt << alt << "\t";
                }
            }
            else
            {
                stream << ref << alt << "\t";
            }
        }
        stream << "\n";
    }
    counter->done();
}

}  // namespace vcfbox
