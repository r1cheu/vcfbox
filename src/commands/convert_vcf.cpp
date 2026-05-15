#include "commands/convert_vcf.h"

#include <cstdint>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>

#include "common/progress.h"
#include "hts/hts_raii.h"

extern "C"
{
#include <htslib/vcf.h>
}

namespace vcfbox
{
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
    auto counter = create_counter("Converting to HapMap format", processd_snp);

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
