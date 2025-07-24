#include "CLI11.hpp"
#include "utils.h"
#include "vcf.h"

int main(int argc, char** argv)
{
    CLI::App app{"A simple cli tool contains some operation on vcf"};
    argv = app.ensure_utf8(argv);
    app.require_subcommand(1);
    std::string vcf;
    std::string paired_sample;
    std::string output;
    bool keep_old_samples = false;

    auto* combine = app.add_subcommand(
        "combine", "Combine genotypes from paired samples in a VCF file");
    auto* count
        = app.add_subcommand("count", "Count the number of snps in a VCF file");

    count->add_option("-v,--vcf", vcf, "Path to input VCF file")->required();
    combine->add_option("-v,--vcf", vcf, "Path to input VCF file")->required();
    combine
        ->add_option(
            "-p,--paired-sample",
            paired_sample,
            "Path to file with paired samples, one pair per line, separated by "
            "space.")
        ->required();
    combine
        ->add_option(
            "-o,--output",
            output,
            "Path to output VCF file, if not provided, will be the same as "
            "input "
            "VCF file.")
        ->default_str("output.vcf");
    combine->add_flag(
        "-k,--keep-old-samples",
        keep_old_samples,
        "Keep old samples in the output VCF file, default is false.");
    CLI11_PARSE(app, argc, argv);

    if (*combine)
    {
        std::string mode = vcfbox::parse_mode(output);
        try
        {
            auto pairs = vcfbox::parse_sample_pairs(paired_sample);
            vcfbox::combine_genotypes(
                vcf, pairs, keep_old_samples, output, mode);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        }
    }
    else if (*count)
    {
        vcfbox::count_records(vcf);
    }

    return 0;
}
