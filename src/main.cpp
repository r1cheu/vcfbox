#include "CLI11.hpp"
#include "utils.h"
#include "vcf.h"

int main(int argc, char** argv)
{
    CLI::App app{
        "A simple cli tool to combine genotypes according to paired sample "
        "list from vcf."};
    argv = app.ensure_utf8(argv);

    std::string vcf;
    std::string paired_sample;
    std::string output;
    bool keep_old_samples = false;
    app.add_option("-v,--vcf", vcf, "Path to input VCF file")->required();
    app.add_option(
           "-p,--paired-sample",
           paired_sample,
           "Path to file with paired samples, one pair per line, separated by "
           "space.")
        ->required();
    app.add_option(
           "-o,--output",
           output,
           "Path to output VCF file, if not provided, will be the same as "
           "input "
           "VCF file.")
        ->default_str("output.vcf");
    app.add_flag(
        "-k,--keep-old-samples",
        keep_old_samples,
        "Keep old samples in the output VCF file, default is false.");
    CLI11_PARSE(app, argc, argv);

    std::string mode = parse_mode(output);
    try
    {
        auto pairs = parse_sample_pairs(paired_sample);
        combine_genotypes(vcf, pairs, keep_old_samples, output, mode);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
