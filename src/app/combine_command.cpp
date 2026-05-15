#include "app/subcommands.h"

#include <CLI/CLI.hpp>

#include "combine/combine_genotypes.h"

namespace vcfbox
{
CLI::App* add_combine_command(CLI::App& app, CombineOptions& options)
{
    auto* command = app.add_subcommand(
        "combine", "Combine genotypes from paired samples in a VCF file");

    command
        ->add_option("-v,--vcf", options.vcf_path, "Path to input VCF file")
        ->required();
    command
        ->add_option(
            "-p,--paired-sample",
            options.paired_sample_path,
            "Path to file with paired samples, one pair per line, separated "
            "by space.")
        ->required();
    command
        ->add_option("-o,--output", options.output_path, "Path to output VCF file")
        ->capture_default_str();
    command->add_flag(
        "-k,--keep-old-samples",
        options.keep_old_samples,
        "Keep old samples in the output VCF file.");

    return command;
}

void run_combine_command(const CombineOptions& options)
{
    auto pairs = parse_sample_pairs(options.paired_sample_path);
    combine_genotypes(
        options.vcf_path,
        pairs,
        options.keep_old_samples,
        options.output_path);
}
}  // namespace vcfbox
