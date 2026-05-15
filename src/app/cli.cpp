#include "cli.h"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include <CLI/CLI.hpp>
#include "commands/combine_genotypes.h"
#include "commands/convert_vcf.h"
#include "common/path.h"
#include "parentage/parentage.h"

namespace
{
struct CombineOptions
{
    std::string vcf_path;
    std::string paired_sample_path;
    std::string output_path = "output.vcf";
    bool keep_old_samples = false;
};

struct ConvertOptions
{
    std::string vcf_path;
    std::string output_path = "output.hmp";
};

template <typename Function>
int run_command(Function run)
{
    try
    {
        run();
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}

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
        ->add_option(
            "-o,--output",
            options.output_path,
            "Path to output VCF file, if not provided, will be the same as "
            "input VCF file.")
        ->default_str(options.output_path);
    command->add_flag(
        "-k,--keep-old-samples",
        options.keep_old_samples,
        "Keep old samples in the output VCF file, default is false.");

    return command;
}

CLI::App* add_convert_command(CLI::App& app, ConvertOptions& options)
{
    auto* command = app.add_subcommand(
        "convert", "Convert VCF to various formats (e.g., HapMap)");
    command
        ->add_option("-v,--vcf", options.vcf_path, "Path to input VCF file")
        ->required();
    command
        ->add_option(
            "-o,--output",
            options.output_path,
            "Path to output file, if not provided, will be the same as input "
            "VCF file.")
        ->default_str(options.output_path);
    return command;
}

CLI::App* add_parentage_command(
    CLI::App& app,
    vcfbox::ParentageMatrixOptions& matrix_options,
    vcfbox::ParentageTestOptions& test_options,
    CLI::App*& build_matrix_command,
    CLI::App*& test_command)
{
    auto* command = app.add_subcommand(
        "parentage", "Build parentage matrices and test parentage from BAMs");
    command->require_subcommand(1);

    build_matrix_command = command->add_subcommand(
        "build-matrix", "Build parentage indicator matrices from parent VCF");
    build_matrix_command
        ->add_option(
            "--parents",
            matrix_options.parents_path,
            "Path to parent genotype VCF/BCF file")
        ->required();
    build_matrix_command
        ->add_option(
            "--maternal",
            matrix_options.maternal_patterns,
            "Regex patterns matching maternal samples")
        ->expected(1, -1)
        ->required();
    build_matrix_command
        ->add_option(
            "--paternal",
            matrix_options.paternal_patterns,
            "Regex patterns matching paternal samples")
        ->expected(1, -1)
        ->required();
    build_matrix_command
        ->add_option(
            "--matrix-prefix",
            matrix_options.matrix_prefix,
            "Output prefix for parentage matrices")
        ->required();

    test_command = command->add_subcommand(
        "test", "Test parentage from BAM list and parentage matrices");
    test_command
        ->add_option("--bam", test_options.bam_list_path, "Path to BAM list")
        ->required();
    test_command
        ->add_option(
            "--matrix-prefix",
            test_options.matrix_prefix,
            "Input prefix for parentage matrices")
        ->required();
    test_command
        ->add_option(
            "-o,--output", test_options.output_path, "Path to output result")
        ->required();

    return command;
}

void run_combine(const CombineOptions& options)
{
    std::string mode = vcfbox::parse_mode(options.output_path);
    auto pairs = vcfbox::parse_sample_pairs(options.paired_sample_path);
    vcfbox::combine_genotypes(
        options.vcf_path,
        pairs,
        options.keep_old_samples,
        options.output_path,
        mode);
}

void run_convert(const ConvertOptions& options)
{
    if (options.output_path.substr(options.output_path.find_last_of('.') + 1)
        == "hmp")
    {
        vcfbox::to_hapmap(options.vcf_path, options.output_path);
        return;
    }

    throw std::runtime_error("Unsupported format: " + options.output_path);
}
}  // namespace

namespace vcfbox
{
int run_cli(int argc, char** argv)
{
    CLI::App app{"A simple cli tool contains some operation on vcf"};
    argv = app.ensure_utf8(argv);
    app.require_subcommand(1);

    CombineOptions combine_options;
    ConvertOptions convert_options;
    ParentageMatrixOptions parentage_matrix_options;
    ParentageTestOptions parentage_test_options;

    auto* combine_command = add_combine_command(app, combine_options);
    auto* convert_command = add_convert_command(app, convert_options);
    CLI::App* build_matrix_command = nullptr;
    CLI::App* parentage_test_command = nullptr;
    add_parentage_command(
        app,
        parentage_matrix_options,
        parentage_test_options,
        build_matrix_command,
        parentage_test_command);

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

    if (*combine_command)
    {
        return run_command([&] { run_combine(combine_options); });
    }
    if (*convert_command)
    {
        return run_command([&] { run_convert(convert_options); });
    }
    if (*build_matrix_command)
    {
        return run_command(
            [&] { build_parentage_matrices(parentage_matrix_options); });
    }
    if (*parentage_test_command)
    {
        return run_command([&] { test_parentage(parentage_test_options); });
    }

    return 0;
}
}  // namespace vcfbox
