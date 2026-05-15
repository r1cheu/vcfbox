#include "app/subcommands.h"

#include <CLI/CLI.hpp>

#include "parentage/parentage.h"

namespace vcfbox
{
ParentageCommands add_parentage_command(
    CLI::App& app,
    ParentageMatrixOptions& matrix_options,
    ParentageTestOptions& test_options)
{
    auto* parentage = app.add_subcommand(
        "parentage", "Build parentage matrices and test parentage from BAMs");
    parentage->require_subcommand(1);

    auto* build_matrix = parentage->add_subcommand(
        "build-matrix", "Build parentage indicator matrices from parent VCF");
    build_matrix
        ->add_option(
            "--parents",
            matrix_options.parents_path,
            "Path to parent genotype VCF/BCF file")
        ->required();
    build_matrix
        ->add_option(
            "--maternal",
            matrix_options.maternal_patterns,
            "Regex patterns matching maternal samples")
        ->expected(1, -1)
        ->required();
    build_matrix
        ->add_option(
            "--paternal",
            matrix_options.paternal_patterns,
            "Regex patterns matching paternal samples")
        ->expected(1, -1)
        ->required();
    build_matrix
        ->add_option(
            "--matrix-prefix",
            matrix_options.matrix_prefix,
            "Output prefix for parentage matrices")
        ->required();

    auto* test = parentage->add_subcommand(
        "test", "Test parentage from BAM list and parentage matrices");
    test
        ->add_option("--bam", test_options.bam_list_path, "Path to BAM list")
        ->required();
    test
        ->add_option(
            "--matrix-prefix",
            test_options.matrix_prefix,
            "Input prefix for parentage matrices")
        ->required();
    test
        ->add_option(
            "-o,--output", test_options.output_path, "Path to output result")
        ->required();

    return {parentage, build_matrix, test};
}
}  // namespace vcfbox
