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
            "-v,--vcf",
            matrix_options.vcf_path,
            "Path to parent genotype VCF/BCF file")
        ->required();
    build_matrix
        ->add_option(
            "-m,--maternal",
            matrix_options.maternal_patterns,
            "Regex patterns matching maternal samples")
        ->expected(1, -1)
        ->required();
    build_matrix
        ->add_option(
            "-p,--paternal",
            matrix_options.paternal_patterns,
            "Regex patterns matching paternal samples")
        ->expected(1, -1)
        ->required();
    build_matrix
        ->add_option(
            "-o,--prefix",
            matrix_options.prefix,
            "Output prefix for parentage matrices")
        ->required();

    auto* test = parentage->add_subcommand(
        "test", "Test parentage from BAM list and parentage matrices");
    test->add_option("-b,--bam", test_options.bam_list_path, "Path to BAM list")
        ->required();
    test
        ->add_option(
            "-x,--prefix",
            test_options.prefix,
            "Input prefix for parentage matrices")
        ->required();
    test
        ->add_option(
            "-o,--output",
            test_options.output_path,
            "Path to top-K result TSV")
        ->required();
    test
        ->add_option(
            "-s,--summary",
            test_options.summary_path,
            "Path to per-sample summary TSV")
        ->required();
    test->add_option(
        "--raw",
        test_options.raw_path,
        "Optional path to full K_m x K_p loglik TSV");
    test->add_option(
        "-k,--top-k",
        test_options.top_k,
        "Number of top candidate pairs per sample");
    test->add_option(
        "--threshold",
        test_options.threshold,
        "Posterior threshold for accept/partial/reject call");
    test->add_option(
        "-e,--error-rate",
        test_options.error_rate,
        "Per-base sequencing error rate");
    test->add_option(
        "-q,--min-mapq", test_options.min_mapq, "Minimum read mapping quality");
    test->add_option(
        "-Q,--min-baseq", test_options.min_baseq, "Minimum base quality");
    test->add_option(
        "-t,--threads", test_options.threads, "Number of worker threads");

    return {parentage, build_matrix, test};
}
}  // namespace vcfbox
