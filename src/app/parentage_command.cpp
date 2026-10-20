#include "app/subcommands.h"

#include <CLI/CLI.hpp>

#include "parentage/parentage.h"

namespace vcfbox
{
auto add_parentage_command(CLI::App& app, ParentageOptions& options)
    -> CLI::App*
{
    auto* command = app.add_subcommand(
        "parentage",
        "Assign F1 samples to parent pairs from a homozygous parent panel");

    command
        ->add_option(
            "-v,--vcf",
            options.vcf_path,
            "Panel VCF of bi-allelic SNP genotypes")
        ->required();
    command->add_option(
        "--exclude-sites",
        options.exclude_sites_path,
        "File of `contig<TAB>pos` (1-based) sites to drop from the panel");
    command
        ->add_option(
            "--duplicate-max-discordance",
            options.duplicate_max_discordance,
            "Warn about parent pairs differing on at most this fraction of "
            "sites")
        ->capture_default_str();
    command->add_option(
        "-b,--bams",
        options.bams_path,
        "List of BAMs: `sample<TAB>path` or `path` per line");
    command->add_option(
        "-o,--output",
        options.output_prefix,
        "Output prefix of the per-status .<status>.tsv, .topk.tsv and "
        ".summary.txt "
        "files; required with --bams");
    command->add_option(
        "--expected",
        options.expected_path,
        "Expected pairs: `sample<TAB>mother<TAB>father` per line");
    command->add_option(
        "--counts-dir",
        options.counts_dir,
        "Directory of per-sample counts files; existing ones skip the pileup");
    command->add_option("--min-mapq", options.min_mapq, "Minimum read MAPQ")
        ->capture_default_str();
    command
        ->add_option("--min-baseq", options.min_baseq, "Minimum base quality")
        ->capture_default_str();

    command->add_flag(
        "--keep-overlaps",
        options.keep_overlaps,
        "Count both mates of a pair at a site they both cover");
    command
        ->add_option(
            "--error-rate",
            options.error_rate,
            "Per-read probability of observing the wrong allele")
        ->capture_default_str();
    command
        ->add_option("--top", options.top, "Parent pairs reported per sample")
        ->capture_default_str();
    command
        ->add_option(
            "--candidates",
            options.candidates,
            "Parents with the best single-parent lod among which pairs are "
            "formed; 0 pairs every parent")
        ->capture_default_str();
    command
        ->add_option("--threads", options.threads, "Samples processed at once")
        ->capture_default_str();
    command
        ->add_option(
            "--max-hom-mismatch",
            options.max_hom_mismatch,
            "fit FAIL above this best-pair hom_mismatch_rate")
        ->capture_default_str();
    command
        ->add_option(
            "--min-lod",
            options.min_lod,
            "AMBIGUOUS when a model or parent margin is below this LOD")
        ->capture_default_str();

    return command;
}

auto run_parentage_command(const ParentageOptions& options) -> void
{
    parentage::RunOptions run_options;
    run_options.vcf_path = options.vcf_path;
    run_options.exclude_sites_path = options.exclude_sites_path;
    run_options.duplicate_max_discordance = options.duplicate_max_discordance;
    run_options.bams_path = options.bams_path;
    run_options.output_prefix = options.output_prefix;
    run_options.expected_path = options.expected_path;
    run_options.counts_dir = options.counts_dir;
    run_options.min_mapq = options.min_mapq;
    run_options.min_baseq = options.min_baseq;
    run_options.keep_overlaps = options.keep_overlaps;
    run_options.error_rate = options.error_rate;
    run_options.top = options.top;
    run_options.candidates = options.candidates;
    run_options.threads = options.threads;
    run_options.call.max_hom_mismatch = options.max_hom_mismatch;
    run_options.call.min_lod = options.min_lod;
    parentage::run(run_options);
}
}  // namespace vcfbox
