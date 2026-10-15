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
            "Panel VCF (homozygous, fully called)")
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

    return command;
}

auto run_parentage_command(const ParentageOptions& options) -> void
{
    parentage::RunOptions run_options;
    run_options.vcf_path = options.vcf_path;
    run_options.exclude_sites_path = options.exclude_sites_path;
    run_options.duplicate_max_discordance = options.duplicate_max_discordance;
    parentage::run(run_options);
}
}  // namespace vcfbox
