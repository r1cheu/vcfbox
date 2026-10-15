#include "parentage/parentage.h"

#include <format>
#include <iostream>

#include "parentage/panel.h"

namespace vcfbox::parentage
{
namespace
{
auto print_panel_summary(const Panel& panel) -> void
{
    std::cerr << std::format(
        "Panel: {} parents, {} sites on {} contigs{}\n",
        panel.n_parents(),
        panel.n_sites(),
        panel.sites().contigs().size(),
        panel.n_excluded_sites() > 0
            ? std::format(" ({} excluded)", panel.n_excluded_sites())
            : "");

    for (const auto& dup : panel.near_duplicates())
    {
        std::cerr << std::format(
            "Warning: parents {} and {} differ at only {} of {} sites; "
            "samples derived from either will be AMBIGUOUS\n",
            panel.names()[dup.parent_a],
            panel.names()[dup.parent_b],
            dup.n_discordant,
            panel.n_sites());
    }
}
}  // namespace

auto run(const RunOptions& options) -> void
{
    PanelOptions panel_options;
    panel_options.exclude_sites_path = options.exclude_sites_path;
    panel_options.duplicate_max_discordance = options.duplicate_max_discordance;

    const Panel panel = load_panel(options.vcf_path, panel_options);
    print_panel_summary(panel);
}
}  // namespace vcfbox::parentage
