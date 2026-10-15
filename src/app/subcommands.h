#pragma once

#include <string>

namespace CLI
{
class App;
}

namespace vcfbox
{
struct CombineOptions
{
    std::string vcf_path;
    std::string paired_sample_path;
    std::string output_path = "output.vcf";
    bool keep_old_samples = false;
};

CLI::App* add_combine_command(CLI::App& app, CombineOptions& options);
void run_combine_command(const CombineOptions& options);

struct ParentageOptions
{
    std::string vcf_path;
    std::string exclude_sites_path;
    double duplicate_max_discordance = 0.001;
    std::string bams_path;
    std::string counts_dir;
    int min_mapq = 20;
    int min_baseq = 20;
    bool keep_overlaps = false;
};

CLI::App* add_parentage_command(CLI::App& app, ParentageOptions& options);
void run_parentage_command(const ParentageOptions& options);
}  // namespace vcfbox
