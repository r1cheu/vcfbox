#pragma once

#include <cstddef>
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
    std::string output_prefix;
    std::string expected_path;
    std::string counts_dir;
    int min_mapq = 20;
    int min_baseq = 20;
    bool keep_overlaps = false;
    double error_rate = 0.01;
    double site_error_rate = 0.02;
    size_t top = 10;
    size_t candidates = 20;
    size_t threads = 1;
    double max_unmapped_rate = 0.8;
    double max_hom_mismatch = 0.03;
    double min_lod = 50;
};

CLI::App* add_parentage_command(CLI::App& app, ParentageOptions& options);
void run_parentage_command(const ParentageOptions& options);
}  // namespace vcfbox
