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
}  // namespace vcfbox
