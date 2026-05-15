#pragma once

#include <string>

#include "parentage/parentage.h"

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

struct ParentageCommands
{
    CLI::App* parentage;
    CLI::App* build_matrix;
    CLI::App* test;
};

CLI::App* add_combine_command(CLI::App& app, CombineOptions& options);
void run_combine_command(const CombineOptions& options);

ParentageCommands add_parentage_command(
    CLI::App& app,
    ParentageMatrixOptions& matrix_options,
    ParentageTestOptions& test_options);
}  // namespace vcfbox
