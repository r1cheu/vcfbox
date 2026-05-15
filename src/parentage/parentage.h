#pragma once

#include <string>
#include <vector>

namespace vcfbox
{
struct ParentageMatrixOptions
{
    std::string vcf_path;
    std::vector<std::string> maternal_patterns;
    std::vector<std::string> paternal_patterns;
    std::string prefix;
};

struct ParentageTestOptions
{
    std::string bam_list_path;
    std::string prefix;
    std::string output_path;
    double error_rate = 0.01;
    int min_mapq = 20;
    int min_baseq = 13;
    int threads = 1;
};

void build_parentage_matrices(const ParentageMatrixOptions& options);
void test_parentage(const ParentageTestOptions& options);
}  // namespace vcfbox
