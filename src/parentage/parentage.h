#pragma once

#include <string>
#include <vector>

namespace vcfbox
{
struct ParentageMatrixOptions
{
    std::string parents_path;
    std::vector<std::string> maternal_patterns;
    std::vector<std::string> paternal_patterns;
    std::string matrix_prefix;
};

struct ParentageTestOptions
{
    std::string bam_list_path;
    std::string matrix_prefix;
    std::string output_path;
    double error_rate = 0.01;
    int min_mapq = 20;
    int min_baseq = 13;
};

void build_parentage_matrices(const ParentageMatrixOptions& options);
void test_parentage(const ParentageTestOptions& options);
}  // namespace vcfbox
