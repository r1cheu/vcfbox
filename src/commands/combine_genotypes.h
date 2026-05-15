#pragma once
#include <string>
#include <utility>
#include <vector>

namespace vcfbox
{
using SamplePair = std::pair<std::string, std::string>;

std::vector<SamplePair> parse_sample_pairs(const std::string& file_path);

void combine_genotypes(
    const std::string& vcf_path,
    const std::vector<SamplePair>& sample_pairs,
    bool keep_old_samples,
    const std::string& out_path,
    const std::string& mode = "w");
}  // namespace vcfbox
