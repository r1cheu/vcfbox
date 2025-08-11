#pragma once
#include <cstdlib>
#include <string>
#include "utils.h"

namespace vcfbox
{
std::vector<detail::SamplePair> parse_sample_pairs(
    const std::string& file_path);

void combine_genotypes(
    const std::string& vcf_path,
    const std::vector<detail::SamplePair>& sample_pairs,
    bool keep_old_samples,
    const std::string& out_path,
    const std::string& mode = "w");

void to_hapmap(const std::string& vcf_path, const std::string& out_path);

}  // namespace vcfbox
