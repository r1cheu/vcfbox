#pragma once
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C"
{
#include <htslib/vcf.h>
}

namespace detail
{
using SamplePair = std::pair<std::string, std::string>;

void check_sample_consistence(
    std::string_view vcf_path,
    const std::vector<SamplePair>& sample_pairs);

bcf_hdr_t* init_bcf_head(
    bcf_hdr_t* header,
    const std::vector<SamplePair>& sample_pairs,
    bool keep_old_samples);

void copy_rec_info(
    bcf_hdr_t* header,
    bcf_hdr_t* output_header,
    bcf1_t* in_rec,
    bcf1_t* out_rec);

std::vector<int32_t> concat_gt(
    const std::vector<std::pair<std::string, std::string>>& sample_pairs,
    const std::unordered_map<std::string, int>& sample_to_idx,
    const int32_t* gt_arr,
    bool keep_old_samples,
    int n_gt);
}  // namespace detail
//

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

}  // namespace vcfbox
