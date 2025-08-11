#include <memory>
#include <string>
#include <string_view>
#include "barkeep.h"

extern "C"
{
#include <htslib/vcf.h>
}

namespace vcfbox
{
std::string parse_mode(std::string_view file_path);
size_t count_records(std::string_view vcf_path);

}  // namespace vcfbox

namespace detail
{
std::shared_ptr<barkeep::CompositeDisplay> create_progress(
    size_t total,
    size_t& progress_counters);

std::shared_ptr<barkeep::CompositeDisplay> create_counter(
    const std::string& message,
    size_t& progress_counters);

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
