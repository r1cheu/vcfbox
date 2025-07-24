#include <memory>
#include <string>
#include <string_view>
#include "barkeep.h"

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
}  // namespace detail
