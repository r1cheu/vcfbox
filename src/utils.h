#include <memory>
#include <string>
#include <string_view>
#include "barkeep.h"

std::string parse_mode(std::string_view file_path);
namespace detail
{
std::shared_ptr<barkeep::CompositeDisplay> create_progress(
    size_t total,
    size_t& progress_counters);
}  // namespace detail
