#include "common/path.h"

#include <string>

namespace vcfbox
{
std::string parse_mode(std::string_view file_path)
{
    auto ext = file_path.substr(file_path.find_last_of('.') + 1);
    if (ext == "bam")
    {
        return "wb";
    }
    if (ext == "sam")
    {
        return "w";
    }
    if (ext == "cram")
    {
        return "wc";
    }
    if (ext == "bcf")
    {
        return "wb";
    }
    if (ext == "vcf")
    {
        return "w";
    }
    if (ext == "gz")
    {
        auto base = file_path.substr(0, file_path.find_last_of('.'));
        auto base_ext = base.substr(base.find_last_of('.') + 1);
        if (base_ext == "vcf")
        {
            return "wz";
        }
        if (base_ext == "bcf")
        {
            return "wb";
        }
    }
    return "w";
}
}  // namespace vcfbox
