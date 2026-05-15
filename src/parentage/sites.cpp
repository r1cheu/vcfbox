#include "parentage/sites.h"

#include <charconv>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "common/strings.h"

namespace vcfbox
{
namespace
{
Site parse_line(std::string_view line, const std::string& path, size_t lineno)
{
    const auto fields = split(line, '\t');
    if (fields.size() < 5 || fields[0].empty() || fields[1].empty()
        || fields[3].size() != 1 || fields[4].size() != 1)
    {
        throw std::runtime_error(
            "Malformed or non-SNP BED line at " + path + ":"
            + std::to_string(lineno));
    }

    int64_t pos = 0;
    const auto& pos_field = fields[1];
    const auto* first = pos_field.data();
    const auto* last = first + pos_field.size();
    const auto [ptr, ec] = std::from_chars(first, last, pos);
    if (ec != std::errc{} || ptr != last)
    {
        throw std::runtime_error(
            "Bad pos in BED line at " + path + ":" + std::to_string(lineno));
    }

    return Site{
        std::string(fields[0]), pos, fields[3].front(), fields[4].front()};
}
}  // namespace

std::vector<Site> load_sites(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
    {
        throw std::runtime_error("Cannot open sites bed: " + path);
    }

    std::vector<Site> sites;
    std::string line;
    size_t lineno = 0;
    while (std::getline(in, line))
    {
        ++lineno;
        if (line.empty())
        {
            continue;
        }
        sites.push_back(parse_line(line, path, lineno));
    }
    return sites;
}
}  // namespace vcfbox
