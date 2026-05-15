#include "common/lines.h"

#include <fstream>
#include <stdexcept>
#include <utility>

namespace vcfbox
{
std::vector<std::string> read_lines(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
    {
        throw std::runtime_error("Cannot open file: " + path);
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line))
    {
        lines.push_back(std::move(line));
    }
    return lines;
}
}  // namespace vcfbox
