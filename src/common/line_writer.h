#pragma once

#include <fstream>
#include <string>
#include <string_view>

namespace vcfbox
{
class LineWriter
{
public:
    explicit LineWriter(std::string path);
    void write_line(std::string_view line);

private:
    std::string path_;
    std::ofstream out_;
};
}  // namespace vcfbox
