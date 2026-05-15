#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace vcfbox
{
class LineWriter
{
public:
    explicit LineWriter(std::string path);
    void write_line(std::string_view line);

    template <typename First, typename... Rest>
    void write_fields(const First& first, const Rest&... rest)
    {
        out_ << first;
        ((out_ << '\t' << rest), ...);
        out_ << '\n';
        if (!out_)
        {
            throw std::runtime_error("Failed writing: " + path_);
        }
    }

private:
    std::string path_;
    std::ofstream out_;
};
}  // namespace vcfbox
