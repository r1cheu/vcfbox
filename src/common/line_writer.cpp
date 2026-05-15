#include "common/line_writer.h"

#include <stdexcept>
#include <utility>

namespace vcfbox
{
LineWriter::LineWriter(std::string path)
    : path_(std::move(path)), out_(path_)
{
    if (!out_)
    {
        throw std::runtime_error("Cannot open output file: " + path_);
    }
}

void LineWriter::write_line(std::string_view line)
{
    out_ << line << '\n';
    if (!out_)
    {
        throw std::runtime_error("Failed writing: " + path_);
    }
}
}  // namespace vcfbox
