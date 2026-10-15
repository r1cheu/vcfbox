#include "file_fixture.h"

#include <atomic>
#include <format>
#include <fstream>
#include <stdexcept>
#include <system_error>

#include <unistd.h>

namespace fs = std::filesystem;

namespace vcfbox::test
{
FileFixture::FileFixture()
{
    static std::atomic<int> instance{0};
    dir_ = fs::temp_directory_path()
           / std::format("vcfbox_test_{}_{}", ::getpid(), instance++);
    fs::create_directories(dir_);
}

FileFixture::~FileFixture()
{
    std::error_code ec;
    fs::remove_all(dir_, ec);
}

auto FileFixture::create_text_file(
    std::string_view content,
    std::string_view suffix) -> fs::path
{
    const fs::path path = dir_ / std::format("file_{}{}", counter_++, suffix);
    std::ofstream out(path, std::ios::binary);
    if (!out)
    {
        throw std::runtime_error(
            std::format("Failed to create {}", path.string()));
    }
    out << content;
    return path;
}

auto FileFixture::missing_path(std::string_view suffix) const -> fs::path
{
    return dir_ / std::format("missing{}", suffix);
}
}  // namespace vcfbox::test
