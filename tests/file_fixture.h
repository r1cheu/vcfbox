#pragma once

#include <filesystem>
#include <string_view>

namespace vcfbox::test
{
// A temporary directory that lives as long as the fixture; every file
// created through it is removed with the directory.
class FileFixture
{
   public:
    FileFixture();
    ~FileFixture();

    FileFixture(const FileFixture&) = delete;
    auto operator=(const FileFixture&) -> FileFixture& = delete;
    FileFixture(FileFixture&&) = delete;
    auto operator=(FileFixture&&) -> FileFixture& = delete;

    // Writes `content` to a fresh file ending in `suffix` and returns its
    // path.
    [[nodiscard]] auto create_text_file(
        std::string_view content,
        std::string_view suffix = ".txt") -> std::filesystem::path;

    // A path inside the directory that no file has been created at.
    [[nodiscard]] auto missing_path(std::string_view suffix = ".txt") const
        -> std::filesystem::path;

    auto dir() const -> const std::filesystem::path& { return dir_; }

   private:
    std::filesystem::path dir_;
    int counter_ = 0;
};
}  // namespace vcfbox::test
