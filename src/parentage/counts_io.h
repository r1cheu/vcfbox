#pragma once

#include <cstddef>
#include <string>

#include "parentage/counts.h"

namespace vcfbox::parentage
{
// Counts file (bgzf): magic `VBXCNT01`, sample name, site count, stats, then
// the n0 / n1 / n_other arrays. The site count is the only tie to the panel.

// Writes `<path>.partial` and renames it to `path`, so an interrupted write
// never leaves a file at `path`.
auto write_counts(const SampleCounts& counts, const std::string& path) -> void;

// Throws when the file is not a counts file, is truncated, or holds a
// different number of sites than `expected_n_sites`.
[[nodiscard]] auto read_counts(const std::string& path, size_t expected_n_sites)
    -> SampleCounts;
}  // namespace vcfbox::parentage
