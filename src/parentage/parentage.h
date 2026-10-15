#pragma once

#include <string>
#include <vector>

namespace vcfbox::parentage
{
struct RunOptions
{
    std::string vcf_path;
    std::string exclude_sites_path;
    double duplicate_max_discordance = 0.001;
    std::string bams_path;   // empty: load the panel only
    std::string counts_dir;  // empty: never read or write counts files
    int min_mapq = 20;
    int min_baseq = 20;
    bool keep_overlaps = false;  // count both mates over the same site
};

struct BamEntry
{
    std::string sample;
    std::string path;
};

// Parses a BAM list: one `sample<TAB>path` or bare `path` per line (sample =
// file name without `.bam`); blank lines and `#` comments are skipped.
// Throws on a line with more than two fields or a repeated sample name.
[[nodiscard]] auto read_bam_list(const std::string& path)
    -> std::vector<BamEntry>;

// Loads the panel, then counts alleles for every listed BAM (or reads the
// counts back from `counts_dir`). Scoring stages are added on top of this
// entry point as they are implemented.
auto run(const RunOptions& options) -> void;
}  // namespace vcfbox::parentage
