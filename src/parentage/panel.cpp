#include "parentage/panel.h"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <format>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/lines.h"
#include "common/progress.h"
#include "common/strings.h"
#include "hts/hts_raii.h"

extern "C"
{
#include <htslib/vcf.h>
}

namespace vcfbox::parentage
{
namespace
{
using ExcludedSites
    = std::unordered_map<std::string, std::unordered_set<int32_t>>;

auto read_excluded_sites(const std::string& path) -> ExcludedSites
{
    ExcludedSites excluded;
    if (path.empty())
    {
        return excluded;
    }
    size_t line_no = 0;
    for (const auto& line : read_lines(path))
    {
        ++line_no;
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        const auto fields = split(line, '\t');
        int32_t pos = 0;
        if (fields.size() < 2
            || std::from_chars(
                   fields[1].data(), fields[1].data() + fields[1].size(), pos)
                       .ec
                   != std::errc{})
        {
            throw std::runtime_error(
                std::format("{}:{}: expected `contig<TAB>pos`", path, line_no));
        }
        excluded[std::string(fields[0])].insert(pos - 1);
    }
    return excluded;
}

auto gt_to_string(int32_t a0, int32_t a1) -> std::string
{
    auto one = [](int32_t a) -> std::string
    {
        if (a == bcf_int32_vector_end)
        {
            return "";
        }
        return bcf_gt_is_missing(a) ? "." : std::format("{}", bcf_gt_allele(a));
    };
    const auto s1 = one(a1);
    return s1.empty() ? one(a0) : std::format("{}/{}", one(a0), s1);
}

auto site_label(const bcf_hdr_t* hdr, const bcf1_t* rec) -> std::string
{
    return std::format("{}:{}", bcf_hdr_id2name(hdr, rec->rid), rec->pos + 1);
}

// Per-parent staging: 64 sites are accumulated into one word per parent,
// then appended to that parent's row chunk list.
class RowBuilder
{
   public:
    explicit RowBuilder(size_t n_parents)
        : chunks_(n_parents), current_(n_parents, 0)
    {
    }

    auto set(size_t parent) -> void { current_[parent] |= uint64_t{1} << bit_; }

    auto next_site() -> void
    {
        if (++bit_ == 64)
        {
            flush();
        }
    }

    auto finish(size_t n_sites) -> BitRows
    {
        if (bit_ != 0)
        {
            flush();
        }
        BitRows rows(chunks_.size(), n_sites);
        for (size_t k = 0; k < chunks_.size(); ++k)
        {
            std::ranges::copy(chunks_[k], rows.row(k).begin());
        }
        return rows;
    }

   private:
    auto flush() -> void
    {
        for (size_t k = 0; k < chunks_.size(); ++k)
        {
            chunks_[k].push_back(current_[k]);
            current_[k] = 0;
        }
        bit_ = 0;
    }

    std::vector<std::vector<uint64_t>> chunks_;
    std::vector<uint64_t> current_;
    unsigned bit_ = 0;
};

auto check_alleles(const bcf_hdr_t* hdr, const bcf1_t* rec) -> void
{
    if (rec->n_allele != 2)
    {
        throw std::runtime_error(
            std::format(
                "{}: panel sites must be bi-allelic, found {} allele(s)",
                site_label(hdr, rec),
                rec->n_allele));
    }
    for (int i = 0; i < 2; ++i)
    {
        const char* a = rec->d.allele[i];
        if (std::strlen(a) != 1 || encode_base(a[0]) == no_base)
        {
            throw std::runtime_error(
                std::format(
                    "{}: panel sites must be SNPs with A/C/G/T alleles, found "
                    "{}",
                    site_label(hdr, rec),
                    a));
        }
    }
}

struct OpenedVcf
{
    hts_file_ptr file;
    bcf_hdr_ptr hdr;
};

auto open_panel_vcf(const std::string& vcf_path) -> OpenedVcf
{
    hts_file_ptr file(hts_open(vcf_path.c_str(), "r"));
    if (!file)
    {
        throw std::runtime_error(
            std::format("Failed to open panel VCF: {}", vcf_path));
    }
    bcf_hdr_ptr hdr(bcf_hdr_read(file.get()));
    if (!hdr)
    {
        throw std::runtime_error(
            std::format("Failed to read VCF header: {}", vcf_path));
    }
    if (bcf_hdr_nsamples(hdr.get()) == 0)
    {
        throw std::runtime_error(
            std::format("Panel VCF has no samples: {}", vcf_path));
    }
    return {std::move(file), std::move(hdr)};
}

auto sample_names(const bcf_hdr_t* hdr) -> std::vector<std::string>
{
    const std::span<char*> samples(
        hdr->samples, static_cast<size_t>(bcf_hdr_nsamples(hdr)));
    return {samples.begin(), samples.end()};
}

// Contig name and length from the `##contig` header line of `rid`.
auto contig_from_header(const bcf_hdr_t* hdr, int rid) -> ContigInfo
{
    ContigInfo info;
    info.name = bcf_hdr_id2name(hdr, rid);
    info.length = static_cast<int64_t>(hdr->id[BCF_DT_CTG][rid].val->info[0]);
    if (info.length <= 0)
    {
        throw std::runtime_error(
            std::format(
                "Panel VCF header has no length for contig {}; a "
                "`##contig=<ID=...,length=...>` line is required",
                info.name));
    }
    return info;
}

// Allele1 dosage (0, 1 or 2) of a diploid GT; nullopt when either allele is
// missing. Throws for haploid calls.
auto allele1_dosage(
    int32_t a0,
    int32_t a1,
    const bcf_hdr_t* hdr,
    const bcf1_t* rec,
    const std::string& parent) -> std::optional<int>
{
    if (a1 == bcf_int32_vector_end)
    {
        throw std::runtime_error(
            std::format(
                "{} {}: haploid genotype {}; parents must be diploid "
                "homozygous",
                site_label(hdr, rec),
                parent,
                gt_to_string(a0, a1)));
    }
    if (bcf_gt_is_missing(a0) || bcf_gt_is_missing(a1))
    {
        return std::nullopt;
    }
    return bcf_gt_allele(a0) + bcf_gt_allele(a1);
}

// Bit planes being built for every parent.
struct PlaneBuilders
{
    RowBuilder geno;
    RowBuilder het;
    RowBuilder missing;
};

// Reads the GT of every parent at `rec`, sets its bits in `planes` and
// advances them to the next site. Returns the allele1 frequency over the
// called parents; throws when no parent is called.
auto read_site_genotypes(
    const bcf_hdr_t* hdr,
    bcf1_t* rec,
    genotype_buffer& gt,
    std::span<const std::string> names,
    PlaneBuilders& planes) -> float
{
    const int n_gt = read_genotypes(hdr, rec, gt);
    if (n_gt != 2 * static_cast<int>(names.size()))
    {
        throw std::runtime_error(
            std::format(
                "{}: expected diploid GT for every parent",
                site_label(hdr, rec)));
    }
    const std::span<const int32_t> g(gt.data.get(), static_cast<size_t>(n_gt));
    size_t dosage = 0;
    size_t n_called = 0;
    for (size_t k = 0; k < names.size(); ++k)
    {
        const std::optional<int> d
            = allele1_dosage(g[2 * k], g[(2 * k) + 1], hdr, rec, names[k]);
        if (!d)
        {
            planes.missing.set(k);
            continue;
        }
        if (*d == 2)
        {
            planes.geno.set(k);
        }
        else if (*d == 1)
        {
            planes.het.set(k);
        }
        dosage += static_cast<size_t>(*d);
        ++n_called;
    }
    if (n_called == 0)
    {
        throw std::runtime_error(
            std::format("{}: no parent is called", site_label(hdr, rec)));
    }
    planes.geno.next_site();
    planes.het.next_site();
    planes.missing.next_site();
    return static_cast<float>(dosage) / static_cast<float>(2 * n_called);
}
}  // namespace

auto load_panel(const std::string& vcf_path, const PanelOptions& options)
    -> Panel
{
    const auto excluded = read_excluded_sites(options.exclude_sites_path);
    auto [file, hdr] = open_panel_vcf(vcf_path);
    std::vector<std::string> names = sample_names(hdr.get());

    PlaneBuilders planes{
        RowBuilder(names.size()),
        RowBuilder(names.size()),
        RowBuilder(names.size())};
    SiteTable sites;
    size_t n_excluded_sites = 0;
    int current_rid = -1;
    const std::unordered_set<int32_t>* current_excluded = nullptr;

    bcf1_ptr rec(bcf_init());
    genotype_buffer gt;
    size_t processed = 0;
    auto counter = create_counter("Loading panel", processed, "site/s");
    counter->show();

    while (bcf_read(file.get(), hdr.get(), rec.get()) == 0)
    {
        ++processed;
        if (rec->rid != current_rid)
        {
            current_rid = rec->rid;
            ContigInfo info = contig_from_header(hdr.get(), rec->rid);
            const auto it = excluded.find(info.name);
            current_excluded = it == excluded.end() ? nullptr : &it->second;
            sites.add_contig(std::move(info));
        }
        if (current_excluded != nullptr
            && current_excluded->contains(static_cast<int32_t>(rec->pos)))
        {
            ++n_excluded_sites;
            continue;
        }

        bcf_unpack(rec.get(), BCF_UN_STR);
        check_alleles(hdr.get(), rec.get());
        const float freq1
            = read_site_genotypes(hdr.get(), rec.get(), gt, names, planes);
        sites.add_site(
            rec->pos,
            encode_base(rec->d.allele[0][0]),
            encode_base(rec->d.allele[1][0]),
            freq1);
    }
    counter->done();
    sites.finish();

    BitRows geno = planes.geno.finish(sites.size());
    BitRows het = planes.het.finish(sites.size());
    BitRows missing = planes.missing.finish(sites.size());
    return Panel(
        std::move(sites),
        std::move(names),
        std::move(geno),
        std::move(het),
        std::move(missing),
        options.duplicate_max_discordance,
        n_excluded_sites);
}

Panel::Panel(
    SiteTable sites,
    std::vector<std::string> names,
    BitRows geno,
    BitRows het,
    BitRows missing,
    double duplicate_max_discordance,
    size_t n_excluded_sites)
    : sites_(std::move(sites)),
      names_(std::move(names)),
      geno_(std::move(geno)),
      het_(std::move(het)),
      missing_(std::move(missing)),
      n_excluded_sites_(n_excluded_sites)
{
    if (geno_.rows() != names_.size())
    {
        throw std::logic_error(
            std::format(
                "Panel: {} parent names but {} genotype rows",
                names_.size(),
                geno_.rows()));
    }
    if (geno_.n_bits() != sites_.size())
    {
        throw std::logic_error(
            std::format(
                "Panel: {} sites but {} genotype columns",
                sites_.size(),
                geno_.n_bits()));
    }
    if (het_.rows() != geno_.rows() || het_.n_bits() != geno_.n_bits())
    {
        throw std::logic_error(
            std::format(
                "Panel: genotype matrix is {} x {} but heterozygosity matrix "
                "is {} x {}",
                geno_.rows(),
                geno_.n_bits(),
                het_.rows(),
                het_.n_bits()));
    }
    if (missing_.rows() != geno_.rows() || missing_.n_bits() != geno_.n_bits())
    {
        throw std::logic_error(
            std::format(
                "Panel: genotype matrix is {} x {} but missing matrix is {} x "
                "{}",
                geno_.rows(),
                geno_.n_bits(),
                missing_.rows(),
                missing_.n_bits()));
    }
    for (size_t i = 0; i < names_.size(); ++i)
    {
        n_het_calls_ += popcount_rows(het_.row(i));
        n_missing_calls_ += popcount_rows(missing_.row(i));
    }

    const auto max_diff = static_cast<size_t>(
        duplicate_max_discordance * static_cast<double>(sites_.size()));
    for (size_t i = 0; i < names_.size(); ++i)
    {
        for (size_t j = i + 1; j < names_.size(); ++j)
        {
            const size_t n_diff
                = popcount_discordant(parent_rows(i), parent_rows(j));
            if (n_diff <= max_diff)
            {
                near_duplicates_.push_back(
                    {static_cast<uint32_t>(i),
                     static_cast<uint32_t>(j),
                     n_diff});
            }
        }
    }
}
}  // namespace vcfbox::parentage
