#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <string>
#include <vector>

#include "file_fixture.h"
#include "parentage/counts.h"
#include "parentage/pileup.h"
#include "parentage/sites.h"
#include "sam_fixture.h"

using Catch::Matchers::ContainsSubstring;
using vcfbox::parentage::encode_base;
using vcfbox::parentage::pileup_bam;
using vcfbox::parentage::PileupOptions;
using vcfbox::parentage::SampleCounts;
using vcfbox::parentage::SiteTable;
using vcfbox::test::FileFixture;
using vcfbox::test::repeated_reads;
using vcfbox::test::SamFile;
using vcfbox::test::SamRead;

namespace
{
// Five sites: chr01 at 1-based 10 A>G, 20 C>T, 30 G>A, 35 T>C; chr02 at 5
// A>G. Indices 0..4.
auto make_sites() -> SiteTable
{
    SiteTable sites;
    sites.add_contig({"chr01", 1000});
    sites.add_site(9, encode_base('A'), encode_base('G'), 0.5F);
    sites.add_site(19, encode_base('C'), encode_base('T'), 0.5F);
    sites.add_site(29, encode_base('G'), encode_base('A'), 0.5F);
    sites.add_site(34, encode_base('T'), encode_base('C'), 0.5F);
    sites.add_contig({"chr02", 500});
    sites.add_site(4, encode_base('A'), encode_base('G'), 0.5F);
    sites.finish();
    return sites;
}

auto run(FileFixture& files, const SamFile& sam, PileupOptions options = {})
    -> SampleCounts
{
    return pileup_bam(
        files.create_text_file(sam.text(), ".sam"), "S", make_sites(), options);
}

auto n0_of(const SampleCounts& counts) -> std::vector<uint8_t>
{
    return {counts.n0().begin(), counts.n0().end()};
}

auto n1_of(const SampleCounts& counts) -> std::vector<uint8_t>
{
    return {counts.n1().begin(), counts.n1().end()};
}

const std::string q30(30, 'I');  // base quality 40 at every position

// 30 bases from position 1: allele0 at sites 10, 20, 30; N elsewhere.
const SamRead allele0_read{
    .cigar = "30M",
    .seq = "NNNNNNNNNANNNNNNNNNCNNNNNNNNNG",
    .qual = q30};
}  // namespace

TEST_CASE(
    "pileup_bam counts the base under each panel site by allele",
    "[parentage][pileup]")
{
    FileFixture files;
    SamFile sam;

    SECTION("allele0 bases")
    {
        sam.reads = {allele0_read};
        const SampleCounts counts = run(files, sam);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{1, 1, 1, 0, 0});
        REQUIRE(n1_of(counts) == std::vector<uint8_t>{0, 0, 0, 0, 0});
        REQUIRE(counts.sites_covered() == 3);
        REQUIRE(counts.bases_used() == 3);
        REQUIRE(counts.stats().reads_total == 1);
        REQUIRE(counts.stats().reads_used == 1);
    }
    SECTION("allele1 bases")
    {
        sam.reads
            = {{.cigar = "30M",
                .seq = "NNNNNNNNNGNNNNNNNNNTNNNNNNNNNA",
                .qual = q30}};
        const SampleCounts counts = run(files, sam);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{0, 0, 0, 0, 0});
        REQUIRE(n1_of(counts) == std::vector<uint8_t>{1, 1, 1, 0, 0});
    }
    SECTION("a third base goes to n_other and N is ignored")
    {
        sam.reads
            = {{.cigar = "30M",
                .seq = "NNNNNNNNNTNNNNNNNNNNNNNNNNNNNN",
                .qual = q30}};
        const SampleCounts counts = run(files, sam);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{0, 0, 0, 0, 0});
        REQUIRE(n1_of(counts) == std::vector<uint8_t>{0, 0, 0, 0, 0});
        REQUIRE(counts.n_other()[0] == 1);
        REQUIRE(counts.n_other()[1] == 0);
        REQUIRE(counts.sites_covered() == 0);
    }
    SECTION("a read beyond the last site counts nothing but is used")
    {
        SamRead read = allele0_read;
        read.pos = 900;
        sam.reads = {read};
        const SampleCounts counts = run(files, sam);
        REQUIRE(counts.bases_used() == 0);
        REQUIRE(counts.stats().reads_used == 1);
    }
    SECTION("the site cursor restarts on the next contig")
    {
        sam.reads
            = {allele0_read,
               {.contig = "chr02",
                .cigar = "10M",
                .seq = "NNNNANNNNN",
                .qual = std::string(10, 'I')}};
        const SampleCounts counts = run(files, sam);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{1, 1, 1, 0, 1});
    }
    SECTION("no reads at all")
    {
        const SampleCounts counts = run(files, sam);
        REQUIRE(counts.bases_used() == 0);
        REQUIRE(counts.stats().reads_total == 0);
    }
}

TEST_CASE("pileup_bam skips reads by flag and MAPQ", "[parentage][pileup]")
{
    FileFixture files;
    SamFile sam;
    SamRead read = allele0_read;

    SECTION("unmapped")
    {
        read.flag = 0x4;
        sam.reads = {read};
        const SampleCounts counts = run(files, sam);
        REQUIRE(counts.bases_used() == 0);
        REQUIRE(counts.stats().reads_skipped_flag == 1);
    }
    SECTION("secondary")
    {
        read.flag = 0x100;
        sam.reads = {read};
        REQUIRE(run(files, sam).stats().reads_skipped_flag == 1);
    }
    SECTION("supplementary")
    {
        read.flag = 0x800;
        sam.reads = {read};
        REQUIRE(run(files, sam).stats().reads_skipped_flag == 1);
    }
    SECTION("QC fail")
    {
        read.flag = 0x200;
        sam.reads = {read};
        REQUIRE(run(files, sam).stats().reads_skipped_flag == 1);
    }
    SECTION("duplicate")
    {
        read.flag = 0x400;
        sam.reads = {read};
        REQUIRE(run(files, sam).stats().reads_skipped_flag == 1);
    }
    SECTION("paired but not a proper pair")
    {
        read.flag = 0x41;
        sam.reads = {read};
        const SampleCounts counts = run(files, sam);
        REQUIRE(counts.bases_used() == 0);
        REQUIRE(counts.stats().reads_skipped_improper == 1);
    }
    SECTION("proper pair whose mate is elsewhere")
    {
        read.flag = 0x43;
        read.mate_contig = "chr02";
        read.mpos = 1;
        sam.reads = {read};
        REQUIRE(run(files, sam).bases_used() == 3);
    }
    SECTION("MAPQ below the threshold")
    {
        read.mapq = 19;
        sam.reads = {read};
        const SampleCounts counts = run(files, sam);
        REQUIRE(counts.bases_used() == 0);
        REQUIRE(counts.stats().reads_skipped_mapq == 1);
    }
    SECTION("MAPQ at the threshold")
    {
        read.mapq = 20;
        sam.reads = {read};
        REQUIRE(run(files, sam).bases_used() == 3);
    }
    SECTION("MAPQ 255 means unavailable")
    {
        read.mapq = 255;
        sam.reads = {read};
        const SampleCounts counts = run(files, sam);
        REQUIRE(counts.bases_used() == 0);
        REQUIRE(counts.stats().reads_skipped_mapq == 1);
    }
}

TEST_CASE(
    "pileup_bam skips bases below the base quality threshold",
    "[parentage][pileup]")
{
    FileFixture files;
    SamFile sam;
    SamRead read = allele0_read;

    SECTION("quality 19 at the first site")
    {
        read.qual = "IIIIIIIII4IIIIIIIIIIIIIIIIIIII";
        sam.reads = {read};
        const SampleCounts counts = run(files, sam);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{0, 1, 1, 0, 0});
        REQUIRE(counts.stats().bases_skipped_baseq == 1);
    }
    SECTION("quality 20 at the first site")
    {
        read.qual = "IIIIIIIII5IIIIIIIIIIIIIIIIIIII";
        sam.reads = {read};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{1, 1, 1, 0, 0});
    }
}

TEST_CASE(
    "pileup_bam follows the CIGAR to find the base under a site",
    "[parentage][pileup]")
{
    FileFixture files;
    SamFile sam;

    SECTION("deletion over a site")
    {
        sam.reads
            = {{.cigar = "5M10D15M",
                .seq = "NNNNNNNNNCNNNNNNNNNG",
                .qual = std::string(20, 'I')}};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{0, 1, 1, 0, 0});
    }
    SECTION("insertion before the sites")
    {
        sam.reads
            = {{.cigar = "9M2I21M",
                .seq = "NNNNNNNNNNNANNNNNNNNNCNNNNNNNNNG",
                .qual = std::string(32, 'I')}};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{1, 1, 1, 0, 0});
    }
    SECTION("soft clip before the sites")
    {
        sam.reads
            = {{.cigar = "3S30M",
                .seq = "NNNNNNNNNNNNANNNNNNNNNCNNNNNNNNNG",
                .qual = std::string(33, 'I')}};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{1, 1, 1, 0, 0});
    }
    SECTION("hard clip is not in the sequence")
    {
        SamRead read = allele0_read;
        read.cigar = "3H30M";
        sam.reads = {read};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{1, 1, 1, 0, 0});
    }
    SECTION("reference skip over a site")
    {
        sam.reads
            = {{.cigar = "5M10N15M",
                .seq = "NNNNNNNNNCNNNNNNNNNG",
                .qual = std::string(20, 'I')}};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{0, 1, 1, 0, 0});
    }
}

TEST_CASE(
    "pileup_bam counts a site once when both mates cover it",
    "[parentage][pileup]")
{
    FileFixture files;
    SamFile sam;
    // First mate: positions 1-30. Second mate: positions 11-40 carrying
    // allele0 at sites 20, 30 and 35.
    const SamRead first{
        .qname = "p",
        .flag = 0x63,
        .pos = 1,
        .cigar = "30M",
        .mate_contig = "=",
        .mpos = 11,
        .tlen = 40,
        .seq = "NNNNNNNNNANNNNNNNNNCNNNNNNNNNG",
        .qual = q30};
    const SamRead second{
        .qname = "p",
        .flag = 0x93,
        .pos = 11,
        .cigar = "30M",
        .mate_contig = "=",
        .mpos = 1,
        .tlen = -40,
        .seq = "NNNNNNNNNCNNNNNNNNNGNNNNTNNNNN",
        .qual = q30};

    SECTION("agreeing mates count once")
    {
        sam.reads = {first, second};
        const SampleCounts counts = run(files, sam);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{1, 1, 1, 1, 0});
        REQUIRE(counts.stats().bases_skipped_baseq == 2);
    }
    SECTION("dedup disabled counts both mates")
    {
        sam.reads = {first, second};
        PileupOptions options;
        options.dedup_overlaps = false;
        const SampleCounts counts = run(files, sam, options);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{1, 2, 2, 1, 0});
        REQUIRE(counts.stats().bases_skipped_baseq == 0);
    }
    SECTION("disagreeing mates keep the higher-quality base")
    {
        SamRead other_allele = second;
        other_allele.seq = "NNNNNNNNNTNNNNNNNNNGNNNNTNNNNN";   // T at site 20
        other_allele.qual = "IIIIIIIII5IIIIIIIIIIIIIIIIIIII";  // quality 20
        sam.reads = {first, other_allele};
        const SampleCounts counts = run(files, sam);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{1, 1, 1, 1, 0});
        REQUIRE(n1_of(counts) == std::vector<uint8_t>{0, 0, 0, 0, 0});
    }
    SECTION("disagreeing mates drop the lower-quality base even on mate 1")
    {
        SamRead weak_first = first;
        weak_first.qual = "IIIIIIIIIIIIIIIIIII5IIIIIIIIII";  // 20 at site 20
        SamRead other_allele = second;
        other_allele.seq = "NNNNNNNNNTNNNNNNNNNGNNNNTNNNNN";
        sam.reads = {weak_first, other_allele};
        const SampleCounts counts = run(files, sam);
        REQUIRE(n0_of(counts) == std::vector<uint8_t>{1, 0, 1, 1, 0});
        REQUIRE(n1_of(counts) == std::vector<uint8_t>{0, 1, 0, 0, 0});
    }
    SECTION("second mate contained in the first counts once")
    {
        SamRead contained = second;
        contained.cigar = "10M";
        contained.tlen = -30;
        contained.seq = "NNNNNNNNNC";
        contained.qual = std::string(10, 'I');
        SamRead left = first;
        left.tlen = 30;
        sam.reads = {left, contained};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{1, 1, 1, 0, 0});
    }
    SECTION("mates starting at the same position count once")
    {
        SamRead left = first;
        left.mpos = 1;
        left.tlen = 30;
        SamRead right = first;
        right.flag = 0x93;
        right.mpos = 1;
        right.tlen = -30;
        sam.reads = {left, right};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{1, 1, 1, 0, 0});
    }
    SECTION("unpaired reads at the same position both count")
    {
        sam.reads = {allele0_read, allele0_read};
        REQUIRE(n0_of(run(files, sam)) == std::vector<uint8_t>{2, 2, 2, 0, 0});
    }
}

TEST_CASE("pileup_bam saturates a site count at 255", "[parentage][pileup]")
{
    FileFixture files;
    SamFile sam;
    sam.reads = repeated_reads(allele0_read, 300);

    const SampleCounts counts = run(files, sam);

    REQUIRE(counts.n0()[0] == 255);
    REQUIRE(counts.stats().reads_used == 300);
}

TEST_CASE(
    "pileup_bam checks the BAM contigs against the panel",
    "[parentage][pileup]")
{
    FileFixture files;
    SamFile sam;

    SECTION("reads on a contig without panel sites are counted as off-panel")
    {
        sam.contigs = {{"chr01", 1000}, {"chr02", 500}, {"chrM", 100}};
        SamRead read = allele0_read;
        read.contig = "chrM";
        sam.reads = {read};
        const SampleCounts counts = run(files, sam);
        REQUIRE(counts.bases_used() == 0);
        REQUIRE(counts.stats().reads_off_panel == 1);
        REQUIRE(counts.stats().reads_used == 0);
    }
    SECTION("panel contig missing from the BAM")
    {
        sam.contigs = {{"chr01", 1000}};
        REQUIRE_THROWS_WITH(
            run(files, sam),
            ContainsSubstring("Panel contig chr02 is not present in"));
    }
    SECTION("contig length differs")
    {
        sam.contigs = {{"chr01", 1000}, {"chr02", 501}};
        REQUIRE_THROWS_WITH(
            run(files, sam),
            ContainsSubstring(
                "Contig chr02 has length 500 in the panel but 501"));
    }
}

TEST_CASE(
    "pileup_bam rejects a BAM that is not coordinate sorted",
    "[parentage][pileup]")
{
    FileFixture files;
    SamFile sam;
    SamRead later = allele0_read;
    later.pos = 20;

    SECTION("position goes backwards")
    {
        sam.reads = {later, allele0_read};
        REQUIRE_THROWS_WITH(
            run(files, sam),
            ContainsSubstring(
                "not coordinate sorted (chr01:1 follows chr01:20)"));
    }
    SECTION("contig goes backwards")
    {
        SamRead other = allele0_read;
        other.contig = "chr02";
        sam.reads = {other, allele0_read};
        REQUIRE_THROWS_WITH(
            run(files, sam), ContainsSubstring("not coordinate sorted"));
    }
    SECTION("file does not exist")
    {
        REQUIRE_THROWS_WITH(
            pileup_bam(files.missing_path(".bam"), "S", make_sites(), {}),
            ContainsSubstring("Failed to open BAM"));
    }
}
