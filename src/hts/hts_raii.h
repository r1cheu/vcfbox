#pragma once

#include <cstdlib>
#include <memory>
extern "C"
{
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
}

struct hts_file_deleter
{
    void operator()(htsFile* p) const
    {
        if (p != nullptr)
        {
            hts_close(p);
        }
    }
};

struct hts_idx_deleter
{
    void operator()(hts_idx_t* p) const
    {
        if (p != nullptr)
        {
            hts_idx_destroy(p);
        }
    }
};

struct hts_itr_deleter
{
    void operator()(hts_itr_t* p) const
    {
        if (p != nullptr)
        {
            hts_itr_destroy(p);
        }
    }
};

struct bcf_hdr_deleter
{
    void operator()(bcf_hdr_t* p) const
    {
        if (p != nullptr)
        {
            bcf_hdr_destroy(p);
        }
    }
};

struct bcf1_deleter
{
    void operator()(bcf1_t* p) const
    {
        if (p != nullptr)
        {
            bcf_destroy(p);
        }
    }
};

struct sam_hdr_deleter
{
    void operator()(sam_hdr_t* p) const
    {
        if (p != nullptr)
        {
            sam_hdr_destroy(p);
        }
    }
};

struct bam1_deleter
{
    void operator()(bam1_t* p) const
    {
        if (p != nullptr)
        {
            bam_destroy1(p);
        }
    }
};

struct hts_free_deleter
{
    void operator()(void* p) const { std::free(p); }
};

using hts_file_ptr = std::unique_ptr<htsFile, hts_file_deleter>;
using hts_idx_ptr = std::unique_ptr<hts_idx_t, hts_idx_deleter>;
using hts_itr_ptr = std::unique_ptr<hts_itr_t, hts_itr_deleter>;
using bcf_hdr_ptr = std::unique_ptr<bcf_hdr_t, bcf_hdr_deleter>;
using bcf1_ptr = std::unique_ptr<bcf1_t, bcf1_deleter>;
using sam_hdr_ptr = std::unique_ptr<sam_hdr_t, sam_hdr_deleter>;
using bam1_ptr = std::unique_ptr<bam1_t, bam1_deleter>;

struct genotype_buffer
{
    std::unique_ptr<int32_t, hts_free_deleter> data;
    int capacity = 0;
};

inline int read_genotypes(
    const bcf_hdr_t* hdr, bcf1_t* rec, genotype_buffer& buf)
{
    int32_t* raw = buf.data.release();
    const int ret = bcf_get_genotypes(hdr, rec, &raw, &buf.capacity);
    buf.data.reset(raw);
    return ret;
}
