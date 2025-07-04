#pragma once

#include <memory>
extern "C"
{
#include <htslib/vcf.h>
}

struct HtsFileDeleter
{
    void operator()(htsFile* p) const
    {
        if (p != nullptr)
        {
            hts_close(p);
        }
    }
};

struct BcfHdrDeleter
{
    void operator()(bcf_hdr_t* p) const
    {
        if (p != nullptr)
        {
            bcf_hdr_destroy(p);
        }
    }
};

struct BcfRecDeleter
{
    void operator()(bcf1_t* p) const
    {
        if (p != nullptr)
        {
            bcf_destroy(p);
        }
    }
};

class Genotypes
{
   public:
    int32_t* p_ = nullptr;
    int n_ = 0;

    ~Genotypes()
    {
        if (p_ != nullptr)
        {
            free(p_);
        }
    }
    // 禁止拷贝，允许移动
    Genotypes() = default;
    Genotypes(const Genotypes&) = delete;
    Genotypes& operator=(const Genotypes&) = delete;
    Genotypes(Genotypes&&) = default;
    Genotypes& operator=(Genotypes&&) = default;
};
using HtsFile = std::unique_ptr<htsFile, HtsFileDeleter>;
using BcfHdr = std::unique_ptr<bcf_hdr_t, BcfHdrDeleter>;
using BcfRec = std::unique_ptr<bcf1_t, BcfRecDeleter>;
