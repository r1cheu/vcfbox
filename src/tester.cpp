#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C"
{
#include <htslib/vcf.h>
}

// ----------------------------------------------------------------------------
//  要测试的核心代码
// ----------------------------------------------------------------------------

using SamplePair = std::pair<std::string, std::string>;

std::vector<SamplePair> parse_sample_pairs(const std::string& file_path)
{
    std::vector<SamplePair> pairs;
    std::ifstream file(file_path);
    if (!file)
    {
        throw std::runtime_error("Cannot open sample pairs file: " + file_path);
    }
    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        std::istringstream iss(line);
        std::string s1, s2;
        if (iss >> s1 >> s2)
        {
            pairs.emplace_back(s1, s2);
        }
    }
    return pairs;
}

void check_sample_consistence(
    std::string_view vcf_path,
    const std::vector<SamplePair>& sample_pairs)
{
    htsFile* vcf_file = bcf_open(vcf_path.data(), "r");
    if (vcf_file == nullptr)
    {
        throw std::runtime_error(
            "Could not open VCF file: " + std::string(vcf_path));
    }
    bcf_hdr_t* header = bcf_hdr_read(vcf_file);
    if (header == nullptr)
    {
        bcf_close(vcf_file);
        throw std::runtime_error(
            "Could not read VCF header from: " + std::string(vcf_path));
    }

    std::set<std::string> all_samples;
    for (int i = 0; i < bcf_hdr_nsamples(header); ++i)
    {
        all_samples.insert(header->samples[i]);
    }
    bcf_hdr_destroy(header);
    bcf_close(vcf_file);

    std::vector<std::string> missing_samples;
    for (const auto& pair : sample_pairs)
    {
        if (all_samples.find(pair.first) == all_samples.end())
        {
            missing_samples.push_back(pair.first);
        }
        if (all_samples.find(pair.second) == all_samples.end())
        {
            missing_samples.push_back(pair.second);
        }
    }

    if (!missing_samples.empty())
    {
        std::string missing_str;
        for (const auto& s : missing_samples)
        {
            missing_str += s + ", ";
        }
        throw std::runtime_error("Samples not found in VCF: " + missing_str);
    }
}

void extract_genotypes(
    const std::string& vcf_path,
    const std::vector<SamplePair>& sample_pairs,
    const std::string& out_path)
{
    check_sample_consistence(vcf_path, sample_pairs);

    htsFile* vcf_file = bcf_open(vcf_path.c_str(), "r");
    if (vcf_file == nullptr)
    {
        throw std::runtime_error("Could not open VCF file: " + vcf_path);
    }
    bcf_hdr_t* header = bcf_hdr_read(vcf_file);
    if (header == nullptr)
    {
        bcf_close(vcf_file);
        throw std::runtime_error("Could not read VCF header from: " + vcf_path);
    }

    std::unordered_map<std::string, int> sample_to_idx;
    for (int i = 0; i < bcf_hdr_nsamples(header); ++i)
    {
        sample_to_idx[header->samples[i]] = i;
    }

    htsFile* output_file = hts_open(out_path.c_str(), "w");
    if (output_file == nullptr)
    {
        bcf_hdr_destroy(header);
        bcf_close(vcf_file);
        throw std::runtime_error("Could not open output file: " + out_path);
    }

    bcf_hdr_t* output_header = bcf_hdr_init("w");
    for (int i = 0; i < header->nhrec; ++i)
    {
        if (header->hrec[i]->type == BCF_HL_CTG)
        {
            bcf_hdr_add_hrec(output_header, bcf_hrec_dup(header->hrec[i]));
        }
    }
    bcf_hrec_t* gt_hrec
        = bcf_hdr_get_hrec(header, BCF_HL_FMT, "ID", "GT", NULL);
    if (gt_hrec)
    {
        bcf_hdr_add_hrec(output_header, bcf_hrec_dup(gt_hrec));
    }
    else
    {
        bcf_hdr_append(
            output_header,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    }
    for (const auto& pair : sample_pairs)
    {
        bcf_hdr_add_sample(
            output_header, (pair.first + "_" + pair.second).c_str());
    }
    bcf_hdr_add_sample(output_header, NULL);

    if (bcf_hdr_write(output_file, output_header) != 0)
    {
        throw std::runtime_error("Failed to write output header");
    }

    bcf1_t* in_rec = bcf_init();
    bcf1_t* out_rec = bcf_init();
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;

    while (bcf_read(vcf_file, header, in_rec) == 0)
    {
        bcf_unpack(in_rec, BCF_UN_ALL);

        if (in_rec->n_allele > 2)
        {
            continue;
        }

        if (bcf_get_genotypes(header, in_rec, &gt_arr, &ngt_arr) <= 0)
        {
            continue;
        }

        std::vector<int32_t> out_gts;
        out_gts.reserve(sample_pairs.size() * 2);

        for (const auto& pair : sample_pairs)
        {
            int f_idx = sample_to_idx[pair.first];
            int m_idx = sample_to_idx[pair.second];
            int32_t* f_gt = &gt_arr[f_idx * 2];
            int32_t* m_gt = &gt_arr[m_idx * 2];

            if (bcf_gt_is_missing(f_gt[0]) || bcf_gt_is_missing(f_gt[1])
                || bcf_gt_is_missing(m_gt[0]) || bcf_gt_is_missing(m_gt[1])
                || bcf_gt_allele(f_gt[0]) != bcf_gt_allele(f_gt[1])
                || bcf_gt_allele(m_gt[0]) != bcf_gt_allele(m_gt[1]))
            {
                out_gts.push_back(bcf_gt_missing);
                out_gts.push_back(bcf_gt_missing);
            }
            else
            {
                out_gts.push_back(bcf_gt_unphased(bcf_gt_allele(f_gt[0])));
                out_gts.push_back(bcf_gt_unphased(bcf_gt_allele(m_gt[0])));
            }
        }

        bcf_clear(out_rec);
        out_rec->rid = bcf_hdr_name2id(
            output_header, bcf_hdr_id2name(header, in_rec->rid));
        out_rec->pos = in_rec->pos;
        bcf_update_id(output_header, out_rec, in_rec->d.id);
        bcf_update_alleles(
            output_header,
            out_rec,
            (const char**)in_rec->d.allele,
            in_rec->n_allele);
        out_rec->qual = in_rec->qual;

        bcf_update_genotypes(
            output_header, out_rec, out_gts.data(), out_gts.size());

        if (bcf_write(output_file, output_header, out_rec) != 0)
        {
            throw std::runtime_error("Failed to write VCF record");
        }
    }

    free(gt_arr);
    bcf_destroy(in_rec);
    bcf_destroy(out_rec);
    bcf_hdr_destroy(header);
    bcf_hdr_destroy(output_header);
    bcf_close(vcf_file);
    if (hts_close(output_file) != 0)
    {
        // Handle error
    }
}

// ----------------------------------------------------------------------------
//  测试框架
// ----------------------------------------------------------------------------

/**
 * @brief 创建测试所需的VCF和样本对文件
 */
void create_test_files()
{
    // 1. 创建 test.vcf
    std::ofstream vcf_file("test.vcf");
    vcf_file
        << "##fileformat=VCFv4.2\n"
        << "##contig=<ID=chr1,length=1000>\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\tP3\t"
           "P4\tP5\n"
        << "chr1\t100\trs1\tA\tC\t.\tPASS\t.\tGT\t0/0\t1/1\t./.\t0/1\t0/0\n"
        << "chr1\t200\trs2\tA\tC,G\t.\tPASS\t.\tGT\t0/0\t1/1\t1/2\t0/1\t0/0\n"
        << "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT\t1/1\t0/0\t./.\t0/1\t1/1\n";
    vcf_file.close();

    // 2. 创建 pairs.txt
    std::ofstream pairs_file("pairs.txt");
    pairs_file << "P1\tP2\n"
               << "P1\tP3\n"
               << "P1\tP4\n"
               << "P1\tP5\n";
    pairs_file.close();
}

/**
 * @brief 读取文件全部内容到字符串
 */
std::string read_file_content(const std::string& path)
{
    std::ifstream file(path);
    if (!file)
    {
        return "";
    }
    return std::string(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());
}

int main(int argc, char* argv[])
{
    try
    {
        // 步骤 1: 创建测试文件
        std::cout << "1. 创建测试文件 (test.vcf, pairs.txt)..." << std::endl;
        create_test_files();

        // 步骤 2: 定义期望的输出内容
        std::cout << "2. 定义期望的输出..." << std::endl;
        const std::string expected_output
            = "##fileformat=VCFv4.2\n"
              "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
              "##contig=<ID=chr1,length=1000>\n"
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1_"
              "P2\tP1_P3\tP1_P4\tP1_P5\n"
              "chr1\t100\trs1\tA\tC\t.\t.\t.\tGT\t0/1\t./.\t./.\t0/0\n"
              "chr1\t300\trs3\tG\tT\t.\t.\t.\tGT\t1/0\t./.\t./.\t1/1\n";

        std::cout << "3. 执行 extract_genotypes..." << "\n";
        auto pairs = parse_sample_pairs("pairs.txt");
        extract_genotypes("test.vcf", pairs, "output.vcf");

        std::cout << "4. 比对实际输出与期望输出..." << '\n';
        std::string actual_output = read_file_content("output.vcf");

        std::string cleaned_actual;
        std::stringstream ss(actual_output);
        std::string line;
        while (std::getline(ss, line, '\n'))
        {
            if (!line.empty() && line.back() == '\r')
            {
                line.pop_back();
            }
            cleaned_actual += line + "\n";
        }

        if (cleaned_actual == expected_output)
        {
            std::cout << "\n✅ 测试通过！输出与期望完全一致。\n" << std::endl;
        }
        else
        {
            std::cout << "\n❌ 测试失败！输出与期望不符。\n" << std::endl;
            std::cout << "--- 期望输出 ---\n"
                      << expected_output << "-------------------\n"
                      << '\n';
            std::cout << "--- 实际输出 ---\n"
                      << actual_output << "-------------------\n"
                      << '\n';
            return 1;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "测试过程中发生错误: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
