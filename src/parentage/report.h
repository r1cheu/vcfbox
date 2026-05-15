#pragma once

#include <Eigen/Core>
#include <string_view>
#include <vector>

namespace vcfbox
{
struct PairScore
{
    Eigen::Index m;
    Eigen::Index f;
    double loglik;
    double posterior;
};

struct SampleReport
{
    std::vector<PairScore> top_k;
    PairScore best;
    double m_marg_post;
    double f_marg_post;
};

SampleReport build_report(const Eigen::MatrixXd& ll, int top_k);

enum class Call
{
    Accept,
    Partial,
    Reject,
};

Call classify_call(const SampleReport& report, double threshold);
std::string_view to_string(Call call);
}  // namespace vcfbox
