#include "parentage/report.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace vcfbox
{
namespace
{
std::vector<PairScore> select_top_k(
    const Eigen::MatrixXd& ll, const Eigen::MatrixXd& posterior, int top_k)
{
    const auto n_pairs = ll.size();
    const auto k = std::min<Eigen::Index>(top_k, n_pairs);
    std::vector<Eigen::Index> order(static_cast<size_t>(n_pairs));
    for (Eigen::Index i = 0; i < n_pairs; ++i)
    {
        order[static_cast<size_t>(i)] = i;
    }
    std::partial_sort(
        order.begin(),
        order.begin() + k,
        order.end(),
        [&](Eigen::Index a, Eigen::Index b) {
            return *(ll.data() + a) > *(ll.data() + b);
        });

    std::vector<PairScore> out;
    out.reserve(static_cast<size_t>(k));
    const Eigen::Index rows = ll.rows();
    for (Eigen::Index r = 0; r < k; ++r)
    {
        const Eigen::Index idx = order[static_cast<size_t>(r)];
        const Eigen::Index col = idx / rows;
        const Eigen::Index row = idx % rows;
        out.push_back(
            {row, col, ll(row, col), posterior(row, col)});
    }
    return out;
}
}  // namespace

SampleReport build_report(const Eigen::MatrixXd& ll, int top_k)
{
    if (ll.size() == 0)
    {
        throw std::runtime_error("build_report: empty likelihood matrix");
    }
    const double ll_max = ll.maxCoeff();
    const auto shifted = (ll.array() - ll_max).exp();
    const double z = shifted.sum();
    const Eigen::MatrixXd posterior = shifted / z;

    SampleReport report;
    report.top_k = select_top_k(ll, posterior, std::max(1, top_k));
    report.best = report.top_k.front();
    report.m_marg_post = posterior.row(report.best.m).sum();
    report.f_marg_post = posterior.col(report.best.f).sum();
    return report;
}

Call classify_call(const SampleReport& report, double threshold)
{
    if (report.best.posterior >= threshold)
    {
        return Call::Accept;
    }
    if (std::max(report.m_marg_post, report.f_marg_post) >= threshold)
    {
        return Call::Partial;
    }
    return Call::Reject;
}

std::string_view to_string(Call call)
{
    switch (call)
    {
    case Call::Accept:
        return "accept";
    case Call::Partial:
        return "partial";
    case Call::Reject:
        return "reject";
    }
    return "unknown";
}
}  // namespace vcfbox
