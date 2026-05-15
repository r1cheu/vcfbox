#include "parentage/likelihood.h"

#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace vcfbox
{
namespace
{
static_assert(
    sizeof(AlleleCount) == 2 * sizeof(uint32_t),
    "AlleleCount must be a packed pair of uint32_t for zero-copy Map");

using CountRowM = Eigen::Matrix<uint32_t, Eigen::Dynamic, 2, Eigen::RowMajor>;

struct PerSiteLogLik
{
    Eigen::VectorXd rr;
    Eigen::VectorXd aa;
    Eigen::VectorXd ra;
};

PerSiteLogLik build_per_site(
    std::span<const AlleleCount> counts,
    double error_rate)
{
    const auto n = static_cast<Eigen::Index>(counts.size());
    const Eigen::Map<const CountRowM> mat(
        reinterpret_cast<const uint32_t*>(counts.data()), n, 2);

    const Eigen::VectorXd n_ref = mat.col(0).cast<double>();
    const Eigen::VectorXd n_alt = mat.col(1).cast<double>();

    const double log_correct = std::log(1.0 - error_rate);
    const double log_error = std::log(error_rate / 3.0);
    const double log_het
        = std::log((0.5 * (1.0 - error_rate)) + (error_rate / 6.0));

    PerSiteLogLik out;
    out.rr = (n_ref * log_correct) + (n_alt * log_error);
    out.aa = (n_alt * log_correct) + (n_ref * log_error);
    out.ra = (n_ref + n_alt) * log_het;
    return out;
}

void check_dims(
    Eigen::Index n,
    const Eigen::MatrixXd& m0,
    const Eigen::MatrixXd& m1,
    const Eigen::MatrixXd& p0,
    const Eigen::MatrixXd& p1)
{
    if (m0.rows() != n || m1.rows() != n || p0.rows() != n || p1.rows() != n)
    {
        throw std::runtime_error(
            "likelihood: bitmatrix row count does not match site count");
    }
    if (m0.cols() != m1.cols())
    {
        throw std::runtime_error("likelihood: M0/M1 column count mismatch");
    }
    if (p0.cols() != p1.cols())
    {
        throw std::runtime_error("likelihood: P0/P1 column count mismatch");
    }
}

}  // namespace

Eigen::MatrixXd compute_pair_loglik(
    std::span<const AlleleCount> counts,
    const Eigen::MatrixXd& m0,
    const Eigen::MatrixXd& m1,
    const Eigen::MatrixXd& p0,
    const Eigen::MatrixXd& p1,
    double error_rate)
{
    const auto n = static_cast<Eigen::Index>(counts.size());
    check_dims(n, m0, m1, p0, p1);

    const auto per_site = build_per_site(counts, error_rate);

    const auto p_m0
        = (per_site.rr.asDiagonal() * p0) + (per_site.ra.asDiagonal() * p1);
    const auto p_m1
        = (per_site.aa.asDiagonal() * p1) + (per_site.ra.asDiagonal() * p0);

    return (m0.transpose() * p_m0) + (m1.transpose() * p_m1);
}
}  // namespace vcfbox
