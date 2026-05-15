#pragma once

#include <Eigen/Core>
#include <span>

#include "parentage/pileup.h"

namespace vcfbox
{
Eigen::MatrixXd compute_pair_loglik(
    std::span<const AlleleCount> counts,
    const Eigen::MatrixXd& m0,
    const Eigen::MatrixXd& m1,
    const Eigen::MatrixXd& p0,
    const Eigen::MatrixXd& p1,
    double error_rate);
}  // namespace vcfbox
