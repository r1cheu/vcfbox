#include "parentage/parentage.h"

#include <stdexcept>

namespace vcfbox
{
void build_parentage_matrices(const ParentageMatrixOptions& /*options*/)
{
    throw std::runtime_error("parentage build-matrix is not implemented yet");
}

void test_parentage(const ParentageTestOptions& /*options*/)
{
    throw std::runtime_error("parentage test is not implemented yet");
}
}  // namespace vcfbox
