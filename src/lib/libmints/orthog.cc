#include <cstdlib>
#include <cmath>

#include "matrix.h"
#include "vector.h"
#include "orthog.h"

namespace psi {

namespace {

struct max_abs {
  bool operator() (const double& i, const double& j)
    { return std::abs(i)<std::abs(j); }
};

}

OverlapOrthog::OverlapOrthog(OrthogMethod method,
                             SharedMatrix overlap,
                             double lindep_tolerance,
                             int debug)
    : nlindep_(0), min_orthog_res_(0.0), max_orthog_res_(0.0)
{
    orthog_method_ = method;
    overlap_ = overlap;
    lindep_tol_ = lindep_tolerance;
    debug_ = debug;
    dim_ = overlap_->rowspi();
}

void OverlapOrthog::compute_overlap_eig(SharedMatrix overlap_eigvec, SharedVector isqrt_eigval, SharedVector sqrt_eigval)
{
    SharedMatrix U(new Matrix("U", overlap_->rowspi(), overlap_->colspi()));
    SharedVector m(new Vector(overlap_->colspi()));

    overlap_->diagonalize(U, m);

    double maxabs = *std::max_element(m->begin(), m->end(), max_abs());
    double s_tol = lindep_tol_ * maxabs;
    double minabs = maxabs;

    std::vector<double> m_sqrt(dim_.sum());
    std::vector<double> m_isqrt(dim_.sum());
    std::vector<int> m_index(dim_.sum());
    std::vector<int> nfunc(dim_.n());
    int nfunctotal = 0;

    nlindep_ = 0;
    for (Vector::const_iterator iter=m->begin();
         iter != m->end(); ++iter) {
        if (*iter > s_tol) {
            if (*iter < minabs)
                minabs = *iter;
            m_sqrt[nfunctotal] = sqrt(*iter);
            m_sqrt[nfunctotal] = 1.0/m_sqrt[nfunctotal];
            m_index[nfunctotal] = std::distance(m->begin(), iter);
            nfunc
        }
    }
}
}
