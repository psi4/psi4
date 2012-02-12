#include <cstdlib>
#include <cmath>

#include "matrix.h"
#include "vector.h"
#include "orthog.h"

namespace psi {

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

//    overlap_->diagonalize(U, m);


}
}
