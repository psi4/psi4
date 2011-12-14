
#include "opimpl.h"
#include "laplace.h"
#include "class.h"
#include "tensor.h"
#include "tensorblock.h"
#include "quadratures.h"
#include "env.h"

using namespace yeti;
using namespace std;


#ifdef redefine_size_t
#define size_t custom_size_t
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
//  LaplaceTransformOp class
////////////////////////////////////////////////////////////////////////////////////////////////////


//-------------------------------------------------
//  Constructors and Initializers
//-------------------------------------------------

LaplaceTransformOp::LaplaceTransformOp(
    LaplaceTransform* lt,
    const double* evals,
    usi total_indices
) : evals_(evals),
    lt_(lt),
    total_indices_(total_indices)
{

}


LaplaceTransformOp::LaplaceTransformOp(
    LaplaceTransform* lt,
    const double* evals,
    usi total_indices,
    double efermi
) : evals_(evals),
    lt_(lt),
    efermi_(efermi),
    total_indices_(total_indices)
{

}


LaplaceTransformOp::~LaplaceTransformOp() {
}


//-------------------------------------------------
//  Virtual methods from ElementOp
//-------------------------------------------------

void
LaplaceTransformOp::retrieve(TensorBlock* block) const
{
    block->retrieve_verbatim();
}

void
LaplaceTransformOp::release(TensorBlock* block) const
{
    block->release_verbatim();
}

void
LaplaceTransformOp::configure(usi alpha) {
    alpha_ = alpha;
}

//-------------------------------------------------
//  element_op implementations
//-------------------------------------------------

void
LaplaceTransformOp::element_op(
    uli nblock,
    const uli* index_starts,
    const uli* sizes,
    double* data
)
{
    uli istart = index_starts[total_indices_ - 1];
    uli size = sizes[total_indices_ - 1];
    uli istop = istart + size;
    uli ntrans = nblock / size;
    double* dataptr = data;
    for(int x = 0; x < ntrans; ++x) {
        for(int p = istart; p < istop; ++p, ++dataptr) {
            double deltae = -fabs(evals_[p]);
            (*dataptr) *= exp(deltae * lt_->t[alpha_]);
        }
    }
}



////////////////////////////////////////////////////////////////////////////////////////////////////
//  LaplaceTransform class
////////////////////////////////////////////////////////////////////////////////////////////////////

LaplaceTransform::LaplaceTransform()
{
    nquad_ = 13;

    t = new double[nquad_];
    w = new double[nquad_];
    get_quadrature(nquad_, t, w);
}


LaplaceTransform::~LaplaceTransform()
{
    delete[] w;
    delete[] t;
}

