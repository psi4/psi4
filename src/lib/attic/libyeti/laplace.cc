/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */


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

