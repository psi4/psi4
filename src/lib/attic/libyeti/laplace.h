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

#ifndef LAPLACE_H
#define LAPLACE_H

#include "elementop.h"
#include <math.h>


namespace yeti {

class LaplaceTransform
{
    protected:

        int nquad_;

        double max_err_;

    public:

        /// Constructor that uses some default values for testing purposes
        LaplaceTransform();

        ~LaplaceTransform();

        /// Returns the number of quadrature points
        int nquad() { return nquad_; }

        /// The Gaussian quadrature weights
        double* w;

        /// The Gaussian quadrature abcissa
        double* t;


};

class LaplaceTransformOp :
        public ElementOp
{
    protected:
        const double* evals_;

        /// The number of indices in the tensor
        usi total_indices_;

        /// The underlying LaplaceTransform object
        LaplaceTransform* lt_;

        /// The current quadrature point to compute
        usi alpha_;

        double efermi_;

    public:
        using ElementOp::element_op;

        /* OLD COMMENT
          Constructs a LaplaceTransformOp object that uses evals and is capable of operating on a
          tensor with total_indices indices.  The constructed object only transforms the first trans_indices
          indices of the tensor.  If trans_indices = 0, all of the indices are transformed (this is the default
          if no value is specified for this parameter.

          @param evals  The eigenvalues (orbital energies) used in the transformation
          @param lt A pointer to the LaplaceTransform object to use
          @param total_indices The operator is constructed to act only on tensors with the specified number of indices
          @param trans_indices The operator only transforms this number of indices.  If 0 or not given, transform all
                               indices (i. e. trans_indices = total_indices)
         */


        LaplaceTransformOp(
            LaplaceTransform* lt,
            const double* evals,
            usi total_indices
        );

        LaplaceTransformOp(
            LaplaceTransform* lt,
            const double* evals,
            usi total_indices,
            double efermi
        );


        ~LaplaceTransformOp();

        void element_op(
            uli nblock,
            const uli* index_starts,
            const uli* sizes,
            double* data
        );

        void configure(usi alpha);

        void retrieve(TensorBlock* block) const;

        void release(TensorBlock* block) const;

};








} // End namespace yeti


#endif // LAPLACE_H
