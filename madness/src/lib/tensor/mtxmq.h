/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
#ifndef MADNESS_TENSOR_MTXMQ_H__INCLUDED
#define MADNESS_TENSOR_MTXMQ_H__INCLUDED

#include <madness_config.h>

#ifdef HAVE_IBMBGP
#include <linalg/cblas.h>
#endif

namespace madness {

    /// Matrix = Matrix transpose * matrix ... reference implementation

    /// Does \c C=AT*B whereas mTxm does C=C+AT*B.  It also supposed
    /// to be fast which it achieves thru restrictions
    ///   * All dimensions even
    ///   * All pointers aligned
    /// \code
    ///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
    /// \endcode
    template <typename aT, typename bT, typename cT>
    void mTxmq(long dimi, long dimj, long dimk,
               cT* restrict c, const aT* a, const bT* b) {

        //std::cout << "IN GENERIC mTxmq " << tensor_type_names[TensorTypeData<aT>::id] << " " << tensor_type_names[TensorTypeData<bT>::id] << " " << tensor_type_names[TensorTypeData<cT>::id] << "\n";

        for (long i=0; i<dimi; ++i,c+=dimj,++a) {
            for (long j=0; j<dimj; ++j) c[j] = 0.0;
            const aT *aik_ptr = a;
            for (long k=0; k<dimk; ++k,aik_ptr+=dimi) {
                aT aki = *aik_ptr;
                for (long j=0; j<dimj; ++j) {
                    c[j] += aki*b[k*dimj+j];
                }
            }
        }

    }

#ifdef HAVE_IBMBGP
    template <>
    inline void mTxmq(long ni, long nj, long nk, double* restrict c, const double* a, const double* b) {
	double one=1.0;
	double zero=0.0;
	dgemm_("n","t",&nj,&ni,&nk,&one,b,&nj,a,&ni,&zero,c,&nj,1,1);
    }

    template <>
    inline void mTxmq(long ni, long nj, long nk, double_complex* restrict c, const double_complex* a, const double_complex* b) {
	double_complex one=1.0;
	double_complex zero=0.0;
	zgemm_("n","t",&nj,&ni,&nk,&one,b,&nj,a,&ni,&zero,c,&nj,1,1);
    }

#elif defined(X86_64) && !defined(DISABLE_SSE3)
    template <>
    void mTxmq(long dimi, long dimj, long dimk,
               double* restrict c, const double* a, const double* b);

    template <>
    void mTxmq(long dimi, long dimj, long dimk,
               double_complex* restrict c, const double_complex* a, const double_complex* b);

#ifndef __INTEL_COMPILER
    template <>
    void mTxmq(long dimi, long dimj, long dimk,
               double_complex* restrict c, const double_complex* a, const double* b);
#endif

#elif defined(X86_32)
    template <>
    void mTxmq(long dimi, long dimj, long dimk,
               double* restrict c, const double* a, const double* b);
#endif

}

#endif // MADNESS_TENSOR_MTXMQ_H__INCLUDED
