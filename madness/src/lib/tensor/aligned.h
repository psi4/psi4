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
#ifndef MADNESS_TENSOR_ALIGNED_H__INCLUDED
#define MADNESS_TENSOR_ALIGNED_H__INCLUDED

/*!
  \file tensor/aligned.h
  \brief Provides routines for internal use optimized for aligned data

  This stuff used to be implemented in assembly but it is too much
  effort keeping that working especially for multiple compilers.
*/

#include <madness_config.h>
#include <tensor/tensor.h>
#include <cstring>

namespace madness {

    template <typename T>
    static
    inline
    void aligned_zero(long n, T* a) {
#ifdef HAVE_MEMSET
        // A hand coded SSE2 loop is faster only for data in the L1 cache
        std::memset((void *) a, 0, n*sizeof(T));
#else
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        for (long i=0; i<n4; i+=4,a+=4) {
            a[0] = 0;
            a[1] = 0;
            a[2] = 0;
            a[3] = 0;
        }
        for (long i=0; i<rem; ++i) *a++ = 0;
#endif
    }

    template <typename T, typename Q>
    static
    inline
    void aligned_axpy(long n, T* restrict a, const T* restrict b, Q s) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        for (long i=0; i<n4; i+=4,a+=4,b+=4) {
            a[0] += s*b[0];
            a[1] += s*b[1];
            a[2] += s*b[2];
            a[3] += s*b[3];
        }
        for (long i=0; i<rem; ++i) *a++ += s * *b++;
    }

    template <typename T, typename Q>
    static
    inline
    void aligned_add(long n, T* restrict a, const Q* restrict b) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        for (long i=0; i<n4; i+=4,a+=4,b+=4) {
            a[0] += b[0];
            a[1] += b[1];
            a[2] += b[2];
            a[3] += b[3];
        }
        for (long i=0; i<rem; ++i) *a++ += *b++;
    }

    template <typename T, typename Q>
    static
    inline
    void aligned_sub(long n, T* restrict a, const Q* restrict b) {
        long n4 = (n>>2)<<2;
        long rem = n-n4;
        for (long i=0; i<n4; i+=4,a+=4,b+=4) {
            a[0] -= b[0];
            a[1] -= b[1];
            a[2] -= b[2];
            a[3] -= b[3];
        }
        for (long i=0; i<rem; ++i) *a++ -= *b++;
    }
}

#endif // MADNESS_TENSOR_ALIGNED_H__INCLUDED
