/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \defgroup QT libqt: The Quantum-Trio Miscellaneous Library */

/*!
 \file
 \ingroup QT
 \brief The PSI3 BLAS1 interface routines

 Interface to the BLAS routines

 C. David Sherrill
 Anna I. Krylov

 May 1998

 NOTE: Refactored by Rob Parrish on 1/24/2010
 This file now contains all relevant BLAS1
 routines, with provisions made for >2^31
 elements (size_t sizes).

 All BLAS2 and BLAS3 routines are now wrapped
 and are in blas_inftc23.cc

*/

#include "blas_level1.h"

#include <cstddef>
#include <climits>
#include <cmath>

#ifdef USING_LAPACK_MKL
#include "mkl_cblas.h"
#else
#include "cblas.h"
#endif

namespace psi {
void C_DSWAP(size_t length, double *x, int inc_x, double *y, int inc_y) {
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *x_s = &x[block * inc_x * (size_t)INT_MAX];
        double *y_s = &y[block * inc_y * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        cblas_dswap(length_s, x_s, inc_x, y_s, inc_y);
    }
}

void C_DAXPY(size_t length, double a, double *x, int inc_x, double *y, int inc_y) {
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *x_s = &x[block * inc_x * (size_t)INT_MAX];
        double *y_s = &y[block * inc_y * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        cblas_daxpy(length_s, a, x_s, inc_x, y_s, inc_y);
    }
}

void C_DCOPY(size_t length, double *x, int inc_x, double *y, int inc_y) {
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *x_s = &x[block * inc_x * (size_t)INT_MAX];
        double *y_s = &y[block * inc_y * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        cblas_dcopy(length_s, x_s, inc_x, y_s, inc_y);
    }
}

void C_DSCAL(size_t length, double alpha, double *vec, int inc) {
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *vec_s = &vec[block * inc * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        cblas_dscal(length_s, alpha, vec_s, inc);
    }
}

void C_DROT(size_t length, double *x, int inc_x, double *y, int inc_y, double costheta, double sintheta) {
    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *x_s = &x[block * inc_x * (size_t)INT_MAX];
        double *y_s = &y[block * inc_y * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        cblas_drot(length_s, x_s, inc_x, y_s, inc_y, costheta, sintheta);
    }
}

double C_DDOT(size_t length, double *x, int inc_x, double *y, int inc_y) {
    if (length == 0) return 0.0;

    double reg = 0.0;

    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *x_s = &x[block * inc_x * (size_t)INT_MAX];
        double *y_s = &y[block * inc_y * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        reg += cblas_ddot(length_s, x_s, inc_x, y_s, inc_y);
    }

    return reg;
}

double C_DNRM2(size_t length, double *x, int inc_x) {
    if (length == 0) return 0.0;

    double reg = 0.0;

    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *x_s = &x[block * inc_x * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        reg += cblas_dnrm2(length_s, x_s, inc_x);
    }

    return reg;
}

double C_DASUM(size_t length, double *x, int inc_x) {
    if (length == 0) return 0.0;

    double reg = 0.0;

    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *x_s = &x[block * inc_x * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        reg += cblas_dasum(length_s, x_s, inc_x);
    }

    return reg;
}

size_t C_IDAMAX(size_t length, double *x, int inc_x) {
    if (length == 0) return 0L;

    size_t reg = 0L;
    size_t reg2 = 0L;

    int big_blocks = (int)(length / INT_MAX);
    int small_size = (int)(length % INT_MAX);
    for (int block = 0; block <= big_blocks; block++) {
        double *x_s = &x[block * inc_x * (size_t)INT_MAX];
        signed int length_s = (block == big_blocks) ? small_size : INT_MAX;
        reg2 = cblas_idamax(length_s, x_s, inc_x) + block * inc_x * (size_t)INT_MAX;
        if (std::fabs(x[reg]) > std::fabs(x[reg2])) reg = reg2;
    }

    return reg;
}
}  // namespace psi
