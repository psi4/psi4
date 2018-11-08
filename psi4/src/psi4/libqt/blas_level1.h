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

/*!
 * \file
 * \brief Function signatures for wrappers to BLAS level 1 subroutines
 * \ingroup QT
 */

#pragma once

#include <cstddef>

#include "psi4/pragma.h"

namespace psi {
PSI_API
void C_DROT(size_t ntot, double* x, int incx, double* y, int incy, double costheta, double sintheta);
PSI_API
void C_DSWAP(size_t length, double* x, int incx, double* y, int inc_y);
PSI_API
void C_DSCAL(size_t len, double alpha, double* vec, int inc);
PSI_API
void C_DCOPY(size_t length, double* x, int inc_x, double* y, int inc_y);
PSI_API
void C_DAXPY(size_t length, double a, double* x, int inc_x, double* y, int inc_y);
PSI_API
double C_DDOT(size_t n, double* X, int inc_x, double* Y, int inc_y);
PSI_API
double C_DNRM2(size_t n, double* X, int inc_x);
PSI_API
double C_DASUM(size_t n, double* X, int inc_x);
PSI_API
size_t C_IDAMAX(size_t n, double* X, int inc_x);
}  // namespace psi
