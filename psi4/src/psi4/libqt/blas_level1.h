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
/**
 * Swaps a vector with another vector.
 *
 * @param length Specifies the number of elements in vectors x and y.
 * @param x Array, DIMENSION at least (1 + (n-1)*abs(inc_x)).
 * @param inc_x Specifies the increment for the elements of x.
 * @param y Array, DIMENSION at least (1 + (n-1)*abs(inc_y)).
 * @param inc_y Specifies the increment for the elements of y.
 *
 * @ingroup QT
 */
PSI_API
void C_DSWAP(size_t length, double* x, int inc_x, double* y, int inc_y);

/*!
 * This function performs y = a * x + y.
 *
 * Steps every inc_x in x and every inc_y in y (normally both 1).
 *
 * \param length   length of arrays
 * \param a        scalar a to multiply vector x
 * \param x        vector x
 * \param inc_x    how many places to skip to get to next element in x
 * \param y        vector y
 * \param inc_y    how many places to skip to get to next element in y
 *
 * \ingroup QT
 */
PSI_API
void C_DAXPY(size_t length, double a, double* x, int inc_x, double* y, int inc_y);

/*!
 * This function copies x into y.
 *
 * Steps every inc_x in x and every inc_y in y (normally both 1).
 *
 * \param length  = length of array
 * \param x       = vector x
 * \param inc_x   = how many places to skip to get to next element in x
 * \param y       = vector y
 * \param inc_y   = how many places to skip to get to next element in y
 *
 * \ingroup QT
 */
PSI_API
void C_DCOPY(size_t length, double* x, int inc_x, double* y, int inc_y);

/*!
 * This function scales a vector by a real scalar.
 *
 * \param len   length of array
 * \param alpha scale factor
 * \param vec   vector to scale
 * \param inc   how many places to skip to get to next element in vec
 *
 * \ingroup QT
 */
PSI_API
void C_DSCAL(size_t len, double alpha, double* vec, int inc);

/*!
 *Calculates a plane Givens rotation for vectors x, y and
 * angle theta.  x = x*cos + y*sin, y = -x*sin + y*cos.
 *
 * \param ntot   length of x,y
 * \param x      vector x
 * \param inc_x  how many places to skip to get to the next element of x
 * \param y      vector Y
 * \param inc_y  how many places to skip to get to the next element of y
 * \param costheta cosine
 * \param sintheta sine
 *
 * \ingroup QT
 */
PSI_API
void C_DROT(size_t ntot, double* x, int inc_x, double* y, int inc_y, double costheta, double sintheta);

/*!
 * This function returns the dot product of two vectors, x and y.
 *
 * \param n      Number of elements in x and y.
 * \param x      A pointer to the beginning of the data in x.
 *               Must be of at least length (1+(N-1)*abs(inc_x).
 * \param inc_x  how many places to skip to get to next element in x
 * \param y      A pointer to the beginning of the data in y.
 * \param inc_y  how many places to skip to get to next element in y
 *
 * @returns the dot product
 *
 * \ingroup QT
 */
PSI_API
double C_DDOT(size_t n, double* x, int inc_x, double* y, int inc_y);

/*!
 * This function returns the square of the norm of this vector.
 *
 * \param n      Number of elements in x.
 * \param x      A pointer to the beginning of the data in x.
 *               Must be of at least length (1+(N-1)*abs(inc_x).
 * \param inc_x  how many places to skip to get to next element in x
 *
 * @returns the norm squared product
 *
 * \ingroup QT
 */
PSI_API
double C_DNRM2(size_t n, double* x, int inc_x);

/*!
 * This function returns the sum of the absolute value of this vector.
 *
 * \param n      Number of elements in x.
 * \param x      A pointer to the beginning of the data in x.
 *               Must be of at least length (1+(N-1)*abs(inc_x).
 * \param inc_x  how many places to skip to get to next element in x
 *
 * @returns the sum of the absolute value
 *
 * \ingroup QT
 */
PSI_API
double C_DASUM(size_t n, double* x, int inc_x);

/*!
 * This function returns the index of the largest absolute value compoment of this vector.
 *
 * \param n      Number of elements in x.
 * \param x      A pointer to the beginning of the data in x.
 *               Must be of at least length (1+(N-1)*abs(inc_x).
 * \param inc_x  how many places to skip to get to next element in x
 *
 * @returns the index of the largest absolute value
 *
 * \ingroup QT
 */
PSI_API
size_t C_IDAMAX(size_t n, double* x, int inc_x);
}  // namespace psi
