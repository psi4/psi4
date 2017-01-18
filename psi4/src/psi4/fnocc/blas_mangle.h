/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef BLAS_MANGLE_H
#define BLAS_MANGLE_H

/**
 * Declare mangling for fortran-ordered blas routines
 */

#ifndef FC_SYMBOL
#define FC_SYMBOL 2
#endif

#ifdef USE_FCMANGLE_H
#include "FCMangle.h"
#define dgemv  FC_GLOBAL(dgemv , DGEMV)
#define dgemm  FC_GLOBAL(dgemm , DGEMM)
#define dcopy  FC_GLOBAL(dcopy , DCOPY)
#define daxpy  FC_GLOBAL(daxpy , DAXPY)
#define dnrm2  FC_GLOBAL(drnm2 , DRNM2)
#define dgesv  FC_GLOBAL(dgesv , DGESV)
#define ddot   FC_GLOBAL(ddot  , DDOT)
#define dsyev  FC_GLOBAL(dsyev , DSYEV)
#define dspev  FC_GLOBAL(dspev , DSPEV)
#define dgesvd FC_GLOBAL(dgesvd, DGESVD)
#else // USE_FCMANGLE_H
#if FC_SYMBOL==2
#define dgemv  dgemv_
#define dgemm  dgemm_
#define dcopy  dcopy_
#define daxpy  daxpy_
#define dnrm2  drnm2_
#define dgesv  dgesv_
#define ddot   ddot_  
#define dsyev  dsyev_ 
#define dspev  dspev_
#define dgesvd dgesvd_
#elif FC_SYMBOL==1
#define dgemv  dgemv
#define dgemm  dgemm
#define dcopy  dcopy
#define daxpy  daxpy
#define dnrm2  drnm2
#define dgesv  dgesv
#define ddot   ddot 
#define dsyev  dsyev
#define dspev  dspev
#define dgesvd dgesvd
#elif FC_SYMBOL==3
#define dgemv  DGEMV
#define dgemm  DGEMM
#define dcopy  DCOPY
#define daxpy  DAXPY
#define dnrm2  DRNM2
#define dgesv  DGESV
#define ddot   DDOT
#define dsyev  DSYEV
#define dspev  DSPEV
#define dgesvd DGESVD
#elif FC_SYMBOL==4
#define dgemv  DGEMV_
#define dgemm  DGEMM_
#define dcopy  DCOPY_
#define daxpy  DAXPY_
#define dnrm2  DRNM2_
#define dgesv  DGESV_
#define ddot   DDOT_
#define dsyev  DSYEV_
#define dspev  DSPEV_
#define dgesvd DGESVD_
#endif // FC_SYMBOL
#endif // USE_FCMANGLE_H

#endif