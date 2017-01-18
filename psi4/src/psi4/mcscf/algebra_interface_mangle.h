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

#ifndef _psi_src_bin_mcscf_algebra_interface_mangle_h_
#define _psi_src_bin_mcscf_algebra_interface_mangle_h_

#ifndef FC_SYMBOL
#define FC_SYMBOL 2
#endif

#ifdef USE_FCMANGLE_H
#include "FCMangle.h"
#define F_DAXPY  FC_GLOBAL(daxpy,  DAXPY) 
#define F_DCOPY  FC_GLOBAL(dcopy,  DCOPY) 
#define F_DGEMM  FC_GLOBAL(dgemm,  DGEMM)
#define F_DROT   FC_GLOBAL(drot,   DROT)
#define F_DSCAL  FC_GLOBAL(dscal,  DSCAL)
#define F_DGEMV  FC_GLOBAL(dgemv,  DGEMV)
#define F_DSPMV  FC_GLOBAL(dfpmv,  DSPMV)
#define F_DDOT   FC_GLOBAL(ddot,   DDOT)
#define F_DGEEV  FC_GLOBAL(dgeev,  DGEEV)
#define F_DGESV  FC_GLOBAL(dgesv,  DGESV)
#define F_DGETRF FC_GLOBAL(dgetrf, DGETRF)
#define F_DGETRI FC_GLOBAL(dgetri, DGETRI)
#define F_DGESVD FC_GLOBAL(dgesvd, DGESVD)
#define F_DSYEV  FC_GLOBAL(dsyev,  DSYEV)
#else // USE_FCMANGLE_H
#if FC_SYMBOL==2
#define F_DAXPY daxpy_
#define F_DCOPY dcopy_
#define F_DGEMM dgemm_
#define F_DROT drot_
#define F_DSCAL dscal_
#define F_DGEMV dgemv_
#define F_DSPMV dspmv_
#define F_DDOT  ddot_
#define F_DGEEV dgeev_
#define F_DGESV dgesv_
#define F_DGETRF dgetrf_
#define F_DGETRI dgetri_
#define F_DGESVD dgesvd_
#define F_DSYEV dsyev_
#elif FC_SYMBOL==1
#define F_DAXPY daxpy
#define F_DCOPY dcopy
#define F_DGEMM dgemm
#define F_DROT drot
#define F_DSCAL dscal
#define F_DGEMV dgemv
#define F_DSPMV dspmv
#define F_DDOT  ddot
#define F_DGEEV dgeev
#define F_DGESV dgesv
#define F_DGETRF dgetrf
#define F_DGETRI dgetri
#define F_DGESVD dgesvd
#define F_DSYEV dsyev
#elif FC_SYMBOL==3
#define F_DAXPY DAXPY
#define F_DCOPY DCOPY
#define F_DGEMM DGEMM
#define F_DROT DROT
#define F_DSCAL DSCAL
#define F_DGEMV DGEMV
#define F_DSPMV DSPMV
#define F_DDOT  DDOT
#define F_DGEEV DGEEV
#define F_DGESV DGESV
#define F_DGETRF DGETRF
#define F_DGETRI DGETRI
#define F_DGESVD DGESVD
#define F_DSYEV DSYEV
#elif FC_SYMBOL==4
#define F_DAXPY DAXPY_
#define F_DCOPY DCOPY_
#define F_DGEMM DGEMM_
#define F_DROT DROT_
#define F_DSCAL DSCAL_
#define F_DGEMV DGEMV_
#define F_DSPMV DSPMV_
#define F_DDOT  DDOT_
#define F_DGEEV DGEEV_
#define F_DGESV DGESV_
#define F_DGETRF DGETRF_
#define F_DGETRI DGETRI_
#define F_DGESVD DGESVD_
#define F_DSYEV DSYEV_
#endif // FC_SYMBOL
#endif // USE_FCMANGLE_H

#endif // _psi_src_bin_mcscf_algebra_interface_mangle_h_