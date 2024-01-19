/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#ifndef _psi_src_lib_libqt_blas_intfc_mangle_h_
#define _psi_src_lib_libqt_blas_intfc_mangle_h_

/*! \defgroup QT libqt: The Quantum-Trio Miscellaneous Library */

/*!
 \file
 \ingroup QT
 \brief The PSI3 BLAS1 interface routines

 Declares mangling for BLAS1 interface routines

*/

#ifdef USE_FCMANGLE_H
#include "FCMangle.h"
#define F_DSWAP FC_GLOBAL(dswap, DSWAP)
#define F_DAXPY FC_GLOBAL(daxpy, DAXPY)
#define F_DAXPBY FC_GLOBAL(daxpby, DAXPBY)
#define F_DCOPY FC_GLOBAL(dcopy, DCOPY)
#define F_DROT FC_GLOBAL(drot, DROT)
#define F_DSCAL FC_GLOBAL(dscal, DSCAL)
#define F_DDOT FC_GLOBAL(ddot, DDOT)
#define F_DASUM FC_GLOBAL(dasum, DASUM)
#define F_DNRM2 FC_GLOBAL(dnrm2, DNRM2)
#define F_IDAMAX FC_GLOBAL(idamax, IDAMAX)
#else  // USE_FCMANGLE_H
#if FC_SYMBOL == 2
#define F_DSWAP dswap_
#define F_DAXPY daxpy_
#define F_DAXPBY daxpby_
#define F_DCOPY dcopy_
#define F_DROT drot_
#define F_DSCAL dscal_
#define F_DDOT ddot_
#define F_DASUM dasum_
#define F_DNRM2 dnrm2_
#define F_IDAMAX idamax_
#elif FC_SYMBOL == 1
#define F_DSWAP dswap
#define F_DAXPY daxpy
#define F_DAXPBY daxpby
#define F_DCOPY dcopy
#define F_DROT drot
#define F_DSCAL dscal
#define F_DDOT ddot
#define F_DASUM dasum
#define F_DNRM2 dnrm2
#define F_IDAMAX idamax
#elif FC_SYMBOL == 3
#define F_DSWAP DSWAP
#define F_DAXPY DAXPY
#define F_DAXPBY DAXPBY
#define F_DCOPY DCOPY
#define F_DROT DROT
#define F_DSCAL DSCAL
#define F_DDOT DDOT
#define F_DASUM DASUM
#define F_DNRM2 DNRM2
#define F_IDAMAX IDAMAX
#elif FC_SYMBOL == 4
#define F_DSWAP DSWAP_
#define F_DAXPY DAXPY_
#define F_DAXPBY DAXPBY_
#define F_DCOPY DCOPY_
#define F_DROT DROT_
#define F_DSCAL DSCAL_
#define F_DDOT DDOT_
#define F_DASUM DASUM_
#define F_DNRM2 DNRM2_
#define F_IDAMAX IDAMAX_
#endif
#endif

#endif
