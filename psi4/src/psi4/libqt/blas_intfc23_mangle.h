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

#ifndef _psi_src_lib_libqt_blas_intfc23_mangle_h_
#define _psi_src_lib_libqt_blas_intfc23_mangle_h_

/*! \defgroup QT libqt: The Quantum-Trio Miscellaneous Library */

/*!
 \file
 \ingroup QT
 \brief The PSI3 BLAS2 and BLAS3 interface routines

 Declares mangling for BLAS2 and BLAS3 interface routines

*/

#ifdef USE_FCMANGLE_H
#include "FCMangle.h"
#define F_DGBMV FC_GLOBAL(dgbmv, DGBMV)
#define F_DGEMM FC_GLOBAL(dgemm, DGEMM)
#define F_DGEMV FC_GLOBAL(dgemv, DGEMV)
#define F_DGER FC_GLOBAL(dger, DGER)
#define F_DSBMV FC_GLOBAL(dsbmv, DSBMV)
#define F_DSPMV FC_GLOBAL(dspmv, DSPMV)
#define F_DSPR FC_GLOBAL(dspr, DSPR)
#define F_DSPR2 FC_GLOBAL(dspr2, DSPR2)
#define F_DSYMM FC_GLOBAL(dsymm, DSYMM)
#define F_DSYMV FC_GLOBAL(dsymv, DSYMV)
#define F_DSYR FC_GLOBAL(dsyr, DSYR)
#define F_DSYR2 FC_GLOBAL(dsyr2, DSYR2)
#define F_DSYR2K FC_GLOBAL(dsyr2k, DSYR2K)
#define F_DSYRK FC_GLOBAL(dsyrk, DSYRK)
#define F_DTBMV FC_GLOBAL(dtbmv, DTBMV)
#define F_DTBSV FC_GLOBAL(dtbsv, DTBSV)
#define F_DTPMV FC_GLOBAL(dtpmv, DTPMV)
#define F_DTPSV FC_GLOBAL(dtpsv, DTPSV)
#define F_DTRMM FC_GLOBAL(dtrmm, DTRMM)
#define F_DTRMV FC_GLOBAL(dtrmv, DTRMV)
#define F_DTRSM FC_GLOBAL(dtrsm, DTRSM)
#define F_DTRSV FC_GLOBAL(dtrsv, DTRSV)
#else  // USE_FCMANGLE_H
#if FC_SYMBOL == 2
#define F_DGBMV dgbmv_
#define F_DGEMM dgemm_
#define F_DGEMV dgemv_
#define F_DGER dger_
#define F_DSBMV dsbmv_
#define F_DSPMV dspmv_
#define F_DSPR dspr_
#define F_DSPR2 dspr2_
#define F_DSYMM dsymm_
#define F_DSYMV dsymv_
#define F_DSYR dsyr_
#define F_DSYR2 dsyr2_
#define F_DSYR2K dsyr2k_
#define F_DSYRK dsyrk_
#define F_DTBMV dtbmv_
#define F_DTBSV dtbsv_
#define F_DTPMV dtpmv_
#define F_DTPSV dtpsv_
#define F_DTRMM dtrmm_
#define F_DTRMV dtrmv_
#define F_DTRSM dtrsm_
#define F_DTRSV dtrsv_
#elif FC_SYMBOL == 1
#define F_DGBMV dgbmv
#define F_DGEMM dgemm
#define F_DGEMV dgemv
#define F_DGER dger
#define F_DSBMV dsbmv
#define F_DSPMV dspmv
#define F_DSPR dspr
#define F_DSPR2 dspr2
#define F_DSYMM dsymm
#define F_DSYMV dsymv
#define F_DSYR dsyr
#define F_DSYR2 dsyr2
#define F_DSYR2K dsyr2k
#define F_DSYRK dsyrk
#define F_DTBMV dtbmv
#define F_DTBSV dtbsv
#define F_DTPMV dtpmv
#define F_DTPSV dtpsv
#define F_DTRMM dtrmm
#define F_DTRMV dtrmv
#define F_DTRSM dtrsm
#define F_DTRSV dtrsv
#elif FC_SYMBOL == 3
#define F_DGBMV DGBMV
#define F_DGEMM DGEMM
#define F_DGEMV DGEMV
#define F_DGER DGER
#define F_DSBMV DSBMV
#define F_DSPMV DSPMV
#define F_DSPR DSPR
#define F_DSPR2 DSPR2
#define F_DSYMM DSYMM
#define F_DSYMV DSYMV
#define F_DSYR DSYR
#define F_DSYR2 DSYR2
#define F_DSYR2K DSYR2K
#define F_DSYRK DSYRK
#define F_DTBMV DTBMV
#define F_DTBSV DTBSV
#define F_DTPMV DTPMV
#define F_DTPSV DTPSV
#define F_DTRMM DTRMM
#define F_DTRMV DTRMV
#define F_DTRSM DTRSM
#define F_DTRSV DTRSV
#elif FC_SYMBOL == 4
#define F_DGBMV DGBMV_
#define F_DGEMM DGEMM_
#define F_DGEMV DGEMV_
#define F_DGER DGER_
#define F_DSBMV DSBMV_
#define F_DSPMV DSPMV_
#define F_DSPR DSPR_
#define F_DSPR2 DSPR2_
#define F_DSYMM DSYMM_
#define F_DSYMV DSYMV_
#define F_DSYR DSYR_
#define F_DSYR2 DSYR2_
#define F_DSYR2K DSYR2K_
#define F_DSYRK DSYRK_
#define F_DTBMV DTBMV_
#define F_DTBSV DTBSV_
#define F_DTPMV DTPMV_
#define F_DTPSV DTPSV_
#define F_DTRMM DTRMM_
#define F_DTRMV DTRMV_
#define F_DTRSM DTRSM_
#define F_DTRSV DTRSV_
#endif
#endif

#endif
