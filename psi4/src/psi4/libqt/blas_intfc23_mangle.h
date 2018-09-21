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
#define F_DGER  FC_GLOBAL(dger,  DGER)
#define F_DSBMV FC_GLOBAL(dsbmv, DSBMV)
#define F_DSPMV FC_GLOBAL(dspmv, DSPMV)
#define F_DSPR  FC_GLOBAL(dspr,  DSPR)
#define F_DSPR2 FC_GLOBAL(dspr2, DSPR2)
#define F_DSYMM FC_GLOBAL(dsymm, DSYMM)
#define F_DSYMV FC_GLOBAL(dsymv, DSYMV)
#define F_DSYR  FC_GLOBAL(dsyr,  DSYR)
#define F_DSYR2 FC_GLOBAL(dsyr2, DSYR2)
#define F_DSYR2K FC_GLOBAL(dsyr2k DSYR2K)
#define F_DSYRK FC_GLOBAL(dsyrk, DSYRK)
#define F_DTBMV FC_GLOBAL(dtbmv, DTBMV)
#define F_DTBSV FC_GLOBAL(dtbsv, DTBSV)
#define F_DTPMV FC_GLOBAL(dtpmv, DTPMV)
#define F_DTPSV FC_GLOBAL(dtpsv, DTPSV)
#define F_DTRMM FC_GLOBAL(dtrmm, DTRMM)
#define F_DTRMV FC_GLOBAL(dtrmv, DTRMV)
#define F_DTRSM FC_GLOBAL(dtrsm, DTRSM)
#define F_DTRSV FC_GLOBAL(dtrsv, DTRSV)
#define F_SGBMV FC_GLOBAL(sgbmv, SGBMV)
#define F_SGEMM FC_GLOBAL(sgemm, SGEMM)
#define F_SGEMV FC_GLOBAL(sgemv, SGEMV)
#define F_SGER  FC_GLOBAL(sger,  SGER)
#define F_SSBMV FC_GLOBAL(ssbmv, SSBMV)
#define F_SSPMV FC_GLOBAL(sspmv, SSPMV)
#define F_SSPR  FC_GLOBAL(sspr,  SSPR)
#define F_SSPR2 FC_GLOBAL(sspr2, SSPR2)
#define F_SSYMM FC_GLOBAL(ssymm, SSYMM)
#define F_SSYMV FC_GLOBAL(ssymv, SSYMV)
#define F_SSYR  FC_GLOBAL(ssyr,  SSYR)
#define F_SSYR2 FC_GLOBAL(ssyr2, SSYR2)
#define F_SSYR2K FC_GLOBAL(ssyr2k SSYR2K)
#define F_SSYRK FC_GLOBAL(ssyrk, SSYRK)
#define F_STBMV FC_GLOBAL(stbmv, STBMV)
#define F_STBSV FC_GLOBAL(stbsv, STBSV)
#define F_STPMV FC_GLOBAL(stpmv, STPMV)
#define F_STPSV FC_GLOBAL(stpsv, STPSV)
#define F_STRMM FC_GLOBAL(strmm, STRMM)
#define F_STRMV FC_GLOBAL(strmv, STRMV)
#define F_STRSM FC_GLOBAL(strsm, STRSM)
#define F_STRSV FC_GLOBAL(strsv, STRSV)
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
#define F_SGBMV sgbmv_
#define F_SGEMM sgemm_
#define F_SGEMV sgemv_
#define F_SGER  sger_
#define F_SSBMV ssbmv_
#define F_SSPMV sspmv_
#define F_SSPR  sspr_
#define F_SSPR2 sspr2_
#define F_SSYMM ssymm_
#define F_SSYMV ssymv_
#define F_SSYR  ssyr_
#define F_SSYR2 ssyr2_
#define F_SSYR2K ssyr2k_
#define F_SSYRK ssyrk_
#define F_STBMV stbmv_
#define F_STBSV stbsv_
#define F_STPMV stpmv_
#define F_STPSV stpsv_
#define F_STRMM strmm_
#define F_STRMV strmv_
#define F_STRSM strsm_
#define F_STRSV strsv_
#elif FC_SYMBOL == 1
#define F_DGBMV dgbmv
#define F_DGEMM dgemm
#define F_DGEMV dgemv
#define F_DGER  dger
#define F_DSBMV dsbmv
#define F_DSPMV dspmv
#define F_DSPR  dspr
#define F_DSPR2 dspr2
#define F_DSYMM dsymm
#define F_DSYMV dsymv
#define F_DSYR  dsyr
#define F_DSYR2 dsyr2
#define F_DSYR2Kdsyr2k
#define F_DSYRK dsyrk
#define F_DTBMV dtbmv
#define F_DTBSV dtbsv
#define F_DTPMV dtpmv
#define F_DTPSV dtpsv
#define F_DTRMM dtrmm
#define F_DTRMV dtrmv
#define F_DTRSM dtrsm
#define F_DTRSV dtrsv
#define F_SGBMV  sgbmv
#define F_SGEMM  sgemm
#define F_SGEMV  sgemv
#define F_SGER   sger
#define F_SSBMV  ssbmv
#define F_SSPMV  sspmv
#define F_SSPR   sspr
#define F_SSPR2  sspr2
#define F_SSYMM  ssymm
#define F_SSYMV  ssymv
#define F_SSYR   ssyr
#define F_SSYR2  ssyr2
#define F_SSYR2K ssyr2k
#define F_SSYRK  ssyrk
#define F_STBMV  stbmv
#define F_STBSV  stbsv
#define F_STPMV  stpmv
#define F_STPSV  stpsv
#define F_STRMM  strmm
#define F_STRMV  strmv
#define F_STRSM  strsm
#define F_STRSV  strsv
#elif FC_SYMBOL == 3
#define F_DGBMV DGBMV
#define F_DGEMM DGEMM
#define F_DGEMV DGEMV
#define F_DGER  DGER
#define F_DSBMV DSBMV
#define F_DSPMV DSPMV
#define F_DSPR  DSPR
#define F_DSPR2 DSPR2
#define F_DSYMM DSYMM
#define F_DSYMV DSYMV
#define F_DSYR  DSYR
#define F_DSYR2 DSYR2
#define F_DSYR2KDSYR2K
#define F_DSYRK DSYRK
#define F_DTBMV DTBMV
#define F_DTBSV DTBSV
#define F_DTPMV DTPMV
#define F_DTPSV DTPSV
#define F_DTRMM DTRMM
#define F_DTRMV DTRMV
#define F_DTRSM DTRSM
#define F_DTRSV DTRSV
#define F_SGBMV SGBMV
#define F_SGEMM SGEMM
#define F_SGEMV SGEMV
#define F_SGER  SGER
#define F_SSBMV SSBMV
#define F_SSPMV SSPMV
#define F_SSPR  SSPR
#define F_SSPR2 SSPR2
#define F_SSYMM SSYMM
#define F_SSYMV SSYMV
#define F_SSYR  SSYR
#define F_SSYR2 SSYR2
#define F_SSYR2K SSYR2K
#define F_SSYRK SSYRK
#define F_STBMV STBMV
#define F_STBSV STBSV
#define F_STPMV STPMV
#define F_STPSV STPSV
#define F_STRMM STRMM
#define F_STRMV STRMV
#define F_STRSM STRSM
#define F_STRSV STRSV
#elif FC_SYMBOL == 4
#define F_DGBMV DGBMV_
#define F_DGEMM DGEMM_
#define F_DGEMV DGEMV_
#define F_DGER  DGER_
#define F_DSBMV DSBMV_
#define F_DSPMV DSPMV_
#define F_DSPR  DSPR_
#define F_DSPR2 DSPR2_
#define F_DSYMM DSYMM_
#define F_DSYMV DSYMV_
#define F_DSYR  DSYR_
#define F_DSYR2 DSYR2_
#define F_DSYR2KDSYR2K_
#define F_DSYRK DSYRK_
#define F_DTBMV DTBMV_
#define F_DTBSV DTBSV_
#define F_DTPMV DTPMV_
#define F_DTPSV DTPSV_
#define F_DTRMM DTRMM_
#define F_DTRMV DTRMV_
#define F_DTRSM DTRSM_
#define F_DTRSV DTRSV_
#define F_SGBMV  SGBMV_
#define F_SGEMM  SGEMM_
#define F_SGEMV  SGEMV_
#define F_SGER   SGER_
#define F_SSBMV  SSBMV_
#define F_SSPMV  SSPMV_
#define F_SSPR   SSPR_
#define F_SSPR2  SSPR2_
#define F_SSYMM  SSYMM_
#define F_SSYMV  SSYMV_
#define F_SSYR   SSYR_
#define F_SSYR2  SSYR2_
#define F_SSYR2K SSYR2K_
#define F_SSYRK  SSYRK_
#define F_STBMV  STBMV_
#define F_STBSV  STBSV_
#define F_STPMV  STPMV_
#define F_STPSV  STPSV_
#define F_STRMM  STRMM_
#define F_STRMV  STRMV_
#define F_STRSM  STRSM_
#define F_STRSV  STRSV_
#endif
#endif

#endif
