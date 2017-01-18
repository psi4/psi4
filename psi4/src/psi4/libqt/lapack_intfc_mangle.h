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

#ifndef _psi_src_lib_libqt_lapack_intfc_mangle_h_
#define _psi_src_lib_libqt_lapack_intfc_mangle_h_

/*! \defgroup QT libqt: The Quantum-Trio Miscellaneous Library */

/*!
 \file
 \ingroup QT
 \brief The PSI3 LAPACK interface routines

 Declares mangling for LAPACK interface routines

*/

#ifdef USE_FCMANGLE_H
#include "FCMangle.h"
#define F_DGEEV  FC_GLOBAL(dgeev,  DGEEV )
#define F_DGESV  FC_GLOBAL(dgesv,  DGESV )
#define F_DGETRF FC_GLOBAL(dgetrf, DGETRF)
#define F_DGETRI FC_GLOBAL(dgetri, DGETRI)
#define F_DPOTRF FC_GLOBAL(dpotrf, DPOTRF)
#define F_DPOTRI FC_GLOBAL(dpotri, DPOTRI)
#define F_DPOTRS FC_GLOBAL(dpotrs, DPOTRS)
#define F_DGESVD FC_GLOBAL(dgesvd, DGESVD)
#define F_DSYEV  FC_GLOBAL(dsyev,  DSYEV )
#define F_DBDSDC FC_GLOBAL(dbdsdc, DBDSDC)
#define F_DBDSQR FC_GLOBAL(dbdsqr, DBDSQR)
#define F_DDISNA FC_GLOBAL(ddisna, DDISNA)
#define F_DGBBRD FC_GLOBAL(dgbbrd, DGBBRD)
#define F_DGBCON FC_GLOBAL(dgbcon, DGBCON)
#define F_DGBEQU FC_GLOBAL(dgbequ, DGBEQU)
#define F_DGBRFS FC_GLOBAL(dgbrfs, DGBRFS)
#define F_DGBSV  FC_GLOBAL(dgbsv,  DGBSV )
#define F_DGBSVX FC_GLOBAL(dgbsvx, DGBSVX)
#define F_DGBTRF FC_GLOBAL(dgbtrf, DGBTRF)
#define F_DGBTRS FC_GLOBAL(dgbtrs, DGBTRS)
#define F_DGEBAK FC_GLOBAL(dgebak, DGEBAK)
#define F_DGEBAL FC_GLOBAL(dgebal, DGEBAL)
#define F_DGEBRD FC_GLOBAL(dgebrd, DGEBRD)
#define F_DGECON FC_GLOBAL(dgecon, DGECON)
#define F_DGEEQU FC_GLOBAL(dgeequ, DGEEQU)
#define F_DGEES  FC_GLOBAL(dgees,  DGEES )
#define F_DGEESX FC_GLOBAL(dgeesx, DGEESX)
#define F_DGEEV  FC_GLOBAL(dgeev,  DGEEV )
#define F_DGEEVX FC_GLOBAL(dgeevx, DGEEVX)
#define F_DGEGS  FC_GLOBAL(dgegs,  DGEGS )
#define F_DGEGV  FC_GLOBAL(dgegv,  DGEGV )
#define F_DGEHRD FC_GLOBAL(dgehrd, DGEHRD)
#define F_DGELQF FC_GLOBAL(dgelqf, DGELQF)
#define F_DGELS  FC_GLOBAL(dgels,  DGELS )
#define F_DGELSD FC_GLOBAL(dgelsd, DGELSD)
#define F_DGELSS FC_GLOBAL(dgelss, DGELSS)
#define F_DGELSX FC_GLOBAL(dgelsx, DGELSX)
#define F_DGELSY FC_GLOBAL(dgelsy, DGELSY)
#define F_DGEQLF FC_GLOBAL(dgeqlf, DGEQLF)
#define F_DGEQP3 FC_GLOBAL(dgeqp3, DGEQP3)
#define F_DGEQPF FC_GLOBAL(dgeqpf, DGEQPF)
#define F_DGEQRF FC_GLOBAL(dgeqrf, DGEQRF)
#define F_DGERFS FC_GLOBAL(dgerfs, DGERFS)
#define F_DGERQF FC_GLOBAL(dgerqf, DGERQF)
#define F_DGESDD FC_GLOBAL(dgesdd, DGESDD)
#define F_DGESV  FC_GLOBAL(dgesv,  DGESV )
#define F_DGESVX FC_GLOBAL(dgesvx, DGESVX)
#define F_DGETRF FC_GLOBAL(dgetrf, DGETRF)
#define F_DGETRI FC_GLOBAL(dgetri, DGETRI)
#define F_DGETRS FC_GLOBAL(dgetrs, DGETRS)
#define F_DGGBAK FC_GLOBAL(dggbak, DGGBAK)
#define F_DGGBAL FC_GLOBAL(dggbal, DGGBAL)
#define F_DGGES  FC_GLOBAL(dgges,  DGGES )
#define F_DGGESX FC_GLOBAL(dggesx, DGGESX)
#define F_DGGEV  FC_GLOBAL(dggev,  DGGEV )
#define F_DGGEVX FC_GLOBAL(dggevx, DGGEVX)
#define F_DGGGLM FC_GLOBAL(dggglm, DGGGLM)
#define F_DGGHRD FC_GLOBAL(dgghrd, DGGHRD)
#define F_DGGLSE FC_GLOBAL(dgglse, DGGLSE)
#define F_DGGQRF FC_GLOBAL(dggqrf, DGGQRF)
#define F_DGGRQF FC_GLOBAL(dggrqf, DGGRQF)
#define F_DGGSVD FC_GLOBAL(dggsvd, DGGSVD)
#define F_DGGSVP FC_GLOBAL(dggsvp, DGGSVP)
#define F_DGTCON FC_GLOBAL(dgtcon, DGTCON)
#define F_DGTRFS FC_GLOBAL(dgtrfs, DGTRFS)
#define F_DGTSV  FC_GLOBAL(dgtsv,  DGTSV )
#define F_DGTSVX FC_GLOBAL(dgtsvx, DGTSVX)
#define F_DGTTRF FC_GLOBAL(dgttrf, DGTTRF)
#define F_DGTTRS FC_GLOBAL(dgttrs, DGTTRS)
#define F_DHGEQZ FC_GLOBAL(dhgeqz, DHGEQZ)
#define F_DHSEIN FC_GLOBAL(dhsein, DHSEIN)
#define F_DHSEQR FC_GLOBAL(dhseqr, DHSEQR)
#define F_DOPGTR FC_GLOBAL(dopgtr, DOPGTR)
#define F_DOPMTR FC_GLOBAL(dopmtr, DOPMTR)
#define F_DORGBR FC_GLOBAL(dorgbr, DORGBR)
#define F_DORGHR FC_GLOBAL(dorghr, DORGHR)
#define F_DORGLQ FC_GLOBAL(dorglq, DORGLQ)
#define F_DORGQL FC_GLOBAL(dorgql, DORGQL)
#define F_DORGQR FC_GLOBAL(dorgqr, DORGQR)
#define F_DORGRQ FC_GLOBAL(dorgrq, DORGRQ)
#define F_DORGTR FC_GLOBAL(dorgtr, DORGTR)
#define F_DORMBR FC_GLOBAL(dormbr, DORMBR)
#define F_DORMHR FC_GLOBAL(dormhr, DORMHR)
#define F_DORMLQ FC_GLOBAL(dormlq, DORMLQ)
#define F_DORMQL FC_GLOBAL(dormql, DORMQL)
#define F_DORMQR FC_GLOBAL(dormqr, DORMQR)
#define F_DORMR3 FC_GLOBAL(dormr3, DORMR3)
#define F_DORMRQ FC_GLOBAL(dormrq, DORMRQ)
#define F_DORMRZ FC_GLOBAL(dormrz, DORMRZ)
#define F_DORMTR FC_GLOBAL(dormtr, DORMTR)
#define F_DPBCON FC_GLOBAL(dpbcon, DPBCON)
#define F_DPBEQU FC_GLOBAL(dpbequ, DPBEQU)
#define F_DPBRFS FC_GLOBAL(dpbrfs, DPBRFS)
#define F_DPBSTF FC_GLOBAL(dpbstf, DPBSTF)
#define F_DPBSV  FC_GLOBAL(dpbsv,  DPBSV )
#define F_DPBSVX FC_GLOBAL(dpbsvx, DPBSVX)
#define F_DPBTRF FC_GLOBAL(dpbtrf, DPBTRF)
#define F_DPBTRS FC_GLOBAL(dpbtrs, DPBTRS)
#define F_DPOCON FC_GLOBAL(dpocon, DPOCON)
#define F_DPOEQU FC_GLOBAL(dpoequ, DPOEQU)
#define F_DPORFS FC_GLOBAL(dporfs, DPORFS)
#define F_DPOSV  FC_GLOBAL(dposv,  DPOSV )
#define F_DPOSVX FC_GLOBAL(dposvx, DPOSVX)
#define F_DPOTRF FC_GLOBAL(dpotrf, DPOTRF)
#define F_DPOTRI FC_GLOBAL(dpotri, DPOTRI)
#define F_DPOTRS FC_GLOBAL(dpotrs, DPOTRS)
#define F_DPPCON FC_GLOBAL(dppcon, DPPCON)
#define F_DPPEQU FC_GLOBAL(dppequ, DPPEQU)
#define F_DPPRFS FC_GLOBAL(dpprfs, DPPRFS)
#define F_DPPSV  FC_GLOBAL(dppsv,  DPPSV )
#define F_DPPSVX FC_GLOBAL(dppsvx, DPPSVX)
#define F_DPPTRF FC_GLOBAL(dpptrf, DPPTRF)
#define F_DPPTRI FC_GLOBAL(dpptri, DPPTRI)
#define F_DPPTRS FC_GLOBAL(dpptrs, DPPTRS)
#define F_DPTCON FC_GLOBAL(dptcon, DPTCON)
#define F_DPTEQR FC_GLOBAL(dpteqr, DPTEQR)
#define F_DPTRFS FC_GLOBAL(dptrfs, DPTRFS)
#define F_DPTSV  FC_GLOBAL(dptsv,  DPTSV )
#define F_DPTSVX FC_GLOBAL(dptsvx, DPTSVX)
#define F_DPTTRF FC_GLOBAL(dpttrf, DPTTRF)
#define F_DPTTRS FC_GLOBAL(dpttrs, DPTTRS)
#define F_DSBEV  FC_GLOBAL(dsbev,  DSBEV )
#define F_DSBEVD FC_GLOBAL(dsbevd, DSBEVD)
#define F_DSBEVX FC_GLOBAL(dsbevx, DSBEVX)
#define F_DSBGST FC_GLOBAL(dsbgst, DSBGST)
#define F_DSBGV  FC_GLOBAL(dsbgv,  DSBGV )
#define F_DSBGVD FC_GLOBAL(dsbgvd, DSBGVD)
#define F_DSBGVX FC_GLOBAL(dsbgvx, DSBGVX)
#define F_DSBTRD FC_GLOBAL(dsbtrd, DSBTRD)
#define F_DSGESV FC_GLOBAL(dsgesv, DSGESV)
#define F_DSPCON FC_GLOBAL(dspcon, DSPCON)
#define F_DSPEV  FC_GLOBAL(dspev,  DSPEV )
#define F_DSPEVD FC_GLOBAL(dspevd, DSPEVD)
#define F_DSPEVX FC_GLOBAL(dspevx, DSPEVX)
#define F_DSPGST FC_GLOBAL(dspgst, DSPGST)
#define F_DSPGV  FC_GLOBAL(dspgv,  DSPGV )
#define F_DSPGVD FC_GLOBAL(dspgvd, DSPGVD)
#define F_DSPGVX FC_GLOBAL(dspgvx, DSPGVX)
#define F_DSPRFS FC_GLOBAL(dsprfs, DSPRFS)
#define F_DSPSV  FC_GLOBAL(dspsv,  DSPSV )
#define F_DSPSVX FC_GLOBAL(dspsvx, DSPSVX)
#define F_DSPTRD FC_GLOBAL(dsptrd, DSPTRD)
#define F_DSPTRF FC_GLOBAL(dsptrf, DSPTRF)
#define F_DSPTRI FC_GLOBAL(dsptri, DSPTRI)
#define F_DSPTRS FC_GLOBAL(dsptrs, DSPTRS)
#define F_DSTEBZ FC_GLOBAL(dstebz, DSTEBZ)
#define F_DSTEDC FC_GLOBAL(dstedc, DSTEDC)
#define F_DSTEGR FC_GLOBAL(dstegr, DSTEGR)
#define F_DSTEIN FC_GLOBAL(dstein, DSTEIN)
#define F_DSTEQR FC_GLOBAL(dsteqr, DSTEQR)
#define F_DSTERF FC_GLOBAL(dsterf, DSTERF)
#define F_DSTEV  FC_GLOBAL(dstev,  DSTEV )
#define F_DSTEVD FC_GLOBAL(dstevd, DSTEVD)
#define F_DSTEVR FC_GLOBAL(dstevr, DSTEVR)
#define F_DSTEVX FC_GLOBAL(dstevx, DSTEVX)
#define F_DSYCON FC_GLOBAL(dsycon, DSYCON)
#define F_DSYEV  FC_GLOBAL(dsyev,  DSYEV )
#define F_DSYEVD FC_GLOBAL(dsyevd, DSYEVD)
#define F_DSYEVR FC_GLOBAL(dsyevr, DSYEVR)
#define F_DSYEVX FC_GLOBAL(dsyevx, DSYEVX)
#define F_DSYGST FC_GLOBAL(dsygst, DSYGST)
#define F_DSYGV  FC_GLOBAL(dsygv,  DSYGV )
#define F_DSYGVD FC_GLOBAL(dsygvd, DSYGVD)
#define F_DSYGVX FC_GLOBAL(dsygvx, DSYGVX)
#define F_DSYRFS FC_GLOBAL(dsyrfs, DSYRFS)
#define F_DSYSV  FC_GLOBAL(dsysv,  DSYSV )
#define F_DSYSVX FC_GLOBAL(dsysvx, DSYSVX)
#define F_DSYTRD FC_GLOBAL(dsytrd, DSYTRD)
#define F_DSYTRF FC_GLOBAL(dsytrf, DSYTRF)
#define F_DSYTRI FC_GLOBAL(dsytri, DSYTRI)
#define F_DSYTRS FC_GLOBAL(dsytrs, DSYTRS)
#define F_DTBCON FC_GLOBAL(dtbcon, DTBCON)
#define F_DTBRFS FC_GLOBAL(dtbrfs, DTBRFS)
#define F_DTBTRS FC_GLOBAL(dtbtrs, DTBTRS)
#define F_DTGEVC FC_GLOBAL(dtgevc, DTGEVC)
#define F_DTGEXC FC_GLOBAL(dtgexc, DTGEXC)
#define F_DTGSEN FC_GLOBAL(dtgsen, DTGSEN)
#define F_DTGSJA FC_GLOBAL(dtgsja, DTGSJA)
#define F_DTGSNA FC_GLOBAL(dtgsna, DTGSNA)
#define F_DTGSYL FC_GLOBAL(dtgsyl, DTGSYL)
#define F_DTPCON FC_GLOBAL(dtpcon, DTPCON)
#define F_DTPRFS FC_GLOBAL(dtprfs, DTPRFS)
#define F_DTPTRI FC_GLOBAL(dtptri, DTPTRI)
#define F_DTPTRS FC_GLOBAL(dtptrs, DTPTRS)
#define F_DTRCON FC_GLOBAL(dtrcon, DTRCON)
#define F_DTREVC FC_GLOBAL(dtrevc, DTREVC)
#define F_DTREXC FC_GLOBAL(dtrexc, DTREXC)
#define F_DTRRFS FC_GLOBAL(dtrrfs, DTRRFS)
#define F_DTRSEN FC_GLOBAL(dtrsen, DTRSEN)
#define F_DTRSNA FC_GLOBAL(dtrsna, DTRSNA)
#define F_DTRSYL FC_GLOBAL(dtrsyl, DTRSYL)
#define F_DTRTRI FC_GLOBAL(dtrtri, DTRTRI)
#define F_DTRTRS FC_GLOBAL(dtrtrs, DTRTRS)
#define F_DTZRQF FC_GLOBAL(dtzrqf, DTZRQF)
#define F_DTZRZF FC_GLOBAL(dtzrzf, DTZRZF)
#else // USE_FCMANGLE_H
#if FC_SYMBOL==2
#define F_DGEEV dgeev_
#define F_DGESV dgesv_
#define F_DGETRF dgetrf_
#define F_DGETRI dgetri_
#define F_DPOTRF dpotrf_
#define F_DPOTRI dpotri_
#define F_DPOTRS dpotrs_
#define F_DGESVD dgesvd_
#define F_DSYEV dsyev_
#define F_DBDSDC dbdsdc_
#define F_DBDSQR dbdsqr_
#define F_DDISNA ddisna_
#define F_DGBBRD dgbbrd_
#define F_DGBCON dgbcon_
#define F_DGBEQU dgbequ_
#define F_DGBRFS dgbrfs_
#define F_DGBSV dgbsv_
#define F_DGBSVX dgbsvx_
#define F_DGBTRF dgbtrf_
#define F_DGBTRS dgbtrs_
#define F_DGEBAK dgebak_
#define F_DGEBAL dgebal_
#define F_DGEBRD dgebrd_
#define F_DGECON dgecon_
#define F_DGEEQU dgeequ_
#define F_DGEES dgees_
#define F_DGEESX dgeesx_
#define F_DGEEV dgeev_
#define F_DGEEVX dgeevx_
#define F_DGEGS dgegs_
#define F_DGEGV dgegv_
#define F_DGEHRD dgehrd_
#define F_DGELQF dgelqf_
#define F_DGELS dgels_
#define F_DGELSD dgelsd_
#define F_DGELSS dgelss_
#define F_DGELSX dgelsx_
#define F_DGELSY dgelsy_
#define F_DGEQLF dgeqlf_
#define F_DGEQP3 dgeqp3_
#define F_DGEQPF dgeqpf_
#define F_DGEQRF dgeqrf_
#define F_DGERFS dgerfs_
#define F_DGERQF dgerqf_
#define F_DGESDD dgesdd_
#define F_DGESV dgesv_
#define F_DGESVX dgesvx_
#define F_DGETRF dgetrf_
#define F_DGETRI dgetri_
#define F_DGETRS dgetrs_
#define F_DGGBAK dggbak_
#define F_DGGBAL dggbal_
#define F_DGGES dgges_
#define F_DGGESX dggesx_
#define F_DGGEV dggev_
#define F_DGGEVX dggevx_
#define F_DGGGLM dggglm_
#define F_DGGHRD dgghrd_
#define F_DGGLSE dgglse_
#define F_DGGQRF dggqrf_
#define F_DGGRQF dggrqf_
#define F_DGGSVD dggsvd_
#define F_DGGSVP dggsvp_
#define F_DGTCON dgtcon_
#define F_DGTRFS dgtrfs_
#define F_DGTSV dgtsv_
#define F_DGTSVX dgtsvx_
#define F_DGTTRF dgttrf_
#define F_DGTTRS dgttrs_
#define F_DHGEQZ dhgeqz_
#define F_DHSEIN dhsein_
#define F_DHSEQR dhseqr_
#define F_DOPGTR dopgtr_
#define F_DOPMTR dopmtr_
#define F_DORGBR dorgbr_
#define F_DORGHR dorghr_
#define F_DORGLQ dorglq_
#define F_DORGQL dorgql_
#define F_DORGQR dorgqr_
#define F_DORGRQ dorgrq_
#define F_DORGTR dorgtr_
#define F_DORMBR dormbr_
#define F_DORMHR dormhr_
#define F_DORMLQ dormlq_
#define F_DORMQL dormql_
#define F_DORMQR dormqr_
#define F_DORMR3 dormr3_
#define F_DORMRQ dormrq_
#define F_DORMRZ dormrz_
#define F_DORMTR dormtr_
#define F_DPBCON dpbcon_
#define F_DPBEQU dpbequ_
#define F_DPBRFS dpbrfs_
#define F_DPBSTF dpbstf_
#define F_DPBSV dpbsv_
#define F_DPBSVX dpbsvx_
#define F_DPBTRF dpbtrf_
#define F_DPBTRS dpbtrs_
#define F_DPOCON dpocon_
#define F_DPOEQU dpoequ_
#define F_DPORFS dporfs_
#define F_DPOSV dposv_
#define F_DPOSVX dposvx_
#define F_DPOTRF dpotrf_
#define F_DPOTRI dpotri_
#define F_DPOTRS dpotrs_
#define F_DPPCON dppcon_
#define F_DPPEQU dppequ_
#define F_DPPRFS dpprfs_
#define F_DPPSV dppsv_
#define F_DPPSVX dppsvx_
#define F_DPPTRF dpptrf_
#define F_DPPTRI dpptri_
#define F_DPPTRS dpptrs_
#define F_DPTCON dptcon_
#define F_DPTEQR dpteqr_
#define F_DPTRFS dptrfs_
#define F_DPTSV dptsv_
#define F_DPTSVX dptsvx_
#define F_DPTTRF dpttrf_
#define F_DPTTRS dpttrs_
#define F_DSBEV dsbev_
#define F_DSBEVD dsbevd_
#define F_DSBEVX dsbevx_
#define F_DSBGST dsbgst_
#define F_DSBGV dsbgv_
#define F_DSBGVD dsbgvd_
#define F_DSBGVX dsbgvx_
#define F_DSBTRD dsbtrd_
#define F_DSGESV dsgesv_
#define F_DSPCON dspcon_
#define F_DSPEV dspev_
#define F_DSPEVD dspevd_
#define F_DSPEVX dspevx_
#define F_DSPGST dspgst_
#define F_DSPGV dspgv_
#define F_DSPGVD dspgvd_
#define F_DSPGVX dspgvx_
#define F_DSPRFS dsprfs_
#define F_DSPSV dspsv_
#define F_DSPSVX dspsvx_
#define F_DSPTRD dsptrd_
#define F_DSPTRF dsptrf_
#define F_DSPTRI dsptri_
#define F_DSPTRS dsptrs_
#define F_DSTEBZ dstebz_
#define F_DSTEDC dstedc_
#define F_DSTEGR dstegr_
#define F_DSTEIN dstein_
#define F_DSTEQR dsteqr_
#define F_DSTERF dsterf_
#define F_DSTEV dstev_
#define F_DSTEVD dstevd_
#define F_DSTEVR dstevr_
#define F_DSTEVX dstevx_
#define F_DSYCON dsycon_
#define F_DSYEV dsyev_
#define F_DSYEVD dsyevd_
#define F_DSYEVR dsyevr_
#define F_DSYEVX dsyevx_
#define F_DSYGST dsygst_
#define F_DSYGV dsygv_
#define F_DSYGVD dsygvd_
#define F_DSYGVX dsygvx_
#define F_DSYRFS dsyrfs_
#define F_DSYSV dsysv_
#define F_DSYSVX dsysvx_
#define F_DSYTRD dsytrd_
#define F_DSYTRF dsytrf_
#define F_DSYTRI dsytri_
#define F_DSYTRS dsytrs_
#define F_DTBCON dtbcon_
#define F_DTBRFS dtbrfs_
#define F_DTBTRS dtbtrs_
#define F_DTGEVC dtgevc_
#define F_DTGEXC dtgexc_
#define F_DTGSEN dtgsen_
#define F_DTGSJA dtgsja_
#define F_DTGSNA dtgsna_
#define F_DTGSYL dtgsyl_
#define F_DTPCON dtpcon_
#define F_DTPRFS dtprfs_
#define F_DTPTRI dtptri_
#define F_DTPTRS dtptrs_
#define F_DTRCON dtrcon_
#define F_DTREVC dtrevc_
#define F_DTREXC dtrexc_
#define F_DTRRFS dtrrfs_
#define F_DTRSEN dtrsen_
#define F_DTRSNA dtrsna_
#define F_DTRSYL dtrsyl_
#define F_DTRTRI dtrtri_
#define F_DTRTRS dtrtrs_
#define F_DTZRQF dtzrqf_
#define F_DTZRZF dtzrzf_
#elif FC_SYMBOL==1
#define F_DGEEV dgeev
#define F_DGESV dgesv
#define F_DGETRF dgetrf
#define F_DGETRI dgetri
#define F_DPOTRF dpotrf
#define F_DPOTRI dpotri
#define F_DPOTRS dpotrs
#define F_DGESVD dgesvd
#define F_DSYEV dsyev
#define F_DBDSDC dbdsdc
#define F_DBDSQR dbdsqr
#define F_DDISNA ddisna
#define F_DGBBRD dgbbrd
#define F_DGBCON dgbcon
#define F_DGBEQU dgbequ
#define F_DGBRFS dgbrfs
#define F_DGBSV dgbsv
#define F_DGBSVX dgbsvx
#define F_DGBTRF dgbtrf
#define F_DGBTRS dgbtrs
#define F_DGEBAK dgebak
#define F_DGEBAL dgebal
#define F_DGEBRD dgebrd
#define F_DGECON dgecon
#define F_DGEEQU dgeequ
#define F_DGEES dgees
#define F_DGEESX dgeesx
#define F_DGEEV dgeev
#define F_DGEEVX dgeevx
#define F_DGEGS dgegs
#define F_DGEGV dgegv
#define F_DGEHRD dgehrd
#define F_DGELQF dgelqf
#define F_DGELS dgels
#define F_DGELSD dgelsd
#define F_DGELSS dgelss
#define F_DGELSX dgelsx
#define F_DGELSY dgelsy
#define F_DGEQLF dgeqlf
#define F_DGEQP3 dgeqp3
#define F_DGEQPF dgeqpf
#define F_DGEQRF dgeqrf
#define F_DGERFS dgerfs
#define F_DGERQF dgerqf
#define F_DGESDD dgesdd
#define F_DGESV dgesv
#define F_DGESVX dgesvx
#define F_DGETRF dgetrf
#define F_DGETRI dgetri
#define F_DGETRS dgetrs
#define F_DGGBAK dggbak
#define F_DGGBAL dggbal
#define F_DGGES dgges
#define F_DGGESX dggesx
#define F_DGGEV dggev
#define F_DGGEVX dggevx
#define F_DGGGLM dggglm
#define F_DGGHRD dgghrd
#define F_DGGLSE dgglse
#define F_DGGQRF dggqrf
#define F_DGGRQF dggrqf
#define F_DGGSVD dggsvd
#define F_DGGSVP dggsvp
#define F_DGTCON dgtcon
#define F_DGTRFS dgtrfs
#define F_DGTSV dgtsv
#define F_DGTSVX dgtsvx
#define F_DGTTRF dgttrf
#define F_DGTTRS dgttrs
#define F_DHGEQZ dhgeqz
#define F_DHSEIN dhsein
#define F_DHSEQR dhseqr
#define F_DOPGTR dopgtr
#define F_DOPMTR dopmtr
#define F_DORGBR dorgbr
#define F_DORGHR dorghr
#define F_DORGLQ dorglq
#define F_DORGQL dorgql
#define F_DORGQR dorgqr
#define F_DORGRQ dorgrq
#define F_DORGTR dorgtr
#define F_DORMBR dormbr
#define F_DORMHR dormhr
#define F_DORMLQ dormlq
#define F_DORMQL dormql
#define F_DORMQR dormqr
#define F_DORMR3 dormr3
#define F_DORMRQ dormrq
#define F_DORMRZ dormrz
#define F_DORMTR dormtr
#define F_DPBCON dpbcon
#define F_DPBEQU dpbequ
#define F_DPBRFS dpbrfs
#define F_DPBSTF dpbstf
#define F_DPBSV dpbsv
#define F_DPBSVX dpbsvx
#define F_DPBTRF dpbtrf
#define F_DPBTRS dpbtrs
#define F_DPOCON dpocon
#define F_DPOEQU dpoequ
#define F_DPORFS dporfs
#define F_DPOSV dposv
#define F_DPOSVX dposvx
#define F_DPOTRF dpotrf
#define F_DPOTRI dpotri
#define F_DPOTRS dpotrs
#define F_DPPCON dppcon
#define F_DPPEQU dppequ
#define F_DPPRFS dpprfs
#define F_DPPSV dppsv
#define F_DPPSVX dppsvx
#define F_DPPTRF dpptrf
#define F_DPPTRI dpptri
#define F_DPPTRS dpptrs
#define F_DPTCON dptcon
#define F_DPTEQR dpteqr
#define F_DPTRFS dptrfs
#define F_DPTSV dptsv
#define F_DPTSVX dptsvx
#define F_DPTTRF dpttrf
#define F_DPTTRS dpttrs
#define F_DSBEV dsbev
#define F_DSBEVD dsbevd
#define F_DSBEVX dsbevx
#define F_DSBGST dsbgst
#define F_DSBGV dsbgv
#define F_DSBGVD dsbgvd
#define F_DSBGVX dsbgvx
#define F_DSBTRD dsbtrd
#define F_DSGESV dsgesv
#define F_DSPCON dspcon
#define F_DSPEV dspev
#define F_DSPEVD dspevd
#define F_DSPEVX dspevx
#define F_DSPGST dspgst
#define F_DSPGV dspgv
#define F_DSPGVD dspgvd
#define F_DSPGVX dspgvx
#define F_DSPRFS dsprfs
#define F_DSPSV dspsv
#define F_DSPSVX dspsvx
#define F_DSPTRD dsptrd
#define F_DSPTRF dsptrf
#define F_DSPTRI dsptri
#define F_DSPTRS dsptrs
#define F_DSTEBZ dstebz
#define F_DSTEDC dstedc
#define F_DSTEGR dstegr
#define F_DSTEIN dstein
#define F_DSTEQR dsteqr
#define F_DSTERF dsterf
#define F_DSTEV dstev
#define F_DSTEVD dstevd
#define F_DSTEVR dstevr
#define F_DSTEVX dstevx
#define F_DSYCON dsycon
#define F_DSYEV dsyev
#define F_DSYEVD dsyevd
#define F_DSYEVR dsyevr
#define F_DSYEVX dsyevx
#define F_DSYGST dsygst
#define F_DSYGV dsygv
#define F_DSYGVD dsygvd
#define F_DSYGVX dsygvx
#define F_DSYRFS dsyrfs
#define F_DSYSV dsysv
#define F_DSYSVX dsysvx
#define F_DSYTRD dsytrd
#define F_DSYTRF dsytrf
#define F_DSYTRI dsytri
#define F_DSYTRS dsytrs
#define F_DTBCON dtbcon
#define F_DTBRFS dtbrfs
#define F_DTBTRS dtbtrs
#define F_DTGEVC dtgevc
#define F_DTGEXC dtgexc
#define F_DTGSEN dtgsen
#define F_DTGSJA dtgsja
#define F_DTGSNA dtgsna
#define F_DTGSYL dtgsyl
#define F_DTPCON dtpcon
#define F_DTPRFS dtprfs
#define F_DTPTRI dtptri
#define F_DTPTRS dtptrs
#define F_DTRCON dtrcon
#define F_DTREVC dtrevc
#define F_DTREXC dtrexc
#define F_DTRRFS dtrrfs
#define F_DTRSEN dtrsen
#define F_DTRSNA dtrsna
#define F_DTRSYL dtrsyl
#define F_DTRTRI dtrtri
#define F_DTRTRS dtrtrs
#define F_DTZRQF dtzrqf
#define F_DTZRZF dtzrzf
#elif FC_SYMBOL==3
#define F_DGEEV DGEEV
#define F_DGESV DGESV
#define F_DGETRF DGETRF
#define F_DGETRI DGETRI
#define F_DPOTRF DPOTRF
#define F_DPOTRI DPOTRI
#define F_DPOTRS DPOTRS
#define F_DGESVD DGESVD
#define F_DSYEV DSYEV
#define F_DBDSDC DBDSDC
#define F_DBDSQR DBDSQR
#define F_DDISNA DDISNA
#define F_DGBBRD DGBBRD
#define F_DGBCON DGBCON
#define F_DGBEQU DGBEQU
#define F_DGBRFS DGBRFS
#define F_DGBSV DGBSV
#define F_DGBSVX DGBSVX
#define F_DGBTRF DGBTRF
#define F_DGBTRS DGBTRS
#define F_DGEBAK DGEBAK
#define F_DGEBAL DGEBAL
#define F_DGEBRD DGEBRD
#define F_DGECON DGECON
#define F_DGEEQU DGEEQU
#define F_DGEES DGEES
#define F_DGEESX DGEESX
#define F_DGEEV DGEEV
#define F_DGEEVX DGEEVX
#define F_DGEGS DGEGS
#define F_DGEGV DGEGV
#define F_DGEHRD DGEHRD
#define F_DGELQF DGELQF
#define F_DGELS DGELS
#define F_DGELSD DGELSD
#define F_DGELSS DGELSS
#define F_DGELSX DGELSX
#define F_DGELSY DGELSY
#define F_DGEQLF DGEQLF
#define F_DGEQP3 DGEQP3
#define F_DGEQPF DGEQPF
#define F_DGEQRF DGEQRF
#define F_DGERFS DGERFS
#define F_DGERQF DGERQF
#define F_DGESDD DGESDD
#define F_DGESV DGESV
#define F_DGESVX DGESVX
#define F_DGETRF DGETRF
#define F_DGETRI DGETRI
#define F_DGETRS DGETRS
#define F_DGGBAK DGGBAK
#define F_DGGBAL DGGBAL
#define F_DGGES DGGES
#define F_DGGESX DGGESX
#define F_DGGEV DGGEV
#define F_DGGEVX DGGEVX
#define F_DGGGLM DGGGLM
#define F_DGGHRD DGGHRD
#define F_DGGLSE DGGLSE
#define F_DGGQRF DGGQRF
#define F_DGGRQF DGGRQF
#define F_DGGSVD DGGSVD
#define F_DGGSVP DGGSVP
#define F_DGTCON DGTCON
#define F_DGTRFS DGTRFS
#define F_DGTSV DGTSV
#define F_DGTSVX DGTSVX
#define F_DGTTRF DGTTRF
#define F_DGTTRS DGTTRS
#define F_DHGEQZ DHGEQZ
#define F_DHSEIN DHSEIN
#define F_DHSEQR DHSEQR
#define F_DOPGTR DOPGTR
#define F_DOPMTR DOPMTR
#define F_DORGBR DORGBR
#define F_DORGHR DORGHR
#define F_DORGLQ DORGLQ
#define F_DORGQL DORGQL
#define F_DORGQR DORGQR
#define F_DORGRQ DORGRQ
#define F_DORGTR DORGTR
#define F_DORMBR DORMBR
#define F_DORMHR DORMHR
#define F_DORMLQ DORMLQ
#define F_DORMQL DORMQL
#define F_DORMQR DORMQR
#define F_DORMR3 DORMR3
#define F_DORMRQ DORMRQ
#define F_DORMRZ DORMRZ
#define F_DORMTR DORMTR
#define F_DPBCON DPBCON
#define F_DPBEQU DPBEQU
#define F_DPBRFS DPBRFS
#define F_DPBSTF DPBSTF
#define F_DPBSV DPBSV
#define F_DPBSVX DPBSVX
#define F_DPBTRF DPBTRF
#define F_DPBTRS DPBTRS
#define F_DPOCON DPOCON
#define F_DPOEQU DPOEQU
#define F_DPORFS DPORFS
#define F_DPOSV DPOSV
#define F_DPOSVX DPOSVX
#define F_DPOTRF DPOTRF
#define F_DPOTRI DPOTRI
#define F_DPOTRS DPOTRS
#define F_DPPCON DPPCON
#define F_DPPEQU DPPEQU
#define F_DPPRFS DPPRFS
#define F_DPPSV DPPSV
#define F_DPPSVX DPPSVX
#define F_DPPTRF DPPTRF
#define F_DPPTRI DPPTRI
#define F_DPPTRS DPPTRS
#define F_DPTCON DPTCON
#define F_DPTEQR DPTEQR
#define F_DPTRFS DPTRFS
#define F_DPTSV DPTSV
#define F_DPTSVX DPTSVX
#define F_DPTTRF DPTTRF
#define F_DPTTRS DPTTRS
#define F_DSBEV DSBEV
#define F_DSBEVD DSBEVD
#define F_DSBEVX DSBEVX
#define F_DSBGST DSBGST
#define F_DSBGV DSBGV
#define F_DSBGVD DSBGVD
#define F_DSBGVX DSBGVX
#define F_DSBTRD DSBTRD
#define F_DSGESV DSGESV
#define F_DSPCON DSPCON
#define F_DSPEV DSPEV
#define F_DSPEVD DSPEVD
#define F_DSPEVX DSPEVX
#define F_DSPGST DSPGST
#define F_DSPGV DSPGV
#define F_DSPGVD DSPGVD
#define F_DSPGVX DSPGVX
#define F_DSPRFS DSPRFS
#define F_DSPSV DSPSV
#define F_DSPSVX DSPSVX
#define F_DSPTRD DSPTRD
#define F_DSPTRF DSPTRF
#define F_DSPTRI DSPTRI
#define F_DSPTRS DSPTRS
#define F_DSTEBZ DSTEBZ
#define F_DSTEDC DSTEDC
#define F_DSTEGR DSTEGR
#define F_DSTEIN DSTEIN
#define F_DSTEQR DSTEQR
#define F_DSTERF DSTERF
#define F_DSTEV DSTEV
#define F_DSTEVD DSTEVD
#define F_DSTEVR DSTEVR
#define F_DSTEVX DSTEVX
#define F_DSYCON DSYCON
#define F_DSYEV DSYEV
#define F_DSYEVD DSYEVD
#define F_DSYEVR DSYEVR
#define F_DSYEVX DSYEVX
#define F_DSYGST DSYGST
#define F_DSYGV DSYGV
#define F_DSYGVD DSYGVD
#define F_DSYGVX DSYGVX
#define F_DSYRFS DSYRFS
#define F_DSYSV DSYSV
#define F_DSYSVX DSYSVX
#define F_DSYTRD DSYTRD
#define F_DSYTRF DSYTRF
#define F_DSYTRI DSYTRI
#define F_DSYTRS DSYTRS
#define F_DTBCON DTBCON
#define F_DTBRFS DTBRFS
#define F_DTBTRS DTBTRS
#define F_DTGEVC DTGEVC
#define F_DTGEXC DTGEXC
#define F_DTGSEN DTGSEN
#define F_DTGSJA DTGSJA
#define F_DTGSNA DTGSNA
#define F_DTGSYL DTGSYL
#define F_DTPCON DTPCON
#define F_DTPRFS DTPRFS
#define F_DTPTRI DTPTRI
#define F_DTPTRS DTPTRS
#define F_DTRCON DTRCON
#define F_DTREVC DTREVC
#define F_DTREXC DTREXC
#define F_DTRRFS DTRRFS
#define F_DTRSEN DTRSEN
#define F_DTRSNA DTRSNA
#define F_DTRSYL DTRSYL
#define F_DTRTRI DTRTRI
#define F_DTRTRS DTRTRS
#define F_DTZRQF DTZRQF
#define F_DTZRZF DTZRZF
#elif FC_SYMBOL==4
#define F_DGEEV DGEEV_
#define F_DGESV DGESV_
#define F_DGETRF DGETRF_
#define F_DGETRI DGETRI_
#define F_DPOTRF DPOTRF_
#define F_DPOTRI DPOTRI_
#define F_DPOTRS DPOTRS_
#define F_DGESVD DGESVD_
#define F_DSYEV DSYEV_
#define F_DBDSDC DBDSDC_
#define F_DBDSQR DBDSQR_
#define F_DDISNA DDISNA_
#define F_DGBBRD DGBBRD_
#define F_DGBCON DGBCON_
#define F_DGBEQU DGBEQU_
#define F_DGBRFS DGBRFS_
#define F_DGBSV DGBSV_
#define F_DGBSVX DGBSVX_
#define F_DGBTRF DGBTRF_
#define F_DGBTRS DGBTRS_
#define F_DGEBAK DGEBAK_
#define F_DGEBAL DGEBAL_
#define F_DGEBRD DGEBRD_
#define F_DGECON DGECON_
#define F_DGEEQU DGEEQU_
#define F_DGEES DGEES_
#define F_DGEESX DGEESX_
#define F_DGEEV DGEEV_
#define F_DGEEVX DGEEVX_
#define F_DGEGS DGEGS_
#define F_DGEGV DGEGV_
#define F_DGEHRD DGEHRD_
#define F_DGELQF DGELQF_
#define F_DGELS DGELS_
#define F_DGELSD DGELSD_
#define F_DGELSS DGELSS_
#define F_DGELSX DGELSX_
#define F_DGELSY DGELSY_
#define F_DGEQLF DGEQLF_
#define F_DGEQP3 DGEQP3_
#define F_DGEQPF DGEQPF_
#define F_DGEQRF DGEQRF_
#define F_DGERFS DGERFS_
#define F_DGERQF DGERQF_
#define F_DGESDD DGESDD_
#define F_DGESV DGESV_
#define F_DGESVX DGESVX_
#define F_DGETRF DGETRF_
#define F_DGETRI DGETRI_
#define F_DGETRS DGETRS_
#define F_DGGBAK DGGBAK_
#define F_DGGBAL DGGBAL_
#define F_DGGES DGGES_
#define F_DGGESX DGGESX_
#define F_DGGEV DGGEV_
#define F_DGGEVX DGGEVX_
#define F_DGGGLM DGGGLM_
#define F_DGGHRD DGGHRD_
#define F_DGGLSE DGGLSE_
#define F_DGGQRF DGGQRF_
#define F_DGGRQF DGGRQF_
#define F_DGGSVD DGGSVD_
#define F_DGGSVP DGGSVP_
#define F_DGTCON DGTCON_
#define F_DGTRFS DGTRFS_
#define F_DGTSV DGTSV_
#define F_DGTSVX DGTSVX_
#define F_DGTTRF DGTTRF_
#define F_DGTTRS DGTTRS_
#define F_DHGEQZ DHGEQZ_
#define F_DHSEIN DHSEIN_
#define F_DHSEQR DHSEQR_
#define F_DOPGTR DOPGTR_
#define F_DOPMTR DOPMTR_
#define F_DORGBR DORGBR_
#define F_DORGHR DORGHR_
#define F_DORGLQ DORGLQ_
#define F_DORGQL DORGQL_
#define F_DORGQR DORGQR_
#define F_DORGRQ DORGRQ_
#define F_DORGTR DORGTR_
#define F_DORMBR DORMBR_
#define F_DORMHR DORMHR_
#define F_DORMLQ DORMLQ_
#define F_DORMQL DORMQL_
#define F_DORMQR DORMQR_
#define F_DORMR3 DORMR3_
#define F_DORMRQ DORMRQ_
#define F_DORMRZ DORMRZ_
#define F_DORMTR DORMTR_
#define F_DPBCON DPBCON_
#define F_DPBEQU DPBEQU_
#define F_DPBRFS DPBRFS_
#define F_DPBSTF DPBSTF_
#define F_DPBSV DPBSV_
#define F_DPBSVX DPBSVX_
#define F_DPBTRF DPBTRF_
#define F_DPBTRS DPBTRS_
#define F_DPOCON DPOCON_
#define F_DPOEQU DPOEQU_
#define F_DPORFS DPORFS_
#define F_DPOSV DPOSV_
#define F_DPOSVX DPOSVX_
#define F_DPOTRF DPOTRF_
#define F_DPOTRI DPOTRI_
#define F_DPOTRS DPOTRS_
#define F_DPPCON DPPCON_
#define F_DPPEQU DPPEQU_
#define F_DPPRFS DPPRFS_
#define F_DPPSV DPPSV_
#define F_DPPSVX DPPSVX_
#define F_DPPTRF DPPTRF_
#define F_DPPTRI DPPTRI_
#define F_DPPTRS DPPTRS_
#define F_DPTCON DPTCON_
#define F_DPTEQR DPTEQR_
#define F_DPTRFS DPTRFS_
#define F_DPTSV DPTSV_
#define F_DPTSVX DPTSVX_
#define F_DPTTRF DPTTRF_
#define F_DPTTRS DPTTRS_
#define F_DSBEV DSBEV_
#define F_DSBEVD DSBEVD_
#define F_DSBEVX DSBEVX_
#define F_DSBGST DSBGST_
#define F_DSBGV DSBGV_
#define F_DSBGVD DSBGVD_
#define F_DSBGVX DSBGVX_
#define F_DSBTRD DSBTRD_
#define F_DSGESV DSGESV_
#define F_DSPCON DSPCON_
#define F_DSPEV DSPEV_
#define F_DSPEVD DSPEVD_
#define F_DSPEVX DSPEVX_
#define F_DSPGST DSPGST_
#define F_DSPGV DSPGV_
#define F_DSPGVD DSPGVD_
#define F_DSPGVX DSPGVX_
#define F_DSPRFS DSPRFS_
#define F_DSPSV DSPSV_
#define F_DSPSVX DSPSVX_
#define F_DSPTRD DSPTRD_
#define F_DSPTRF DSPTRF_
#define F_DSPTRI DSPTRI_
#define F_DSPTRS DSPTRS_
#define F_DSTEBZ DSTEBZ_
#define F_DSTEDC DSTEDC_
#define F_DSTEGR DSTEGR_
#define F_DSTEIN DSTEIN_
#define F_DSTEQR DSTEQR_
#define F_DSTERF DSTERF_
#define F_DSTEV DSTEV_
#define F_DSTEVD DSTEVD_
#define F_DSTEVR DSTEVR_
#define F_DSTEVX DSTEVX_
#define F_DSYCON DSYCON_
#define F_DSYEV DSYEV_
#define F_DSYEVD DSYEVD_
#define F_DSYEVR DSYEVR_
#define F_DSYEVX DSYEVX_
#define F_DSYGST DSYGST_
#define F_DSYGV DSYGV_
#define F_DSYGVD DSYGVD_
#define F_DSYGVX DSYGVX_
#define F_DSYRFS DSYRFS_
#define F_DSYSV DSYSV_
#define F_DSYSVX DSYSVX_
#define F_DSYTRD DSYTRD_
#define F_DSYTRF DSYTRF_
#define F_DSYTRI DSYTRI_
#define F_DSYTRS DSYTRS_
#define F_DTBCON DTBCON_
#define F_DTBRFS DTBRFS_
#define F_DTBTRS DTBTRS_
#define F_DTGEVC DTGEVC_
#define F_DTGEXC DTGEXC_
#define F_DTGSEN DTGSEN_
#define F_DTGSJA DTGSJA_
#define F_DTGSNA DTGSNA_
#define F_DTGSYL DTGSYL_
#define F_DTPCON DTPCON_
#define F_DTPRFS DTPRFS_
#define F_DTPTRI DTPTRI_
#define F_DTPTRS DTPTRS_
#define F_DTRCON DTRCON_
#define F_DTREVC DTREVC_
#define F_DTREXC DTREXC_
#define F_DTRRFS DTRRFS_
#define F_DTRSEN DTRSEN_
#define F_DTRSNA DTRSNA_
#define F_DTRSYL DTRSYL_
#define F_DTRTRI DTRTRI_
#define F_DTRTRS DTRTRS_
#define F_DTZRQF DTZRQF_
#define F_DTZRZF DTZRZF_
#endif

#endif

#endif 