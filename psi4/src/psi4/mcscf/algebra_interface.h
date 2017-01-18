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

#ifndef _psi_src_bin_mcscf_algebra_interface_h_
#define _psi_src_bin_mcscf_algebra_interface_h_

#include "algebra_interface_mangle.h"

namespace psi{ namespace mcscf{

extern "C" void F_DAXPY(int *length, double *a, double *x, int *inc_x,
                    double *y, int *inc_y);
extern "C" void F_DCOPY(int *length, double *x, int *inc_x,
                    double *y, int *inc_y);
extern "C" void F_DGEMM(const char *transa, const char *transb, int *m, int *n, int *k,
                    double *alpha, double *A, int *lda, double *B, int *ldb,
                    double *beta, double *C, int *ldc);
extern "C" void F_DROT(int *ntot,double *x, int *incx,double *y, int *incy,
                  double *cotheta,double *sintheta);
extern "C" void F_DSCAL(int *n, double *alpha, double *vec, int *inc);
extern "C" void F_DGEMV(char *transa, int *m, int *n, double *alpha, double *A,
                    int *lda, double *X, int *inc_x, double *beta,
                    double *Y, int *inc_y);
extern "C" double F_DDOT(int *n, double *x, int *incx, double *y, int *incy);

void C_DGEMM_12(int m, int n, int k, double alpha,double *A, int nra,
                double *B, int ncb, double beta, double *C, int ncc);
void C_DGEMM_22(int m, int n, int k, double alpha,double *A, int nca,
                double *B, int ncb, double beta, double *C, int ncc);

// void C_DGEMM_11(int m, int n, int k, double alpha,double *A, int nca,
//                 double *B, int ncb, double beta, double *C, int ncc);
// void C_DGEMM_21(int m, int n, int k, double alpha,double *A, int nca,
//                 double *B, int ncb, double beta, double *C, int ncc);


extern "C" void F_DGEEV(const char *jobvl, const char *jobvr, int *n, double *a, int *lda,
                    double *wr, double *wi, double *vl, int *ldvl, double *vr,
                    int *ldvr, double *work, int *lwork, int *info);
extern "C" void F_DGESV(int *n, int *nrhs, double *A, int *lda, int *ipiv,
                    double *B, int *ldb, int *info);

extern "C" void F_DSYEV(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO );

}}

#endif // _psi_src_bin_mcscf_algebra_interface_h_