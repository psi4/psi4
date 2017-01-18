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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* relax_I_RHF(): Add the ROHF orbital-response contributions from
** the one-electron density matrix to the I(I,J) and I(I,A) blocks of
** the Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases. */

void relax_I_RHF(void)
{
  dpdfile2 I, D, f, I1, I2, I3;
  dpdbuf4 E, A, C, D4;
  int h, nirreps, i, j, e, *occpi, *virtpi, *openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;

  /* I(I,A) = I'(I,A) - sum_M f(I,M) D(orb)(A,M) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,A)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->contract222(&f, &D, &I, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&I);
  /* Add the contributions of dependent pairs (i,j) and (a,b) to the Lagrangian i
   * due to the use of canonical perturbed orbitals  */
  /*  Iji -= fjj * I'ij   */ 
  if (params.wfn == "CCSD_T" && params.dertype == 1){
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
  global_dpd_->file2_init(&I2, PSIF_CC_TMP, 0, 0, 0, "delta_I/delta_f_IJ");
  global_dpd_->file2_init(&I3, PSIF_CC_TMP, 0, 1, 1, "delta_I/delta_f_AB");
  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->contract222(&f, &I2, &I, 0, 0, -1, 1);
  global_dpd_->file2_close(&f);

  global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  global_dpd_->buf4_scmcopy(&A, PSIF_CC_MISC, "4 <kj|li> - <kj|il> - <ki|jl>", 4);
  global_dpd_->buf4_sort_axpy(&A, PSIF_CC_MISC, pqsr, 0, 0, "4 <kj|li> - <kj|il> - <ki|jl>", -1);
  global_dpd_->buf4_sort_axpy(&A, PSIF_CC_MISC, prsq, 0, 0, "4 <kj|li> - <kj|il> - <ki|jl>", -1);
  global_dpd_->buf4_close(&A);

  /*  Iji -= 1/2 sum_kl (4<kj|li> - <kj|il> - <ki|jl>)(I'kl - I'lk)/(fkk - fll)   */ 
  global_dpd_->buf4_init(&A, PSIF_CC_MISC, 0, 0, 0, 0, 0, 0, "4 <kj|li> - <kj|il> - <ki|jl>");
  global_dpd_->dot13(&I2, &A, &I, 0, 0, -0.5, 1.0);
  global_dpd_->buf4_close(&A);
  global_dpd_->file2_close(&I2);


  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 11, 11, 11, 11, 0, "C <ai|bj>");
  global_dpd_->buf4_scmcopy(&C, PSIF_CC_MISC, "4 <aj|bi> - <aj|ib> - <ai|jb>", 4);
  global_dpd_->buf4_init(&D4, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_sort_axpy(&D4, PSIF_CC_MISC, rqsp, 11, 11, "4 <aj|bi> - <aj|ib> - <ai|jb>", -1);
  global_dpd_->buf4_sort_axpy(&D4, PSIF_CC_MISC, rpsq, 11, 11, "4 <aj|bi> - <aj|ib> - <ai|jb>", -1);
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_close(&D4);

  /*  Iji -= 1/2 sum_ab (4<aj|bi> - <aj|ib> - <ai|jb>)(I'ab - I'ba)/(faa - fbb)   */ 
  global_dpd_->buf4_init(&C, PSIF_CC_MISC, 0, 11, 11, 11, 11, 0, "4 <aj|bi> - <aj|ib> - <ai|jb>");
  global_dpd_->dot13(&I3, &C, &I, 0, 0, -0.5, 1.0);
  global_dpd_->buf4_close(&C);
  global_dpd_->file2_close(&I3);
  global_dpd_->file2_close(&I);

  /*  Iab -= faa * I'ab   */ 
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 1, 1, "I'AB");
  global_dpd_->file2_init(&I1, PSIF_CC_TMP, 0, 1, 1, "delta_I/delta_f_AB");
  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->contract222(&f, &I1, &I, 0, 1, -1, 1);
  global_dpd_->file2_close(&I);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_close(&f);
   }
  /* RHF Case: I(i,j) = I'(i,j) - D(orb)(e,c) [4 <ei|mj> - <ei|jm> - <ej|im>] */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,J)");
  global_dpd_->file2_close(&I);

  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  global_dpd_->buf4_scmcopy(&E, PSIF_CC_EINTS, "4 <ei|mj> - <ei|jm> - <ej|im>", 4);
  global_dpd_->buf4_sort_axpy(&E, PSIF_CC_EINTS, pqsr, 11, 0, "4 <ei|mj> - <ei|jm> - <ej|im>", -1);
  global_dpd_->buf4_sort_axpy(&E, PSIF_CC_EINTS, psqr, 11, 0, "4 <ei|mj> - <ei|jm> - <ej|im>", -1);
  global_dpd_->buf4_close(&E);

  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "4 <ei|mj> - <ei|jm> - <ej|im>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_close(&I);
}

}} // namespace psi::ccdensity
