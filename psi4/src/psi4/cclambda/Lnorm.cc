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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

extern double pseudoenergy(struct L_Params L_params);

void Lnorm(struct L_Params L_params)
{
  dpdfile2 R1, L1, LIA, Lia, RIA, Ria;
  dpdbuf4 R2, L2, LIJAB, Lijab, LIjAb, RIJAB, Rijab, RIjAb;
  double tval, overlap, overlap0, overlap1, overlap2, L0;
  char R1A_lbl[32], R1B_lbl[32], R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];
  int L_irr;
  L_irr = L_params.irrep;

  if (L_params.ground)
    L0 = 1.0;
  else
    L0 = 0.0;

  sprintf(R1A_lbl, "RIA %d %d", L_irr, L_params.root);
  sprintf(R1B_lbl, "Ria %d %d", L_irr, L_params.root);
  sprintf(R2AA_lbl, "RIJAB %d %d", L_irr, L_params.root);
  sprintf(R2BB_lbl, "Rijab %d %d", L_irr, L_params.root);
  sprintf(R2AB_lbl, "RIjAb %d %d", L_irr, L_params.root);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    overlap0 = L0 * L_params.R0;
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");

    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    overlap1 = global_dpd_->file2_dot(&LIA, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1B_lbl);
    overlap1 += global_dpd_->file2_dot(&Lia, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    overlap2 = global_dpd_->buf4_dot(&LIJAB, &R2);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2BB_lbl);
    overlap2 += global_dpd_->buf4_dot(&Lijab, &R2);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 0, 5, 0, 5, 0, R2AB_lbl);
    overlap2 += global_dpd_->buf4_dot(&LIjAb, &R2);
    global_dpd_->buf4_close(&R2);
  }
  else {
    overlap0 = L0 * L_params.R0;
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");

    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 0, 1, R1A_lbl);
    overlap1 = global_dpd_->file2_dot(&LIA, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, L_irr, 2, 3, R1B_lbl);
    overlap1 += global_dpd_->file2_dot(&Lia, &R1);
    global_dpd_->file2_close(&R1);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 2, 7, 2, 7, 0, R2AA_lbl);
    overlap2 = global_dpd_->buf4_dot(&LIJAB, &R2);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 12, 17, 12, 17, 0, R2BB_lbl);
    overlap2 += global_dpd_->buf4_dot(&Lijab, &R2);
    global_dpd_->buf4_close(&R2);
    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, L_irr, 22, 28, 22, 28, 0, R2AB_lbl);
    overlap2 += global_dpd_->buf4_dot(&LIjAb, &R2);
    global_dpd_->buf4_close(&R2);
  }

  overlap = overlap0 + overlap1 + overlap2;

  outfile->Printf("\n\tInitial  <L|R>  =     %15.10lf\n", overlap);

  global_dpd_->file2_scm(&LIA, 1.0/overlap);
  global_dpd_->file2_scm(&Lia, 1.0/overlap);
  global_dpd_->buf4_scm(&LIJAB, 1.0/overlap);
  global_dpd_->buf4_scm(&Lijab, 1.0/overlap);
  global_dpd_->buf4_scm(&LIjAb, 1.0/overlap);

  outfile->Printf("\tNormalizing L...\n");
  outfile->Printf("\tL0 * R0 =     %15.10lf\n", overlap0/overlap);
  outfile->Printf("\tL1 * R1 =     %15.10lf\n", overlap1/overlap);
  outfile->Printf("\tL2 * R2 =     %15.10lf\n", overlap2/overlap);
  outfile->Printf("\t <L|R>  =     %15.10lf\n", overlap/overlap);

  global_dpd_->file2_close(&LIA);
  global_dpd_->file2_close(&Lia);
  global_dpd_->buf4_close(&LIJAB);
  global_dpd_->buf4_close(&Lijab);
  global_dpd_->buf4_close(&LIjAb);

  tval = pseudoenergy(L_params);
  outfile->Printf("\tPseudoenergy or Norm of normalized L = %20.15lf\n",tval);

  return;
}

}} // namespace psi::cclambda
