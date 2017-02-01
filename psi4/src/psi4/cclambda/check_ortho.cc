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
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

double LR_overlap_ROHF(int IRR, int L_index, int R_index);
double LR_overlap_RHF(int IRR, int L_index, int R_index);

void check_ortho(struct L_Params *pL_params) {
  int L_state_index, root_L_irr, L_irr;
  int R_state_index, root_R_irr, R_irr;
  double **O, tval;
  int L,R;


  if (params.ref <= 1) {
    O = block_matrix(params.nstates,params.nstates);
    for (L=0;L<params.nstates;++L) {
      L_irr = pL_params[L].irrep;
      root_L_irr = pL_params[L].root;

      for (R=0;R<params.nstates;++R) {
        R_irr = pL_params[R].irrep;
        root_R_irr = pL_params[R].root;

        if (L_irr == R_irr) {
          tval = LR_overlap_ROHF(L_irr, root_L_irr, root_R_irr);
          if (pL_params[L].ground)
            tval += 1.0 * pL_params[R].R0;
        }
        else
          tval = -99.0;

        O[L][R] = tval;
      }
    }
    outfile->Printf("\t<L|R> overlap matrix with ROHF quantities (-99 => 0 by symmetry)\n");
    print_mat(O, params.nstates, params.nstates, "outfile");
    free_block(O);
  }

  if (params.ref == 0) { /* test RHF quantities */
    O = block_matrix(params.nstates, params.nstates);
    for (L=0; L<params.nstates;++L) {
      L_irr = pL_params[L].irrep;
      root_L_irr = pL_params[L].root;

      for (R=0;R<params.nstates;++R) {
        R_irr = pL_params[R].irrep;
        root_R_irr = pL_params[R].root;

        if (L_irr == R_irr) {
          tval = LR_overlap_RHF(L_irr, root_L_irr, root_R_irr);
          if (pL_params[L].ground)
            tval += 1.0 * pL_params[R].R0;
        }
        else
          tval = -99.0;

        O[L][R] = tval;
      }
    }
    outfile->Printf("\t<L|R> overlap matrix with RHF quantities (-99 => 0 by symmetry)\n");
    print_mat(O, params.nstates, params.nstates, "outfile");
    free_block(O);
  }
  return;
}

double LR_overlap_ROHF(int IRR, int L_index, int R_index) {
  double overlap;
  dpdfile2 R1, L1;
  dpdbuf4 R2, L2;
  char R1A_lbl[32], R1B_lbl[32], R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];
  char L1A_lbl[32], L1B_lbl[32], L2AA_lbl[32], L2BB_lbl[32], L2AB_lbl[32];

  sprintf(R1A_lbl, "RIA %d %d", IRR, R_index);
  sprintf(R1B_lbl, "Ria %d %d", IRR, R_index);
  sprintf(R2AA_lbl, "RIJAB %d %d", IRR, R_index);
  sprintf(R2BB_lbl, "Rijab %d %d", IRR, R_index);
  sprintf(R2AB_lbl, "RIjAb %d %d", IRR, R_index);

  sprintf(L1A_lbl, "LIA %d %d", IRR, L_index);
  sprintf(L1B_lbl, "Lia %d %d", IRR, L_index);
  sprintf(L2AA_lbl, "LIJAB %d %d", IRR, L_index);
  sprintf(L2BB_lbl, "Lijab %d %d", IRR, L_index);
  sprintf(L2AB_lbl, "LIjAb %d %d", IRR, L_index);

  global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, IRR, 0, 1, R1A_lbl);
  global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, IRR, 0, 1, L1A_lbl);
  overlap = global_dpd_->file2_dot(&L1, &R1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&L1);

  global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, IRR, 0, 1, R1B_lbl);
  global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, IRR, 0, 1, L1B_lbl);
  overlap += global_dpd_->file2_dot(&L1, &R1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&L1);

  global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 2, 7, 2, 7, 0, R2AA_lbl);
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2AA_lbl);
  overlap += global_dpd_->buf4_dot(&L2, &R2);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&L2);

  global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 2, 7, 2, 7, 0, R2BB_lbl);
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2BB_lbl);
  overlap += global_dpd_->buf4_dot(&L2, &R2);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&L2);

  global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 0, 5, 0, 5, 0, R2AB_lbl);
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl);
  overlap += global_dpd_->buf4_dot(&L2, &R2);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&L2);

  return overlap;
}

double LR_overlap_RHF(int IRR, int L_index, int R_index) {
  dpdfile2 R1, L1;
  dpdbuf4 R2, L2;
  double overlap, overlap2, overlap3;
  char L1A_lbl[32], R1A_lbl[32], lbl[32];

  sprintf(L1A_lbl, "LIA %d %d", IRR, L_index);
  sprintf(R1A_lbl, "RIA %d %d", IRR, R_index);

  global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, IRR, 0, 1, R1A_lbl);
  global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, IRR, 0, 1, L1A_lbl);
  overlap = 2.0 * global_dpd_->file2_dot(&L1, &R1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&L1);

  sprintf(lbl, "2RIjAb - RIjbA %d %d", IRR, R_index);
  global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 0, 5, 0, 5, 0, lbl);

  sprintf(lbl, "LIjAb %d %d", IRR, L_index);
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, lbl);
  overlap2 = global_dpd_->buf4_dot(&L2, &R2);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&R2);

  sprintf(lbl, "2LIjAb - LIjbA %d %d", IRR, L_index);
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, lbl);

  sprintf(lbl, "RIjAb %d %d", IRR, R_index);
  global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 0, 5, 0, 5, 0, lbl);
  overlap3 = global_dpd_->buf4_dot(&L2, &R2);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&L2);

  if (fabs(overlap2 - overlap3) > 1E-14) {
    outfile->Printf("Bad anti-symmetry detected in RHF quantities\n");
    outfile->Printf("error: %15.10lf\n",overlap2-overlap3);
  }

  overlap += overlap2;
  return overlap;
}

}} // namespace psi::cclambda
