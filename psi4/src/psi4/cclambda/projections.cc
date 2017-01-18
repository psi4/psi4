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

/* projections() computes the projection of the EOM CCSD wavefunction onto the
reference, singles and doubles spaces, respectively.  These values may be
compared to those output by ACES2.  Also, an approximate excitation level (AEL)
may be defined as <Psi|S><S|Psi> + 2<Psi|D><D|Psi>.  For symmetric excitations,
this AEL theoretically might have a value slightly less than 1.  The sum of
these projections must equal 1 if R and L are normalized on this space.

<Psi|0><0|Psi> = <0|Le^(-T)|0><0|Re^T|0> =
    R0 * (-dot_L1T1 - dot_L2T2 + 0.5*dot_L2T1T1);
<Psi|S><S|Psi> = <0|Le^(-T)|S><S|Re^T|S> =
    R0 * (+dot_L1T1 - dot_L2T1T1) + dot_L1R1 - dot_L2T1R1;
<Psi|D><D|Psi> = <0|Le^(-T)|D><D|Re^T|D> =
    R0 * (+dot_L2T2 + 0.5*dot_L2T1T1) + dot_L2T1R1 + dot_L2R2;

where dot_L1T1 = Lme Tme               dot_L2T2 = Lmnef Tmnef
      dot_L2T1T1 = Lmnef Tme Tnf       dot_L2T1R1 = Lmnef Tme Rnf
      dot_L1R1 = Lme Rme               dot_L2R2 = Lmnef Rmnef
*/

void projections(struct L_Params *pL_params) {
  int IRR;
  int i,j, root;
  double R0, projection_0, projection_S, projection_D, projection_tot, ael;
  double dot_L1T1=0, dot_L2T2=0, dot_L2T1T1=0, dot_L1R1=0, dot_L2T1R1=0;
  double dot_L2R2=0;
  char R1A_lbl[32], R1B_lbl[32], R2AA_lbl[32], R2BB_lbl[32], R2AB_lbl[32];
  char L1A_lbl[32], L1B_lbl[32], L2AA_lbl[32], L2BB_lbl[32], L2AB_lbl[32], L2AB_lbl2[32];
  dpdfile2 L1A, L1B, T1A, T1B, I1, R1A, R1B;
  dpdbuf4 T2, L2, R2;
  /* assumes that the excited-Rs are available */

  if (params.ref == 0) {
    global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  }
  else if (params.ref == 1) {
    global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  }
  else if (params.ref == 2) {
    global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  }

  for (i=1; i<params.nstates;++i) {
    psio_close(PSIF_CC_TMP, 0);
    psio_open(PSIF_CC_TMP, PSIO_OPEN_NEW);

    IRR = pL_params[i].irrep;
    R0 = pL_params[i].R0;
    root = pL_params[i].root;
    sprintf(R1A_lbl, "RIA %d %d", IRR, root);
    sprintf(R1B_lbl, "Ria %d %d", IRR, root);
    sprintf(R2AA_lbl, "RIJAB %d %d", IRR, root);
    sprintf(R2BB_lbl, "Rijab %d %d", IRR, root);
    sprintf(R2AB_lbl, "RIjAb %d %d", IRR, root);
    sprintf(L1A_lbl, "LIA %d %d", IRR, root);
    sprintf(L1B_lbl, "Lia %d %d", IRR, root);
    sprintf(L2AA_lbl, "LIJAB %d %d", IRR, root);
    sprintf(L2BB_lbl, "Lijab %d %d", IRR, root);
    sprintf(L2AB_lbl, "LIjAb %d %d", IRR, root);
    sprintf(L2AB_lbl2, "2LIjAb - LIjbA %d %d", IRR, root);

    if (params.ref == 0) {
      global_dpd_->file2_init(&L1A, PSIF_CC_LAMPS, IRR, 0, 1, L1A_lbl);
    }
    else if (params.ref == 1) {
      global_dpd_->file2_init(&L1A, PSIF_CC_LAMPS, IRR, 0, 1, L1A_lbl);
      global_dpd_->file2_init(&L1B, PSIF_CC_LAMPS, IRR, 0, 1, L1B_lbl);
    }
    else if (params.ref == 2) {
      global_dpd_->file2_init(&L1A, PSIF_CC_LAMPS, IRR, 0, 1, L1A_lbl);
      global_dpd_->file2_init(&L1B, PSIF_CC_LAMPS, IRR, 2, 3, L1B_lbl);
    }

    if (IRR == 0) {
      /* dot_L1T1 = <0|Lme Tme|0>, assuming L0=0 (excited state) */
      if (params.ref == 0) {
        dot_L1T1 = 2.0 * global_dpd_->file2_dot(&L1A, &T1A);
      }
      else if (params.ref >= 1) {
        dot_L1T1 = global_dpd_->file2_dot(&L1A, &T1A);
        dot_L1T1 += global_dpd_->file2_dot(&L1B, &T1B);
      }

      /* dot_L2T2 = <0|Lmnef Tmnef|0> */
      if (params.ref == 0) {
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
        dot_L2T2 = global_dpd_->buf4_dot(&L2, &T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&L2);
      }
      else if (params.ref == 1) {
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2AA_lbl);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
        dot_L2T2 = global_dpd_->buf4_dot(&L2, &T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&L2);
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2BB_lbl);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
        dot_L2T2 += global_dpd_->buf4_dot(&L2, &T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&L2);
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        dot_L2T2 += global_dpd_->buf4_dot(&L2, &T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&L2);
      }
      else if (params.ref == 2) {
       global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2AA_lbl);
       global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
       dot_L2T2 = global_dpd_->buf4_dot(&L2, &T2);
       global_dpd_->buf4_close(&T2);
       global_dpd_->buf4_close(&L2);
       global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 12, 17, 12, 17, 0, L2BB_lbl);
       global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
       dot_L2T2 += global_dpd_->buf4_dot(&L2, &T2);
       global_dpd_->buf4_close(&T2);
       global_dpd_->buf4_close(&L2);
       global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 22, 28, 22, 28, 0, L2AB_lbl);
       global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
       dot_L2T2 += global_dpd_->buf4_dot(&L2, &T2);
       global_dpd_->buf4_close(&T2);
       global_dpd_->buf4_close(&L2);
     }
     /* dot_L2T1T1 = TNF (TME LMNEF) + Tnf (Tme Lmnef) + 2*Tnf(TME LMnEf)  */
     if (params.ref == 0) {
        global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X(N,F)");
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl2);
        global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&L2);
        dot_L2T1T1 = 2.0 * global_dpd_->file2_dot(&T1A, &I1);
        global_dpd_->file2_close(&I1);
     }
     else if (params.ref == 1) {
        global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X(N,F)");
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 2, 7, 0, L2AA_lbl);
        global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&L2);
        dot_L2T1T1 = global_dpd_->file2_dot(&T1A, &I1);
        global_dpd_->file2_close(&I1);
        global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X(n,f)");
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 2, 7, 0, L2BB_lbl);
        global_dpd_->dot13(&T1B, &L2, &I1, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&L2);
        dot_L2T1T1 += global_dpd_->file2_dot(&T1B, &I1);
        global_dpd_->file2_close(&I1);
        global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X(n,f)");
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl);
        global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&L2);
        dot_L2T1T1 += 2.0 * global_dpd_->file2_dot(&T1B, &I1);
        global_dpd_->file2_close(&I1);
      }
      else if (params.ref == 2) {
        global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X(N,F)");
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 2, 7, 0, L2AA_lbl);
        global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&L2);
        dot_L2T1T1 = global_dpd_->file2_dot(&T1A, &I1);
        global_dpd_->file2_close(&I1);
        global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 2, 3, "X(n,f)");
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 10, 15, 12, 17, 0, L2BB_lbl);
        global_dpd_->dot13(&T1B, &L2, &I1, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&L2);
        dot_L2T1T1 += global_dpd_->file2_dot(&T1B, &I1);
        global_dpd_->file2_close(&I1);
        global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 2, 3, "X(n,f)");
        global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 22, 28, 22, 28, 0, L2AB_lbl);
        global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&L2);
        dot_L2T1T1 += 2.0 * global_dpd_->file2_dot(&T1B, &I1);
        global_dpd_->file2_close(&I1);
      }
    }

    /* open R files */
    if (params.ref == 0) {
      global_dpd_->file2_init(&R1A, PSIF_CC_RAMPS, IRR, 0, 1, R1A_lbl);
    }
    else if (params.ref == 1) {
      global_dpd_->file2_init(&R1A, PSIF_CC_RAMPS, IRR, 0, 1, R1A_lbl);
      global_dpd_->file2_init(&R1B, PSIF_CC_RAMPS, IRR, 0, 1, R1B_lbl);
    }
    else if (params.ref == 2) {
      global_dpd_->file2_init(&R1A, PSIF_CC_RAMPS, IRR, 0, 1, R1A_lbl);
      global_dpd_->file2_init(&R1B, PSIF_CC_RAMPS, IRR, 2, 3, R1B_lbl);
    }

    /* dot_L1R1 = Lme Rme */
    if (params.ref == 0) {
      dot_L1R1 = 2.0 * global_dpd_->file2_dot(&L1A,&R1A);
    }
    else if (params.ref >= 1) {
      dot_L1R1 = global_dpd_->file2_dot(&L1A,&R1A);
      dot_L1R1 += global_dpd_->file2_dot(&L1B,&R1B);
    }

    /* dot_L2R2 = <0|Lmnef Rmnef|0> */
    if (params.ref == 0) {
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl2);
      global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 0, 5, 0, 5, 0, R2AB_lbl);
      dot_L2R2 = global_dpd_->buf4_dot(&L2, &R2);
      global_dpd_->buf4_close(&R2);
      global_dpd_->buf4_close(&L2);
    }
    else if (params.ref == 1) {
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2AA_lbl);
      global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 2, 7, 2, 7, 0, R2AA_lbl);
      dot_L2R2 = global_dpd_->buf4_dot(&L2, &R2);
      global_dpd_->buf4_close(&R2);
      global_dpd_->buf4_close(&L2);
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2BB_lbl);
      global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 2, 7, 2, 7, 0, R2BB_lbl);
      dot_L2R2 += global_dpd_->buf4_dot(&L2, &R2);
      global_dpd_->buf4_close(&R2);
      global_dpd_->buf4_close(&L2);
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl);
      global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 0, 5, 0, 5, 0, R2AB_lbl);
      dot_L2R2 += global_dpd_->buf4_dot(&L2, &R2);
      global_dpd_->buf4_close(&R2);
      global_dpd_->buf4_close(&L2);
    }
    else if (params.ref == 2) {
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 2, 7, 2, 7, 0, L2AA_lbl);
      global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 2, 7, 2, 7, 0, R2AA_lbl);
      dot_L2R2 = global_dpd_->buf4_dot(&L2, &R2);
      global_dpd_->buf4_close(&R2);
      global_dpd_->buf4_close(&L2);
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 12, 17, 12, 17, 0, L2BB_lbl);
      global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 12, 17, 12, 17, 0, R2BB_lbl);
      dot_L2R2 += global_dpd_->buf4_dot(&L2, &R2);
      global_dpd_->buf4_close(&R2);
      global_dpd_->buf4_close(&L2);
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 22, 28, 22, 28, 0, L2AB_lbl);
      global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, IRR, 22, 28, 22, 28, 0, R2AB_lbl);
      dot_L2R2 += global_dpd_->buf4_dot(&L2, &R2);
      global_dpd_->buf4_close(&R2);
      global_dpd_->buf4_close(&L2);
    }

    /* dot_L2T1R1 = RNF (TME LMNEF) + Rnf (Tme Lmnef) */
    /*            + RNF (Tme LNmFe) + Rnf (TME LMnEf) */
    if (params.ref == 0) {
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl2);
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X2(N,F)");
      global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 = 2.0 * global_dpd_->file2_dot(&R1A, &I1);
      global_dpd_->file2_close(&I1);
    }
    else if (params.ref == 1) {
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X2(N,F)");
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 2, 7, 0, L2AA_lbl);
      global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 = global_dpd_->file2_dot(&R1A, &I1);
      global_dpd_->file2_close(&I1);
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X2(n,f)");
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 2, 7, 0, L2BB_lbl);
      global_dpd_->dot13(&T1B, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 += global_dpd_->file2_dot(&R1B, &I1);
      global_dpd_->file2_close(&I1);
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl);
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X2(N,F)");
      global_dpd_->dot24(&T1B, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 += global_dpd_->file2_dot(&R1A, &I1);
      global_dpd_->file2_close(&I1);
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X2(n,f)");
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 0, 5, 0, L2AB_lbl);
      global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 += global_dpd_->file2_dot(&R1B, &I1);
      global_dpd_->file2_close(&I1);
    }
    else if (params.ref == 2) {
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X2(N,F)");
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 0, 5, 2, 7, 0, L2AA_lbl);
      global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 = global_dpd_->file2_dot(&R1A, &I1);
      global_dpd_->file2_close(&I1);
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 2, 3, "X2(n,f)");
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 10, 15, 12, 17, 0, L2BB_lbl);
      global_dpd_->dot13(&T1B, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 += global_dpd_->file2_dot(&R1B, &I1);
      global_dpd_->file2_close(&I1);
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 22, 28, 22, 28, 0, L2AB_lbl);
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 0, 1, "X2(N,F)");
      global_dpd_->dot24(&T1B, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 += global_dpd_->file2_dot(&R1A, &I1);
      global_dpd_->file2_close(&I1);
      global_dpd_->file2_init(&I1, PSIF_CC_TMP, IRR, 2, 3, "X2(n,f)");
      global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, IRR, 22, 28, 22, 28, 0, L2AB_lbl);
      global_dpd_->dot13(&T1A, &L2, &I1, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&L2);
      dot_L2T1R1 += global_dpd_->file2_dot(&R1B, &I1);
      global_dpd_->file2_close(&I1);
    }

    /* close L1, R1 */
    if (params.ref == 0) {
      global_dpd_->file2_close(&R1A); global_dpd_->file2_close(&L1A);
    }
    else if (params.ref >= 1) {
      global_dpd_->file2_close(&R1A); global_dpd_->file2_close(&R1B);
      global_dpd_->file2_close(&L1A); global_dpd_->file2_close(&L1B);
    }

    projection_0 = R0 * (-dot_L1T1 - dot_L2T2 + 0.5*dot_L2T1T1);
    projection_S = R0 * (+dot_L1T1 - dot_L2T1T1) + dot_L1R1 - dot_L2T1R1;
    projection_D = R0 * (+dot_L2T2 + 0.5*dot_L2T1T1) + dot_L2T1R1 + dot_L2R2;
    projection_tot = projection_0 + projection_S + projection_D;
    ael = projection_S + 2.0 * projection_D;

    outfile->Printf("\n\tProjections for excited state, irrep %s, root %d:\n", moinfo.labels[0], root);
    outfile->Printf("\t<0|Le^(-T)|0><0|Re^T|0>  = %15.10lf\n", projection_0);
    outfile->Printf("\t<0|Le^(-T)|S><S|Re^T|0>  = %15.10lf\n", projection_S);
    outfile->Printf("\t<0|Le^(-T)|D><D|Re^T|0>  = %15.10lf\n", projection_D);
    outfile->Printf("\tSum of above             = %15.10lf\n", projection_tot);
    outfile->Printf("\tApprox. excitation level = %15.10lf\n", ael);
  } /* sum over states */

  /* close T1 */
  global_dpd_->file2_close(&T1A);
  if (params.ref >= 1)
    global_dpd_->file2_close(&T1B);
}

}} // namespace psi::cclambda
