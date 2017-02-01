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
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void ex_tdensity_rohf(struct TD_Params S, struct TD_Params U);
void ex_tdensity_uhf(struct TD_Params S, struct TD_Params U);
void ex_tdensity_intermediates(struct TD_Params S, struct TD_Params U);
void ex_sort_td_rohf(char hand, int Tirrep);
void ex_sort_td_uhf(char hand, int Tirrep);

void ex_tdensity_rohf(struct TD_Params S, struct TD_Params U)
{
  dpdfile2 DAI, Dai, DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, Int, XIJ, Xij, R1;
  dpdbuf4 T2, L2, R2, I2;
  dpdfile2 D, T1, L1, Z;

  // For generalization, may want to put an if-check here that
  // allows me to set the LHS DPD file number to either ground
  // or whichever excited state.

  /* Set LHS to Excited State */
  int LHS = PSIF_CC_GL;
  int Lirrep = S.irrep;
  int Rirrep = U.irrep;
  int Tirrep = Lirrep^Rirrep;
  /* Auxiliary Irrep Handles */
     // -> This should always be Lirrep
  int Airrep = Lirrep;

/*
  /// These have no R-dependence, so don't include in XTD's.
  if(Tirrep == 0) {
    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 0, "LTDIJ");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 0, "LTDij");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 1, 1, "LTDAB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 1, 1, "LTDab");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);
  }
*/

/*
  // Double Check LHS and RHS
  printf("*** Are the LHS and RHS ok?\n");
  global_dpd_->file2_init(&L1, PSIF_CC_GL, Lirrep, 0, 1, "LIA");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
  global_dpd_->file2_print(&L1, stdout);
  global_dpd_->file2_print(&R1, stdout);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&R1);
*/


  /* R_I^A */
     //-> Exclude this term for excited -> excited transition densities.
     //   This the R-only term (no L, so doesnt' contribute).
/*
  global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
  global_dpd_->file2_copy(&R1, PSIF_CC_TMP, "LTDIA");
  global_dpd_->file2_close(&R1);
*/

/*
  if(Tirrep == 0) {
    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 1, "LTDIA");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);

    //  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A)

    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);

    // T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E)

    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(A,E)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);
  }
*/

  /* R_i^a */
     //-> Exclude this term for excited -> excited transition densities.
     //   This the R-only term (no L, so doesnt' contribute).
/*
  global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "Ria");
  global_dpd_->file2_copy(&R1, PSIF_CC_TMP, "LTDia");
  global_dpd_->file2_close(&R1);
*/

/*
  if(Tirrep == 0) {
    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 1, "LTDia");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(i,m)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(i,m)");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(a,e)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);

    // Note that these blocks are still stored occ/vir
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_TMP, "LTDAI");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_TMP, "LTDai");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->file2_scm(&L1, (1/S.R0));
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "Lia");
    global_dpd_->file2_scm(&L1, (1/S.R0));
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_scm(&L2, (1/S.R0));
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_scm(&L2, (1/S.R0));
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_scm(&L2, (1/S.R0));
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_scm(&L2, (1/S.R0));
    global_dpd_->buf4_close(&L2);
  }
*/

  ex_tdensity_intermediates(S,U);

  global_dpd_->file2_init(&TIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&Tia, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_init(&RIA, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
  global_dpd_->file2_init(&Ria, PSIF_CC_GR, Rirrep, 0, 1, "Ria");
  global_dpd_->file2_init(&LIA, LHS, Lirrep, 0, 1, "LIA");
  global_dpd_->file2_init(&Lia, LHS, Lirrep, 0, 1, "Lia");

  /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */

  global_dpd_->file2_init(&DIJ, PSIF_CC_TMP, Tirrep, 0, 0, "LTDIJ");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 0, "LR_OO");
  global_dpd_->file2_axpy(&Int, &DIJ, -1.0, 1);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->contract222(&TIA, &Int, &DIJ, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&DIJ);

  global_dpd_->file2_init(&Dij, PSIF_CC_TMP, Tirrep, 0, 0, "LTDij");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 0, "LR_oo");
  global_dpd_->file2_axpy(&Int, &Dij, -1.0, 1);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_ov");
  global_dpd_->contract222(&Tia, &Int, &Dij, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&Dij);

  /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */

  global_dpd_->file2_init(&DAB, PSIF_CC_TMP, Tirrep, 1, 1, "LTDAB");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 1, 1, "LR_VV");
  global_dpd_->file2_axpy(&Int, &DAB, 1.0, 0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->contract222(&Int, &TIA, &DAB, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&DAB);

  global_dpd_->file2_init(&Dab, PSIF_CC_TMP, Tirrep, 1, 1, "LTDab");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 1, 1, "LR_vv");
  global_dpd_->file2_axpy(&Int, &Dab, 1.0, 0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_ov");
  global_dpd_->contract222(&Int, &Tia, &Dab, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&Dab);

  /* D[a][i] = +L2R1_ov[i][a] */

  global_dpd_->file2_init(&DAI, PSIF_CC_TMP, Tirrep, 0, 1, "LTDAI");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->file2_axpy(&Int, &DAI, 1.0, 0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&DAI);

  global_dpd_->file2_init(&Dai, PSIF_CC_TMP, Tirrep, 0, 1, "LTDai");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_ov");
  global_dpd_->file2_axpy(&Int, &Dai, 1.0, 0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&Dai);

  global_dpd_->file2_init(&DIA, PSIF_CC_TMP, Tirrep, 0, 1, "LTDIA");
  global_dpd_->file2_init(&Dia, PSIF_CC_TMP, Tirrep, 0, 1, "LTDia");

  /* D[i][a] = L1R2_ov[i][a] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L1R2_OV");
  global_dpd_->file2_axpy(&Int, &DIA, 1.0, 0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L1R2_ov");
  global_dpd_->file2_axpy(&Int, &Dia, 1.0, 0);
  global_dpd_->file2_close(&Int);

  /* - LR_OO[M][I] * t1[M][A] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 0, "LR_OO");
  global_dpd_->contract222(&Int, &TIA, &DIA, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 0, "LR_oo");
  global_dpd_->contract222(&Int, &Tia, &Dia, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  /* - t1[I][E] * LR_vv[E][A] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 1, 1, "LR_VV");
  global_dpd_->contract222(&TIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 1, 1, "LR_vv");
  global_dpd_->contract222(&Tia, &Int, &Dia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  /* - LT2_OO[M][I] * r1[M][A] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Airrep, 0, 0, "LT2_OO");
  global_dpd_->contract222(&Int, &RIA, &DIA, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Airrep, 0, 0, "LT2_oo");
  global_dpd_->contract222(&Int, &Ria, &Dia, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  /* - r1[I][E] * LT2_VV[E][A] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Airrep, 1, 1, "LT2_VV");
  global_dpd_->contract222(&RIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Airrep, 1, 1, "LT2_vv");
  global_dpd_->contract222(&Ria, &Int, &Dia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  /* + L2R1_ov[M][E] * t2[i][m][a][e] */

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_ov");
  global_dpd_->dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_ov");
  global_dpd_->dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->buf4_close(&T2);

  /* - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->file2_init(&XIJ, PSIF_EOM_TMP, Tirrep, 0, 0, "XIJ");
  global_dpd_->contract222(&TIA, &Int, &XIJ, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&XIJ, PSIF_EOM_TMP, Tirrep, 0, 0, "XIJ");
  global_dpd_->contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&XIJ);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_ov");
  global_dpd_->file2_init(&Xij, PSIF_EOM_TMP, Tirrep, 0, 0, "Xij");
  global_dpd_->contract222(&Tia, &Int, &Xij, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Xij, PSIF_EOM_TMP, Tirrep, 0, 0, "Xij");
  global_dpd_->contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Xij);

  global_dpd_->file2_close(&DIA);
  global_dpd_->file2_close(&Dia);

  global_dpd_->file2_close(&TIA);
  global_dpd_->file2_close(&Tia);
  global_dpd_->file2_close(&RIA);
  global_dpd_->file2_close(&Ria);
  global_dpd_->file2_close(&LIA);
  global_dpd_->file2_close(&Lia);

  return;
}

void ex_tdensity_uhf(struct TD_Params S, struct TD_Params U)
{
  dpdfile2 DAI, Dai, DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, Int, XIJ, Xij, R1;
  dpdbuf4 T2, L2, R2, I2;
  dpdfile2 D, T1, L1, Z;

  /* Set LHS to Excited State */
  int LHS = PSIF_CC_GL;
  int Lirrep = S.irrep;
  int Rirrep = U.irrep;
  int Tirrep = Lirrep^Rirrep;
  /* Auxiliary Irrep Handles */
  int Airrep = Lirrep;

/*
  if(Tirrep == 0 ) { // Symmetric Transition

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 0, "LTDIJ");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 2, 2, "LTDij");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 1, 1, "LTDAB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 3, 3, "LTDab");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&D);
  }
*/

  /* R_I^A */
     //-> Exclude this term for excited -> excited transition densities.
     //   This the R-only term (no L, so doesnt' contribute).
/*
  global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
  global_dpd_->file2_copy(&R1, PSIF_CC_TMP, "LTDIA");
  global_dpd_->file2_close(&R1);
*/

/*
  if(Tirrep == 0) { // Symmetric Transitions

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 1, "LTDIA");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);

    //  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A)

    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    // T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E)
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(A,E)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);
  }
*/

  /* R_i^a */
     //-> Exclude this term for excited -> excited transition densities.
     //   This the R-only term (no L, so doesnt' contribute).
/*
  global_dpd_->file2_init(&R1, PSIF_CC_GR, Rirrep, 2, 3, "Ria");
  global_dpd_->file2_copy(&R1, PSIF_CC_TMP, "LTDia");
  global_dpd_->file2_close(&R1);
*/

/*
  if(Tirrep == 0) { // Symmetric Transitions

    global_dpd_->file2_init(&D, PSIF_CC_TMP, 0, 2, 3, "LTDia");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 2, 2, "Z(i,m)");
    global_dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 2, 2, "Z(i,m)");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 3, 3, "Z(a,e)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&Z);
    global_dpd_->file2_close(&D);

    // Note that these blocks are still stored occ/vir
    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_CC_TMP, "LTDAI");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_CC_TMP, "LTDai");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, LHS, Lirrep, 0, 1, "LIA");
    global_dpd_->file2_scm(&L1, (1/S.R0));
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&L1, LHS, Lirrep, 2, 3, "Lia");
    global_dpd_->file2_scm(&L1, (1/S.R0));
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_scm(&L2, (1/S.R0));
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_scm(&L2, (1/S.R0));
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_scm(&L2, (1/S.R0));
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, LHS, Lirrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_scm(&L2, (1/S.R0));
    global_dpd_->buf4_close(&L2);
  }
*/

  ex_tdensity_intermediates(S,U);

  global_dpd_->file2_init(&TIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&Tia, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_init(&RIA, PSIF_CC_GR, Rirrep, 0, 1, "RIA");
  global_dpd_->file2_init(&Ria, PSIF_CC_GR, Rirrep, 2, 3, "Ria");
  global_dpd_->file2_init(&LIA, LHS, Lirrep, 0, 1, "LIA");
  global_dpd_->file2_init(&Lia, LHS, Lirrep, 2, 3, "Lia");

  /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */

  global_dpd_->file2_init(&DIJ, PSIF_CC_TMP, Tirrep, 0, 0, "LTDIJ");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 0, "LR_OO");
  global_dpd_->file2_axpy(&Int, &DIJ, -1.0, 1);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->contract222(&TIA, &Int, &DIJ, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&DIJ);

  global_dpd_->file2_init(&Dij, PSIF_CC_TMP, Tirrep, 2, 2, "LTDij");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 2, "LR_oo");
  global_dpd_->file2_axpy(&Int, &Dij, -1.0, 1);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 3, "L2R1_ov");
  global_dpd_->contract222(&Tia, &Int, &Dij, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&Dij);

  /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */

  global_dpd_->file2_init(&DAB, PSIF_CC_TMP, Tirrep, 1, 1, "LTDAB");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 1, 1, "LR_VV");
  global_dpd_->file2_axpy(&Int, &DAB, 1.0, 0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->contract222(&Int, &TIA, &DAB, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&DAB);

  global_dpd_->file2_init(&Dab, PSIF_CC_TMP, Tirrep, 3, 3, "LTDab");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 3, 3, "LR_vv");
  global_dpd_->file2_axpy(&Int, &Dab, 1.0, 0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 3, "L2R1_ov");
  global_dpd_->contract222(&Int, &Tia, &Dab, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&Dab);

  /* D[a][i] = +L2R1_ov[i][a] */

  global_dpd_->file2_init(&DAI, PSIF_CC_TMP, Tirrep, 0, 1, "LTDAI");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->file2_axpy(&Int, &DAI, 1.0, 0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&DAI);

  global_dpd_->file2_init(&Dai, PSIF_CC_TMP, Tirrep, 2, 3, "LTDai");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 3, "L2R1_ov");
  global_dpd_->file2_axpy(&Int, &Dai, 1.0, 0);
  global_dpd_->file2_close(&Int);
  global_dpd_->file2_close(&Dai);

  global_dpd_->file2_init(&DIA, PSIF_CC_TMP, Tirrep, 0, 1, "LTDIA");

  global_dpd_->file2_init(&Dia, PSIF_CC_TMP, Tirrep, 2, 3, "LTDia");

  /* D[i][a] = L1R2_ov[i][a] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L1R2_OV");
  global_dpd_->file2_axpy(&Int, &DIA, 1.0, 0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 3, "L1R2_ov");
  global_dpd_->file2_axpy(&Int, &Dia, 1.0, 0);
  global_dpd_->file2_close(&Int);

  /* - LR_OO[M][I] * t1[M][A] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 0, "LR_OO");
  global_dpd_->contract222(&Int, &TIA, &DIA, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 2, "LR_oo");
  global_dpd_->contract222(&Int, &Tia, &Dia, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  /* - t1[I][E] * LR_vv[E][A] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 1, 1, "LR_VV");
  global_dpd_->contract222(&TIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 3, 3, "LR_vv");
  global_dpd_->contract222(&Tia, &Int, &Dia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  /* - LT2_OO[M][I] * r1[M][A] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Airrep, 0, 0, "LT2_OO");
  global_dpd_->contract222(&Int, &RIA, &DIA, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Airrep, 2, 2, "LT2_oo");
  global_dpd_->contract222(&Int, &Ria, &Dia, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  /* - r1[I][E] * LT2_vv[E][A] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Airrep, 1, 1, "LT2_VV");
  global_dpd_->contract222(&RIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Airrep, 3, 3, "LT2_vv");
  global_dpd_->contract222(&Ria, &Int, &Dia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Int);

  /* + L2R1_ov[M][E] * t2[i][m][a][e] */

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 3, "L2R1_ov");
  global_dpd_->dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 3, "L2R1_ov");
  global_dpd_->dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&Int);
  global_dpd_->buf4_close(&T2);

  /* - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 0, 1, "L2R1_OV");
  global_dpd_->file2_init(&XIJ, PSIF_EOM_TMP, Tirrep, 0, 0, "XIJ");
  global_dpd_->contract222(&TIA, &Int, &XIJ, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&XIJ, PSIF_EOM_TMP, Tirrep, 0, 0, "XIJ");
  global_dpd_->contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&XIJ);

  global_dpd_->file2_init(&Int, PSIF_EOM_TMP, Tirrep, 2, 3, "L2R1_ov");
  global_dpd_->file2_init(&Xij, PSIF_EOM_TMP, Tirrep, 2, 2, "Xij");
  global_dpd_->contract222(&Tia, &Int, &Xij, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&Int);

  global_dpd_->file2_init(&Xij, PSIF_EOM_TMP, Tirrep, 2, 2, "Xij");
  global_dpd_->contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&Xij);

  global_dpd_->file2_close(&DIA);
  global_dpd_->file2_close(&Dia);

  global_dpd_->file2_close(&TIA);
  global_dpd_->file2_close(&Tia);
  global_dpd_->file2_close(&RIA);
  global_dpd_->file2_close(&Ria);
  global_dpd_->file2_close(&LIA);
  global_dpd_->file2_close(&Lia);

  return;
}

}} // namespace psi::ccdensity
