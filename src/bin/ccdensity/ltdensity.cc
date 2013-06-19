/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void ltdensity_rohf(struct TD_Params S);
void ltdensity_uhf(struct TD_Params S);
void ltdensity_intermediates(struct TD_Params S);
void sort_ltd_rohf(struct TD_Params S);
void sort_ltd_uhf(struct TD_Params S);

void ltdensity_rohf(struct TD_Params S)
{
  dpdfile2 DAI, Dai, DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, Int, XIJ, Xij, R1;
  dpdbuf4 T2, L2, R2, I2;
  dpdfile2 D, T1, L1, Z;

  if(S.irrep == 0) {
    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 0, "LTDIJ");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2); 
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2); 
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&D);

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 0, "LTDij");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&D);

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 1, 1, "LTDAB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&D);

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 1, 1, "LTDab");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&D);
  }

  /* R_I^A */
  dpd_->file2_init(&R1, PSIF_CC_GR, S.irrep, 0, 1, "RIA");
  dpd_->file2_copy(&R1, PSIF_CC_TMP, "LTDIA");
  dpd_->file2_close(&R1);

  if(S.irrep == 0) {
    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 1, "LTDIA");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&L1);
    dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&Z);

    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */

    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&Z);
    dpd_->file2_close(&T1);

    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */

    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(A,E)");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&Z);
    dpd_->file2_close(&D);
  }

  /* R_i^a */
  dpd_->file2_init(&R1, PSIF_CC_GR, S.irrep, 0, 1, "Ria");
  dpd_->file2_copy(&R1, PSIF_CC_TMP, "LTDia");
  dpd_->file2_close(&R1);

  if(S.irrep == 0) {
    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 1, "LTDia");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(i,m)");
    dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&L1);
    dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&Z);
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(i,m)");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&Z);
    dpd_->file2_close(&T1);
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(a,e)");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&Z);
    dpd_->file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->file2_copy(&L1, PSIF_CC_TMP, "LTDAI");
    dpd_->file2_close(&L1);

    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    dpd_->file2_copy(&L1, PSIF_CC_TMP, "LTDai");
    dpd_->file2_close(&L1);

    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->file2_scm(&L1, (1/S.R0));
    dpd_->file2_close(&L1);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    dpd_->file2_scm(&L1, (1/S.R0));
    dpd_->file2_close(&L1);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_scm(&L2, (1/S.R0));
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 7, 2, 7, 0, "Lijab");
    dpd_->buf4_scm(&L2, (1/S.R0));
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_scm(&L2, (1/S.R0));
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_->buf4_scm(&L2, (1/S.R0));
    dpd_->buf4_close(&L2);
  }

  ltdensity_intermediates(S);

  dpd_->file2_init(&TIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  dpd_->file2_init(&Tia, PSIF_CC_OEI, 0, 0, 1, "tia");
  dpd_->file2_init(&RIA, PSIF_CC_GR, S.irrep, 0, 1, "RIA");
  dpd_->file2_init(&Ria, PSIF_CC_GR, S.irrep, 0, 1, "Ria");
  dpd_->file2_init(&LIA, PSIF_CC_GLG, 0, 0, 1, "LIA");
  dpd_->file2_init(&Lia, PSIF_CC_GLG, 0, 0, 1, "Lia");

  /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */

  dpd_->file2_init(&DIJ, PSIF_CC_TMP, S.irrep, 0, 0, "LTDIJ");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_->file2_axpy(&Int, &DIJ, -1.0, 1);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->contract222(&TIA, &Int, &DIJ, 0, 0, -1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&DIJ);

  dpd_->file2_init(&Dij, PSIF_CC_TMP, S.irrep, 0, 0, "LTDij");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 0, "LR_oo");
  dpd_->file2_axpy(&Int, &Dij, -1.0, 1);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_->contract222(&Tia, &Int, &Dij, 0, 0, -1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&Dij);

  /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */

  dpd_->file2_init(&DAB, PSIF_CC_TMP, S.irrep, 1, 1, "LTDAB");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_->file2_axpy(&Int, &DAB, 1.0, 0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->contract222(&Int, &TIA, &DAB, 1, 1, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&DAB);

  dpd_->file2_init(&Dab, PSIF_CC_TMP, S.irrep, 1, 1, "LTDab");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 1, 1, "LR_vv");
  dpd_->file2_axpy(&Int, &Dab, 1.0, 0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_->contract222(&Int, &Tia, &Dab, 1, 1, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&Dab);

  /* D[a][i] = +L2R1_ov[i][a] */

  dpd_->file2_init(&DAI, PSIF_CC_TMP, S.irrep, 0, 1, "LTDAI");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->file2_axpy(&Int, &DAI, 1.0, 0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&DAI);

  dpd_->file2_init(&Dai, PSIF_CC_TMP, S.irrep, 0, 1, "LTDai");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_->file2_axpy(&Int, &Dai, 1.0, 0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&Dai);

  dpd_->file2_init(&DIA, PSIF_CC_TMP, S.irrep, 0, 1, "LTDIA");
  dpd_->file2_init(&Dia, PSIF_CC_TMP, S.irrep, 0, 1, "LTDia");

  /* D[i][a] = L1R2_ov[i][a] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L1R2_OV");
  dpd_->file2_axpy(&Int, &DIA, 1.0, 0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L1R2_ov");
  dpd_->file2_axpy(&Int, &Dia, 1.0, 0);
  dpd_->file2_close(&Int);

  /* - LR_OO[M][I] * t1[M][A] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_->contract222(&Int, &TIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 0, "LR_oo");
  dpd_->contract222(&Int, &Tia, &Dia, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  /* - t1[I][E] * LR_vv[E][A] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_->contract222(&TIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 1, 1, "LR_vv");
  dpd_->contract222(&Tia, &Int, &Dia, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  /* - LT2_OO[M][I] * r1[M][A] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, 0, 0, 0, "LT2_OO");
  dpd_->contract222(&Int, &RIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, 0, 0, 0, "LT2_oo");
  dpd_->contract222(&Int, &Ria, &Dia, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  /* - r1[I][E] * LT2_VV[E][A] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, 0, 1, 1, "LT2_VV");
  dpd_->contract222(&RIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, 0, 1, 1, "LT2_vv");
  dpd_->contract222(&Ria, &Int, &Dia, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  /* + L2R1_ov[M][E] * t2[i][m][a][e] */

  dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB"); 
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->buf4_close(&T2);

  dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb"); 
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_->dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->buf4_close(&T2);

  dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab"); 
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_->dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->buf4_close(&T2);

  dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB"); 
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->buf4_close(&T2);
    
  /* - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->file2_init(&XIJ, PSIF_EOM_TMP, S.irrep, 0, 0, "XIJ");
  dpd_->contract222(&TIA, &Int, &XIJ, 0, 0, 1.0, 0.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&XIJ, PSIF_EOM_TMP, S.irrep, 0, 0, "XIJ");
  dpd_->contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&XIJ);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_ov");
  dpd_->file2_init(&Xij, PSIF_EOM_TMP, S.irrep, 0, 0, "Xij");
  dpd_->contract222(&Tia, &Int, &Xij, 0, 0, 1.0, 0.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Xij, PSIF_EOM_TMP, S.irrep, 0, 0, "Xij");
  dpd_->contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Xij);

  dpd_->file2_close(&DIA);
  dpd_->file2_close(&Dia);

  dpd_->file2_close(&TIA);
  dpd_->file2_close(&Tia);
  dpd_->file2_close(&RIA);
  dpd_->file2_close(&Ria);
  dpd_->file2_close(&LIA);
  dpd_->file2_close(&Lia);

  return;
}

void ltdensity_uhf(struct TD_Params S)
{
  dpdfile2 DAI, Dai, DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, Int, XIJ, Xij, R1;
  dpdbuf4 T2, L2, R2, I2;
  dpdfile2 D, T1, L1, Z;

  if(S.irrep == 0 ) { /* Symmetric Transition */

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 0, "LTDIJ");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2); 
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2); 
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&D);

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 2, 2, "LTDij");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_->contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    dpd_->contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&D);

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 1, 1, "LTDAB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&D);

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 3, 3, "LTDab");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&D);
  }

  /* R_I^A */
    
  dpd_->file2_init(&R1, PSIF_CC_GR, S.irrep, 0, 1, "RIA");
  dpd_->file2_copy(&R1, PSIF_CC_TMP, "LTDIA");
  dpd_->file2_close(&R1); 

  if(S.irrep == 0) { // Symmetric Transitions

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 0, 1, "LTDIA");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&L1);
    dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&Z);

    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */

    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&Z);
    dpd_->file2_close(&T1);
    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 1, 1, "Z(A,E)");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&Z);
    dpd_->file2_close(&D);
  }

  /* R_i^a */
    
  dpd_->file2_init(&R1, PSIF_CC_GR, S.irrep, 2, 3, "Ria");
  dpd_->file2_copy(&R1, PSIF_CC_TMP, "LTDia"); 
  dpd_->file2_close(&R1);

  if(S.irrep == 0) { /* Symmetric Transitions */

    dpd_->file2_init(&D, PSIF_CC_TMP, 0, 2, 3, "LTDia");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&L1);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 2, 2, "Z(i,m)");
    dpd_->contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&L1);
    dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&Z);
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 2, 2, "Z(i,m)");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_->contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&T2);
    dpd_->buf4_close(&L2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&Z);
    dpd_->file2_close(&T1);
    dpd_->file2_init(&Z, PSIF_CC_TMP0, 0, 3, 3, "Z(a,e)");
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_->contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&T2);
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->file2_close(&Z);
    dpd_->file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->file2_copy(&L1, PSIF_CC_TMP, "LTDAI");
    dpd_->file2_close(&L1);

    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    dpd_->file2_copy(&L1, PSIF_CC_TMP, "LTDai");
    dpd_->file2_close(&L1);

    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    dpd_->file2_scm(&L1, (1/S.R0));
    dpd_->file2_close(&L1);
    dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    dpd_->file2_scm(&L1, (1/S.R0));
    dpd_->file2_close(&L1);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_scm(&L2, (1/S.R0));
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_->buf4_scm(&L2, (1/S.R0));
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->buf4_scm(&L2, (1/S.R0));
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_->buf4_scm(&L2, (1/S.R0));
    dpd_->buf4_close(&L2);
  }

  ltdensity_intermediates(S);

  dpd_->file2_init(&TIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  dpd_->file2_init(&Tia, PSIF_CC_OEI, 0, 2, 3, "tia");
  dpd_->file2_init(&RIA, PSIF_CC_GR, S.irrep, 0, 1, "RIA");
  dpd_->file2_init(&Ria, PSIF_CC_GR, S.irrep, 2, 3, "Ria");
  dpd_->file2_init(&LIA, PSIF_CC_GL, S.irrep, 0, 1, "LIA");
  dpd_->file2_init(&Lia, PSIF_CC_GL, S.irrep, 2, 3, "Lia");

  /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */

  dpd_->file2_init(&DIJ, PSIF_CC_TMP, S.irrep, 0, 0, "LTDIJ");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_->file2_axpy(&Int, &DIJ, -1.0, 1);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->contract222(&TIA, &Int, &DIJ, 0, 0, -1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&DIJ);

  dpd_->file2_init(&Dij, PSIF_CC_TMP, S.irrep, 2, 2, "LTDij");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 2, "LR_oo");
  dpd_->file2_axpy(&Int, &Dij, -1.0, 1);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_->contract222(&Tia, &Int, &Dij, 0, 0, -1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&Dij);

  /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */

  dpd_->file2_init(&DAB, PSIF_CC_TMP, S.irrep, 1, 1, "LTDAB");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_->file2_axpy(&Int, &DAB, 1.0, 0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->contract222(&Int, &TIA, &DAB, 1, 1, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&DAB);

  dpd_->file2_init(&Dab, PSIF_CC_TMP, S.irrep, 3, 3, "LTDab");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 3, 3, "LR_vv");
  dpd_->file2_axpy(&Int, &Dab, 1.0, 0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_->contract222(&Int, &Tia, &Dab, 1, 1, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&Dab);

  /* D[a][i] = +L2R1_ov[i][a] */

  dpd_->file2_init(&DAI, PSIF_CC_TMP, S.irrep, 0, 1, "LTDAI");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->file2_axpy(&Int, &DAI, 1.0, 0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&DAI);

  dpd_->file2_init(&Dai, PSIF_CC_TMP, S.irrep, 2, 3, "LTDai");
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_->file2_axpy(&Int, &Dai, 1.0, 0);
  dpd_->file2_close(&Int);
  dpd_->file2_close(&Dai);

  dpd_->file2_init(&DIA, PSIF_CC_TMP, S.irrep, 0, 1, "LTDIA");

  dpd_->file2_init(&Dia, PSIF_CC_TMP, S.irrep, 2, 3, "LTDia");

  /* D[i][a] = L1R2_ov[i][a] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L1R2_OV");
  dpd_->file2_axpy(&Int, &DIA, 1.0, 0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 3, "L1R2_ov");
  dpd_->file2_axpy(&Int, &Dia, 1.0, 0);
  dpd_->file2_close(&Int);

  /* - LR_OO[M][I] * t1[M][A] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 0, "LR_OO");
  dpd_->contract222(&Int, &TIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 2, "LR_oo");
  dpd_->contract222(&Int, &Tia, &Dia, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  /* - t1[I][E] * LR_vv[E][A] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 1, 1, "LR_VV");
  dpd_->contract222(&TIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 3, 3, "LR_vv");
  dpd_->contract222(&Tia, &Int, &Dia, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  /* - LT2_OO[M][I] * r1[M][A] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, 0, 0, 0, "LT2_OO");
  dpd_->contract222(&Int, &RIA, &DIA, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, 0, 2, 2, "LT2_oo");
  dpd_->contract222(&Int, &Ria, &Dia, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  /* - r1[I][E] * LT2_vv[E][A] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, 0, 1, 1, "LT2_VV");
  dpd_->contract222(&RIA, &Int, &DIA, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, 0, 3, 3, "LT2_vv");
  dpd_->contract222(&Ria, &Int, &Dia, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Int);

  /* + L2R1_ov[M][E] * t2[i][m][a][e] */

  dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB"); 
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->buf4_close(&T2);

  dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb"); 
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_->dot24(&Int, &T2, &DIA, 0, 0, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->buf4_close(&T2);

  dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab"); 
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_->dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->buf4_close(&T2);

  dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB"); 
  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->dot24(&Int, &T2, &Dia, 0, 0, 1.0, 1.0);
  dpd_->file2_close(&Int);
  dpd_->buf4_close(&T2);
    
  /* - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 0, 1, "L2R1_OV");
  dpd_->file2_init(&XIJ, PSIF_EOM_TMP, S.irrep, 0, 0, "XIJ");
  dpd_->contract222(&TIA, &Int, &XIJ, 0, 0, 1.0, 0.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&XIJ, PSIF_EOM_TMP, S.irrep, 0, 0, "XIJ");
  dpd_->contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&XIJ);

  dpd_->file2_init(&Int, PSIF_EOM_TMP, S.irrep, 2, 3, "L2R1_ov");
  dpd_->file2_init(&Xij, PSIF_EOM_TMP, S.irrep, 2, 2, "Xij");
  dpd_->contract222(&Tia, &Int, &Xij, 0, 0, 1.0, 0.0);
  dpd_->file2_close(&Int);

  dpd_->file2_init(&Xij, PSIF_EOM_TMP, S.irrep, 2, 2, "Xij");
  dpd_->contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&Xij);

  dpd_->file2_close(&DIA);
  dpd_->file2_close(&Dia);

  dpd_->file2_close(&TIA);
  dpd_->file2_close(&Tia);
  dpd_->file2_close(&RIA);
  dpd_->file2_close(&Ria);
  dpd_->file2_close(&LIA);
  dpd_->file2_close(&Lia);

  return;
}

}} // namespace psi::ccdensity
