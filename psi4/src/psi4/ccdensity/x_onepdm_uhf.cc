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

/* onepdm_uhf(): Computes the non-R0 parts of the unrelaxed EOM 1pdm
* intermediates are defined in x_oe_intermediates.c
*
*   D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f]
*
*   D[a][b] = +LR_vv[a][b] + t1[n][b] * L2R1_ov[n][a]
*
*   D[a][i] = +L2R1_ov[i][a]
*
*   D[i][a] = + L1R2_ov[i][a]
*            - t1[m][a] * LR_oo[m][i]
*            - t1[i][e] * LR_vv[e][a]
*            - r1[m][a] * LT2_oo[m][i]
*            - r1[i][e] * LT2_vv[e][a]
*            + L2R1_ov[M][E] * (t2[i][m][a][e] - t1[i][e] * t1[m][a])
*
* RAK 2003
*/

void x_onepdm_uhf(struct RHO_Params rho_params)
{
  dpdfile2 DAI, Dai, DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, I, XIJ, Xij;
  dpdbuf4 T2, L2, R2, I2;
  int L_irr, R_irr, G_irr;
  double dot_IA, dot_ia, dot_AI, dot_ai;
  L_irr = rho_params.L_irr;
  R_irr = rho_params.R_irr;
  G_irr = rho_params.G_irr;

  global_dpd_->file2_init(&TIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&Tia, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_init(&RIA, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&Ria, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->file2_init(&LIA, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  global_dpd_->file2_init(&Lia, PSIF_CC_GL, L_irr, 2, 3, "Lia");

  /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */
  global_dpd_->file2_init(&DIJ, PSIF_CC_OEI, G_irr, 0, 0, rho_params.DIJ_lbl);
  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->file2_axpy(&I, &DIJ, -1.0, 1);
  global_dpd_->file2_close(&I);

  if (!params.connect_xi) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract222(&TIA, &I, &DIJ, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&I);
  }
  global_dpd_->file2_close(&DIJ);

  global_dpd_->file2_init(&Dij, PSIF_CC_OEI, G_irr, 2, 2, rho_params.Dij_lbl);
  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 2, "LR_oo");
  global_dpd_->file2_axpy(&I, &Dij, -1.0, 1);
  global_dpd_->file2_close(&I);

  if (!params.connect_xi) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->contract222(&Tia, &I, &Dij, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&I);
  }
  global_dpd_->file2_close(&Dij);

  /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */
  global_dpd_->file2_init(&DAB, PSIF_CC_OEI, G_irr, 1, 1, rho_params.DAB_lbl);
  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 1, 1, "LR_VV");
  global_dpd_->file2_axpy(&I, &DAB, 1.0, 0);
  global_dpd_->file2_close(&I);

  if (!params.connect_xi) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract222(&I, &TIA, &DAB, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&I);
  }
  global_dpd_->file2_close(&DAB);

  global_dpd_->file2_init(&Dab, PSIF_CC_OEI, G_irr, 3, 3, rho_params.Dab_lbl);
  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 3, 3, "LR_vv");
  global_dpd_->file2_axpy(&I, &Dab, 1.0, 0);
  global_dpd_->file2_close(&I);

  if (!params.connect_xi) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->contract222(&I, &Tia, &Dab, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&I);
  }
  global_dpd_->file2_close(&Dab);

  /* D[a][i] = +L2R1_ov[i][a] */

  if (!params.connect_xi) {
    global_dpd_->file2_init(&DAI, PSIF_CC_OEI, G_irr, 0, 1, rho_params.DAI_lbl);
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_axpy(&I, &DAI, 1.0, 0);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_close(&DAI);

    global_dpd_->file2_init(&Dai, PSIF_CC_OEI, G_irr, 2, 3, rho_params.Dai_lbl);
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->file2_axpy(&I, &Dai, 1.0, 0);
    global_dpd_->file2_close(&I);
    global_dpd_->file2_close(&Dai);
  }

  /*
     D[I][A] = (1-R0)*tIA + L1R2_OV[I][A]
               - LR_OO[M][I] * t1[M][A]
               - t1[I][E] * LR_vv[E][A]
               - LT2_OO[M][I] * r1[M][A]
               - r1[I][E] * LT2_vv[E][A]
               + L2R1_ov[M][E] * (t2[i][m][a][e] - t1[i][e] * t1[m][a])
  */

  /* (1-R0) * tIA */
  global_dpd_->file2_init(&DIA, PSIF_CC_OEI, G_irr, 0, 1, rho_params.DIA_lbl);

  if (G_irr == 0) {
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_axpy(&I, &DIA, 1.0, 0);
    global_dpd_->file2_close(&I);
  }

  global_dpd_->file2_init(&Dia, PSIF_CC_OEI, G_irr, 2, 3, rho_params.Dia_lbl);
  if (G_irr == 0) {
    global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_axpy(&I, &Dia, 1.0, 0);
    global_dpd_->file2_close(&I);
  }

  /* D[i][a] = L1R2_ov[i][a] */
  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  global_dpd_->file2_axpy(&I, &DIA, 1.0, 0);
  global_dpd_->file2_close(&I);

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 3, "L1R2_ov");
  global_dpd_->file2_axpy(&I, &Dia, 1.0, 0);
  global_dpd_->file2_close(&I);

  /* - LR_OO[M][I] * t1[M][A] */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->contract222(&I, &TIA, &DIA, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I);

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 2, "LR_oo");
  global_dpd_->contract222(&I, &Tia, &Dia, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I);

  /* - t1[I][E] * LR_vv[E][A] */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 1, 1, "LR_VV");
  global_dpd_->contract222(&TIA, &I, &DIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I);

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 3, 3, "LR_vv");
  global_dpd_->contract222(&Tia, &I, &Dia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I);

  /* - LT2_OO[M][I] * r1[M][A] */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, L_irr, 0, 0, "LT2_OO");
  global_dpd_->contract222(&I, &RIA, &DIA, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I);

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, L_irr, 2, 2, "LT2_oo");
  global_dpd_->contract222(&I, &Ria, &Dia, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I);

  /* - r1[I][E] * LT2_vv[E][A] */

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, L_irr, 1, 1, "LT2_VV");
  global_dpd_->contract222(&RIA, &I, &DIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I);

  global_dpd_->file2_init(&I, PSIF_EOM_TMP, L_irr, 3, 3, "LT2_vv");
  global_dpd_->contract222(&Ria, &I, &Dia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I);

  /* term 6, + L2R1_ov[M][E] * t2[i][m][a][e] */

  if (!params.connect_xi) {
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->dot24(&I, &T2, &DIA, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&I);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->dot24(&I, &T2, &DIA, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&I);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->dot24(&I, &T2, &Dia, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&I);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->dot24(&I, &T2, &Dia, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&I);
    global_dpd_->buf4_close(&T2);
  }

  /* term 7, - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

  if (!params.connect_xi) {
    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&XIJ, PSIF_EOM_TMP, G_irr, 0, 0, "XIJ");
    global_dpd_->contract222(&TIA, &I, &XIJ, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&XIJ, PSIF_EOM_TMP, G_irr, 0, 0, "XIJ");
    global_dpd_->contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&XIJ);

    global_dpd_->file2_init(&I, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->file2_init(&Xij, PSIF_EOM_TMP, G_irr, 2, 2, "Xij");
    global_dpd_->contract222(&Tia, &I, &Xij, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&I);

    global_dpd_->file2_init(&Xij, PSIF_EOM_TMP, G_irr, 2, 2, "Xij");
    global_dpd_->contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&Xij);
  }

  global_dpd_->file2_close(&DIA);
  global_dpd_->file2_close(&Dia);

  global_dpd_->file2_close(&TIA);
  global_dpd_->file2_close(&Tia);
  global_dpd_->file2_close(&RIA);
  global_dpd_->file2_close(&Ria);
  global_dpd_->file2_close(&LIA);
  global_dpd_->file2_close(&Lia);

  /* compute overlaps */
  global_dpd_->file2_init(&DIA, PSIF_CC_OEI, G_irr, 0, 1, rho_params.DIA_lbl);
  dot_IA = global_dpd_->file2_dot_self(&DIA);
  global_dpd_->file2_close(&DIA);
  global_dpd_->file2_init(&Dia, PSIF_CC_OEI, G_irr, 2, 3, rho_params.Dia_lbl);
  dot_ia = global_dpd_->file2_dot_self(&Dia);
  global_dpd_->file2_close(&Dia);
  global_dpd_->file2_init(&DAI, PSIF_CC_OEI, G_irr, 0, 1, rho_params.DAI_lbl);
  dot_AI = global_dpd_->file2_dot_self(&DAI);
  global_dpd_->file2_close(&DAI);
  global_dpd_->file2_init(&Dai, PSIF_CC_OEI, G_irr, 2, 3, rho_params.Dai_lbl);
  dot_ai = global_dpd_->file2_dot_self(&Dai);
  global_dpd_->file2_close(&Dai);
    outfile->Printf("\tOverlaps of onepdm after excited-state parts added.\n");
    outfile->Printf("\t<DIA|DIA> = %15.10lf     <Dia|Dia> = %15.10lf\n", dot_IA, dot_ia);
    outfile->Printf("\t<DAI|DAI> = %15.10lf     <Dai|Dai> = %15.10lf\n", dot_AI, dot_ai);
    outfile->Printf("\t<Dpq|Dqp> = %15.10lf\n", dot_IA+dot_ia+dot_AI+dot_ai);
  return;
}

}} // namespace psi::ccdensity
