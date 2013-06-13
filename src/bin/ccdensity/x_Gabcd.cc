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

/* x_Gabcd(): computes non-R0 parts of Gijkl 2pdm */
/* Gabcd = 0.5 Lmnab * Rmncd + P(cd) Lmnab * Rmc * Tnd */
/* Gabcd = 0.5 Lmnab * Rmncd + Lmnab (Rmc * Tnd + tmc * Rnd) */

void x_Gabcd(void)
{
  dpdfile2 R1, T1;
  dpdbuf4 L2, I2, I3, R2, GABCD, Gabcd, GAbCd;
  int L_irr, R_irr, G_irr;
  double value;
  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

  if (params.ref == 0 || params.ref == 1) {
    /* GABCD += 0.5 * Lmnab * Rmncd */
    dpd_->buf4_init(&GABCD, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_->contract444(&L2, &R2, &GABCD, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&GABCD);

    dpd_->buf4_init(&Gabcd, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
    dpd_->contract444(&L2, &R2, &Gabcd, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&Gabcd);

    dpd_->buf4_init(&GAbCd, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_->contract444(&L2, &R2, &GAbCd, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&GAbCd);
  }
  else {
    /* GABCD += 0.5 * Lmnab * Rmncd */
    dpd_->buf4_init(&GABCD, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_->contract444(&L2, &R2, &GABCD, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&GABCD);

    dpd_->buf4_init(&Gabcd, PSIF_CC_GAMMA, G_irr, 17, 17, 17, 17, 0, "Gabcd");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 12, 17, 12, 17, 0, "Rijab");
    dpd_->contract444(&L2, &R2, &Gabcd, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&Gabcd);

    dpd_->buf4_init(&GAbCd, PSIF_CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_->contract444(&L2, &R2, &GAbCd, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&L2);
    dpd_->buf4_close(&GAbCd);
  }

  if (params.ref == 0 || params.ref == 1) {
    /* GABCD = LMNAB (RMC * TND + tMC * RND) */
    dpd_->buf4_init(&GABCD, PSIF_CC_GAMMA, G_irr, 7, 5, 7, 7, 0, "GABCD");
    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 7, 11, 7, 11, 0, "L2R1_VVOV(pqsr)");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract424(&I2, &T1, &GABCD, 3, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&I2);
  
    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_VVOV");
    dpd_->contract244(&T1, &I2, &GABCD, 0, 2, 1, 1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->buf4_close(&I2);
    dpd_->buf4_close(&GABCD);

    /* Gabcd = Lmnab (Rmc * Tnd + tmc * Rnd) */
    dpd_->buf4_init(&Gabcd, PSIF_CC_GAMMA, G_irr, 7, 5, 7, 7, 0, "Gabcd");
    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 7, 11, 7, 11, 0, "L2R1_vvov(pqsr)");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->contract424(&I2, &T1, &Gabcd, 3, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&I2);

    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_vvov");
    dpd_->contract244(&T1, &I2, &Gabcd, 0, 2, 1, 1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->buf4_close(&I2);
    dpd_->buf4_close(&Gabcd);

    /* GAbCd = LMnAb (RMC * Tnd + TMC * Rnd) */
    dpd_->buf4_init(&GAbCd, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 5, 11, 5, 11, 0, "L2R1_VvoV(pqsr)");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->contract424(&I2, &T1, &GAbCd, 3, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->buf4_close(&I2);

    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvOv");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract244(&T1, &I2, &GAbCd, 0, 2, 1, 1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->buf4_close(&I2);
    dpd_->buf4_close(&GAbCd);
  }
  else {
    /* GABCD = LMNAB (RMC * TND + tMC * RND) */
    dpd_->buf4_init(&GABCD, PSIF_CC_GAMMA, G_irr, 7, 5, 7, 7, 0, "GABCD");
    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 7, 21, 7, 21, 0, "L2R1_VVOV(pqsr)");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract424(&I2, &T1, &GABCD, 3, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&I2);
  
    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 7, 20, 7, 20, 0, "L2R1_VVOV");
    dpd_->contract244(&T1, &I2, &GABCD, 0, 2, 1, 1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->buf4_close(&I2);
    dpd_->buf4_close(&GABCD);

    /* Gabcd = Lmnab (Rmc * Tnd + tmc * Rnd) */
    dpd_->buf4_init(&Gabcd, PSIF_CC_GAMMA, G_irr, 17, 15, 17, 17, 0, "Gabcd");
    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 17, 31, 17, 31, 0, "L2R1_vvov(pqsr)");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract424(&I2, &T1, &Gabcd, 3, 0, 0, -1.0, 1.0);
    dpd_->buf4_close(&I2);

    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 17, 30, 17, 30, 0, "L2R1_vvov");
    dpd_->contract244(&T1, &I2, &Gabcd, 0, 2, 1, 1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->buf4_close(&I2);
    dpd_->buf4_close(&Gabcd);

    /* GAbCd = LMnAb (RMC * Tnd + TMC * Rnd) */
    dpd_->buf4_init(&GAbCd, PSIF_CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 28, 26, 28, 26, 0, "L2R1_VvoV(pqsr)");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->contract424(&I2, &T1, &GAbCd, 3, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->buf4_close(&I2);

    dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 28, 24, 28, 24, 0, "L2R1_VvOv");
    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->contract244(&T1, &I2, &GAbCd, 0, 2, 1, 1.0, 1.0);
    dpd_->file2_close(&T1);
    dpd_->buf4_close(&I2);
    dpd_->buf4_close(&GAbCd);
  }

  /* bra-ket symmetrize */
  if (params.ref == 0 || params.ref == 1) {
    dpd_->buf4_init(&GABCD, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_->buf4_symm(&GABCD);
    dpd_->buf4_close(&GABCD);
    dpd_->buf4_init(&Gabcd, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
    dpd_->buf4_symm(&Gabcd);
    dpd_->buf4_close(&Gabcd);
    dpd_->buf4_init(&GAbCd, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_->buf4_symm(&GAbCd);
    dpd_->buf4_close(&GAbCd);
  }
  else {
    dpd_->buf4_init(&GABCD, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_->buf4_symm(&GABCD);
    dpd_->buf4_close(&GABCD);
    dpd_->buf4_init(&Gabcd, PSIF_CC_GAMMA, G_irr, 17, 17, 17, 17, 0, "Gabcd");
    dpd_->buf4_symm(&Gabcd);
    dpd_->buf4_close(&Gabcd);
    dpd_->buf4_init(&GAbCd, PSIF_CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    dpd_->buf4_symm(&GAbCd);
    dpd_->buf4_close(&GAbCd);
  }
  return;
}

}} // namespace psi::ccdensity
