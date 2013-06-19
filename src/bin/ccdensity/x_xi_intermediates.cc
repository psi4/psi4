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
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* build one-electron intermediates for the construction of the xi amplitudes */

void x_xi_intermediates(void)
{
  dpdfile2 L1, R1, T1, I, LR1, LR2, LT1, LT2;
  dpdbuf4 L2, T2, R2, D, H2, I2;
  int L_irr, R_irr, G_irr;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

  if (params.ref == 0) {

    /* RD_OO(I,J) = (2RImEf - RImFe) <jm|ef> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);

    /* RD_VV(E,A) = (2RMnFe - RMnEf) <mn|fa> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);

    /* R2Wamef_OV(n,g) = Rnmef Wgmef */
    dpd_->file2_init(&I, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);

    /* R2Wamef_OOOV = Rmnfg Weifg */
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(Mn,eI)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(Mn,eI)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 0, 10, "R2Wamef_OoOv");
    dpd_->buf4_close(&I2);

    /* R1Wamef_VV (a,e) = Rmf Wamef */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->file2_close(&I);

    /* R1Wmnie_OO (M,I) = Rne Wmnie */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->file2_close(&I);

  }
  else if (params.ref == 1) {
    /* RD_OO(I,J)  =  RIMEF * <JM||EF> + RImEf <Jm|Ef> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);
  
    /* RD_oo(i,j)  =  Rimef * <jm||ef> + RiMeF <jM|eF> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);

    /* RD_VV(E,A)  =  RMNFE * <MN||FA> + RmNfE <mN|fA> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);
  
    /* RD_vv(e,a)  =  Rmnfe * <mn||fa> + RMnFe <Mn|Fa> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_vv");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "Rijab");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);
  
    /* R2Wamef_OOVO = Rmnfg Weifg -> OOOV (mn,ie) */
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 2, 11, 2, 11, 0, "Z(M>N,EI)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 2, 11, 2, 11, 0, "Z(M>N,EI)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 2, 10, "R2Wamef_OOOV"); 
    dpd_->buf4_close(&I2);

    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(Mn,eI)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(Mn,eI)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 0, 10, "R2Wamef_OoOv");
    dpd_->buf4_close(&I2);
  
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(mN,Ei)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJAb");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 0, 11, 0, 11, 0, "Z(mN,Ei)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 0, 10, "R2Wamef_oOoV");
    dpd_->buf4_close(&I2);

    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 2, 11, 2, 11, 0, "Z(m>n,ei)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 2, 11, 2, 11, 0, "Z(m>n,ei)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 2, 10, "R2Wamef_ooov");
    dpd_->buf4_close(&I2);
  
    /* R2Wamef_OV(n,g) = Rnmef Wgmef */
    dpd_->file2_init(&I, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);
  
    dpd_->file2_init(&I, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_ov");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);

    /* R1Wmnie_OO (M,I) = Rne Wmnie */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 2, 10, 0, "WMNIE");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->file2_close(&I);

    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_oo");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 2, 10, 0, "Wmnie");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->file2_close(&I);

    /* R1Wamef_VV (a,e) = Rmf Wamef */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->file2_close(&I);

    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_vv");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
  
    dpd_->file2_close(&I);
  }
  else { /* UHF */
    /* RD_OO(I,J)  =  RIMEF * <JM||EF> + RImEf <Jm|Ef> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);
  
    /* RD_oo(i,j)  =  Rimef * <jm||ef> + RiMeF <jM|eF> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 2, 2, "RD_oo");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 10, 17, 12, 17, 0, "Rijab");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_->contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);

    /* RD_VV(E,A)  =  RMNFE * <MN||FA> + RmNfE <mN|fA> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);
  
    /* RD_vv(e,a)  =  Rmnfe * <mn||fa> + RMnFe <Mn|Fa> */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 3, 3, "RD_vv");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 12, 15, 12, 17, 0, "Rijab");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&I);
  
    /* R2Wamef_OOVO = Rmnfg W(eifg) -> Z(mn,ei) -> OOOV (mn,ie) */
    /* make using (AM,EF) so we can avoid making the (MA,EF) at all*/
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 2, 21, 2, 21, 0, "Z(M>N,EI)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 2, 21, 2, 21, 0, "Z(M>N,EI)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 2, 20, "R2Wamef_OOOV"); 
    dpd_->buf4_close(&I2);
  
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 22, 25, 22, 25, 0, "Z(Mn,eI)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 22, 29, 22, 29, 0, "RIjaB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 22, 25, 22, 25, 0, "Z(Mn,eI)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 22, 24, "R2Wamef_OoOv");
    dpd_->buf4_close(&I2);
  
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 23, 26, 23, 26, 0, "Z(mN,Ei)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 23, 28, 23, 28, 0, "RiJAb");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 23, 26, 23, 26, 0, "Z(mN,Ei)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 23, 27, "R2Wamef_oOoV");
    dpd_->buf4_close(&I2);
  
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 12, 31, 12, 31, 0, "Z(m>n,ei)");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 12, 17, 12, 17, 0, "Rijab");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    dpd_->contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_close(&I2);
    dpd_->buf4_init(&I2, PSIF_EOM_TMP1, R_irr, 12, 31, 12, 31, 0, "Z(m>n,ei)");
    dpd_->buf4_sort(&I2, PSIF_EOM_TMP, pqsr, 12, 30, "R2Wamef_ooov");
    dpd_->buf4_close(&I2);

    psio_close(PSIF_EOM_TMP1,0);
    psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);
  
    /* R2Wamef_OV(n,g) = Rnmef Wgmef */
    dpd_->file2_init(&I, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
  
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
  
    dpd_->file2_close(&I);
  
    dpd_->file2_init(&I, PSIF_EOM_TMP, R_irr, 2, 3, "R2Wamef_ov");
  
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 10, 17, 12, 17, 0, "Rijab");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_->contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&H2);
    dpd_->buf4_close(&R2);

    dpd_->file2_close(&I);

    /* R1Wmnie_OO (M,I) = Rne Wmnie */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");

    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 20, 2, 20, 0, "WMNIE");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 22, 24, 22, 24, 0, "WMnIe");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);

    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 2, 2, "R1Wmnie_oo");

    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 10, 30, 12, 30, 0, "Wmnie");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 23, 27, 23, 27, 0, "WmNiE");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);

    dpd_->file2_close(&I);

    /* R1Wamef_VV (a,e) = Rmf Wamef */
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");

    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 21, 5, 21, 7, 0, "WAMEF");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);

    dpd_->file2_close(&I);
    dpd_->file2_init(&I, PSIF_EOM_TMP_XI, R_irr, 3, 3, "R1Wamef_vv");

    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 31, 15, 31, 17, 0, "Wamef");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
    dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    dpd_->dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
    dpd_->file2_close(&R1);
    dpd_->buf4_close(&H2);
  
    dpd_->file2_close(&I);
  }
  return;
}

}} // namespace psi::ccdensity
