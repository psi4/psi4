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

extern void x_xi_check(char *term_lbl);
extern void x_xi1_connected(void);

void x_xi1_uhf(void)
{
  dpdfile2 L1, XIA, Xia, I1, R1, F1, Z1A, Z1B;
  int L_irr, R_irr, G_irr;
  dpdbuf4 D, R2, L2, H2, I2;
  double tval;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

#ifdef DEBUG_XI
x_xi_check("begin xi1");
#endif
  /* term 1, XIA += 0.25 LIA Rmnef <mn||ef> */
  if ((R_irr == 0)  && (!params.connect_xi)) {
    dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    params.RD_overlap = 0.5 * dpd_->file2_trace(&I1);
    dpd_->file2_close(&I1);
    dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 2, 2, "RD_oo");
    params.RD_overlap += 0.5 * dpd_->file2_trace(&I1);
    dpd_->file2_close(&I1);
   /*fprintf(outfile,"RD overlap %15.10lf\n", params.RD_overlap);*/

    dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    dpd_->file2_copy(&L1, PSIF_EOM_XI, "XIA");
    dpd_->file2_close(&L1);
    dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    dpd_->file2_scm(&XIA, params.RD_overlap);
    dpd_->file2_close(&XIA);
    dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
    dpd_->file2_copy(&L1, PSIF_EOM_XI, "Xia");
    dpd_->file2_close(&L1);
    dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");
    dpd_->file2_scm(&Xia, params.RD_overlap);
    dpd_->file2_close(&Xia);
  }
  /* unnecessary 
  else {
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_scm(&XIA, 0.0);
    dpd_file2_close(&XIA);
    dpd_file2_init(&Xia, EOM_XI, G_irr, 2, 3, "Xia");
    dpd_file2_scm(&Xia, 0.0);
    dpd_file2_close(&Xia);
  }
  */
#ifdef DEBUG_XI
x_xi_check("term 1");
#endif

  /* term 2, Xia -= (Rmnef <in||ef>) * Lma */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  dpd_->contract222(&I1, &L1, &XIA, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&L1);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 2, 2, "RD_oo");
  dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
  dpd_->contract222(&I1, &L1, &Xia, 1, 1, -1.0, 1.0);
  dpd_->file2_close(&L1);
  dpd_->file2_close(&I1);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 2");
#endif

  /* term 3, XIA -= 0.5 LIE (Rmnfe <mn||fa>) */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
  dpd_->contract222(&L1, &I1, &XIA, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&I1);
  dpd_->file2_close(&L1);

  dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
  dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 3, 3, "RD_vv");
  dpd_->contract222(&L1, &I1, &Xia, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&I1);
  dpd_->file2_close(&L1);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 3");
#endif

  /* term 4, XIA += (Lme Rmnef) <in||af> */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  /* = OV(N,f) <IN||AF> + ov(nf) <In||Af> */
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  dpd_->dot24(&I1, &D, &XIA, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&D);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L1R2_ov");
  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_->dot24(&I1, &D, &XIA, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&D);
  dpd_->file2_close(&I1);

  /* = ov(n,f) <in||af> + OV(NF) <iN||aF> */
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L1R2_ov");
  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  dpd_->dot24(&I1, &D, &Xia, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&D);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
  dpd_->dot24(&I1, &D, &Xia, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&D);
  dpd_->file2_close(&I1);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 4");
#endif

  /* XIA += (Lmnef * Rmnef) FIA */ 
  dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  params.overlap1 = dpd_->file2_dot(&R1, &L1);
  dpd_->file2_close(&R1);
  dpd_->file2_close(&L1);
  dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 2, 3, "Lia");
  params.overlap1 += tval = dpd_->file2_dot(&R1, &L1);
  dpd_->file2_close(&R1);
  dpd_->file2_close(&L1);
  params.overlap2 = 1.0e0 - params.overlap1;

  /* explanation in xi1_connected and ROHF code */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  dpd_->file2_axpy(&F1, &XIA, params.overlap2, 0);
  dpd_->file2_close(&F1);
  dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 2, 3, "Fme");
  dpd_->file2_axpy(&F1, &Xia, params.overlap2, 0);
  dpd_->file2_close(&F1);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 5");
#endif

    /* term 6, -0.5 (Linef Rmnef) Fma */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_OO");
  dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  dpd_->contract222(&I1, &F1, &XIA, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&F1);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 2, "LR2_oo");
  dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 2, 3, "Fme");
  dpd_->contract222(&I1, &F1, &Xia, 0, 1, -1.0, 1.0);
  dpd_->file2_close(&F1);
  dpd_->file2_close(&I1);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 6");
#endif
  /* term 7, -0.5 (Lmnaf Rmnef) Fie */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_VV");
  dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  dpd_->contract222(&F1, &I1, &XIA, 0, 0, -1.0, 1.0);
  dpd_->file2_close(&F1);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 3, 3, "LR2_vv");
  dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 2, 3, "Fme");
  dpd_->contract222(&F1, &I1, &Xia, 0, 0, -1.0, 1.0);
  dpd_->file2_close(&F1);
  dpd_->file2_close(&I1);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

  if (!params.connect_xi) {
    /* term 8, (Fme Rmnef) Linaf) */
    dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(N,F)");
    dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 2, 7, 0, "RIJAB");
    dpd_->dot13(&F1, &R2, &I1, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&F1);
    dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 2, 3, "Fme");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
    dpd_->dot13(&F1, &R2, &I1, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&F1);
    dpd_->file2_close(&I1);

    dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 2, 3, "Z(n,f)");
    dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 2, 3, "Fme");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 10, 15, 12, 17, 0, "Rijab");
    dpd_->dot13(&F1, &R2, &I1, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&F1);
    dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_->dot13(&F1, &R2, &I1, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&R2);
    dpd_->file2_close(&F1);
    dpd_->file2_close(&I1);

    dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

    dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(N,F)");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_->dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_->dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->file2_close(&I1);

    dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 2, 3, "Z(n,f)");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 15, 12, 17, 0, "Lijab");
    dpd_->dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->file2_close(&I1);

    dpd_->file2_close(&XIA);
    dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 8");
#endif
  }

  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_OO");
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 21, 2, 21, 0, "WMNIE (M>N,EI)");
  dpd_->dot24(&I1, &H2, &XIA, 1, 0, 1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 2, "LR2_oo");
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
  dpd_->dot14(&I1, &H2, &XIA, 1, 0, -1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->file2_close(&I1);

  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 2, "LR2_oo");
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 10, 31, 12, 31, 0, "Wmnie (m>n,ei)");
  dpd_->dot24(&I1, &H2, &Xia, 1, 0, 1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_OO");
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
  dpd_->dot14(&I1, &H2, &Xia, 1, 0, -1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->file2_close(&I1);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 9");
#endif

/*  term 11 XIA -= Rmnef Lmoea Winof */
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  /* this would be easier if it would work but 13 and 31 shifts are
     incompatible when symmetry is on
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_contract442(&H2, &I2, &XIA, 0, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  */
  /* if I could do a 442(0,3) I could avoid these sorts */
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 21, 2, 21, 0, "WMNIE (M>N,EI)");
  dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, qrsp, 20, 0, "W (NF,OI)");
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 10, 31, 12, 31, 0, "Wmnie (m>n,ei)");
  dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, qrsp, 30, 10, "W (nf,oi)");
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
  dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, qrsp, 30, 0, "WMnIe qrsp");
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
  dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, qrsp, 20, 10, "WmNiE qrsp");
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
  dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, prsq, 27, 23, "WmNiE prsq");
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
  dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, prsq, 24, 22, "WMnIe prsq");
  dpd_->buf4_close(&H2);

  dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 20, 0, 20, 0, 0, "W (NF,OI)");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  dpd_->contract442(&H2, &I2, &XIA, 3, 3, -1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 30, 0, 30, 0, 0, "WMnIe qrsp");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  dpd_->contract442(&H2, &I2, &XIA, 3, 3, -1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 27, 23, 27, 23, 0, "WmNiE prsq");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
  dpd_->contract442(&H2, &I2, &XIA, 3, 3, 1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 30, 10, 30, 10, 0, "W (nf,oi)");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  dpd_->contract442(&H2, &I2, &Xia, 3, 3, -1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 20, 10, 20, 10, 0, "WmNiE qrsp");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  dpd_->contract442(&H2, &I2, &Xia, 3, 3, -1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 24, 22, 24, 22, 0, "WMnIe prsq");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
  dpd_->contract442(&H2, &I2, &Xia, 3, 3, 1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 11");
#endif

  /* term 14, +0.25 * (Rmnef Loief) * Wmnoa */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 20, 2, 20, 0, "WMNIE");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_OOOO");
  dpd_->contract442(&I2, &H2, &XIA, 3, 3, 1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 23, 27, 23, 27, 0, "WmNiE");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 23, 23, 23, 23, 0, "R2L2_oOoO");
  dpd_->contract442(&I2, &H2, &XIA, 3, 3, 1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 12, 30, 12, 30, 0, "Wmnie");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 12, 10, 12, 12, 0, "R2L2_oooo");
  dpd_->contract442(&I2, &H2, &Xia, 3, 3, 1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 22, 24, 22, 24, 0, "WMnIe");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 22, 22, 22, 22, 0, "R2L2_OoOo");
  dpd_->contract442(&I2, &H2, &Xia, 3, 3, 1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 14");
#endif

  if (!params.connect_xi) {
  /*  term 16 XIA += 0.5 Lioaf (Rmnef Wmnoe) */
  dpd_->file2_init(&Z1A, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(O,F)");
  dpd_->file2_init(&Z1B, PSIF_EOM_TMP1, R_irr, 2, 3, "Z(o,f)");

  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 21, 2, 21, 0, "WMNIE (M>N,EI)");
  dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
  dpd_->contract442(&H2, &R2, &Z1A, 3, 3, 1.0, 0.0);
  dpd_->buf4_close(&R2);
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
  dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 22, 29, 22, 29, 0, "RIjaB");
  dpd_->contract442(&H2, &R2, &Z1A, 3, 3, -1.0, 1.0);
  dpd_->buf4_close(&R2);
  dpd_->buf4_close(&H2);

  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 12, 31, 12, 31, 0, "Wmnie (m>n,ei)");
  dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 12, 15, 12, 17, 0, "Rijab");
  dpd_->contract442(&H2, &R2, &Z1B, 3, 3, 1.0, 0.0);
  dpd_->buf4_close(&R2);
  dpd_->buf4_close(&H2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
  dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 23, 28, 23, 28, 0, "RiJAb");
  dpd_->contract442(&H2, &R2, &Z1B, 3, 3, -1.0, 1.0);
  dpd_->buf4_close(&R2);
  dpd_->buf4_close(&H2);

  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
  dpd_->dot24(&Z1A, &L2, &XIA, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&L2);
  dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  dpd_->dot24(&Z1B, &L2, &XIA, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&L2);

  dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 15, 12, 17, 0, "Lijab");
  dpd_->dot24(&Z1B, &L2, &Xia, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&L2);
  dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
  dpd_->dot24(&Z1A, &L2, &Xia, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&L2);

  dpd_->file2_close(&Z1A);
  dpd_->file2_close(&Z1B);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 16");
#endif
  }

/*  term 10 XIA += 0.5 (Lmnef Rmneg) = VV(f,g) Wfiga */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_VV");
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 21, 5, 21, 7, 0, "WAMEF");
  dpd_->dot13(&I1, &H2, &XIA, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 3, 3, "LR2_vv");
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
  dpd_->dot13(&I1, &H2, &XIA, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->file2_close(&I1);

  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 3, 3, "LR2_vv");
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 31, 15, 31, 17, 0, "Wamef");
  dpd_->dot13(&I1, &H2, &Xia, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->file2_close(&I1);
  dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_VV");
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
  dpd_->dot13(&I1, &H2, &Xia, 0, 0, 1.0, 1.0);
  dpd_->buf4_close(&H2);
  dpd_->file2_close(&I1);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 10");
#endif

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

/* term 12, + (Rmnef Lmieg) Wgnaf = OVOV(nf,ig) W(gn,af) */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  /* + OVOV(NF,IG) WGNAF = OVOV(GN,IF) WGNAF(GN,AF) */
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, sprq, 21, 20, "Z (GN,IF)"); 
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 21, 5, 21, 7, 0, "WAMEF");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 21, 20, 21, 20, 0, "Z (GN,IF)");
  dpd_->contract442(&I2, &H2, &XIA, 2, 2, 1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  /* + OVOV(nf,IG) WGnAf = +OVOV(Gn,If) WGnIf(Gn,Af) */
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, sprq, 26, 24, "Z (Gn,If)"); 
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 26, 24, 26, 24, 0, "Z (Gn,If)");
  dpd_->contract442(&I2, &H2, &XIA, 2, 2, 1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  /* + OVOV(Nf,Ig) WgNIf = - OVOV(Nf,Ig) WgNfI = -OVOV(gN,If) WaMeF(gN,Af) */
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
  dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, spqr, 25, 25, "Z (gN,fI)"); 
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 25, 25, 25, 25, 0, "Z (gN,fI)");
  dpd_->contract442(&I2, &H2, &XIA, 3, 3, -1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  /* + OVOV(nf,ig) Wgnif = +OVOV(gn,if) W(gn,af) */
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, sprq, 31, 30, "Z (gn,if)"); 
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 31, 15, 31, 17, 0, "Wamef");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 31, 30, 31, 30, 0, "Z (gn,if)");
  dpd_->contract442(&I2, &H2, &Xia, 2, 2, 1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  /* + OVOV(NF,ig) WgNaF = +OVOV(gN,iF) W(gN,aF) */
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, sprq, 25, 27, "Z (gN,iF)"); 
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 25, 27, 25, 27, 0, "Z (gN,iF)");
  dpd_->contract442(&I2, &H2, &Xia, 2, 2, 1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  /* + OVOV(nF,iG) WGnaF (-WGnFi) = -OVOV(Gn,Fi) WAmEf(Gn,Fi) */
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
  dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, spqr, 26, 26, "Z (Gn,Fi)"); 
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
  dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 26, 26, 26, 26, 0, "Z (Gn,Fi)");
  dpd_->contract442(&I2, &H2, &Xia, 3, 3, -1.0, 1.0);
  dpd_->buf4_close(&I2);
  dpd_->buf4_close(&H2);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 12");
#endif

/* term 13 -0.25 (Rmnfg Weifg) Lmnea = +OOOV(MN,IE) L(MN,EA) */
  dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

  dpd_->buf4_init(&I2, PSIF_EOM_TMP, R_irr, 2, 20, 2, 20, 0, "R2Wamef_OOOV");
  dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_->contract442(&I2, &L2, &XIA, 2, 2, 1.0, 1.0);
  dpd_->buf4_close(&L2);
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, R_irr, 22, 24, 22, 24, 0, "R2Wamef_OoOv");
  dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  dpd_->contract442(&I2, &L2, &XIA, 2, 2, 1.0, 1.0);
  dpd_->buf4_close(&L2);
  dpd_->buf4_close(&I2);

  dpd_->buf4_init(&I2, PSIF_EOM_TMP, R_irr, 12, 30, 12, 30, 0, "R2Wamef_ooov");
  dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 12, 15, 12, 17, 0, "Lijab");
  dpd_->contract442(&I2, &L2, &Xia, 2, 2, 1.0, 1.0);
  dpd_->buf4_close(&L2);
  dpd_->buf4_close(&I2);
  dpd_->buf4_init(&I2, PSIF_EOM_TMP, R_irr, 23, 27, 23, 27, 0, "R2Wamef_oOoV");
  dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
  dpd_->contract442(&I2, &L2, &Xia, 2, 2, 1.0, 1.0);
  dpd_->buf4_close(&L2);
  dpd_->buf4_close(&I2);

  dpd_->file2_close(&XIA);
  dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 13");
#endif

  if (!params.connect_xi) {
    /* term 15 Linag (Rnmef Wgmef) */
    dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");

    dpd_->file2_init(&I1, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_->dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->file2_close(&I1);
    dpd_->file2_init(&I1, PSIF_EOM_TMP, R_irr, 2, 3, "R2Wamef_ov");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->file2_close(&I1);

    dpd_->file2_init(&I1, PSIF_EOM_TMP, R_irr, 2, 3, "R2Wamef_ov");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 15, 12, 17, 0, "Lijab");
    dpd_->dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->file2_close(&I1);
    dpd_->file2_init(&I1, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_->dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L2);
    dpd_->file2_close(&I1);

    dpd_->file2_close(&XIA);
    dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 15");
#endif
  }

  if (params.connect_xi) x_xi1_connected();

#ifdef DEBUG_XI
x_xi_check("extra doubles terms");
#endif
  
  return;
}

}} // namespace psi::ccdensity
