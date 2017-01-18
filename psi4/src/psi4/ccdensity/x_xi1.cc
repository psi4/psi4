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

void x_xi1_rohf(void);
extern void x_xi_check(char *term_lbl);
extern void x_xi1_connected(void);
extern void x_xi1_uhf(void);
extern void x_xi1_rhf(void);

/* compute xi_1 amplitudes for zeta equations */

void x_xi1(void) {
  if (params.ref == 0)
    x_xi1_rhf();
  else if (params.ref == 1)
    x_xi1_rohf();
  else
    x_xi1_uhf();
}

void x_xi1_rohf(void)
{
  dpdfile2 L1, XIA, Xia, I1, R1, F1, Z1A, Z1B;
  int L_irr, R_irr, G_irr;
  dpdbuf4 D, R2, L2, H2, I2;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

/*  dpd_buf4_init(&H2, CC_HBAR, 0, 2, 10, 2, 10, 0, "WMNIE");
  dpd_buf4_print(&H2,outfile,1);
  dpd_buf4_close(&H2); */

#ifdef DEBUG_XI
x_xi_check("begin xi1");
#endif
  /* term 1, XIA += 0.25 LIA Rmnef <mn||ef> */
  if ((R_irr == 0)  && (!params.connect_xi)) {
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    params.RD_overlap = 0.5 * global_dpd_->file2_trace(&I1);
    global_dpd_->file2_close(&I1);
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
    params.RD_overlap += 0.5 * global_dpd_->file2_trace(&I1);
    global_dpd_->file2_close(&I1);

    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
    global_dpd_->file2_copy(&L1, PSIF_EOM_XI, "XIA");
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    global_dpd_->file2_scm(&XIA, params.RD_overlap);
    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
    global_dpd_->file2_copy(&L1, PSIF_EOM_XI, "Xia");
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
    global_dpd_->file2_scm(&Xia, params.RD_overlap);
    global_dpd_->file2_close(&Xia);
  }
#ifdef DEBUG_XI
x_xi_check("term 1");
#endif
  /* term 2, Xia -= (Rmnef <in||ef>) * Lma */
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  global_dpd_->contract222(&I1, &L1, &XIA, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
  global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
  global_dpd_->contract222(&I1, &L1, &Xia, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 2");
#endif

  /* term 3, XIA -= 0.5 LIE (Rmnfe <mn||fa>) */
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");

  global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
  global_dpd_->contract222(&L1, &I1, &XIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_close(&L1);

  global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_vv");
  global_dpd_->contract222(&L1, &I1, &Xia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_close(&L1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 3");
#endif

  /* term 4, XIA += (Lme Rmnef) <in||af> */
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  global_dpd_->dot24(&I1, &D, &XIA, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L1R2_ov");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->dot24(&I1, &D, &XIA, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L1R2_ov");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  global_dpd_->dot24(&I1, &D, &Xia, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->dot24(&I1, &D, &Xia, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 4");
#endif

  /* XIA += (Lmnef * Rmnef) FIA */
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  params.overlap1 = global_dpd_->file2_dot(&R1, &L1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->file2_init(&L1, PSIF_CC_GL, L_irr, 0, 1, "Lia");
  params.overlap1 += global_dpd_->file2_dot(&R1, &L1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&L1);
  params.overlap2 = 1.0e0 - params.overlap1 - (params.R0 * params.L0);

  /* When (connect_xi), we still include the following term, even though Hbar
     is not connected to R.  The <Rmnef|Lmnef> Fia term here along with the
     <Rme|Lme> Fia term which is _not_ substracted out in xi_connected add up
     to (1)*Fia.  This constant term causes cclambda to be solving the
     ground-state lambda equations implicitly as well. */
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");

  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->file2_axpy(&F1, &XIA, params.overlap2, 0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->file2_axpy(&F1, &Xia, params.overlap2, 0);
  global_dpd_->file2_close(&F1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 5");
#endif

    /* term 6, -0.5 (Linef Rmnef) Fma */
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_OO");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->contract222(&I1, &F1, &XIA, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_oo");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->contract222(&I1, &F1, &Xia, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 6");
#endif
  /* term 7, -0.5 (Lmnaf Rmnef) Fie */
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_VV");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->contract222(&F1, &I1, &XIA, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_vv");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->contract222(&F1, &I1, &Xia, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

  if (!params.connect_xi) {
  /* term 8, (Fme Rmnef) Linaf) */
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(N,F)");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 2, 7, 0, "RIJAB");
  global_dpd_->dot13(&F1, &R2, &I1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  global_dpd_->dot13(&F1, &R2, &I1, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(n,f)");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 2, 7, 0, "Rijab");
  global_dpd_->dot13(&F1, &R2, &I1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  global_dpd_->dot13(&F1, &R2, &I1, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(N,F)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
  global_dpd_->dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  global_dpd_->dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(n,f)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
  global_dpd_->dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 8");
#endif
    }

  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_OO");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  global_dpd_->dot24(&I1, &H2, &XIA, 1, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_oo");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  global_dpd_->dot14(&I1, &H2, &XIA, 1, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_oo");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
  global_dpd_->dot24(&I1, &H2, &Xia, 1, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR2_OO");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  global_dpd_->dot14(&I1, &H2, &Xia, 1, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 9");
#endif

/*  term 11 XIA -= Rmnef Lmoea Winof */
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  /* this would be easier if it would work but 13 and 31 shifts are
     incompatible when symmetry is on
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_contract442(&H2, &I2, &XIA, 0, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  */
  /* if I could do a 442(0,3) I could avoid these sorts */
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, qrsp, 10, 0, "W (NF,OI)");
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
  global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, qrsp, 10, 0, "W (nf,oi)");
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, qrsp, 10, 0, "WMnIe qrsp");
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, qrsp, 10, 0, "WmNiE qrsp");
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, prsq, 10, 0, "WmNiE prsq");
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  global_dpd_->buf4_sort(&H2, PSIF_EOM_TMP1, prsq, 10, 0, "WMnIe prsq");
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 10, 0, 10, 0, 0, "W (NF,OI)");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  global_dpd_->contract442(&H2, &I2, &XIA, 3, 3, -1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 10, 0, 10, 0, 0, "WMnIe qrsp");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  global_dpd_->contract442(&H2, &I2, &XIA, 3, 3, -1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 10, 0, 10, 0, 0, "WmNiE prsq");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  global_dpd_->contract442(&H2, &I2, &XIA, 3, 3, 1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 10, 0, 10, 0, 0, "W (nf,oi)");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  global_dpd_->contract442(&H2, &I2, &Xia, 3, 3, -1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 10, 0, 10, 0, 0, "WmNiE qrsp");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  global_dpd_->contract442(&H2, &I2, &Xia, 3, 3, -1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_EOM_TMP1, 0, 10, 0, 10, 0, 0, "WMnIe prsq");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  global_dpd_->contract442(&H2, &I2, &Xia, 3, 3, 1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 11");
#endif

  /* term 14, +0.25 * (Rmnef Loief) * Wmnoa */
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 10, 2, 10, 0, "WMNIE");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_OOOO");
  global_dpd_->contract442(&I2, &H2, &XIA, 3, 3, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_oOoO");
  global_dpd_->contract442(&I2, &H2, &XIA, 3, 3, 1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 10, 2, 10, 0, "Wmnie");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_oooo");
  global_dpd_->contract442(&I2, &H2, &Xia, 3, 3, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  global_dpd_->contract442(&I2, &H2, &Xia, 3, 3, 1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 14");
#endif

  if (!params.connect_xi) {
  /*  term 16 XIA += 0.5 Lioaf (Rmnef Wmnoe) */
  global_dpd_->file2_init(&Z1A, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(O,F)");
  global_dpd_->file2_init(&Z1B, PSIF_EOM_TMP1, R_irr, 0, 1, "Z(o,f)");

  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
  global_dpd_->contract442(&H2, &R2, &Z1A, 3, 3, 1.0, 0.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
  global_dpd_->contract442(&H2, &R2, &Z1A, 3, 3, -1.0, 1.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie (m>n,ei)");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "Rijab");
  global_dpd_->contract442(&H2, &R2, &Z1B, 3, 3, 1.0, 0.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJAb");
  global_dpd_->contract442(&H2, &R2, &Z1B, 3, 3, -1.0, 1.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
  global_dpd_->dot24(&Z1A, &L2, &XIA, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->dot24(&Z1B, &L2, &XIA, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);

  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
  global_dpd_->dot24(&Z1B, &L2, &Xia, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  global_dpd_->dot24(&Z1A, &L2, &Xia, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);

  global_dpd_->file2_close(&Z1A);
  global_dpd_->file2_close(&Z1B);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 16");
#endif
  }

/*  term 10 XIA += 0.5 (Lmnef Rmneg) Wfiga */
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_VV");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
  global_dpd_->dot13(&I1, &H2, &XIA, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_vv");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  global_dpd_->dot13(&I1, &H2, &XIA, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_vv");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
  global_dpd_->dot13(&I1, &H2, &Xia, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR2_VV");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->dot13(&I1, &H2, &Xia, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&I1);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 10");
#endif

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

/* term 12, + (Rmnef Lmieg) Wgnaf */
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  global_dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, sprq, 11, 10, "Z (GN,IF)");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z (GN,IF)");
  global_dpd_->contract442(&I2, &H2, &XIA, 2, 2, 1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  global_dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, sprq, 11, 10, "Z (Gn,If)");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z (Gn,If)");
  global_dpd_->contract442(&I2, &H2, &XIA, 2, 2, 1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  global_dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, spqr, 11, 11, "Z (gN,fI)");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "Z (gN,fI)");
  global_dpd_->contract442(&I2, &H2, &XIA, 3, 3, -1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  global_dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, sprq, 11, 10, "Z (gn,if)");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z (gn,if)");
  global_dpd_->contract442(&I2, &H2, &Xia, 2, 2, 1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  global_dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, sprq, 11, 10, "Z (gN,iF)");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z (gN,iF)");
  global_dpd_->contract442(&I2, &H2, &Xia, 2, 2, 1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  global_dpd_->buf4_sort(&I2, PSIF_EOM_TMP1, spqr, 11, 11, "Z (Gn,Fi)");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "Z (Gn,Fi)");
  global_dpd_->contract442(&I2, &H2, &Xia, 3, 3, -1.0, 1.0);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&H2);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 12");
#endif

/* term 13 -0.25 (Rmnfg Weifg) Lmnea */
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, R_irr, 2, 10, 2, 10, 0, "R2Wamef_OOOV");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  global_dpd_->contract442(&I2, &L2, &XIA, 2, 2, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, R_irr, 0, 10, 0, 10, 0, "R2Wamef_OoOv");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->contract442(&I2, &L2, &XIA, 2, 2, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&I2);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, R_irr, 2, 10, 2, 10, 0, "R2Wamef_ooov");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  global_dpd_->contract442(&I2, &L2, &Xia, 2, 2, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, R_irr, 0, 10, 0, 10, 0, "R2Wamef_oOoV");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  global_dpd_->contract442(&I2, &L2, &Xia, 2, 2, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&I2);

  global_dpd_->file2_close(&XIA);
  global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 13");
#endif

  if (!params.connect_xi) {
    /* term 15 Linag (Rnmef Wgmef) */
    global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I1);
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_ov");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I1);

    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_ov");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
    global_dpd_->dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I1);
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&I1);

    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 15");
#endif
  }

  if (params.connect_xi) x_xi1_connected();

  return;
}



void x_xi_zero(void) {
  dpdfile2 XIA, Xia;
  dpdbuf4 XIJAB, Xijab, XIjAb;
  int G_irr;
  G_irr = params.G_irr;

  if (params.ref == 0) {
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    global_dpd_->file2_scm(&XIA, 0.0);
    global_dpd_->file2_close(&XIA);
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    global_dpd_->buf4_scm(&XIjAb, 0.0);
    global_dpd_->buf4_close(&XIjAb);
  }
  else if (params.ref == 1) {
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    global_dpd_->file2_scm(&XIA, 0.0);
    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
    global_dpd_->file2_scm(&Xia, 0.0);
    global_dpd_->file2_close(&Xia);
    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    global_dpd_->buf4_scm(&XIJAB, 0.0);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
    global_dpd_->buf4_scm(&Xijab, 0.0);
    global_dpd_->buf4_close(&Xijab);
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    global_dpd_->buf4_scm(&XIjAb, 0.0);
    global_dpd_->buf4_close(&XIjAb);
  }
  else {
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
    global_dpd_->file2_scm(&XIA, 0.0);
    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 2, 3, "Xia");
    global_dpd_->file2_scm(&Xia, 0.0);
    global_dpd_->file2_close(&Xia);
    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    global_dpd_->buf4_scm(&XIJAB, 0.0);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 12, 17, 12, 17, 0, "Xijab");
    global_dpd_->buf4_scm(&Xijab, 0.0);
    global_dpd_->buf4_close(&Xijab);
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
    global_dpd_->buf4_scm(&XIjAb, 0.0);
    global_dpd_->buf4_close(&XIjAb);
  }
  return;
}

}} // namespace psi::ccdensity
