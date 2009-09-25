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

/* compute xi_1 amplitudes for zeta equations */

void x_xi1_rhf(void)
{
  dpdfile2 L1, XIA, Xia, I1, R1, F1, Z1A, Z1B;
  int L_irr, R_irr, G_irr;
  dpdbuf4 D, R2, L2, H2, I2, Z2;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

#ifdef DEBUG_XI
x_xi_check("begin xi1");
#endif
  /* term 1, XIA += 0.25 LIA Rmnef <mn||ef> */
  if ((R_irr == 0)  && (!params.connect_xi)) {
    dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    params.RD_overlap = dpd_file2_trace(&I1);
    dpd_file2_close(&I1);
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, EOM_XI, "XIA");
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_scm(&L1, params.RD_overlap);
    dpd_file2_close(&L1);
  }
#ifdef DEBUG_XI
x_xi_check("term 1");
#endif

  /* term 2, Xia -= (Rmnef <in||ef>) * Lma */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
  dpd_contract222(&I1, &L1, &XIA, 1, 1, -1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_file2_close(&I1);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 2");
#endif

  /* term 3, XIA -= 0.5 LIE (Rmnfe <mn||fa>) */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
  dpd_contract222(&L1, &I1, &XIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_file2_close(&L1);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 3");
#endif

  /* term 4, XIA += (Lme Rmnef) <in||af> */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_dot24(&I1, &D, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_file2_close(&I1);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 4");
#endif

  /* term 5, XIA += (Lmnef * Rmnef) FIA */ 
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
  params.overlap1 = 2.0 * dpd_file2_dot(&R1, &L1);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);
  params.overlap2 = 1.0e0 - params.overlap1 - (params.R0 * params.L0);

  /* When (connect_xi), we still include the following term, even though Hbar
     is not connected to R.  The <Rmnef|Lmnef> Fia term here along with the
     <Rme|Lme> Fia term which is _not_ substracted out in xi_connected add up
     to (1)*Fia.  This constant term causes cclambda to be solving the
     ground-state lambda equations implicitly when it solves for zeta. */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_axpy(&F1, &XIA, params.overlap2, 0);
  dpd_file2_close(&F1);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 5");
#endif

  /* term 6, XIA -= (0.5 Linef Rmnef) Fma */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR2_OO");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&I1, &F1, &XIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&I1);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 6");
#endif

  /* term 7, XIA -= (0.5 Lmnaf Rmnef) Fie */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR2_VV");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&F1, &I1, &XIA, 0, 0, -1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&I1);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

  if (!params.connect_xi) {
    /* term 8, XIA += (Fme Rmnef) Linaf) */
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 1, "Z(N,F)");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
    dpd_dot13(&F1, &R2, &I1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_file2_close(&F1);
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
    dpd_dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&XIA);
    dpd_file2_close(&I1);
#ifdef DEBUG_XI
x_xi_check("term 8");
#endif
  }

  /* term 9, XIA -= (0.5 Lmnef Rmoef) Woina */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR2_OO");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
  dpd_dot14(&I1, &H2, &XIA, 1, 0, -1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 9");
#endif

/*  term 10 XIA += (0.5 Lmnef Rmneg) Wfiga */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR2_VV");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
  dpd_dot13(&I1, &H2, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 10");
#endif

/*  term 11 XIA -= Rmnef Lmoea Winof */
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
  dpd_buf4_sort(&H2, EOM_TMP1, qrsp, 10, 0, "WMnIe - 2WnMIe qrsp");
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
  dpd_buf4_sort(&H2, EOM_TMP1, qrsp, 10, 0, "2WMnIe - WnMIe qrsp");
  dpd_buf4_close(&H2);

  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_buf4_init(&H2, EOM_TMP1, 0, 10, 0, 10, 0, 0, "2WMnIe - WnMIe qrsp");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_contract442(&H2, &I2, &XIA, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, EOM_TMP1, 0, 10, 0, 10, 0, 0, "WMnIe - 2WnMIe qrsp");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_contract442(&H2, &I2, &XIA, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 11");
#endif

/* term 12, + (Rmnef Lmieg) Wgnaf */
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_copy(&I2, EOM_TMP1, "R2L2 2OVov + OvOv");
  dpd_buf4_close(&I2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "R2L2 2OVov + OvOv");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_axpy(&I2, &Z2, 2.0);
  dpd_buf4_close(&I2);
  dpd_buf4_sort(&Z2, EOM_TMP1, sprq, 11, 10, "2OVov + OvOv (Gn,If)"); 
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_copy(&I2, EOM_TMP1, "R2L2 OVov + 2OvOv");
  dpd_buf4_close(&I2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "R2L2 OVov + 2OvOv");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_axpy(&I2, &Z2, 2.0);
  dpd_buf4_close(&I2);
  dpd_buf4_sort(&Z2, EOM_TMP1, spqr, 11, 11, "OVov + 2OvOv (Gn,fI)"); 
  dpd_buf4_close(&Z2);

  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "2OVov + OvOv (Gn,If)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_contract442(&I2, &H2, &XIA, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "OVov + 2OvOv (Gn,fI)");
  dpd_contract442(&I2, &H2, &XIA, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 12");
#endif

/* term 13 -0.25 (Rmnfg Weifg) Lmnea */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_buf4_init(&I2, EOM_TMP, R_irr, 0, 10, 0, 10, 0, "R2Wamef_OoOv");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
  dpd_contract442(&I2, &L2, &XIA, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&I2);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 13");
#endif

  /* term 14, +0.25 * (Rmnef Loief) * Wmnoa */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  dpd_contract442(&I2, &H2, &XIA, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 14");
#endif

  /* term 15 Linag (Rnmef Wgmef) */
  if (!params.connect_xi) {
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_init(&I1, EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
    dpd_dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_file2_close(&XIA);
#ifdef DEBUG_XI
x_xi_check("term 15");
#endif
  }

  /*  term 16 XIA += 0.5 Lioaf (Rmnef Wmnoe) */
  if (!params.connect_xi) {
    dpd_file2_init(&Z1A, EOM_TMP1, R_irr, 0, 1, "Z(O,F)");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_contract442(&H2, &R2, &Z1A, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&R2);
    dpd_buf4_close(&H2);
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
    dpd_dot24(&Z1A, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&XIA);
    dpd_file2_close(&Z1A);
#ifdef DEBUG_XI
x_xi_check("term 16");
#endif
  }

  if (params.connect_xi) x_xi1_connected();
  
  return;
}

}} // namespace psi::ccdensity
