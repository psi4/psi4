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
    dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    params.RD_overlap = 0.5 * dpd_file2_trace(&I1);
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
    params.RD_overlap += 0.5 * dpd_file2_trace(&I1);
    dpd_file2_close(&I1);

    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, EOM_XI, "XIA");
    dpd_file2_close(&L1);
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_scm(&XIA, params.RD_overlap);
    dpd_file2_close(&XIA);
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
    dpd_file2_copy(&L1, EOM_XI, "Xia");
    dpd_file2_close(&L1);
    dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
    dpd_file2_scm(&Xia, params.RD_overlap);
    dpd_file2_close(&Xia);
  }
#ifdef DEBUG_XI
x_xi_check("term 1");
#endif
  /* term 2, Xia -= (Rmnef <in||ef>) * Lma */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");

  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
  dpd_contract222(&I1, &L1, &XIA, 1, 1, -1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
  dpd_contract222(&I1, &L1, &Xia, 1, 1, -1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_file2_close(&I1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 2");
#endif

  /* term 3, XIA -= 0.5 LIE (Rmnfe <mn||fa>) */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");

  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
  dpd_contract222(&L1, &I1, &XIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_file2_close(&L1);

  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 1, 1, "RD_vv");
  dpd_contract222(&L1, &I1, &Xia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_file2_close(&L1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 3");
#endif

  /* term 4, XIA += (Lme Rmnef) <in||af> */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");

  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_dot24(&I1, &D, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L1R2_ov");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_dot24(&I1, &D, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_file2_close(&I1);

  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L1R2_ov");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_dot24(&I1, &D, &Xia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_dot24(&I1, &D, &Xia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_file2_close(&I1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 4");
#endif

  /* XIA += (Lmnef * Rmnef) FIA */ 
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
  params.overlap1 = dpd_file2_dot(&R1, &L1);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
  params.overlap1 += dpd_file2_dot(&R1, &L1);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);
  params.overlap2 = 1.0e0 - params.overlap1 - (params.R0 * params.L0);

  /* When (connect_xi), we still include the following term, even though Hbar
     is not connected to R.  The <Rmnef|Lmnef> Fia term here along with the
     <Rme|Lme> Fia term which is _not_ substracted out in xi_connected add up
     to (1)*Fia.  This constant term causes cclambda to be solving the
     ground-state lambda equations implicitly as well. */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");

  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_axpy(&F1, &XIA, params.overlap2, 0);
  dpd_file2_close(&F1);
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
  dpd_file2_axpy(&F1, &Xia, params.overlap2, 0);
  dpd_file2_close(&F1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 5");
#endif

    /* term 6, -0.5 (Linef Rmnef) Fma */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");

  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR2_OO");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&I1, &F1, &XIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR2_oo");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
  dpd_contract222(&I1, &F1, &Xia, 0, 1, -1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&I1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 6");
#endif
  /* term 7, -0.5 (Lmnaf Rmnef) Fie */
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");

  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR2_VV");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&F1, &I1, &XIA, 0, 0, -1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR2_vv");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
  dpd_contract222(&F1, &I1, &Xia, 0, 0, -1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&I1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

  if (!params.connect_xi) {
  /* term 8, (Fme Rmnef) Linaf) */
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 1, "Z(N,F)");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 2, 7, 0, "RIJAB");
  dpd_dot13(&F1, &R2, &I1, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&R2);
  dpd_file2_close(&F1);
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_dot13(&F1, &R2, &I1, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&R2);
  dpd_file2_close(&F1);
  dpd_file2_close(&I1);

  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 1, "Z(n,f)");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 2, 7, 0, "Rijab");
  dpd_dot13(&F1, &R2, &I1, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&R2);
  dpd_file2_close(&F1);
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_dot13(&F1, &R2, &I1, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&R2);
  dpd_file2_close(&F1);
  dpd_file2_close(&I1);

  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 1, "Z(N,F)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
  dpd_dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);

  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 1, "Z(n,f)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
  dpd_dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 8");
#endif
    }

  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR2_OO");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  dpd_dot24(&I1, &H2, &XIA, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR2_oo");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  dpd_dot14(&I1, &H2, &XIA, 1, 0, -1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);

  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR2_oo");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
  dpd_dot24(&I1, &H2, &Xia, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR2_OO");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_dot14(&I1, &H2, &Xia, 1, 0, -1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 9");
#endif

/*  term 11 XIA -= Rmnef Lmoea Winof */
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

  /* this would be easier if it would work but 13 and 31 shifts are
     incompatible when symmetry is on
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_contract442(&H2, &I2, &XIA, 0, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  */
  /* if I could do a 442(0,3) I could avoid these sorts */
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  dpd_buf4_sort(&H2, EOM_TMP1, qrsp, 10, 0, "W (NF,OI)");
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
  dpd_buf4_sort(&H2, EOM_TMP1, qrsp, 10, 0, "W (nf,oi)");
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_buf4_sort(&H2, EOM_TMP1, qrsp, 10, 0, "WMnIe qrsp");
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  dpd_buf4_sort(&H2, EOM_TMP1, qrsp, 10, 0, "WmNiE qrsp");
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  dpd_buf4_sort(&H2, EOM_TMP1, prsq, 10, 0, "WmNiE prsq");
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_buf4_sort(&H2, EOM_TMP1, prsq, 10, 0, "WMnIe prsq");
  dpd_buf4_close(&H2);

  dpd_buf4_init(&H2, EOM_TMP1, 0, 10, 0, 10, 0, 0, "W (NF,OI)");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_contract442(&H2, &I2, &XIA, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, EOM_TMP1, 0, 10, 0, 10, 0, 0, "WMnIe qrsp");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_contract442(&H2, &I2, &XIA, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, EOM_TMP1, 0, 10, 0, 10, 0, 0, "WmNiE prsq");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  dpd_contract442(&H2, &I2, &XIA, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_buf4_init(&H2, EOM_TMP1, 0, 10, 0, 10, 0, 0, "W (nf,oi)");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_contract442(&H2, &I2, &Xia, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, EOM_TMP1, 0, 10, 0, 10, 0, 0, "WmNiE qrsp");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_contract442(&H2, &I2, &Xia, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, EOM_TMP1, 0, 10, 0, 10, 0, 0, "WMnIe prsq");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_contract442(&H2, &I2, &Xia, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 11");
#endif

  /* term 14, +0.25 * (Rmnef Loief) * Wmnoa */
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

  dpd_buf4_init(&H2, CC_HBAR, 0, 2, 10, 2, 10, 0, "WMNIE");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_OOOO");
  dpd_contract442(&I2, &H2, &XIA, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_oOoO");
  dpd_contract442(&I2, &H2, &XIA, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_buf4_init(&H2, CC_HBAR, 0, 2, 10, 2, 10, 0, "Wmnie");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_oooo");
  dpd_contract442(&I2, &H2, &Xia, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  dpd_contract442(&I2, &H2, &Xia, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 14");
#endif

  if (!params.connect_xi) {
  /*  term 16 XIA += 0.5 Lioaf (Rmnef Wmnoe) */
  dpd_file2_init(&Z1A, EOM_TMP1, R_irr, 0, 1, "Z(O,F)");
  dpd_file2_init(&Z1B, EOM_TMP1, R_irr, 0, 1, "Z(o,f)");

  dpd_buf4_init(&H2, CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  dpd_buf4_init(&R2, CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
  dpd_contract442(&H2, &R2, &Z1A, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
  dpd_contract442(&H2, &R2, &Z1A, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&H2);

  dpd_buf4_init(&H2, CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie (m>n,ei)");
  dpd_buf4_init(&R2, CC_GR, R_irr, 2, 5, 2, 7, 0, "Rijab");
  dpd_contract442(&H2, &R2, &Z1B, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJAb");
  dpd_contract442(&H2, &R2, &Z1B, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&H2);

  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
  dpd_dot24(&Z1A, &L2, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_dot24(&Z1B, &L2, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
  dpd_dot24(&Z1B, &L2, &Xia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_dot24(&Z1A, &L2, &Xia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);

  dpd_file2_close(&Z1A);
  dpd_file2_close(&Z1B);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 16");
#endif
  }

/*  term 10 XIA += 0.5 (Lmnef Rmneg) Wfiga */
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR2_VV");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
  dpd_dot13(&I1, &H2, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR2_vv");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  dpd_dot13(&I1, &H2, &XIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);

  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR2_vv");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
  dpd_dot13(&I1, &H2, &Xia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR2_VV");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_dot13(&I1, &H2, &Xia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&I1);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 10");
#endif

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

/* term 12, + (Rmnef Lmieg) Wgnaf */
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_buf4_sort(&I2, EOM_TMP1, sprq, 11, 10, "Z (GN,IF)"); 
  dpd_buf4_close(&I2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z (GN,IF)");
  dpd_contract442(&I2, &H2, &XIA, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_buf4_sort(&I2, EOM_TMP1, sprq, 11, 10, "Z (Gn,If)"); 
  dpd_buf4_close(&I2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z (Gn,If)");
  dpd_contract442(&I2, &H2, &XIA, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_sort(&I2, EOM_TMP1, spqr, 11, 11, "Z (gN,fI)"); 
  dpd_buf4_close(&I2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "Z (gN,fI)");
  dpd_contract442(&I2, &H2, &XIA, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_buf4_sort(&I2, EOM_TMP1, sprq, 11, 10, "Z (gn,if)"); 
  dpd_buf4_close(&I2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z (gn,if)");
  dpd_contract442(&I2, &H2, &Xia, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_sort(&I2, EOM_TMP1, sprq, 11, 10, "Z (gN,iF)"); 
  dpd_buf4_close(&I2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z (gN,iF)");
  dpd_contract442(&I2, &H2, &Xia, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  dpd_buf4_sort(&I2, EOM_TMP1, spqr, 11, 11, "Z (Gn,Fi)"); 
  dpd_buf4_close(&I2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "Z (Gn,Fi)");
  dpd_contract442(&I2, &H2, &Xia, 3, 3, -1.0, 1.0);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&H2);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 12");
#endif

/* term 13 -0.25 (Rmnfg Weifg) Lmnea */
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

  dpd_buf4_init(&I2, EOM_TMP, R_irr, 2, 10, 2, 10, 0, "R2Wamef_OOOV");
  dpd_buf4_init(&L2, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_contract442(&I2, &L2, &XIA, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&I2);
  dpd_buf4_init(&I2, EOM_TMP, R_irr, 0, 10, 0, 10, 0, "R2Wamef_OoOv");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract442(&I2, &L2, &XIA, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP, R_irr, 2, 10, 2, 10, 0, "R2Wamef_ooov");
  dpd_buf4_init(&L2, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_contract442(&I2, &L2, &Xia, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&I2);
  dpd_buf4_init(&I2, EOM_TMP, R_irr, 0, 10, 0, 10, 0, "R2Wamef_oOoV");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract442(&I2, &L2, &Xia, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&I2);

  dpd_file2_close(&XIA);
  dpd_file2_close(&Xia);
#ifdef DEBUG_XI
x_xi_check("term 13");
#endif

  if (!params.connect_xi) {
    /* term 15 Linag (Rnmef Wgmef) */
    dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");

    dpd_file2_init(&I1, EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP, R_irr, 0, 1, "R2Wamef_ov");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_dot24(&I1, &L2, &XIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);

    dpd_file2_init(&I1, EOM_TMP, R_irr, 0, 1, "R2Wamef_ov");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
    dpd_dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_dot24(&I1, &L2, &Xia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);

    dpd_file2_close(&XIA);
    dpd_file2_close(&Xia);
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
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_scm(&XIA, 0.0);
    dpd_file2_close(&XIA);
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_scm(&XIjAb, 0.0);
    dpd_buf4_close(&XIjAb);
  } 
  else if (params.ref == 1) {
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_scm(&XIA, 0.0);
    dpd_file2_close(&XIA);
    dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
    dpd_file2_scm(&Xia, 0.0);
    dpd_file2_close(&Xia);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_scm(&XIJAB, 0.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
    dpd_buf4_scm(&Xijab, 0.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_scm(&XIjAb, 0.0);
    dpd_buf4_close(&XIjAb);
  }
  else {
    dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
    dpd_file2_scm(&XIA, 0.0);
    dpd_file2_close(&XIA);
    dpd_file2_init(&Xia, EOM_XI, G_irr, 2, 3, "Xia");
    dpd_file2_scm(&Xia, 0.0);
    dpd_file2_close(&Xia);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_scm(&XIJAB, 0.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 17, 12, 17, 0, "Xijab");
    dpd_buf4_scm(&Xijab, 0.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
    dpd_buf4_scm(&XIjAb, 0.0);
    dpd_buf4_close(&XIjAb);
  }
  return;
}

}} // namespace psi::ccdensity
