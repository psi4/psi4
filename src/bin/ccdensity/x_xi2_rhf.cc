/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "math.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void x_xi2_4_rhf(void);
extern void x_xi2_14(void);
extern void x_xi_check(char *term_lbl);
extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);

/* compute xi_2 amplitudes for RHF wavefunctions for zeta equations */

void x_xi2_rhf(void)
{
  dpdfile2 L1, XIA, Xia, I1, R1, F1, Z1A, Z1B;
  int L_irr, R_irr, G_irr;
  double tval;
  dpdbuf4 D2, R2, L2, H2, I2, Z, Z2, XIJAB, Xijab, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

#ifdef DEBUG_XI
  x_xi_check("reset");
#endif

  /* terms 1 and 5, Xijab += (Lme Rme + 0.25 Lmnef Rmnef) <ij||eb> */
  /* overlaps in params are assigned in x_xi1.c */
  /* see comments in xi1_connected.c */
  dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_scmcopy(&D2, EOM_XI, "XIjAb", params.overlap1+params.overlap2);
  dpd_buf4_close(&D2);
#ifdef DEBUG_XI
x_xi_check("terms 1 and 5");
#endif

  /* terms 2 and 9, Xijab -= P(ab) (Lma Rme + Lmnfa Rmnfe) <ij||eb */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract424(&D2, &I1, &Z2, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&D2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_axpy(&Z2, &XIjAb, -1.0);
  dpd_buf4_sort(&Z2, EOM_TMP1, qpsr, 0, 5, "X (jI,bA)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (jI,bA)");
  dpd_buf4_axpy(&Z2, &XIjAb, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 2 and 9");
#endif

  /* terms 3 and 10, Xijab -= P(ij) (Lie Rme + 0.5 Linef Rmnef) <mj||ab> */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_OO");
  dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract244(&I1, &D2, &Z2, 1, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&D2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_axpy(&Z2, &XIjAb, -1.0);
  dpd_buf4_sort(&Z2, EOM_TMP1, qpsr, 0, 5, "X (jI,bA)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (jI,bA)");
  dpd_buf4_axpy(&Z2, &XIjAb, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 3 and 10");
#endif

  x_xi2_4_rhf();
#ifdef DEBUG_XI
x_xi_check("terms 4 and 6");
#endif

  /* term 7, Xijab += 0.25 Lmnab Rmnef <ij||ef> */
  dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_init(&Z2, EOM_TMP1, R_irr, 0, 0, 0, 0, 0, "Z (Ij,Mn)");
  dpd_contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0); 
  dpd_buf4_close(&R2);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract444(&Z2, &L2, &XIjAb, 0, 1, 1.0, 1.0); 
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

  /* term 8, Xijab += 0.25 Rmnef Lijef <mn||ab> */
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract444(&I2, &D2, &XIjAb, 1, 1, 1.0, 1.0); 
  dpd_buf4_close(&D2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 8");
#endif

  /* term 11, Xijab -= 0.5 P(ab) Lijfb (Rmnef <mn||ea>) */
  /* term 17        -=     P(ab) Lijfb (Rmf Fma) */
  /* term 20        +=     P(ab) Lijfb (Rme Wfmae) */
  /* build 1-e intermediates to include term 17 */
      /* for term 11: */
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
  dpd_file2_copy(&I1, EOM_TMP1, "X1 (F,A)");
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "X1 (F,A)");
      /* for term 20: */
  dpd_file2_init(&Z1A, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");
  dpd_file2_axpy(&Z1A, &I1, -1.0, 0);
  dpd_file2_close(&Z1A);
      /* for term 17: */
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&R1, &F1, &I1, 1, 1, 1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&R1);
  dpd_file2_close(&I1);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "X1 (F,A)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract244(&I1, &L2, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_axpy(&Z2, &XIjAb, -1.0);
  dpd_buf4_close(&XIjAb);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 0, 5, "XIjAb", -1.0);
  dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("terms 11, 17 and 20");
#endif

  /* term 12, Xijab -= 0.5 P(ij) Lmjab (Rmnef <in||ef>) */
  /* term 16,       -=     P(ij) Lmjab (Rme Fie) */
  /* term 21,       -=     P(ij) Lmjab (Rne Winme) */
  /* make 1-electron intermediates to include terms 16 and 21 as well */
      /* for term 12: */
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  dpd_file2_copy(&I1, EOM_TMP1, "X1 (M,I)");
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "X1 (M,I)");
       /* for term 21 */
  dpd_file2_init(&Z1A, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");
  dpd_file2_axpy(&Z1A, &I1, 1.0, 1);
  dpd_file2_close(&Z1A);
      /* for term 16: */
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&R1, &F1, &I1, 0, 0, 1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&R1);
  dpd_file2_close(&I1);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "X1 (M,I)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract244(&I1, &L2, &Z2, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_axpy(&Z2, &XIjAb, -1.0);
  dpd_buf4_close(&XIjAb);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 0, 5, "XIjAb", -1.0);
  dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 12, 16 and 21");
#endif

  /* term 13 + 15, (0.25 Rmnef <mn||ef> + Rme Fme) Lijab */
  if ( (L_irr == 0) && (!params.connect_xi)) {
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
    tval = 2.0 * dpd_file2_dot(&R1, &F1);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    tval += params.RD_overlap;

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_axpy(&L2, &XIjAb, tval);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 13 (ijab) and 15 (Fme)");
#endif
  }

  /* term 14, +P(ij) P(ab) Lmjeb Rme Fia */
  if (!params.connect_xi) {
    x_xi2_14();
#ifdef DEBUG_XI
x_xi_check("term 14 (Fme)");
#endif
  }

  /* term 22, +P(ij) (Limfe Rme) Wfjab */
  if (!params.connect_xi) {
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_axpy(&Z2, &XIjAb, 1.0);
    dpd_buf4_close(&XIjAb);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
    dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 22 (Wamef)");
#endif
  }

  /* term 23, -P(ab) (Lmnea Rme) Wijnb */
  if (!params.connect_xi) {
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    dpd_contract244(&I1, &H2, &Z2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_axpy(&Z2, &XIjAb, -1.0);
    dpd_buf4_close(&XIjAb);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 0, 5, "XIjAb", -1.0);
    dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 23 (Wmnie)");
#endif
  }

  /* term 25, Xijab += (Lnmab Rme) Wijne */
  dpd_buf4_init(&Z2, EOM_TMP1, R_irr, 0, 0, 0, 0, 0, "X (Ij,Nm)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&H2, &R1, &Z2, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract444(&Z2, &L2, &Z, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_axpy(&Z, &XIjAb, 1.0);
  dpd_buf4_close(&XIjAb);
  dpd_buf4_sort_axpy(&Z, EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
  dpd_buf4_close(&Z);
#ifdef DEBUG_XI
x_xi_check("term 25 (Wmnie)");
#endif

  /* term 24, Xijab -= (Lijfe Rme) Wfmab */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_contract444(&I2, &H2, &Z2, 0, 1, 1.0, 0.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_axpy(&Z2, &XIjAb, -1.0);
  dpd_buf4_close(&XIjAb);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 0, 5, "XIjAb", -1.0);
  dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 24 (Wamef)");
#endif

  /* terms 18, 19: Xijab -= P(ij) P(ab) Linae (Rme Wmjnb + Rnf Wejbf) */
  /* construct Z(JB,NE) = RME WMJNB + RNF WEJBF */
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 11, 10, 11, 10, 0, "Z (Ej,Nb)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_contract244(&R1, &H2, &Z, 1, 2, 1, -1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&R1);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 10, 10, "Z (jE,Nb)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jE,Nb)");
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 10, 10, "Z (jb,NE)");
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Je,Nb)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&H2, &R1, &Z, 1, 0, 1, -1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 11, 11, 11, 11, 0, "Z (eJ,bN)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_sort_axpy(&Z, EOM_TMP1, qpsr, 10, 10, "Z (Je,Nb)", 1.0);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Je,Nb)");
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 10, 10, "Z (Jb,Ne)");
  dpd_buf4_close(&Z);

  /* XIjAb += - L(IA,NE) Z(jb,NE) - L(IA,ne) Z(jb,ne)  */
  /*          - L(jb,ne) Z(IA,ne) - L(jb,NE) Z(IA,NE)  transpose of line above */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "XIjAb (IA,jb)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "2LIjAb - LIjbA (IA,jb)");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,NE)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Jb,Ne)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);

  dpd_buf4_sort(&Z2, EOM_TMP1, rspq, 10, 10, "XIjAb (jb,IA)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "XIjAb (jb,IA)");
  dpd_buf4_axpy(&Z, &Z2, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, prqs, 0, 5, "XIjAb", 1.0);
  dpd_buf4_close(&Z2);

  /* XIjAb += + L(jA,Ne) Z(Ib,Ne) + L(Ib,nE) Z(jA,nE)  */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (Ib,jA)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Jb,Ne)");
  dpd_contract444(&Z, &L2, &Z2, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);

  dpd_buf4_sort(&Z2, EOM_TMP1, rspq, 10, 10, "Z2 (jA,Ib)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (jA,Ib)");
  dpd_buf4_axpy(&Z, &Z2, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, prsq, 0, 5, "XIjAb", 1.0);
  dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("terms 18, 19 (Wmnie, Wamef)");
#endif

  /* Write irrep of XI amplitudes to CC_INFO */
  psio_write_entry(CC_INFO, "XI Irrep", (char *) &G_irr,sizeof(int));

  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  tval = 2.0 * dpd_file2_dot_self(&XIA);
  dpd_file2_close(&XIA);
  fprintf(outfile,"XI_IA amplitudes: Norm=%15.10lf, Dot=%15.10lf\n", sqrt(tval), tval );

  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_sort(&XIjAb, EOM_TMP1, pqsr, 0, 5, "XIjbA");
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjbA");
  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  tval = norm_C_rhf(&XIA, &XIjAb, &Z2);
  dpd_file2_close(&XIA);
  dpd_buf4_close(&XIjAb);
  dpd_buf4_close(&Z2);
  fprintf(outfile,"XI amplitudes   : Norm=%15.10lf, Dot=%15.10lf\n", sqrt(tval), tval );

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);
  return;
}



/* compute terms 4 and 6 of xi2 amplitudes */
/* Xijab += P(ij) P(ab) (Rme Lia + Rmnef Linaf) <mj||eb> */
void x_xi2_4_rhf(void)
{
  dpdfile2 RIA, LIA;
  int L_irr, R_irr, G_irr, nirreps;
  int I, A, M, E, i, a, m, e, h, row, col, Isym, Esym, Asym, Msym;
  dpdbuf4 D, R2, L2, H2, I2, Z, Z2, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* construct RL = Rme Lia + Rmnef Linaf */
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_copy(&I2, EOM_TMP1, "RL_OVov");
  dpd_buf4_close(&I2);

  /* RL_OVOV(me,ia) += Rme Lia */
  dpd_file2_init(&RIA, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&LIA, CC_GL, L_irr, 0, 1, "LIA");

  dpd_file2_mat_init(&RIA);
  dpd_file2_mat_init(&LIA);
  dpd_file2_mat_rd(&RIA);
  dpd_file2_mat_rd(&LIA);

  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&I2, h);
    dpd_buf4_mat_irrep_rd(&I2, h);
    for(row=0; row < I2.params->rowtot[h]; row++) {
      m = I2.params->roworb[h][row][0];
      e = I2.params->roworb[h][row][1];
      M = RIA.params->rowidx[m]; Msym = RIA.params->psym[m];
      E = RIA.params->colidx[e]; Esym = RIA.params->qsym[e];
      for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
        i = I2.params->colorb[h^G_irr][col][0];
        a = I2.params->colorb[h^G_irr][col][1];
        I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
        A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
        if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
          I2.matrix[h][row][col] += RIA.matrix[Msym][M][E] * LIA.matrix[Isym][I][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&I2, h);
    dpd_buf4_mat_irrep_close(&I2, h);
  }
  dpd_buf4_close(&I2);

  /* XIjAb += RL_OVov(me,IA) * (2<Mj|Eb> - <Mj|Be>) for the ZOVOV alpha case
   and the reverse of this multiplication for the Zovov beta case */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "X2 (IA,jb)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
  dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);
  /* XIjAb += RL_OvOv(Me,Ia) * (<Mj|Eb> (ME,jb) and the reverse of this
   to finish the all beta and all alpha spin orbital components */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_contract444(&D, &Z, &Z2, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "XIjAb");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_axpy(&Z2, &XIjAb, 1.0);
  dpd_buf4_close(&XIjAb);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
  dpd_buf4_close(&Z2);

  /* Now do Z alpha beta parts */
  /* XIjAb += ZmEjA <mI|bE> (mE,Ib) + ZMeIb <Mj|Ae> (Me,jA) */
  dpd_buf4_init(&Z2, EOM_TMP1, 0, 10, 10, 10, 10, 0, "X2 (Ib,jA)");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);

  dpd_buf4_sort(&Z2, EOM_TMP1, prsq, 0, 5, "XIjAb");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  dpd_buf4_axpy(&Z2, &XIjAb, 1.0);
  dpd_buf4_close(&XIjAb);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
  dpd_buf4_close(&Z2);

  return;
}


}} // namespace psi::ccdensity
