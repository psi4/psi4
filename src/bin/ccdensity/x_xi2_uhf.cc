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

extern void x_xi_check(char *term_lbl);
extern void x_xi2_14(void); /* in x_xi2.c */
void x_xi2_4_uhf(void); /* in x_xi2.c */

/* compute UHF xi_2 amplitudes for zeta equations */

void x_xi2_uhf(void)
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
  dpd_buf4_init(&D2, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  dpd_buf4_scmcopy(&D2, EOM_XI, "XIJAB", params.overlap1+params.overlap2);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&D2, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_scmcopy(&D2, EOM_XI, "Xijab", params.overlap1+params.overlap2);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&D2, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_buf4_scmcopy(&D2, EOM_XI, "XIjAb", params.overlap1+params.overlap2);
  dpd_buf4_close(&D2);
#ifdef DEBUG_XI
x_xi_check("terms 1 and 5");
#endif

  /* terms 2 and 9, XIJAB -= P(AB) (LMA RME + LMNFA RMNFE) <IJ||EB> */
  dpd_buf4_init(&Z2, EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z (I>J,AB)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_buf4_init(&D2, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
  dpd_contract244(&I1, &D2, &Z2, 1, 2, 1, 1.0, 0.0);
  dpd_buf4_close(&D2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
  dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z2, EOM_TMP1, 0, 12, 15, 12, 15, 0, "Z (i>j,ab)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 3, 3, "LR_vv");
  dpd_buf4_init(&D2, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
  dpd_contract244(&I1, &D2, &Z2, 1, 2, 1, 1.0, 0.0);
  dpd_buf4_close(&D2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 15, 12, 17, 0, "Xijab");
  dpd_buf4_axpy(&Z2, &Xijab, -1.0);
  dpd_buf4_close(&Xijab);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 12, 17, "Xijab", 1.0);
  dpd_buf4_close(&Z2);

  /* -LR_VV(A,E) <Ij||Eb> + LR_vv(b,e) <Ij||eA> */
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_buf4_init(&D2, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_contract244(&I1, &D2, &XIjAb, 1, 2, 1, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 3, 3, "LR_vv");
  dpd_contract424(&D2, &I1, &XIjAb, 3, 1, 0, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&D2);
  dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 2 and 9");
#endif

  /* terms 3 and 10, Xijab -= P(ij) (Lie Rme + 0.5 Linef Rmnef) <mj||ab> */
  dpd_buf4_init(&Z2, EOM_TMP1, 0, 0, 7, 0, 7, 0, "Z (IJ,A>B)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_OO");
  dpd_buf4_init(&D2, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
  dpd_contract244(&I1, &D2, &Z2, 1, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&D2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
  dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z2, EOM_TMP1, 0, 10, 17, 10, 17, 0, "Z (ij,a>b)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 2, "LR_oo");
  dpd_buf4_init(&D2, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
  dpd_contract244(&I1, &D2, &Z2, 1, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&D2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 10, 17, 12, 17, 0, "Xijab");
  dpd_buf4_axpy(&Z2, &Xijab, -1.0);
  dpd_buf4_close(&Xijab);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 12, 17, "Xijab", 1.0);
  dpd_buf4_close(&Z2);

  /*  - LR_OO(I,M) <Mj||Ab> + LR_OO(j,m) <mI||Ab> */
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_OO");
  dpd_buf4_init(&D2, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_contract244(&I1, &D2, &XIjAb, 1, 0, 0, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 2, "LR_oo");
  dpd_contract424(&D2, &I1, &XIjAb, 1, 1, 1, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&D2);
  dpd_buf4_close(&XIjAb);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("terms 3 and 10");
#endif

  x_xi2_4_uhf();

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("terms 4 and 6");
#endif

  /* term 7, Xijab += (0.25 Lmnab Rmnef) <ij||ef> */
  dpd_buf4_init(&D2, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  dpd_buf4_init(&R2, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
  dpd_buf4_init(&Z2, EOM_TMP1, R_irr, 2, 2, 2, 2, 0, "Z (I>J,M>N)");
  dpd_contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0); 
  dpd_buf4_close(&R2);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  dpd_buf4_init(&L2, CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
  dpd_contract444(&Z2, &L2, &XIJAB, 0, 1, 1.0, 1.0); 
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&XIJAB);

  dpd_buf4_init(&D2, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_init(&R2, CC_GR, R_irr, 12, 17, 12, 17, 0, "Rijab");
  dpd_buf4_init(&Z2, EOM_TMP1, R_irr, 12, 12, 12, 12, 0, "Z (i>j,m>n)");
  dpd_contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0); 
  dpd_buf4_close(&R2);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 17, 12, 17, 0, "Xijab");
  dpd_buf4_init(&L2, CC_GL, L_irr, 12, 17, 12, 17, 0, "Lijab");
  dpd_contract444(&Z2, &L2, &Xijab, 0, 1, 1.0, 1.0); 
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Xijab);

  dpd_buf4_init(&D2, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_buf4_init(&R2, CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
  dpd_buf4_init(&Z2, EOM_TMP1, R_irr, 22, 22, 22, 22, 0, "Z (Ij,Mn)");
  dpd_contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0); 
  dpd_buf4_close(&R2);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  dpd_buf4_init(&L2, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  dpd_contract444(&Z2, &L2, &XIjAb, 0, 1, 1.0, 1.0); 
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&XIjAb);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

  /* term 8, Xijab += (0.25 Rmnef Lijef) <mn||ab> */
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_OOOO");
  dpd_buf4_init(&D2, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  dpd_contract444(&I2, &D2, &XIJAB, 1, 1, 1.0, 1.0); 
  dpd_buf4_close(&D2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 17, 12, 17, 0, "Xijab");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 12, 12, 12, 12, 0, "R2L2_oooo");
  dpd_buf4_init(&D2, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  dpd_contract444(&I2, &D2, &Xijab, 1, 1, 1.0, 1.0); 
  dpd_buf4_close(&D2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&Xijab);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 22, 22, 22, 22, 0, "R2L2_OoOo");
  dpd_buf4_init(&D2, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
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
  /* build 1-e intermediate in parentheses above */
      /* for term 11: */
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
  dpd_file2_copy(&I1, EOM_TMP1, "Z (F,A)");
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
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

      /* for term 11: */
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 3, 3, "RD_vv");
  dpd_file2_copy(&I1, EOM_TMP1, "Z (f,a)");
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 3, 3, "Z (f,a)");
      /* for term 20: */
  dpd_file2_init(&Z1B, EOM_TMP_XI, R_irr, 3, 3, "R1Wamef_vv");
  dpd_file2_axpy(&Z1B, &I1, -1.0, 0);
  dpd_file2_close(&Z1B);
      /* for term 17: */
  dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
  dpd_file2_init(&F1, CC_OEI, 0, 2, 3, "Fme");
  dpd_contract222(&R1, &F1, &I1, 1, 1, 1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&R1);
  dpd_file2_close(&I1);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (I>J,AB)"); 
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_contract244(&I1, &L2, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
  dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z (i>j,ab)"); 
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 3, 3, "Z (f,a)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 12, 15, 12, 17, 0, "Lijab");
  dpd_contract244(&I1, &L2, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 15, 12, 17, 0, "Xijab");
  dpd_buf4_axpy(&Z2, &Xijab, -1.0);
  dpd_buf4_close(&Xijab);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 12, 17, "Xijab", 1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  dpd_contract244(&I1, &L2, &XIjAb, 0, 2, 1, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 3, 3, "Z (f,a)");
  dpd_contract424(&L2, &I1, &XIjAb, 3, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);
  dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 11, 17 and 20");
#endif

  /* term 12, Xijab -= 0.5 P(ij) Lmjab (Rmnef <in||ef>) */
  /* term 16,       -=     P(ij) Lmjab (Rme Fie) */
  /* term 21,       -=     P(ij) Lmjab (Rne Winme) */
  /* make 1-electron intermediate in parentheses above */
      /* for term 12: */
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  dpd_file2_copy(&I1, EOM_TMP1, "Z (M,I)");
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
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

      /* for term 12 */
  dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 2, 2, "RD_oo");
  dpd_file2_copy(&I1, EOM_TMP1, "Z (m,i)");
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 2, 2, "Z (m,i)");
      /* for term 21 */
  dpd_file2_init(&Z1B, EOM_TMP_XI, R_irr, 2, 2, "R1Wmnie_oo");
  dpd_file2_axpy(&Z1B, &I1, 1.0, 1);
  dpd_file2_close(&Z1B);
      /* for term 16 */
  dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
  dpd_file2_init(&F1, CC_OEI, 0, 2, 3, "Fme");
  dpd_contract222(&R1, &F1, &I1, 0, 0, 1.0, 1.0);
  dpd_file2_close(&F1);
  dpd_file2_close(&R1);
  dpd_file2_close(&I1);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (IJ,A>B)"); 
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 7, 2, 7, 0, "LIJAB");
  dpd_contract244(&I1, &L2, &Z2, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
  dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 17, 10, 17, 0, "Z (ij,a>b)"); 
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 2, 2, "Z (m,i)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 10, 17, 12, 17, 0, "Lijab");
  dpd_contract244(&I1, &L2, &Z2, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 10, 17, 12, 17, 0, "Xijab");
  dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 12, 17, "Xijab", 1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  dpd_contract244(&I1, &L2, &XIjAb, 0, 0, 0, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP1, R_irr, 2, 2, "Z (m,i)");
  dpd_contract424(&L2, &I1, &XIjAb, 1, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&I1);
  dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 12, 16 and 21");
#endif

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

    /* term 13 + 15, (0.25 Rmnef <mn||ef> + Rme Fme) Lijab */
  if ( (L_irr == 0) && (!params.connect_xi)) {
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
    tval = dpd_file2_dot(&R1, &F1);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
    dpd_file2_init(&F1, CC_OEI, 0, 2, 3, "Fme");
    tval += dpd_file2_dot(&R1, &F1);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    tval += params.RD_overlap;

    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_init(&L2, CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_axpy(&L2, &XIJAB, tval);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 17, 12, 17, 0, "Xijab");
    dpd_buf4_init(&L2, CC_GL, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_axpy(&L2, &Xijab, tval);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Xijab);
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
    dpd_buf4_init(&L2, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
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
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

    /* term 22, +P(ij) (Limfe Rme) Wfjab = +P(ij) L2R1(I,F) WFJAB */
  if (!params.connect_xi) {
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (IJ,A>B)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    dpd_contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, 1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", -1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 17, 10, 17, 0, "Z (ij,a>b)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    dpd_contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 10, 17, 12, 17, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &Xijab, 1.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 12, 17, "Xijab", -1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_contract244(&I1, &H2, &XIjAb, 1, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_close(&XIjAb);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 23, 29, 23, 29, 0, "Z (jI,bA)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 22, 28, "XIjAb", 1.0);
    dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 22 (Wamef)");
#endif
    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);
  }

    /* term 23, -P(ab) (Lmnea Rme) Wijnb */
  if (!params.connect_xi) {
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (I>J,AB)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 2, 20, 2, 20, 0, "WMNIE");
    dpd_contract244(&I1, &H2, &Z2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z (i>j,ab)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 12, 30, 12, 30, 0, "Wmnie");
    dpd_contract244(&I1, &H2, &Z2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 12, 15, 12, 17, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 12, 17, "Xijab", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 22, 24, 22, 24, 0, "WMnIe");
    dpd_contract244(&I1, &H2, &XIjAb, 0, 2, 1, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_close(&XIjAb);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 23, 28, 23, 28, 0, "Z (jI,Ab)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
    dpd_contract424(&H2, &I1, &Z2, 3, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 22, 28, "XIjAb", -1.0);
    dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 23 (Wmnie)");
#endif
  }

  /* term 25, Xijab += (Lnmab Rme) Wijne */
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 7, 20, 7, 20, 0, "L2R1_VVOV");
  dpd_buf4_init(&H2, CC_HBAR, 0, 2, 20, 2, 20, 0, "WMNIE");
  dpd_contract444(&H2, &I2, &XIJAB, 0, 0, 1.0, 1.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&XIJAB);

  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 17, 12, 17, 0, "Xijab");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 17, 30, 17, 30, 0, "L2R1_vvov");
  dpd_buf4_init(&H2, CC_HBAR, 0, 12, 30, 12, 30, 0, "Wmnie");
  dpd_contract444(&H2, &I2, &Xijab, 0, 0, 1.0, 1.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&Xijab);

  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 28, 24, 28, 24, 0, "L2R1_VvOv");
  dpd_buf4_init(&H2, CC_HBAR, 0, 22, 24, 22, 24, 0, "WMnIe");
  dpd_contract444(&H2, &I2, &XIjAb, 0, 0, 1.0, 1.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&XIjAb);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 23, 28, 23, 28, 0, "Z (jI,Ab)");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 28, 27, 28, 27, 0, "L2R1_VvoV");
  dpd_buf4_init(&H2, CC_HBAR, 0, 23, 27, 23, 27, 0, "WmNiE");
  dpd_contract444(&H2, &I2, &Z2, 0, 0, 1.0, 0.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 22, 28, "XIjAb", 1.0);
  dpd_buf4_close(&Z2);

#ifdef DEBUG_XI
x_xi_check("term 25 (Wmnie)");
#endif
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

    /* term 24, Xijab -= (Lijfe Rme) Wfmab */
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 21, 2, 21, 0, "L2R1_OOVO");
  dpd_buf4_init(&H2, CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
  dpd_contract444(&I2, &H2, &XIJAB, 0, 1, -1.0, 1.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&XIJAB);

  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 17, 12, 17, 0, "Xijab");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 12, 31, 12, 31, 0, "L2R1_oovo");
  dpd_buf4_init(&H2, CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
  dpd_contract444(&I2, &H2, &Xijab, 0, 1, -1.0, 1.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&Xijab);

  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 22, 26, 22, 26, 0, "L2R1_OoVo");
  dpd_buf4_init(&H2, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
  dpd_contract444(&I2, &H2, &XIjAb, 0, 1, -1.0, 1.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_close(&XIjAb);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 22, 29, 22, 29, 0, "Z (Ij,bA)");
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 22, 25, 22, 25, 0, "L2R1_OovO");
  dpd_buf4_init(&H2, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
  dpd_contract444(&I2, &H2, &Z2, 0, 1, 1.0, 0.0); 
  dpd_buf4_close(&H2);
  dpd_buf4_close(&I2);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 22, 28, "XIjAb", -1.0);
  dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 24 (Wamef)");
#endif
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* terms 18, 19: Xijab -= P(ij) P(ab) Linae (Rme Wmjnb + Rnf Wejbf) */
  /* construct Z(JB,NE) = RME WMJNB + RNF WEJBF */
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 21, 21, 21, 21, 0, "Z (EJ,BN)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 21, 2, 21, 0, "WMNIE (M>N,EI)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 21, 5, 21, 7, 0, "WAMEF");
  dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&R1);
  dpd_buf4_sort(&Z, EOM_TMP1, qpsr, 20, 20, "Z (JE,NB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 20, 20, 20, 20, 0, "Z (JE,NB)");
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 20, 20, "Z (JB,NE)");
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 25, 27, 25, 27, 0, "Z (eJ,nB)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 23, 27, 23, 27, 0, "WmNiE");
  dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
  dpd_contract244(&R1, &H2, &Z, 1, 2, 1, -1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&R1);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 24, 27, "Z (Je,nB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 24, 27, 24, 27, 0, "Z (Je,nB)");
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 20, 30, "Z (JB,ne)");
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 31, 31, 31, 31, 0, "Z (ej,bn)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 10, 31, 12, 31, 0, "Wmnie (m>n,ei)");
  dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 31, 15, 31, 17, 0, "Wamef");
  dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&R1);
  dpd_buf4_sort(&Z, EOM_TMP1, qpsr, 30, 30, "Z (je,nb)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 30, 30, 30, 30, 0, "Z (je,nb)");
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 30, 30, "Z (jb,ne)");
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 26, 24, 26, 24, 0, "Z (Ej,Nb)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 22, 24, 22, 24, 0, "WMnIe");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
  dpd_contract244(&R1, &H2, &Z, 1, 2, 1, -1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_file2_close(&R1);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 27, 24, "Z (jE,Nb)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 27, 24, 27, 24, 0, "Z (jE,Nb)");
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 30, 20, "Z (jb,NE)");
  dpd_buf4_close(&Z);

  /* construct Z(Jb,Ne) <= Z(Je,bN) = - Rme WJmNb + WeJbF RNF */
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 24, 25, 24, 25, 0, "Z (Je,bN)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
  dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract424(&H2, &R1, &Z, 1, 0, 1, -1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 25, 25, 25, 25, 0, "Z (eJ,bN)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_sort_axpy(&Z, EOM_TMP1, qprs, 24, 25, "Z (Je,bN)", 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 24, 25, 24, 25, 0, "Z (Je,bN)");
  dpd_buf4_sort(&Z, EOM_TMP1, prsq, 24, 24, "Z (Jb,Ne)");
  dpd_buf4_close(&Z);

  /* construct Z(jB,nE) <= Z(jE,Bn) = - RME WjMnB + WEjBf Rnf */
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 27, 26, 27, 26, 0, "Z (jE,Bn)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&H2, &R1, &Z, 1, 0, 1, -1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 26, 26, 26, 26, 0, "Z (Ej,Bn)");
  dpd_buf4_init(&H2, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
  dpd_file2_init(&R1, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_sort_axpy(&Z, EOM_TMP1, qprs, 27, 26, "Z (jE,Bn)", 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 27, 26, 27, 26, 0, "Z (jE,Bn)");
  dpd_buf4_sort(&Z, EOM_TMP1, prsq, 27, 27, "Z (jB,nE)");
  dpd_buf4_close(&Z);

  /* XIJAB -= P(IJ) P(AB) L(IA,NE) Z(NE,JB) */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 20, 20, 20, 20, 0, "Z2 (IA,JB)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 20, 20, 20, 20, 0, "Z (JB,NE)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 20, 30, 20, 30, 0, "Z (JB,ne)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);
  dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "Z2 (IJ,AB)");
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (IJ,AB)");
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 5, 2, 7, 0, "XIJAB");
  dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 2, 7, "XIJAB", -1.0);
  dpd_buf4_close(&Z2);

  /* Xijab -= P(ij) P(ab) L(ia,ne) Z(ne,jb) */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 30, 30, 30, 30, 0, "Z2 (ia,jb)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 30, 30, 30, 30, 0, "Z (jb,ne)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_GL, L_irr, 30, 20, 30, 20, 0, "LiaJB");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 30, 20, 30, 20, 0, "Z (jb,NE)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);
  dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 10, 15, "Z2 (ij,ab)");
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 15, 10, 15, 0, "Z2 (ij,ab)");
  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 10, 15, 12, 17, 0, "Xijab");
  dpd_buf4_axpy(&Z2, &Xijab, -1.0);
  dpd_buf4_close(&Xijab);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 12, 17, "Xijab", 1.0);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 12, 17, "Xijab", 1.0);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 12, 17, "Xijab", -1.0);
  dpd_buf4_close(&Z2);

  /* XIjAb += - L(IA,NE) Z(jb,NE) - L(IA,ne) Z(jb,ne)  */
  /*          - L(jb,ne) Z(IA,ne) - L(jb,NE) Z(IA,NE)  */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 20, 30, 20, 30, 0, "Z2 (IA,jb)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 30, 20, 30, 20, 0, "Z (jb,NE)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 30, 30, 30, 30, 0, "Z (jb,ne)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 20, 30, 20, 30, 0, "Z (JB,ne)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
  dpd_contract444(&Z, &L2, &Z2, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 20, 20, 20, 20, 0, "Z (JB,NE)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 30, 20, 30, 20, 0, "LiaJB");
  dpd_contract444(&Z, &L2, &Z2, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, prqs, 22, 28, "XIjAb", -1.0);
  dpd_buf4_close(&Z2);

  /* XIjAb += + L(jA,Ne) Z(Ib,Ne) + L(Ib,nE) Z(jA,nE)  */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 24, 27, 24, 27, 0, "Z2 (Ib,jA)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 27, 24, 27, 24, 0, "LjAIb");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 24, 24, 24, 24, 0, "Z (Jb,Ne)");
  dpd_contract444(&Z, &L2, &Z2, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_GL, L_irr, 24, 27, 24, 27, 0, "LIbjA");
  dpd_buf4_init(&Z, EOM_TMP1, R_irr, 27, 27, 27, 27, 0, "Z (jB,nE)");
  dpd_contract444(&L2, &Z, &Z2, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&L2);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, prsq, 22, 28, "XIjAb", 1.0);
  dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("terms 18, 19 (Wmnie, Wamef)");
#endif
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* Write irrep of XI amplitudes to CC_INFO */
  psio_write_entry(CC_INFO, "XI Irrep", (char *) &G_irr,sizeof(int));

  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  tval = dpd_file2_dot_self(&XIA);
  dpd_file2_close(&XIA);
  fprintf(outfile,"XIA amplitudes: norm=%20.15lf dot=%20.15lf\n", sqrt(tval), tval );
  dpd_file2_init(&Xia, EOM_XI, G_irr, 2, 3, "Xia");
  tval += dpd_file2_dot_self(&Xia);
  fprintf(outfile,"X1 amplitudes:  norm=%20.15lf dot=%20.15lf\n", sqrt(tval), tval );
  dpd_file2_close(&Xia);
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  tval += dpd_buf4_dot_self(&XIJAB);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 12, 17, 12, 17, 0, "Xijab");
  tval += dpd_buf4_dot_self(&Xijab);
  dpd_buf4_close(&Xijab);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
  tval += dpd_buf4_dot_self(&XIjAb);
  dpd_buf4_close(&XIjAb);
  fprintf(outfile,"Norm of Xi: %20.15lf\n", sqrt(tval) );
  return;
}



/* compute terms 4 and 6 of xi2 amplitudes */
/* Xijab += P(ij) P(ab) (Rme Lia + Rmnef Linaf) <mj||eb> */
/* Xijab += P(ij) P(ab) Z2(me,ia) <mj||eb>(me,jb) */
void x_xi2_4_uhf(void)
{
  dpdfile2 RIA, Ria, LIA, Lia;
  int L_irr, R_irr, G_irr, nirreps;
  int I, A, M, E, i, a, m, e, h, row, col, Isym, Esym, Asym, Msym;
  dpdbuf4 D, R2, L2, H2, I2, Z, Z2, XIJAB, Xijab, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  dpd_buf4_init(&I2, EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  dpd_buf4_copy(&I2, EOM_TMP1, "RL_OVOV");
  dpd_buf4_close(&I2);
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  dpd_buf4_copy(&I2, EOM_TMP1, "RL_ovov");
  dpd_buf4_close(&I2);
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  dpd_buf4_copy(&I2, EOM_TMP1, "RL_OVov");
  dpd_buf4_close(&I2);
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  dpd_buf4_copy(&I2, EOM_TMP1, "RL_ovOV");
  dpd_buf4_close(&I2);
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
  dpd_buf4_copy(&I2, EOM_TMP1, "RL_oVoV");
  dpd_buf4_close(&I2);
  dpd_buf4_init(&I2, EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
  dpd_buf4_copy(&I2, EOM_TMP1, "RL_OvOv");
  dpd_buf4_close(&I2);

  /* RL_OVOV(me,ia) += Rme Lia */
  dpd_file2_init(&RIA, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&LIA, CC_GL, L_irr, 0, 1, "LIA");
  dpd_file2_init(&Ria, CC_GR, R_irr, 2, 3, "Ria");
  dpd_file2_init(&Lia, CC_GL, L_irr, 2, 3, "Lia");

  dpd_file2_mat_init(&RIA); dpd_file2_mat_init(&Ria);
  dpd_file2_mat_init(&LIA); dpd_file2_mat_init(&Lia);
  dpd_file2_mat_rd(&RIA); dpd_file2_mat_rd(&Ria);
  dpd_file2_mat_rd(&LIA); dpd_file2_mat_rd(&Lia);

  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 20, 20, 20, 20, 0, "RL_OVOV");
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

  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 30, 30, 30, 30, 0, "RL_ovov");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&I2, h);
    dpd_buf4_mat_irrep_rd(&I2, h);
    for(row=0; row < I2.params->rowtot[h]; row++) {
      m = I2.params->roworb[h][row][0];
      e = I2.params->roworb[h][row][1];
      M = Ria.params->rowidx[m]; Msym = Ria.params->psym[m];
      E = Ria.params->colidx[e]; Esym = Ria.params->qsym[e];
      for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
        i = I2.params->colorb[h^G_irr][col][0];
        a = I2.params->colorb[h^G_irr][col][1];
        I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
        A = Lia.params->colidx[a]; Asym = Lia.params->qsym[a];
        if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
          I2.matrix[h][row][col] += Ria.matrix[Msym][M][E] * Lia.matrix[Isym][I][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&I2, h);
    dpd_buf4_mat_irrep_close(&I2, h);
  }
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 20, 30, 20, 30, 0, "RL_OVov");
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
        I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
        A = Lia.params->colidx[a]; Asym = Lia.params->qsym[a];
        if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
          I2.matrix[h][row][col] += RIA.matrix[Msym][M][E] * Lia.matrix[Isym][I][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&I2, h);
    dpd_buf4_mat_irrep_close(&I2, h);
  }
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP1, G_irr, 30, 20, 30, 20, 0, "RL_ovOV");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&I2, h);
    dpd_buf4_mat_irrep_rd(&I2, h);
    for(row=0; row < I2.params->rowtot[h]; row++) {
      m = I2.params->roworb[h][row][0];
      e = I2.params->roworb[h][row][1];
      M = Ria.params->rowidx[m]; Msym = Ria.params->psym[m];
      E = Ria.params->colidx[e]; Esym = Ria.params->qsym[e];
      for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
        i = I2.params->colorb[h^G_irr][col][0];
        a = I2.params->colorb[h^G_irr][col][1];
        I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
        A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
        if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
          I2.matrix[h][row][col] += Ria.matrix[Msym][M][E] * LIA.matrix[Isym][I][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&I2, h);
    dpd_buf4_mat_irrep_close(&I2, h);
  }
  dpd_buf4_close(&I2);

  dpd_file2_mat_close(&RIA); dpd_file2_mat_close(&Ria);
  dpd_file2_mat_close(&LIA); dpd_file2_mat_close(&Lia);
  dpd_file2_close(&RIA); dpd_file2_close(&Ria);
  dpd_file2_close(&LIA); dpd_file2_close(&Lia);

  /* Z2(IA,JB) = RL(ME,IA) <MJ|EB>(ME,JB) + Z2(me,IA) <mJ|eB>(me,JB) */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 20, 20, 20, 20, 0, "Z2 (IA,JB)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 20, 20, 20, 20, 0, "RL_OVOV");
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 30, 20, 30, 20, 0, "RL_ovOV");
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  dpd_contract444(&Z, &D, &Z2, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);
  dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "Z2 (IJ,AB)");
  dpd_buf4_close(&Z2);

  /* XIJAB += P(IJ) P(AB) Z2(IJ,AB) */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (IJ,AB)");
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 5, 2, 7, 0, "XIJAB");
  dpd_buf4_axpy(&Z2, &XIJAB, 1.0);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", -1.0);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", -1.0);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 2, 7, "XIJAB", 1.0);
  dpd_buf4_close(&Z2);

  /* Z2(ia,jb) = RL(me,ia) <mj|eb>(me,jb) + Z2(ME,ia) <Mj|Eb>(ME,jb) */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 30, 30, 30, 30, 0, "Z2 (ia,jb)");
  dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 30, 30, 30, 30, 0, "RL_ovov");
  dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 20, 30, 20, 30, 0, "RL_OVov");
  dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);
  dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 10, 15, "Z2 (ij,ab)");
  dpd_buf4_close(&Z2);

  /* Xijab += P(ij) P(ab) Z2(ij,ab) */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 15, 10, 15, 0, "Z2 (ij,ab)");
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 10, 15, 12, 17, 0, "Xijab");
  dpd_buf4_axpy(&Z2, &XIJAB, 1.0);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 12, 17, "Xijab", -1.0);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 12, 17, "Xijab", -1.0);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 12, 17, "Xijab", 1.0);
  dpd_buf4_close(&Z2);
    
  /* Z2 (IA,jb) += RL(ME,IA) <Mj|Eb>(ME,jb) + Z2(me,IA) <mj|eb>(me,jb) */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 20, 30, 20, 30, 0, "Z2 (IA,jb)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 20, 20, 20, 20, 0, "RL_OVOV");
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 30, 20, 30, 20, 0, "RL_ovOV");
  dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);

  /* Z2 (IA,jb) +=  RL(ME,jb) <MI||EA>(ME,IA) + <mI|eA>(me,IA) RL(me,jb) */
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 20, 30, 20, 30, 0, "RL_OVov");
  dpd_contract444(&D, &Z, &Z2, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 30, 30, 30, 30, 0, "RL_ovov");
  dpd_contract444(&D, &Z, &Z2, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, prqs, 22, 28, "XIjAb", 1.0);
  dpd_buf4_close(&Z2);

  /* Z2 (Ib,jA) += - <mI||Eb>(mE,Ib) RL(mE,jA) = + <mI|bE>(mE,Ib) RL(mE,jA) */
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 24, 27, 24, 27, 0, "Z2 (Ib,jA)");
  dpd_buf4_init(&D, CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 27, 27, 27, 27, 0, "RL_oVoV");
  dpd_contract444(&D, &Z, &Z2, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);

  /* Z2 (Ib,jA) += - <Mj||eA>(Me,jA) RL(Me,Ia) = + <Mj|Ae>(Me,jA) RL(Me,Ia)*/
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 24, 24, 24, 24, 0, "RL_OvOv");
  dpd_buf4_init(&D, CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
  dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&D);
  dpd_buf4_sort_axpy(&Z2, EOM_XI, prsq, 22, 28, "XIjAb", 1.0);
  dpd_buf4_close(&Z2);

  return;
}

}} // namespace psi::ccdensity
