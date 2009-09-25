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

void x_Gibja_rohf(void);
extern void x_Gibja_uhf(void);

void x_Gibja(void) {
  if (params.ref == 0 || params.ref == 1)
    x_Gibja_rohf();
  else
    x_Gibja_uhf();
  return;
}

/* x_Gibja(): computes non-R0 parts of Gibja 2pdm
   really Dajib = Djabi, then * -1 to get Djaib
   and arranged as G(ia,jb) until final sort
*/

void x_Gibja_rohf(void)
{
  int h, nirreps, row, col, L_irr, R_irr, G_irr;
  int i, j, a, b, I, J, A, B, Isym, Jsym, Asym, Bsym;
  dpdfile2 L1, T1A, T1B, L1A, L1B, R1A, R1B, I1A, I1B;
  dpdbuf4 I2, L2, R2, T2, Z, Z1, V2, G, GIBJA, Gibja, GIbJa, GiBjA, GIbjA, GiBJa;
  double value;
  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* Gajib = lia rjb + limae (rjmbe - rmb tje - rje tmb + rme tjb) */

  /* term 3  Gibja += L(im,ae) R(jm,be) */
  dpd_buf4_init(&V2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_buf4_sort(&V2, EOM_TMP0, rspq, 10, 10, "GIAJB");
  dpd_buf4_close(&V2);
  dpd_buf4_init(&V2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_buf4_sort(&V2, EOM_TMP0, rspq, 10, 10, "Giajb");
  dpd_buf4_close(&V2);
  dpd_buf4_init(&V2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_sort(&V2, EOM_TMP0, rspq, 10, 10, "GIaJb");
  dpd_buf4_close(&V2);
  dpd_buf4_init(&V2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  dpd_buf4_sort(&V2, EOM_TMP0, rspq, 10, 10, "GiAjB");
  dpd_buf4_close(&V2);
  dpd_buf4_init(&V2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_buf4_sort(&V2, EOM_TMP0, rspq, 10, 10, "GIAjb");
  dpd_buf4_close(&V2);
  dpd_buf4_init(&V2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_sort(&V2, EOM_TMP0, rspq, 10, 10, "GiaJB");
  dpd_buf4_close(&V2);

  /* term 4, G(IA,JB) <-- - L(IM,AE) T(J,E) R(M,B) */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(IM,AJ)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&L2, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IB,AJ)");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z, &R1A, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&R1A);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(IA,BJ)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,BJ)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(IA,JB)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(IA,JB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAJB");
  dpd_buf4_axpy(&Z1, &G, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 4, G(ia,jb) <-- - L(im,ae) T(j,e) R(m,b) */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(im,aj)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&L2, &T1B, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ib,aj)");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z, &R1B, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&R1B);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(ia,bj)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ia,bj)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(ia,jb)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(ia,jb)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Giajb");
  dpd_buf4_axpy(&Z1, &G, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 4, G(Ia,Jb) <-- - L(Im,aE) T(J,E) R(m,b) */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(Im,aJ)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjaB");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&L2, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,aJ)");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z, &R1B, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&R1B);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(Ia,bJ)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ia,bJ)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(Ia,Jb)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(Ia,Jb)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIaJb");
  dpd_buf4_axpy(&Z1, &G, 1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 4, G(iA,jB) <-- - L(iM,Ae) T(j,e) R(M,B) */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(iM,Aj)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJAb");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&L2, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iB,Aj)");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z, &R1A, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&R1A);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(iA,Bj)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iA,Bj)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(iA,jB)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(iA,jB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiAjB");
  dpd_buf4_axpy(&Z1, &G, 1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 4, G(IA,jb) <-- - L2(Im,Ae) T(j,e) R(m,b) */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(Im,Aj)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&L2, &T1B, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&T1B);
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
  dpd_contract424(&Z, &R1B, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_file2_close(&R1B);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(IA,bj)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,bj)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(IA,jb)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(IA,jb)");
  dpd_buf4_init(&G, EOM_TMP0, 0, 10, 10, 10, 10, 0, "GIAjb");
  dpd_buf4_axpy(&Z1, &G, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 4, G(ia,JB) <-- - L(iM,aE) T(J,E) R(M,B) */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(iM,aJ)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&L2, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&T1A);
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iB,aJ)");
  dpd_contract424(&Z, &R1A, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_file2_close(&R1A);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(ia,BJ)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ia,BJ)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(ia,JB)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(ia,JB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiaJB");
  dpd_buf4_axpy(&Z1, &G, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1, 0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  /* term 5, G(IA,JB) <-- - (L(IM,AE)*R(J,E)) * T(M,B) */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 11, 2, 11, 0, "L2R1_OOVO");
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IB,AJ)");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(IA,BJ)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,BJ)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(IA,JB)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(IA,JB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAJB");
  dpd_buf4_axpy(&Z1, &G, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 5, G(ia,jb) <-- - (L(im,ae)*R(j,e)) * T(m,b) */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 11, 2, 11, 0, "L2R1_oovo");
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ib,aj)");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1B, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(ia,bj)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ia,bj)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(ia,jb)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(ia,jb)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Giajb");
  dpd_buf4_axpy(&Z1, &G, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 5, G(Ia,Jb) <-- - (L2(Im,Ea)*R(J,E)) * T(m,b) */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO");
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,aJ)");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1B, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(Ia,bJ)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ia,bJ)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(Ia,Jb)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(Ia,Jb)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIaJb");
  dpd_buf4_axpy(&Z1, &G, 1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 5, G(iA,jB) <-- - (L2(iM,eA)*R(j,e)) * T(M,B) */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 11, 0, 11, 0, "Z(iM,Aj)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJAb");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&L2, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1A);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iB,Aj)");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1B, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(iA,Bj)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iA,Bj)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(iA,jB)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(iA,jB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiAjB");
  dpd_buf4_axpy(&Z1, &G, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z1);

  /* term 5, G(IA,jb) <-- - L2(Im,Ae) R(j,e) T(m,b) */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 11, 0, 11, 0, "Z(Im,Aj)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&L2, &R1B, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&R1B);
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
  dpd_contract424(&Z, &T1B, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_file2_close(&T1B);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(IA,bj)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,bj)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(IA,jb)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(IA,jb)");
  dpd_buf4_init(&G, EOM_TMP0, 0, 10, 10, 10, 10, 0, "GIAjb");
  dpd_buf4_axpy(&Z1, &G, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  /* term 5, G(ia,JB) <-- - L(iM,aE) R(J,E) T(M,B) */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 11, 0, 11, 0, "Z(iM,aJ)");
  dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&L2, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_file2_close(&R1A);
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iB,aJ)");
  dpd_contract424(&Z, &T1A, &Z1, 1, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&Z);
  dpd_file2_close(&T1A);
  dpd_buf4_sort(&Z1, EOM_TMP1, prqs, 10, 11, "Z(ia,BJ)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ia,BJ)");
  dpd_buf4_sort(&Z1, EOM_TMP1, pqsr, 10, 10, "Z(ia,JB)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(ia,JB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiaJB");
  dpd_buf4_axpy(&Z1, &G, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1, 0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_mat_init(&R1A);
  dpd_file2_mat_rd(&R1A);
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_file2_mat_init(&R1B);
  dpd_file2_mat_rd(&R1B);
  dpd_file2_init(&L1A, CC_GL, L_irr, 0, 1, "LIA");
  dpd_file2_mat_init(&L1A);
  dpd_file2_mat_rd(&L1A);
  dpd_file2_init(&L1B, CC_GL, L_irr, 0, 1, "Lia");
  dpd_file2_mat_init(&L1B);
  dpd_file2_mat_rd(&L1B);
  if (!params.connect_xi) {
    dpd_file2_init(&I1A, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_mat_init(&I1A);
    dpd_file2_mat_rd(&I1A);
    dpd_file2_init(&I1B, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_file2_mat_init(&I1B);
    dpd_file2_mat_rd(&I1B);
  }
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1A);
  dpd_file2_mat_rd(&T1A);
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_mat_init(&T1B);
  dpd_file2_mat_rd(&T1B);

  /* term 1, G(IA,JB) <-- +  L(I,A) R(J,B) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAJB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h); 
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      I = L1A.params->rowidx[i]; Isym = L1A.params->psym[i];
      a = G.params->roworb[h][row][1];
      A = L1A.params->colidx[a]; Asym = L1A.params->qsym[a];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        j = G.params->colorb[h^G_irr][col][0];
        J = R1A.params->rowidx[j]; Jsym = R1A.params->psym[j];
        b = G.params->colorb[h^G_irr][col][1];
        B = R1A.params->colidx[b]; Bsym = R1A.params->qsym[b];
        if( ((Isym^Asym)==L_irr) && ((Jsym^Bsym)==R_irr))
          G.matrix[h][row][col] += L1A.matrix[Isym][I][A] * R1A.matrix[Jsym][J][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  /* term 2, G(IA,JB) <-- + (L(IM,AE) R(M,E))  T(J,B) = L2R1_OV(I,A) * T(J,B) */
  if (!params.connect_xi) {
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h); 
      dpd_buf4_mat_irrep_rd(&G, h);
      for(row=0; row < G.params->rowtot[h]; row++) {
        i = G.params->roworb[h][row][0];
        I = I1A.params->rowidx[i]; Isym = I1A.params->psym[i];
        a = G.params->roworb[h][row][1];
        A = I1A.params->colidx[a]; Asym = I1A.params->qsym[a];
        for(col=0; col < G.params->coltot[h^G_irr]; col++) {
          j = G.params->colorb[h^G_irr][col][0];
          J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
          b = G.params->colorb[h^G_irr][col][1];
          B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
          if( ((Isym^Asym)==G_irr) && (Jsym==Bsym) )
            G.matrix[h][row][col] += I1A.matrix[Isym][I][A] * T1A.matrix[Jsym][J][B];
	    }
      }
      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
  }
  dpd_buf4_scm(&G, -1.0);
  dpd_buf4_close(&G);

  /* term 1, G(ia,jb) <-- +  L(i,a) R(j,b) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Giajb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h); 
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      I = L1B.params->rowidx[i]; Isym = L1B.params->psym[i];
      a = G.params->roworb[h][row][1];
      A = L1B.params->colidx[a]; Asym = L1B.params->qsym[a];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        j = G.params->colorb[h^G_irr][col][0];
        J = R1B.params->rowidx[j]; Jsym = R1B.params->psym[j];
        b = G.params->colorb[h^G_irr][col][1];
        B = R1B.params->colidx[b]; Bsym = R1B.params->qsym[b];
        if( ((Isym^Asym)==L_irr) && ((Jsym^Bsym)==R_irr))
          G.matrix[h][row][col] += L1B.matrix[Isym][I][A] * R1B.matrix[Jsym][J][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  /* term 2, G(ia,jb) <-- + (L(im,ae) R(m,e))*T(j,b) = L2R1_ov(i,a) * T(j,b) */
  if (!params.connect_xi) {
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h); 
      dpd_buf4_mat_irrep_rd(&G, h);
      for(row=0; row < G.params->rowtot[h]; row++) {
        i = G.params->roworb[h][row][0];
        I = I1B.params->rowidx[i]; Isym = I1B.params->psym[i];
        a = G.params->roworb[h][row][1];
        A = I1B.params->colidx[a]; Asym = I1B.params->qsym[a];
        for(col=0; col < G.params->coltot[h^G_irr]; col++) {
          j = G.params->colorb[h^G_irr][col][0];
          J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
          b = G.params->colorb[h^G_irr][col][1];
          B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
          if( ((Isym^Asym)==G_irr) && (Jsym==Bsym))
            G.matrix[h][row][col] += I1B.matrix[Isym][I][A] * T1B.matrix[Jsym][J][B];
	    }
      }
      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
  }
  dpd_buf4_scm(&G, -1.0);
  dpd_buf4_close(&G);

    /* term 1, G(IA,jb) <-- L(I,A) R(j,b) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAjb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      I = L1A.params->rowidx[i]; Isym = L1A.params->psym[i];
      a = G.params->roworb[h][row][1];
      A = L1A.params->colidx[a]; Asym = L1A.params->qsym[a];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        j = G.params->colorb[h^G_irr][col][0];
        J = R1B.params->rowidx[j]; Jsym = R1B.params->psym[j];
        b = G.params->colorb[h^G_irr][col][1];
        B = R1B.params->colidx[b]; Bsym = R1B.params->qsym[b];
        if( ((Isym^Asym)==L_irr) && ((Jsym^Bsym)==R_irr) )
          G.matrix[h][row][col] += L1A.matrix[Isym][I][A] * R1B.matrix[Jsym][J][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  /* term 2, G(IA,jb) <-- L2R1_OV(I,A) *  T(j,b) */
  if (!params.connect_xi) {
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);
      for(row=0; row < G.params->rowtot[h]; row++) {
        i = G.params->roworb[h][row][0];
        I = I1A.params->rowidx[i]; Isym = I1A.params->psym[i];
        a = G.params->roworb[h][row][1];
        A = I1A.params->colidx[a]; Asym = I1A.params->qsym[a];
        for(col=0; col < G.params->coltot[h^G_irr]; col++) {
          j = G.params->colorb[h^G_irr][col][0];
          J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
          b = G.params->colorb[h^G_irr][col][1];
          B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
          if( ((Isym^Asym)==G_irr) && (Jsym==Bsym))
            G.matrix[h][row][col] += I1A.matrix[Isym][I][A] * T1B.matrix[Jsym][J][B];
        }
      }
      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
  }
  dpd_buf4_scm(&G, -1.0);
  dpd_buf4_close(&G);

  /* term 1, G(ia,JB) <-- L(i,a) R(J,B) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiaJB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      I = L1B.params->rowidx[i]; Isym = L1B.params->psym[i];
      a = G.params->roworb[h][row][1];
      A = L1B.params->colidx[a]; Asym = L1B.params->qsym[a];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        j = G.params->colorb[h^G_irr][col][0];
        J = R1A.params->rowidx[j]; Jsym = R1A.params->psym[j];
        b = G.params->colorb[h^G_irr][col][1];
        B = R1A.params->colidx[b]; Bsym = R1A.params->qsym[b];
        if( ((Isym^Asym)==L_irr) && ((Jsym^Bsym)==R_irr))
          G.matrix[h][row][col] += L1B.matrix[Isym][I][A] * R1A.matrix[Jsym][J][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  /* term 2, G(ia,JB) <-- L2R1_ov(i,a) T(J,B) */
  if (!params.connect_xi) {
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);
      for(row=0; row < G.params->rowtot[h]; row++) {
        i = G.params->roworb[h][row][0];
        I = I1B.params->rowidx[i]; Isym = I1B.params->psym[i];
        a = G.params->roworb[h][row][1];
        A = I1B.params->colidx[a]; Asym = I1B.params->qsym[a];
        for(col=0; col < G.params->coltot[h^G_irr]; col++) {
          j = G.params->colorb[h^G_irr][col][0];
          J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
          b = G.params->colorb[h^G_irr][col][1];
          B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
          if( ((Isym^Asym)==G_irr) && (Jsym==Bsym))
            G.matrix[h][row][col] += I1B.matrix[Isym][I][A] * T1A.matrix[Jsym][J][B];
        }
      }
      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
  }
  dpd_buf4_scm(&G, -1.0);
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&R1A);
  dpd_file2_close(&R1A);
  dpd_file2_mat_close(&R1B);
  dpd_file2_close(&R1B);
  dpd_file2_mat_close(&L1A);
  dpd_file2_close(&L1A);
  dpd_file2_mat_close(&L1B);
  dpd_file2_close(&L1B);
  if (!params.connect_xi) {
    dpd_file2_mat_close(&I1A);
    dpd_file2_close(&I1A);
    dpd_file2_mat_close(&I1B);
    dpd_file2_close(&I1B);
  }
  dpd_file2_mat_close(&T1A);
  dpd_file2_close(&T1A);
  dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1B);

    /* Sort all spin cases to correct ordering (ib,ja) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAJB");
  dpd_buf4_sort(&G, EOM_TMP0, psrq, 10, 10, "GIBJA");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Giajb");
  dpd_buf4_sort(&G, EOM_TMP0, psrq, 10, 10, "Gibja");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIaJb");
  dpd_buf4_scm(&G,-1.0);
  dpd_buf4_sort(&G, EOM_TMP0, psrq, 10, 10, "GIbJa");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiAjB");
  dpd_buf4_scm(&G,-1.0);
  dpd_buf4_sort(&G, EOM_TMP0, psrq, 10, 10, "GiBjA");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAjb");
  dpd_buf4_sort(&G, EOM_TMP0, psrq, 10, 10, "GIbjA");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiaJB");
  dpd_buf4_sort(&G, EOM_TMP0, psrq, 10, 10, "GiBJa");
  dpd_buf4_close(&G);

  /* Add to ground state terms in CC_GAMMA */
  dpd_buf4_init(&GIBJA, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_axpy(&GIBJA, &G, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&GIBJA);
  dpd_buf4_init(&Gibja, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_axpy(&Gibja, &G, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Gibja);
  dpd_buf4_init(&GIbJa, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  dpd_buf4_axpy(&GIbJa, &G, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&GIbJa);
  dpd_buf4_init(&GiBjA, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  dpd_buf4_axpy(&GiBjA, &G, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&GiBjA);
  dpd_buf4_init(&GIbjA, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  dpd_buf4_axpy(&GIbjA, &G, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&GIbjA);
  dpd_buf4_init(&GiBJa, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  dpd_buf4_axpy(&GiBJa, &G, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&GiBJa);

  /* symmetrize after adding to CC_GAMMA */
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_symm(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_symm(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  dpd_buf4_symm(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  dpd_buf4_symm(&G);
  dpd_buf4_close(&G);
  dpd_buf4_init(&GIbjA, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  dpd_buf4_init(&GiBJa, CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  dpd_buf4_symm2(&GIbjA, &GiBJa);
  dpd_buf4_close(&GiBJa);
  dpd_buf4_sort(&GIbjA, CC_GAMMA, rspq, 10, 10, "GiBJa");
  dpd_buf4_close(&GIbjA);

  psio_close(EOM_TMP0, 0);
  psio_open(EOM_TMP0,PSIO_OPEN_NEW);
  psio_close(EOM_TMP1, 0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);
  return;
}

}} // namespace psi::ccdensity
