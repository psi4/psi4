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
  global_dpd_->buf4_init(&V2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  global_dpd_->buf4_sort(&V2, PSIF_EOM_TMP0, rspq, 10, 10, "GIAJB");
  global_dpd_->buf4_close(&V2);
  global_dpd_->buf4_init(&V2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  global_dpd_->buf4_sort(&V2, PSIF_EOM_TMP0, rspq, 10, 10, "Giajb");
  global_dpd_->buf4_close(&V2);
  global_dpd_->buf4_init(&V2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  global_dpd_->buf4_sort(&V2, PSIF_EOM_TMP0, rspq, 10, 10, "GIaJb");
  global_dpd_->buf4_close(&V2);
  global_dpd_->buf4_init(&V2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  global_dpd_->buf4_sort(&V2, PSIF_EOM_TMP0, rspq, 10, 10, "GiAjB");
  global_dpd_->buf4_close(&V2);
  global_dpd_->buf4_init(&V2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  global_dpd_->buf4_sort(&V2, PSIF_EOM_TMP0, rspq, 10, 10, "GIAjb");
  global_dpd_->buf4_close(&V2);
  global_dpd_->buf4_init(&V2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  global_dpd_->buf4_sort(&V2, PSIF_EOM_TMP0, rspq, 10, 10, "GiaJB");
  global_dpd_->buf4_close(&V2);

  /* term 4, G(IA,JB) <-- - L(IM,AE) T(J,E) R(M,B) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(IM,AJ)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "LIJAB");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&L2, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IB,AJ)");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z, &R1A, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(IA,BJ)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,BJ)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(IA,JB)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(IA,JB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAJB");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 4, G(ia,jb) <-- - L(im,ae) T(j,e) R(m,b) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(im,aj)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 2, 7, 0, "Lijab");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&L2, &T1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ib,aj)");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&Z, &R1B, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(ia,bj)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ia,bj)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(ia,jb)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(ia,jb)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Giajb");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 4, G(Ia,Jb) <-- - L(Im,aE) T(J,E) R(m,b) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(Im,aJ)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjaB");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&L2, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,aJ)");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&Z, &R1B, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(Ia,bJ)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ia,bJ)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(Ia,Jb)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(Ia,Jb)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIaJb");
  global_dpd_->buf4_axpy(&Z1, &G, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 4, G(iA,jB) <-- - L(iM,Ae) T(j,e) R(M,B) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(iM,Aj)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJAb");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&L2, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iB,Aj)");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z, &R1A, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(iA,Bj)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iA,Bj)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(iA,jB)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(iA,jB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiAjB");
  global_dpd_->buf4_axpy(&Z1, &G, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 4, G(IA,jb) <-- - L2(Im,Ae) T(j,e) R(m,b) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(Im,Aj)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&L2, &T1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&T1B);
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
  global_dpd_->contract424(&Z, &R1B, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(IA,bj)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,bj)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(IA,jb)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(IA,jb)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, 0, 10, 10, 10, 10, 0, "GIAjb");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 4, G(ia,JB) <-- - L(iM,aE) T(J,E) R(M,B) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(iM,aJ)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&L2, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iB,aJ)");
  global_dpd_->contract424(&Z, &R1A, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(ia,BJ)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ia,BJ)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(ia,JB)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(ia,JB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiaJB");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* term 5, G(IA,JB) <-- - (L(IM,AE)*R(J,E)) * T(M,B) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 11, 2, 11, 0, "L2R1_OOVO");
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IB,AJ)");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1A, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(IA,BJ)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,BJ)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(IA,JB)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(IA,JB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAJB");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 5, G(ia,jb) <-- - (L(im,ae)*R(j,e)) * T(m,b) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 11, 2, 11, 0, "L2R1_oovo");
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ib,aj)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&Z, &T1B, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(ia,bj)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ia,bj)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(ia,jb)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(ia,jb)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Giajb");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 5, G(Ia,Jb) <-- - (L2(Im,Ea)*R(J,E)) * T(m,b) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO");
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,aJ)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&Z, &T1B, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(Ia,bJ)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ia,bJ)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(Ia,Jb)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(Ia,Jb)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIaJb");
  global_dpd_->buf4_axpy(&Z1, &G, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 5, G(iA,jB) <-- - (L2(iM,eA)*R(j,e)) * T(M,B) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 11, 0, 11, 0, "Z(iM,Aj)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJAb");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&L2, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iB,Aj)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1B, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(iA,Bj)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iA,Bj)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(iA,jB)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(iA,jB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiAjB");
  global_dpd_->buf4_axpy(&Z1, &G, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z1);

  /* term 5, G(IA,jb) <-- - L2(Im,Ae) R(j,e) T(m,b) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 11, 0, 11, 0, "Z(Im,Aj)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&L2, &R1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&R1B);
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
  global_dpd_->contract424(&Z, &T1B, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(IA,bj)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,bj)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(IA,jb)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(IA,jb)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, 0, 10, 10, 10, 10, 0, "GIAjb");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  /* term 5, G(ia,JB) <-- - L(iM,aE) R(J,E) T(M,B) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 11, 0, 11, 0, "Z(iM,aJ)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&L2, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&R1A);
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(iB,aJ)");
  global_dpd_->contract424(&Z, &T1A, &Z1, 1, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, prqs, 10, 11, "Z(ia,BJ)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(ia,BJ)");
  global_dpd_->buf4_sort(&Z1, PSIF_EOM_TMP1, pqsr, 10, 10, "Z(ia,JB)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(ia,JB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiaJB");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_mat_init(&R1A);
  global_dpd_->file2_mat_rd(&R1A);
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->file2_mat_init(&R1B);
  global_dpd_->file2_mat_rd(&R1B);
  global_dpd_->file2_init(&L1A, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  global_dpd_->file2_mat_init(&L1A);
  global_dpd_->file2_mat_rd(&L1A);
  global_dpd_->file2_init(&L1B, PSIF_CC_GL, L_irr, 0, 1, "Lia");
  global_dpd_->file2_mat_init(&L1B);
  global_dpd_->file2_mat_rd(&L1B);
  if (!params.connect_xi) {
    global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_mat_init(&I1A);
    global_dpd_->file2_mat_rd(&I1A);
    global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->file2_mat_init(&I1B);
    global_dpd_->file2_mat_rd(&I1B);
  }
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1A);
  global_dpd_->file2_mat_rd(&T1A);
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_mat_init(&T1B);
  global_dpd_->file2_mat_rd(&T1B);

  /* term 1, G(IA,JB) <-- +  L(I,A) R(J,B) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAJB");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  /* term 2, G(IA,JB) <-- + (L(IM,AE) R(M,E))  T(J,B) = L2R1_OV(I,A) * T(J,B) */
  if (!params.connect_xi) {
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);
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
      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
  }
  global_dpd_->buf4_scm(&G, -1.0);
  global_dpd_->buf4_close(&G);

  /* term 1, G(ia,jb) <-- +  L(i,a) R(j,b) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Giajb");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  /* term 2, G(ia,jb) <-- + (L(im,ae) R(m,e))*T(j,b) = L2R1_ov(i,a) * T(j,b) */
  if (!params.connect_xi) {
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);
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
      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
  }
  global_dpd_->buf4_scm(&G, -1.0);
  global_dpd_->buf4_close(&G);

    /* term 1, G(IA,jb) <-- L(I,A) R(j,b) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAjb");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  /* term 2, G(IA,jb) <-- L2R1_OV(I,A) *  T(j,b) */
  if (!params.connect_xi) {
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);
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
      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
  }
  global_dpd_->buf4_scm(&G, -1.0);
  global_dpd_->buf4_close(&G);

  /* term 1, G(ia,JB) <-- L(i,a) R(J,B) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiaJB");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  /* term 2, G(ia,JB) <-- L2R1_ov(i,a) T(J,B) */
  if (!params.connect_xi) {
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);
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
      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
  }
  global_dpd_->buf4_scm(&G, -1.0);
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&R1A);
  global_dpd_->file2_close(&R1A);
  global_dpd_->file2_mat_close(&R1B);
  global_dpd_->file2_close(&R1B);
  global_dpd_->file2_mat_close(&L1A);
  global_dpd_->file2_close(&L1A);
  global_dpd_->file2_mat_close(&L1B);
  global_dpd_->file2_close(&L1B);
  if (!params.connect_xi) {
    global_dpd_->file2_mat_close(&I1A);
    global_dpd_->file2_close(&I1A);
    global_dpd_->file2_mat_close(&I1B);
    global_dpd_->file2_close(&I1B);
  }
  global_dpd_->file2_mat_close(&T1A);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_mat_close(&T1B);
  global_dpd_->file2_close(&T1B);

    /* Sort all spin cases to correct ordering (ib,ja) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAJB");
  global_dpd_->buf4_sort(&G, PSIF_EOM_TMP0, psrq, 10, 10, "GIBJA");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Giajb");
  global_dpd_->buf4_sort(&G, PSIF_EOM_TMP0, psrq, 10, 10, "Gibja");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIaJb");
  global_dpd_->buf4_scm(&G,-1.0);
  global_dpd_->buf4_sort(&G, PSIF_EOM_TMP0, psrq, 10, 10, "GIbJa");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiAjB");
  global_dpd_->buf4_scm(&G,-1.0);
  global_dpd_->buf4_sort(&G, PSIF_EOM_TMP0, psrq, 10, 10, "GiBjA");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIAjb");
  global_dpd_->buf4_sort(&G, PSIF_EOM_TMP0, psrq, 10, 10, "GIbjA");
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiaJB");
  global_dpd_->buf4_sort(&G, PSIF_EOM_TMP0, psrq, 10, 10, "GiBJa");
  global_dpd_->buf4_close(&G);

  /* Add to ground state terms in CC_GAMMA */
  global_dpd_->buf4_init(&GIBJA, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  global_dpd_->buf4_axpy(&GIBJA, &G, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&GIBJA);
  global_dpd_->buf4_init(&Gibja, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Gibja");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
  global_dpd_->buf4_axpy(&Gibja, &G, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Gibja);
  global_dpd_->buf4_init(&GIbJa, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  global_dpd_->buf4_axpy(&GIbJa, &G, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&GIbJa);
  global_dpd_->buf4_init(&GiBjA, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  global_dpd_->buf4_axpy(&GiBjA, &G, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&GiBjA);
  global_dpd_->buf4_init(&GIbjA, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  global_dpd_->buf4_axpy(&GIbjA, &G, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&GIbjA);
  global_dpd_->buf4_init(&GiBJa, PSIF_EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  global_dpd_->buf4_axpy(&GiBJa, &G, 1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&GiBJa);

  /* symmetrize after adding to CC_GAMMA */
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIBJA");
  global_dpd_->buf4_symm(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "Gibja");
  global_dpd_->buf4_symm(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbJa");
  global_dpd_->buf4_symm(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBjA");
  global_dpd_->buf4_symm(&G);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&GIbjA, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GIbjA");
  global_dpd_->buf4_init(&GiBJa, PSIF_CC_GAMMA, G_irr, 10, 10, 10, 10, 0, "GiBJa");
  global_dpd_->buf4_symm2(&GIbjA, &GiBJa);
  global_dpd_->buf4_close(&GiBJa);
  global_dpd_->buf4_sort(&GIbjA, PSIF_CC_GAMMA, rspq, 10, 10, "GiBJa");
  global_dpd_->buf4_close(&GIbjA);

  psio_close(PSIF_EOM_TMP0, 0);
  psio_open(PSIF_EOM_TMP0,PSIO_OPEN_NEW);
  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);
  return;
}

}} // namespace psi::ccdensity
