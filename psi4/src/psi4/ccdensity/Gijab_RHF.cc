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
#include <stdio.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void Gijab_RHF(void)
{
  int h, nirreps, i, a, m, e, I, A, M, E, Isym, Asym, Msym, Esym, row, col;
  int j, b, J, B, Jsym, Bsym;
  double value;
  dpdfile2 T1, L1, g, ZZ, ZZ2, T1A, T1B;
  dpdbuf4 G, L, T, V, Z, Z1, Z2;
  bool T2_L2_V = true;

  /*  T2 * L2 * V is absent in CC2 Lagrangian */
  if (params.wfn == "CC2" && params.dertype ==1) T2_L2_V = false;

  nirreps = moinfo.nirreps;

  if (T2_L2_V){
  /* ( g(I,M) + L(M,E) T(I,E) ) --> Z(I,M)(TMP0)  */
  global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 0, 0, "GMI");
  global_dpd_->file2_copy(&g, PSIF_CC_TMP0, "Z(I,M)");
  global_dpd_->file2_close(&g);
  }
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  if (T2_L2_V) 
     global_dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 1.0);
  else  
     global_dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&ZZ);

  if (T2_L2_V){
  /* ( g(i,m) + L(m,e) T(i,e) ) --> Z(i,m)(TMP1)  */
  global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 0, 0, "Gmi");
  global_dpd_->file2_copy(&g, PSIF_CC_TMP1, "Z(i,m)");
  global_dpd_->file2_close(&g);
  }
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP1, 0, 0, 0, "Z(i,m)");
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  if(T2_L2_V) 
    global_dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 1.0);
  else 
    global_dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&ZZ);

  if (T2_L2_V){
  /* ( g(E,A) - L(M,E) T(M,A) ) --> Z(E,A)(TMP2) */
  global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "GAE");
  global_dpd_->file2_copy(&g, PSIF_CC_TMP2, "Z(E,A)");
  global_dpd_->file2_close(&g);
  }
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP2, 0, 1, 1, "Z(E,A)");
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  if (T2_L2_V)
     global_dpd_->contract222(&L1, &T1, &ZZ, 1, 1, -1.0, 1.0);
  else 
     global_dpd_->contract222(&L1, &T1, &ZZ, 1, 1, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&ZZ);

  if (T2_L2_V){
  /* ( g(e,a) - L(m,e) T(m,a) ) --> Z(e,a)(TMP3) */
  global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "Gae");
  global_dpd_->file2_copy(&g, PSIF_CC_TMP3, "Z(e,a)");
  global_dpd_->file2_close(&g);
  }
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP3, 0, 1, 1, "Z(e,a)");
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  if (T2_L2_V) 
     global_dpd_->contract222(&L1, &T1, &ZZ, 1, 1, -1.0, 1.0);
  else 
     global_dpd_->contract222(&L1, &T1, &ZZ, 1, 1, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&ZZ);

  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1A);
  global_dpd_->file2_mat_rd(&T1A);
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_mat_init(&T1B);
  global_dpd_->file2_mat_rd(&T1B);

  /* ( - T(IA,ME) + 2 * T(I,E) T(M,A) ) --> Z(IA,ME) */
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  global_dpd_->buf4_copy(&T, PSIF_CC_TMP4, "Z(IA,ME)");
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 10, 10, 10, 0, "Z(IA,ME)");
  global_dpd_->buf4_scm(&Z, -1.0);
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row=0; row < Z.params->rowtot[h]; row++) {
      i = Z.params->roworb[h][row][0];
      a = Z.params->roworb[h][row][1];
      I = T1A.params->rowidx[i];  Isym = T1A.params->psym[i];
      A = T1A.params->colidx[a];  Asym = T1A.params->qsym[a];

      for(col=0; col < Z.params->coltot[h]; col++) {
    m = Z.params->colorb[h][col][0];
    e = Z.params->colorb[h][col][1];
    M = T1A.params->rowidx[m];  Msym = T1A.params->psym[m];
    E = T1A.params->colidx[e];  Esym = T1A.params->qsym[e];

    if((Isym==Esym) && (Msym==Asym))
      Z.matrix[h][row][col] += (2* T1A.matrix[Isym][I][E] *
                    T1A.matrix[Msym][M][A]);
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->buf4_close(&Z);

  /* ( - T(ia,me) + 2 * T(i,e) T(m,a) ) --> Z(ia,me) */
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  global_dpd_->buf4_copy(&T, PSIF_CC_TMP5, "Z(ia,me)");
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP5, 0, 10, 10, 10, 10, 0, "Z(ia,me)");
  global_dpd_->buf4_scm(&Z, -1.0);
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row=0; row < Z.params->rowtot[h]; row++) {
      i = Z.params->roworb[h][row][0];
      a = Z.params->roworb[h][row][1];
      I = T1B.params->rowidx[i];  Isym = T1B.params->psym[i];
      A = T1B.params->colidx[a];  Asym = T1B.params->qsym[a];

      for(col=0; col < Z.params->coltot[h]; col++) {
    m = Z.params->colorb[h][col][0];
    e = Z.params->colorb[h][col][1];
    M = T1B.params->rowidx[m];  Msym = T1B.params->psym[m];
    E = T1B.params->colidx[e];  Esym = T1B.params->qsym[e];

    if((Isym==Esym) && (Msym==Asym))
      Z.matrix[h][row][col] += (2* T1B.matrix[Isym][I][E] *
                    T1B.matrix[Msym][M][A]);
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->buf4_close(&Z);

  /* ( - T(iA,Me) + 2 * T(i,e) T(M,A) ) --> Z(iA,Me) */
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
  global_dpd_->buf4_copy(&T, PSIF_CC_TMP6, "Z(iA,Me)");
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP6, 0, 10, 10, 10, 10, 0, "Z(iA,Me)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row=0; row < Z.params->rowtot[h]; row++) {
      i = Z.params->roworb[h][row][0];
      a = Z.params->roworb[h][row][1];
      I = T1B.params->rowidx[i];  Isym = T1B.params->psym[i];
      A = T1A.params->colidx[a];  Asym = T1A.params->qsym[a];

      for(col=0; col < Z.params->coltot[h]; col++) {
    m = Z.params->colorb[h][col][0];
    e = Z.params->colorb[h][col][1];
    M = T1A.params->rowidx[m];  Msym = T1A.params->psym[m];
    E = T1B.params->colidx[e];  Esym = T1B.params->qsym[e];

    if((Isym==Esym) && (Msym==Asym))
      Z.matrix[h][row][col] += (2* T1B.matrix[Isym][I][E] *
                    T1A.matrix[Msym][M][A]);
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->buf4_close(&Z);

  /* ( - T(Ia,mE) + 2 * T(I,E) T(m,a) ) --> Z(Ia,mE) */
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  global_dpd_->buf4_copy(&T, PSIF_CC_TMP7, "Z(Ia,mE)");
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP7, 0, 10, 10, 10, 10, 0, "Z(Ia,mE)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row=0; row < Z.params->rowtot[h]; row++) {
      i = Z.params->roworb[h][row][0];
      a = Z.params->roworb[h][row][1];
      I = T1A.params->rowidx[i];  Isym = T1A.params->psym[i];
      A = T1B.params->colidx[a];  Asym = T1B.params->qsym[a];

      for(col=0; col < Z.params->coltot[h]; col++) {
    m = Z.params->colorb[h][col][0];
    e = Z.params->colorb[h][col][1];
    M = T1B.params->rowidx[m];  Msym = T1B.params->psym[m];
    E = T1A.params->colidx[e];  Esym = T1A.params->qsym[e];

    if((Isym==Esym) && (Msym==Asym))
      Z.matrix[h][row][col] += (2* T1A.matrix[Isym][I][E] *
                    T1B.matrix[Msym][M][A]);
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->buf4_close(&Z);

  global_dpd_->file2_mat_close(&T1A);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_mat_close(&T1B);
  global_dpd_->file2_close(&T1B);


  /* L(Ij,Ab) */
  global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->buf4_copy(&L, PSIF_CC_GAMMA, "GIjAb");
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  /* add the T3 contributions to CCSD(T) tpdm calculated in cctriples*/
  if (params.wfn == "CCSD_T"){
      global_dpd_->buf4_init(&V, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "GIjAb(T)");
      global_dpd_->buf4_axpy(&V, &G, 1.0);
      global_dpd_->buf4_close(&V);
    }

  /* Tau(Ij,Ab) * (L0*R0 = 1, ground or 0, excited */
  if (params.ground) {
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_axpy(&T, &G, 1.0);    
    global_dpd_->buf4_close(&T);
  }
  /* V(Ij,Mn) Tau(Mn,Ab) */
  if (T2_L2_V)
     global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  else
     global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_IjAb");
  global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 0, 0, 0, 0, 0, "VMnIj");
  global_dpd_->contract444(&V, &T, &G, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&G);
  /* - ( Z(I,M) Tau(Mj,Ab) - Z(j,m) Tau(mI,bA) ) */
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP8, 0, 0, 5, 0, 5, 0, "Z1(Ij,Ab)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
  global_dpd_->contract244(&ZZ, &T, &Z1, 1, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&ZZ);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 0, 5, 0, 5, 0, "Z2(jI,bA)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP1, 0, 0, 0, "Z(i,m)");
  global_dpd_->contract244(&ZZ, &T, &Z2, 1, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&ZZ);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP10, qprs, 0, 5, "Z2(Ij,bA)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP10, 0, 0, 5, 0, 5, 0, "Z2(Ij,bA)");
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP9, pqsr, 0, 5, "Z2(Ij,Ab)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 0, 5, 0, 5, 0, "Z2(Ij,Ab)");
  global_dpd_->buf4_axpy(&Z2, &Z1, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  global_dpd_->buf4_axpy(&Z1, &G, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);
  /* - ( Z(E,A) Tau(Ij,bE) - Z(e,b) Tau(Ij,Ae) ) */
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP8, 0, 0, 5, 0, 5, 0, "ZZ1(Ij,Ab)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP3, 0, 1, 1, "Z(e,a)");
  global_dpd_->contract424(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&ZZ);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 0, 5, 0, 5, 0, "Z2(jI,bA)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP2, 0, 1, 1, "Z(E,A)");
  global_dpd_->contract424(&T, &ZZ, &Z2, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&ZZ);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP10, qprs, 0, 5, "Z2(Ij,bA)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP10, 0, 0, 5, 0, 5, 0, "Z2(Ij,bA)");
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP9, pqsr, 0, 5, "Z2(Ij,Ab)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 0, 5, 0, 5, 0, "Z2(Ij,Ab)");
  global_dpd_->buf4_axpy(&Z2, &Z1, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  global_dpd_->buf4_axpy(&Z1, &G, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&G);
  
  if(T2_L2_V){
  /* - P(Ij) P(Ab) ( T'(IA,me) (T2) V(jb,me) + T'(IA,ME) (TMP4) V(jb,ME) ) */
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP8, 0, 10, 10, 10, 10, 0, "Z(IA,jb)");
  global_dpd_->buf4_init(&T, PSIF_CC_TMP4, 0, 10, 10, 10, 10, 0, "Z(IA,ME)");
  global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "ViaJB");
  global_dpd_->contract444(&T, &V, &Z, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "Viajb");
  global_dpd_->contract444(&T, &V, &Z, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&T);
  /* T'(jA,Me) V(Ib,Me) */
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP9, 0, 10, 10, 10, 10, 0, "Z(jA,Ib)");
  global_dpd_->buf4_init(&T, PSIF_CC_TMP6, 0, 10, 10, 10, 10, 0, "Z(iA,Me)");
  global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIaJb");
  global_dpd_->contract444(&T, &V, &Z1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP10, rqps, 10, 10, "Z(IA,jb)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP10, 0, 10, 10, 10, 10, 0, "Z(IA,jb)");
  global_dpd_->buf4_axpy(&Z1, &Z, -1.0);
  global_dpd_->buf4_close(&Z1);
  /* T'(Ib,mE) V(jA,mE) */
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP9, 0, 10, 10, 10, 10, 0, "Z(Ib,jA)");
  global_dpd_->buf4_init(&T, PSIF_CC_TMP7, 0, 10, 10, 10, 10, 0, "Z(Ia,mE)");
  global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "ViAjB");
  global_dpd_->contract444(&T, &V, &Z1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP10, psrq, 10, 10, "Z(IA,jb)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP10, 0, 10, 10, 10, 10, 0, "Z(IA,jb)");
  global_dpd_->buf4_axpy(&Z1, &Z, -1.0);
  global_dpd_->buf4_close(&Z1);
  /* T'(jb,ME) (T2) V(IA,ME) + T'(jb,me) (TMP5) V(IA,me) */
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP9, 0, 10, 10, 10, 10, 0, "Z(IA,jb)");
  global_dpd_->buf4_init(&T, PSIF_CC_TMP5, 0, 10, 10, 10, 10, 0, "Z(ia,me)");
  global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIAjb");
  global_dpd_->contract444(&V, &T, &Z1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIAJB");
  global_dpd_->contract444(&V, &T, &Z1, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_axpy(&Z1, &Z, 1.0);
  global_dpd_->buf4_close(&Z1);
  /* - Z(IA,jb) --> G(Ij,Ab) */
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP9, prqs, 0, 5, "Z(Ij,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP9, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  global_dpd_->buf4_axpy(&Z, &G, -0.5);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  }
  /* T'(IA,me) (T2) L(m,e) + T'(IA,ME) (TMP4) L(M,E) --> ZZ(I,A) */
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 1, "ZZ(I,A)");
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
  global_dpd_->buf4_init(&T, PSIF_CC_TMP4, 0, 10, 10, 10, 10, 0, "Z(IA,ME)");
  global_dpd_->contract422(&T, &L1, &ZZ, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  global_dpd_->contract422(&T, &L1, &ZZ, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&ZZ);
  /* - ZZ(I,A) T(j,b) --> G(Ij,Ab) */
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 1, "ZZ(I,A)");
  global_dpd_->file2_mat_init(&ZZ);
  global_dpd_->file2_mat_rd(&ZZ);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];

      for(col=0; col < G.params->coltot[h]; col++) {
    a = G.params->colorb[h][col][0];
    b = G.params->colorb[h][col][1];

    value = 0.0;

    I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
    J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
    A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
    B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

    if((Isym==Asym) && (Jsym==Bsym))
      value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

    G.matrix[h][row][col] -= value;

      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_mat_close(&ZZ);
  global_dpd_->file2_close(&ZZ);

  /* T'(jb,ME) (T2) L(M,E) + T'(jb,me) (TMP5) L(m,e) --> ZZ(j,b) */
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 1, "ZZ(j,b)");
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  global_dpd_->contract422(&T, &L1, &ZZ, 0, 0, -1.0, 0.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
  global_dpd_->buf4_init(&T, PSIF_CC_TMP5, 0, 10, 10, 10, 10, 0, "Z(ia,me)");
  global_dpd_->contract422(&T, &L1, &ZZ, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&ZZ);
  /* - ZZ(j,b) T(I,A) --> G(Ij,Ab) */
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 1, "ZZ(j,b)");
  global_dpd_->file2_mat_init(&ZZ);
  global_dpd_->file2_mat_rd(&ZZ);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];

      for(col=0; col < G.params->coltot[h]; col++) {
    a = G.params->colorb[h][col][0];
    b = G.params->colorb[h][col][1];

    value = 0.0;

    I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
    J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
    A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
    B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

    if((Isym==Asym) && (Jsym==Bsym))
      value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

    G.matrix[h][row][col] -= value;

      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_mat_close(&ZZ);
  global_dpd_->file2_close(&ZZ);

  /* T(j,e) L(m,e) --> ZZ(j,m) */
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 0, "Z(j,m)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
  global_dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&T1);
  /* ZZ(j,m) T(m,b) --> ZZ2(j,b) */
  global_dpd_->file2_init(&ZZ2, PSIF_CC_TMP9, 0, 0, 1, "ZZ2(j,b)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&ZZ);
  global_dpd_->file2_close(&ZZ2);

  /* 3 T(I,A) ZZ(j,b) --> G(Ij,Ab) */
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP9, 0, 0, 1, "ZZ2(j,b)");
  global_dpd_->file2_mat_init(&ZZ);
  global_dpd_->file2_mat_rd(&ZZ);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];

      for(col=0; col < G.params->coltot[h]; col++) {
    a = G.params->colorb[h][col][0];
    b = G.params->colorb[h][col][1];

    value = 0.0;

    I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
    J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
    A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
    B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

    if((Isym==Asym) && (Jsym==Bsym))
      value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

    G.matrix[h][row][col] += 3.0 * value;

      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_mat_close(&ZZ);
  global_dpd_->file2_close(&ZZ);

  /* T(I,E) L(M,E) --> ZZ(I,M) */
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 0, "Z(I,M)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
  global_dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&L1);
  global_dpd_->file2_close(&T1);

  /* ZZ(I,M) T(M,A) --> ZZ2(I,A) */
  global_dpd_->file2_init(&ZZ2, PSIF_CC_TMP9, 0, 0, 1, "ZZ2(I,A)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&ZZ);
  global_dpd_->file2_close(&ZZ2);

  /* 3 T(j,b) ZZ(I,A) --> G(Ij,Ab) */
  global_dpd_->file2_init(&ZZ, PSIF_CC_TMP9, 0, 0, 1, "ZZ2(I,A)");
  global_dpd_->file2_mat_init(&ZZ);
  global_dpd_->file2_mat_rd(&ZZ);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];

      for(col=0; col < G.params->coltot[h]; col++) {
    a = G.params->colorb[h][col][0];
    b = G.params->colorb[h][col][1];

    value = 0.0;

    I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
    J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
    A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
    B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

    if((Isym==Asym) && (Jsym==Bsym))
      value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

    G.matrix[h][row][col] += 3.0 * value;

      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_scm(&G, 0.5);
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_mat_close(&ZZ);
  global_dpd_->file2_close(&ZZ);
}

}} // namespace psi::ccdensity
