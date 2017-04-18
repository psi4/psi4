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
#include <strings.h>
#include <string.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void Gciab(void)
{
  int h, nirreps, a, b, c, i, A, B, C, I, Asym, Bsym, Csym, Isym, row, col;
  double value;
  dpdfile2 L1, T1, g;
  dpdbuf4 G, L, T, Z, Z1, Z2, V;
  double factor=0.0;
  bool T2_L2_V = true;

  /*  T2 * L2 * V is absent in CC2 Lagrangian */
  if (params.wfn == "CC2" && params.dertype ==1) T2_L2_V = false;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
    /* t(M,C) L(Mi,Ab) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(M,C) Tau(Mi,Ab) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(i,e) L(Mn,Ce) --> Z(Mn,Ci) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(Mn,Ci)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(Mn,Ci) Tau(Mn,Ab) --> G(Ci,Ab) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
    if(T2_L2_V) 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    else 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_IjAb");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);

    if(T2_L2_V){
    /* - V(iA,mC) T(m,b) --> Z(iA,bC) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(iA,bC)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "ViAjB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 5, "Z(iC,bA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(iC,bA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 11, 5, "Z(Ci,bA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 11, 5, 11, 5, 0, "Z(Ci,bA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(Ci,Ab)");
    global_dpd_->buf4_close(&Z);
    /* V(ib,MC) T(M,A) --> Z(ib,AC) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(ib,AC)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "ViaJB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, psrq, 10, 5, "Z(iC,Ab)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "Z(iC,Ab)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(Ci,Ab)");
    global_dpd_->buf4_close(&Z);
    /* Z1(Ci,AB) + Z1(Ci,AB) --> G(Ci,AB) */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(Ci,Ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(Ci,Ab)");
    global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
    global_dpd_->buf4_axpy(&Z2, &G, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&G);

    /* g(C,A) T(i,b) --> G(Ci,Ab) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "GAE");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
      global_dpd_->file2_mat_close(&g);
      global_dpd_->file2_close(&g);
      global_dpd_->file2_mat_close(&T1);
      global_dpd_->file2_close(&T1);
  }
   /* add the T3 contributions to CCSD(T) tpdm calculated in cctriples*/
    if (params.wfn == "CCSD_T"){
        global_dpd_->buf4_init(&V, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "GCiAb(T)");
        global_dpd_->buf4_axpy(&V, &G, 1.0);
        global_dpd_->buf4_close(&V);
    }
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
    /* t(M,C) L(MI,AB) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(M,C) Tau(MI,AB) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(I,E) L(MN,CE) --> Z(MN,CI) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 11, 2, 11, 0, "Z(MN,CI)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(MN,CI) Tau(MN,AB) --> G(CI,AB) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
    if(T2_L2_V) 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    else 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_IJAB");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);
    
    if(T2_L2_V){
    /* - V(IA,MC) T(M,B) --> Z(IA,BC) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(IA,BC)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIAJB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 5, "Z(IC,BA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(IC,BA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 11, 5, "Z(CI,BA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 11, 5, 11, 5, 0, "Z(CI,BA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(CI,AB)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(CI,AB)");
    global_dpd_->buf4_axpy(&Z, &Z1, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 7, 0, "GCIAB");
    global_dpd_->buf4_axpy(&Z1, &G, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&G);

    /* - ( g(C,A) T(I,B) - g(C,B) T(I,A) ) --> G(CI,AB) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "GAE");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      B = g.params->colidx[b];  A = T1.params->colidx[a];
      Bsym = g.params->qsym[b];  Asym = T1.params->qsym[a];

      if((Csym==Bsym) && (Isym==Asym))
        value -= g.matrix[Csym][C][B] * T1.matrix[Isym][I][A];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
      global_dpd_->file2_mat_close(&g);
      global_dpd_->file2_close(&g);
      global_dpd_->file2_mat_close(&T1);
      global_dpd_->file2_close(&T1);
   } 
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/


    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
    /* t(m,c) L(mi,ab) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(m,c) Tau(mi,ab) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(i,e) L(mn,ce) --> Z(mn,ci) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 11, 2, 11, 0, "Z(mn,ci)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(mn,ci) Tau(mn,ab) --> G(ci,ab) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
    if(T2_L2_V) 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    else 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_ijab");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);

    if(T2_L2_V){
    /* - V(ia,mc) T(m,b) --> Z(ia,bc) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(ia,bc)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "Viajb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 5, "Z(ic,ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(ic,ba)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 11, 5, "Z(ci,ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 11, 5, 11, 5, 0, "Z(ci,ba)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(ci,ab)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(ci,ab)");
    global_dpd_->buf4_axpy(&Z, &Z1, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 7, 0, "Gciab");
    global_dpd_->buf4_axpy(&Z1, &G, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&G);

    /* - ( g(c,a) T(i,b) - g(c,b) T(i,a) ) --> G(ci,ab) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "Gae");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      B = g.params->colidx[b];  A = T1.params->colidx[a];
      Bsym = g.params->qsym[b];  Asym = T1.params->qsym[a];

      if((Csym==Bsym) && (Isym==Asym))
        value -= g.matrix[Csym][C][B] * T1.matrix[Isym][I][A];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
      global_dpd_->file2_mat_close(&g);
      global_dpd_->file2_close(&g);
      global_dpd_->file2_mat_close(&T1);
      global_dpd_->file2_close(&T1);
  }
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
    /* t(M,C) L(Mi,Ab) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(M,C) Tau(Mi,Ab) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(i,e) L(Mn,Ce) --> Z(Mn,Ci) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(Mn,Ci)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(Mn,Ci) Tau(Mn,Ab) --> G(Ci,Ab) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
    if(T2_L2_V) 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    else 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_IjAb");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);

    if(T2_L2_V){
    /* - V(iA,mC) T(m,b) --> Z(iA,bC) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(iA,bC)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "ViAjB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 5, "Z(iC,bA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(iC,bA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 11, 5, "Z(Ci,bA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 11, 5, 11, 5, 0, "Z(Ci,bA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(Ci,Ab)");
    global_dpd_->buf4_close(&Z);
    /* V(ib,MC) T(M,A) --> Z(ib,AC) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(ib,AC)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "ViaJB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, psrq, 10, 5, "Z(iC,Ab)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "Z(iC,Ab)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(Ci,Ab)");
    global_dpd_->buf4_close(&Z);
    /* Z1(Ci,AB) + Z1(Ci,AB) --> G(Ci,AB) */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(Ci,Ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(Ci,Ab)");
    global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
    global_dpd_->buf4_axpy(&Z2, &G, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&G);

    /* g(C,A) T(i,b) --> G(Ci,Ab) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "GAE");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);

  }
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/



    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
    /* t(m,c) L(mI,aB) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(m,c) Tau(mI,aB) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(I,E) L(mN,cE) --> Z(mN,cI) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(mN,cI)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(mN,cI) Tau(mN,aB) --> G(cI,aB) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
    if(T2_L2_V) 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
    else 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "t1_iJaB");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);

    if(T2_L2_V) {
    /* - V(Ia,Mc) T(M,B) --> Z(Ia,Bc) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(Ia,Bc)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIaJb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 5, "Z(Ic,Ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(Ic,Ba)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 11, 5, "Z(cI,Ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 11, 5, 11, 5, 0, "Z(cI,Ba)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(cI,aB)");
    global_dpd_->buf4_close(&Z);
    /* V(IB,mc) T(m,a) --> Z(IB,ac) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(IB,ac)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIAjb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, psrq, 10, 5, "Z(Ic,aB)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "Z(Ic,aB)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(cI,aB)");
    global_dpd_->buf4_close(&Z);
    /* Z1(cI,aB) + Z2(cI,aB) --> G(cI,aB) */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(cI,aB)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(cI,aB)");
    global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
    global_dpd_->buf4_axpy(&Z2, &G, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&G);

    /* g(c,a) T(I,B) --> G(cI,aB) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "Gae");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
  }  
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/
  }
  else if(params.ref == 2) { /** UHF **/

    if(params.wfn == "CCSD_T" && params.dertype==1) {
      /* For CCSD(T) gradients, some density contributions are
     calculated in cctriples */
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 7, 20, 7, 0, "GIDAB");
      global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, qprs, 21, 7, "GCIAB");
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 17, 30, 17, 0, "Gidab");
      global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, qprs, 31, 17, "Gciab");
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 28, 24, 28, 0, "GIdAb");
      global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, qpsr, 25, 29, "GcIaB");
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 29, 27, 29, 0, "GiDaB");
      global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, qpsr, 26, 28, "GCiAb");
      global_dpd_->buf4_close(&G);

      factor = 1.0;
    }

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");
    /* t(M,C) L(MI,AB) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, factor);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(M,C) Tau(MI,AB) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(I,E) L(MN,CE) --> Z(MN,CI) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 21, 2, 21, 0, "Z(MN,CI)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(MN,CI) Tau(MN,AB) --> G(CI,AB) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");
    if(T2_L2_V)
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    else 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "t1_IJAB");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);

    if(T2_L2_V){
    /* - V(IA,MC) T(M,B) --> Z(IA,BC) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 5, 20, 5, 0, "Z(IA,BC)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 20, 20, 20, 20, 0, "VIAJB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 20, 5, "Z(IC,BA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 20, 5, 20, 5, 0, "Z(IC,BA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 21, 5, "Z(CI,BA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 21, 5, 21, 5, 0, "Z(CI,BA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 21, 5, "Z(CI,AB)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z(CI,AB)");
    global_dpd_->buf4_axpy(&Z, &Z1, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 5, 21, 7, 0, "GCIAB");
    global_dpd_->buf4_axpy(&Z1, &G, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&G);

    /* - ( g(C,A) T(I,B) - g(C,B) T(I,A) ) --> G(CI,AB) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "GAE");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 21, 7, 21, 7, 0, "GCIAB");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      B = g.params->colidx[b];  A = T1.params->colidx[a];
      Bsym = g.params->qsym[b];  Asym = T1.params->qsym[a];

      if((Csym==Bsym) && (Isym==Asym))
        value -= g.matrix[Csym][C][B] * T1.matrix[Isym][I][A];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
      global_dpd_->file2_mat_close(&g);
      global_dpd_->file2_close(&g);
      global_dpd_->file2_mat_close(&T1);
      global_dpd_->file2_close(&T1);
  }
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/


    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");
    /* t(m,c) L(mi,ab) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, factor);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(m,c) Tau(mi,ab) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tauijab");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(i,e) L(mn,ce) --> Z(mn,ci) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 12, 31, 12, 31, 0, "Z(mn,ci)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(mn,ci) Tau(mn,ab) --> G(ci,ab) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");
    if(T2_L2_V) 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    else 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "t1_ijab");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);

    if(T2_L2_V) {
    /* - V(ia,mc) T(m,b) --> Z(ia,bc) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 15, 30, 15, 0, "Z(ia,bc)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 30, 30, 30, 30, 0, "Viajb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 30, 15, "Z(ic,ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 30, 15, 30, 15, 0, "Z(ic,ba)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 31, 15, "Z(ci,ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 31, 15, 31, 15, 0, "Z(ci,ba)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 31, 15, "Z(ci,ab)");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 31, 15, 31, 15, 0, "Z(ci,ab)");
    global_dpd_->buf4_axpy(&Z, &Z1, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 15, 31, 17, 0, "Gciab");
    global_dpd_->buf4_axpy(&Z1, &G, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&G);

    /* - ( g(c,a) T(i,b) - g(c,b) T(i,a) ) --> G(ci,ab) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 3, 3, "Gae");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 31, 17, 31, 17, 0, "Gciab");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      B = g.params->colidx[b];  A = T1.params->colidx[a];
      Bsym = g.params->qsym[b];  Asym = T1.params->qsym[a];

      if((Csym==Bsym) && (Isym==Asym))
        value -= g.matrix[Csym][C][B] * T1.matrix[Isym][I][A];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
      global_dpd_->file2_mat_close(&g);
      global_dpd_->file2_close(&g);
      global_dpd_->file2_mat_close(&T1);
      global_dpd_->file2_close(&T1);
  }
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/


    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
    /* t(M,C) L(Mi,Ab) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, factor);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(M,C) Tau(Mi,Ab) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(i,e) L(Mn,Ce) --> Z(Mn,Ci) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 22, 26, 22, 26, 0, "Z(Mn,Ci)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(Mn,Ci) Tau(Mn,Ab) --> G(Ci,Ab) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
    if(T2_L2_V) 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    else 
      global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "t1_IjAb");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);

    if(T2_L2_V) {
    /* - V(iA,mC) T(m,b) --> Z(iA,bC) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 29, 27, 29, 0, "Z(iA,bC)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 27, 27, 27, 27, 0, "ViAjB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 27, 29, "Z(iC,bA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 27, 29, 27, 29, 0, "Z(iC,bA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 26, 29, "Z(Ci,bA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 26, 29, 26, 29, 0, "Z(Ci,bA)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 26, 28, "Z(Ci,Ab)");
    global_dpd_->buf4_close(&Z);
    /* V(ib,MC) T(M,A) --> Z(ib,AC) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 30, 5, 30, 5, 0, "Z(ib,AC)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 30, 20, 30, 20, 0, "ViaJB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, psrq, 27, 28, "Z(iC,Ab)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 27, 28, 27, 28, 0, "Z(iC,Ab)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 26, 28, "Z(Ci,Ab)");
    global_dpd_->buf4_close(&Z);
    /* Z1(Ci,AB) + Z1(Ci,AB) --> G(Ci,AB) */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 26, 28, 26, 28, 0, "Z(Ci,Ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 26, 28, 26, 28, 0, "Z(Ci,Ab)");
    global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");
    global_dpd_->buf4_axpy(&Z2, &G, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&G);

    /* g(C,A) T(i,b) --> G(Ci,Ab) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "GAE");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 26, 28, 26, 28, 0, "GCiAb");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
      global_dpd_->file2_mat_close(&g);
      global_dpd_->file2_close(&g);
      global_dpd_->file2_mat_close(&T1);
      global_dpd_->file2_close(&T1);
  }
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/



    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
    /* t(m,c) L(mI,aB) */
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract244(&T1, &L, &G, 0, 0, 0, 1.0, factor);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* l(m,c) Tau(mI,aB) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
    global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->contract244(&L1, &T, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&L1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&G);
    /* t(I,E) L(mN,cE) --> Z(mN,cI) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(mN,cI)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&L);
    /* -Z(mN,cI) Tau(mN,aB) --> G(cI,aB) */
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
    if (T2_L2_V) 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
    else 
       global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "t1_iJaB");
    global_dpd_->contract444(&Z, &T, &G, 1, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&G);

    if (T2_L2_V) {
    /* - V(Ia,Mc) T(M,B) --> Z(Ia,Bc) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 28, 24, 28, 0, "Z(Ia,Bc)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 24, 24, 24, 24, 0, "VIaJb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, -1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 24, 28, "Z(Ic,Ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 24, 28, 24, 28, 0, "Z(Ic,Ba)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, qprs, 25, 28, "Z(cI,Ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 25, 28, 25, 28, 0, "Z(cI,Ba)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 25, 29, "Z(cI,aB)");
    global_dpd_->buf4_close(&Z);
    /* V(IB,mc) T(m,a) --> Z(IB,ac) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 20, 15, 20, 15, 0, "Z(IB,ac)");
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 20, 30, 20, 30, 0, "VIAjb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract244(&T1, &V, &Z, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, psrq, 24, 29, "Z(Ic,aB)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 24, 29, 24, 29, 0, "Z(Ic,aB)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 25, 29, "Z(cI,aB)");
    global_dpd_->buf4_close(&Z);
    /* Z1(cI,aB) + Z2(cI,aB) --> G(cI,aB) */
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 25, 29, 25, 29, 0, "Z(cI,aB)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 25, 29, 25, 29, 0, "Z(cI,aB)");
    global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");
    global_dpd_->buf4_axpy(&Z2, &G, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&G);

    /* g(c,a) T(I,B) --> G(cI,aB) */
    global_dpd_->file2_init(&g, PSIF_CC_GLG, 0, 3, 3, "Gae");
    global_dpd_->file2_mat_init(&g);
    global_dpd_->file2_mat_rd(&g);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 25, 29, 25, 29, 0, "GcIaB");

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    c = G.params->roworb[h][row][0];
    i = G.params->roworb[h][row][1];
    for(col=0; col < G.params->coltot[h]; col++) {
      a = G.params->colorb[h][col][0];
      b = G.params->colorb[h][col][1];

      value = 0.0;

      C = g.params->rowidx[c];  I = T1.params->rowidx[i];
      Csym = g.params->psym[c]; Isym = T1.params->psym[i];
      A = g.params->colidx[a];  B = T1.params->colidx[b];
      Asym = g.params->qsym[a];  Bsym = T1.params->qsym[b];

      if((Csym==Asym) && (Isym==Bsym))
        value += g.matrix[Csym][C][A] * T1.matrix[Isym][I][B];

      G.matrix[h][row][col] -= value;
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
      global_dpd_->file2_mat_close(&g);
      global_dpd_->file2_close(&g);
      global_dpd_->file2_mat_close(&T1);
      global_dpd_->file2_close(&T1);
  }
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_close(&G);

    /*global_dpd_->file2_mat_close(&g);
    global_dpd_->file2_close(&g);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);*/

  }
}

}} // namespace psi::ccdensity
