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

/*
** Gibja(): Constructs the ibja block of the two-electron density
** matrix, defined based on the energy contribution:
**
**  E(TWO) <-- sum_ibja G(ib,ja) <ja||ib>
**
** or by spin-case:
**
** E(AA) <-- G(IB,JA) <JA||IB>
** E(BB) <-- G(ib,ja) <ja||ib>
** E(AB) <-- G(Ib,Ja) <Ja|Ib> + G(iB,jA) <jA|iB> - G(Ib,jA) <jA|bI> - G(iB,Ja) <Ja|Bi>
**
** The (spin-orbital) equation for G(ib,ja) is:
**
** G(ib,ja) = -L(i,a) T(j,b) - L(im,ae) [ T(jm,be) - T(j,e) T(m,b) ]
**
** I actually build the negative of all the terms above for each spin
** case and then multiply each of these by -1.  Hence the dpd_buf4_scm()
** calls you see in the code.
**
** Each Gibja spin-case is built here with the storage G(ia,jb) in
** CC_MISC, but finally resorted to a "proper" G(ib,ja) ordering at
** the end of this routine and placed in CC_GAMMA.  In addition, all
** blocks are "bra-ket symmetrized" (as required by the
** backtransformation) after all contributions have been accounted
** for.  */

void Gibja(void)
{
  int h, nirreps, row, col;
  int i, j, a, b, I, J, A, B, Isym, Jsym, Asym, Bsym;
  dpdfile2 T1, L1, T1A, T1B, L1A, L1B;
  dpdbuf4 G, L, T, Z, Z1, V, G1, G2;
  bool T2_L2_V = true;

  /*  T2 * L2 * V is absent in CC2 Lagrangian */
  if (params.wfn == "CC2" && params.dertype ==1) T2_L2_V = false;

  nirreps = moinfo.nirreps;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    if (T2_L2_V){
    /* G(ia,jb) <-- L(im,ae) T(jm,be) */
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIAJB");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 10, 10, "GIAJB");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "Viajb");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 10, 10, "Giajb");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "ViaJB");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 10, 10, "GIAjb");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIAjb");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 10, 10, "GiaJB");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "VIaJb");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 10, 10, "GIaJb");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "ViAjB");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 10, 10, "GiAjB");
    global_dpd_->buf4_close(&V);
    }

    /* G(IA,JB) <-- - L(IM,AE) T(J,E) T(M,B) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(IM,AJ)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(IB,AJ)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 10, 11, "Z(IA,BJ)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(IA,BJ)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 10, 10, "Z(IA,JB)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(IA,JB)");
    if (T2_L2_V){
       global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GIAJB");
       global_dpd_->buf4_axpy(&Z1, &G, -1.0);
       global_dpd_->buf4_close(&G);
    }
    else
    global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GIAJB", -1.0);
      global_dpd_->buf4_close(&Z1);

    /* G(ia,jb) <-- - L(im,ae) T(j,e) T(m,b) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(im,aj)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 2, 7, 0, "Lijab");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(ib,aj)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 10, 11, "Z(ia,bj)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(ia,bj)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 10, 10, "Z(ia,jb)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(ia,jb)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "Giajb");
      global_dpd_->buf4_axpy(&Z1, &G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "Giajb", -1.0);
    global_dpd_->buf4_close(&Z1);


    /* G(IA,jb) <-- - L(Im,Ae) T(j,e) T(m,b) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(Im,Aj)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 10, 11, "Z(IA,bj)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(IA,bj)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 10, 10, "Z(IA,jb)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(IA,jb)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GIAjb");
      global_dpd_->buf4_axpy(&Z1, &G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GIAjb", -1.0);
    global_dpd_->buf4_close(&Z1);

    /* G(ia,JB) <-- - L(iM,aE) T(J,E) T(M,B) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(iM,aJ)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(iB,aJ)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 10, 11, "Z(ia,BJ)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(ia,BJ)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 10, 10, "Z(ia,JB)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(ia,JB)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GiaJB");
      global_dpd_->buf4_axpy(&Z1, &G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GiaJB", -1.0);
    global_dpd_->buf4_close(&Z1);

    /* G(Ia,Jb) <-- - L(Im,Ea) T(J,E) T(m,b) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(Im,Ja)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(Ib,Ja)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, psrq, 10, 10, "Z(Ia,Jb)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Ia,Jb)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GIaJb");
      global_dpd_->buf4_axpy(&Z1, &G, 1.0);
      global_dpd_->buf4_scm(&G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GIaJb", -1.0);
    global_dpd_->buf4_close(&Z1);

    /* G(iA,jB) <-- - L(iM,eA) T(j,e) T(M,B) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(iM,jA)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract244(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(iB,jA)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, psrq, 10, 10, "Z(iA,jB)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(iA,jB)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GiAjB");
      global_dpd_->buf4_axpy(&Z1, &G, 1.0);
      global_dpd_->buf4_scm(&G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GiAjB", -1.0);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1A);
    global_dpd_->file2_mat_rd(&T1A);
    global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_mat_init(&T1B);
    global_dpd_->file2_mat_rd(&T1B);
    global_dpd_->file2_init(&L1A, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->file2_mat_init(&L1A);
    global_dpd_->file2_mat_rd(&L1A);
    global_dpd_->file2_init(&L1B, PSIF_CC_GLG, 0, 0, 1, "Lia");
    global_dpd_->file2_mat_init(&L1B);
    global_dpd_->file2_mat_rd(&L1B);

    /* G(IA,JB) <-- L(I,A) T(J,B) */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GIAJB");
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    i = G.params->roworb[h][row][0];
    I = L1A.params->rowidx[i]; Isym = L1A.params->psym[i];
    a = G.params->roworb[h][row][1];
    A = L1A.params->colidx[a]; Asym = L1A.params->qsym[a];

    for(col=0; col < G.params->coltot[h]; col++) {
      j = G.params->colorb[h][col][0];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      b = G.params->colorb[h][col][1];
      B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];

      if((Isym==Asym) && (Jsym==Bsym))
        G.matrix[h][row][col] += L1A.matrix[Isym][I][A] *
          T1A.matrix[Jsym][J][B];

    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);


    /* G(ia,jb) <-- L(i,a) T(j,b) */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "Giajb");
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    i = G.params->roworb[h][row][0];
    I = L1B.params->rowidx[i]; Isym = L1B.params->psym[i];
    a = G.params->roworb[h][row][1];
    A = L1B.params->colidx[a]; Asym = L1B.params->qsym[a];

    for(col=0; col < G.params->coltot[h]; col++) {
      j = G.params->colorb[h][col][0];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      b = G.params->colorb[h][col][1];
      B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];

      if((Isym==Asym) && (Jsym==Bsym))
        G.matrix[h][row][col] += L1B.matrix[Isym][I][A] *
          T1B.matrix[Jsym][J][B];
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);

    /* G(IA,jb) <-- L(I,A) T(j,b) */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GIAjb");
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    i = G.params->roworb[h][row][0];
    I = L1A.params->rowidx[i]; Isym = L1A.params->psym[i];
    a = G.params->roworb[h][row][1];
    A = L1A.params->colidx[a]; Asym = L1A.params->qsym[a];

    for(col=0; col < G.params->coltot[h]; col++) {
      j = G.params->colorb[h][col][0];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      b = G.params->colorb[h][col][1];
      B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];

      if((Isym==Asym) && (Jsym==Bsym))
        G.matrix[h][row][col] += L1A.matrix[Isym][I][A] *
          T1B.matrix[Jsym][J][B];
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);

    /* G(ia,JB) <-- L(i,a) T(J,B) */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GiaJB");
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    i = G.params->roworb[h][row][0];
    I = L1B.params->rowidx[i]; Isym = L1B.params->psym[i];
    a = G.params->roworb[h][row][1];
    A = L1B.params->colidx[a]; Asym = L1B.params->qsym[a];

    for(col=0; col < G.params->coltot[h]; col++) {
      j = G.params->colorb[h][col][0];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      b = G.params->colorb[h][col][1];
      B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];

      if((Isym==Asym) && (Jsym==Bsym))
        G.matrix[h][row][col] += L1B.matrix[Isym][I][A] *
          T1A.matrix[Jsym][J][B];
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->file2_mat_close(&L1A);
    global_dpd_->file2_close(&L1A);
    global_dpd_->file2_mat_close(&L1B);
    global_dpd_->file2_close(&L1B);
    global_dpd_->file2_mat_close(&T1A);
    global_dpd_->file2_close(&T1A);
    global_dpd_->file2_mat_close(&T1B);
    global_dpd_->file2_close(&T1B);

    /* Sort all spin cases to correct ordering */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GIAJB");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 10, 10, "GIBJA");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "Giajb");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 10, 10, "Gibja");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GIAjb");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 10, 10, "GIbjA");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GiaJB");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 10, 10, "GiBJa");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GIaJb");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 10, 10, "GIbJa");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 10, 10, 10, 10, 0, "GiAjB");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 10, 10, "GiBjA");
    global_dpd_->buf4_close(&G);

    if (params.ground) { /* otherwise, sort in x_Gibja */
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      global_dpd_->buf4_symm(&G);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
      global_dpd_->buf4_symm(&G);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
      global_dpd_->buf4_symm(&G);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
      global_dpd_->buf4_symm(&G);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
      global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
      global_dpd_->buf4_symm2(&G1, &G2);
      global_dpd_->buf4_close(&G2);
      global_dpd_->buf4_sort(&G1, PSIF_CC_GAMMA, rspq, 10, 10, "GiBJa");
      global_dpd_->buf4_close(&G1);
    }
  }
  else if(params.ref == 2) { /** UHF **/

    if (T2_L2_V){
    /* G(ia,jb) <-- L(im,ae) T(jm,be) */
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 20, 20, 20, 20, 0, "VIAJB");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 20, 20, "GIAJB");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 30, 30, 30, 30, 0, "Viajb");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 30, 30, "Giajb");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 30, 20, 30, 20, 0, "ViaJB");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 20, 30, "GIAjb");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 20, 30, 20, 30, 0, "VIAjb");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 30, 20, "GiaJB");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 24, 24, 24, 24, 0, "VIaJb");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 24, 24, "GIaJb");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 27, 27, 27, 27, 0, "ViAjB");
    global_dpd_->buf4_sort(&V, PSIF_CC_MISC, rspq, 27, 27, "GiAjB");
    global_dpd_->buf4_close(&V);
    }

    /* G(IA,JB) <-- - L(IM,AE) T(J,E) T(M,B) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 21, 0, 21, 0, "Z(IM,AJ)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 20, 21, 20, 21, 0, "Z(IB,AJ)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 20, 21, "Z(IA,BJ)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(IA,BJ)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 20, 20, "Z(IA,JB)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 20, 20, 20, 20, 0, "Z(IA,JB)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 20, 20, 20, 20, 0, "GIAJB");
      global_dpd_->buf4_axpy(&Z1, &G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GIAJB", -1.0);
    global_dpd_->buf4_close(&Z1);

    /* G(ia,jb) <-- - L(im,ae) T(j,e) T(m,b) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 31, 10, 31, 0, "Z(im,aj)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 10, 15, 12, 17, 0, "Lijab");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 30, 31, 30, 31, 0, "Z(ib,aj)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 30, 31, "Z(ia,bj)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(ia,bj)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 30, 30, "Z(ia,jb)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 30, 30, 30, 30, 0, "Z(ia,jb)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 30, 30, 30, 30, 0, "Giajb");
      global_dpd_->buf4_axpy(&Z1, &G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "Giajb", -1.0);
    global_dpd_->buf4_close(&Z1);

    /* G(IA,jb) <-- - L(Im,Ae) T(j,e) T(m,b) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 22, 26, 22, 26, 0, "Z(Im,Aj)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 24, 26, 24, 26, 0, "Z(Ib,Aj)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 20, 31, "Z(IA,bj)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 31, 20, 31, 0, "Z(IA,bj)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 20, 30, "Z(IA,jb)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 20, 30, 20, 30, 0, "GIAjb");
      global_dpd_->buf4_axpy(&Z1, &G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GIAjb", -1.0);
    global_dpd_->buf4_close(&Z1);


    /* G(ia,JB) <-- - L(iM,aE) T(J,E) T(M,B) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(iM,aJ)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 27, 25, 27, 25, 0, "Z(iB,aJ)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 30, 21, "Z(ia,BJ)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 30, 21, 30, 21, 0, "Z(ia,BJ)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, pqsr, 30, 20, "Z(ia,JB)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 30, 20, 30, 20, 0, "Z(ia,JB)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 30, 20, 30, 20, 0, "GiaJB");
      global_dpd_->buf4_axpy(&Z1, &G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GiaJB", -1.0);
    global_dpd_->buf4_close(&Z1);

    /* G(Ia,Jb) <-- - L(Im,Ea) T(J,E) T(m,b) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 22, 24, 22, 24, 0, "Z(Im,Ja)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 24, 24, 24, 24, 0, "Z(Ib,Ja)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, psrq, 24, 24, "Z(Ia,Jb)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Ia,Jb)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 24, 24, 24, 24, 0, "GIaJb");
      global_dpd_->buf4_axpy(&Z1, &G, 1.0);
      global_dpd_->buf4_scm(&G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GIaJb", -1.0);
    global_dpd_->buf4_close(&Z1);

    /* G(iA,jB) <-- - L(iM,eA) T(j,e) T(M,B) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 27, 23, 27, 0, "Z(iM,jA)");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract244(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 27, 27, 27, 27, 0, "Z(iB,jA)");
    global_dpd_->contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, psrq, 27, 27, "Z(iA,jB)");
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(iA,jB)");
    if (T2_L2_V){
      global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 27, 27, 27, 27, 0, "GiAjB");
      global_dpd_->buf4_axpy(&Z1, &G, 1.0);
      global_dpd_->buf4_scm(&G, -1.0);
      global_dpd_->buf4_close(&G);
    }
    else
      global_dpd_->buf4_scmcopy(&Z1, PSIF_CC_MISC, "GiAjB", -1.0);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1A);
    global_dpd_->file2_mat_rd(&T1A);
    global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_mat_init(&T1B);
    global_dpd_->file2_mat_rd(&T1B);
    global_dpd_->file2_init(&L1A, PSIF_CC_GLG, 0, 0, 1, "LIA");
    global_dpd_->file2_mat_init(&L1A);
    global_dpd_->file2_mat_rd(&L1A);
    global_dpd_->file2_init(&L1B, PSIF_CC_GLG, 0, 2, 3, "Lia");
    global_dpd_->file2_mat_init(&L1B);
    global_dpd_->file2_mat_rd(&L1B);

    /* G(IA,JB) <-- L(I,A) T(J,B) */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 20, 20, 20, 20, 0, "GIAJB");
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    i = G.params->roworb[h][row][0];
    I = L1A.params->rowidx[i]; Isym = L1A.params->psym[i];
    a = G.params->roworb[h][row][1];
    A = L1A.params->colidx[a]; Asym = L1A.params->qsym[a];

    for(col=0; col < G.params->coltot[h]; col++) {
      j = G.params->colorb[h][col][0];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      b = G.params->colorb[h][col][1];
      B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];

      if((Isym==Asym) && (Jsym==Bsym))
        G.matrix[h][row][col] += L1A.matrix[Isym][I][A] *
          T1A.matrix[Jsym][J][B];

    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);


    /* G(ia,jb) <-- L(i,a) T(j,b) */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 30, 30, 30, 30, 0, "Giajb");
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    i = G.params->roworb[h][row][0];
    I = L1B.params->rowidx[i]; Isym = L1B.params->psym[i];
    a = G.params->roworb[h][row][1];
    A = L1B.params->colidx[a]; Asym = L1B.params->qsym[a];

    for(col=0; col < G.params->coltot[h]; col++) {
      j = G.params->colorb[h][col][0];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      b = G.params->colorb[h][col][1];
      B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];

      if((Isym==Asym) && (Jsym==Bsym))
        G.matrix[h][row][col] += L1B.matrix[Isym][I][A] *
          T1B.matrix[Jsym][J][B];
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);

    /* G(IA,jb) <-- L(I,A) T(j,b) */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 20, 30, 20, 30, 0, "GIAjb");
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    i = G.params->roworb[h][row][0];
    I = L1A.params->rowidx[i]; Isym = L1A.params->psym[i];
    a = G.params->roworb[h][row][1];
    A = L1A.params->colidx[a]; Asym = L1A.params->qsym[a];

    for(col=0; col < G.params->coltot[h]; col++) {
      j = G.params->colorb[h][col][0];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      b = G.params->colorb[h][col][1];
      B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];

      if((Isym==Asym) && (Jsym==Bsym))
        G.matrix[h][row][col] += L1A.matrix[Isym][I][A] *
          T1B.matrix[Jsym][J][B];
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);

    /* G(ia,JB) <-- L(i,a) T(J,B) */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 30, 20, 30, 20, 0, "GiaJB");
    for(h=0; h < nirreps; h++) {

      global_dpd_->buf4_mat_irrep_init(&G, h);
      global_dpd_->buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
    i = G.params->roworb[h][row][0];
    I = L1B.params->rowidx[i]; Isym = L1B.params->psym[i];
    a = G.params->roworb[h][row][1];
    A = L1B.params->colidx[a]; Asym = L1B.params->qsym[a];

    for(col=0; col < G.params->coltot[h]; col++) {
      j = G.params->colorb[h][col][0];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      b = G.params->colorb[h][col][1];
      B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];

      if((Isym==Asym) && (Jsym==Bsym))
        G.matrix[h][row][col] += L1B.matrix[Isym][I][A] *
          T1A.matrix[Jsym][J][B];
    }
      }

      global_dpd_->buf4_mat_irrep_wrt(&G, h);
      global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->file2_mat_close(&L1A);
    global_dpd_->file2_close(&L1A);
    global_dpd_->file2_mat_close(&L1B);
    global_dpd_->file2_close(&L1B);
    global_dpd_->file2_mat_close(&T1A);
    global_dpd_->file2_close(&T1A);
    global_dpd_->file2_mat_close(&T1B);
    global_dpd_->file2_close(&T1B);

    /* Sort all spin cases to correct ordering */
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 20, 20, 20, 20, 0, "GIAJB");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 20, 20, "GIBJA");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 30, 30, 30, 30, 0, "Giajb");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 30, 30, "Gibja");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 20, 30, 20, 30, 0, "GIAjb");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 24, 27, "GIbjA");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 30, 20, 30, 20, 0, "GiaJB");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 27, 24, "GiBJa");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 24, 24, 24, 24, 0, "GIaJb");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 24, 24, "GIbJa");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_MISC, 0, 27, 27, 27, 27, 0, "GiAjB");
    global_dpd_->buf4_sort(&G, PSIF_CC_GAMMA, psrq, 27, 27, "GiBjA");
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
    global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
    global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
    global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
    global_dpd_->buf4_symm(&G);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
    global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 27, 24, 27, 24, 0, "GiBJa");
    global_dpd_->buf4_symm2(&G1, &G2);
    global_dpd_->buf4_close(&G2);
    global_dpd_->buf4_sort(&G1, PSIF_CC_GAMMA, rspq, 27, 24, "GiBJa");
    global_dpd_->buf4_close(&G1);
  }
}

}} // namespace psi::ccdensity
