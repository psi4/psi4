/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <strings.h>
#include <string.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    void Gijab_UHF(void)
    {
      int h, nirreps, i, a, m, e, I, A, M, E, Isym, Asym, Msym, Esym, row, col;
      int j, b, J, B, Jsym, Bsym;
      double value;
      dpdfile2 T1, L1, g, ZZ, ZZ2, T1A, T1B;
      dpdbuf4 G, L, T, V, Z, Z1, Z2;

      nirreps = moinfo.nirreps;

      /* ( g(I,M) + L(M,E) T(I,E) ) --> Z(I,M)(TMP0)  */
      dpd_->file2_init(&g, PSIF_CC_GLG, 0, 0, 0, "GMI");
      dpd_->file2_copy(&g, PSIF_CC_TMP0, "Z(I,M)");
      dpd_->file2_close(&g);
      dpd_->file2_init(&ZZ, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 1.0);
      dpd_->file2_close(&T1);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&ZZ);

      /* ( g(i,m) + L(m,e) T(i,e) ) --> Z(i,m)(TMP1)  */
      dpd_->file2_init(&g, PSIF_CC_GLG, 0, 2, 2, "Gmi");
      dpd_->file2_copy(&g, PSIF_CC_TMP1, "Z(i,m)");
      dpd_->file2_close(&g);
      dpd_->file2_init(&ZZ, PSIF_CC_TMP1, 0, 2, 2, "Z(i,m)");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 1.0);
      dpd_->file2_close(&T1);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&ZZ);

      /* ( g(E,A) - L(M,E) T(M,A) ) --> Z(E,A)(TMP2) */
      dpd_->file2_init(&g, PSIF_CC_GLG, 0, 1, 1, "GAE");
      dpd_->file2_copy(&g, PSIF_CC_TMP2, "Z(E,A)");
      dpd_->file2_close(&g);
      dpd_->file2_init(&ZZ, PSIF_CC_TMP2, 0, 1, 1, "Z(E,A)");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->contract222(&L1, &T1, &ZZ, 1, 1, -1.0, 1.0);
      dpd_->file2_close(&T1);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&ZZ);

      /* ( g(e,a) - L(m,e) T(m,a) ) --> Z(e,a)(TMP3) */
      dpd_->file2_init(&g, PSIF_CC_GLG, 0, 3, 3, "Gae");
      dpd_->file2_copy(&g, PSIF_CC_TMP3, "Z(e,a)");
      dpd_->file2_close(&g);
      dpd_->file2_init(&ZZ, PSIF_CC_TMP3, 0, 3, 3, "Z(e,a)");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->contract222(&L1, &T1, &ZZ, 1, 1, -1.0, 1.0);
      dpd_->file2_close(&T1);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&ZZ);

      dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_mat_init(&T1A);
      dpd_->file2_mat_rd(&T1A);
      dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->file2_mat_init(&T1B);
      dpd_->file2_mat_rd(&T1B);

      /* ( - T(IA,ME) + 2 * T(I,E) T(M,A) ) --> Z(IA,ME) */
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
      dpd_->buf4_copy(&T, PSIF_CC_TMP4, "Z(IA,ME)");
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_->buf4_scm(&Z, -1.0);
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&Z, h);
	dpd_->buf4_mat_irrep_rd(&Z, h);
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
	dpd_->buf4_mat_irrep_wrt(&Z, h);
	dpd_->buf4_mat_irrep_close(&Z, h);
      }
      dpd_->buf4_close(&Z);

      /* ( - T(ia,me) + 2 * T(i,e) T(m,a) ) --> Z(ia,me) */
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
      dpd_->buf4_copy(&T, PSIF_CC_TMP5, "Z(ia,me)");
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&Z, PSIF_CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_->buf4_scm(&Z, -1.0);
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&Z, h);
	dpd_->buf4_mat_irrep_rd(&Z, h);
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
	dpd_->buf4_mat_irrep_wrt(&Z, h);
	dpd_->buf4_mat_irrep_close(&Z, h);
      }
      dpd_->buf4_close(&Z);

      /* ( - T(iA,Me) + 2 * T(i,e) T(M,A) ) --> Z(iA,Me) */
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
      dpd_->buf4_copy(&T, PSIF_CC_TMP6, "Z(iA,Me)");
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&Z, PSIF_CC_TMP6, 0, 27, 24, 27, 24, 0, "Z(iA,Me)");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&Z, h);
	dpd_->buf4_mat_irrep_rd(&Z, h);
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
	dpd_->buf4_mat_irrep_wrt(&Z, h);
	dpd_->buf4_mat_irrep_close(&Z, h);
      }
      dpd_->buf4_close(&Z);

      /* ( - T(Ia,mE) + 2 * T(I,E) T(m,a) ) --> Z(Ia,mE) */
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
      dpd_->buf4_copy(&T, PSIF_CC_TMP7, "Z(Ia,mE)");
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&Z, PSIF_CC_TMP7, 0, 24, 27, 24, 27, 0, "Z(Ia,mE)");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&Z, h);
	dpd_->buf4_mat_irrep_rd(&Z, h);
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
	dpd_->buf4_mat_irrep_wrt(&Z, h);
	dpd_->buf4_mat_irrep_close(&Z, h);
      }
      dpd_->buf4_close(&Z);

      dpd_->file2_mat_close(&T1A);
      dpd_->file2_close(&T1A);
      dpd_->file2_mat_close(&T1B);
      dpd_->file2_close(&T1B);

      /* L(IJ,AB) */
      if(params.wfn == "CCSD_T" && params.dertype==1) {
	/* For CCSD(T) gradients, some density contributions are
	   calculated in cctriples */
	dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
	dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
	dpd_->buf4_axpy(&L, &G, 1.0);
	dpd_->buf4_close(&L);
	dpd_->buf4_close(&G);
      }
      else {
	dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
	dpd_->buf4_copy(&L, PSIF_CC_GAMMA, "GIJAB");
	dpd_->buf4_close(&L);
      }
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
      /* Tau(IJ,AB) * (L0*R0 = 1, ground or 0, excited */
      if (params.ground) {
	dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
	dpd_->buf4_axpy(&T, &G, 1.0);
	dpd_->buf4_close(&T);
      }
      /* V(IJ,MN) Tau(MN,AB) */
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 2, 2, 2, 2, 0, "VMNIJ");
      dpd_->contract444(&V, &T, &G, 0, 1, 1.0, 1.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_close(&G);
      /* - ( Z(I,M) Tau(MJ,AB) - Z(J,M) Tau(MI,AB) ) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP8, 0, 0, 7, 0, 7, 0, "Z1(IJ,AB)");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
      dpd_->file2_init(&ZZ, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
      dpd_->contract244(&ZZ, &T, &Z1, 1, 0, 0, -1.0, 0.0);
      dpd_->file2_close(&ZZ);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z1, PSIF_CC_TMP9, qprs, 0, 7, "Z2(JI,AB)");
      dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 0, 7, 0, 7, 0, "Z2(JI,AB)");
      dpd_->buf4_axpy(&Z2, &Z1, -1.0);
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
      dpd_->buf4_axpy(&Z1, &G, 1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_close(&G);
      /* - ( Z(E,A) Tau(IJ,BE) - Z(E,B) Tau(IJ,AE) ) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP8, 0, 2, 5, 2, 5, 0, "ZZ1(IJ,AB)");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");
      dpd_->file2_init(&ZZ, PSIF_CC_TMP2, 0, 1, 1, "Z(E,A)");
      dpd_->contract424(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&ZZ);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z1, PSIF_CC_TMP9, pqsr, 2, 5, "Z2(IJ,BA)");
      dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 2, 5, 2, 5, 0, "Z2(IJ,BA)");
      dpd_->buf4_axpy(&Z2, &Z1, -1.0);
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
      dpd_->buf4_axpy(&Z1, &G, 1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_close(&G);
      /* - P(IJ) P(AB) ( T'(IA,ME) (TMP4) V(JB,ME) + T(IA,me) (T2) V(JB,me) ) */
      dpd_->buf4_init(&Z, PSIF_CC_TMP8, 0, 20, 20, 20, 20, 0, "Z(IA,JB)");
      dpd_->buf4_init(&T, PSIF_CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 20, 20, 20, 20, 0, "VIAJB");
      dpd_->contract444(&T, &V, &Z, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 20, 30, 20, 30, 0, "VIAjb");
      dpd_->contract444(&T, &V, &Z, 0, 0, -1.0, 1.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z, PSIF_CC_TMP9, rqps, 20, 20, "Z(JA,IB)");
      dpd_->buf4_sort(&Z, PSIF_CC_TMP10, psrq, 20, 20, "Z(IB,JA)");
      dpd_->buf4_sort(&Z, PSIF_CC_TMP11, rspq, 20, 20, "Z(JB,IA)");
      dpd_->buf4_init(&Z1, PSIF_CC_TMP9, 0, 20, 20, 20, 20, 0, "Z(JA,IB)");
      dpd_->buf4_axpy(&Z1, &Z, -1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_init(&Z1, PSIF_CC_TMP10, 0, 20, 20, 20, 20, 0, "Z(IB,JA)");
      dpd_->buf4_axpy(&Z1, &Z, -1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_init(&Z1, PSIF_CC_TMP11, 0, 20, 20, 20, 20, 0, "Z(JB,IA)");
      dpd_->buf4_axpy(&Z1, &Z, 1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_sort(&Z, PSIF_CC_TMP9, prqs, 0, 5, "Z(IJ,AB)");
      dpd_->buf4_close(&Z);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 2, 7, 0, "GIJAB");
      dpd_->buf4_init(&Z, PSIF_CC_TMP9, 0, 0, 5, 0, 5, 0, "Z(IJ,AB)");
      /* I don't understand this factor of 1/2 that shows up here */
      dpd_->buf4_axpy(&Z, &G, -0.5);
      dpd_->buf4_close(&Z);
      dpd_->buf4_close(&G);
      /* T'(IA,ME) (TMP4) L(M,E) + T'(IA,me) (T2) L(m,e) --> ZZ(I,A) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 1, "ZZ(I,A)");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      dpd_->buf4_init(&T, PSIF_CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_->contract422(&T, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&T);
      dpd_->file2_close(&L1);
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_->contract422(&T, &L1, &ZZ, 0, 0, -1.0, 1.0);
      dpd_->buf4_close(&T);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&ZZ);
      /* - P(IJ) P(AB) ZZ(I,A) T(J,B) --> G(IJ,AB) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 1, "ZZ(I,A)");
      dpd_->file2_mat_init(&ZZ);
      dpd_->file2_mat_rd(&ZZ);
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_mat_init(&T1);
      dpd_->file2_mat_rd(&T1);

      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&G, h);
	dpd_->buf4_mat_irrep_rd(&G, h);
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

	    I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	    J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];

	    if((Jsym==Asym) && (Isym==Bsym))
	      value -= ZZ.matrix[Jsym][J][A] * T1.matrix[Isym][I][B];

	    I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	    J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
	    A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
	    B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

	    if((Isym==Bsym) && (Jsym==Asym))
	      value -= ZZ.matrix[Isym][I][B] * T1.matrix[Jsym][J][A];

	    I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	    J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];

	    if((Isym==Asym) && (Jsym==Bsym))
	      value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

	    G.matrix[h][row][col] -= value;
	      
	  }
	}
	dpd_->buf4_mat_irrep_wrt(&G, h);
	dpd_->buf4_mat_irrep_close(&G, h);
      }
      dpd_->buf4_close(&G);

      dpd_->file2_mat_close(&T1);
      dpd_->file2_close(&T1);
      dpd_->file2_mat_close(&ZZ);
      dpd_->file2_close(&ZZ);

      /* T(J,E) L(M,E) --> ZZ(J,M) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 0, "Z(J,M)");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&T1);
      /* ZZ(J,M) T(M,B) --> ZZ2(J,B) */
      dpd_->file2_init(&ZZ2, PSIF_CC_TMP9, 0, 0, 1, "ZZ2(J,B)");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
      dpd_->file2_close(&T1);
      dpd_->file2_close(&ZZ);
      dpd_->file2_close(&ZZ2);

      /* 3 P(IJ) P(AB) T(I,A) ZZ(J,B) --> G(IJ,AB) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP9, 0, 0, 1, "ZZ2(J,B)");
      dpd_->file2_mat_init(&ZZ);
      dpd_->file2_mat_rd(&ZZ);
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_mat_init(&T1);
      dpd_->file2_mat_rd(&T1);

      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&G, h);
	dpd_->buf4_mat_irrep_rd(&G, h);
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

	    I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	    J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];

	    if((Jsym==Asym) && (Isym==Bsym))
	      value -= T1.matrix[Jsym][J][A] * ZZ.matrix[Isym][I][B];

	    I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	    J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
	    A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
	    B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

	    if((Isym==Bsym) && (Jsym==Asym))
	      value -= T1.matrix[Isym][I][B] * ZZ.matrix[Jsym][J][A];

	    I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	    J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];

	    if((Isym==Asym) && (Jsym==Bsym))
	      value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

	    G.matrix[h][row][col] += 3.0 * value;
	      
	  }
	}
	dpd_->buf4_mat_irrep_wrt(&G, h);
	dpd_->buf4_mat_irrep_close(&G, h);
      }
      dpd_->buf4_scm(&G, 0.5);
      dpd_->buf4_close(&G);

      dpd_->file2_mat_close(&T1);
      dpd_->file2_close(&T1);
      dpd_->file2_mat_close(&ZZ);
      dpd_->file2_close(&ZZ);



      /* L(ij,ab) */
      if(params.wfn == "CCSD_T" && params.dertype==1) {
	/* For CCSD(T) gradients, some density contributions are
	   calculated in cctriples */
	dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
	dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 12, 17, 12, 17, 0, "Lijab");
	dpd_->buf4_axpy(&L, &G, 1.0);
	dpd_->buf4_close(&L);
	dpd_->buf4_close(&G);
      }
      else {
	dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 12, 17, 12, 17, 0, "Lijab");
	dpd_->buf4_copy(&L, PSIF_CC_GAMMA, "Gijab");
	dpd_->buf4_close(&L);
      }
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
      /* Tau(ij,ab) * (L0*R0 = 1, ground or 0, excited */
      if (params.ground) {
	dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
	dpd_->buf4_axpy(&T, &G, 1.0);
	dpd_->buf4_close(&T);
      }
      /* V(ij,mn) Tau(mn,ab) */
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 12, 12, 12, 12, 0, "Vmnij");
      dpd_->contract444(&V, &T, &G, 0, 1, 1.0, 1.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_close(&G);
      /* - ( Z(i,m) Tau(mj,ab) - Z(j,m) Tau(mi,ab) ) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP8, 0, 10, 17, 10, 17, 0, "Z1(ij,ab)");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tauijab");
      dpd_->file2_init(&ZZ, PSIF_CC_TMP1, 0, 2, 2, "Z(i,m)");
      dpd_->contract244(&ZZ, &T, &Z1, 1, 0, 0, -1.0, 0.0);
      dpd_->file2_close(&ZZ);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z1, PSIF_CC_TMP9, qprs, 10, 17, "Z2(ji,ab)");
      dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 10, 17, 10, 17, 0, "Z2(ji,ab)");
      dpd_->buf4_axpy(&Z2, &Z1, -1.0);
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 17, 12, 17, 0, "Gijab");
      dpd_->buf4_axpy(&Z1, &G, 1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_close(&G);
      /* - ( Z(e,a) Tau(ij,be) - Z(e,b) Tau(ij,ae) ) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP8, 0, 12, 15, 12, 15, 0, "ZZ1(ij,ab)");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tauijab");
      dpd_->file2_init(&ZZ, PSIF_CC_TMP3, 0, 3, 3, "Z(e,a)");
      dpd_->contract424(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&ZZ);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z1, PSIF_CC_TMP9, pqsr, 12, 15, "Z2(ij,ba)");
      dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 12, 15, 12, 15, 0, "Z2(ij,ba)");
      dpd_->buf4_axpy(&Z2, &Z1, -1.0);
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 15, 12, 17, 0, "Gijab");
      dpd_->buf4_axpy(&Z1, &G, 1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_close(&G);
      /* - P(ij) P(ab) ( T'(ia,me) (TMP5) V(jb,me) + T(ia,ME) (T2) V(jb,ME) ) */
      dpd_->buf4_init(&Z, PSIF_CC_TMP8, 0, 30, 30, 30, 30, 0, "Z(ia,jb)");
      dpd_->buf4_init(&T, PSIF_CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 30, 30, 30, 30, 0, "Viajb");
      dpd_->contract444(&T, &V, &Z, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 30, 20, 30, 20, 0, "ViaJB");
      dpd_->contract444(&T, &V, &Z, 0, 0, -1.0, 1.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z, PSIF_CC_TMP9, rqps, 30, 30, "Z(ja,ib)");
      dpd_->buf4_sort(&Z, PSIF_CC_TMP10, psrq, 30, 30, "Z(ib,ja)");
      dpd_->buf4_sort(&Z, PSIF_CC_TMP11, rspq, 30, 30, "Z(jb,ia)");
      dpd_->buf4_init(&Z1, PSIF_CC_TMP9, 0, 30, 30, 30, 30, 0, "Z(ja,ib)");
      dpd_->buf4_axpy(&Z1, &Z, -1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_init(&Z1, PSIF_CC_TMP10, 0, 30, 30, 30, 30, 0, "Z(ib,ja)");
      dpd_->buf4_axpy(&Z1, &Z, -1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_init(&Z1, PSIF_CC_TMP11, 0, 30, 30, 30, 30, 0, "Z(jb,ia)");
      dpd_->buf4_axpy(&Z1, &Z, 1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_sort(&Z, PSIF_CC_TMP9, prqs, 10, 15, "Z(ij,ab)");
      dpd_->buf4_close(&Z);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 15, 12, 17, 0, "Gijab");
      dpd_->buf4_init(&Z, PSIF_CC_TMP9, 0, 10, 15, 10, 15, 0, "Z(ij,ab)");
      /* I don't understand this factor of 1/2 that shows up here */
      dpd_->buf4_axpy(&Z, &G, -0.5);
      dpd_->buf4_close(&Z);
      dpd_->buf4_close(&G);
      /* T'(ia,me) (TMP5) L(m,e) + T'(ia,ME) (T2) L(M,E) --> ZZ(i,a) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 2, 3, "ZZ(i,a)");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      dpd_->buf4_init(&T, PSIF_CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_->contract422(&T, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&T);
      dpd_->file2_close(&L1);
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_->contract422(&T, &L1, &ZZ, 0, 0, -1.0, 1.0);
      dpd_->buf4_close(&T);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&ZZ);
      /* - P(ij) P(ab) ZZ(i,a) T(j,b) --> G(ij,ab) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 2, 3, "ZZ(i,a)");
      dpd_->file2_mat_init(&ZZ);
      dpd_->file2_mat_rd(&ZZ);
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->file2_mat_init(&T1);
      dpd_->file2_mat_rd(&T1);

      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&G, h);
	dpd_->buf4_mat_irrep_rd(&G, h);
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

	    I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	    J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];

	    if((Jsym==Asym) && (Isym==Bsym))
	      value -= ZZ.matrix[Jsym][J][A] * T1.matrix[Isym][I][B];

	    I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	    J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
	    A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
	    B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

	    if((Isym==Bsym) && (Jsym==Asym))
	      value -= ZZ.matrix[Isym][I][B] * T1.matrix[Jsym][J][A];

	    I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	    J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];

	    if((Isym==Asym) && (Jsym==Bsym))
	      value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

	    G.matrix[h][row][col] -= value;
	      
	  }
	}
	dpd_->buf4_mat_irrep_wrt(&G, h);
	dpd_->buf4_mat_irrep_close(&G, h);
      }
      dpd_->buf4_close(&G);

      dpd_->file2_mat_close(&T1);
      dpd_->file2_close(&T1);
      dpd_->file2_mat_close(&ZZ);
      dpd_->file2_close(&ZZ);

      /* T(j,e) L(m,e) --> ZZ(j,m) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 2, 2, "Z(j,m)");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&T1);
      /* ZZ(j,m) T(m,b) --> ZZ2(j,b) */
      dpd_->file2_init(&ZZ2, PSIF_CC_TMP9, 0, 2, 3, "ZZ2(j,b)");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
      dpd_->file2_close(&T1);
      dpd_->file2_close(&ZZ);
      dpd_->file2_close(&ZZ2);

      /* 3 P(ij) P(ab) T(i,a) ZZ(j,b) --> G(ij,ab) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP9, 0, 2, 3, "ZZ2(j,b)");
      dpd_->file2_mat_init(&ZZ);
      dpd_->file2_mat_rd(&ZZ);
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->file2_mat_init(&T1);
      dpd_->file2_mat_rd(&T1);

      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&G, h);
	dpd_->buf4_mat_irrep_rd(&G, h);
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

	    I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	    J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];

	    if((Jsym==Asym) && (Isym==Bsym))
	      value -= T1.matrix[Jsym][J][A] * ZZ.matrix[Isym][I][B];

	    I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	    J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
	    A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
	    B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

	    if((Isym==Bsym) && (Jsym==Asym))
	      value -= T1.matrix[Isym][I][B] * ZZ.matrix[Jsym][J][A];

	    I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	    J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];

	    if((Isym==Asym) && (Jsym==Bsym))
	      value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

	    G.matrix[h][row][col] += 3.0 * value;
	      
	  }
	}
	dpd_->buf4_mat_irrep_wrt(&G, h);
	dpd_->buf4_mat_irrep_close(&G, h);
      }
      dpd_->buf4_scm(&G, 0.5);
      dpd_->buf4_close(&G);

      dpd_->file2_mat_close(&T1);
      dpd_->file2_close(&T1);
      dpd_->file2_mat_close(&ZZ);
      dpd_->file2_close(&ZZ);



      /* L(Ij,Ab) */
      if(params.wfn == "CCSD_T" && params.dertype==1) {
	/* For CCSD(T) gradients, some density contributions are
	   calculated in cctriples */
	dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
	dpd_->buf4_axpy(&L, &G, 1.0);
	dpd_->buf4_close(&L);
	dpd_->buf4_close(&G);
      }
      else {
	dpd_->buf4_init(&L, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
	dpd_->buf4_copy(&L, PSIF_CC_GAMMA, "GIjAb");
	dpd_->buf4_close(&L);
      }
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      /* Tau(Ij,Ab) * (L0*R0 = 1, ground or 0, excited */
      if (params.ground) {
	dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
	dpd_->buf4_axpy(&T, &G, 1.0);
	dpd_->buf4_close(&T);
      }
      /* V(Ij,Mn) Tau(Mn,Ab) */
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 22, 22, 22, 22, 0, "VMnIj");
      dpd_->contract444(&V, &T, &G, 0, 1, 1.0, 1.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_close(&G);
      /* - ( Z(I,M) Tau(Mj,Ab) - Z(j,m) Tau(mI,bA) ) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP8, 0, 22, 28, 22, 28, 0, "Z1(Ij,Ab)");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
      dpd_->file2_init(&ZZ, PSIF_CC_TMP0, 0, 0, 0, "Z(I,M)");
      dpd_->contract244(&ZZ, &T, &Z1, 1, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&ZZ);
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 23, 29, 23, 29, 0, "Z2(jI,bA)");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
      dpd_->file2_init(&ZZ, PSIF_CC_TMP1, 0, 2, 2, "Z(i,m)");
      dpd_->contract244(&ZZ, &T, &Z2, 1, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&ZZ);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z2, PSIF_CC_TMP10, qprs, 22, 29, "Z2(Ij,bA)");
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&Z2, PSIF_CC_TMP10, 0, 22, 29, 22, 29, 0, "Z2(Ij,bA)");
      dpd_->buf4_sort(&Z2, PSIF_CC_TMP9, pqsr, 22, 28, "Z2(Ij,Ab)");
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 22, 28, 22, 28, 0, "Z2(Ij,Ab)");
      dpd_->buf4_axpy(&Z2, &Z1, 1.0);
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      dpd_->buf4_axpy(&Z1, &G, -1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_close(&G);
      /* - ( Z(E,A) Tau(Ij,bE) - Z(e,b) Tau(Ij,Ae) ) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP8, 0, 22, 28, 22, 28, 0, "ZZ1(Ij,Ab)");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
      dpd_->file2_init(&ZZ, PSIF_CC_TMP3, 0, 3, 3, "Z(e,a)");
      dpd_->contract424(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&ZZ);
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 23, 29, 23, 29, 0, "Z2(jI,bA)");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
      dpd_->file2_init(&ZZ, PSIF_CC_TMP2, 0, 1, 1, "Z(E,A)");
      dpd_->contract424(&T, &ZZ, &Z2, 3, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&ZZ);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z2, PSIF_CC_TMP10, qprs, 22, 29, "Z2(Ij,bA)");
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&Z2, PSIF_CC_TMP10, 0, 22, 29, 22, 29, 0, "Z2(Ij,bA)");
      dpd_->buf4_sort(&Z2, PSIF_CC_TMP9, pqsr, 22, 28, "Z2(Ij,Ab)");
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&Z2, PSIF_CC_TMP9, 0, 22, 28, 22, 28, 0, "Z2(Ij,Ab)");
      dpd_->buf4_axpy(&Z2, &Z1, 1.0);
      dpd_->buf4_close(&Z2);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      dpd_->buf4_axpy(&Z1, &G, 1.0);
      dpd_->buf4_close(&Z1);
      dpd_->buf4_close(&G);
      /* - P(Ij) P(Ab) ( T'(IA,me) (T2) V(jb,me) + T'(IA,ME) (TMP4) V(jb,ME) ) */
      dpd_->buf4_init(&Z, PSIF_CC_TMP8, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
      dpd_->buf4_init(&T, PSIF_CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 30, 20, 30, 20, 0, "ViaJB");
      dpd_->contract444(&T, &V, &Z, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 30, 30, 30, 30, 0, "Viajb");
      dpd_->contract444(&T, &V, &Z, 0, 0, -1.0, 1.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      /* T'(jA,Me) V(Ib,Me) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP9, 0, 27, 24, 27, 24, 0, "Z(jA,Ib)");
      dpd_->buf4_init(&T, PSIF_CC_TMP6, 0, 27, 24, 27, 24, 0, "Z(iA,Me)");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 24, 24, 24, 24, 0, "VIaJb");
      dpd_->contract444(&T, &V, &Z1, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z1, PSIF_CC_TMP10, rqps, 20, 30, "Z(IA,jb)");
      dpd_->buf4_close(&Z1);
      dpd_->buf4_init(&Z1, PSIF_CC_TMP10, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
      dpd_->buf4_axpy(&Z1, &Z, -1.0);
      dpd_->buf4_close(&Z1);
      /* T'(Ib,mE) V(jA,mE) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP9, 0, 24, 27, 24, 27, 0, "Z(Ib,jA)");
      dpd_->buf4_init(&T, PSIF_CC_TMP7, 0, 24, 27, 24, 27, 0, "Z(Ia,mE)");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 27, 27, 27, 27, 0, "ViAjB");
      dpd_->contract444(&T, &V, &Z1, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_sort(&Z1, PSIF_CC_TMP10, psrq, 20, 30, "Z(IA,jb)");
      dpd_->buf4_close(&Z1);
      dpd_->buf4_init(&Z1, PSIF_CC_TMP10, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
      dpd_->buf4_axpy(&Z1, &Z, -1.0);
      dpd_->buf4_close(&Z1);
      /* T'(jb,ME) (T2) V(IA,ME) + T'(jb,me) (TMP5) V(IA,me) */
      dpd_->buf4_init(&Z1, PSIF_CC_TMP9, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
      dpd_->buf4_init(&T, PSIF_CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 20, 30, 20, 30, 0, "VIAjb");
      dpd_->contract444(&V, &T, &Z1, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_->buf4_init(&V, PSIF_CC_MISC, 0, 20, 20, 20, 20, 0, "VIAJB");
      dpd_->contract444(&V, &T, &Z1, 0, 0, -1.0, 1.0);
      dpd_->buf4_close(&V);
      dpd_->buf4_close(&T);
      dpd_->buf4_axpy(&Z1, &Z, 1.0);
      dpd_->buf4_close(&Z1);
      /* - Z(IA,jb) --> G(Ij,Ab) */
      dpd_->buf4_sort(&Z, PSIF_CC_TMP9, prqs, 22, 28, "Z(Ij,Ab)");
      dpd_->buf4_close(&Z);
      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      dpd_->buf4_init(&Z, PSIF_CC_TMP9, 0, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
      dpd_->buf4_axpy(&Z, &G, -0.5);
      dpd_->buf4_close(&Z);
      dpd_->buf4_close(&G);
      /* T'(IA,me) (T2) L(m,e) + T'(IA,ME) (TMP4) L(M,E) --> ZZ(I,A) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 1, "ZZ(I,A)");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      dpd_->buf4_init(&T, PSIF_CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_->contract422(&T, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_->buf4_close(&T);
      dpd_->file2_close(&L1);
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_->contract422(&T, &L1, &ZZ, 0, 0, -1.0, 1.0);
      dpd_->buf4_close(&T);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&ZZ);
      /* - ZZ(I,A) T(j,b) --> G(Ij,Ab) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 1, "ZZ(I,A)");
      dpd_->file2_mat_init(&ZZ);
      dpd_->file2_mat_rd(&ZZ);
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->file2_mat_init(&T1);
      dpd_->file2_mat_rd(&T1);

      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&G, h);
	dpd_->buf4_mat_irrep_rd(&G, h);
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
	dpd_->buf4_mat_irrep_wrt(&G, h);
	dpd_->buf4_mat_irrep_close(&G, h);
      }
      dpd_->buf4_close(&G);

      dpd_->file2_mat_close(&T1);
      dpd_->file2_close(&T1);
      dpd_->file2_mat_close(&ZZ);
      dpd_->file2_close(&ZZ);

      /* T'(jb,ME) (T2) L(M,E) + T'(jb,me) (TMP5) L(m,e) --> ZZ(j,b) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 2, 3, "ZZ(j,b)");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_->contract422(&T, &L1, &ZZ, 0, 0, -1.0, 0.0);
      dpd_->buf4_close(&T);
      dpd_->file2_close(&L1);
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      dpd_->buf4_init(&T, PSIF_CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_->contract422(&T, &L1, &ZZ, 0, 0, 1.0, 1.0);
      dpd_->buf4_close(&T);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&ZZ);
      /* - ZZ(j,b) T(I,A) --> G(Ij,Ab) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 2, 3, "ZZ(j,b)");
      dpd_->file2_mat_init(&ZZ);
      dpd_->file2_mat_rd(&ZZ);
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_mat_init(&T1);
      dpd_->file2_mat_rd(&T1);

      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&G, h);
	dpd_->buf4_mat_irrep_rd(&G, h);
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
	dpd_->buf4_mat_irrep_wrt(&G, h);
	dpd_->buf4_mat_irrep_close(&G, h);
      }
      dpd_->buf4_close(&G);

      dpd_->file2_mat_close(&T1);
      dpd_->file2_close(&T1);
      dpd_->file2_mat_close(&ZZ);
      dpd_->file2_close(&ZZ);

      /* T(j,e) L(m,e) --> ZZ(j,m) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 2, 2, "Z(j,m)");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&T1);
      /* ZZ(j,m) T(m,b) --> ZZ2(j,b) */
      dpd_->file2_init(&ZZ2, PSIF_CC_TMP9, 0, 2, 3, "ZZ2(j,b)");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
      dpd_->file2_close(&T1);
      dpd_->file2_close(&ZZ);
      dpd_->file2_close(&ZZ2);

      /* 3 T(I,A) ZZ(j,b) --> G(Ij,Ab) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP9, 0, 2, 3, "ZZ2(j,b)");
      dpd_->file2_mat_init(&ZZ);
      dpd_->file2_mat_rd(&ZZ);
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_mat_init(&T1);
      dpd_->file2_mat_rd(&T1);

      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&G, h);
	dpd_->buf4_mat_irrep_rd(&G, h);
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
	dpd_->buf4_mat_irrep_wrt(&G, h);
	dpd_->buf4_mat_irrep_close(&G, h);
      }
      dpd_->buf4_close(&G);

      dpd_->file2_mat_close(&T1);
      dpd_->file2_close(&T1);
      dpd_->file2_mat_close(&ZZ);
      dpd_->file2_close(&ZZ);

      /* T(I,E) L(M,E) --> ZZ(I,M) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP8, 0, 0, 0, "Z(I,M)");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      dpd_->contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_->file2_close(&L1);
      dpd_->file2_close(&T1);

      /* ZZ(I,M) T(M,A) --> ZZ2(I,A) */
      dpd_->file2_init(&ZZ2, PSIF_CC_TMP9, 0, 0, 1, "ZZ2(I,A)");
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      dpd_->contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
      dpd_->file2_close(&T1);
      dpd_->file2_close(&ZZ);
      dpd_->file2_close(&ZZ2);

      /* 3 T(j,b) ZZ(I,A) --> G(Ij,Ab) */
      dpd_->file2_init(&ZZ, PSIF_CC_TMP9, 0, 0, 1, "ZZ2(I,A)");
      dpd_->file2_mat_init(&ZZ);
      dpd_->file2_mat_rd(&ZZ);
      dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
      dpd_->file2_mat_init(&T1);
      dpd_->file2_mat_rd(&T1);

      dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	dpd_->buf4_mat_irrep_init(&G, h);
	dpd_->buf4_mat_irrep_rd(&G, h);
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
	dpd_->buf4_mat_irrep_wrt(&G, h);
	dpd_->buf4_mat_irrep_close(&G, h);
      }
      dpd_->buf4_scm(&G, 0.5);
      dpd_->buf4_close(&G);

      dpd_->file2_mat_close(&T1);
      dpd_->file2_close(&T1);
      dpd_->file2_mat_close(&ZZ);
      dpd_->file2_close(&ZZ);
    }

  }} // namespace psi::ccdensity
