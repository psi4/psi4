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
      dpd_file2_init(&g, CC_GLG, 0, 0, 0, "GMI");
      dpd_file2_copy(&g, CC_TMP0, "Z(I,M)");
      dpd_file2_close(&g);
      dpd_file2_init(&ZZ, CC_TMP0, 0, 0, 0, "Z(I,M)");
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 1.0);
      dpd_file2_close(&T1);
      dpd_file2_close(&L1);
      dpd_file2_close(&ZZ);

      /* ( g(i,m) + L(m,e) T(i,e) ) --> Z(i,m)(TMP1)  */
      dpd_file2_init(&g, CC_GLG, 0, 2, 2, "Gmi");
      dpd_file2_copy(&g, CC_TMP1, "Z(i,m)");
      dpd_file2_close(&g);
      dpd_file2_init(&ZZ, CC_TMP1, 0, 2, 2, "Z(i,m)");
      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 1.0);
      dpd_file2_close(&T1);
      dpd_file2_close(&L1);
      dpd_file2_close(&ZZ);

      /* ( g(E,A) - L(M,E) T(M,A) ) --> Z(E,A)(TMP2) */
      dpd_file2_init(&g, CC_GLG, 0, 1, 1, "GAE");
      dpd_file2_copy(&g, CC_TMP2, "Z(E,A)");
      dpd_file2_close(&g);
      dpd_file2_init(&ZZ, CC_TMP2, 0, 1, 1, "Z(E,A)");
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_contract222(&L1, &T1, &ZZ, 1, 1, -1.0, 1.0);
      dpd_file2_close(&T1);
      dpd_file2_close(&L1);
      dpd_file2_close(&ZZ);

      /* ( g(e,a) - L(m,e) T(m,a) ) --> Z(e,a)(TMP3) */
      dpd_file2_init(&g, CC_GLG, 0, 3, 3, "Gae");
      dpd_file2_copy(&g, CC_TMP3, "Z(e,a)");
      dpd_file2_close(&g);
      dpd_file2_init(&ZZ, CC_TMP3, 0, 3, 3, "Z(e,a)");
      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_contract222(&L1, &T1, &ZZ, 1, 1, -1.0, 1.0);
      dpd_file2_close(&T1);
      dpd_file2_close(&L1);
      dpd_file2_close(&ZZ);

      dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_mat_init(&T1A);
      dpd_file2_mat_rd(&T1A);
      dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_mat_init(&T1B);
      dpd_file2_mat_rd(&T1B);

      /* ( - T(IA,ME) + 2 * T(I,E) T(M,A) ) --> Z(IA,ME) */
      dpd_buf4_init(&T, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
      dpd_buf4_copy(&T, CC_TMP4, "Z(IA,ME)");
      dpd_buf4_close(&T);
      dpd_buf4_init(&Z, CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_buf4_scm(&Z, -1.0);
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&Z, h);
	dpd_buf4_mat_irrep_rd(&Z, h);
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
	dpd_buf4_mat_irrep_wrt(&Z, h);
	dpd_buf4_mat_irrep_close(&Z, h);
      }
      dpd_buf4_close(&Z);

      /* ( - T(ia,me) + 2 * T(i,e) T(m,a) ) --> Z(ia,me) */
      dpd_buf4_init(&T, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
      dpd_buf4_copy(&T, CC_TMP5, "Z(ia,me)");
      dpd_buf4_close(&T);
      dpd_buf4_init(&Z, CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_buf4_scm(&Z, -1.0);
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&Z, h);
	dpd_buf4_mat_irrep_rd(&Z, h);
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
	dpd_buf4_mat_irrep_wrt(&Z, h);
	dpd_buf4_mat_irrep_close(&Z, h);
      }
      dpd_buf4_close(&Z);

      /* ( - T(iA,Me) + 2 * T(i,e) T(M,A) ) --> Z(iA,Me) */
      dpd_buf4_init(&T, CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
      dpd_buf4_copy(&T, CC_TMP6, "Z(iA,Me)");
      dpd_buf4_close(&T);
      dpd_buf4_init(&Z, CC_TMP6, 0, 27, 24, 27, 24, 0, "Z(iA,Me)");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&Z, h);
	dpd_buf4_mat_irrep_rd(&Z, h);
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
	dpd_buf4_mat_irrep_wrt(&Z, h);
	dpd_buf4_mat_irrep_close(&Z, h);
      }
      dpd_buf4_close(&Z);

      /* ( - T(Ia,mE) + 2 * T(I,E) T(m,a) ) --> Z(Ia,mE) */
      dpd_buf4_init(&T, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
      dpd_buf4_copy(&T, CC_TMP7, "Z(Ia,mE)");
      dpd_buf4_close(&T);
      dpd_buf4_init(&Z, CC_TMP7, 0, 24, 27, 24, 27, 0, "Z(Ia,mE)");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&Z, h);
	dpd_buf4_mat_irrep_rd(&Z, h);
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
	dpd_buf4_mat_irrep_wrt(&Z, h);
	dpd_buf4_mat_irrep_close(&Z, h);
      }
      dpd_buf4_close(&Z);

      dpd_file2_mat_close(&T1A);
      dpd_file2_close(&T1A);
      dpd_file2_mat_close(&T1B);
      dpd_file2_close(&T1B);

      /* L(IJ,AB) */
      if(params.wfn == "CCSD_T" && params.dertype==1) {
	/* For CCSD(T) gradients, some density contributions are
	   calculated in cctriples */
	dpd_buf4_init(&G, CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
	dpd_buf4_init(&L, CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
	dpd_buf4_axpy(&L, &G, 1.0);
	dpd_buf4_close(&L);
	dpd_buf4_close(&G);
      }
      else {
	dpd_buf4_init(&L, CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
	dpd_buf4_copy(&L, CC_GAMMA, "GIJAB");
	dpd_buf4_close(&L);
      }
      dpd_buf4_init(&G, CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
      /* Tau(IJ,AB) * (L0*R0 = 1, ground or 0, excited */
      if (params.ground) {
	dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
	dpd_buf4_axpy(&T, &G, 1.0);
	dpd_buf4_close(&T);
      }
      /* V(IJ,MN) Tau(MN,AB) */
      dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
      dpd_buf4_init(&V, CC_MISC, 0, 2, 2, 2, 2, 0, "VMNIJ");
      dpd_contract444(&V, &T, &G, 0, 1, 1.0, 1.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_close(&G);
      /* - ( Z(I,M) Tau(MJ,AB) - Z(J,M) Tau(MI,AB) ) */
      dpd_buf4_init(&Z1, CC_TMP8, 0, 0, 7, 0, 7, 0, "Z1(IJ,AB)");
      dpd_buf4_init(&T, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
      dpd_file2_init(&ZZ, CC_TMP0, 0, 0, 0, "Z(I,M)");
      dpd_contract244(&ZZ, &T, &Z1, 1, 0, 0, -1.0, 0.0);
      dpd_file2_close(&ZZ);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z1, CC_TMP9, qprs, 0, 7, "Z2(JI,AB)");
      dpd_buf4_init(&Z2, CC_TMP9, 0, 0, 7, 0, 7, 0, "Z2(JI,AB)");
      dpd_buf4_axpy(&Z2, &Z1, -1.0);
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
      dpd_buf4_axpy(&Z1, &G, 1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_close(&G);
      /* - ( Z(E,A) Tau(IJ,BE) - Z(E,B) Tau(IJ,AE) ) */
      dpd_buf4_init(&Z1, CC_TMP8, 0, 2, 5, 2, 5, 0, "ZZ1(IJ,AB)");
      dpd_buf4_init(&T, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");
      dpd_file2_init(&ZZ, CC_TMP2, 0, 1, 1, "Z(E,A)");
      dpd_contract424(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0);
      dpd_file2_close(&ZZ);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z1, CC_TMP9, pqsr, 2, 5, "Z2(IJ,BA)");
      dpd_buf4_init(&Z2, CC_TMP9, 0, 2, 5, 2, 5, 0, "Z2(IJ,BA)");
      dpd_buf4_axpy(&Z2, &Z1, -1.0);
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
      dpd_buf4_axpy(&Z1, &G, 1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_close(&G);
      /* - P(IJ) P(AB) ( T'(IA,ME) (TMP4) V(JB,ME) + T(IA,me) (T2) V(JB,me) ) */
      dpd_buf4_init(&Z, CC_TMP8, 0, 20, 20, 20, 20, 0, "Z(IA,JB)");
      dpd_buf4_init(&T, CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_buf4_init(&V, CC_MISC, 0, 20, 20, 20, 20, 0, "VIAJB");
      dpd_contract444(&T, &V, &Z, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_init(&T, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_buf4_init(&V, CC_MISC, 0, 20, 30, 20, 30, 0, "VIAjb");
      dpd_contract444(&T, &V, &Z, 0, 0, -1.0, 1.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z, CC_TMP9, rqps, 20, 20, "Z(JA,IB)");
      dpd_buf4_sort(&Z, CC_TMP10, psrq, 20, 20, "Z(IB,JA)");
      dpd_buf4_sort(&Z, CC_TMP11, rspq, 20, 20, "Z(JB,IA)");
      dpd_buf4_init(&Z1, CC_TMP9, 0, 20, 20, 20, 20, 0, "Z(JA,IB)");
      dpd_buf4_axpy(&Z1, &Z, -1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_init(&Z1, CC_TMP10, 0, 20, 20, 20, 20, 0, "Z(IB,JA)");
      dpd_buf4_axpy(&Z1, &Z, -1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_init(&Z1, CC_TMP11, 0, 20, 20, 20, 20, 0, "Z(JB,IA)");
      dpd_buf4_axpy(&Z1, &Z, 1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_sort(&Z, CC_TMP9, prqs, 0, 5, "Z(IJ,AB)");
      dpd_buf4_close(&Z);
      dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 2, 7, 0, "GIJAB");
      dpd_buf4_init(&Z, CC_TMP9, 0, 0, 5, 0, 5, 0, "Z(IJ,AB)");
      /* I don't understand this factor of 1/2 that shows up here */
      dpd_buf4_axpy(&Z, &G, -0.5);
      dpd_buf4_close(&Z);
      dpd_buf4_close(&G);
      /* T'(IA,ME) (TMP4) L(M,E) + T'(IA,me) (T2) L(m,e) --> ZZ(I,A) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 0, 1, "ZZ(I,A)");
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_buf4_init(&T, CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_contract422(&T, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&T);
      dpd_file2_close(&L1);
      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_buf4_init(&T, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_contract422(&T, &L1, &ZZ, 0, 0, -1.0, 1.0);
      dpd_buf4_close(&T);
      dpd_file2_close(&L1);
      dpd_file2_close(&ZZ);
      /* - P(IJ) P(AB) ZZ(I,A) T(J,B) --> G(IJ,AB) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 0, 1, "ZZ(I,A)");
      dpd_file2_mat_init(&ZZ);
      dpd_file2_mat_rd(&ZZ);
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_buf4_init(&G, CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
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
	dpd_buf4_mat_irrep_wrt(&G, h);
	dpd_buf4_mat_irrep_close(&G, h);
      }
      dpd_buf4_close(&G);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);
      dpd_file2_mat_close(&ZZ);
      dpd_file2_close(&ZZ);

      /* T(J,E) L(M,E) --> ZZ(J,M) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 0, 0, "Z(J,M)");
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_file2_close(&L1);
      dpd_file2_close(&T1);
      /* ZZ(J,M) T(M,B) --> ZZ2(J,B) */
      dpd_file2_init(&ZZ2, CC_TMP9, 0, 0, 1, "ZZ2(J,B)");
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
      dpd_file2_close(&T1);
      dpd_file2_close(&ZZ);
      dpd_file2_close(&ZZ2);

      /* 3 P(IJ) P(AB) T(I,A) ZZ(J,B) --> G(IJ,AB) */
      dpd_file2_init(&ZZ, CC_TMP9, 0, 0, 1, "ZZ2(J,B)");
      dpd_file2_mat_init(&ZZ);
      dpd_file2_mat_rd(&ZZ);
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_buf4_init(&G, CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
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
	dpd_buf4_mat_irrep_wrt(&G, h);
	dpd_buf4_mat_irrep_close(&G, h);
      }
      dpd_buf4_scm(&G, 0.5);
      dpd_buf4_close(&G);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);
      dpd_file2_mat_close(&ZZ);
      dpd_file2_close(&ZZ);



      /* L(ij,ab) */
      if(params.wfn == "CCSD_T" && params.dertype==1) {
	/* For CCSD(T) gradients, some density contributions are
	   calculated in cctriples */
	dpd_buf4_init(&G, CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
	dpd_buf4_init(&L, CC_GLG, 0, 12, 17, 12, 17, 0, "Lijab");
	dpd_buf4_axpy(&L, &G, 1.0);
	dpd_buf4_close(&L);
	dpd_buf4_close(&G);
      }
      else {
	dpd_buf4_init(&L, CC_GLG, 0, 12, 17, 12, 17, 0, "Lijab");
	dpd_buf4_copy(&L, CC_GAMMA, "Gijab");
	dpd_buf4_close(&L);
      }
      dpd_buf4_init(&G, CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
      /* Tau(ij,ab) * (L0*R0 = 1, ground or 0, excited */
      if (params.ground) {
	dpd_buf4_init(&T, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
	dpd_buf4_axpy(&T, &G, 1.0);
	dpd_buf4_close(&T);
      }
      /* V(ij,mn) Tau(mn,ab) */
      dpd_buf4_init(&T, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
      dpd_buf4_init(&V, CC_MISC, 0, 12, 12, 12, 12, 0, "Vmnij");
      dpd_contract444(&V, &T, &G, 0, 1, 1.0, 1.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_close(&G);
      /* - ( Z(i,m) Tau(mj,ab) - Z(j,m) Tau(mi,ab) ) */
      dpd_buf4_init(&Z1, CC_TMP8, 0, 10, 17, 10, 17, 0, "Z1(ij,ab)");
      dpd_buf4_init(&T, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tauijab");
      dpd_file2_init(&ZZ, CC_TMP1, 0, 2, 2, "Z(i,m)");
      dpd_contract244(&ZZ, &T, &Z1, 1, 0, 0, -1.0, 0.0);
      dpd_file2_close(&ZZ);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z1, CC_TMP9, qprs, 10, 17, "Z2(ji,ab)");
      dpd_buf4_init(&Z2, CC_TMP9, 0, 10, 17, 10, 17, 0, "Z2(ji,ab)");
      dpd_buf4_axpy(&Z2, &Z1, -1.0);
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 17, 12, 17, 0, "Gijab");
      dpd_buf4_axpy(&Z1, &G, 1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_close(&G);
      /* - ( Z(e,a) Tau(ij,be) - Z(e,b) Tau(ij,ae) ) */
      dpd_buf4_init(&Z1, CC_TMP8, 0, 12, 15, 12, 15, 0, "ZZ1(ij,ab)");
      dpd_buf4_init(&T, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tauijab");
      dpd_file2_init(&ZZ, CC_TMP3, 0, 3, 3, "Z(e,a)");
      dpd_contract424(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0);
      dpd_file2_close(&ZZ);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z1, CC_TMP9, pqsr, 12, 15, "Z2(ij,ba)");
      dpd_buf4_init(&Z2, CC_TMP9, 0, 12, 15, 12, 15, 0, "Z2(ij,ba)");
      dpd_buf4_axpy(&Z2, &Z1, -1.0);
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&G, CC_GAMMA, 0, 12, 15, 12, 17, 0, "Gijab");
      dpd_buf4_axpy(&Z1, &G, 1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_close(&G);
      /* - P(ij) P(ab) ( T'(ia,me) (TMP5) V(jb,me) + T(ia,ME) (T2) V(jb,ME) ) */
      dpd_buf4_init(&Z, CC_TMP8, 0, 30, 30, 30, 30, 0, "Z(ia,jb)");
      dpd_buf4_init(&T, CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_buf4_init(&V, CC_MISC, 0, 30, 30, 30, 30, 0, "Viajb");
      dpd_contract444(&T, &V, &Z, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_init(&T, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_buf4_init(&V, CC_MISC, 0, 30, 20, 30, 20, 0, "ViaJB");
      dpd_contract444(&T, &V, &Z, 0, 0, -1.0, 1.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z, CC_TMP9, rqps, 30, 30, "Z(ja,ib)");
      dpd_buf4_sort(&Z, CC_TMP10, psrq, 30, 30, "Z(ib,ja)");
      dpd_buf4_sort(&Z, CC_TMP11, rspq, 30, 30, "Z(jb,ia)");
      dpd_buf4_init(&Z1, CC_TMP9, 0, 30, 30, 30, 30, 0, "Z(ja,ib)");
      dpd_buf4_axpy(&Z1, &Z, -1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_init(&Z1, CC_TMP10, 0, 30, 30, 30, 30, 0, "Z(ib,ja)");
      dpd_buf4_axpy(&Z1, &Z, -1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_init(&Z1, CC_TMP11, 0, 30, 30, 30, 30, 0, "Z(jb,ia)");
      dpd_buf4_axpy(&Z1, &Z, 1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_sort(&Z, CC_TMP9, prqs, 10, 15, "Z(ij,ab)");
      dpd_buf4_close(&Z);
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 15, 12, 17, 0, "Gijab");
      dpd_buf4_init(&Z, CC_TMP9, 0, 10, 15, 10, 15, 0, "Z(ij,ab)");
      /* I don't understand this factor of 1/2 that shows up here */
      dpd_buf4_axpy(&Z, &G, -0.5);
      dpd_buf4_close(&Z);
      dpd_buf4_close(&G);
      /* T'(ia,me) (TMP5) L(m,e) + T'(ia,ME) (T2) L(M,E) --> ZZ(i,a) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 2, 3, "ZZ(i,a)");
      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_buf4_init(&T, CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_contract422(&T, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&T);
      dpd_file2_close(&L1);
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_buf4_init(&T, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_contract422(&T, &L1, &ZZ, 0, 0, -1.0, 1.0);
      dpd_buf4_close(&T);
      dpd_file2_close(&L1);
      dpd_file2_close(&ZZ);
      /* - P(ij) P(ab) ZZ(i,a) T(j,b) --> G(ij,ab) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 2, 3, "ZZ(i,a)");
      dpd_file2_mat_init(&ZZ);
      dpd_file2_mat_rd(&ZZ);
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_buf4_init(&G, CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
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
	dpd_buf4_mat_irrep_wrt(&G, h);
	dpd_buf4_mat_irrep_close(&G, h);
      }
      dpd_buf4_close(&G);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);
      dpd_file2_mat_close(&ZZ);
      dpd_file2_close(&ZZ);

      /* T(j,e) L(m,e) --> ZZ(j,m) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 2, 2, "Z(j,m)");
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_file2_close(&L1);
      dpd_file2_close(&T1);
      /* ZZ(j,m) T(m,b) --> ZZ2(j,b) */
      dpd_file2_init(&ZZ2, CC_TMP9, 0, 2, 3, "ZZ2(j,b)");
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
      dpd_file2_close(&T1);
      dpd_file2_close(&ZZ);
      dpd_file2_close(&ZZ2);

      /* 3 P(ij) P(ab) T(i,a) ZZ(j,b) --> G(ij,ab) */
      dpd_file2_init(&ZZ, CC_TMP9, 0, 2, 3, "ZZ2(j,b)");
      dpd_file2_mat_init(&ZZ);
      dpd_file2_mat_rd(&ZZ);
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_buf4_init(&G, CC_GAMMA, 0, 12, 17, 12, 17, 0, "Gijab");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
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
	dpd_buf4_mat_irrep_wrt(&G, h);
	dpd_buf4_mat_irrep_close(&G, h);
      }
      dpd_buf4_scm(&G, 0.5);
      dpd_buf4_close(&G);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);
      dpd_file2_mat_close(&ZZ);
      dpd_file2_close(&ZZ);



      /* L(Ij,Ab) */
      if(params.wfn == "CCSD_T" && params.dertype==1) {
	/* For CCSD(T) gradients, some density contributions are
	   calculated in cctriples */
	dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
	dpd_buf4_init(&L, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
	dpd_buf4_axpy(&L, &G, 1.0);
	dpd_buf4_close(&L);
	dpd_buf4_close(&G);
      }
      else {
	dpd_buf4_init(&L, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
	dpd_buf4_copy(&L, CC_GAMMA, "GIjAb");
	dpd_buf4_close(&L);
      }
      dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      /* Tau(Ij,Ab) * (L0*R0 = 1, ground or 0, excited */
      if (params.ground) {
	dpd_buf4_init(&T, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
	dpd_buf4_axpy(&T, &G, 1.0);
	dpd_buf4_close(&T);
      }
      /* V(Ij,Mn) Tau(Mn,Ab) */
      dpd_buf4_init(&T, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
      dpd_buf4_init(&V, CC_MISC, 0, 22, 22, 22, 22, 0, "VMnIj");
      dpd_contract444(&V, &T, &G, 0, 1, 1.0, 1.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_close(&G);
      /* - ( Z(I,M) Tau(Mj,Ab) - Z(j,m) Tau(mI,bA) ) */
      dpd_buf4_init(&Z1, CC_TMP8, 0, 22, 28, 22, 28, 0, "Z1(Ij,Ab)");
      dpd_buf4_init(&T, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
      dpd_file2_init(&ZZ, CC_TMP0, 0, 0, 0, "Z(I,M)");
      dpd_contract244(&ZZ, &T, &Z1, 1, 0, 0, 1.0, 0.0);
      dpd_file2_close(&ZZ);
      dpd_buf4_close(&T);
      dpd_buf4_init(&Z2, CC_TMP9, 0, 23, 29, 23, 29, 0, "Z2(jI,bA)");
      dpd_buf4_init(&T, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
      dpd_file2_init(&ZZ, CC_TMP1, 0, 2, 2, "Z(i,m)");
      dpd_contract244(&ZZ, &T, &Z2, 1, 0, 0, 1.0, 0.0);
      dpd_file2_close(&ZZ);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z2, CC_TMP10, qprs, 22, 29, "Z2(Ij,bA)");
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&Z2, CC_TMP10, 0, 22, 29, 22, 29, 0, "Z2(Ij,bA)");
      dpd_buf4_sort(&Z2, CC_TMP9, pqsr, 22, 28, "Z2(Ij,Ab)");
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&Z2, CC_TMP9, 0, 22, 28, 22, 28, 0, "Z2(Ij,Ab)");
      dpd_buf4_axpy(&Z2, &Z1, 1.0);
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      dpd_buf4_axpy(&Z1, &G, -1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_close(&G);
      /* - ( Z(E,A) Tau(Ij,bE) - Z(e,b) Tau(Ij,Ae) ) */
      dpd_buf4_init(&Z1, CC_TMP8, 0, 22, 28, 22, 28, 0, "ZZ1(Ij,Ab)");
      dpd_buf4_init(&T, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
      dpd_file2_init(&ZZ, CC_TMP3, 0, 3, 3, "Z(e,a)");
      dpd_contract424(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0);
      dpd_file2_close(&ZZ);
      dpd_buf4_close(&T);
      dpd_buf4_init(&Z2, CC_TMP9, 0, 23, 29, 23, 29, 0, "Z2(jI,bA)");
      dpd_buf4_init(&T, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
      dpd_file2_init(&ZZ, CC_TMP2, 0, 1, 1, "Z(E,A)");
      dpd_contract424(&T, &ZZ, &Z2, 3, 0, 0, 1.0, 0.0);
      dpd_file2_close(&ZZ);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z2, CC_TMP10, qprs, 22, 29, "Z2(Ij,bA)");
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&Z2, CC_TMP10, 0, 22, 29, 22, 29, 0, "Z2(Ij,bA)");
      dpd_buf4_sort(&Z2, CC_TMP9, pqsr, 22, 28, "Z2(Ij,Ab)");
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&Z2, CC_TMP9, 0, 22, 28, 22, 28, 0, "Z2(Ij,Ab)");
      dpd_buf4_axpy(&Z2, &Z1, 1.0);
      dpd_buf4_close(&Z2);
      dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      dpd_buf4_axpy(&Z1, &G, 1.0);
      dpd_buf4_close(&Z1);
      dpd_buf4_close(&G);
      /* - P(Ij) P(Ab) ( T'(IA,me) (T2) V(jb,me) + T'(IA,ME) (TMP4) V(jb,ME) ) */
      dpd_buf4_init(&Z, CC_TMP8, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
      dpd_buf4_init(&T, CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_buf4_init(&V, CC_MISC, 0, 30, 20, 30, 20, 0, "ViaJB");
      dpd_contract444(&T, &V, &Z, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_init(&T, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_buf4_init(&V, CC_MISC, 0, 30, 30, 30, 30, 0, "Viajb");
      dpd_contract444(&T, &V, &Z, 0, 0, -1.0, 1.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      /* T'(jA,Me) V(Ib,Me) */
      dpd_buf4_init(&Z1, CC_TMP9, 0, 27, 24, 27, 24, 0, "Z(jA,Ib)");
      dpd_buf4_init(&T, CC_TMP6, 0, 27, 24, 27, 24, 0, "Z(iA,Me)");
      dpd_buf4_init(&V, CC_MISC, 0, 24, 24, 24, 24, 0, "VIaJb");
      dpd_contract444(&T, &V, &Z1, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z1, CC_TMP10, rqps, 20, 30, "Z(IA,jb)");
      dpd_buf4_close(&Z1);
      dpd_buf4_init(&Z1, CC_TMP10, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
      dpd_buf4_axpy(&Z1, &Z, -1.0);
      dpd_buf4_close(&Z1);
      /* T'(Ib,mE) V(jA,mE) */
      dpd_buf4_init(&Z1, CC_TMP9, 0, 24, 27, 24, 27, 0, "Z(Ib,jA)");
      dpd_buf4_init(&T, CC_TMP7, 0, 24, 27, 24, 27, 0, "Z(Ia,mE)");
      dpd_buf4_init(&V, CC_MISC, 0, 27, 27, 27, 27, 0, "ViAjB");
      dpd_contract444(&T, &V, &Z1, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_sort(&Z1, CC_TMP10, psrq, 20, 30, "Z(IA,jb)");
      dpd_buf4_close(&Z1);
      dpd_buf4_init(&Z1, CC_TMP10, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
      dpd_buf4_axpy(&Z1, &Z, -1.0);
      dpd_buf4_close(&Z1);
      /* T'(jb,ME) (T2) V(IA,ME) + T'(jb,me) (TMP5) V(IA,me) */
      dpd_buf4_init(&Z1, CC_TMP9, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
      dpd_buf4_init(&T, CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_buf4_init(&V, CC_MISC, 0, 20, 30, 20, 30, 0, "VIAjb");
      dpd_contract444(&V, &T, &Z1, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_init(&T, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_buf4_init(&V, CC_MISC, 0, 20, 20, 20, 20, 0, "VIAJB");
      dpd_contract444(&V, &T, &Z1, 0, 0, -1.0, 1.0);
      dpd_buf4_close(&V);
      dpd_buf4_close(&T);
      dpd_buf4_axpy(&Z1, &Z, 1.0);
      dpd_buf4_close(&Z1);
      /* - Z(IA,jb) --> G(Ij,Ab) */
      dpd_buf4_sort(&Z, CC_TMP9, prqs, 22, 28, "Z(Ij,Ab)");
      dpd_buf4_close(&Z);
      dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      dpd_buf4_init(&Z, CC_TMP9, 0, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
      dpd_buf4_axpy(&Z, &G, -0.5);
      dpd_buf4_close(&Z);
      dpd_buf4_close(&G);
      /* T'(IA,me) (T2) L(m,e) + T'(IA,ME) (TMP4) L(M,E) --> ZZ(I,A) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 0, 1, "ZZ(I,A)");
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_buf4_init(&T, CC_TMP4, 0, 20, 20, 20, 20, 0, "Z(IA,ME)");
      dpd_contract422(&T, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&T);
      dpd_file2_close(&L1);
      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_buf4_init(&T, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      dpd_contract422(&T, &L1, &ZZ, 0, 0, -1.0, 1.0);
      dpd_buf4_close(&T);
      dpd_file2_close(&L1);
      dpd_file2_close(&ZZ);
      /* - ZZ(I,A) T(j,b) --> G(Ij,Ab) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 0, 1, "ZZ(I,A)");
      dpd_file2_mat_init(&ZZ);
      dpd_file2_mat_rd(&ZZ);
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
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
	dpd_buf4_mat_irrep_wrt(&G, h);
	dpd_buf4_mat_irrep_close(&G, h);
      }
      dpd_buf4_close(&G);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);
      dpd_file2_mat_close(&ZZ);
      dpd_file2_close(&ZZ);

      /* T'(jb,ME) (T2) L(M,E) + T'(jb,me) (TMP5) L(m,e) --> ZZ(j,b) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 2, 3, "ZZ(j,b)");
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_buf4_init(&T, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      dpd_contract422(&T, &L1, &ZZ, 0, 0, -1.0, 0.0);
      dpd_buf4_close(&T);
      dpd_file2_close(&L1);
      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_buf4_init(&T, CC_TMP5, 0, 30, 30, 30, 30, 0, "Z(ia,me)");
      dpd_contract422(&T, &L1, &ZZ, 0, 0, 1.0, 1.0);
      dpd_buf4_close(&T);
      dpd_file2_close(&L1);
      dpd_file2_close(&ZZ);
      /* - ZZ(j,b) T(I,A) --> G(Ij,Ab) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 2, 3, "ZZ(j,b)");
      dpd_file2_mat_init(&ZZ);
      dpd_file2_mat_rd(&ZZ);
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
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
	dpd_buf4_mat_irrep_wrt(&G, h);
	dpd_buf4_mat_irrep_close(&G, h);
      }
      dpd_buf4_close(&G);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);
      dpd_file2_mat_close(&ZZ);
      dpd_file2_close(&ZZ);

      /* T(j,e) L(m,e) --> ZZ(j,m) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 2, 2, "Z(j,m)");
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
      dpd_contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_file2_close(&L1);
      dpd_file2_close(&T1);
      /* ZZ(j,m) T(m,b) --> ZZ2(j,b) */
      dpd_file2_init(&ZZ2, CC_TMP9, 0, 2, 3, "ZZ2(j,b)");
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
      dpd_file2_close(&T1);
      dpd_file2_close(&ZZ);
      dpd_file2_close(&ZZ2);

      /* 3 T(I,A) ZZ(j,b) --> G(Ij,Ab) */
      dpd_file2_init(&ZZ, CC_TMP9, 0, 2, 3, "ZZ2(j,b)");
      dpd_file2_mat_init(&ZZ);
      dpd_file2_mat_rd(&ZZ);
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
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
	dpd_buf4_mat_irrep_wrt(&G, h);
	dpd_buf4_mat_irrep_close(&G, h);
      }
      dpd_buf4_close(&G);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);
      dpd_file2_mat_close(&ZZ);
      dpd_file2_close(&ZZ);

      /* T(I,E) L(M,E) --> ZZ(I,M) */
      dpd_file2_init(&ZZ, CC_TMP8, 0, 0, 0, "Z(I,M)");
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
      dpd_contract222(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0);
      dpd_file2_close(&L1);
      dpd_file2_close(&T1);

      /* ZZ(I,M) T(M,A) --> ZZ2(I,A) */
      dpd_file2_init(&ZZ2, CC_TMP9, 0, 0, 1, "ZZ2(I,A)");
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_contract222(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0);
      dpd_file2_close(&T1);
      dpd_file2_close(&ZZ);
      dpd_file2_close(&ZZ2);

      /* 3 T(j,b) ZZ(I,A) --> G(Ij,Ab) */
      dpd_file2_init(&ZZ, CC_TMP9, 0, 0, 1, "ZZ2(I,A)");
      dpd_file2_mat_init(&ZZ);
      dpd_file2_mat_rd(&ZZ);
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_buf4_init(&G, CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&G, h);
	dpd_buf4_mat_irrep_rd(&G, h);
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
	dpd_buf4_mat_irrep_wrt(&G, h);
	dpd_buf4_mat_irrep_close(&G, h);
      }
      dpd_buf4_scm(&G, 0.5);
      dpd_buf4_close(&G);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);
      dpd_file2_mat_close(&ZZ);
      dpd_file2_close(&ZZ);
    }

  }} // namespace psi::ccdensity
