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

  nirreps = moinfo.nirreps;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    /* G(ia,jb) <-- L(im,ae) T(jm,be) */
    dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "VIAJB");
    dpd_buf4_sort(&V, CC_MISC, rspq, 10, 10, "GIAJB");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "Viajb");
    dpd_buf4_sort(&V, CC_MISC, rspq, 10, 10, "Giajb");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "ViaJB");
    dpd_buf4_sort(&V, CC_MISC, rspq, 10, 10, "GIAjb");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "VIAjb");
    dpd_buf4_sort(&V, CC_MISC, rspq, 10, 10, "GiaJB");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "VIaJb");
    dpd_buf4_sort(&V, CC_MISC, rspq, 10, 10, "GIaJb");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 10, 10, 10, 10, 0, "ViAjB");
    dpd_buf4_sort(&V, CC_MISC, rspq, 10, 10, "GiAjB");
    dpd_buf4_close(&V);

    /* G(IA,JB) <-- - L(IM,AE) T(J,E) T(M,B) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(IM,AJ)");
    dpd_buf4_init(&L, CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(IB,AJ)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, prqs, 10, 11, "Z(IA,BJ)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(IA,BJ)");
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 10, 10, "Z(IA,JB)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(IA,JB)");
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GIAJB");
    dpd_buf4_axpy(&Z1, &G, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&G);

    /* G(ia,jb) <-- - L(im,ae) T(j,e) T(m,b) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(im,aj)");
    dpd_buf4_init(&L, CC_GLG, 0, 0, 5, 2, 7, 0, "Lijab");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(ib,aj)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, prqs, 10, 11, "Z(ia,bj)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(ia,bj)");
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 10, 10, "Z(ia,jb)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(ia,jb)");
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "Giajb");
    dpd_buf4_axpy(&Z1, &G, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&G);

    /* G(IA,jb) <-- - L(Im,Ae) T(j,e) T(m,b) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(Im,Aj)");
    dpd_buf4_init(&L, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, prqs, 10, 11, "Z(IA,bj)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(IA,bj)");
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 10, 10, "Z(IA,jb)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(IA,jb)");
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GIAjb");
    dpd_buf4_axpy(&Z1, &G, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&G);

    /* G(ia,JB) <-- - L(iM,aE) T(J,E) T(M,B) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(iM,aJ)");
    dpd_buf4_init(&L, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(iB,aJ)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, prqs, 10, 11, "Z(ia,BJ)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(ia,BJ)");
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 10, 10, "Z(ia,JB)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(ia,JB)");
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GiaJB");
    dpd_buf4_axpy(&Z1, &G, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&G);

    /* G(Ia,Jb) <-- - L(Im,Ea) T(J,E) T(m,b) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(Im,Ja)");
    dpd_buf4_init(&L, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_file2_close(&T1);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(Ib,Ja)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, psrq, 10, 10, "Z(Ia,Jb)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Ia,Jb)");
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GIaJb");
    dpd_buf4_axpy(&Z1, &G, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    /* G(iA,jB) <-- - L(iM,eA) T(j,e) T(M,B) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(iM,jA)");
    dpd_buf4_init(&L, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract244(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_file2_close(&T1);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(iB,jA)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, psrq, 10, 10, "Z(iA,jB)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(iA,jB)");
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GiAjB");
    dpd_buf4_axpy(&Z1, &G, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1A);
    dpd_file2_mat_rd(&T1A);
    dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_mat_init(&T1B);
    dpd_file2_mat_rd(&T1B);
    dpd_file2_init(&L1A, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_mat_init(&L1A);
    dpd_file2_mat_rd(&L1A);
    dpd_file2_init(&L1B, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_mat_init(&L1B);
    dpd_file2_mat_rd(&L1B);

    /* G(IA,JB) <-- L(I,A) T(J,B) */
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GIAJB");
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&G, h); 0,
					dpd_buf4_mat_irrep_rd(&G, h);

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

      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);


    /* G(ia,jb) <-- L(i,a) T(j,b) */
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "Giajb");
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&G, h); 0,
					dpd_buf4_mat_irrep_rd(&G, h);

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

      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    /* G(IA,jb) <-- L(I,A) T(j,b) */
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GIAjb");
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&G, h); 0,
					dpd_buf4_mat_irrep_rd(&G, h);

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

      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    /* G(ia,JB) <-- L(i,a) T(J,B) */
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GiaJB");
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&G, h); 0,
					dpd_buf4_mat_irrep_rd(&G, h);

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

      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    dpd_file2_mat_close(&L1A);
    dpd_file2_close(&L1A);
    dpd_file2_mat_close(&L1B);
    dpd_file2_close(&L1B);
    dpd_file2_mat_close(&T1A);
    dpd_file2_close(&T1A);
    dpd_file2_mat_close(&T1B);
    dpd_file2_close(&T1B);

    /* Sort all spin cases to correct ordering */
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GIAJB");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 10, 10, "GIBJA");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "Giajb");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 10, 10, "Gibja");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GIAjb");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 10, 10, "GIbjA");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GiaJB");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 10, 10, "GiBJa");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GIaJb");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 10, 10, "GIbJa");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 10, 10, 10, 10, 0, "GiAjB");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 10, 10, "GiBjA");
    dpd_buf4_close(&G);

    if (params.ground) { /* otherwise, sort in x_Gibja */
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
      dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
      dpd_buf4_symm(&G);
      dpd_buf4_close(&G);
      dpd_buf4_init(&G1, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
      dpd_buf4_init(&G2, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
      dpd_buf4_symm2(&G1, &G2);
      dpd_buf4_close(&G2);
      dpd_buf4_sort(&G1, CC_GAMMA, rspq, 10, 10, "GiBJa");
      dpd_buf4_close(&G1);
    }
  }
  else if(params.ref == 2) { /** UHF **/

    /* G(ia,jb) <-- L(im,ae) T(jm,be) */
    dpd_buf4_init(&V, CC_MISC, 0, 20, 20, 20, 20, 0, "VIAJB");
    dpd_buf4_sort(&V, CC_MISC, rspq, 20, 20, "GIAJB");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 30, 30, 30, 30, 0, "Viajb");
    dpd_buf4_sort(&V, CC_MISC, rspq, 30, 30, "Giajb");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 30, 20, 30, 20, 0, "ViaJB");
    dpd_buf4_sort(&V, CC_MISC, rspq, 20, 30, "GIAjb");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 20, 30, 20, 30, 0, "VIAjb");
    dpd_buf4_sort(&V, CC_MISC, rspq, 30, 20, "GiaJB");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 24, 24, 24, 24, 0, "VIaJb");
    dpd_buf4_sort(&V, CC_MISC, rspq, 24, 24, "GIaJb");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, CC_MISC, 0, 27, 27, 27, 27, 0, "ViAjB");
    dpd_buf4_sort(&V, CC_MISC, rspq, 27, 27, "GiAjB");
    dpd_buf4_close(&V);

    /* G(IA,JB) <-- - L(IM,AE) T(J,E) T(M,B) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 21, 0, 21, 0, "Z(IM,AJ)");
    dpd_buf4_init(&L, CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 20, 21, 20, 21, 0, "Z(IB,AJ)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, prqs, 20, 21, "Z(IA,BJ)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(IA,BJ)");
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 20, 20, "Z(IA,JB)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 20, 20, 20, 20, 0, "Z(IA,JB)");
    dpd_buf4_init(&G, CC_MISC, 0, 20, 20, 20, 20, 0, "GIAJB");
    dpd_buf4_axpy(&Z1, &G, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&G);

    /* G(ia,jb) <-- - L(im,ae) T(j,e) T(m,b) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 31, 10, 31, 0, "Z(im,aj)");
    dpd_buf4_init(&L, CC_GLG, 0, 10, 15, 12, 17, 0, "Lijab");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 30, 31, 30, 31, 0, "Z(ib,aj)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, prqs, 30, 31, "Z(ia,bj)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(ia,bj)");
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 30, 30, "Z(ia,jb)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 30, 30, 30, 30, 0, "Z(ia,jb)");
    dpd_buf4_init(&G, CC_MISC, 0, 30, 30, 30, 30, 0, "Giajb");
    dpd_buf4_axpy(&Z1, &G, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&G);

    /* G(IA,jb) <-- - L(Im,Ae) T(j,e) T(m,b) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 22, 26, 22, 26, 0, "Z(Im,Aj)");
    dpd_buf4_init(&L, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 24, 26, 24, 26, 0, "Z(Ib,Aj)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, prqs, 20, 31, "Z(IA,bj)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 20, 31, 20, 31, 0, "Z(IA,bj)");
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 20, 30, "Z(IA,jb)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 20, 30, 20, 30, 0, "Z(IA,jb)");
    dpd_buf4_init(&G, CC_MISC, 0, 20, 30, 20, 30, 0, "GIAjb");
    dpd_buf4_axpy(&Z1, &G, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&G);

    /* G(ia,JB) <-- - L(iM,aE) T(J,E) T(M,B) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(iM,aJ)");
    dpd_buf4_init(&L, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 27, 25, 27, 25, 0, "Z(iB,aJ)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, prqs, 30, 21, "Z(ia,BJ)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 30, 21, 30, 21, 0, "Z(ia,BJ)");
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 30, 20, "Z(ia,JB)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 30, 20, 30, 20, 0, "Z(ia,JB)");
    dpd_buf4_init(&G, CC_MISC, 0, 30, 20, 30, 20, 0, "GiaJB");
    dpd_buf4_axpy(&Z1, &G, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&G);

    /* G(Ia,Jb) <-- - L(Im,Ea) T(J,E) T(m,b) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 22, 24, 22, 24, 0, "Z(Im,Ja)");
    dpd_buf4_init(&L, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_file2_close(&T1);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_buf4_init(&Z1, CC_TMP1, 0, 24, 24, 24, 24, 0, "Z(Ib,Ja)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, psrq, 24, 24, "Z(Ia,Jb)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Ia,Jb)");
    dpd_buf4_init(&G, CC_MISC, 0, 24, 24, 24, 24, 0, "GIaJb");
    dpd_buf4_axpy(&Z1, &G, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    /* G(iA,jB) <-- - L(iM,eA) T(j,e) T(M,B) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 23, 27, 23, 27, 0, "Z(iM,jA)");
    dpd_buf4_init(&L, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract244(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_file2_close(&T1);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z1, CC_TMP1, 0, 27, 27, 27, 27, 0, "Z(iB,jA)");
    dpd_contract424(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_file2_close(&T1);
    dpd_buf4_sort(&Z1, CC_TMP0, psrq, 27, 27, "Z(iA,jB)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(iA,jB)");
    dpd_buf4_init(&G, CC_MISC, 0, 27, 27, 27, 27, 0, "GiAjB");
    dpd_buf4_axpy(&Z1, &G, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1A);
    dpd_file2_mat_rd(&T1A);
    dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_mat_init(&T1B);
    dpd_file2_mat_rd(&T1B);
    dpd_file2_init(&L1A, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_mat_init(&L1A);
    dpd_file2_mat_rd(&L1A);
    dpd_file2_init(&L1B, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_mat_init(&L1B);
    dpd_file2_mat_rd(&L1B);

    /* G(IA,JB) <-- L(I,A) T(J,B) */
    dpd_buf4_init(&G, CC_MISC, 0, 20, 20, 20, 20, 0, "GIAJB");
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

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

      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);


    /* G(ia,jb) <-- L(i,a) T(j,b) */
    dpd_buf4_init(&G, CC_MISC, 0, 30, 30, 30, 30, 0, "Giajb");
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

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

      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    /* G(IA,jb) <-- L(I,A) T(j,b) */
    dpd_buf4_init(&G, CC_MISC, 0, 20, 30, 20, 30, 0, "GIAjb");
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

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

      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    /* G(ia,JB) <-- L(i,a) T(J,B) */
    dpd_buf4_init(&G, CC_MISC, 0, 30, 20, 30, 20, 0, "GiaJB");
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

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

      dpd_buf4_mat_irrep_wrt(&G, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_scm(&G, -1.0);
    dpd_buf4_close(&G);

    dpd_file2_mat_close(&L1A);
    dpd_file2_close(&L1A);
    dpd_file2_mat_close(&L1B);
    dpd_file2_close(&L1B);
    dpd_file2_mat_close(&T1A);
    dpd_file2_close(&T1A);
    dpd_file2_mat_close(&T1B);
    dpd_file2_close(&T1B);

    /* Sort all spin cases to correct ordering */
    dpd_buf4_init(&G, CC_MISC, 0, 20, 20, 20, 20, 0, "GIAJB");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 20, 20, "GIBJA");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 30, 30, 30, 30, 0, "Giajb");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 30, 30, "Gibja");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 20, 30, 20, 30, 0, "GIAjb");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 24, 27, "GIbjA");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 30, 20, 30, 20, 0, "GiaJB");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 27, 24, "GiBJa");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 24, 24, 24, 24, 0, "GIaJb");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 24, 24, "GIbJa");
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_MISC, 0, 27, 27, 27, 27, 0, "GiAjB");
    dpd_buf4_sort(&G, CC_GAMMA, psrq, 27, 27, "GiBjA");
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G1, CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
    dpd_buf4_init(&G2, CC_GAMMA, 0, 27, 24, 27, 24, 0, "GiBJa");
    dpd_buf4_symm2(&G1, &G2);
    dpd_buf4_close(&G2);
    dpd_buf4_sort(&G1, CC_GAMMA, rspq, 27, 24, "GiBJa");
    dpd_buf4_close(&G1);
  }
}

}} // namespace psi::ccdensity
