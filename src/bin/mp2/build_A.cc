/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_build_A(void);
void uhf_build_A(void);
void rhf_sf_build_A(void);
void uhf_sf_build_A(void);

void build_A(void)
{
  if(params.gradient) {
    if(params.ref == 0) rhf_sf_build_A();
    else if(params.ref == 2) uhf_sf_build_A();
  }
  else {
    if(params.ref == 0) rhf_build_A();
    else if(params.ref == 2) uhf_build_A();
  }
}

void rhf_build_A(void)
{
  dpdbuf4 Amat;  /* MO Hessian */
  dpdbuf4 D;     /* Two-electron integral */
  dpdbuf4 C;     /* Two-electron integral */
  dpdfile2 fIJ;  /* Occ-Occ Fock matrix */
  dpdfile2 fAB;  /* Vir-Vir Fock matrix */
  int h, nirreps;
  int i, j, a, b, ai, bj;
  int I, J, A, B;
  int Asym, Bsym, Isym, Jsym;
  
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&D, CC_MISC, rpsq, 11, 11, "A(AI,BJ)");
  dpd_buf4_close(&D);

  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_buf4_scm(&Amat, 4.0);
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort_axpy(&D, CC_MISC, sprq, 11, 11, "A(AI,BJ)", -1.0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_sort_axpy(&C, CC_MISC, qpsr, 11, 11, "A(AI,BJ)", -1.0);
  dpd_buf4_close(&C);
  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_buf4_close(&Amat);

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fAB);

  nirreps = mo.nirreps;

  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Amat, h);
    dpd_buf4_mat_irrep_rd(&Amat, h);
    for(ai=0; ai < Amat.params->rowtot[h]; ai++) {
      a = Amat.params->roworb[h][ai][0];
      i = Amat.params->roworb[h][ai][1];
      A = fAB.params->rowidx[a];
      I = fIJ.params->rowidx[i];
      Asym = fAB.params->psym[a];
      Isym = fIJ.params->psym[i];
      for(bj=0; bj < Amat.params->coltot[h]; bj++) {
        b = Amat.params->colorb[h][bj][0];
        j = Amat.params->colorb[h][bj][1];
        B = fAB.params->colidx[b];
        J = fIJ.params->colidx[j];
        Bsym = fAB.params->qsym[b];
        Jsym = fIJ.params->qsym[j];
        if((I==J)&&(Asym==Bsym)) Amat.matrix[h][ai][bj] += fAB.matrix[Asym][A][B];
        if((A==B)&&(Isym==Jsym)) Amat.matrix[h][ai][bj] -= fIJ.matrix[Isym][I][J];
      }
    }
    dpd_buf4_mat_irrep_wrt(&Amat, h);
    dpd_buf4_mat_irrep_close(&Amat, h);
  }
  dpd_buf4_close(&Amat);

  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fAB);
  dpd_file2_mat_close(&fIJ);
  dpd_file2_close(&fIJ);
}

void uhf_build_A(void)
{

}

void rhf_sf_build_A(void)
{
  int h, nirreps, e, m, a, i, em, ai, E, M, A, I;
  int Esym, Msym, Asym, Isym;
  int *virtpi, *occpi, *occ_off, *vir_off;
  int *qt_occ, *qt_vir; /* Spatial orbital translators */
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia;
  dpdbuf4 Amat, Amat2, D, C;

  nirreps = mo.nirreps;
  occpi = mo.occpi; 
  virtpi = mo.virtpi;
  occ_off = mo.occ_off; 
  vir_off = mo.vir_off;
  qt_occ = mo.qt_occ; 
  qt_vir = mo.qt_vir;

  /* Two-electron integral contributions */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&D, CC_MISC, rpsq, 11, 11, "A(EM,AI)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(EM,AI)");
  dpd_buf4_sort(&Amat, CC_TMP0, psrq, 11, 11, "D <im|ea> (ei,am)");
  dpd_buf4_scm(&Amat, 2.0);
  dpd_buf4_copy(&Amat, CC_TMP0, "A(EM,ai)");
  dpd_buf4_init(&D, CC_TMP0, 0, 11, 11, 11, 11, 0, "D <im|ea> (ei,am)");
  dpd_buf4_axpy(&D, &Amat, -1.0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_sort(&C, CC_TMP0, qpsr, 11, 11, "C <ai|bj>");
  dpd_buf4_close(&C);
  dpd_buf4_init(&C, CC_TMP0, 0, 11, 11, 11, 11, 0, "C <ai|bj>");
  dpd_buf4_axpy(&C, &Amat, -1.0);
  dpd_buf4_close(&C);
  dpd_buf4_copy(&Amat, CC_TMP0, "A(em,ai)");
  dpd_buf4_close(&Amat);

  /* Fock matrix contributions */
  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
  dpd_file2_mat_init(&fab);
  dpd_file2_mat_rd(&fab);
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_mat_init(&fIA);
  dpd_file2_mat_rd(&fIA);
  dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
  dpd_file2_mat_init(&fia);
  dpd_file2_mat_rd(&fia);

  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(EM,AI)");
  
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Amat, h);
    dpd_buf4_mat_irrep_rd(&Amat, h);

    for(em=0; em < Amat.params->rowtot[h]; em++) {
      e = Amat.params->roworb[h][em][0];
      m = Amat.params->roworb[h][em][1];
      E = fAB.params->rowidx[e];
      M = fIJ.params->rowidx[m];
      Esym = fAB.params->psym[e];
      Msym = fIJ.params->psym[m];
      for(ai=0; ai < Amat.params->coltot[h]; ai++) {
	a = Amat.params->colorb[h][ai][0];
	i = Amat.params->colorb[h][ai][1];
	A = fAB.params->colidx[a];
	I = fIJ.params->colidx[i];
	Asym = fAB.params->qsym[a];
	Isym = fIJ.params->qsym[i];

	if((M==I) && (Esym==Asym))
	  Amat.matrix[h][em][ai] += fAB.matrix[Esym][E][A];
	if((E==A) && (Msym==Isym))
	  Amat.matrix[h][em][ai] -= fIJ.matrix[Msym][M][I];
      }
    }

    dpd_buf4_mat_irrep_wrt(&Amat, h);
    dpd_buf4_mat_irrep_close(&Amat, h);
  }

  dpd_buf4_close(&Amat);

  dpd_buf4_init(&Amat, CC_TMP0, 0, 11, 11, 11, 11, 0, "A(em,ai)");
  
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Amat, h);
    dpd_buf4_mat_irrep_rd(&Amat, h);

    for(em=0; em < Amat.params->rowtot[h]; em++) {
      e = Amat.params->roworb[h][em][0];
      m = Amat.params->roworb[h][em][1];
      E = fab.params->rowidx[e];
      M = fij.params->rowidx[m];
      Esym = fab.params->psym[e];
      Msym = fij.params->psym[m];
      for(ai=0; ai < Amat.params->coltot[h]; ai++) {
	a = Amat.params->colorb[h][ai][0];
	i = Amat.params->colorb[h][ai][1];
	A = fab.params->colidx[a];
	I = fij.params->colidx[i];
	Asym = fab.params->qsym[a];
	Isym = fij.params->qsym[i];

	if((M==I) && (Esym==Asym))
	  Amat.matrix[h][em][ai] += fab.matrix[Esym][E][A];
	if((E==A) && (Msym==Isym))
	  Amat.matrix[h][em][ai] -= fij.matrix[Msym][M][I];
      }
    }

    dpd_buf4_mat_irrep_wrt(&Amat, h);
    dpd_buf4_mat_irrep_close(&Amat, h);
  }

  dpd_buf4_close(&Amat);

  dpd_buf4_init(&Amat, CC_TMP0, 0, 11, 11, 11, 11, 0, "A(EM,ai)");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Amat, h);
    dpd_buf4_mat_irrep_rd(&Amat, h);

    for(em=0; em < Amat.params->rowtot[h]; em++) {
      e = Amat.params->roworb[h][em][0];
      m = Amat.params->roworb[h][em][1];
      Esym = Amat.params->psym[e];
      Msym = Amat.params->qsym[m];
      E = e - vir_off[Esym];
      M = m - occ_off[Msym];
      for(ai=0; ai < Amat.params->coltot[h]; ai++) {
	a = Amat.params->colorb[h][ai][0];
	i = Amat.params->colorb[h][ai][1];
	Asym = Amat.params->rsym[a];
	Isym = Amat.params->ssym[i];
	A = a - vir_off[Asym];
	I = i - occ_off[Isym];

	/* This comparison is somewhat tricky.  The algebraic
	   expression for the Fock matrix shift here is:

	   A(EM,ai) += delta(M,a) f(E,i)(beta)

	   The Kronecker Delta is actually a comparison between
	   the *spatial* orbitals associated with M, and A.
	   Hence we have to compare the spatial orbital
	   translation of the the two absolute orbital indices. */
	if((qt_occ[m] == qt_vir[a]) && (Esym==Isym))
	  Amat.matrix[h][em][ai] += fia.matrix[Isym][I][E];
      }
    }

    dpd_buf4_mat_irrep_wrt(&Amat, h);
    dpd_buf4_mat_irrep_close(&Amat, h);
  }
  dpd_buf4_sort(&Amat, CC_TMP0, rspq, 11, 11, "A(em,AI)");
  dpd_buf4_close(&Amat);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_close(&fij);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fab);
  dpd_file2_mat_close(&fIA);
  dpd_file2_close(&fIA);
  dpd_file2_mat_close(&fia);
  dpd_file2_close(&fia);

  /* Now sum all three A-matrix components and divide by 2 */
  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(EM,AI)");
  dpd_buf4_init(&Amat2, CC_TMP0, 0, 11, 11, 11, 11, 0, "A(em,ai)");
  dpd_buf4_axpy(&Amat2, &Amat, 1.0);
  dpd_buf4_close(&Amat2);
  dpd_buf4_init(&Amat2, CC_TMP0, 0, 11, 11, 11, 11, 0, "A(EM,ai)");
  dpd_buf4_axpy(&Amat2, &Amat, 1.0);
  dpd_buf4_close(&Amat2);
  dpd_buf4_init(&Amat2, CC_TMP0, 0, 11, 11, 11, 11, 0, "A(em,AI)");
  dpd_buf4_axpy(&Amat2, &Amat, 1.0);
  dpd_buf4_close(&Amat2);
  dpd_buf4_scm(&Amat, 0.5);
  dpd_buf4_close(&Amat);
}

void uhf_sf_build_A(void)
{

}

}} /* End namespace */
